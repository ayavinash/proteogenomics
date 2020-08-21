import pandas as pd
from argparse import ArgumentParser
import os,sys,math
from useful_functions import get_data_from_fasta,chunks,assign_cores,fastaToDict,map_refseq_chr_accessions,complement_data
from multiprocessing import Process,Manager
from peptide_to_protein_mapping import translate_on_all_cores

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq

def get_exon_skip_proteins(database):
    fasta_data = get_data_from_fasta(database)
    id_split = fasta_data.fasta_id.str.split("_",expand=True)
    id_split.rename(columns ={0:"trans_id"},inplace=True)
    del id_split[1]
    del id_split[2] 
    fasta_data = pd.concat([fasta_data,id_split],axis=1)
    fasta_data["skipped_exon"] = fasta_data.description.str.split(" ").apply(lambda x: x[-1])
    fasta_data.skipped_exon = fasta_data.skipped_exon.apply(int)
    return fasta_data    

def remove_row(skip_data,trans_data):
    skip_cds_list = skip_data.cds_no.tolist()
    frames = {x:trans_data[trans_data.cds_no!=x] for x in skip_cds_list}
    return pd.concat(frames.values(),ignore_index=True)


def remove_each_row(trans_data,cds_no):
    trans_data = trans_data[trans_data.cds_no!=cds_no]
    trans_data.skipped_exon = cds_no
    return trans_data

def modify_gff(gff_file_path,cores):
    data=pd.read_csv(gff_file_path,header=None,comment="#",sep="\t",\
		names=["chr","source","feature","start","end","score","strand","frame","attributes"])
    
    data = data[data.feature=="CDS"]
    data = data[data.chr.str.startswith("chr")]
    data.attributes= data.attributes+";"
    data["parent_id"] = data.attributes.str.extract('Parent=(.*?);')
    data["feature_id"] = data.attributes.str.extract('ID=(.*?);')
    data["trans_id"] = data.attributes.str.extract(';transcript_id=(.*?);')
    data["exon_no"] = data.attributes.str.extract(';exon_number=(.*?);')
    data["gene_name"] = data.attributes.str.extract('gene_name=(.*?);')
    del data["attributes"]
    data.exon_no = data.exon_no.apply(int)
    data["cds_no"] =1
    data.cds_no = data.groupby("parent_id")["cds_no"].cumsum()
    data.set_index("parent_id",inplace=True)
    data["total_cds"] = data.groupby("parent_id")["trans_id"].count()
    data = data[data.total_cds>2]
    data.reset_index(inplace=True)
    data["frame_below"] = data.frame.shift(-1)
    data["parent_id_below"] = data.parent_id.shift(-1)
    data["frame_same"] = data.frame==data.frame_below
    data["parent_same"] = data.parent_id==data.parent_id_below
    data["to_skip"] = data.parent_same==data.frame_same
    #data["skipped_exon"] = np.nan
    tmp = data[(data.to_skip) & (data.parent_same) & (data.frame_same) & (data.cds_no!=1)]
    skip_data = tmp.copy(deep=True)
    skip_data["loop_no"] =1
    #skip_data.set_index("parent_id",inplace=True)
    skip_data.loop_no = skip_data.groupby("parent_id")["loop_no"].cumsum()
    #skip_data.reset_index(inplace=True)
    loops = set(skip_data.loop_no)
    frames = {}
    for loop in loops:
        #print("Looping....",loop)
        tmp = skip_data[skip_data.loop_no==loop]
        transcripts = tmp.parent_id
        
        #print(transcripts.head())
        tmp_data = data.drop(transcripts.index)
        tmp_data = tmp_data[tmp_data.parent_id.isin(transcripts)]
        print("Looping....",loop,len(tmp_data),len(transcripts))
        tmp_data["loop_no"] = loop
        tmp_data["skipped_exon"] = tmp_data.parent_id.map(tmp.set_index("parent_id")["cds_no"].to_dict())
        frames[loop] = tmp_data
    new_skip = pd.concat(frames.values(),ignore_index=True)
    #new_skip["skipped_exon"] = skip_data["cds_no"]
    #new_skip = new_skip[new_skip.skipped_exon.notnull()]

    return new_skip



    #data["exons_to_skip"] = 

    grouped = data.set_index("parent_id").groupby("parent_id")
#   skip_grouped = skip_data.groupby("parent_id")

    data_series = skip_data.apply(lambda x: remove_each_row(grouped.get_group(x["parent_id"]),x["cds_no"]),axis=1)
    skip = pd.concat(data_series.tolist())

    index_list = skip_data.index.tolist()
    index_lists = list(chunks(index_list, math.ceil(len(index_list)/cores)))

    all_args_list = []
    manager = Manager()
    return_dict = manager.dict()    
    for index_list in index_lists:
        print("Distributing....")
        tmp = skip_data.loc[index_list]
        cpu_skip_data = tmp.copy(deep=True)
        tmp = data[data.parent_id.isin(set(cpu_skip_data.parent_id))]
        cpu_data = tmp.copy(deep=True)
        arg_list =(cpu_skip_data,cpu_data)
        all_args_list.append(arg_list)
    p_list = []
    for cpu in range(len(all_args_list)):
        p = Process(target = modify_gff_each_core,args =(all_args_list[cpu][0],all_args_list[cpu][1],return_dict))
        p_list.append(p)
    for p in p_list:
        p.start()
    for p in p_list:
        p.join()
    #print(return_dict)
    if len(return_dict)>1:
        new_gff = pd.concat(return_dict.values(),ignore_index=True)   
    elif len(return_dict)==1:
        new_gff = return_dict.values()[0]
    else:
        print("No hits....")
        sys.exit()
    return new_gff
    

def modify_gff_each_core(skip_data,data,return_dict):
    grouped = data.groupby("parent_id")
    print("Generating-----")
    data_series = skip_data.apply(lambda x: remove_each_row(grouped.get_group(x["parent_id"]),x["cds_no"]),axis=1)
    new_gff = pd.concat(data_series.tolist())
    return_dict[str(os.getpid())] = new_gff

def change_ids(cds_data):
    cds_data.skipped_exon = cds_data.skipped_exon.apply(int).apply(str)
    cds_data.feature_id = "ID="+cds_data.feature_id
    cds_data.parent_id = "Parent="+cds_data.parent_id+"_NovIso_"+cds_data.skipped_exon
    cds_data.trans_id = "transcript_id="+cds_data.trans_id+"_NovIso_"+cds_data.skipped_exon
    cds_data.gene_name = "gene_name="+cds_data.gene_name
    cds_data["attributes"] = cds_data.feature_id.str.cat([cds_data.parent_id,cds_data.trans_id,cds_data.gene_name],sep=";")
    cds_data = cds_data[["chr","source","feature","start","end","score","strand","frame","attributes"]]

    return cds_data

def generate_fasta_from_gff(data,genome_file_path,frames,cpus,primary_accessions,annotation_source):
    print("Reading gtf file..")
#    data=pd.read_csv(gff_file_path,header=None,comment="#",sep="\t",\
#		names=["chr","source","feature","start","end","score","strand","frame","attributes"])
#    if annotation_source=="gencode":
    data = data[data.feature == "CDS"]
    #data["sequence_id"]= data.attributes.str.split(";").apply(lambda x: x[0].split("ID=")[1])
    print("Loading genome fasta...")
    genome_dict = fastaToDict(genome_file_path)
    print("Done loading...")
    ##remove ids that are not present in fasta file..
    chrs_in_gff = set(data.chr)
    chrs_in_genome = set(genome_dict.keys())
    chrs_not_present = chrs_in_gff.difference(chrs_in_genome)
    if chrs_not_present:
        print("Following gff chromsome ids could not be found in genome file...")
        print(chrs_not_present)
        data = data[data.chr.isin(set(genome_dict.keys()))]
    data["sequence"] = data.apply(lambda x: genome_dict[x["chr"]][x["start"]-1:x["end"]].seq,axis=1)
    data.sequence[data.strand=="-"] =data.sequence[data.strand=="-"].apply(lambda x: x.reverse_complement())
    if frames == 6:
        data = complement_data(data)
    else:
        data["seq_type"] = "original"
    #data = generate_frames_data(data)
    if annotation_source in ("refseq","gnomon"):
        data = map_refseq_chr_accessions(data,primary_accessions)
    return data


def get_transcript_from_exon_or_cds(data):
    data["parent_id"] = data.attributes.str.split(";").apply(lambda x: x[1].split("Parent=")[1])
    data.sequence= data.sequence.apply(lambda x: str(x))
    tmp = data.drop_duplicates("parent_id").set_index("parent_id")
    grouped = data.groupby("parent_id")
    trans_series = grouped["sequence"].apply(lambda x: x.str.cat())
    trans_data = pd.DataFrame({"sequence":trans_series,"chr":tmp.chr,"frame":tmp.frame},columns=["sequence","chr","frame"])
    trans_data["frame"] = trans_data.frame.apply(int)
    trans_data.sequence = trans_data.sequence.apply(lambda x: Seq(x))
    trans_data["sequence_id"] = trans_data.index
    return trans_data

def get_arguments():
    parser= ArgumentParser()
    #parser.add_argument("-peptides",required=True,help="full path to peptides file")
    #parser.add_argument("-fasta", help="full path to exon skip fasta file")
    parser.add_argument("-gff", help="full path to exon skip database file")
    parser.add_argument("-out",required=True,help ="full path to output dir")
    parser.add_argument("-cores",default=1,type=int,help ="full path output gff file")
    parser.add_argument("-genome",help="path to genome fasta file.")
    parser.add_argument("-frames",default=3,type=int,choices =[3,6], help ="full path output gff file")
    parser.add_argument("-accessions",help="path to accessions file.")
    parser.add_argument("-annotation",required=True,choices= ["refseq","gencode","gnomon"],help="annotation source")
    args= parser.parse_args()
    return args

def main():
    args = get_arguments()
    print("Reading fasta...")
    #exon_skipped_fasta = get_exon_skip_proteins(args.fasta)
    #exon_skipped_fasta = exon_skipped_fasta
    print("reading gff...")
    gff_base_name = os.path.basename(args.gff)
    gff_name,ext = os.path.splitext(gff_base_name)
    out_gff_name = gff_name+"_single_exon_skipped.gff3"
    out_fasta_name = gff_name+"_single_exon_skipped_proteins.fasta"
    out_gff_file_path = os.path.join(args.out,out_gff_name)
    out_fasta_file_path = os.path.join(args.out,out_fasta_name)
    print(out_gff_file_path)
    print(out_fasta_file_path)
    cores =assign_cores(args.cores)
    if not os.path.exists(out_gff_file_path):
        exonskiped_gff = modify_gff(args.gff,cores)
        exonskiped_gff  = change_ids(exonskiped_gff)
        exonskiped_gff.to_csv(out_gff_file_path,sep="\t",header=False,index=False)    
    print("Loading exon skipped gff file..")
    exonskiped_gff=pd.read_csv(out_gff_file_path,header=None,comment="#",sep="\t",\
		    names=["chr","source","feature","start","end","score","strand","frame","attributes"])
    #exonskiped_gff = exonskiped_gff.loc[:1000]

    #
    data= generate_fasta_from_gff(exonskiped_gff,args.genome,args.frames,args.cores,args.accessions,args.annotation)
    trans_data = get_transcript_from_exon_or_cds(data)
    #print(trans_data.sequence.loc["ENST00000263741.11_NovIso_3"])
    protein_data = translate_on_all_cores(trans_data,cores)
    protein_data.sequence = protein_data.apply(lambda x: SeqRecord(Seq(x["sequence"]),id=x["sequence_id"]+"|"+x["sequence_id"]),axis=1)
    #print(protein_data.sequence.iloc[0])
    #print(protein_data.sequence)
    #print(protein_data.sequence.loc["ENST00000263741.11_NovIso_3"])
    recs = protein_data.sequence.to_dict().values()
    SeqIO.write(recs, out_fasta_file_path, "fasta")


if __name__=="__main__":
    main()

