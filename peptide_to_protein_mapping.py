
### Takes peptides and finds all hits in all proteins in the database##
### Author: Avinash Yadav
### PhD student, Scuola Normale Superiore
### contact: ayavinash@gmail.com

from argparse import ArgumentParser
import pandas as pd
import os,sys,math,time
from multiprocessing import Process,cpu_count,Manager
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from useful_functions import *

def generate_data_from_fasta(database_file_path,annotation_source,db_type,gff_file_path,primary_accessions,frames,cpus):
    split_point =get_split_point(annotation_source,db_type)
    fasta_data = get_data_from_fasta(database_file_path)
    fasta_data["sequence_id"] = fasta_data.fasta_id.apply(lambda x: x.split("|")[split_point])
    del fasta_data["fasta_id"]
    del fasta_data["description"]

    if (db_type=="protein") and (annotation_source=="gencode"):
        fasta_data["sequence"] = fasta_data.sequence.apply(lambda x: x[1:] if x.startswith("X") else x)
    if (db_type=="transcript"):
        fasta_data = get_strand_chr_for_transcripts(fasta_data,gff_file_path,annotation_source,primary_accessions,frames,cpus)
    #print(fasta_data.sequence.head())
    return fasta_data

def generate_data_from_genome(genome_file_path,annotation_source,gff_file_path,primary_accessions,frames,cpus):
    print("Reading gtf file..")
    data=pd.read_csv(gff_file_path,header=None,comment="#",sep="\t",\
		names=["chr","source","feature","start","end","score","strand","frame","attributes"])
#    if annotation_source=="gencode":
    data = data[data.feature == "gene"]
    data["sequence_id"]= data.attributes.str.split(";").apply(lambda x: x[0].split("ID=")[1])
    genome_dict = fastaToDict(genome_file_path)
    ##remove ids that are not present in fasta file..
    chrs_in_gff = set(data.chr)
    chrs_in_genome = set(genome_dict.keys())
    chrs_not_present = chrs_in_gff.difference(chrs_in_genome)
    if chrs_not_present:
        print("Following gff chromsome ids could not be found in genome file...")
        print(chrs_not_present)
        data = data[data.chr.isin(set(genome_dict.keys()))]
    data["sequence"] = data.apply(lambda x: genome_dict[x["chr"]][x["start"]-1:x["end"]],axis=1)
    data.sequence[data.strand=="-"] =data.sequence[data.strand=="-"].apply(lambda x: x.reverse_complement())
    del_list =["source","feature","score","frame","attributes"]
    for col in del_list:
        del data[col]
    if frames == 6:
        data = complement_data(data)
    else:
        data["seq_type"] = "original"
    data = generate_frames_data(data)
    if annotation_source in ("refseq","gnomon"):
        data = map_refseq_chr_accessions(data,primary_accessions)
    data = translate_on_all_cores(data,cpus)
    return data    


def get_strand_chr_for_transcripts(fasta_data,gff_file_path,annotation_source,primary_accessions,frames,cpus):
    trans_data = get_trans_data(gff_file_path,annotation_source)
   
    fasta_data["strand"] = fasta_data.sequence_id.map(trans_data.strand.to_dict())
    fasta_data = fasta_data[fasta_data.strand.notnull()]
    fasta_data["chr"] = fasta_data.sequence_id.map(trans_data.chr.to_dict())

    del trans_data
    
    if frames == 6:
        fasta_data = complement_data(fasta_data)
    else:
        fasta_data["seq_type"] = "original"

    fasta_data = generate_frames_data(fasta_data)
    if annotation_source in ("refseq","gnomon"):
        fasta_data = map_refseq_chr_accessions(fasta_data,primary_accessions)
    fasta_data = translate_on_all_cores(fasta_data,cpus)
    #fasta_data = translate(fasta_data)
    
    return fasta_data

def translate_on_all_cores(data,cpus):
    manager = Manager()
    return_dict = manager.dict()    
    index_list = data.index
    index_lists = list(chunks(index_list, math.ceil(len(index_list)/cpus)))
    cpu_data_frame_list = []
    for cpu in range(cpus):
        #print(len(pep_index_lists[cpu]))
        small_data = data.loc[index_lists[cpu]]
        cpu_data_frame = small_data.copy(deep=True)
        cpu_data_frame_list.append(cpu_data_frame)
    p_list = []
    #print(pep_index_lists)
    for cpu in range(cpus):
        p = Process(target = translate,args =(cpu_data_frame_list[cpu],return_dict))
        p_list.append(p)
    for p in p_list:
        p.start()
    for p in p_list:
        p.join()
    translated_data = pd.concat(return_dict.values(),ignore_index=True)   
    return translated_data    



def translate(fasta_data,return_dict):
    #print(fasta_data)
    #print(fasta_data.iloc[0])
    print("Translating on process...",os.getpid())
    chr_data = fasta_data[fasta_data.chr=="chrM"]
    other_data =fasta_data[fasta_data.chr!="chrM"]
    #print(chr_data)
    frames =[]
    if not chr_data.empty:
        
        chr_data["sequence"] = chr_data.apply(lambda x: x["sequence"][x["frame"]:].translate(table =2),axis = 1)
        frames.append(chr_data)
    if not other_data.empty:
        print("Translating....")
        other_data["sequence"] = other_data.apply(lambda x: x["sequence"][x["frame"]:].translate(table =1),axis = 1)
        frames.append(other_data)
    if not frames:
        print("No sequences to map peptides to..on process..",os.getpid())
        fasta_data = pd.DataFrame([])

    else:
        print("Length of frames...is",len(frames))
        fasta_data = pd.concat(frames,ignore_index=True)
    if not fasta_data.empty:
        fasta_data.sequence = fasta_data.sequence.apply(str)
    return_dict[str(os.getpid())] = fasta_data


def get_trans_data(gff_file_path,annotation_source):
    print("Reading gtf file..")
    data=pd.read_csv(gff_file_path,header=None,comment="#",sep="\t",\
		names=["chr","source","feature","start","end","score","strand","frame","attributes"])
    if annotation_source=="gencode":
        data = data[data.feature == "transcript"]
        data["trans_id"]= data.attributes.str.split(";").apply(lambda x: x[0].split("ID=")[1])
    elif annotation_source in ("refseq","gnomon"):
        data = data[(data.feature.str.contains("RNA")) | (data.feature.str.contains("transcript")) ]
        data["trans_id"] = data.attributes.str.extract('transcript_id=(.*)')
        data = data[data.trans_id.notnull()]
        if annotation_source=="gnomon":
            data["trans_id"] = data.trans_id.str.split("GNOMON:").apply(lambda x: x[1])
    else:
        print(annotation_source,"is not implemented...Exiting...")
        sys.exit()
    data.set_index("trans_id",inplace=True)
    return data

def substring_indexes(substring, string):
    last_found = -1  
    while True:
        last_found = string.find(substring, last_found + 1)
        if last_found == -1:  
            break  
        yield last_found

def find_peptide_hits(data, pep_data, out_file_path,return_dict):
    filename, ext = os.path.splitext(out_file_path) 
    out_file_path = filename+"_"+str(os.getpid())+ext
    #print(filename)
    #print(out_file_path)
    
    out = pd.concat(list(pep_data.apply(lambda x:row_function(data,x["pep_seq"],x["peptide"]),axis=1)),ignore_index=True)

    if not out.empty:   
        out["pep_hits"] = out.apply(lambda x: list(substring_indexes(x["pep_seq"],x["sequence"])),axis=1)
        out["pep_len"] = out.pep_seq.apply(len)
        out["aa_res_before"] = out.apply(lambda x:[x["sequence"][y-1] if y else "-" for y in x["pep_hits"]],axis=1)
        out["aa_res_after"] = out.apply(lambda x:[x["sequence"][y+x["pep_len"]] if len(x["sequence"])>y+x["pep_len"] else "-" for y in x["pep_hits"]],axis=1)

        #out.to_csv(out_file_path,index=False)
     
    else:
        print("Could not find any peptide hits....")
    del out["sequence"]
    return_dict[str(os.getpid())]= out



def row_function(data,pep,original_pep):
    pep_data = data[data.sequence.str.contains(pep,regex=False)]
    pep_data["pep_seq"] = pep
    pep_data["peptide"] = original_pep
    return pep_data


def check_arguments(args):
    if not (os.path.isfile(args.peptides)):
        sys.exit("Input peptides file not valid...!Exiting...")
    if args.database:
        if not ((os.path.isfile(args.database)) and (os.path.splitext(args.database)[1] in (".fa",".fasta",".fsa"))):
            sys.exit("Input database not valid ...! Exiting...")
    if args.gff:
        if not ((os.path.isfile(args.gff)) and (os.path.splitext(args.gff)[1] == ".gff3")):
            sys.exit("Input GFF3 file not valid ...! Exiting...")
    if args.accessions:
        if not os.path.isfile(args.accessions):
            sys.exit("Input Accessions file does not exists ...! Exiting...")
    if args.exDirs:
        for exDir in args.exDirs:
            if not os.path.isdir(exDir):
                print("directory...",exDir,"...to exclude peptides does not exists..! exiting...")
                sys.exit()
    if not os.path.exists(args.out):
        os.makedirs(args.out)

#    os.makedirs(args.out)
    
    log_file_path = os.path.join(args.out,"log.log")
    log_file = open(log_file_path,"a")
    log_file.writelines(sys.executable+" "+" ".join(sys.argv)+"\n")
    log_file.close()
    #print(args.exDirs)
    #sys.exit()

    if not args.outfile:
        pep_file_name = os.path.splitext(os.path.basename(os.path.basename(args.peptides)))[0]
        if args.dbtype=="gene":
            db_name=args.genome
        else:
            db_name = args.database

        database_file_name = os.path.splitext(os.path.basename(os.path.basename(db_name)))[0]
        out_file_name = pep_file_name+"_"+database_file_name+".csv"
        out_file_path = os.path.join(os.path.basename(args.out),out_file_name)
    else:
        out_file_path = os.path.join(args.out,args.outfile)
    if os.path.exists(out_file_path):
        print("output file",out_file_path,"...exists...Exiting.")
    return out_file_path


def get_peptides_from_fasta(pep_file_path):
    pep_data = get_data_from_fasta(pep_file_path)
    del pep_data["fasta_id"]
    del pep_data["description"]
    pep_data.rename(columns ={"sequence":"pep_seq"},inplace=True)
    pep_data.pep_seq = pep_data.pep_seq.apply(str)
    return pep_data

def get_peptides(pep_file_path,header):
    ext = os.path.splitext(pep_file_path)[1]
    fasta_exts = (".fa",".fsa",".fasta")
    if ext in fasta_exts:
        pep_data = get_peptides_from_fasta(pep_file_path)
    elif ext ==".csv":
        if not header:
            pep_data = pd.read_csv(pep_file_path,header=None,names=["pep_seq"])
        else:
            pep_data = pd.read_csv(pep_file_path)
            pep_data = pep_data[["pep_seq"]]
    #print(pep_data.pep_seq)
    pep_data = pep_data.drop_duplicates("pep_seq")
    pep_data["peptide"] = pep_data.pep_seq
    pep_data.pep_seq = pep_data.pep_seq.str.replace("L","I")
    #pep_data = pep_data.iloc[:1000,:]
    print(pep_data)
    return pep_data


def distribute_peptides(cores,pep_data,data_from_fasta,out_file_path):
    print("Finding peptide hits...")
    manager = Manager()
    return_dict = manager.dict()    
    cpus = cores
    pep_index_list = pep_data.index
    print("Total peptides to map...",len(pep_index_list))
    pep_index_lists = list(chunks(pep_index_list, math.ceil(len(pep_index_list)/cpus)))
    #print(pep_index_lists)
    cpu_data_frame_list = []
    for cpu in range(cpus):
        #print(len(pep_index_lists[cpu]))
        cpu_data_frame = data_from_fasta.copy(deep=True)
        cpu_data_frame_list.append(cpu_data_frame)
    p_list = []
    #print(pep_index_lists)
    for cpu in range(cpus):
        print(cpu)
        print(len(cpu_data_frame_list))
        p = Process(target = find_peptide_hits,args =(cpu_data_frame_list[cpu],pep_data.loc[pep_index_lists[cpu],:],out_file_path,return_dict))
        p_list.append(p)
    for p in p_list:
        p.start()
    for p in p_list:
        p.join()
    map_data = pd.concat(return_dict.values(),ignore_index=True)
    map_data.pep_hits = map_data.pep_hits.apply(lambda x: ";".join([str(y) for y in x]))
    map_data.aa_res_before = map_data.aa_res_before.apply(lambda x: ";".join(x))
    map_data.aa_res_after = map_data.aa_res_after.apply(lambda x: ";".join(x))
    map_data.to_csv(out_file_path,index=False)
    


def distribute_proteins(cores,pep_data,data_from_fasta,out_file_path):
    print("Finding peptide hits...")
    cpus = cores
    data_index_list = data_from_fasta.index
    data_index_lists = list(chunks(data_index_list, math.ceil(len(data_index_list)/cpus)))
    all_args_list =[]
    for cpu in range(cpus):
        args_list =[]
        small_data = data_from_fasta.loc[data_index_lists[cpu]]
        passing_fasta_data = small_data.copy(deep=True)
        passing_pep_data = pep_data.copy(deep=True)
        args_list =[passing_fasta_data,passing_pep_data,out_file_path]
        all_args_list.append(tuple(args_list))
    del data_from_fasta
    p_list = []
    #print(pep_index_lists)
    for cpu in range(cpus):
        p = Process(target = find_peptide_hits,args =all_args_list[cpu])
        p_list.append(p)
    for p in p_list:
        p.start()

def get_peptides_to_exclude(exclude_dir):
    exclude_files = os.listdir(exclude_dir)
    #print(exclude_files)
    exclude_files = [x for x in exclude_files if os.path.splitext(x)[1]==".csv"]
    frames = {}
    data = pd.DataFrame([])
    for exclude_file in exclude_files:
        exclude_file_path = os.path.join(exclude_dir,exclude_file)
        data = pd.read_csv(exclude_file_path)
        frames[exclude_file] = data
    
    if len(frames)>1:
        data = pd.concat(frames.values(),ignore_index=True)
    
    return data

def walk_in_exclude_dirs(exclude_dirs):
    frames ={}
    exclude_dir_paths = exclude_dirs
    #print(exclude_dir_paths)

    for exclude_dir_path in exclude_dir_paths:
        current_exclude_dir= exclude_dir_path
        #print(current_exclude_dir)
        data = get_peptides_to_exclude(current_exclude_dir)
        #print(data)
        if not data.empty:
            frames[current_exclude_dir] = data
            print("Excluding..",len(set(data.peptide)),".from..",current_exclude_dir)
    if len(frames)>1:
        data = pd.concat(frames.values(),ignore_index=True)
    return data


def main():
    args = get_arguments()
    out_file_path = check_arguments(args)
    #os.makedirs(args.out)
    pep_data = get_peptides(args.peptides,args.header)
    print(len(pep_data))
    if args.exDirs:
        exclude_data = walk_in_exclude_dirs(args.exDirs)
        #exclude_data = get_peptides_to_exclude(args.exclude)
        #print(exclude_data.columns)
        print("Excluding...",len(set(exclude_data.pep_seq)),"..peptides from main file..!")
        pep_data = pep_data[~pep_data.peptide.isin(set(exclude_data.peptide))]

    #print(pep_data)
    if not args.cores:
        cores = cpu_count()-2
    else:
        cores = assign_cores(args.cores)
    if args.dbtype=="gene":
        data_from_fasta = generate_data_from_genome(args.genome,args.annotation,args.gff,args.accessions,args.frames,cores)
    else:
        data_from_fasta = generate_data_from_fasta(args.database,args.annotation,args.dbtype,args.gff,args.accessions,args.frames,cores)
    data_from_fasta.sequence=data_from_fasta.sequence.apply(str)
    #print(data_from_fasta.sequence.head())
    data_from_fasta["sequence"] = data_from_fasta.sequence.str.replace("L","I")
    
    #print(data_from_fasta)
    #sys.exit()
    #print(data_from_fasta.head())
#    if not args.dbtype=="gene":
#        distribute_proteins(cores,pep_data,data_from_fasta,out_file_path)
#    else:
    distribute_peptides(cores,pep_data,data_from_fasta,out_file_path)

def get_arguments():
    parser= ArgumentParser()
    parser.add_argument("-peptides",required=True,help="full path to peptides file")
    parser.add_argument("-database", help="full path to database file")
    parser.add_argument("-out",required=True,help ="full path of output directory")
    parser.add_argument("-outfile",required=True,help ="output file name")
    parser.add_argument("-cores",type=int,help ="Number of cores to utilize..")
    #parser.add_argument("-exRoot",help="Root exclude folder")
    parser.add_argument("-exDirs",nargs="+",help="Name of directories to exclude...")
    parser.add_argument("-annotation",choices=["gencode","refseq","gnomon","custom","uniprot"],help="annotation generated by..")
    parser.add_argument("-dbtype",choices=["transcript","protein","gene"],required=True,help="Nucleic acid or amino acid")
    parser.add_argument("-frames",choices =[3,6], default=6,type=int,help="Frames to map peptides..")
    parser.add_argument("-gff",help="full path to gff3 file path")
    parser.add_argument("-accessions",help="full path to accessions file path")
    parser.add_argument("-genome",help="full path to genome file path")
    parser.add_argument("-header",type=bool,default=False,help="does peptides file has header")
    args = parser.parse_args()
    if args.dbtype in ("transcript","gene"):
        if not args.gff:
            print("for transcript mapping provide the gff file...with option -gff")
            sys.exit()
        if args.dbtype=="gene":
            if not args.genome:
                print("for gene mapping provide the genome file...with option -genome")
                sys.exit()

        if args.annotation in ("refseq","gnomon"):
            if not args.accessions:
                print("For refseq or gnomon mapping provide the accessions file...with option -accessions")
                sys.exit()

    return args

if __name__=="__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    print("total time for finding peptides..",str(round((end_time-start_time)/60.0,2)),"..mins..")
    print("Done....!")