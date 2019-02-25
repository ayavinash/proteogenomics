### Takes peptides hits in proteins and convertes to genomic coordinate##
### Author: Avinash Yadav
### PhD student, Scuola Normale Superiore
### contact: ayavinash@gmail.com

import pandas as pd
import json
import numpy as np
import sys
import os,math,shutil
from multiprocessing import Process,cpu_count
from argparse import ArgumentParser


def flatten_csv_file(data):
    row_list = []
    for index in data.index:
        row = data.loc[index,]
        for pep_hit in row["pep_hits"]:
            new_row = row.copy()
            new_row["pep_hits"] = pep_hit
            row_list.append(new_row)
    new_data = pd.concat(row_list, axis = 1)
    new_data = new_data.transpose()
    new_data.index = range(len(row_list))
    new_data["protein_id"] = new_data.prot_id.str.split(".").apply(lambda x: x[0])
    new_data["trans_id"] = new_data.protein_id
    return new_data

def get_pep_start_on_transcript(data,gtf):
    first_cds_gtf = gtf.drop_duplicates("trans_id")
    first_cds_gtf = first_cds_gtf.set_index("trans_id")
    first_cds_gtf.frame =first_cds_gtf.frame.apply(int) 
    trans_frame_dict = first_cds_gtf.frame.to_dict()
    data["frame"] = data.protein_id.map(trans_frame_dict)
    data.frame[data.frame.isnull()] = 0
    data["pep_start_on_trans"] = data.pep_hits*3 + data.frame + 1
    data["pep_end_on_trans"] = data["pep_start_on_trans"] + 3*data.peptide.str.len() -1
    return data

def read_gtf(gtf_file_path,primary_accessions):
    print("Reading gtf file..")
    data=pd.read_csv(gtf_file_path,header=None,comment="#",sep="\t",\
		names=["chr","source","feature","start","end","score","strand","frame","attributes"])
    data = data[data.feature == "CDS"]
    tmp = pd.read_csv(primary_accessions,sep="\t")
    tmp = tmp.set_index("RefSeq Accession.version")
    tmp["#Chromosome"] =tmp["#Chromosome"].apply(lambda x: "chr"+x)
    tmp["#Chromosome"][tmp["#Chromosome"]=="chrMT"]="chrM"
    refseq_chr_dict = tmp["#Chromosome"].to_dict()

    data.chr = data.chr.map(refseq_chr_dict)
    data = data[data.chr.notnull()]
    data.index = range(len(data))
    data = data.assign(gene_name = data.attributes.str.extract(';gene=(.*?);'))
    data = data.assign(trans_id = data.attributes.str.extract(';Name=(.*?);'))
    data = data[data.trans_id.notnull()]
    data.trans_id = data.trans_id.str.split(".").apply(lambda x:x[0])
    trans_exon_dict = data.groupby("trans_id")["trans_id"].count().to_dict()
    for trans in trans_exon_dict:
        count = trans_exon_dict[trans]
        last_val = count+2
        trans_exon_dict[trans] = list(range(1,last_val))
    data["exon_no"] = data.trans_id.apply(lambda x: trans_exon_dict[x].pop(0))
    data.exon_no = data.exon_no.apply(lambda x: int(x))
    data["exon_len"] = data.end - data.start + 1
    data["exon_cum_len"] = data.groupby("trans_id")["exon_len"].cumsum()
    data["seq_type"] ="original"
    return data

def make_bed(data, gtf,bed_file_path,out_file_path):
    tmp_gtf = gtf.set_index("trans_id")
    prot_strand_dict = tmp_gtf.strand.to_dict()
    prot_chr_dict = tmp_gtf.chr.to_dict()
    prot_gene_dict = tmp_gtf.gene_name.to_dict()
    data = get_pep_start_on_transcript(data,gtf)
    data["strand"] = data.protein_id.map(prot_strand_dict)
    data["chr"] = data.protein_id.map(prot_chr_dict)
    data["gene_name"] = data.protein_id.map(prot_gene_dict)
    data["seq_type"] ="original"
    data = data[data.protein_id.isin(set(gtf.trans_id))]
    print("@@@@@@@@@@@@@")
    print(data)
    data["exon_start_index"] = data.apply(lambda x: x["exon_cum_len"][x["pep_start_on_trans"]<=x["exon_cum_len"]].first_valid_index(),axis=1)
    data["exon_end_index"] = data.apply(lambda x: x["exon_cum_len"][x["pep_end_on_trans"]<=x["exon_cum_len"]].first_valid_index(),axis=1)
    data= data[(data.exon_start_index.notnull()) & (data.exon_end_index.notnull())]
    data = data[data.exon_start_index<=data.exon_end_index]
    data.exon_start_index = data.exon_start_index.apply(int)
    data.exon_end_index = data.exon_end_index.apply(int)
    data["start_exon_no"] = data.exon_start_index.apply(lambda x: gtf.exon_no.loc[x])
    data["end_exon_no"] = data.exon_end_index.apply(lambda x: gtf.exon_no.loc[x])
    data["total_exons"] = data.end_exon_no - data.start_exon_no +1
    data["previous_exon_cum_len"] = data.apply(lambda x: 0 if x["start_exon_no"]==1 else gtf.loc[x["exon_start_index"]-1,"exon_cum_len"],axis=1)
    data["previous_exon_cum_len_end"] = data.apply(lambda x: 0 if x["end_exon_no"]==1 else gtf.loc[x["exon_end_index"]-1,"exon_cum_len"],axis=1)
    data["cds_start_of_start_exon"] = data.exon_start_index.apply(lambda x: gtf.loc[x,"start"])
    data["cds_end_of_start_exon"] = data.exon_start_index.apply(lambda x: gtf.loc[x,"end"])
    data["cds_start_of_end_exon"] = data.exon_end_index.apply(lambda x: gtf.loc[x,"start"])
    data["cds_end_of_end_exon"] = data.exon_end_index.apply(lambda x: gtf.loc[x,"end"])
    pos_data = data[data.strand=="+"]
    neg_data = data[data.strand=="-"]
    pos_data["pep_start_on_chr"] = pos_data.cds_start_of_start_exon+pos_data.pep_start_on_trans-pos_data.previous_exon_cum_len-1
    pos_data["pep_end_on_chr"] = pos_data.cds_start_of_end_exon+pos_data.pep_end_on_trans-pos_data.previous_exon_cum_len_end-1
    neg_data["pep_start_on_chr"] = neg_data.cds_end_of_start_exon-(neg_data.pep_start_on_trans-neg_data.previous_exon_cum_len)+1
    neg_data["pep_end_on_chr"] = neg_data.cds_end_of_end_exon-(neg_data.pep_end_on_trans-neg_data.previous_exon_cum_len_end)+1
    data = pd.concat([pos_data,neg_data])
    data.pep_start_on_chr = data.pep_start_on_chr.apply(int)
    data.pep_end_on_chr = data.pep_end_on_chr.apply(int)
    data["blocks_starts"] = 0
    data.blocks_starts= data.blocks_starts.apply(lambda x:[int(x)])
    data["blocks_lengths"] = data.peptide.apply(len)*3
    data.blocks_lengths = data.blocks_lengths.apply(lambda x: [x])
    data["passing_lengths"] = data.apply(lambda x: [],axis=1)
    data["passing_starts"] = data.apply(lambda x: [],axis=1)
    data["passing_exon_index"] = data.apply(lambda x: list(range(x["exon_start_index"],x["exon_end_index"]+1)),axis=1)
    print(data.passing_exon_index)

    if len(data[data.total_exons>1]):
            
        data.blocks_lengths[data.total_exons>1]=data[data.total_exons>1].apply(lambda x: [x["cds_end_of_start_exon"]-x["pep_start_on_chr"]+1, x["pep_end_on_chr"]-x["cds_start_of_end_exon"]+1] if x["strand"]=="+" else 
                                [x["cds_end_of_end_exon"]-x["pep_end_on_chr"]+1,x["pep_start_on_chr"]-x["cds_start_of_start_exon"]+1],axis=1)
        data.blocks_starts[data.total_exons>1] = data[data.total_exons>1].apply(lambda x: [0,x["cds_start_of_end_exon"]-x["pep_start_on_chr"]] if x["strand"]=="+" else
                                    [0,x["cds_start_of_start_exon"]-x["pep_end_on_chr"]],axis=1)
#    data.passing_lengths[data.total_exons>2]=data[data.total_exons>2].apply(lambda x: (
#        ((gtf.end.loc[x["passing_exon_index"][1]:x["passing_exon_index"][-2]]).sub(
#            (gtf.start.loc[x["passing_exon_index"][1]:x["passing_exon_index"][-2]])))+1).tolist(),axis=1)
    if len(data[data.total_exons>2]):
            
        data.passing_lengths[data.total_exons>2]=data[data.total_exons>2].apply(lambda x: (
            ((gtf.end.loc[x["passing_exon_index"][1]:x["passing_exon_index"][-2]]).sub(
                (gtf.start.loc[x["passing_exon_index"][1]:x["passing_exon_index"][-2]])))+1).tolist(),axis=1)


        data.passing_starts[data.total_exons>2] = data[data.total_exons>2].apply(lambda x: [gtf.start.loc[y]-x["pep_start_on_chr"] for y in x["passing_exon_index"][1:-1]] if x["strand"]=="+" else
                                    [gtf.start.loc[y]-x["pep_end_on_chr"] for y in x["passing_exon_index"][1:-1]],axis=1)
        data.blocks_lengths[data.total_exons>2] = data.apply(lambda x: [x["blocks_lengths"][0]] + x["passing_lengths"] + [x["blocks_lengths"][-1]],axis=1)
        data.blocks_starts[data.total_exons>2] = data.apply(lambda x:  [0]+ x["passing_starts"]+[x["blocks_starts"][-1]],axis=1)
    data["chromStart"]=data.apply(lambda x: x["pep_start_on_chr"]-1 if x["strand"]=="+" else x["pep_end_on_chr"]-1,axis=1)
    data["chromEnd"]= data.apply(lambda x: x["pep_end_on_chr"] if x["strand"]=="+" else x["pep_start_on_chr"],axis=1)
    data.chromStart = data.chromStart.apply(int)
    data.chromEnd = data.chromEnd.apply(int)
    data["name"]= data.peptide.str.cat([data.gene_name,data.trans_id],sep="#")
    data["itemRGB"]="255,0,0"
    data["score"]="1000"
    data["thickStart"]=data.chromStart
    data["thickEnd"]=data.chromEnd
    data.blocks_starts= data.blocks_starts.apply(lambda x: ",".join([str(y) for y in x]))
    data.blocks_lengths= data.blocks_lengths.apply(lambda x: ",".join([str(y) for y in x])) 
    bed_data = data[["chr","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRGB","total_exons","blocks_lengths","blocks_starts"]]
    bed_data.rename(columns={"chr":"#chr"},inplace=True)
    bed_file_path = bed_file_path.replace(".bed","_"+str(os.getpid())+".bed")
    out_file_path = out_file_path.replace(".csv","_"+str(os.getpid())+".csv")
    bed_data.to_csv(bed_file_path,sep="\t",index=False)
    del data["exon_cum_len"]
    print(bed_data)

    data.to_csv(out_file_path,index=False)

def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i+n]



def make_dirs(out_dir):
    if os.path.exists(out_dir):
        try:
            shutil.rmtree(out_dir)
        except:
            sys.exit("Could not clean output directory..",out_dir,"! Exiting....")
    if not os.path.exists(out_dir):
        try:
            os.makedirs(out_dir)
        except:
            sys.exit("Could not create output directory..",out_dir,"! Exiting....")
    return out_dir

def assign_cores(cores):
    if cores >=cpu_count():
        cores = cpu_count()-2
    return cores


def initialize(args):
 
    out_dir = make_dirs(args.out)
    bed_dir = make_dirs(os.path.join(out_dir,"bed"))
    csv_dir = make_dirs(os.path.join(out_dir,"csv"))
    csv_files = os.listdir(args.maps)
    base_filename = os.path.commonprefix(os.listdir(args.maps)).strip("_")
    bed_file_path = os.path.join(bed_dir,base_filename+".bed")
    out_file_path = os.path.join(csv_dir,base_filename+".csv")
  

    frames ={}
    for csv_file in csv_files:
        csv_file_path = os.path.join(args.maps,csv_file)
        csv_data = pd.read_csv(csv_file_path)
        frames[csv_file] = csv_data
    if len(frames)>1:
        csv_data = pd.concat(frames.values(),ignore_index=True)

    #csv_data = csv_data.iloc[:1000,]

    csv_data.pep_hits = csv_data.pep_hits.apply(lambda x: json.loads(x))
    csv_data["pep_hits_len"] = csv_data.pep_hits.apply(len)
    print(csv_data.pep_hits_len)
    csv_data1= csv_data[csv_data.pep_hits_len==1]
    csv_data1["protein_id"] = csv_data1.prot_id.str.split(".").apply(lambda x: x[0])
    csv_data1["trans_id"] = csv_data1.protein_id    
    csv_data2= csv_data[csv_data.pep_hits_len>1]
    print("flattening...")
    if len(csv_data2)>0:
        csv_data2=flatten_csv_file(csv_data2)
    csv_data1.pep_hits= csv_data1.pep_hits.apply(lambda x: x[0])
    csv_data = pd.concat([csv_data1,csv_data2],ignore_index=True)
    print(csv_data)
    #sys.exit()
    gtf = read_gtf(args.gff,args.accessions)
    grouped = gtf.groupby("trans_id")
    print(len(csv_data1),len(csv_data2))
    #sys.exit()
    print("grouping")
    cum_len_dict = {x[0]:x[1]["exon_cum_len"] for x in grouped}
    csv_data["exon_cum_len"] = csv_data.protein_id.map(cum_len_dict)
    index_list = csv_data.index.tolist()
    cpus = assign_cores(args.cores)
    index_lists = list(chunks(index_list, math.ceil(len(index_list)/cpus)))
    print(index_lists)
    cpu_csv_data_list = []
    gtf_list = []
    for index_list in index_lists:
        print(len(index_list))
        cpu_csv_data = csv_data.copy(deep=True)
        cpu_csv_data = cpu_csv_data.loc[index_list,:]
        cpu_csv_data_list.append(cpu_csv_data)
        gtf_list.append(gtf.copy(deep=True))
    p_list = []
    for cpu in range(cpus):
        p = Process(target = make_bed,args =(cpu_csv_data_list[cpu],gtf_list[cpu],bed_file_path,out_file_path))
        p_list.append(p)
    for p in p_list:
        p.start()
    #for p in p_list:
    #    p.join()
#    generate_single_file(single_dir,bed_dir,out_dir)

def get_arguments():
    parser= ArgumentParser()
    parser.add_argument("-maps",help="full path to peptide to protein mapped folder")
    parser.add_argument("-gff",help="full path to gff/gtf file path")
    parser.add_argument("-accessions",help="full path to accessions file path")
    parser.add_argument("-out",help ="full path of output directory")
    parser.add_argument("-cores",type=int,help ="Number of cores to utilize..")
    
    args = parser.parse_args()
    return args

def main():
    args = get_arguments()
    initialize(args)

if __name__=="__main__":
    main()

