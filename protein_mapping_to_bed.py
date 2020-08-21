### Takes peptides hits in proteins and convertes to genomic coordinate##
### Author: Avinash Yadav
### PhD student, Scuola Normale Superiore
### contact: ayavinash@gmail.com

import pandas as pd
#import json
import numpy as np
import sys
import os,math
from multiprocessing import Process,cpu_count,Manager
from argparse import ArgumentParser
from useful_functions import map_refseq_chr_accessions,chunks,assign_cores,get_data_from_fasta,get_split_point


def flatten_csv_file(data):
    row_list = []
    for index in data.index:
        row = data.loc[index,]
        for index,pep_hit in enumerate(row["pep_hits"]):
            new_row = row.copy()
            new_row["pep_hits"] = pep_hit
            new_row["aa_res_before"] = row["aa_res_before"][index]
            new_row["aa_res_after"] = row["aa_res_after"][index]
            row_list.append(new_row)
    new_data = pd.concat(row_list, axis = 1)
    new_data = new_data.transpose()
    new_data.index = range(len(row_list))
    return new_data


def get_pep_start_on_transcript(data,gtf,map_type):
    if map_type=="transcript":

        data.set_index(["sequence_id","strand","chr","seq_type"],inplace=True)
        data["parent_id"] = gtf.drop_duplicates(["trans_id","strand","chr","seq_type"]).set_index(["trans_id","strand","chr","seq_type"])["parent_id"]
    elif map_type=="protein":
        data["seq_type"] ="original"
        data.set_index(["sequence_id","seq_type"],inplace=True)
        data["parent_id"] = gtf.drop_duplicates(["trans_id","seq_type"]).set_index(["trans_id","seq_type"])["parent_id"]
    else:
        print("This is not implemented...")
        sys.exit()

    data = data[data.parent_id.notnull()]
    data.reset_index(inplace=True)
    data.set_index(["parent_id","seq_type"],inplace=True)
    first_cds_gtf = gtf.drop_duplicates(["parent_id","seq_type"]).set_index(["parent_id","seq_type"])
  
    if "frame" not in data.columns:
        data["frame"] = first_cds_gtf.frame.apply(int)
    if "strand" not in data.columns:
        data["strand"] = first_cds_gtf.strand
    if "chr" not in data.columns:
        data["chr"] = first_cds_gtf.chr
    

    
    data["gene_name"] = first_cds_gtf.gene_name
    data.frame[data.frame.isnull()] = 0
    data["pep_start_on_trans"] = data.pep_hits*3 + data.frame + 1
    data["pep_end_on_trans"] = data["pep_start_on_trans"] + 3*data.peptide.str.len() -1
    #data["exon_cum_len"] = gtf.groupby(["parent_id","seq_type"])["exon_cum_len"]
    data.reset_index(inplace=True)
    return data

def shorten_gtf(data,map_type):
    if map_type=="protein":
        data = data[data.feature == "CDS"]
    elif map_type=="transcript":
        data = data[data.feature=="exon"]
    elif map_type=="gene":
        data = data[data.feature=="gene"]
    else:
        print("This is not implemmented....")
        sys.exit()
    return data



def read_gtf(gtf_file_path,primary_accessions,annotation,map_type,frames):
    print("Reading gtf file..")
    gff=pd.read_csv(gtf_file_path,header=None,comment="#",sep="\t",\
		names=["chr","source","feature","start","end","score","strand","frame","attributes"])
    data = gff.copy(deep=True)
    data = shorten_gtf(data,map_type)
    data["seq_type"] ="original"
    data["parent_id"] = data.attributes.str.extract('Parent=(.*?);')
    data["feature_id"] = data.attributes.str.extract('ID=(.*?);')
    if annotation in ("refseq","gnomon"):
        data = map_refseq_chr_accessions(data,primary_accessions)
  
    data.attributes = data.attributes+";"
    data = get_protein_data_mapping_type(data,annotation,map_type)
    data = make_cum_len_columns(data)
    if (frames==6) & (map_type!="protein"):
        data = generate_complement_annotations(data)
    
    data = data[data.trans_id.notnull()]
    return data

def generate_complement_annotations(data):
    comp_data = data.copy(deep=True)
    comp_data["seq_type"] ="complement"
    #comp_data.strand = comp_data.strand.map({"+":"-","-":"+"})
    comp_data.sort_values(["parent_id","seq_type"],ascending = [0,0],inplace=True)
    data = pd.concat([data,comp_data],ignore_index=True)
    return data    


def get_protein_data_mapping_type(data,annotation,map_type):
    if annotation =="refseq":
        data["gene_name"] = data.attributes.str.extract('gene=(.*?);')
        if map_type=="transcript":
            data["trans_id"] = data.attributes.str.extract(';transcript_id=(.*?);')
        elif map_type=="protein":
            data["trans_id"] = data.attributes.str.extract(';protein_id=(.*?);')
        else:
            print("#"*20, "this is unexpected..")

    elif annotation =="gencode":
        data = data.assign(gene_name = data.attributes.str.extract('gene_name=(.*?);'))
        data["trans_id"] = data.attributes.str.extract('transcript_id=(.*?);')
    elif annotation=="gnomon":
        data = data.assign(gene_name = data.attributes.str.extract(';gene=gene.(.*?);'))
        data = data.assign(trans_id = data.attributes.str.extract(';Name=GNOMON:(.*?);'))
    else:
        print("This annotation type is not implemented...Exiting..!")
        sys.exit()
    #data = data[data.trans_id.notnull()]
    return data

    #data.trans_id = data.trans_id.str.split(".").apply(lambda x:x[0])

def make_cum_len_columns_old(data):
    trans_exon_dict = data.groupby("trans_id")["trans_id"].count().to_dict()
    for trans in trans_exon_dict:
        count = trans_exon_dict[trans]
        last_val = count+2
        trans_exon_dict[trans] = list(range(1,last_val))
    data["exon_no"] = data.trans_id.apply(lambda x: trans_exon_dict[x].pop(0))
    data.exon_no = data.exon_no.apply(lambda x: int(x))
    data["exon_len"] = data.end - data.start + 1
    data["exon_cum_len"] = data.groupby("trans_id")["exon_len"].cumsum()
    
    #print("#"*20)
    #print(data)
    return data

def make_cum_len_columns(data):
    
    data.set_index(["parent_id","seq_type"],inplace=True)
    grouped = data.groupby(["parent_id","seq_type"])
    data["exon_no"] =1
    data["exon_no"] = grouped["exon_no"].cumsum()
    data["exon_len"] = data.end - data.start + 1
    data["exon_cum_len"] = grouped["exon_len"].cumsum()
    data.reset_index(inplace=True)
    return data




def get_exon_start_or_end_index(grouped,pep_start_or_end_on_trans,parent_id,seq_type):
    exon_cum_len = grouped.get_group((parent_id,seq_type))
    exon_start_or_end_index = exon_cum_len[pep_start_or_end_on_trans<=exon_cum_len].first_valid_index()
    return exon_start_or_end_index

def make_bed(data, gtf,return_dict):
    grouped = gtf.groupby(["parent_id","seq_type"])["exon_cum_len"]
    #print(data.columns)
    #print(gtf.columns)
    data = data[(data.parent_id.isin(set(gtf.parent_id)))]
    data["exon_start_index"] = data.apply(lambda x: get_exon_start_or_end_index(grouped,x["pep_start_on_trans"]
                            ,x["parent_id"],x["seq_type"]),axis=1)

    data["exon_end_index"] = data.apply(lambda x: get_exon_start_or_end_index(grouped,x["pep_end_on_trans"]
                            ,x["parent_id"],x["seq_type"]),axis=1)
    #print(data)
    data= data[(data.exon_start_index.notnull()) & (data.exon_end_index.notnull())]
    data.exon_start_index = data.exon_start_index.apply(int)
    data.exon_end_index = data.exon_end_index.apply(int)
    data = data[data.exon_start_index<=data.exon_end_index]

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
    #print(data.passing_exon_index)

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
    data["name"]= data.peptide.str.cat([data.gene_name,data.sequence_id],sep="#")
    data["itemRGB"]="255,0,0"
    data["score"]="1000"
    data["thickStart"]=data.chromStart
    data["thickEnd"]=data.chromEnd
    data.blocks_starts= data.blocks_starts.apply(lambda x: ",".join([str(y) for y in x]))
    data.blocks_lengths= data.blocks_lengths.apply(lambda x: ",".join([str(y) for y in x])) 
    bed_data = data[["chr","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRGB","total_exons","blocks_lengths","blocks_starts"]]
    bed_data.rename(columns={"chr":"#chr"},inplace=True)
    return_dict[str(os.getpid())] = data
    #print(bed_data)
    #return bed_data


def get_csv_data(csv_folder):
    
    csv_files = os.listdir(csv_folder)
    csv_files =[x for x in csv_files if x.endswith(".csv")]
    #print(csv_files)
    frames ={}
    for csv_file in csv_files:
        csv_file_path = os.path.join(csv_folder,csv_file)
        csv_data = pd.read_csv(csv_file_path)
        frames[csv_file] = csv_data
    if len(frames)>1:
        csv_data = pd.concat(frames.values(),ignore_index=True)
    return csv_data



def process_csv_maps(csv_data):
    csv_data.pep_hits = csv_data.pep_hits.apply(lambda x: [int(y) for y in x.split(";")])
    csv_data.aa_res_before = csv_data.aa_res_before.str.split(";")
    csv_data.aa_res_after = csv_data.aa_res_after.str.split(";")
    csv_data["pep_hits_len"] = csv_data.pep_hits.apply(len)
    csv_data1= csv_data[csv_data.pep_hits_len==1]
    csv_data2= csv_data[csv_data.pep_hits_len>1]
    print("flattening...")
    if len(csv_data2)>0:
        csv_data2=flatten_csv_file(csv_data2)
    csv_data1.pep_hits= csv_data1.pep_hits.apply(lambda x: x[0])
    csv_data1.aa_res_before= csv_data1.aa_res_before.apply(lambda x: x[0])
    csv_data1.aa_res_after= csv_data1.aa_res_after.apply(lambda x: x[0])
    csv_data = pd.concat([csv_data1,csv_data2],ignore_index=True)
    #print(csv_data)
    if csv_data.empty:
        sys.exit("No data...")
    return csv_data

def convert_to_protein_maps(annotation_source,dbtype,csv_data,database_file_path):
    split_point =get_split_point(annotation_source,dbtype)
    fasta_data = get_data_from_fasta(database_file_path)
    fasta_data["sequence_id"] = fasta_data.fasta_id.apply(lambda x: x.split("|")[split_point])
    del fasta_data["fasta_id"]
    fasta_data.sequence= fasta_data.sequence.apply(str)
    seq_dict = (fasta_data.set_index("sequence_id")["sequence"]).to_dict()
    csv_data["peptide"] = csv_data.sequence_id.map(seq_dict)
    csv_data["pep_seq"] = csv_data.peptide.str.replace("L","I")
    csv_data["pep_hits"] = 0
    return csv_data


def initialize(args):
    if os.path.isfile(args.maps) & args.maps.endswith(".csv"):
        csv_data = pd.read_csv(args.maps)
    else:
        csv_data = get_csv_data(args.maps)
    csv_data.pep_hits=csv_data.pep_hits.astype(str)    
    csv_data = process_csv_maps(csv_data)
    if args.mapProteins:
        if not args.database:
            print("protein database file is required....")
            sys.exit()
        csv_data = convert_to_protein_maps(args.annotation,args.dbtype,csv_data,args.database)
    print(csv_data)
      #bed_file_path,out_file_path = get_out_file_paths(args.out,args.maps)
    bed_file_path = args.out
    if os.path.exists(args.out):
        print("Output file exists...Exiting..")
        sys.exit()
    
    gtf = read_gtf(args.gff,args.accessions,args.annotation,args.dbtype,args.frames)
    csv_data = get_pep_start_on_transcript(csv_data,gtf,args.dbtype)
    #print(csv_data)
    #print(csv_data.columns)
    index_list = csv_data.index.tolist()
    cpus = assign_cores(args.cores)
    #print("-"*20,cpus)
    #print(index_list)
    #print(csv_data)
    index_lists = list(chunks(index_list, math.ceil(len(index_list)/cpus)))
    #print(index_lists)
    cpu_csv_data_list = []
    gtf_list = []
    manager = Manager()
    return_dict = manager.dict()    

    for index_list in index_lists:
        #print(len(index_list))
        cpu_csv_data = csv_data.copy(deep=True)
        cpu_csv_data = cpu_csv_data.loc[index_list,:]
        cpu_csv_data_list.append(cpu_csv_data)
        gtf_list.append(gtf.copy(deep=True))
    p_list = []
    for cpu in range(cpus):
        p = Process(target = make_bed,args =(cpu_csv_data_list[cpu],gtf_list[cpu],return_dict))
        p_list.append(p)
    for p in p_list:
        p.start()
    for p in p_list:
        p.join()
    #print(return_dict)
    if len(return_dict)>1:
        bed_data = pd.concat(return_dict.values(),ignore_index=True)   
    elif len(return_dict)==1:
        bed_data = return_dict.values()[0]
    else:
        print("No hits....")
        sys.exit()
    if args.mapProteins:
        bed_data.name = bed_data.name.apply(lambda x: "#".join(x.split("#")[1:]))
    bed_data.to_csv(bed_file_path,index=False)
    return bed_data    
    
#   generate_single_file(single_dir,bed_dir,out_dir)

def get_arguments():
    parser= ArgumentParser()
    parser.add_argument("-maps",required=True,help="full path to peptide to protein mapped folder or file..")
    parser.add_argument("-mapProteins",default=False,type=bool,help="to map proteins")
    parser.add_argument("-database",default=None,help="full path to protein fasta..")
    parser.add_argument("-gff",required=True,help="full path to gff/gtf file path")
    parser.add_argument("-accessions",help="full path to accessions file path")
    parser.add_argument("-out",required=True,help ="full path of output file")
    parser.add_argument("-cores",type=int,default=1,help ="Number of cores to utilize..")
    parser.add_argument("-annotation",default=None,required=True,choices=["gencode","refseq","gnomon"],help ="Annotation source..")
    parser.add_argument("-dbtype",default=None,required=True,choices=["protein","transcript","gene","genome"],help ="Map onto..")
    parser.add_argument("-frames",type=int,default=6,choices=[3,6],help ="Number of frames ..")
    

    args = parser.parse_args()
    #print(args)
    return args

def main():
    args = get_arguments()
    initialize(args)

if __name__=="__main__":
    main()

