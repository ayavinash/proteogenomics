

import pandas as pd
import os,sys,math,time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from multiprocessing import cpu_count

def generate_frames_data(data):
    frames = [0, 1, 2]
    frames_data = {}
    for frame in frames:
        frame_data = data.copy(deep=True)
        frame_data = frame_data.assign(frame = frame)
        frames_data[frame] =frame_data
    data = pd.concat(frames_data.values(), ignore_index = True)
    del frames_data
    return data


def complement_data(data):
    print("Complementing....")
    comp_data = data.copy(deep=True)
    data = data.assign(seq_type = "original")
    comp_data = comp_data.assign(sequence = comp_data.sequence.apply(lambda x: x.reverse_complement()))
    comp_data = comp_data.assign(seq_type = "complement")
    comp_data["strand"] =comp_data.strand.map({"+":"-","-":"+"})
    data = pd.concat([data, comp_data], ignore_index = True)
    del comp_data
    return data        


def fastaToDict(fasta_file):
    """
    Convert a fasta file to a python dictionary.
    """

    handle = open(fasta_file, "r")
    record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()
    return record_dict


def get_data_from_fasta(database_file_path):
    trans_dict = fastaToDict(database_file_path)
    key_list, value_list = [], []
    des_list = []
    n = 1
    for fasta_id, value in trans_dict.items():
        key_list.append(fasta_id)
        value_list.append(value.seq)
        des_list.append(value.description)
        #n+=1
        #if n>1000:
        #    break
    del trans_dict
    fasta_data_dict = {"fasta_id" : key_list, "sequence" : value_list,"description":des_list}
    fasta_data = pd.DataFrame(data = fasta_data_dict)
    print(fasta_data.sequence.head())
    return fasta_data


def get_refseq_accessions(primary_accessions):
    tmp = pd.read_csv(primary_accessions,sep="\t")
    tmp = tmp.set_index("RefSeq Accession.version")
    tmp["#Chromosome"] =tmp["#Chromosome"].apply(lambda x: "chr"+x)
    tmp["#Chromosome"][tmp["#Chromosome"]=="chrMT"]="chrM"
    refseq_chr_dict = tmp["#Chromosome"].to_dict()
    return refseq_chr_dict

def map_refseq_chr_accessions(data,primary_accessions):
    all_chrs = {k:k for k in set(data.chr)}
    refseq_chr_dict = get_refseq_accessions(primary_accessions)
    refseq_chr_dict ={**all_chrs,**refseq_chr_dict}
    data.chr = data.chr.map(refseq_chr_dict)
    return data

def get_split_point(annotation_source,db_type):
    
    if (annotation_source=="gencode") & (db_type=="protein"):
        split_point = 1
    elif (annotation_source=="refseq") & (db_type=="protein"): 
        split_point = 3
    elif (annotation_source=="gnomon") & (db_type=="protein"):
        split_point = 2
    elif (annotation_source=="gencode") & (db_type=="transcript"):
        split_point = 0
    elif (annotation_source=="refseq") & (db_type=="transcript"):
        split_point = 3
    elif (annotation_source=="gnomon") & (db_type=="transcript"):
        split_point = 2
    elif (annotation_source=="uniprot") & (db_type=="protein"):
        split_point= 1
    else:
        print(annotation_source,db_type,"is not implemented...Exiting...")
        sys.exit()
    return split_point

def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i+n]

def assign_cores(cores):
    if cores >=cpu_count():
        cores = cpu_count()-2
    return cores
