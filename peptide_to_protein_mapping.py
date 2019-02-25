
### Takes peptides and finds all hits in all proteins in the database##
### Author: Avinash Yadav
### PhD student, Scuola Normale Superiore
### contact: ayavinash@gmail.com

from argparse import ArgumentParser
import pandas as pd
import os,sys,math,time,shutil
from multiprocessing import Process,cpu_count
from Bio import SeqIO

def fastaToDict(fasta_file):
    """
    Convert a fasta file to a python dictionary.
    """

    handle = open(fasta_file, "r")
    record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()
    return record_dict


def generate_data_from_fasta(database_file_path):
    trans_dict = fastaToDict(database_file_path)
    key_list, value_list = [], []
    des_list = []

    for fasta_id, value in trans_dict.items():
         key_list.append(fasta_id)
         value_list.append(str(value.seq))
         des_list.append(value.description)
    del trans_dict

    fasta_data_dict = {"fasta_id" : key_list, "sequence" : value_list,"description":des_list}
    fasta_data = pd.DataFrame(data = fasta_data_dict)
    fasta_data["prot_id"] = fasta_data.fasta_id.apply(lambda x: x.split("|")[3])
    fasta_data["sequence"] = fasta_data.sequence.apply(lambda x: x[1:] if x.startswith("X") else x)
    del fasta_data["fasta_id"]
    del fasta_data["description"]
    return fasta_data


def substring_indexes(substring, string):
    last_found = -1  
    while True:
        last_found = string.find(substring, last_found + 1)
        if last_found == -1:  
            break  
        yield last_found

def find_peptide_hits(data, pep_data, out_file_path):
    filename, ext = os.path.splitext(out_file_path) 
    out_file_path = filename+"_"+str(os.getpid())+ext
    #print(filename)
    #print(out_file_path)

    out = pd.concat(list(pep_data.apply(lambda x:row_function(data,x["pep_seq"],x["peptide"]),axis=1)),ignore_index=True)
    out["pep_hits"] = out.apply(lambda x: list(substring_indexes(x["pep_seq"],x["sequence"])),axis=1)
    del out["sequence"]
    if not out.empty:
        out.to_csv(out_file_path,index=False)
       
    else:
        print("Could not find any peptide hits....")


def row_function(data,pep,original_pep):
    pep_data = data[data.sequence.str.contains(pep,regex=False)]
    pep_data["pep_seq"] = pep
    pep_data["peptide"] = original_pep
    return pep_data


def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i+n]



def check_arguments(args):
    if not ((os.path.isfile(args.peptides)) and (os.path.splitext(args.peptides)[1]==".csv")):
        sys.exit("Input peptides is not a file or a .csv file...!Exiting...")

    if not ((os.path.isfile(args.database)) and (os.path.splitext(args.database)[1] in (".fa",".fasta"))):
        sys.exit("Input database is not fasta ...! Exiting...")

    if args.empty>=1:
        if os.path.exists(args.out):
            shutil.rmtree(args.out)
        
    else:
        if os.path.exists(args.out):
            sys.exit("output directory already exists...! Exiting...")
    if not os.path.exists(args.out):
        os.makedirs(args.out)    

    pep_file_name = os.path.splitext(os.path.basename(os.path.basename(args.peptides)))[0]
    database_file_name = os.path.splitext(os.path.basename(os.path.basename(args.database)))[0]
    out_file_name = pep_file_name+"_"+database_file_name+".csv"
    out_file_path = os.path.join(args.out,out_file_name)
    return out_file_path


def get_peptides(pep_file_path):
    pep_data = pd.read_csv(pep_file_path,header=None,names=["pep_seq"])
    #print(pep_data.pep_seq)
    pep_data = pep_data.drop_duplicates("pep_seq")
    pep_data["peptide"] = pep_data.pep_seq
    pep_data.pep_seq = pep_data.pep_seq.str.replace("L","I")
    #pep_data = pep_data.iloc[:1000,:]
    #print(pep_data)
    return pep_data

    #pep_index_list = pep_data.iloc[:100,:].index

def assign_cores(cores):
    if cores >=cpu_count():
        cores = cpu_count()-2
    return cores
  
def distribute(cores,pep_data,data_from_fasta,out_file_path):
    print("Finding peptide hits...")
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
        p = Process(target = find_peptide_hits,args =(cpu_data_frame_list[cpu],pep_data.loc[pep_index_lists[cpu],:],out_file_path))
        p_list.append(p)
    for p in p_list:
        p.start()


def get_arguments():
    parser= ArgumentParser()
    parser.add_argument("-peptides",help="full path to peptides csv file")
    parser.add_argument("-database", help="full path to protein database file")
    parser.add_argument("-out",help ="full path of output directory")
    parser.add_argument("-cores",type=int,help ="Number of cores to utilize..")
    parser.add_argument("-empty",type=int,default=0,help="empty output dir if exists")
    args = parser.parse_args()
    return args

def main():
    args = get_arguments()
    out_file_path = check_arguments(args)
    #os.makedirs(args.out)
    pep_data = get_peptides(args.peptides)
    #print(pep_data)
    cores = assign_cores(args.cores)
    data_from_fasta = generate_data_from_fasta(args.database)
    data_from_fasta["sequence"] = data_from_fasta.sequence.str.replace("L","I")
    distribute(cores,pep_data,data_from_fasta,out_file_path)


if __name__=="__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    print("total time for finding peptides..",str(end_time-start_time/60.0),"..mins..")

