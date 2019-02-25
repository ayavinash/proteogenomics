### Take bed files and generate unique hits
### Can generate unique peptide to gene mapping


import pandas as pd
import sys
import os,math,shutil
from multiprocessing import Process
from argparse import ArgumentParser
from protein_mapping_to_bed import make_dirs


def gene_specific(data,out_dir,base_filename):

    split_data = data.name.str.split("#",expand=True)
    split_data.rename(columns ={0:"peptide",1:"gene",2:"protein"},inplace=True)
    data = pd.concat([data,split_data],axis=1)
    data["gene_hits"] = data.peptide.map(data.groupby("peptide")["gene"].nunique().to_dict())
    data["protein_hits"] = data.peptide.map(data.groupby("peptide")["protein"].nunique().to_dict())
 
    data.itemRGB[data.gene_hits>1]="0,0,0"
    data.itemRGB[data.protein_hits==1] = "0,255,0"
    data.name[data.gene_hits>=1] = data[data.gene_hits>=1].apply(lambda x: x["peptide"]+"#"+x["gene"],axis=1)
    data.name[data.protein_hits==1] = data[data.protein_hits==1].apply(lambda x: x["peptide"]+"#"+x["gene"]+"#"+x["protein"],axis=1)


    bed_data = data[["#chr","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRGB","total_exons","blocks_lengths","blocks_starts"]]
#   bed_data = bed_data.drop_duplicates(["#chr","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRGB","total_exons","blocks_lengths","blocks_starts"])
    bed_data = bed_data.drop_duplicates(["#chr","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRGB","total_exons","blocks_lengths","blocks_starts"])

    bed_data.to_csv(os.path.join(out_dir,base_filename+"_unique.csv"),index=False)
    bed_data.to_csv(os.path.join(out_dir,base_filename+"_unique.bed"),sep="\t",index=False)


def generate_single_file(bed_dir):
    bed_files = os.listdir(bed_dir)
    base_filename = os.path.commonprefix(bed_files).strip("_")
    print(bed_files)
    bed_data_list= []
    for bed_file in bed_files:
        each_bed_data = pd.read_csv(os.path.join(bed_dir,bed_file),sep="\t")
        bed_data_list.append(each_bed_data)
    if len(bed_data_list)>1:
        bed_data = pd.concat(bed_data_list,ignore_index=True)
    else:
        bed_data = bed_data_list[0]
    
    #bed_data=bed_data.rename(columns={"chr":"#chr"})
    return (bed_data,base_filename)


def get_arguments():
    parser= ArgumentParser()
    parser.add_argument("-bed",help="full path to bed folder")
    parser.add_argument("-out",help ="full path of output folder")
    args = parser.parse_args()
    return args



def main():
    args = get_arguments()
    make_dirs(args.out)
    data,base_filename = generate_single_file(args.bed)
    #out_filepath = os.path.join(args.out,base_filename+".bed")
    gene_specific(data,args.out,base_filename)
    


if __name__=="__main__":
    main()
