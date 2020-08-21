import pandas as pd
from argparse import ArgumentParser
from useful_functions import map_refseq_chr_accessions
from pyteomics import parser



def load_bed_file(bed_file_path):
    gnomon = pd.read_csv(bed_file_path,sep="\t")
    name_split = gnomon.name.str.split("#",expand=True)
    name_split.rename(columns={0:"peptide",1:"gene",2:"protein"},inplace=True)
    gnomon = pd.concat([gnomon,name_split],axis=1)
    return gnomon

def get_feature_overlap(gnomon,gff):
    out = pd.concat(list(gnomon.apply(lambda x:row_function(gff,x["peptide"],x["chr"],
                                    x["strand"],x["chromStart"],x["chromEnd"]),axis=1)),ignore_index=True)
    #out["gene_biotype"] = out.attributes.str.extract('gene_biotype=(.*?);')
    out["gene_name_biotype"] = out.gene_name.str.cat(out.gene_biotype,sep="#")
    gnomon["refseq_gene_name_biotype"] = gnomon.peptide.map(out.groupby("peptide")["gene_name_biotype"].apply(
        lambda x: x.drop_duplicates().str.cat(sep=";")))
    return gnomon


def row_function(gff,peptide,chrom,strand,chromStart,chromEnd):
    tmp = gff[(gff.chr==chrom) & (gff.strand==strand) &
                                    (chromStart+1>=gff.start) &
                                    (chromEnd<gff.end)]    
    tmp["peptide"] = peptide
    return tmp
   
def assign_if_peptide_is_tryptic(data):
    data["before_pep_after"] = data.aa_res_before.str.cat([data.peptide,data.aa_res_after],sep="")
    data.before_pep_after = data.before_pep_after.str.strip("*").str.strip("-")
    data["tryptic_match"] = data.apply(lambda x: "Yes" if x["peptide"] 
                                in parser.cleave(x["before_pep_after"], parser.expasy_rules['trypsin'],2)
                                else "No",axis=1)                                                 
    return data



def read_refseq(gff_file_path):
    gff=pd.read_csv(gff_file_path,header=None,comment="#",sep="\t",\
		names=["chr","source","feature","start","end","score","strand","frame","attributes"])
    gff = gff[gff.feature=="gene"]
    gff.attributes = gff.attributes+";"
    gff["gene_name"] = gff.attributes.str.extract(';gene=(.*?);')
    gff["gene_biotype"] = gff.attributes.str.extract(';gene_biotype=(.*?);')

    return gff


def get_arguments():
    parser= ArgumentParser()
    parser.add_argument("-maps",required=True,help="full path to gnomon bed file")
    parser.add_argument("-gff", required=True,help="full path to refseq gff3 file")
    parser.add_argument("-accessions",required=True,help ="full path to refseq accessions file")
    parser.add_argument("-out",required=True,help ="output file path")
    args = parser.parse_args()
    return args

def assign_if_in_main(data):
    main_peptides = data.peptide[data.chr.str.startswith("chr")]
    data["assigned_in_main"] = data.peptide.isin(main_peptides)
    return data

#def assign_pep_locus(data):
#    data



def main():
    args = get_arguments()
    #gnomon = load_bed_file(args.bed)
    gnomon = pd.read_csv(args.maps)
    gff = read_refseq(args.gff)
    gff = map_refseq_chr_accessions(gff,args.accessions)
    gnomon = get_feature_overlap(gnomon,gff)
    gnomon = assign_if_peptide_is_tryptic(gnomon)
    gnomon = assign_if_in_main(gnomon)
    print(gnomon.refseq_gene_name_biotype)
    print(gnomon.refseq_gene_name_biotype.dtype)
    
    #gnomon.refseq_gene_name_biotype= gnomon.refseq_gene_name_biotype.apply(lambda x: pd.Series(x).astype().str.cat(sep=";"))
    
    gnomon.to_csv(args.out,index=False)


if __name__=="__main__":
    main()

