####
#### Script to filter blast result for amino acid difference 
####
####


import pandas as pd
from argparse import ArgumentParser


def get_arguments():
    parser= ArgumentParser()
    parser.add_argument("-blast", help="full path to blast output")
    parser.add_argument("-out",required=True,help ="full path to output file")
    parser.add_argument("-diff",default=2,type=int,help ="amino acids difference")
    args= parser.parse_args()
    return args


def analyze_blast(blast_out_file_path,diff):
    """
    check for amino acid difference and print a mapped column

    """
    columns=["peptide","protein","per_identity","alignment_length","mismatch","gaps","pep_start","pep_end","prot_start","prot_end","eval","bit_score","qseq","sseq"]
    data = pd.read_csv(blast_out_file_path,names =columns,header=None)
    data["pep_len"] = data.peptide.apply(len)
    data["aa_in_gaps"] = data.qseq.str.count("-").add(data.sseq.str.count("-"))
    data["total_aa_diff"] = data.mismatch.add(data.aa_in_gaps)
    data["alignment_diff"] = data.pep_len.sub(data.alignment_length)
    data["alignment_diff_of_2"] = data.alignment_diff.between(-diff,diff,inclusive=True)
    data["mapped"] = data.apply(lambda x: (x["total_aa_diff"]<=diff) & (x["alignment_diff_of_2"]),axis=1)
    return data



def main():
    args = get_arguments()
    data = analyze_blast(args.blast,args.diff)
    data.to_csv(args.out,index=False)

if __name__=="__main__":
    main()
    print("Done!")
    