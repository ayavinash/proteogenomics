import pandas as pd
import os,sys
import numpy as np
#from useful_functions import get_data_from_fasta

def main():
#    gnomon_prots = get_data_from_fasta("D:/Data/genome/ncbi/ann_release_108/Gnomon_prot.fsa") 

    data = pd.read_csv("D:/Data/Davide/BC/results/proteogenomics/filtered/peptides_modx/proteo_peptides.csv")
    exonskip = pd.read_csv("D:/Data/Davide/BC/results/proteogenomics/filtered/proteogenomics_mapping/exonskip/exonskipped_peptides_bed.csv")
    ref = pd.read_csv("D:/Data/Davide/BC/results/proteogenomics/filtered/proteogenomics_mapping/gnomon/gnomon_ref_bed_refseq.csv")
    alt = pd.read_csv("D:/Data/Davide/BC/results/proteogenomics/filtered/proteogenomics_mapping/gnomon/gnomon_alt_bed_refseq.csv")
    gencode = pd.read_csv("D:/Data/Davide/BC/results/proteogenomics/filtered/proteogenomics_mapping/gencode/gencode_transcripts_bed_pgp.csv")

#    gnomon_prots["is_partial"] = gnomon_prots.description.str.contains("internal stops and/or frameshifts")
#    full_lengths = gnomon_prots[~gnomon_prots.sequence.apply(str).str.contains("X")]
#    full_lengths = full_lengths[~full_lengths.description.str.contains("frameshifts")]

    exonskip["pep_locus"] = exonskip.chr.str.cat(exonskip.pep_start_on_chr.apply(str),sep=":").str.cat(exonskip.pep_end_on_chr.apply(str),sep="-")
    ref["pep_locus"] = ref.chr.str.cat(ref.pep_start_on_chr.apply(str),sep=":").str.cat(ref.pep_end_on_chr.apply(str),sep="-")
    alt["pep_locus"] = alt.chr.str.cat(alt.pep_start_on_chr.apply(str),sep=":").str.cat(alt.pep_end_on_chr.apply(str),sep="-")
    gencode["pep_locus"] = gencode.chr.str.cat(gencode.pep_start_on_chr.apply(str),sep=":").str.cat(gencode.pep_end_on_chr.apply(str),sep="-")
    

### #exonskip 
    ref =ref[(ref.tryptic_match=="Yes") & ~(ref.peptide.isin(exonskip.peptide))]
    alt =alt[(alt.tryptic_match=="Yes") & ~(alt.peptide.isin(exonskip.peptide))]
    gencode = gencode[(gencode.tryptic_match=="Yes") & ~(gencode.peptide.isin(exonskip.peptide))]

    

    data["pgp"] = np.nan
    data["gene"] = np.nan
    #data["pep_locus"] = np.nan
    data["pep_locus"] = data.pep_seq.map(exonskip.groupby("peptide")["pep_locus"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.pgp[data.pep_seq.isin(exonskip.peptide)]="exonskip_pep"
    data["gene"] = data.pep_seq.map(exonskip.groupby("peptide")["gene_name"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    #valid = gencode.copy(deep=True)

### Gnomon full lengths..
    ref["assembly"] = "ref"
    alt["assembly"] = "alt"
    
    gnomon = pd.concat([ref,alt],ignore_index=True)
    
    current = gnomon[gnomon.refseq_gene_name_biotype.notnull()]
    
    current =current[current.refseq_gene_name_biotype.str.contains("protein_coding")]
    current["PGP_type"]= "novel_isoform_PEP"
    current["assigned_in_ref"] = current.peptide.isin(current.peptide[current.assembly=="ref"])
    tmp  = current.groupby("assembly").apply(lambda x:x.to_dict())
    ref = pd.DataFrame(tmp["ref"])
    alt = pd.DataFrame(tmp["alt"])
    alt = alt[~alt.assigned_in_ref]
    current = pd.concat([ref,alt],ignore_index=True)
    
    data["refseq_gene_name"] =np.nan
    data["refseq_gene_biotype"] =np.nan
    data["assembly"] = np.nan

    

    data["current"]=False
    data.current[data.pgp.isnull()] =True
    data.assembly[data.current] = data.pep_seq.map(current.groupby("peptide")["assembly"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.pgp[data.current]= data.pep_seq.map(current.groupby("peptide")["PGP_type"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.gene[data.current]= data.pep_seq.map(current.groupby("peptide")["gene_name"].apply(lambda x: x.drop_duplicates().apply(str).str.cat(sep="#")))
    data.pep_locus[data.current]= data.pep_seq.map(current.groupby("peptide")["pep_locus"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.refseq_gene_name[data.current] = data.pep_seq.map(
        current.groupby("peptide")["refseq_gene_name_biotype"].apply(lambda x: x.drop_duplicates().str.cat(sep=";")).apply(
            lambda x: "#".join([y.split("#")[0] for y in x.split(";")])))

    data.refseq_gene_biotype[data.current] = data.pep_seq.map(
        current.groupby("peptide")["refseq_gene_name_biotype"].apply(lambda x: x.drop_duplicates().str.cat(sep=";")).apply(
            lambda x: "#".join([y.split("#")[1] for y in x.split(";")])))





## Gencode protein coding transripts
    data["current"] = False
    data.current[data.pgp.isnull()]=True
    current= gencode[gencode.has_cds=="Yes"]

    data.pgp[data.current]= data.pep_seq.map(current[current.has_cds=="Yes"].groupby("peptide")["PGP_type"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.gene[data.current]= data.pep_seq.map(current[current.has_cds=="Yes"].groupby("peptide")["gene_name"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.pep_locus[data.current]= data.pep_seq.map(current[current.has_cds=="Yes"].groupby("peptide")["pep_locus"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))


## Gencode non coding transcripts
    data["current"] = False
    data.current[data.pgp.isnull()]=True
    gencode.PGP_type[gencode.PGP_type.isnull()]="non_coding_PEP"
    data.pgp[data.current]= data.pep_seq.map(gencode[gencode.has_cds=="No"].groupby("peptide")["PGP_type"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.gene[data.current]= data.pep_seq.map(gencode[gencode.has_cds=="No"].groupby("peptide")["gene_name"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.pep_locus[data.current]= data.pep_seq.map(gencode[gencode.has_cds=="No"].groupby("peptide")["pep_locus"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))


## Gnomon novel cds
    ref = pd.read_csv("D:/Data/Davide/BC/results/proteogenomics/filtered/proteogenomics_mapping/gnomon/gnomon_ref_bed_refseq.csv")
    alt = pd.read_csv("D:/Data/Davide/BC/results/proteogenomics/filtered/proteogenomics_mapping/gnomon/gnomon_alt_bed_refseq.csv")

    ref["pep_locus"] = ref.chr.str.cat(ref.pep_start_on_chr.apply(str),sep=":").str.cat(ref.pep_end_on_chr.apply(str),sep="-")
    alt["pep_locus"] = alt.chr.str.cat(alt.pep_start_on_chr.apply(str),sep=":").str.cat(alt.pep_end_on_chr.apply(str),sep="-")
    ref =ref[(ref.tryptic_match=="Yes") & ~(ref.peptide.isin(exonskip.peptide))]
    alt =alt[(alt.tryptic_match=="Yes") & ~(alt.peptide.isin(exonskip.peptide))]
    ref["assembly"] = "ref"
    alt["assembly"] = "alt"
    alt = alt[~alt.peptide.isin(ref.peptide)]
    
    gnomon = pd.concat([ref,alt],ignore_index=True)
    gnomon.refseq_gene_name_biotype[gnomon.refseq_gene_name_biotype.isnull()]="##NA##"
    current = gnomon[~gnomon.refseq_gene_name_biotype.str.contains("protein_coding")]
    current["PGP_type"]= current.refseq_gene_name_biotype.apply(lambda x: "non_coding_PEP" if x!="##NA##" else "novel_CDS_PEP")

    data["current"]=False
    data.current[data.pgp.isnull()] =True
    data.assembly[data.current] = data.pep_seq.map(current.groupby("peptide")["assembly"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.pgp[data.current]= data.pep_seq.map(current.groupby("peptide")["PGP_type"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.gene[data.current]= data.pep_seq.map(current.groupby("peptide")["gene_name"].apply(lambda x: x.drop_duplicates().apply(str).str.cat(sep="#")))
    data.pep_locus[data.current]= data.pep_seq.map(current.groupby("peptide")["pep_locus"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))

    data.refseq_gene_name[data.current] = data.pep_seq.map(
        current.groupby("peptide")["refseq_gene_name_biotype"].apply(lambda x: x.drop_duplicates().str.cat(sep=";")).apply(
            lambda x: "#".join([y.split("#")[0] for y in x.split(";")])))

    data.refseq_gene_biotype[data.current] = data.pep_seq.map(
        current.groupby("peptide")["refseq_gene_name_biotype"].apply(lambda x: x.drop_duplicates().str.cat(sep=";")).apply(
            lambda x: "#".join([y.split("#")[1] for y in x.split(";")])))


    
    data.gene[data.refseq_gene_name.notnull()]=data.refseq_gene_name[data.refseq_gene_name.notnull()] 
    data["assigned_pgp"]= data.apply(lambda x: "ambigiuos_PEP" if ("#" in x["pgp"]) or ("#" in x["pep_locus"]) or ("#" in x["gene"]) else x["pgp"],axis=1)
    data.assigned_pgp = data.assigned_pgp.str.replace("uORF_PEP","uORF-pep")
    data.assigned_pgp = data.assigned_pgp.str.replace("dORF_PEP","dORF-pep")
    data.assigned_pgp = data.assigned_pgp.str.replace("altCDS_PEP","altCDS-pep")
    data.assigned_pgp = data.assigned_pgp.str.replace("novel_CDS_PEP","Novel-CDS-pep")
    data.assigned_pgp = data.assigned_pgp.str.replace("novel_isoform_PEP","Novel-isoform-pep")
    data.assigned_pgp = data.assigned_pgp.str.replace("ambigiuos_PEP","Ambiguous-pep")
    data.assigned_pgp = data.assigned_pgp.str.replace("exonskip_pep","Exonskip-pep")
    data.assigned_pgp = data.assigned_pgp.str.replace("non_coding_PEP","Non-coding-pep")
    
    data.assembly[data.assembly.isnull()]= "ref"

    categories=[							"Ambiguous-pep",
											"Novel-CDS-pep",
											"Non-coding-pep",
											"dORF-pep",
											"altCDS-pep",
											"uORF-pep",
											"Novel-isoform-pep",
											"Exonskip-pep"
											]
    data.assigned_pgp =data.assigned_pgp.astype("category",categories=categories,ordered=True) 
    data= data.sort_values("assigned_pgp",ascending=False)
    data.to_csv("D:/Data/Davide/BC/results/proteogenomics/filtered/proteogenomics_mapping/assigned/proteo_peptides_pgp.csv",index=False)



## Generate bed file....
    ref = pd.read_csv("D:/Data/Davide/BC/results/proteogenomics/filtered/proteogenomics_mapping/gnomon/gnomon_ref_bed_refseq.csv")
    alt = pd.read_csv("D:/Data/Davide/BC/results/proteogenomics/filtered/proteogenomics_mapping/gnomon/gnomon_alt_bed_refseq.csv")

    bed_data = ref[["chr","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRGB","total_exons","blocks_lengths","blocks_starts"]]
    bed_data.rename(columns={"chr":"#chr"},inplace=True)
    bed_data.to_csv("D:/Data/Davide/BC/results/proteogenomics/filtered/proteogenomics_mapping/gnomon/gnomon_ref_bed_refseq.bed",sep="\t",index=False)

    bed_data = alt[["chr","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRGB","total_exons","blocks_lengths","blocks_starts"]]
    bed_data.rename(columns={"chr":"#chr"},inplace=True)
    bed_data.to_csv("D:/Data/Davide/BC/results/proteogenomics/filtered/proteogenomics_mapping/gnomon/gnomon_alt_bed_refseq.bed",sep="\t",index=False)

    gencode = pd.read_csv("D:/Data/Davide/BC/results/proteogenomics/filtered/proteogenomics_mapping/gencode/gencode_transcripts_bed_pgp.csv")
    bed_data = gencode[["chr","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRGB","total_exons","blocks_lengths","blocks_starts"]]
    bed_data.rename(columns={"chr":"#chr"},inplace=True)

    bed_data.to_csv("D:/Data/Davide/BC/results/proteogenomics/filtered/proteogenomics_mapping/gencode/gencode_transcripts_bed_pgp.bed",sep="\t",index=False)


    exonskip = pd.read_csv("D:/Data/Davide/BC/results/proteogenomics/filtered/proteogenomics_mapping/exonskip/exonskipped_peptides_bed.csv")
    bed_data = exonskip[["chr","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRGB","total_exons","blocks_lengths","blocks_starts"]]
    bed_data.rename(columns={"chr":"#chr"},inplace=True)

    bed_data.to_csv("D:/Data/Davide/BC/results/proteogenomics/filtered/proteogenomics_mapping/exonskip/exonskip_bed_pgp.bed",sep="\t",index=False)



if __name__=="__main__":
    main()