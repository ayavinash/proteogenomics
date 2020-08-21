import pandas as pd
import os,sys
import numpy as np
from pyteomics import parser


def assign_if_peptide_is_tryptic(data):
    data["before_pep_after"] = data.aa_res_before.str.cat([data.peptide,data.aa_res_after],sep="")
    data.before_pep_after = data.before_pep_after.str.strip("*").str.strip("-")
    data["tryptic_match"] = data.apply(lambda x: "Yes" if x["peptide"] 
                                in parser.cleave(x["before_pep_after"], parser.expasy_rules['trypsin'],2)
                                else "No",axis=1)                                                 
    return data

def main():
#    gnomon_prots = get_data_from_fasta("D:/Data/genome/ncbi/ann_release_108/Gnomon_prot.fsa") 
    
    data = pd.read_csv("D:/Data/CPTAC/BC/Mascot/results/all_identified_peptides.csv",names=["pep_seq"],header=None)
    crap = pd.read_csv("D:/Data/CPTAC/BC/Mascot/results/crap_peptides.csv",header=None,names=["peptide"])
    gencode_proteins = pd.read_csv("D:/Data/CPTAC/BC/Mascot/results/proteogenomics_mapping/gencode/gencode_proteins_bed.csv")
    refseq = pd.read_csv("D:/Data/CPTAC/BC/Mascot/results/proteogenomics_mapping/refseq/refseq_proteins_bed.csv")
    uniprot = pd.read_csv("D:/Data/CPTAC/BC/Mascot/results/proteogenomics_mapping/uniprot/uniprot_proteins.csv")
    uniprot_maps = pd.read_csv("D:/Data/CPTAC/BC/Mascot/results/proteogenomics_mapping/uniprot/gene_mappings.tab.txt",sep="\t")
    uniprot_maps = uniprot_maps.drop_duplicates("From").set_index("From")["To"].to_dict()

    data["pgp"] = ""
    data["gene"] = ""
    data["pep_locus"] = ""
    data["assembly"] = ""

    gencode_peptides = set(gencode_proteins.peptide.drop_duplicates())    
    refseq_peptides = set(refseq.peptide.drop_duplicates())    
    uniprot_peptides = set(uniprot.peptide.drop_duplicates())    
    crap_peptides = set(crap.peptide.drop_duplicates())    

    data["gencode"] = data.pep_seq.isin(gencode_peptides)
    data["refseq"] = data.pep_seq.isin(refseq_peptides)
    data["uniprot"] = data.pep_seq.isin(uniprot_peptides)
    data["crap"] = data.pep_seq.isin(crap_peptides)
    
    data["canonical"] =  data.apply(lambda x: x["gencode"] | x["refseq"] | x["uniprot"] | x["crap"],axis=1)
    data["non_canonical"] = ~data.canonical

###Assing contaminant
    data.pgp[data.pep_seq.isin(crap_peptides)] ="contaminant_pep"

### Assign gencode
    gencode_proteins = assign_if_peptide_is_tryptic(gencode_proteins)
    gencode_proteins = gencode_proteins[gencode_proteins.tryptic_match=="Yes"]
    gencode_proteins["pep_locus"] = gencode_proteins.chr.str.cat(gencode_proteins.pep_start_on_chr.apply(str),sep=":").str.cat(gencode_proteins.pep_end_on_chr.apply(str),sep="-")

    gencode_peptides = set(gencode_proteins.peptide.drop_duplicates())
    
    data.pgp[(data.pgp=="") & (data.pep_seq.isin(gencode_peptides))] ="gencode_pep"
    data.pep_locus[data.pgp=="gencode_pep"] = data.pep_seq[data.pgp=="gencode_pep"].map(gencode_proteins.groupby("peptide")["pep_locus"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.gene[data.pgp=="gencode_pep"] = data.pep_seq[data.pgp=="gencode_pep"].map(gencode_proteins.groupby("peptide")["gene_name"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))

### Assign refseq
    refseq = assign_if_peptide_is_tryptic(refseq)
    refseq = refseq[refseq.tryptic_match=="Yes"]
    refseq["pep_locus"] = refseq.chr.str.cat(refseq.pep_start_on_chr.apply(str),sep=":").str.cat(refseq.pep_end_on_chr.apply(str),sep="-")
    refseq_peptides = set(refseq.peptide.drop_duplicates())
    data.pgp[(data.pgp=="") & (data.pep_seq.isin(set(refseq.peptide)))] ="refseq_pep"

    data.pep_locus[data.pgp=="refseq_pep"] = data.pep_seq[data.pgp=="refseq_pep"].map(refseq.groupby("peptide")["pep_locus"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.gene[data.pgp=="refseq_pep"] = data.pep_seq[data.pgp=="refseq_pep"].map(refseq.groupby("peptide")["gene_name"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))

## Assign uniprot
    uniprot = assign_if_peptide_is_tryptic(uniprot)
    uniprot = uniprot[uniprot.tryptic_match=="Yes"]
    uniprot["gene_name"] = uniprot.sequence_id.map(uniprot_maps)
    uniprot["pep_locus"] = ""

    data.pgp[(data.pgp=="") & (data.pep_seq.isin(set(uniprot.peptide)))] ="uniprot_pep"
    data.gene[data.pgp=="uniprot_pep"] = data.pep_seq[data.pgp=="uniprot_pep"].map(uniprot.groupby("peptide")["gene_name"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.pep_locus[data.pgp=="uniprot_pep"] = data.pep_seq[data.pgp=="uniprot_pep"].map(uniprot.groupby("peptide")["pep_locus"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))


### Assign proteogenomics 


    exonskip = pd.read_csv("D:/Data/CPTAC/BC/Mascot/results/proteogenomics_mapping/exonskip/exonskip_proteins_bed.csv")
    ref_orig = pd.read_csv("D:/Data/CPTAC/BC/Mascot/results/proteogenomics_mapping/gnomon/gnomon_proteins_ref_bed_refseq.csv")
    alt_orig = pd.read_csv("D:/Data/CPTAC/BC/Mascot/results/proteogenomics_mapping/gnomon/gnomon_proteins_alt_bed_refseq.csv")
    gencode = pd.read_csv("D:/Data/CPTAC/BC/Mascot/results/proteogenomics_mapping/gencode_transcripts/gencode_transcripts_bed_pgp.csv")
    blast = pd.read_csv("D:/Data/CPTAC/BC/Mascot/results/proteogenomics_mapping/blast/all_non_canonical_peptides_blast_out_extended_analyzed.csv")
    blast_mapped_peptides = set(blast.peptide[blast.mapped])

    data["blast_mapped"] = data.pep_seq.isin(blast_mapped_peptides)



#    gnomon_prots["is_partial"] = gnomon_prots.description.str.contains("internal stops and/or frameshifts")
#    full_lengths = gnomon_prots[~gnomon_prots.sequence.apply(str).str.contains("X")]
#    full_lengths = full_lengths[~full_lengths.description.str.contains("frameshifts")]

###Assign exonskip
    exonskip= exonskip[~exonskip.peptide.isin(blast_mapped_peptides)]
    exonskip = assign_if_peptide_is_tryptic(exonskip)
    exonskip = exonskip[exonskip.tryptic_match=="Yes"]
    exonskip["pep_locus"] = exonskip.chr.str.cat(exonskip.pep_start_on_chr.apply(str),sep=":").str.cat(exonskip.pep_end_on_chr.apply(str),sep="-")

    data.pgp[(data.pgp=="") & (data.pep_seq.isin(set(exonskip.peptide)))]="exonskip_pep"
    data.pep_locus[data.pgp=="exonskip_pep"] = data.pep_seq[data.pgp=="exonskip_pep"].map(exonskip.groupby("peptide")["pep_locus"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.gene[data.pgp=="exonskip_pep"] = data.pep_seq.map(exonskip.groupby("peptide")["gene_name"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))


### Assign gnomon novel isoform from ref assembly
    ref = ref_orig.copy(deep=True)
    ref = ref[~ref.peptide.isin(blast_mapped_peptides)]
    main =ref[(ref.tryptic_match=="Yes") & (ref.chr.str.startswith("chr"))  & (ref.refseq_gene_name_biotype.notnull())]
    alts =ref[(ref.tryptic_match=="Yes") & ~(ref.assigned_in_main)  & (ref.refseq_gene_name_biotype.notnull())]
    ref = pd.concat([main,alts])
    ref["protein_coding_genes_biotypes"] = ref.refseq_gene_name_biotype.str.split(";")
    ref["protein_coding_genes_biotypes"] = ref.protein_coding_genes_biotypes.apply(lambda x: [y for y in x if "protein_coding" in y])
    ref.protein_coding_genes_biotypes= ref.protein_coding_genes_biotypes.apply(lambda x: np.nan if not x else x)
    ref = ref[ref.protein_coding_genes_biotypes.notnull()]
    ref["pep_locus"] = ref.chr.str.cat(ref.pep_start_on_chr.apply(str),sep=":").str.cat(ref.pep_end_on_chr.apply(str),sep="-")
    ref["refseq_gene"] = ref.protein_coding_genes_biotypes.apply(lambda x: "#".join([y.split("#")[0] for y in x]))

    data.pgp[(data.pgp=="") & (data.pep_seq.isin(set(ref.peptide)))]="novel_isoform_pep"
    data.pep_locus[data.pgp=="novel_isoform_pep"] = data.pep_seq[data.pgp=="novel_isoform_pep"].map(ref.groupby("peptide")["pep_locus"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.gene[data.pgp=="novel_isoform_pep"] = data.pep_seq[data.pgp=="novel_isoform_pep"].map(ref.groupby("peptide")["refseq_gene"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.assembly[data.pgp=="novel_isoform_pep"] = "ref"


### Assign gnomon novel isoform from alt assembly
    alt = alt_orig.copy(deep=True)
    alt = alt[~alt.peptide.isin(blast_mapped_peptides)]
    main =alt[(alt.tryptic_match=="Yes") & (alt.chr.str.startswith("chr"))  & (alt.refseq_gene_name_biotype.notnull())]
    alts =alt[(alt.tryptic_match=="Yes") & ~(alt.assigned_in_main)  & (alt.refseq_gene_name_biotype.notnull())]
    alt = pd.concat([main,alts])
    alt["protein_coding_genes_biotypes"] = alt.refseq_gene_name_biotype.str.split(";")
    alt["protein_coding_genes_biotypes"] = alt.protein_coding_genes_biotypes.apply(lambda x: [y for y in x if "protein_coding" in y])
    alt.protein_coding_genes_biotypes= alt.protein_coding_genes_biotypes.apply(lambda x: np.nan if not x else x)
    alt = alt[alt.protein_coding_genes_biotypes.notnull()]
    alt["pep_locus"] = alt.chr.str.cat(alt.pep_start_on_chr.apply(str),sep=":").str.cat(alt.pep_end_on_chr.apply(str),sep="-")
    alt["refseq_gene"] = alt.protein_coding_genes_biotypes.apply(lambda x: "#".join([y.split("#")[0] for y in x]))

    data.assembly[(data.pgp=="") & (data.pep_seq.isin(set(alt.peptide)))] = "alt"
    data.pgp[(data.pgp=="") & (data.pep_seq.isin(set(alt.peptide)))]="novel_isoform_pep"
    data.pep_locus[(data.pgp=="novel_isoform_pep") & (data.assembly=="alt")] = data.pep_seq[(data.pgp=="novel_isoform_pep") & (data.assembly=="alt")].map(alt.groupby("peptide")["pep_locus"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.gene[(data.pgp=="novel_isoform_pep") & (data.assembly=="alt")] = data.pep_seq[(data.pgp=="novel_isoform_pep") & (data.assembly=="alt")].map(alt.groupby("peptide")["refseq_gene"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))


## Gencode protein coding transripts
    gencode = gencode[~gencode.peptide.isin(blast_mapped_peptides)]    
    gencode = gencode[(gencode.tryptic_match=="Yes")]
    gencode["pep_locus"] = gencode.chr.str.cat(gencode.pep_start_on_chr.apply(str),sep=":").str.cat(gencode.pep_end_on_chr.apply(str),sep="-")

    data["current"] = False
    data.current[data.pgp==""]=True
    current= gencode[gencode.has_cds=="Yes"]

    data.pgp[data.current]= data.pep_seq.map(current[current.has_cds=="Yes"].groupby("peptide")["PGP_type"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.gene[data.current]= data.pep_seq.map(current[current.has_cds=="Yes"].groupby("peptide")["gene_name"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.pep_locus[data.current]= data.pep_seq.map(current[current.has_cds=="Yes"].groupby("peptide")["pep_locus"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.pgp[data.pgp.isnull()]=""


## Gencode non coding transcripts
    data["current"] = False
    data.current[data.pgp==""]=True
    gencode.PGP_type[gencode.PGP_type.isnull()]="non_coding_pep"
    data.pgp[data.current]= data.pep_seq.map(gencode[gencode.has_cds=="No"].groupby("peptide")["PGP_type"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.gene[data.current]= data.pep_seq.map(gencode[gencode.has_cds=="No"].groupby("peptide")["gene_name"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.pep_locus[data.current]= data.pep_seq.map(gencode[gencode.has_cds=="No"].groupby("peptide")["pep_locus"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.pgp[data.pgp.isnull()]=""


### Assign gnomon non coding from ref assembly
    ref = ref_orig.copy(deep=True)
    ref = ref[~ref.peptide.isin(blast_mapped_peptides)]

    main =ref[(ref.tryptic_match=="Yes") & (ref.chr.str.startswith("chr"))  & (ref.refseq_gene_name_biotype.notnull())]
    alts =ref[(ref.tryptic_match=="Yes") & ~(ref.assigned_in_main)  & (ref.refseq_gene_name_biotype.notnull())]
    ref = pd.concat([main,alts])

    ref["protein_coding_genes_biotypes"] = ref.refseq_gene_name_biotype.str.split(";")
    ref["protein_coding_genes_biotypes"] = ref.protein_coding_genes_biotypes.apply(lambda x: [y for y in x if "protein_coding" not in y])
    ref.protein_coding_genes_biotypes= ref.protein_coding_genes_biotypes.apply(lambda x: np.nan if not x else x)
    ref = ref[ref.protein_coding_genes_biotypes.notnull()]
    ref["pep_locus"] = ref.chr.str.cat(ref.pep_start_on_chr.apply(str),sep=":").str.cat(ref.pep_end_on_chr.apply(str),sep="-")
    ref["refseq_gene"] = ref.protein_coding_genes_biotypes.apply(lambda x: "#".join([y.split("#")[0] for y in x]))
    

    data.pgp[(data.pgp=="") & (data.pep_seq.isin(set(ref.peptide)))]="non_coding_pep"
    data.pep_locus[data.pgp=="non_coding_pep"] = data.pep_seq[data.pgp=="non_coding_pep"].map(ref.groupby("peptide")["pep_locus"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.gene[data.pgp=="non_coding_pep"] = data.pep_seq[data.pgp=="non_coding_pep"].map(ref.groupby("peptide")["refseq_gene"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.assembly[data.pgp=="non_coding_pep"] = "ref"
    data.pgp[data.pgp.isnull()]=""


### Assign gnomon non coding from alt assembly 
    alt = alt_orig.copy(deep=True)
    alt = alt[~alt.peptide.isin(blast_mapped_peptides)]

    main =alt[(alt.tryptic_match=="Yes") & (alt.chr.str.startswith("chr"))  & (alt.refseq_gene_name_biotype.notnull())]
    alts =alt[(alt.tryptic_match=="Yes") & ~(alt.assigned_in_main)  & (alt.refseq_gene_name_biotype.notnull())]
    alt = pd.concat([main,alts])

    alt["protein_coding_genes_biotypes"] = alt.refseq_gene_name_biotype.str.split(";")
    alt["protein_coding_genes_biotypes"] = alt.protein_coding_genes_biotypes.apply(lambda x: [y for y in x if "protein_coding" not in y])
    alt.protein_coding_genes_biotypes= alt.protein_coding_genes_biotypes.apply(lambda x: np.nan if not x else x)
    alt = alt[alt.protein_coding_genes_biotypes.notnull()]
    alt["pep_locus"] = alt.chr.str.cat(alt.pep_start_on_chr.apply(str),sep=":").str.cat(alt.pep_end_on_chr.apply(str),sep="-")
    alt["refseq_gene"] = alt.protein_coding_genes_biotypes.apply(lambda x: "#".join([y.split("#")[0] for y in x]))

    data.assembly[(data.pgp=="") & (data.pep_seq.isin(set(alt.peptide)))] = "alt"
    data.pgp[(data.pgp=="") & (data.pep_seq.isin(alt.peptide))]="non_coding_pep"
    data.pep_locus[(data.pgp=="non_coding_pep") & (data.assembly=="alt")] = data.pep_seq[(data.pgp=="non_coding_pep") & (data.assembly=="alt")].map(alt.groupby("peptide")["pep_locus"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.gene[(data.pgp=="non_coding_pep") & (data.assembly=="alt")] = data.pep_seq[(data.pgp=="non_coding_pep") & (data.assembly=="alt")].map(alt.groupby("peptide")["refseq_gene"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.pgp[data.pgp.isnull()]=""


## Gnomon novel cds
    ref = ref_orig.copy(deep=True)
    alt = alt_orig.copy(deep=True)

    ref["pep_locus"] = ref.chr.str.cat(ref.pep_start_on_chr.apply(str),sep=":").str.cat(ref.pep_end_on_chr.apply(str),sep="-")
    alt["pep_locus"] = alt.chr.str.cat(alt.pep_start_on_chr.apply(str),sep=":").str.cat(alt.pep_end_on_chr.apply(str),sep="-")
    ref =ref[(ref.tryptic_match=="Yes")]
    alt =alt[(alt.tryptic_match=="Yes")]
    ref["assembly"] = "ref"
    alt["assembly"] = "alt"
    alt = alt[~alt.peptide.isin(set(ref.peptide))]
    
    gnomon = pd.concat([ref,alt],ignore_index=True)
    gnomon = gnomon[gnomon.peptide.isin(blast_mapped_peptides)]
    gnomon.refseq_gene_name_biotype[gnomon.refseq_gene_name_biotype.isnull()]="##NA##"
    current = gnomon[~gnomon.refseq_gene_name_biotype.str.contains("protein_coding")]
    current["PGP_type"]= current.refseq_gene_name_biotype.apply(lambda x: "non_coding_pep" if x!="##NA##" else "novel_CDS_pep")

    data["current"]=False
    data.current[data.pgp==""] =True
    data.assembly[data.current] = data.pep_seq.map(current.groupby("peptide")["assembly"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.pgp[data.current]= data.pep_seq.map(current.groupby("peptide")["PGP_type"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))
    data.gene[data.current]= data.pep_seq.map(current.groupby("peptide")["gene_name"].apply(lambda x: x.drop_duplicates().apply(str).str.cat(sep="#")))
    data.pep_locus[data.current]= data.pep_seq.map(current.groupby("peptide")["pep_locus"].apply(lambda x: x.drop_duplicates().str.cat(sep="#")))

    data["refseq_gene_name"] = np.nan
    data["refseq_gene_biotype"] = np.nan

    data.refseq_gene_name[data.current] = data.pep_seq.map(
        current.groupby("peptide")["refseq_gene_name_biotype"].apply(lambda x: x.drop_duplicates().str.cat(sep=";")).apply(
            lambda x: "#".join([y.split("#")[0] for y in x.split(";")])))

    data.refseq_gene_biotype[data.current] = data.pep_seq.map(
        current.groupby("peptide")["refseq_gene_name_biotype"].apply(lambda x: x.drop_duplicates().str.cat(sep=";")).apply(
            lambda x: "#".join([y.split("#")[1] for y in x.split(";")])))


    
    data.gene[data.refseq_gene_name.notnull()]=data.refseq_gene_name[data.refseq_gene_name.notnull()]

    data.pgp[data.pgp.isnull()]="Undefined"
    data.gene[data.gene.isnull()]="Undefined"
    data.pep_locus[data.pep_locus.isnull()]="Undefined"

    data["assigned_pgp"]= data.apply(lambda x: "ambigiuos_pep" if ("#" in x["pgp"]) or ("#" in x["pep_locus"]) or ("#" in x["gene"]) else x["pgp"],axis=1)
    data.assigned_pgp = data.assigned_pgp.str.replace("uORF_PEP","uORF-pep")
    data.assigned_pgp = data.assigned_pgp.str.replace("dORF_PEP","dORF-pep")
    data.assigned_pgp = data.assigned_pgp.str.replace("altCDS_PEP","altCDS-pep")
    data.assigned_pgp = data.assigned_pgp.str.replace("novel_CDS_pep","Novel-CDS-pep")
    data.assigned_pgp = data.assigned_pgp.str.replace("novel_isoform_pep","Novel-isoform-pep")
    data.assigned_pgp = data.assigned_pgp.str.replace("ambigiuos_pep","Ambiguous-pep")
    data.assigned_pgp = data.assigned_pgp.str.replace("exonskip_pep","Exonskip-pep")
    data.assigned_pgp = data.assigned_pgp.str.replace("non_coding_pep","Non-coding-pep")
    data.assigned_pgp = data.assigned_pgp.str.replace("contaminant_pep","contaminant-pep")
    data.assigned_pgp = data.assigned_pgp.str.replace("uniprot_pep","reference-pep")
    data.assigned_pgp = data.assigned_pgp.str.replace("refseq_pep","reference-pep")
    data.assigned_pgp = data.assigned_pgp.str.replace("gencode_pep","reference-pep")
    
    #data.assembly[data.assembly.isnull()]= "ref"

    categories=[							"Undefined",
                                            "contaminant-pep",
                                            "reference-pep",
                                            "Ambiguous-pep",
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
    data.to_csv("D:/Data/CPTAC/BC/Mascot/results/proteo_peptides_pgp.csv",index=False)



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