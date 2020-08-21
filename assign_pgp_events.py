import pandas as pd
from argparse import ArgumentParser
import numpy as np
import sys,os
from pyteomics import parser
#from useful_functions import map_refseq_chr_accessions



def assign_novel_isoform(data):
    valid_data = data[(data.gene_biotype=="protein_coding") & (data.PGP_type!="assigned_in_ref_assembly")]
    data.PGP_type[valid_data.index]="novel_isoform_PEP"
    print(valid_data)
    valid_data = data[(data.PGP_type!="assigned_in_ref_assembly") & (data.gene_biotype.isin(set(["lncRNA","pseudogene"])))]
    data.PGP_type[valid_data.index]=valid_data.gene_biotype.apply(lambda x: x+"_PEP")
    valid_data = data[(data.mapped_gene=="NA") & (data.PGP_type!="assigned_in_ref_assembly")]
    data.PGP_type[valid_data.index] = "novel_CDS"

    data["assigned_locus"] = data.apply(lambda x: x["chr"]+":"+str(int(x["pep_start_on_chr"]))+"-"+str(int(x["pep_end_on_chr"])),axis=1) 
    return data


def assign_novel_CDS_from_gencode(data,gencode_gtf_file_path,chr_map_file):
    gencode=pd.read_csv(gencode_gtf_file_path,header=None,comment="#",sep="\t",
                        names=["chr","source","feature","start","end","score","strand","frame","attributes"])
#    gencode=pd.read_csv(gencode_gtf_file_path,header=None,comment="#",sep="\t",
#                        names=["chr","source","feature","start","end","score","strand","frame","attributes"])
    gencode_gene = gencode[gencode.feature=="gene"]
    data["mapped_by"] = "NA"
    refseq_mapped = data[data.mapped_gene!="NA"]
    data.mapped_by[refseq_mapped.index] ="RefSeq_annotation"
    chr_data = pd.read_csv(chr_map_file,sep="\t")
    chr_dict = chr_data.set_index("RefSeq Accession.version")["#Chromosome"].to_dict()
    print(chr_dict)

    valid_data = data[(data.assembly=="ref") & (data.mapped_gene=="NA")]
    valid_data["chr"] = valid_data.chr.map(chr_dict)

    gencode_gene["gene_name_commmon"] = gencode_gene.attributes.str.extract('gene_name \"(.*?)\";')
    gencode_gene["gene_biotype"] = gencode_gene.attributes.str.extract('gene_type \"(.*?)\";')
    gencode_gene["gene_name"] = gencode_gene.attributes.str.extract('gene_id \"(.*?)\";')
    
    gencode_gene.gene_biotype[gencode_gene.gene_biotype.isnull()] = gencode_gene.attributes[gencode_gene.gene_biotype.isnull()].str.extract('gene_biotype=(.*)')

    biotype_dict_1 = gencode_gene.set_index("gene_name").gene_biotype.to_dict()
    gene_name_dict_1 = gencode_gene.set_index("gene_name").gene_name_commmon.to_dict()
    print("@@@@@@@@@@@")
    print("applying methods...")

    valid_data.mapped_gene = valid_data.apply(lambda x: "#".join(str(x) for x in set(gencode_gene[(gencode_gene.chr==x["chr"]) &
                                    (gencode_gene.strand==x["strand"]) &
                                    (gencode_gene.start<x["pep_start_on_chr"]) &
                                    (gencode_gene.end>x["pep_end_on_chr"])].gene_name) 
                                if type(x)!=float),axis=1)
    valid_data.mapped_gene[valid_data.mapped_gene==""]="NA"
    data.mapped_gene[valid_data.index] = valid_data.mapped_gene
    data.mapped_by[valid_data.index] = "GENCODE_annotation"



    valid_data = valid_data[valid_data.mapped_gene=="NA"]
    valid_data.mapped_gene = valid_data.apply(lambda x: "#".join(str(x) for x in set(gencode_gene[(gencode_gene.chr==x["chr"]) &
                                    (gencode_gene.strand==x["strand"]) &
                                    (gencode_gene.start<x["pep_start_on_chr"]) &
                                    (gencode_gene.end>x["pep_start_on_chr"])].gene_name) 
                                if not type(x)==float),axis=1)
    valid_data.mapped_gene[valid_data.mapped_gene==""]="NA"
    data.mapped_gene[valid_data.index] = valid_data.mapped_gene
    data.mapped_by[valid_data.index] = "GENCODE_annotation"

    valid_data = valid_data[valid_data.mapped_gene=="NA"]
    valid_data.mapped_gene = valid_data.apply(lambda x: "#".join(str(x) for x in set(gencode_gene[(gencode_gene.chr==x["chr"]) &
                                    (gencode_gene.strand==x["strand"]) &
                                    (gencode_gene.start<x["pep_end_on_chr"]) &
                                    (gencode_gene.end>x["pep_end_on_chr"])].gene_name) 
                                if not type(x)==float),axis=1)
    valid_data.mapped_gene[valid_data.mapped_gene==""]="NA"
    data.mapped_gene[valid_data.index] = valid_data.mapped_gene
    data.mapped_by[valid_data.index] = "GENCODE_annotation"

    data.gene_biotype[data.mapped_by=="GENCODE_annotation"] = data.mapped_gene[data.mapped_by=="GENCODE_annotation"].map(biotype_dict_1)
    data.gene_name_common[data.mapped_by=="GENCODE_annotation"] = data.gene_name[data.mapped_by=="GENCODE_annotation"].map(gene_name_dict_1)
    return data

def assign_novel_isoform_gencode(data):
    valid_data = data[(data.mapped_by=="GENCODE_annotation")]

#    print(valid_data[valid_data.gene_biotype.isnull()])


    valid_data.PGP_type = valid_data.gene_biotype.apply(lambda x: x+"_PEP" if type(x)!=float else "NA")

    valid_data.PGP_type[valid_data.PGP_type.str.contains("protein_coding")] = "novel_isoform_PEP"
    valid_data.PGP_type[(valid_data.mapped_gene=="NA")]= "novel_CDS"
    data.PGP_type[valid_data.index] = valid_data.PGP_type
    return data

def assign_refseq_gene_and_pgp(data,ref_gtf_file_path,alt_gtf_file_path,gencode_gtf_file_path):
    print("reading refseq....")
    refseq=pd.read_csv(ref_gtf_file_path,header=None,comment="#",sep="\t",
                        names=["chr","source","feature","start","end","score","strand","frame","attributes"])

    alt_refseq=pd.read_csv(alt_gtf_file_path,header=None,comment="#",sep="\t",
                        names=["chr","source","feature","start","end","score","strand","frame","attributes"])
#    gencode=pd.read_csv(gencode_gtf_file_path,header=None,comment="#",sep="\t",
#                        names=["chr","source","feature","start","end","score","strand","frame","attributes"])
    refseq.attributes = refseq.attributes+";"
    refseq_gene = refseq[refseq.feature=="gene"]
    alt_refseq_gene = alt_refseq[alt_refseq.feature=="gene"]
    refseq_gene["gene_name"] = refseq_gene.attributes.str.extract('GeneID:(.*?),')
    refseq_gene["gene_name"][refseq_gene.gene_name.isnull()] = refseq_gene[refseq_gene.gene_name.isnull()].attributes.str.extract('GeneID:(.*?);')
    refseq_gene["gene_name_common"] = refseq_gene.attributes.str.extract('gene=(.*?);')
    refseq_gene["gene_biotype"] = refseq_gene.attributes.str.extract('gene_biotype=(.*?);')
    refseq_gene.gene_biotype[refseq_gene.gene_biotype.isnull()] = refseq_gene.attributes[refseq_gene.gene_biotype.isnull()].str.extract('gene_biotype=(.*)')

    alt_refseq_gene["gene_name"] = alt_refseq_gene.attributes.str.extract('GeneID:(.*?),')
    alt_refseq_gene["gene_name"][alt_refseq_gene.gene_name.isnull()] = alt_refseq_gene[alt_refseq_gene.gene_name.isnull()].attributes.str.extract('GeneID:(.*?);')
    alt_refseq_gene["gene_name_common"] = alt_refseq_gene.attributes.str.extract('gene=(.*?);')
    alt_refseq_gene["gene_biotype"] = alt_refseq_gene.attributes.str.extract('gene_biotype=(.*?);')
    alt_refseq_gene.gene_biotype[alt_refseq_gene.gene_biotype.isnull()] = alt_refseq_gene.attributes[alt_refseq_gene.gene_biotype.isnull()].str.extract('gene_biotype=(.*)')

    biotype_dict_1 = refseq_gene.set_index("gene_name").gene_biotype.to_dict()
    biotype_dict_2 =alt_refseq_gene.set_index("gene_name").gene_biotype.to_dict()
    gene_name_dict_1 = refseq_gene.set_index("gene_name").gene_name_common.to_dict()
    gene_name_dict_2 =alt_refseq_gene.set_index("gene_name").gene_name_common.to_dict()


    print("@@@@@@@@@@@")
    print("applying methods...")
    data["mapped_gene"]="NA"
    data.mapped_gene[data.assembly=="ref"] = data[data.assembly=="ref"].apply(lambda x: "#".join(str(x) for x in set(refseq_gene[(refseq_gene.chr==x["chr"]) &
                                    (refseq_gene.strand==x["strand"]) &
                                    (refseq_gene.start<x["pep_start_on_chr"]) &
                                    (refseq_gene.end>x["pep_end_on_chr"])].gene_name) 
                                if type(x)!=float),axis=1)
    data.mapped_gene[data.mapped_gene==""]="NA"

    data.mapped_gene[(data.assembly=="ref") & (data.mapped_gene=="NA")] = data[(data.assembly=="ref")& (data.mapped_gene=="NA")].apply(lambda x: "#".join(str(x) for x in set(refseq_gene[(refseq_gene.chr==x["chr"]) &
                                    (refseq_gene.strand==x["strand"]) &
                                    (refseq_gene.start<x["pep_start_on_chr"]) &
                                    (refseq_gene.end>x["pep_start_on_chr"])].gene_name) 
                                if not type(x)==float),axis=1)
    data.mapped_gene[data.mapped_gene==""]="NA"

    data.mapped_gene[(data.assembly=="ref") & (data.mapped_gene=="NA")] = data[(data.assembly=="ref")& (data.mapped_gene=="NA")].apply(lambda x: "#".join(str(x) for x in set(refseq_gene[(refseq_gene.chr==x["chr"]) &
                                    (refseq_gene.strand==x["strand"]) &
                                    (refseq_gene.start<x["pep_end_on_chr"]) &
                                    (refseq_gene.end>x["pep_end_on_chr"])].gene_name) 
                                if not type(x)==float),axis=1)
    data.mapped_gene[data.mapped_gene==""]="NA"




    data.mapped_gene[data.assembly=="alt"] = data[data.assembly=="alt"].apply(lambda x: "#".join(str(x) for x in set(alt_refseq_gene[(alt_refseq_gene.chr==x["chr"]) &
                                    (alt_refseq_gene.strand==x["strand"]) &
                                    (alt_refseq_gene.start<x["pep_start_on_chr"]) &
                                    (alt_refseq_gene.end>x["pep_end_on_chr"])].gene_name) 
                                if not type(x)==float),axis=1)

    data.mapped_gene[data.mapped_gene==""]="NA"

    data.mapped_gene[(data.assembly=="alt") & (data.mapped_gene=="NA")] = data[(data.assembly=="alt") & (data.mapped_gene=="NA")].apply(lambda x: "#".join(str(x) for x in set(alt_refseq_gene[(alt_refseq_gene.chr==x["chr"]) &
                                    (alt_refseq_gene.strand==x["strand"]) &
                                    (alt_refseq_gene.start<x["pep_start_on_chr"]) &
                                    (alt_refseq_gene.end>x["pep_start_on_chr"])].gene_name) 
                                if not type(x)==float),axis=1)
   
    data.mapped_gene[data.mapped_gene==""]="NA"
    data.mapped_gene[(data.assembly=="alt") & (data.mapped_gene=="NA")] = data[(data.assembly=="alt") & (data.mapped_gene=="NA")].apply(lambda x: "#".join(str(x) for x in set(alt_refseq_gene[(alt_refseq_gene.chr==x["chr"]) &
                                    (alt_refseq_gene.strand==x["strand"]) &
                                    (alt_refseq_gene.start<x["pep_end_on_chr"]) &
                                    (alt_refseq_gene.end>x["pep_end_on_chr"])].gene_name) 
                                if not type(x)==float),axis=1)
    
    data.mapped_gene[data.mapped_gene==""]="NA"
    data["gene_biotype"] = np.nan
    data["gene_biotype"][data.assembly=="ref"]= data.mapped_gene[data.assembly=="ref"].apply(lambda x: "#".join([biotype_dict_1[y] for y in x.split("#")]) if x !="NA"
                                               else np.nan)
    data["gene_biotype"][data.assembly=="alt"]= data.mapped_gene[data.assembly=="alt"].apply(lambda x: "#".join([biotype_dict_2[y] for y in x.split("#")]) if x !="NA"
                                                else np.nan)

    data["gene_name_common"] = np.nan
    data.gene_name_common[data.assembly=="ref"]=data.mapped_gene[data.assembly=="ref"].apply(lambda x: "#".join([gene_name_dict_1[y] for y in x.split("#")]) if x !="NA" 
                                                else np.nan)
    data.gene_name_common[data.assembly=="alt"]=data.mapped_gene[data.assembly=="alt"].apply(lambda x: "#".join([gene_name_dict_2[y] for y in x.split("#")]) if x !="NA" 
                                                else np.nan)


#    data["gene_biotype"] ="NA"
#    data.gene_biotype[(data.mapped_gene!="NA") & (~data.mapped_gene.str.contains("#"))]= data.mapped_gene[(data.mapped_gene!="NA") & (~data.mapped_gene.str.contains("#"))].map(biotype_dict_1)
#    data.gene_biotype[data.gene_biotype==""]= "NA"
#    data.gene_biotype[(data.mapped_gene!="NA") & (~data.mapped_gene.str.contains("#"))]= data.mapped_gene[(data.mapped_gene!="NA") & (~data.mapped_gene.str.contains("#"))].map(biotype_dict_2)



#    data["gene_name_common"] = "NA"
#    data.gene_name_common[(data.mapped_gene!="NA") & (data.assembly=="ref")]=data.mapped_gene[(data.mapped_gene!="NA") & (data.assembly=="ref")].map(gene_name_dict_1)
#    data.gene_name_common[(data.mapped_gene!="NA") & (data.assembly=="alt")]=data.mapped_gene[(data.mapped_gene!="NA") & (data.assembly=="alt")].map(gene_name_dict_2)

    print(data)
    return data



def assign_if_in_ref(data):
    data["PGP_type"] ="NA"
    data.PGP_type[data.tryptic_match=="No"]= "semi-tryptic"

    original = data[data.assembly=="ref"]
    complement = data[data.assembly=="alt"]
    print(set(original.peptide))
    print(set(complement.peptide))
    excl_complement = complement[~(complement.peptide.isin(set(original.peptide)))]
    original_assigned = complement[(complement.peptide.isin(set(original.peptide)))]
#    complement["PGP_type"] = "novel_CDS"
    #print(complement)

    data.PGP_type[excl_complement.index]="alt_assembly"
    data.PGP_type[original_assigned.index] ="assigned_in_ref_assembly"
    #print(data[data.peptide=="APTAGATASTRPGSEDRAPPSPRPAAAADSPWAAETAR"])
    #sys.exit()
    
    #print data.PGP_type 
    return data

def assign_if_in_main_chroms(data,ref_chrom_file,assembly):
    ref_chr_data = pd.read_csv(ref_chrom_file,sep="\t")
    
    ref_chr_dict = ref_chr_data.set_index("RefSeq Accession.version")["#Chromosome"].to_dict()
    #data["chrom_type"] = "NA"
    data.chrom_type[data.assembly==assembly] = data.chr[data.assembly==assembly].apply(lambda x: "Main"
                                  if x in ref_chr_dict else "Alts")
    main = data[(data.assembly==assembly) & (data.chrom_type=="Main")]
    alts = data[(data.assembly==assembly) & (data.chrom_type=="Alts")]

    #excl_alts = alts[~(alts.peptide.isin(set(main.peptide)))]
    assigned_in_main = alts[(alts.peptide.isin(set(main.peptide)))]

    #data.PGP_type[excl.alts.index] = "Alts_PEP"
    data.PGP_type[assigned_in_main.index] = "assigned_in_Main"
    return data
								


def assign_if_novel(data):
    original = data[data.seq_type=="original"]
    complement = data[(data.seq_type=="complement")]
    excl_complement = complement[~(complement.peptide.isin(set(original.peptide)))]
    original_assigned = complement[(complement.peptide.isin(set(original.peptide)))]
    data.PGP_type[excl_complement.index]="novel_CDS"
    data.PGP_type[original_assigned.index] ="assigned_in_forward_strand"
    return data

def assign_cds_start_end(mapped_data,gtf):
    mapped_data["has_cds"] = "No"
    mapped_data_copy= mapped_data.copy()
    mapped_data_assigned = mapped_data[mapped_data.PGP_type!="NA"]
    mapped_data_unassigned = mapped_data[mapped_data.PGP_type=="NA"]
    mapped_data = mapped_data_unassigned
   
    data = gtf[gtf.feature == "CDS"]
    #print data
    start_data = data.drop_duplicates(subset="trans_id")
    has_cds_trans_set= set(start_data.trans_id)
    end_data = data.drop_duplicates(subset="trans_id",keep="last")
    start_data = start_data.set_index("trans_id")
    end_data =end_data.set_index("trans_id")
    
    
    #print has_cds_trans_set
    mapped_data["has_cds"]= mapped_data.trans_id.apply(lambda x: "Yes" if x in has_cds_trans_set else "No")
    mapped_data = mapped_data[mapped_data.has_cds=="Yes"]
    #print mapped_data

    mapped_data["cds_start_on_chr"] = mapped_data.apply(lambda x: start_data["start"][x["trans_id"]] if x["strand"] =="+" else  start_data["end"][x["trans_id"]],axis=1)
    mapped_data["cds_end_on_chr"] = mapped_data.apply(lambda x: end_data["end"][x["trans_id"]] if x["strand"] =="+" else  end_data["start"][x["trans_id"]],axis=1)

    mapped_data_copy["has_cds"] = mapped_data.has_cds
    mapped_data_copy.has_cds[mapped_data_copy.has_cds.isnull()]="No"
    mapped_data_copy["cds_start_on_chr"] = mapped_data.cds_start_on_chr
    mapped_data_copy["cds_end_on_chr"] = mapped_data.cds_end_on_chr
  
   
    return mapped_data_copy

def assign_pgp_event_if_uorf(data,gtf):
    data_copy = data.copy()
    tmp = data[(data.strand=="+") & (data.has_cds=="Yes") & (data.PGP_type=="NA")]
    tmp.PGP_type = tmp.apply(lambda x: "uORF_PEP" if x["cds_start_on_chr"]>x["pep_start_on_chr"] else
                                                    "NA",axis=1)
    #print data.PGP_type
    data.PGP_type[tmp.index]=tmp.PGP_type
    #print data.PGP_type    

    tmp = data[(data.strand=="-") & (data.has_cds=="Yes") & (data.PGP_type=="NA")]
    tmp.PGP_type = tmp.apply(lambda x: "uORF_PEP" if x["cds_start_on_chr"]<x["pep_start_on_chr"] else "NA",axis=1)


    data.PGP_type[tmp.index] = tmp.PGP_type
    data_copy.PGP_type= data.PGP_type    
    
    return data_copy

def assign_pgp_event_if_dorf(data,gtf):
    data_copy = data.copy()

    tmp = data[(data.strand=="+") & (data.has_cds=="Yes") & (data.PGP_type=="NA")]
    tmp.PGP_type = tmp.apply(lambda x: "dORF_PEP" if x["pep_start_on_chr"]> x["cds_end_on_chr"] else
                                                    "NA",axis=1)
    #print data.PGP_type    
    data.PGP_type[tmp.index]=tmp.PGP_type
    #print data.PGP_type    

    tmp = data[(data.strand=="-") & (data.has_cds=="Yes") & (data.PGP_type=="NA")]
    tmp.PGP_type = tmp.apply(lambda x: "dORF_PEP" if x["pep_start_on_chr"]<x["cds_end_on_chr"] else "NA",axis=1)
    

    data.PGP_type[tmp.index] = tmp.PGP_type
    data_copy.PGP_type= data.PGP_type    
    

    return data_copy


def assign_pgp_event_if_cds_orf(data,gtf):
    data_copy = data.copy()

    tmp = data[(data.strand=="+") & (data.has_cds=="Yes") & (data.PGP_type=="NA")]
    tmp.PGP_type = tmp.apply(lambda x: "altCDS_PEP" if x["pep_start_on_chr"]> x["cds_start_on_chr"] 
                                                       and  x["pep_start_on_chr"]< x["cds_end_on_chr"] else
                                                    "NA",axis=1)
    #print data.PGP_type    
    data.PGP_type[tmp.index]=tmp.PGP_type
    #print data.PGP_type    

    tmp = data[(data.strand=="-") & (data.has_cds=="Yes") & (data.PGP_type=="NA")]
    tmp.PGP_type = tmp.apply(lambda x: "altCDS_PEP" if x["pep_start_on_chr"]<x["cds_start_on_chr"] 
                                                        and x["pep_start_on_chr"]> x["cds_end_on_chr"] else
                                                        "NA",axis=1)
    data.PGP_type[tmp.index] = tmp.PGP_type
    data_copy.PGP_type= data.PGP_type    
    

    return data_copy

def assign_pgp_event_if_exonic(data,gtf):
    data_copy = data.copy()

    assigned = data[~data.PGP_type.isin(set(["NA","assigned_in_forward_strand","semi-tryptic"]))]
    assigned_in_exon = data[(data.PGP_type=="NA") & data.peptide.isin(set(assigned.peptide))]
    data.PGP_type[assigned_in_exon.index] = "assigned_in_exon"
    unassigned = data[data.PGP_type=="NA"]
    #print unassigned
    data.PGP_type[unassigned.index] = unassigned.trans_type.apply(lambda x: x+"_PEP")

    #print data.PGP_type    
    data_copy.PGP_type = data.PGP_type
    print(len(data_copy.index))
    print(len(data.index))
    print("#################")
    return data_copy

def assign_all_events(data):
    #data_copy = data.copy()
    data["all_PGP_events"] = data.apply(lambda x: "|".join(set(data.PGP_type[data.peptide==x["peptide"]])),axis =1)
    return data    

def assign_pep_locus(data):

    data.pep_start_on_chr = data.pep_start_on_chr.apply(int)    
    data.pep_end_on_chr = data.pep_end_on_chr.apply(int)    

    data["all_pep_locus"] = data.apply(lambda x: x["chr"]+":"+str(x["pep_start_on_chr"])+"-"+str(x["pep_end_on_chr"])
        if x["strand"]=="+" else
        x["chr"]+":"+str(x["pep_start_on_chr"])+"-"+str(x["pep_end_on_chr"]),axis =1)
    data["assigned_pep_locus"] = "NA"
    exclude_set= set(["semi-tryptic","assigned_in_exon","assigned_in_forward_strand"]) 
    data.assigned_pep_locus[~data.PGP_type.isin(exclude_set)]=data[~data.PGP_type.isin(exclude_set)].apply(
        lambda x: x["chr"]+":"+str(x["pep_start_on_chr"])+"-"+str(x["pep_end_on_chr"])
        if x["strand"]=="+" else
        x["chr"]+":"+str(x["pep_start_on_chr"])+"-"+str(x["pep_end_on_chr"]),axis =1)
    data["all_locus_hits"] = "NA"
    data.all_locus_hits[data.assigned_pep_locus!="NA"] = data[data.assigned_pep_locus!="NA"].apply(
        lambda x: "#".join(set(data.assigned_pep_locus[(data.peptide==x["peptide"]) & (data.assigned_pep_locus!="NA")])),axis=1)
    data["total_locus_hits"] ="NA"
    data.total_locus_hits[data.all_locus_hits!="NA"] = data.all_locus_hits[data.all_locus_hits!="NA"].str.count("#")+ 1

    data["all_gene_hits"] ="NA"
    data.all_gene_hits[data.assigned_pep_locus!="NA"] = data[data.assigned_pep_locus!="NA"].apply(
        lambda x: "#".join(set(data.gene_name[(data.peptide==x["peptide"]) & (data.assigned_pep_locus!="NA")])),axis=1)

    data["total_gene_hits"] ="NA"
    data.total_gene_hits[data.all_locus_hits!="NA"] = data.all_gene_hits[data.all_locus_hits!="NA"].str.count("#")+1
    return data


def assign_Exonskip(pep_data,exonskip):
    pep_data.pgp[pep_data.database=="MCF7_exonskip1"] = "exonskip"
    return pep_data


def assign_PGP_event(data,mapped_data):
    mapped_data = assign_cds_start_end(mapped_data,data)
    mapped_data
    mapped_data = assign_pgp_event_if_uorf(mapped_data,data)
    mapped_data = assign_pgp_event_if_dorf(mapped_data,data)
    mapped_data = assign_pgp_event_if_cds_orf(mapped_data,data)
    #mapped_data = assign_pgp_event_if_exonic(mapped_data,data)
    mapped_data = assign_all_events(mapped_data)
    mapped_data = assign_pep_locus(mapped_data)
    return mapped_data

def assign_if_peptide_is_tryptic(data):
    data["before_pep_after"] = data.aa_res_before.str.cat([data.peptide,data.aa_res_after],sep="")
    data.before_pep_after = data.before_pep_after.str.strip("*").str.strip("-")
    data["tryptic_match"] = data.apply(lambda x: "Yes" if x["peptide"] 
                                in parser.cleave(x["before_pep_after"], parser.expasy_rules['trypsin'],2)
                                else "No",axis=1)                                                 
    return data

def assign_by_biotype(pep_data,data):
    trans_data = data[data.feature=="transcript"].set_index(["feature_id","strand"])
    pep_data.set_index(["parent_id","strand"],inplace=True)
#    transcripts_with_cds = data[data.feature=="CDS"].drop_duplicates(["parent_id","strand"]).set_index(["parent_id","strand"])
#    pep_data["has_cds"] = pep_data.index.isin(transcripts_with_cds.index)
#    first_cds = data[data.feature=="CDS"].drop_duplicates(["parent_id","strand"]).set_index(["parent_id","strand"])
#    last_cds = data[data.feature=="CDS"].drop_duplicates(["parent_id","strand"],keep="last").set_index(["parent_id","strand"])
#    pep_data["gene_biotype"] = trans_data.gene_biotype
    pep_data["trans_type"] = trans_data.trans_type
#    pep_data["cds_start_on_chr"] = first_cds.start
#    pep_data["cds_end_on_chr"] = first_cds.end
    pep_data["gene_type"] = trans_data.gene_type
    #pep_data["cds_end"] = last_cds.end
    pep_data.reset_index(inplace=True)
    print(pep_data)
    return pep_data



def read_gff(gff_file_path):
    print("reading gff...")
    data=pd.read_csv(gff_file_path,header=None,comment="#",sep="\t",
                        names=["chr","source","feature","start","end","score","strand","frame","attributes"])
    data.attributes= data.attributes+";"
    data["feature_id"] = data.attributes.str.extract('ID=(.*?);')
    data["parent_id"] = data.attributes.str.extract('Parent=(.*?);')
    data["gene_name"] = data.attributes.str.extract('gene_name=(.*?);')
    data["gene_type"] = data.attributes.str.extract(';gene_type=(.*?);')
    data["trans_type"] = data.attributes.str.extract(';transcript_type=(.*?);')
    data["trans_id"] = data.attributes.str.extract(';transcript_id=(.*?);')
    
    return data

def get_gnomon_mapped_bed_files(file_paths):
    frames=[]
    for file_path in file_paths:
        data = pd.read_csv(file_path)
        if file_path.endswith("_ref_bed.csv"):
            data["assembly"] = "ref"
        else:
            data["assembly"] = "alt"
        frames.append(data)
    data = pd.concat(frames,ignore_index=True)
    return data

def process_gnomon(args):
    gnomon_folder = args.gnomon
    ref_gtf_file_path = os.path.join(gnomon_folder,"ref_GRCh38.p7_top_level.gff3")
    alt_gtf_file_path = os.path.join(gnomon_folder,"alt_CHM1_1.1_top_level.gff3")
    gencode_gtf_file_path = "D:/Data/genome/gencode/gencode.v27.chr_patch_hapl_scaff.annotation.gtf"
    chr_map_file = os.path.join(gnomon_folder,"chr_accessions_GRCh38.p7")
    alt_assembly_chr_map_file = os.path.join(gnomon_folder,"chr_accessions_CHM1_1.1")
    mapped_data = get_gnomon_mapped_bed_files(args.maps)
    mapped_data = assign_if_peptide_is_tryptic(mapped_data)
    mapped_data = assign_refseq_gene_and_pgp(mapped_data,ref_gtf_file_path,
                                    alt_gtf_file_path,gencode_gtf_file_path)
    mapped_data = assign_if_in_ref(mapped_data)
    mapped_data["chrom_type"] ="NA"
    mapped_data = assign_if_in_main_chroms(mapped_data,chr_map_file,"ref")
    mapped_data = assign_if_in_main_chroms(mapped_data,alt_assembly_chr_map_file,"alt")

    mapped_data = assign_novel_isoform(mapped_data) 
    mapped_data = assign_novel_CDS_from_gencode(mapped_data,gencode_gtf_file_path,chr_map_file)
    mapped_data = assign_novel_isoform_gencode(mapped_data)

    mapped_data.to_csv(args.out,index=False)


def main():
    args = get_arguments()
    if args.source=="gnomon":
        process_gnomon(args)
    else:

        pep_data  = pd.read_csv(args.maps)
        pep_data["PGP_type"] ="NA"
        pep_data = assign_if_peptide_is_tryptic(pep_data)
        pep_data = assign_if_novel(pep_data)
        pep_data["trans_id"] = pep_data.parent_id
        #sys.exit()
        gff = read_gff(args.gff)
        pep_data = assign_by_biotype(pep_data,gff)
        pep_data = assign_PGP_event(gff,pep_data)
        pep_data.to_csv(args.out,index=False)
                

def get_arguments():
    parser = ArgumentParser()
    parser.add_argument("-maps",required=True,help="full path to bed mapped files..")
    parser.add_argument("-event", choices=["protein","transcript","gene","genome"], help="type of event")
    parser.add_argument("-source", choices=["gencode","refseq","gnomon"], help="source")
    parser.add_argument("-gff",help="full path to gff file..")
    parser.add_argument("-out",help="full path to gff file..")
    parser.add_argument("-gnomon",default=None,help="Path to folder containing annotations..")
   
    args = parser.parse_args()
    return args

if __name__=="__main__":
    main()