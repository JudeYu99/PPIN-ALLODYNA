#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@ author: Yu Zhu

@ Email: yzhu99@stu.suda.edu.cn

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

"""

### Introduction of MI_coevolution_v2.py
#
# @ This part of program is delicated for calculating conservation and co-evolution score based on merged MSAs.
# @ This program has been updated to version 2 now.
# @ Two separate MSA files are required when running this sub-program.
# @ A merged MSA file will be created in fasta format.
# @ Reference: Liu Y, Bahar I. Sequence Evolution Correlates with Structural Dynamics 2012 Mol Biol Evol 29(9):2253-2263; Liu Y, Gierasch LM, Bahar I Role of Hsp70 ATPase domain intrinsic dynamics and sequence evolution in enabling its functional interactions with NEFs 2010 PLoS Comput Biol 6(9)
#
# @ Python package in need: prody, numpy, pandas, re, time, os
#
#############################################

from prody import *
import numpy as np
import pandas as pd
import re
import time
import os


def parse_fasta(fasta_file_name):
    fasta_file = open(fasta_file_name, "r")
    label_list = []
    seq_list = []

    is_first_seq = True
    seq_i = ""
    for line in fasta_file:
        line = line.strip()
        if line == "":
            continue
        if line[0] == ">":
            is_new_label = True
            label_list.append(line)
            if is_first_seq:
                is_first_seq = False
                seq_i = ""
            else:
                seq_list.append(seq_i)
                seq_i = ""
        else:
            seq_i += line
    seq_list.append(seq_i)
    fasta_file.close()
    label_list_list = []
    dict_fasta = dict()

    key_list_2 = ['OS', 'OX', 'GN', 'PE', 'SV']
    for index_i in range(len(label_list)):
        label_one = label_list[index_i]
        label_one_list = label_one[1:].split("|")
        temp_list = re.split(' OS=| OX=| GN=| PE=| SV=', label_one_list[2])
        match_items = re.findall(' OS=| OX=| GN=| PE=| SV=', label_one_list[2])
        match_items = [x[1:3] for x in match_items]
        if key_list_2 != match_items:
            dict_key_temp = dict()
            for key_i in key_list_2:
                dict_key_temp[key_i] = ''
            for i, key_i in enumerate(match_items, start=1):
                dict_key_temp[key_i] = temp_list[i]
            temp_list_use = list(dict_key_temp.values())
        else:
            temp_list_use = temp_list[1:]
        uniprotKB_proteinName_list = temp_list[0].split(" ")
        OX_GN = temp_list_use[1]+'_'+temp_list_use[2]
        OX_GN = OX_GN.upper()
        label_i_list = label_one_list[0:2] + [uniprotKB_proteinName_list[0]] + [" ".join(uniprotKB_proteinName_list[1:])] + temp_list_use + [OX_GN]
        label_list_list.append(label_i_list)
        dict_fasta[label_one_list[1]] = seq_list[index_i]
    
    return [dict_fasta, label_list_list]



def write_aligned_fasta(fasta_name_01, fasta_name_02):
    
    label_key_list = ['db', 'UniProtID', 'UniProtKB', 'proteinName', 'OS', 'OX', 'GN', 'PE', 'SV', 'OX_GN']
    fasta_01 = parse_fasta(fasta_name_01)
    fasta_02 = parse_fasta(fasta_name_02)
    df_label_01 = pd.DataFrame(list(fasta_01[1]), columns=label_key_list)
    df_label_02 = pd.DataFrame(list(fasta_02[1]), columns=label_key_list)
    col_use = 'OX'
    # drop duplicated ones with multiple OX (or OX_GN) values, and keep the first occurrence (sp is preferred: sorted by db)
    df_label_01.sort_values("db", inplace = True) 
    df_label_01.drop_duplicates(col_use, inplace = True)
    df_label_02.sort_values("db", inplace = True) 
    df_label_02.drop_duplicates(col_use, inplace = True)
    # drop rows with empty GN values
    df_label_01 = df_label_01[df_label_01['GN'] != '']
    df_label_02 = df_label_02[df_label_02['GN'] != '']
    OX_01_list = df_label_01[col_use].to_list()
    OX_02_list = df_label_02[col_use].to_list()
    label_01_02 = list(set(OX_01_list) & set(OX_02_list))
    df_label_01_02 = df_label_01.loc[df_label_01["OX"].isin(label_01_02)]
    df_label_01_02_sort = df_label_01_02.sort_values("OX")
    df_label_02_01 = df_label_02.loc[df_label_02["OX"].isin(label_01_02)]
    df_label_02_01_sort = df_label_02_01.sort_values("OX")
    
    clustal_name_01 = "./BLASTS/" + fasta_name_01.split("/")[-1][:6] + "_" + fasta_name_02.split("/")[-1][:6] + "_blast.fasta"
    fasta_file = open(clustal_name_01, "w")
    for index, row_i in df_label_01_02_sort.iterrows():
        fasta_file.write(">" + row_i["OX"] + "\n")
        fasta_file.write(fasta_01[0][row_i["UniProtID"]] + "\n")
    fasta_file.close()
    
    clustal_name_02 = "./BLASTS/" + fasta_name_02.split("/")[-1][:6] + "_" + fasta_name_01.split("/")[-1][:6] + "_blast.fasta"
    fasta_file = open(clustal_name_02, "w")
    for index, row_i in df_label_02_01_sort.iterrows():
        fasta_file.write(">" + row_i["OX"] + "\n")
        fasta_file.write(fasta_02[0][row_i["UniProtID"]] + "\n")
    fasta_file.close()


def MI_coevolution(msa_file_01, msa_file_02):

    # Load MSA files.
    msa_01 = parseMSA(msa_file_01)
    msa_02 = parseMSA(msa_file_02)

    # Merge MSA.
    msa = mergeMSA(msa_01, msa_02)
    
    print("\n> *** Merge MSAs finished! *** <\n")

    # Refine MSAs.
    msa_refine_01 = refineMSA(msa_01, label = "9606", rowocc = 0.4, seqid = 0.98)
    msa_refine_02 = refineMSA(msa_02, label = "9606", rowocc = 0.4, seqid = 0.98)
    msa_refine = refineMSA(msa, label = "9606", rowocc = 0.4, seqid = 0.98)

    # Get sequence length.
    protein_01_len = msa_refine_01.numResidues()
    protein_02_len = msa_refine_02.numResidues()
    protein_len = msa_refine.numResidues()

    # Write merged MSA files.
    # writeMSA('msa_aligned.fasta', msa_refine)

    # Co-evolution calculation with mutual information.
    mutinfo = buildMutinfoMatrix(msa_refine)
    mutinfo_inter_protein = mutinfo[:protein_01_len, protein_01_len:]
    MI_score = np.mean(mutinfo_inter_protein)

    print("\n> *** MI coevolution finished! *** <\n")

    return MI_score


def RUN(UniProtID_01, UniProtID_02):

    fasta_01 = "./FASTAS/" + UniProtID_01 + ".fasta"  
    fasta_02 = "./FASTAS/" + UniProtID_02 + ".fasta"  
    write_aligned_fasta(fasta_01, fasta_02)
    
    time.sleep(2)
    
    temp_blast_name_01 = "./BLASTS/" + UniProtID_01 + "_" + UniProtID_02
    temp_blast_name_02 = "./BLASTS/" + UniProtID_02 + "_" + UniProtID_01
    
    os.system("clustalo -i " + temp_blast_name_01 + "_blast.fasta" + " -o " + temp_blast_name_01 + "_aligned_blast.fasta")
    os.system("clustalo -i " + temp_blast_name_02 + "_blast.fasta" + " -o " + temp_blast_name_02 + "_aligned_blast.fasta")
    
    time.sleep(2)
    
    msa_file_01 = temp_blast_name_01 + "_aligned_blast.fasta"
    msa_file_02 = temp_blast_name_02 + "_aligned_blast.fasta"
            
    MI_score = MI_coevolution(msa_file_01, msa_file_02)
    
    print("The coevolution score for " + UniProtID_01 + " and " + UniProtID_02 + " is: " + str(MI_score))
    
    return (UniProtID_01, UniProtID_02, MI_score)


if __name__ == '__main__':

    ID_01 = []
    ID_02 = []
    MI = []

    PPIN = pd.read_table("UniProt_PPIN.txt", sep = "\t")
    for i in range(PPIN.shape[0]):
        UniProtID_01 = PPIN.iloc[0,:][0]
        UniProtID_02 = PPIN.iloc[0,:][1]

        temp_res = RUN(UniProtID_01, UniProtID_02)
        ID_01.append(UniProtID_01)
        ID_02.append(UniProtID_02)
        MI.append(temp_res[2])

    df = pd.DataFrame()
    df["UniProtID_01"] = ID_01
    df["UniProtID_02"] = ID_02
    df["MI_score"] = MI

    df.to_csv("./PPIN_MI.csv", index = False, sep = "\t")

'''

    UniProtID_01 = "Q9Y3B2"
    UniProtID_02 = "Q9Y3B7"

'''
