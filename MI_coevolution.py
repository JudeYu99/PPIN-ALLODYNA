#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@ author: Yu Zhu

@ Email: yzhu99@stu.suda.edu.cn

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

"""

### Introduction of conservation_coevolution.py
#
# @ This part of program is delicated for calculating conservation and co-evolution score based on merged MSAs.
# @ Two separate MSA files are required when running this sub-program.
# @ A merged MSA file will be created in fasta format.
# @ Reference: Liu Y, Bahar I. Sequence Evolution Correlates with Structural Dynamics 2012 Mol Biol Evol 29(9):2253-2263; Liu Y, Gierasch LM, Bahar I Role of Hsp70 ATPase domain intrinsic dynamics and sequence evolution in enabling its functional interactions with NEFs 2010 PLoS Comput Biol 6(9)
#
# @ Python package in need: prody, numpy, pandas
#
#############################################

from prody import *
import numpy as np
import pandas as pd


def MI_coevolution(msa_file_01, msa_file_02, prefix = "./OUTPUTS/protein_protein"):

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
    writeMSA(prefix + '_msa_aligned.fasta', msa_refine)

    # Co-evolution calculation with mutual information.
    mutinfo = buildMutinfoMatrix(msa_refine)
    mutinfo_inter_protein = mutinfo[:protein_01_len, protein_01_len:]
    MI_score = np.mean(mutinfo_inter_protein)

    print("\n> *** MI coevolution finished! *** <\n")

    return MI_score

