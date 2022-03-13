#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@ author: Yu Zhu & Ziyun Zhou

@ Email: yzhu99@stu.suda.edu.cn

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

"""

##############################################################################
#
#   @> Step 0. Import Python packages in need
#
##############################################################################

import os
import time
import numpy as np
import networkx as nx
import pandas as pd
from enm.Enm import *
from enm.utils import *
from enm.visualize import *
from HitCommute import *
from get_BLAST_from_UniProt import *
from align_fasta import *
from get_aligned_MSA_from_Clustal_Omega import *
from MI_coevolution import *


##############################################################################
#
#   @> Step 1. Generate ENM from PPI netwrok file
#
##############################################################################

network_file = "test_PPIN.txt"

enm = Enm('PPIN')
enm.read_network(network_file, sep = '\t')
enm.gnm_analysis(normalized = False)

os.mkdir("./OUTPUTS")

##############################################################################
#
#   @> Step 2. Hitting Time and Commute Time Calculation
#
##############################################################################

K = enm.gnm.getKirchhoff()

# Pass the Kirchhoff matrix to hit/commute time 
hc = IT_HitCommute(K)

# Report hitting time
H = hc.buildHitTimes(K)
np.savetxt('./OUTPUTS/hit_df.txt', H)

# Report commute time
C = hc.buildCommuteTimes()
np.savetxt('./OUTPUTS/commute_df.txt', H)

# Plot hitting time
plt.rcParams['figure.figsize'] = (10, 10)
plt.rcParams['font.size'] = 12
plt.title("Hitting Time of PPIN")
plt.imshow(H, origin = "lower", cmap = "rainbow")
plt.colorbar(shrink = 0.8)
plt.savefig("./OUTPUTS/hit.png")
plt.show()

# Plot commute time
plt.rcParams['figure.figsize'] = (10, 10)
plt.rcParams['font.size'] = 12
plt.title("Commute Time of PPIN")
plt.imshow(C, origin = "lower", cmap = "rainbow")
plt.colorbar(shrink = 0.8)
plt.savefig("./OUTPUTS/commute.png")
plt.show()

print("\n> *** Hitting time and commute time calculation finished! *** <\n")


##############################################################################
#
#   @> Step 3. PRS Calculation
#
##############################################################################

enm.get_sensor_effector(use_threshold = True)
enm.cluster_matrix(enm.prs_mat)
neighbor_degree = []
for i in enm.graph_gc.nodes:
    neighbor_degree.append(np.average([enm.df.loc[enm.df.orf_name == a,'deg'].values for a in nx.neighbors(enm.graph_gc,i)]))
enm.df['neighbor_degree'] = neighbor_degree

enm.df.to_csv("./OUTPUTS/pcc_df.csv", index=True, index_label = 'orf_name_id')
np.savetxt('./OUTPUTS/prs_df.txt', enm.prs_mat)
np.savetxt('./OUTPUTS/prs_mat_df.txt', enm.prs_mat_df)
enm.heatmap_annotated(save_figure = True)

print("\n> *** PRS calculation finished! *** <\n")


##############################################################################
#
#   @> Step 4. PPI Co-evolution Score
#
##############################################################################

UniProtID_01 = "P04637"
UniProtID_02 = "P24941"

BLAST_01 = BLAST_MSA(UniProtID_01, "01")
BLAST_02 = BLAST_MSA(UniProtID_02, "02")
        
fasta_01 = "./OUTPUTS/01_blast.fasta"  
fasta_02 = "./OUTPUTS/02_blast.fasta"
write_aligned_fasta(fasta_01, fasta_02)
    
Clustal_Omega("./OUTPUTS/01_02_blast.fasta")
Clustal_Omega("./OUTPUTS/02_01_blast.fasta")

time.sleep(3)
msa_file_01 = "./OUTPUTS/01_02_aligned_blast.fasta"
msa_file_02 = "./OUTPUTS/02_01_aligned_blast.fasta"
        
MI_score = MI_coevolution(msa_file_01, msa_file_02, prefix = "./OUTPUTS/protein_protein")

print("The coevolution score for " + UniProtID_01 + " and " + UniProtID_02 + " is: " + str(MI_score))
