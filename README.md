# PPIN-ALLODYNA
## 1. Introduction of PPIN-ALLODYNA (PPIN-Allo-Dynamics)

At the structural level, the computation of biological networks, including protein structural networks and elastic network models developed in our preliminary research, has been widely applied to allosteric signaling within proteins, providing effective tools for structural biologists. In particular, the elastic network model can not only study the local energy and conformational changes of protein structures at the molecular level, but also realize the quantitative description of protein dynamics at the omics level, providing a theoretical approach for high-throughput allostery studies. Therefore, extending the method based on the elastic network model calculation to the protein-protein interaction network level can realize the global dynamics of the network or module with the functional coupling of the signaling pathways inherent to the network topology, and theoretically realize the allosteric effect at the network level.
  
**PPIN-ALLODYNA** (**P**rotein-**P**rotein **I**nteraction **N**etwork-**ALLO**steric **DYNA**mics) is an integrated tool for exploring methods of introducing allosteric effects based on protein structure level into the analysis of protein-protein interaction networks to construct a PPIN-based functional atlas. This tool is expected for mapping biophysical ideas and methods into the study of protein-protein interaction networks, enabling quantitative and dynamic portrayal of signaling pathways, as well as anticipating new interpretations of molecular mechanisms of disease.
  
The tool is currently being further developed and refined, and a new version will be available soon in the future.

## 2. Usage
  1) Python packages in need: os, re, string, sys, time, numpy, scipy, pandas, prody, wget, selenium, networkx
  2) External software to install: [ChromeDriver](https://sites.google.com/a/chromium.org/chromedriver/home)
  3) Python language version: >=3.7
  4) Before running ***Run.py***, please modify line 41 for PPIN input, line 116 and line 117 for two interacting proteins' UniProt IDs.
  ```
  python3 Run.py
  ```

## 3. Examples

- Example PPIN input file can be referred to [***test_PPIN.txt***](https://github.com/JudeYu99/PPIN-ALLODYNA/blob/main/test_PPIN.txt).

- Example output of Hitting Time and Commute Time:  
  <div align=center>
  <img src="https://github.com/JudeYu99/PPIN-ALLODYNA/blob/main/OUTPUTS/hit.png" width="380" height="380"><img src="https://github.com/JudeYu99/PPIN-ALLODYNA/blob/main/OUTPUTS/commute.png" width="380" height="380">
  </div>
  
- Example output of PRS:  
  <div align=center>
  <img src="https://github.com/JudeYu99/PPIN-ALLODYNA/blob/main/OUTPUTS/prs_heatmap.png" width="400" height="400"/>
  </div>

- Example output of co-evolution score ([TP53](https://www.uniprot.org/uniprot/P04637) and [CDK2](https://www.uniprot.org/uniprot/P24941)):  
  ```
  The coevolution score for P04637 and P24941 is: 0.06319592447541351
  ```
  The intermediate files obtained during the calculation are detailed in [**OUTPUTS**](https://github.com/JudeYu99/PPIN-ALLODYNA/tree/main/OUTPUTS) folder.
  
  For questions and data requirements, please contact ***yzhu99@stu.suda.edu.cn***.

  
## 4. Reference
- Stephen F. Altschul, Thomas L. Madden, Alejandro A.Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J.Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.
- Liu, Y., & Bahar, I. (2012). Sequence evolution correlates with structural dynamics. Molecular biology and evolution, 29(9), 2253–2263.
- Liu, Y., Gierasch, L. M., & Bahar, I. (2010). Role of Hsp70 ATPase domain intrinsic dynamics and sequence evolution in enabling its functional interactions with NEFs. PLoS computational biology, 6(9), e1000931.
- Madeira F, Park YM, Lee J, et al. The EMBL-EBI search and sequence analysis tools APIs in 2019. Nucleic Acids Research. 2019 Jul;47(W1):W636-W641. 
- UniProt BLAST (https://www.uniprot.org/blast/)
- Clustal Omega Web Server (https://www.ebi.ac.uk/Tools/msa/clustalo/)
- Harris, C. R., Millman, K. J., van der Walt, S. J., Gommers, R., Virtanen, P., Cournapeau, D., … Oliphant, T. E. (2020). Array programming with NumPy. Nature, 585, 357–362.
- Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T., Cournapeau, D., … SciPy 1.0 Contributors. (2020). SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17, 261–272. 
- McKinney, W., & others. (2010). Data structures for statistical computing in python. In Proceedings of the 9th Python in Science Conference (Vol. 445, pp. 51–56).
- Hagberg, A., Swart, P., & S Chult, D. (2008). Exploring network structure, dynamics, and function using NetworkX.
