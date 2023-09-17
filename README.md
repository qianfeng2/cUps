A probabilistic method for classifying ups group of the malaria parasite var genes
-----------------------
[![Python 3.6](https://img.shields.io/pypi/pyversions/Django)](https://www.python.org/downloads/release/python-360/)
[![License: GPL-3.0](https://img.shields.io/cran/l/devtools)](https://opensource.org/licenses/GPL-3.0)
### About
This algorithm is .


### Required softwares
- R
- Python >=3.5
```
pip install -r requirements.txt
```


### Input and output 
- Input fasta format biological sequences, maximum length for identifiers length is 15 (Zilversmit et al., 2013)
- Patial alignment produced by JHMM (please see [MZmosaic](https://github.com/qianfeng2/detREC_program/tree/master/MZmosaic) sub folder)

Produces a series of files based on various stage of the implementation of recombination detection program, and places them in the directory specified by output.

- temp file folder  
This folder provides all the chunks containing original triple and MAFFT processed fasta files.
- complement_chunks file folder  
This folder provides all the equal-length triples, name of each file indicates chunk index, two adjacent segment indices, and identified bkp in this triple.
- output csv file:  
Each row records the chunk index in partial alignment result, target, db1 and db2 are three sequences ID for each triple, rec is the identified recombinant ID from this specific triple, sv is the support value. For instance:  

| chunk        | target  | db1  | db2  | rec  | sv  |
| ------------|------------|------------|------------|------------|------------|
|2 | seq3|seq5|seq2|seq3|1|
|6 | seq7|seq8|seq1|seq8|0.58|
|... | ... |... |... |... |... |


### Usage

```
cd project_dir

python scripts/generate_llk.py query_data/example.fasta results

Rscript scripts/classify_upsABC.R results
```


### Credits

[Test_files](https://github.com/qianfeng2/detREC_program/tree/master/Test_files) sub folder, as a toy example, provides a test input.fasta and all the middle and final output files.



### Reference
- Zilversmit, M. M., Chase, E. K., Chen, D. S., Awadalla, P., Day, K. P., & McVean, G. (2013). Hypervariable antigen genes in malaria have ancient roots. BMC evolutionary biology, 13(1), 110.
- Fletcher, W., & Yang, Z. (2009). INDELible: a flexible simulator of biological sequence evolution. Molecular biology and evolution, 26(8), 1879-1888.
- Spielman, S. J., & Wilke, C. O. (2015). Pyvolve: a flexible Python module for simulating sequences along phylogenies. PloS one, 10(9), e0139047.
