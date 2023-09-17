A probabilistic method for classifying ups groups of the malaria *var* genes
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
[Test_files](https://github.com/qianfeng2/detREC_program/tree/master/Test_files) sub folder, as a toy example, provides a test input.fasta and all the middle and final output files.

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
git clone https://github.com/qianfeng2/cUps

cd cUps

mkdir your_output_dir 

python scripts/generate_llk.py query_data/example.fasta your_output_dir

Rscript scripts/classify_upsABC.R your_output_dir
```


### Credits

This algorithm is developed by Qian Feng, Heejung Shim and Yao-ban Chan at the University of Melbourne. For any problems, please report an issue in Github or send an email to [Qian Feng](mailto:fengq2@student.unimelb.edu.au)



### Reference

- Thomas S Rask, Daniel A Hansen, Thor G Theander,Anders Gorm Pedersen, and Thomas Lavstsen. *Plasmodium falciparum* erythrocyte membrane protein 1 diversity in seven genomes – divide and conquer. PLoS Comput. Biol., 6(9):e1000933, 2010

- Anders Krogh, Michael Brown, I Saira Mian, Kimmen Sjölander, and David Haussler. Hidden Markov models in computational biology: Applications to protein modeling. J. Mol. Biol., 235(5):1501–1531, 1994

- Richard Durbin, Sean R Eddy, Anders Krogh, and Graeme Mitchison. Biological sequence analysis: probabilistic models of proteins and nucleic acids. Cambridge university press, 1998



