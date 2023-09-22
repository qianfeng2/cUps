A probabilistic method for classifying the malaria *var* genes into ups groups 
-----------------------
[![Python 3.6](https://img.shields.io/pypi/pyversions/Django)](https://www.python.org/downloads/release/python-360/)
[![License: GPL-3.0](https://img.shields.io/cran/l/devtools)](https://opensource.org/licenses/GPL-3.0)

### About
Our **cUps** algorithm is for **c**lassifying the malaria *var* genes into **ups** groups using DBLα sequences. It takes as input a reference database of DBLα sequences with known ups groups and a set of DBLα sequences to be classified, and outputs the probabilities of membership to all three ups groups per sequence.




### Input and output 
- Fasta format DBLα sequences to be classified (refer to query sequences). 


- Reference data which consist of (1) reference DBLα sequences, each sequence is annotated with DBLα subclass and ups group; (2) profile HMM for each reference category (combination of DBLα subclass and ups group). Please see [reference_data](https://github.com/qianfeng2/cUps/tree/main/reference_data) folder for details.


Our algorithm produces four csv files based on various algorithmic stages, and places them in the directory specified by output. 


- Log-likelihood tables (middle output). For each query sequence, its likelihood for being drawn from the profile HMM under each category is computed.


- Classification result table (final output). Each row starts with the sequence ID and shows the probabilities of membership to three ups groups. An example table is as following.  

|         | A  | B  | C  | 
| ------------|------------|------------|------------|
|seq1 | 0.98|0.01|0.01|
|seq2 | 0.01|0.04|0.95|
|... | ... |... |... |


### Required softwares
- R
- Python >=3.5


### Usage

```
git clone https://github.com/qianfeng2/cUps

cd cUps

pip install -r requirements.txt

mkdir your_output_dir 

python scripts/generate_llk.py query_data/example.fasta your_output_dir

Rscript scripts/classify_upsABC.R your_output_dir
```

The python script is to generate the log-likelihood of query sequences for every profile HMM of reference categories. This script outputs three tables (PAD.csv, PBD.csv and PCD.csv) representating the results under upsA, upsB and upsC group. For every table, each row refers to a query sequence, each column refers to a DBLα subclass, the entry is the log-likelihood value.

The R script is to calculate the posterior probabilities for each query sequence. The reference data is replaceable. The final result is shown at classification_result.csv.

As a toy example, [results](https://github.com/qianfeng2/cUps/tree/main/results) folder provides all the middle and final output files for the example.fasta stored in [query_data](https://github.com/qianfeng2/cUps/tree/main/query_data) folder. 


### Run example for large number of sequences (>1000 sequences)

Since our algorithm processes each sequence independenly, splitting the big dataset into subsets is recommended. For each subset, you can run our algorithm with the help of HPC.

Below I show steps how to split an example data with 1000 sequences and run the subsets simultaneously.

```
cd query_data

mkdir example_bigdata_split_files

cd ../

python scripts/split_input.py query_data/example_bigdata.fasta query_data/example_bigdata_split_files
```

```
cd results

mkdir example_bigdata

cd example_bigdata

for run in $(seq 1 1 4);do mkdir run_$run;done
```

I use job array in HPC to run the subsets in parallel.

```
python scripts/generate_llk.py query_data/example_bigdata_split_files/input_run${SLURM_ARRAY_TASK_ID}.fasta results/example_bigdata/run_${SLURM_ARRAY_TASK_ID}

Rscript scripts/classify_upsABC.R results/example_bigdata/run_${SLURM_ARRAY_TASK_ID}
```


### Note

- Our algorithm only accepts the protein sequences as input. If your data is DNA, please translate it.

- For the format of identifier in each protein sequence, althouth the blank space is not allowed inside the identifier, punctuations (e.g., "|", ";", "_", "-") generally do not break our algorithm.

- Please run your own data inside the cUps directory, as our algorithm relies on the information from reference_data directory. 


### Credits

This algorithm is developed by Qian Feng, Heejung Shim and Yao-ban Chan at the University of Melbourne. For any problems, please report an issue in Github or send an email to [Qian Feng](mailto:fengq2@student.unimelb.edu.au).



### Reference

- Thomas S Rask, Daniel A Hansen, Thor G Theander,Anders Gorm Pedersen, and Thomas Lavstsen. *Plasmodium falciparum* erythrocyte membrane protein 1 diversity in seven genomes – divide and conquer. *PLoS Comput. Biol.*, 6(9):e1000933, 2010

- Anders Krogh, Michael Brown, I Saira Mian, Kimmen Sjölander, and David Haussler. Hidden Markov models in computational biology: Applications to protein modeling. *J. Mol. Biol.*, 235(5):1501–1531, 1994

- Richard Durbin, Sean R Eddy, Anders Krogh, and Graeme Mitchison. *Biological sequence analysis: probabilistic models of proteins and nucleic acids.* Cambridge university press, 1998



