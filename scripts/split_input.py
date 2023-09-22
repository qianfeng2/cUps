#######################################################################
# Copyright (C) 2023  Qian Feng

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Usuage instruction
# This is used in Python 3.
# This code is to split your input fasta into smaller subset, each with 300 sequences.
# The first parameter is your query sequence as the input, second parameter is in the dir where your split files locate
# Usage example: cd your_dir \\
#                python scripts/split_input.py query_data/example_bigdata.fasta dir 
#######################################################################

from collections import defaultdict
from Bio import SeqIO
import sys, os
import numpy
import glob


input_fasta=sys.argv[1]#query_data/example_bigdata.fasta
output_dir=sys.argv[2]#query_data/split_files

seqs = {}
seq_list = []

for seq_record in SeqIO.parse(input_fasta, "fasta"):
    seqs[str(seq_record.id)] = str(seq_record.seq)
    seq_list.append(str(seq_record.id))
    
n=300;run=0
for i in range(0, len(seq_list), n):
    targets = seq_list[i:i+n]
    run+=1
    with open(output_dir+"/input_run" + str(run) + ".fasta", 'w') as outfile:
        for t in targets:
            outfile.write(">"+t+"\n"+seqs[t]+"\n")






