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
# This code is to calculate the log-likelihood of the query sequence being drawn from the profile HMM of each category, using the forward algorithm.
# First input parameter is your query sequence
# Second input parameter is your preferred output directory 
# Usage example: cd your_project_dir \\
#                python scripts/generate_llk.py query_data/example.fasta results
#######################################################################


from collections import defaultdict
from pomegranate import *
from Bio import SeqIO
import sys, os
import numpy
import math
import csv
import glob

input_fasta=sys.argv[1]#e.g. query_data/example.fasta
output_dir=sys.argv[2]#e.g results
output_PAD_csv=output_dir+"/PAD.csv"
output_PBD_csv=output_dir+"/PBD.csv"
output_PCD_csv=output_dir+"/PCD.csv"


def read_PHMM( PHMM ):
    """
    This script will take in .hmm file stroring parameters of Profile HMM. 
    This script returns pomegranate.hmm.HiddenMarkovModel object.
    """
    def convert_probability(p: str) -> float:
        if p == '*':
            return 0
        else:
            return math.e ** -float(p)
    def normalize_probability(p: list) -> list:
        if len(p)==4:
            normalized_p=[p[m]/sum(p) for m in range(4) ]
        elif len(p)==20:
            normalized_p=[p[m]/sum(p) for m in range(20) ]            
        else:    
            normalized_p=list(range(7));
            normalized_p[0]=p[0]/(p[0]+p[1]+p[2]);
            normalized_p[1]=p[1]/(p[0]+p[1]+p[2]);
            normalized_p[2]=p[2]/(p[0]+p[1]+p[2]);
            normalized_p[3]=p[3]/(p[3]+p[4]);
            normalized_p[4]=p[4]/(p[3]+p[4]);
            normalized_p[5]=p[5]/(p[5]+p[6]);
            normalized_p[6]=p[6]/(p[5]+p[6]);
        return normalized_p
    metadata_dict={}
    with open(PHMM, 'r') as f:
        Par_startline_index = 0;Par_endline_index=0
        for i, line in enumerate(f.readlines()):
            if i<=1:continue
            line = line.strip().split()
        #print(line)
            if line[0] == 'LENG':
                metadata_dict['length'] = int(line[1])
            elif line[0] == 'ALPH':
                metadata_dict['alphabet_type'] = line[1]        
            elif line[0] == 'HMM':
                metadata_dict['alphabet'] = line[1:]
            elif line[0] == 'COMPO':
                Par_startline_index = i
            elif line[0] == '//':
                Par_endline_index = i            
    match_line_index=list(range(Par_startline_index,Par_endline_index,3))
    insertion_line_index=list(range(Par_startline_index+1,Par_endline_index,3))
    state_switch_line_index=list(range(Par_startline_index+2,Par_endline_index,3))
    with open(PHMM, 'r') as f:
        insertion_line={};match_line={};state_switch_line={}
        temp1=0;temp2=0;temp3=0;
        for i, line in enumerate(f.readlines()):
            line = line.strip().split()
            if i in insertion_line_index:
                insertion_line[temp1]=line[0:len(metadata_dict['alphabet'])]
                temp1=temp1+1
            elif i in state_switch_line_index:
                state_switch_line[temp2]=line[0:7]
                temp2=temp2+1
            elif i in match_line_index:
                match_line[temp3]=line[1:(len(metadata_dict['alphabet'])+1)]
                temp3=temp3+1
                
    
    # Define all hidden states name
    var_length=metadata_dict['length']+1
    b=["i"]*var_length;a=list(range(var_length));B=["I"]*var_length;a=list(range(var_length))
    Insert_states=[m+str(n) for m,n in zip(b,a)];Insert_states_name=[m+str(n) for m,n in zip(B,a)]
    b=["m"]*(var_length-1);a=list(range(1,var_length));B=["M"]*(var_length-1);a=list(range(1,var_length))
    Match_states=[m+str(n) for m,n in zip(b,a)];Match_states_name=[m+str(n) for m,n in zip(B,a)]
    b=["d"]*(var_length-1);a=list(range(1,var_length));B=["D"]*(var_length-1);a=list(range(1,var_length))
    Delete_states=[m+str(n) for m,n in zip(b,a)];Delete_states_name=[m+str(n) for m,n in zip(B,a)];
    
    # Define HMM model
    model = HiddenMarkovModel( "Global Alignment")
    
    P=3*metadata_dict['length']+1
    hidden_states = [i for i in range(P)]
    for k in range(var_length): 
        # Create the insert states
        p_insert_char=dict(zip(metadata_dict['alphabet'],normalize_probability([convert_probability(p) for p in insertion_line[k]])))
        #p_insert_char=dict(zip(metadata_dict['alphabet'],[convert_probability(p) for p in insertion_line[k]]))
        hidden_states[k]= State( DiscreteDistribution(p_insert_char), name=Insert_states_name[k] )
        
    for k in range(var_length-1):    
        # Create the match states
        p_match_char=dict(zip(metadata_dict['alphabet'],normalize_probability([convert_probability(p) for p in match_line[1+k]])))
        #p_match_char=dict(zip(metadata_dict['alphabet'],[convert_probability(p) for p in match_line[1+k]]))
        hidden_states[(1+k+metadata_dict['length'])] = State( DiscreteDistribution(p_match_char) , name=Match_states_name[k] )
        
    for k in range(var_length-1):
        # Create the delete states
        hidden_states[(1+k+2*metadata_dict['length'])] = State( None, name=Delete_states_name[k] )
        
    # Add all the states to the model
    S=[i for i in hidden_states]
    model.add_states(S)
    
    # Create transitions for each position
    for k in range(metadata_dict['length']+1):
        temp=[convert_probability(p) for p in state_switch_line[k]]
        temp=normalize_probability(temp)
        if k==0:
            model.add_transition( model.start, hidden_states[(1+metadata_dict['length'])], temp[0] )#-m
            model.add_transition( model.start, hidden_states[0], temp[1] )#I
            model.add_transition( model.start, hidden_states[(1+2*metadata_dict['length'])], temp[2] )#-d
            model.add_transition( hidden_states[0], hidden_states[(1+metadata_dict['length'])], temp[3] )#i-m
            model.add_transition( hidden_states[0], hidden_states[0], temp[4] )#i-i
        elif k!=0 and k!=metadata_dict['length']:
            model.add_transition( hidden_states[(k+metadata_dict['length'])], hidden_states[(1+k+metadata_dict['length'])], temp[0] )#m-m
            model.add_transition( hidden_states[(k+metadata_dict['length'])], hidden_states[k], temp[1] )#m-i
            model.add_transition( hidden_states[(k+metadata_dict['length'])], hidden_states[(1+k+2*metadata_dict['length'])], temp[2] )#m-d
            model.add_transition( hidden_states[k], hidden_states[(1+k+metadata_dict['length'])], temp[3] )#i-m
            model.add_transition( hidden_states[k], hidden_states[k], temp[4] )#i-i
            model.add_transition( hidden_states[(k+2*metadata_dict['length'])], hidden_states[(1+k+metadata_dict['length'])], temp[5] )#d-m
            model.add_transition( hidden_states[(k+2*metadata_dict['length'])], hidden_states[(1+k+2*metadata_dict['length'])], temp[6] )#d-d            
        elif k==metadata_dict['length']:
            model.add_transition( hidden_states[2*metadata_dict['length']], model.end, temp[0] )#m
            model.add_transition( hidden_states[metadata_dict['length']], model.end, temp[3] )#i
            model.add_transition( hidden_states[3*metadata_dict['length']], model.end, temp[5] )#d
            model.add_transition( hidden_states[2*metadata_dict['length']], hidden_states[metadata_dict['length']], temp[1] )#m-i
            model.add_transition( hidden_states[metadata_dict['length']], hidden_states[metadata_dict['length']], temp[4])#i-i

    # Call bake to finalize the structure of the model.
    model.bake(merge="None")    
    return model





seqs_full="";
for seq_record in SeqIO.parse(input_fasta, "fasta"):
    seqs_full=seqs_full+str(seq_record.seq)#check sequences
    
if "X" in seqs_full:
    print("Error: your sequence contains X which is not general amino acid, please remove!")
    exit()
elif len(seqs_full)==0:
    print("Error: your sequence is empty!")
    exit()
else:
    PHMM_file_name=['DBLa0.1','DBLa0.10','DBLa0.11','DBLa0.12','DBLa0.13','DBLa0.14','DBLa0.15','DBLa0.16','DBLa0.17','DBLa0.18','DBLa0.19','DBLa0.2','DBLa0.20','DBLa0.21','DBLa0.22','DBLa0.23','DBLa0.24','DBLa0.3','DBLa0.4','DBLa0.5','DBLa0.6','DBLa0.7','DBLa0.8','DBLa0.9','DBLa1.1','DBLa1.2','DBLa1.3','DBLa1.4','DBLa1.5','DBLa1.6','DBLa1.7','DBLa1.8','DBLa2']
    header1=["seqID"];header=header1+PHMM_file_name[0:33]
    with open(output_PAD_csv, 'a+') as outfile:
        outfile.write(",".join([str(l) for l in header]) + "\n")
    with open(output_PBD_csv, 'a+') as outfile:
        outfile.write(",".join([str(l) for l in header]) + "\n")
    with open(output_PCD_csv, 'a+') as outfile:
        outfile.write(",".join([str(l) for l in header]) + "\n")    


    for seq_record in SeqIO.parse(input_fasta, "fasta"):
        zero_dict={};
        for i in PHMM_file_name[0:33]: zero_dict[i]=0
        log_prob_XAD = [];log_prob_XAD.append(str(seq_record.id))
        for PHMM_file in glob.glob("reference_data/upsA/*.hmm"):
            subclass=PHMM_file.split("upsA/")[1].split(".hmm")[0] 
            zero_dict[str(subclass)]=read_PHMM( str(PHMM_file) ).log_probability(list(str(seq_record.seq)));   
        log_prob_XAD=log_prob_XAD+list(zero_dict.values())
        with open(output_PAD_csv, 'a+') as outfile:
            outfile.write(",".join([str(l) for l in log_prob_XAD]) + "\n")
    
        zero_dict={};
        for i in PHMM_file_name[0:33]: zero_dict[i]=0
        log_prob_XBD = [];log_prob_XBD.append(str(seq_record.id))
        for PHMM_file in glob.glob("reference_data/upsB/*.hmm"):
            subclass=PHMM_file.split("upsB/")[1].split(".hmm")[0] 
            zero_dict[str(subclass)]=read_PHMM( str(PHMM_file) ).log_probability(list(str(seq_record.seq)));
        log_prob_XBD=log_prob_XBD+list(zero_dict.values())
        with open(output_PBD_csv, 'a+') as outfile:
            outfile.write(",".join([str(l) for l in log_prob_XBD]) + "\n")
    
        zero_dict={};
        for i in PHMM_file_name[0:33]: zero_dict[i]=0
        log_prob_XCD = [];log_prob_XCD.append(str(seq_record.id))
        for PHMM_file in glob.glob("reference_data/upsC/*.hmm"):
            subclass=PHMM_file.split("upsC/")[1].split(".hmm")[0] 
            zero_dict[str(subclass)]=read_PHMM( str(PHMM_file) ).log_probability(list(str(seq_record.seq)));
        log_prob_XCD=log_prob_XCD+list(zero_dict.values())
        with open(output_PCD_csv, 'a+') as outfile:
            outfile.write(",".join([str(l) for l in log_prob_XCD]) + "\n") 









