# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 18:59:28 2022

@author: jderoo
"""

import argparse
import sys

# Create the parser
parser = argparse.ArgumentParser()

# Add an argument
parser.add_argument('--scfv',           type=str,  required=True)
parser.add_argument('--antigen',        type=str,  required=True)
parser.add_argument('--peptide-len',    type=int,  required=False, default=10)
parser.add_argument('--sliding-window', type=int,  required=False, default=2)
parser.add_argument('--seqs-per-file',  type=int,  required=False, default=10)
parser.add_argument('--demo',           type=bool, required=False, default=False)
parser.add_argument('--cache',          type=str,  required=False)
parser.add_argument('--num-gpus',       type=int,  required=False, default=1)

args = parser.parse_args()


##################
### EDIT BELOW ###
##################

scFv_seq      = args.scfv
prot_to_chunk = args.antigen
pep_len       = args.peptide_len
window_len    = args.sliding_window
seqs_per_file = args.seqs_per_file
demo          = args.demo
cache         = args.cache
num_gpus      = args.num_gpus
#pep_len       = 10   # length of the peptide
#window_len    = 2    # how far do we slide to generate a new peptide?
#seqs_per_file = 10   # sequences to write to a single file for AF2
#demo          = True # is this a demonstration or is this a production run?


##################
### EDIT ABOVE ###
##################


# write a function to keep things pretty
def write_fasta(all_peptides, N):
       
    f = open(f'batch_{str(N+1).zfill(3)}.fasta', 'w')
    
    for i, pep in enumerate(all_peptides):
        f.write(f'>batch_{N+1}_peptide_{i+1}\n')
        f.write(f'{scFv_seq.upper()}:{pep}\n\n')
        
    f.close()
    
        
    


import numpy as np
import os

# make a new folder to keep things pretty

if not os.path.exists('sliding_window_split'):
    os.mkdir('sliding_window_split')

home = os.getcwd()    
os.chdir('sliding_window_split')

seq = prot_to_chunk.upper() 
all_peptides = []
num_files    = 0

for i in range(len(seq)):
    
    n = i*window_len
    peptide = seq[n:n+pep_len]
    
    if demo == True:  # demonstration purposes: print the first 3 peptide seq's
        print(peptide)
        
        if i == 2:  
            break
    
    else: # production run
    
        all_peptides.append(peptide)
        
        
        # if we have run out of peptide chunks/slices to make, write current peptides
        if len(peptide) < pep_len:
            
            all_peptides = all_peptides[:-1]
            all_peptides.append(seq[-1*pep_len:])
            
            # write to file w/ current peptides
            write_fasta(all_peptides, num_files)

            if (num_files+1) % 2 == 1:
                device_num = 0
            else:
                if num_gpus > 1:
                    device_num = 1
                else:
                    device_num = 0

            if cache:
                cache_line = f' --use-cached-mmseq-results ../{cache}'
            else:
                cache_line = ''



                # write a 'runit' for every chunk because I am lazy
            f = open(f'run_batch_{num_files+1}.sh', 'w')
            f.write (f'export CUDA_VISIBLE_DEVICES={device_num}\n' \
                     f'colabfold_batch --templates --num-recycle 6' \
                     f'{cache_line}' \
                     f' batch_{str(num_files+1).zfill(3)}.fasta' \
                     f' output_batch_{str(num_files+1).zfill(3)}')
            f.close ()
                
 
            break
        
        else: # still under seqs_per_file peptide per file
            
            # if we reach seqs_per_file peptides in memory, write to file for AF2
            
            if len(all_peptides) == seqs_per_file:
                write_fasta(all_peptides, num_files)
                
                if (num_files+1) % 2 == 1:
                    device_num = 0
                else:
                    if num_gpus > 1:
                        device_num = 1
                    else:
                        device_num = 0

                if cache:
                    cache_line = f' --use-cached-mmseq-results ../{cache}'
                else:
                    cache_line = ''

                # write a 'runit' for every chunk because I am lazy
                f = open(f'run_batch_{num_files+1}.sh', 'w')
                f.write (f'export CUDA_VISIBLE_DEVICES={device_num}\n' \
                         f'colabfold_batch --templates --num-recycle 6' \
                         f'{cache_line}' \
                         f' batch_{str(num_files+1).zfill(3)}.fasta' \
                         f' output_batch_{str(num_files+1).zfill(3)}')
                f.close ()
                
                
                all_peptides = [] # reset peptides in memory
                num_files   += 1  # add 1 to number of files written to disk
                

# come back to primary directory when all done because some things hate not coming back                
os.chdir(home)                
                
