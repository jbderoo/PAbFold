import json
import glob
import numpy as np
import sys
import os
import argparse
import shutil
import re
from collections import defaultdict

parser = argparse.ArgumentParser()

# Add an argument
parser.add_argument('--antigen',     type=str,  required=True)
parser.add_argument('--data-dir',    type=str,  required=False, default='sliding_window_split')
parser.add_argument('--window-len',  type=int,  required=False, default=2)
parser.add_argument('--all-models',  action="store_true", required=False)#, default=False)

args             = parser.parse_args()
default_data_dir = args.data_dir
binder_seq       = args.antigen.upper()
window_len       = args.window_len
all_models       = args.all_models

if not os.path.exists(default_data_dir):
    print('did you change the name of the output dir to something other than sliding_window_split? Cant find my files!')
    sys.exit()

os.chdir(default_data_dir)

if all_models == False:
    search_string = '*/*_rank_%s_*.json'
    
    if len(glob.glob(search_string % 1)) == 0:
        search_string = '*/*_rank_00%s_*.json'
    
    my_range = 1
    print('only requested top ranked model.')

else:
    search_string = '*/*_model_%s_*.json'
    my_range = 5
    print('requested all models.')
    


# CDS window function
def window(sequence, winSize, step=1):
    """Returns a generator that will iterate through
    the defined chunks of input sequence.  Input sequence
    must be iterable."""

    # Verify the inputs
    if not isinstance(sequence, (list, str, tuple)):
        raise Exception("**ERROR** sequence must be iterable.")
    if not (isinstance(winSize, int) and isinstance(step, int)):
        raise Exception("**ERROR** winSize and step must be of type int.")
    if step > winSize:
        raise Exception("**ERROR** step must not be larger than winSize.")
    if winSize > len(sequence):
        raise Exception("**ERROR** winSize must not be larger than sequence length.")

    for i in range(0, len(sequence) - winSize + 1, step):
        yield sequence[i:i+winSize]




##### Check files from ColabFold to help identify if EVERYTHING has run #####

def how_many_files(seq, slide, ss):
    detected    = len(glob.glob(ss % 1))

    # find the peptide length
    output_batches = sorted(glob.glob('*.fasta'))
    opener = sorted([x for x in output_batches if '01.' in x])
    f = open(opener[0])
    pep_len = 0
    for line in f:
        if '>' not in line:
            data = line.split(':')[-1]
            pep_len = len(data) - 1 # rude \n
        if pep_len > 0:
            break

    f.close()
   
    theoretical = int(np.ceil( (len(seq) - pep_len) / slide) + 1)
    return detected, theoretical, pep_len

def check_files(directory, pattern):
    file_list = os.listdir(directory)
    batch_peptide_counts = defaultdict(int)
    missing_files = []

    for filename in file_list:
        match = pattern.match(filename)
        if match:
            batch, peptide = map(int, match.groups())
            batch_peptide_counts[batch] = max(batch_peptide_counts[batch], peptide)

    for batch in sorted(batch_peptide_counts):
        for peptide in range(1, batch_peptide_counts[batch]):
            expected_file = f'batch_{batch}_peptide_{peptide}_unrelaxed_rank_001_*.pdb'
            if not any(pattern.match(f) and int(pattern.match(f).group(2)) == peptide for f in file_list):
                missing_files.append(expected_file)

    return missing_files

def extract_shortest_chain_sequence(file_path, amino_acid_map):
    try:
        with open(file_path, 'r') as file:
            chain_residues = defaultdict(set)
            for line in file:
                if line.startswith('ATOM'):
                    chain_id, residue_number = line[21], line[22:26].strip()
                    chain_residues[chain_id].add(residue_number)

            shortest_chain = min(chain_residues, key=lambda k: len(chain_residues[k]))
            sequence, residue_sequence_numbers = '', set()
            file.seek(0)
            for line in file:
                if line.startswith('ATOM') and line[21] == shortest_chain:
                    residue_name, residue_number = line[17:20].strip(), line[22:26].strip()
                    if residue_name in amino_acid_map and residue_number not in residue_sequence_numbers:
                        sequence += amino_acid_map[residue_name]
                        residue_sequence_numbers.add(residue_number)

            return sequence
    except IOError as e:
        print(f"Error reading file {file_path}: {e}")
        return ''

def build_sequence(directory, pattern, amino_acid_map):
    file_list = sorted(os.listdir(directory), key=lambda f: (int(pattern.match(f).group(1)), int(pattern.match(f).group(2))) if pattern.match(f) else (0, 0))
    sequences = [extract_shortest_chain_sequence(os.path.join(directory, f), amino_acid_map) for f in file_list if pattern.match(f)]

    if not sequences:
        print(f"No sequences extracted from any files in {directory}")
        return ''

    def overlap_sequences(seq1, seq2):
        overlap = min(len(seq1), len(seq2))
        while overlap > 0 and seq1[-overlap:] != seq2[:overlap]:
            overlap -= 1
        return seq1 + seq2[overlap:]

    continuous_sequence = sequences[0]
    for seq in sequences[1:]:
        continuous_sequence = overlap_sequences(continuous_sequence, seq)

    return continuous_sequence

def process_all_batches(parent_directory):
    dir_pattern = re.compile(r'output_batch_(\d{3})')
    pdb_pattern = re.compile(r'batch_(\d+)_peptide_(\d+)_unrelaxed_rank_001_.*\.pdb')
    amino_acid_map = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

    batch_directories = [d for d in os.listdir(parent_directory) if os.path.isdir(os.path.join(parent_directory, d)) and dir_pattern.match(d)]
    batch_directories.sort(key=lambda x: int(dir_pattern.match(x).group(1)))

    if not batch_directories:
        print(f"No 'output_batch_XXX' directories found in {parent_directory}")
        return ''

    full_sequence = ''
    for batch_dir in batch_directories:
        directory_path = os.path.join(parent_directory, batch_dir)
        batch_sequence = build_sequence(directory_path, pdb_pattern, amino_acid_map)
        full_sequence += batch_sequence if not full_sequence else batch_sequence[len(full_sequence):]

    return full_sequence



# Checkpoint 1:
found, should, epitope_len = how_many_files(binder_seq, window_len, search_string)
if found != should:
    print('### EXPECT ERROR ###')
    print(f'With a sequence length of {len(binder_seq)}, window length of {window_len}, and an automatically detected epitope sequence length of {epitope_len}, there should be {should} structure prediction calculations but we found {found}. Attempting to ID missing files with more robust identification...')
else:
    print('Checkpoint 1 passed!')

# Checkpoint 2:
pdb_pattern = re.compile(r'batch_(\d+)_peptide_(\d+)_unrelaxed_rank_001_.*\.pdb')
missing = check_files('./', pdb_pattern)
if missing:
    print("Missing files:")

    for file in missing:
        print(file)
    print('excpect errors!')
else:
    continuous_sequence = process_all_batches('./')
    print("Continuous sequence built! Checkpoint 2 succeeded. Should be smooth sailing")

#############################################################################










# begin analysis:

for jd_counter in range(1, my_range+1):
    ssf = search_string % jd_counter
    
    files = glob.glob(ssf)
    files = sorted(files)

    if all_models:
        m = 'model'
    else:
        m = 'rank'
    print(f'working on {m} {jd_counter}')
    
    output_batches = sorted(glob.glob('*.fasta'))
    
    batch_num = 0
    pep_len_fasta = 0
    opener = sorted([x for x in output_batches if '01.' in x])
    f = open(opener[0])

    for line in f:
        if '>' in line:
            batch_num += 1
    
        else:
            if pep_len_fasta == 0:
                data = line.split(':')[-1]
                pep_len_fasta = len(data) - 1 # rude \n
    f.close()
    
    all_plddt = np.zeros(shape=[len(output_batches), batch_num, pep_len_fasta]) 
    all_seq   = np.zeros(shape=[len(output_batches), batch_num, pep_len_fasta], dtype=str)
    
    for file in files:
    
        with open(file) as f:
            data = json.load(f)

        num               = int(file.split('/')[0][-3:])
        fasta             = glob.glob(f'*{str(num).zfill(3)}.fasta')   
        peptide           = file.split('/')[1].split('_')[:4]
        intra_fasta_batch = int(peptide[1])
        batch             = intra_fasta_batch-1
        pep_num           = int(peptide[-1])
        peptide           = '_'.join([str(x) for x in peptide])
        write_next_line   = False
        written           = False
    
        with open(fasta[0]) as f1:
            for line in f1:
                if write_next_line == True:
                    seq = line
                    write_next_line = False
    
                if f'>{peptide}\n' == line:
                    write_next_line = True
    
                if written == True:
                    break
    
        peptide_seq = seq.split(':')[-1][:-1] #rude \n
        plddt       = data['plddt'][-1*len(peptide_seq):]
     
        if len(peptide_seq) != len(plddt):
            print('WARNING! missmatch in plddt size and peptide_seq length. Do you have all your files? Did you specify good inputs?')
            
        for i in range(len(peptide_seq)):
            if  all_plddt[batch, pep_num-1, i] <= plddt[i]:
                all_plddt[batch, pep_num-1, i]  = plddt[i]
                all_seq  [batch, pep_num-1, i]  = peptide_seq[i]
    
    seq_list         = []
    plddt_matrix     = all_plddt
    seq_matrix       = all_seq
    batches, peps, n = plddt_matrix.shape
    plddt_list       = np.zeros(shape=[len(binder_seq)+0, 150])
    seq_mat_list     = np.zeros(shape=[len(binder_seq)+0, 150], dtype=str)
    master_pep_seq   = ''
    cnt              = 0
    my_end           = False
    
    for i in range(batches):
        for j in range(peps):
            pep_seq_list = seq_matrix[i, j, :]
            plddt        = plddt_matrix[i, j, :]
            
            pep_seq = ''.join([str(x) for x in pep_seq_list]) 
            if pep_seq == '':
                continue       
    
            sub_seq = binder_seq[cnt:cnt+(2*pep_len_fasta)]
    
            if i == 0 and j == 0:
                idx     = sub_seq.find(pep_seq)
                master_pep_seq += pep_seq
                peptide_length_JD = len(pep_seq) 
     
            else:
                if len(sub_seq) < len(pep_seq):
                    sub_seq = binder_seq[-1*len(pep_seq):]
                    idx = sub_seq.find(pep_seq)
                    my_end == True

                else:
                    idx     = sub_seq.find(pep_seq)
                
                master_pep_seq += pep_seq[-2:]
     
            
            if idx == -1:
    
                 if pep_seq == binder_seq[-1*len(pep_seq)]:
                     continue
    
                 else:
                     print('couldnt find pep_seq in binder_seq.... double check stuff? maybe weird fasta')
                     print(f'detected error with peptide sequence: {pep_seq}')
                     print(f'in sub sequence: {binder_seq[cnt:cnt+(2*pep_len_fasta)]}, around line 300')
                     print('This is most likely because you gave assess_files.sh a missmatched sequence from make_files.sh !')
                     sys.exit()
            idx += cnt

            if idx < (len(binder_seq) - len(pep_seq)):
                for k in range(len(pep_seq)):
                    findzeros = plddt_list[idx+k,:]
                    zeroidx   = np.where(findzeros == 0)[0]
                    plddt_list[idx+k, zeroidx[0]]   = plddt[k]
                    seq_mat_list[idx+k, zeroidx[0]] = pep_seq[k]

            else:
                idx = len(binder_seq) - 1
                
                for k in range(len(pep_seq)):
                    findzeros = plddt_list[idx-k,:]
                    zeroidx   = np.where(findzeros == 0)[0]
                    plddt_list[idx-k, zeroidx[0]]   = plddt[(len(plddt))-k-1]
                    seq_mat_list[idx-k, zeroidx[0]] = pep_seq[(len(pep_seq))-1-k]
            cnt += window_len        
    
    master_pep_seq = ''.join([str(x) for x in seq_mat_list[:,0]])
    if binder_seq != master_pep_seq:
        print('after parsing through the structure prediction files, I could not correctly rebuild the antigen sequence. Double check input files and input variables - potentially rerun AF2 predictions!')
    avg_plddt = np.zeros(shape=len(binder_seq)) 
    
    for i in range(len(binder_seq)):
        idx          = np.where(plddt_list[i, :] != 0)[0] 
        avg_plddt[i] = np.mean(plddt_list[i, idx])
    
    idx       = np.where(avg_plddt != 0)[0]
    avg_plddt = avg_plddt[idx]    


    # make seqs.py for easy CDS analysis merge
    seqs_file = open('seqs.py', 'w')
    PL        = len(peptide_seq)
    slide     = window_len 
    seqlist   = [''.join(x) for x in window(binder_seq,PL,slide)]    
    seqs_file.write("seqs = " + str(seqlist) + '\n')
    seqs_file.write(f'numaa = {int(len(binder_seq))}')
    seqs_file.write(f'\nslide = {int(slide)}')
    seqs_file.write(f'\npeplen = {int(PL)}')
    seqs_file.close()
    shutil.move('./seqs.py', '../seqs.py')

    
    if all_models == False:    
        max_plddt = np.max(plddt_list, axis=1)    
        np.save(f'../{jd_counter}_avg_plddt.npy', avg_plddt)
        np.save(f'../{jd_counter}_max_plddt.npy', max_plddt)

        a = plddt_list        
        z = np.pad(a, ((0,0), (0, 1000)))

        rows, cols = a.shape
        b = np.zeros(shape=[rows, 1000])

        start = -rows-2

        cnt = 0
        N   = peptide_length_JD
        # truly gross logic to rebuild the all_models sliding window 3D plot. Basically, bump out every row (pad with zeros to the left) to rebuild the sliding window scheme for easy analysis.
        for i in range(0, rows-N, window_len):
           top_chunk    = z[:i+N, :]
           bottom_chunk = z[i+N:]#, :-i+1]

           bcp    = np.pad(bottom_chunk, ((0,0), (1, 0)))[:, :bottom_chunk.shape[1]]
           z      = np.concatenate([top_chunk, bcp], axis=0)
           cnt   += 1
           start += 2

        idx        = np.where(z[-1,:] != 0)[0][0]
        d          = z[:, :idx+1]
        rows, cols = d.shape

        good_data = np.zeros(shape=[cols, peptide_length_JD])

        for i, pep in enumerate(d.transpose()):
            idx = np.where(pep != 0)[0]
            vals = pep[idx]
            good_data[i, :] = vals
        np.save('../max_sliding_window.npy', good_data)


    else:
        
        ### start paste ###
        
        a = plddt_list
        z = np.pad(a, ((0,0), (0, 1000)))
        
        rows, cols = a.shape
        b = np.zeros(shape=[rows, 1000])
        
        start = -rows-2
        
        cnt = 0
        N = peptide_length_JD

        for i in range(0, rows-N, window_len):
           top_chunk = z[:i+N, :]
           bottom_chunk = z[i+N:]#, :-i+1]
           
           bcp    = np.pad(bottom_chunk, ((0,0), (1, 0)))[:, :bottom_chunk.shape[1]]
           z      = np.concatenate([top_chunk, bcp], axis=0)
           cnt   += 1
           start += 2
           
        idx        = np.where(z[-1,:] != 0)[0][0]
        d          = z[:, :idx+1]
        rows, cols = d.shape

                
        good_data = np.zeros(shape=[cols, peptide_length_JD])
        for i, pep in enumerate(d.transpose()):
            idx  = np.where(pep != 0)[0]
            vals = pep[idx]
            good_data[i, :] = vals
        datadir  = os.getcwd().split('/')[-2]
        
        if not os.path.exists(f'../rowified_data'):
            os.mkdir(f'../rowified_data')

        np.save(f'../rowified_data/model_{jd_counter}_sliding_window_plddt.npy', good_data)
        
os.chdir('..')

