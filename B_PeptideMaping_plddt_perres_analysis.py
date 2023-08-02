import json
import glob
import numpy as np
import sys
import os
import argparse

parser = argparse.ArgumentParser()

# Add an argument
parser.add_argument('--antigen',  type=str,  required=True)
parser.add_argument('--data-dir',    type=str,  required=False, default='sliding_window_split')
parser.add_argument('--window-len',  type=int,  required=False, default=2)
parser.add_argument('--all-models',  action="store_true", required=False)#, default=False)


#binder_seq = 'MSDNGPQNQRNAPRITFGGPSDSTGSNQNGERSGARSKQRRPQGLPNNTASWFTALTQHGKEDLKFPRGQGVPINTNSSPDDQIGYYRRATRRIRGGDGKMKDLSPRWYFYYLGTGPEAGLPYGANKDGIIWVATEGALNTPKDHIGTRNPANNAAIVLQLPQGTTLPKGFYAEGSRGGSQASSRSSSRSRNSSRNSTPGSSRGTSPARMAGNGGDAALALLLLDRLNQLESKMSGKGQQQQGQTVTKKSAAEASKKPRQKRTATKAYNVTQAFGRRGPEQTQGNFGDQELIRQGTDYKHWPQIAQFAPSASAFFGMSRIGMEVTPSGTWLTYTGAIKLDDKDPNFKDQVILLNKHIDAYKTFPPTEPKKDKKKKADETQALPQRQKKQQTVTLLPAADLDDFSKQLQQSMSSADSTQA' # binder seq
#default_data_dir = 'sliding_window_split'
#window_len = 2

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
    
    #files = glob.glob(search_string)
    
    my_range = 1
    print('only requested top ranked model.')

else:
    search_string = '*/*_model_%s_*.json'
    my_range = 5
    print('requested all models.')
    


#rank = np.arange(1, num_ranks+1, 1)

for jd_counter in range(1, my_range+1):
    #jd_counter = 5    
    ssf = search_string % jd_counter
    
    files = glob.glob(ssf)
    #if len(files) == 0:
    #    files = glob.glob(f'*/*_rank_00{cds_rank}_*.json')
    files = sorted(files)

    if all_models:
        m = 'model'
    else:
        m = 'rank'
    print(f'working on {m} {jd_counter}')
    
    #print(files[:5])
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
    
    # 3D data of each batch, then each peptide in each batch,
    # then each plddt-aa in each peptide
    
    all_plddt = np.zeros(shape=[len(output_batches), batch_num, pep_len_fasta]) 
    all_seq   = np.zeros(shape=[len(output_batches), batch_num, pep_len_fasta], dtype=str)
    #sys.exit()
    
    for file in files:
    
        f = open(file)
        data = json.load(f)
        f.close()
        num = int(file.split('/')[0][-3:])
        # batch_005.fasta
        fasta = glob.glob(f'*{str(num).zfill(3)}.fasta')   
        peptide = file.split('/')[1].split('_')[:4]
        intra_fasta_batch   = int(peptide[1])
        batch = intra_fasta_batch-1
        pep_num = int(peptide[-1])
        peptide = '_'.join([str(x) for x in peptide])
        
        write_next_line = False
        written = False
    
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
    
        plddt = data['plddt'][-1*len(peptide_seq):]
     
        # debugging 
        # print(f'in batch {batch} and peptide {pep_num}, I have seq {peptide_seq} and scores: {plddt}')
        # print(f'peptide sequence__{peptide_seq}')
        # 3D data of each batch, then each peptide in each batch,
        # then each plddt-aa in each peptide
    
        assert len(peptide_seq) == len(plddt)
    
         
        for i in range(len(peptide_seq)):
     
            if int(all_plddt[batch, pep_num-1, i]) != 0:
                print('overwriting... uh oh!')
                sys.exit() 
            all_plddt[batch, pep_num-1, i]  = plddt[i]
            all_seq[batch, pep_num-1, i]    = peptide_seq[i]
    
    
    #np.save('plddt_matrix.npy', all_plddt)
    #np.save('seq_matrix.npy', all_seq)
    
    # union of 2 scripts 
    
    seq_list = []
    
    plddt_matrix = all_plddt
    seq_matrix   = all_seq
    
    batches, peps, n = plddt_matrix.shape
    
    plddt_list = np.zeros(shape=[len(binder_seq)+0, 150])
    seq_mat_list = np.zeros(shape=[len(binder_seq)+0, 150], dtype=str)
    master_pep_seq = ''
    cnt = 0
    my_end = False
    
    for i in range(batches):
        for j in range(peps):
            #for k in range(n):
    #        print(i)
    #        print(j)
            pep_seq_list = seq_matrix[i, j, :]
            plddt        = plddt_matrix[i, j, :]
            
            pep_seq = ''.join([str(x) for x in pep_seq_list]) 
            if pep_seq == '':
                continue       
    
            sub_seq = binder_seq[cnt:cnt+14]
    
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
             #       print(sub_seq)
                master_pep_seq += pep_seq[-2:]
           # print(idx)   
     
       
            if idx == -1:
    
                 if pep_seq == binder_seq[-1*len(pep_seq)]:
                     continue
    
                 else:
                     print('couldnt find pep_seq in binder_seq.... double check stuff? maybe weird fasta')
                     print(f'detected error with peptide sequence: {pep_seq}')
                     print(f'in sub sequence: {binder_seq[cnt:cnt+14]}')
                     print('This is most likely because you gave assess_files.sh a missmatched sequence from make_files.sh !')
                     sys.exit()
            #sys.exit() 
            idx += cnt
            #print(f'peptide {pep_seq} starts at {idx} of ^^^')
            if idx < (len(binder_seq) - len(pep_seq)):
                for k in range(len(pep_seq)):
                    findzeros = plddt_list[idx+k,:]
                    zeroidx   = np.where(findzeros == 0)[0]
                    plddt_list[idx+k, zeroidx[0]] = plddt[k]
                    seq_mat_list[idx+k, zeroidx[0]] = pep_seq[k]
            else:
                idx = len(binder_seq) - 1
                
                for k in range(len(pep_seq)):
                    findzeros = plddt_list[idx-k,:]
                    zeroidx   = np.where(findzeros == 0)[0]
                    plddt_list[idx-k, zeroidx[0]] = plddt[(len(plddt))-k-1]
                    seq_mat_list[idx-k, zeroidx[0]] = pep_seq[(len(pep_seq))-1-k]
    #        if j == 7 and i == 2:
    #            sys.exit() 
            cnt += window_len        
    #        if len(master_pep_seq) < 1:
    #            master_pep_seq = pep_seq
    #        else:
    #            master_pep_seq += pep_seq[-2:]
    
    #    sys.exit()
    master_pep_seq = ''.join([str(x) for x in seq_mat_list[:,0]])
    assert binder_seq == master_pep_seq
    avg_plddt = np.zeros(shape=len(binder_seq)) 
    
    
    for i in range(len(binder_seq)):
    
        idx = np.where(plddt_list[i, :] != 0)[0] 
        if max(idx) >= 149:
            print('expect an error: sliding window and peptide length are super wonky. change line 112 to something greater than 150')
        
        avg_plddt[i] = np.mean(plddt_list[i, idx])
    
    idx       = np.where(avg_plddt != 0)[0]
    avg_plddt = avg_plddt[idx]    
    
    if all_models == False:    
        np.save(f'../{jd_counter}_plddt_list.npy', plddt_list)
        max_plddt = np.max(plddt_list, axis=1)    
        np.save(f'../{jd_counter}_avg_plddt.npy', avg_plddt)
        np.save(f'../{jd_counter}_max_plddt.npy', max_plddt)
        
    else:
        
        ### start paste ###
        
        '''
        all_data = [plddt_list]#glob.glob('?_plddt_list.npy')

        for letters in all_data:
        l = letters[0]
        a = np.load(f'{l}_plddt_list.npy')
        '''
        a = plddt_list
        
        '''
        [[33.09  0.    0.    0.    0.  ]
         [45.69  0.    0.    0.    0.  ]
         [52.56 30.58  0.    0.    0.  ]
         [56.78 40.47  0.    0.    0.  ]
         [56.06 45.47 54.12  0.    0.  ]
         [63.75 47.   61.19  0.    0.  ]
         [64.94 47.62 63.69 37.72  0.  ]
         [63.72 46.44 64.5  44.59  0.  ]
         [65.31 46.09 65.31 48.81 32.41]
         [66.19 46.12 66.69 49.78 41.59]
         -------------------------------
         [46.88 68.56 50.88 48.12 27.73]
         [45.84 68.44 53.5  46.94 33.88]
         [66.94 57.31 52.34 35.84 32.12]
         [69.25 55.47 53.25 38.16 42.5 ]
         [54.53 50.44 37.53 48.28 25.73]
         [47.78 45.47 41.53 52.88 34.84]
         [40.72 36.75 52.5  40.09 33.72]
         [31.98 38.66 52.94 47.12 43.91]
         [36.97 52.38 47.22 48.81 35.69]
         [34.09 45.94 43.34 48.44 37.81]]
        
        first = [33.09 45.69 52.56 56.78 56.06 63.75 64.94 63.72 65.31 66.19]
        2nd   = [30.58 40.47 45.47 47.   47.62 46.44 46.09 46.12 46.88 45.84]
        
        d = np.array([[1,2,3], [4,5,6], [7,8,9]])
        np.pad(d, ((0,1), (0, 1)))
        '''
        
        z = np.pad(a, ((0,0), (0, 1000)))
        
        rows, cols = a.shape
        b = np.zeros(shape=[rows, 1000])
        
        start = -rows-2
        
        cnt = 0
        N = peptide_length_JD

        for i in range(0, rows-N, window_len):
           top_chunk = z[:i+N, :]
           bottom_chunk = z[i+N:]#, :-i+1]
           
           bcp = np.pad(bottom_chunk, ((0,0), (1, 0)))[:, :bottom_chunk.shape[1]]
           z = np.concatenate([top_chunk, bcp], axis=0)
           cnt += 1
           start += 2
           #print(f'ive run {cnt} times')
           
           #if i >= 4:
           #    print(z[:20, :10])
           #    sys.exit()
           
        idx = np.where(z[-1,:] != 0)[0][0]
        d = z[:, :idx+1]
        
        rows, cols = d.shape
                
        good_data = np.zeros(shape=[cols, peptide_length_JD])
        #np.save('good_data.npy', good_data)
        #print(f'good data shape: {good_data.shape}')


        #print(f'my pep length is {N}') 
        for i, pep in enumerate(d.transpose()):
            idx = np.where(pep != 0)[0]
            vals = pep[idx]
            #print(f'loop {i}')
            #print('values of {vals}')
            #np.save('good_data.npy', good_data)
            good_data[i, :] = vals
        datadir = os.getcwd().split('/')[-2]
        
        if not os.path.exists(f'../rowified_data'):
            os.mkdir(f'../rowified_data')

        if not os.path.exists(f'../rowified_data/{datadir}'):
            os.mkdir(f'../rowified_data/{datadir}')
        np.save(f'../rowified_data/{datadir}/model_{jd_counter}_sliding_window_plddt.npy', good_data)
        ### end paste   ###
        
        

os.chdir('..')

