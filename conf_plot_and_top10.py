#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 4 13:13:32 2022

@author: jderoo
"""

##################
### EDIT BELOW ###
##################

# data generated from epitope mapping pipeline
data = './myc_2e2_max_plddt.npy'


# known sequence for protein that binds to scFv/antibody
# NOTE: I break it up into multiple lines just to make it easier
# to read. It is not necessary to do this. You can just copy the
# whole sequence and paste it in between quotes without issue.
prot_seq = 'MDFFRVVENQPPATMPLNVSFTNRNYDLDYDSVQPYFYCDEEENFYQQQQQ' \
'SELQPPAPSEDIWKKFELLPTPPLSPSRRSGLCSPSYVAVTPFSLRGDNDGGGGSFSTADQL' \
'EMVTELLGGDMVNQSFICDPDDETFIKNIIIQDCMWSGFSAAAKLVSEKLASYQAARKDSGS' \
'PNPARGHSVCSTSSLYLQDLSAAASECIDPSVVFPYPLNDSSSPKSCASQDSSAFSPSSDSL' \
'LSSTESSPQGSPEPLVLHEETPPTTSSDSEEEQEDEEEIDVVSVEKRQAPGKRSESGSPSAG' \
'GHSKPPHSPLVLKRCHVSTHQHNYAAPPSTRKDYPAAKRVKLDSVRVLRQISNNRKCTSPRS' \
'SDTEENVKRRTHNVLERQRRNELKRSFFALRDQIPELENNEKAPKVVILKKATAYILSVQAE' \
'EQKLISEEDLLRKRREQLKHKLEQLRNSCA' 

# known sequence that does the actual binding event
epitope_seq = 'EQKLISEEDL'

# title of plot
title = 'mBG17 sliding window length 1'

# x axis spacing on the generated plot
tick_len = 20

# write top 10 peptides and their score to excel sheet?
peptide_to_excel = True

# length of peptides to investigate
pep_len = 10

# provide the path to a pdb structure to put the *data* plddts
# on that structure. It should be a structure of prot_seq. Set to
# False if you do not want do to this.
pdb = False


##################
### EDIT ABOVE ###
##################

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


conf = np.load(data)

peps = []
scores = []

for i in range(len(prot_seq)-pep_len):
    peps.append(prot_seq[i:i+pep_len])
    scores.append(round(np.mean(conf[i:i+pep_len]), 2))

scores = np.array(scores)
peps   = np.array(peps)
idx    = np.argsort(-1*scores)
top10  = idx[:10]

scores10 = scores[top10]
seqs10   = peps[top10]


df = pd.DataFrame(columns=['Pep Seq', 'Score'])

df['Pep Seq'] = seqs10
df['Score']   = scores10

name_list = data.split('_')
name = '_'.join([str(x) for x in name_list[:-1]])

if peptide_to_excel:
    df.to_excel(name + '.xlsx', index=False)

resis     = np.arange(1, len(conf)+1, 1)
subseq    = prot_seq.find(epitope_seq)
sub_conf  = conf[subseq:subseq+len(epitope_seq)]
sub_resis = resis[subseq:subseq+len(epitope_seq)]

if len(prot_seq) % tick_len > 0:
    adder = 2*tick_len
else:
    adder = tick_len
    
base = len(prot_seq) // tick_len

resi_tick = np.arange(1, base*tick_len+adder, tick_len)


plt.figure(dpi=600, figsize=[12,4])
plt.plot(resis, conf, 'k')
plt.plot(sub_resis, sub_conf, 'r')
plt.xticks(resi_tick, resi_tick-1)
plt.xlabel('residue number')
plt.ylabel('confidence')
plt.title(f'{title}')
plt.plot(resis, [60]*len(resis), 'm', linewidth=0.75)

pts = []
rs  = []
for i in range(len(conf)):
    if conf[i] > 60:
        pts.append(conf[i])
        rs.append(resis[i])

plt.plot(rs, pts, 'g.', markersize=7)
plt.savefig(f'{title}.svg', dpi=600, transparent=True)
#plt.show()

if pdb != False:
    
    f = open(pdb)
    w = open(data[:-4]+'.pdb', 'w') 
    
    #conf = np.load(data)
    
    start_resn = 0
    n_data = -1
    
    for line in f:
        if line[0:4] == 'ATOM' or line[0:4] == 'HETA':
            
            resn = int(line[23:27])
            if resn != start_resn:
                start_resn = resn
                n_data += 1
                
            b = conf[n_data]
            
            new_line = line[:60] + str(np.round(b, 2)).rjust(6) + line[66:]
            w.write(new_line)
            
        else:
            w.write(line)
            
            
    f.close()
    w.close()
