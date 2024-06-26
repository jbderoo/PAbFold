![Logo of the project](images/Fig1.png)

PAbFold (Peptide-AntiBody AlphaFold2) is a pipeline centered around AlphaFold2 for predicting both how well an antigen:antibody pair might interact with one another, an where the epitope subsequence is in the antigen. I welcome you to read our paper:  https://biorxiv.org/cgi/content/short/2024.04.19.590298v1

Requirements: 
- Linux
- Python 3.9
  - numpy
- AlphaFold2 (with local colabfold recommended, https://github.com/YoshitakaMo/localcolabfold)
- Our localcolabfold patch (allows pickling of MSAs/template .cif's for local use, to both prevent requeuing of the MSA server that the folks at MMseqs (https://github.com/soedinglab/MMseqs2) so graciously have provided for general use, and speeds up inference time (we can use the same MSA/template over and over again without negatively impacting performance. See publication).

An antibody sequence should first be converted to a single chain variable fragment (scFv) to drastically reduce the computational requirements (therefore time) needed (__Fig A__). For assistance on this, please see the scFv-ification repo (coming soon!). The script `A_PeptideMaping_prep_submission_files.py` creates a fasta file of several sequences, where we chop up the antigen sequence (__Fig B__) into several subsequences such that we can test every individual epitope of a specificied length (Figure example: 10) and sliding window (Figure example: 1). An example run/wrapper script is provided as `make_files.sh`. Afterwards every sequence is then paired with the scFv sequence in a fasta file for folding. Here is a toy example of the output fasta file, that has the chopped up antigen sequence paired to the scFv sequence:

```
>batch_1_peptide_1
scFv_sequence:ABCDEFGHIJ

>batch_1_peptide_2
scFv_sequence:CDEFGHIJKL

>batch_1_peptide_3
scFv_sequence:EFGHIJKLMN
```

After all these structures are predicted, we can use the pLDDT (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3799472/) metric to assess how well the epitope binds into the CDR loops of the scFv. We can plot the sliding window pLDDT and per residue max pLDDT (__Fig C__) to identify regions in the antigen that may bind well with the provided scFv. The script `B_PeptideMaping_plddt_perres_analysis.py` does this analysis. An example run/wrapper is provided as `assess_files.sh`. After rank sorting the peptide sequences of length 10, the experimentally derived epitope sequence appears as the number one hit, and appears once more in the top 5 but shifted a residue (Example: ABCDEFGHIJ → BCDEFGHIJK). This analysis is also independent of residue positioning and does not consider whether the residues are buried or not - in the case of the HA trimer, the epitope is actually burried and is not easily accessible but is still detectable with PAbFold. This method out performs just giving AlphaFold2 the scFv sequence and antigen sequence and trying to let it blindly fold/dock both of them simultaneously, as this typically fails (__Fig D__).

After the structure files are created, .npy files are created for you to do whatever analysis you'd like with them. To create plots similar to __Fig C__, you must run `B_PeptideMaping_plddt_perres_analysis.py` with the `--all-models` flag. We also over the script `all_model_analysis.py` to process these files and create .svg's of the graphs. The graphs are technically created twice in a single svg; the left hand graphs in the svg file are the actual plots themselves, and the right hand one is used to get a large and viewable color bar for the sliding window plot. These objects are easily viewable and editable in any program that can build and edit svg's like Inkscape (https://inkscape.org/release/). For a more simple analysis that *only* looks at the top values ever seen (this is the bottom plot in __Fig C__), you do not need to include the `--all-models` flag, and can also use the `conf_plot_and_top10.py` python script for the simple analysis. A csv file is created in tandem that rank sorts the epitope sequences and their corresponding confidences (via pLDDT).
