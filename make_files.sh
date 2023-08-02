SCFV_SEQ=MAEVQLVESGGDLVKPGGSLKLSCAASGFTFSHYGMSWVRQTPDKRLEWVATIGSRGTYTHYPDSVKGRFTISRDNAKNTLYLQMSSLKSEDTAMYYCARRSEFYYYGNTYYYSAMDYWGQGTSVTVSSGGGGSGGGGSGGGGSDIVLTQSPASLAVSLGQRATISCRASESVDNYGFSFMNWYQQKPGQPPKLLIYAISNRGSGIPARFSGSGSGTDFTLNIHPVEEEDAATYYCQQTKEVPWTFGGGTKLEI

BINDER_SEQ=MDFFRVVENQPPATMPLNVSFTNRNYDLDYDSVQPYFYCDEEENFYQQQQQSELQPPAPSEDIWKKFELLPTPPLSPSRRSGLCSPSYVAVTPFSLRGDNDGGGGSFSTADQLEMVTELLGGDMVNQSFICDPDDETFIKNIIIQDCMWSGFSAAAKLVSEKLASYQAARKDSGSPNPARGHSVCSTSSLYLQDLSAAASECIDPSVVFPYPLNDSSSPKSCASQDSSAFSPSSDSLLSSTESSPQGSPEPLVLHEETPPTTSSDSEEEQEDEEEIDVVSVEKRQAPGKRSESGSPSAGGHSKPPHSPLVLKRCHVSTHQHNYAAPPSTRKDYPAAKRVKLDSVRVLRQISNNRKCTSPRSSDTEENVKRRTHNVLERQRRNELKRSFFALRDQIPELENNEKAPKVVILKKATAYILSVQAEEQKLISEEDLLRKRREQLKHKLEQLRNSCA

CACHE=myc-2E2.pkl

NUM_GPU=1


PEP_LEN=9
SL=2
SPF=115

python A_PeptideMaping_prep_submission_files.py --dock-against $SCFV_SEQ --dock-with $BINDER_SEQ --cache $CACHE --peptide-len $PEP_LEN --num-gpus $NUM_GPU --sliding-window $SL --seqs-per-file $SPF