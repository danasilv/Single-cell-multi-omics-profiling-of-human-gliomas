README.md

Dependencies:
- Local installation of java (see website for instructions)
- see chromHMM manual for instructions (http://compbio.mit.edu/ChromHMM/ChromHMM_manual.pdf)

Execution:
bash learnHMM.sh

Required arguments:
-chromhmm : If this flag is set, Spectacle is run exactly the same as ChromHMM as described in the ChromHMM manual (see above).
-inputdir : This is the directory containing the binarized input files. Only file names containing ‘_binary’ are used by default (e.g., /BINARIZE/).
-outputdir : This is the directory where the output files are written (e.g., /HMM)
-numstates : This parameter specifies the number of states to include in the model (e.g., 18). 
-assembly : This parameter specifies the assembly. OverlapEnrichment and NeighborhoodEnrichment will be called
with default parameters using this assembly (e.g., hg38). 