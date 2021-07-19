README.md

Dependencies:
- R/3.6.1
- bedtools/2.27.1

Execution:
1) Place .cov files from Bismark into "covs" folder
2) Within the covs folder, run dmr_pipeline_binarize.sh to binarize the .cov files
3) Within the covs folder, run dmr_pipeline_makeTSS.sh to make the RDS files with promoters info
4) Copy the two RDS file to the "GLM_TSS" folder
5) Run dmr_pipeline.sh to perform DMR between two groups 