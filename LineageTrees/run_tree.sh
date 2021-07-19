#!/usr/bin/env bash
#$ -cwd
#$ -N JOB_NAME
#$ -pe smp 4 
#$ -l h_rt=504:00:00
#$ -l h_rss=30G
#$ -m abe
#$ -j y 
#$ -R y
#$ -M your_email
#$ -p -10

/iqtree-1.6.9-Linux/bin/iqtree -s ./my_sample_name.phy -st BIN -m TESTNEW -mset GTR2 -nt AUTO -bb 1000 -wbtl -abayes -alrt 1000 -nstop 500 -wsr -nm 2000 -v -alninfo
