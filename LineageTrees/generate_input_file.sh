#!/usr/bin/env bash
#$ -cwd
#$ -N JOB_NAME
#$ -pe smp 1
#$ -l h_rt=1:00:00
#$ -l h_rss=20G
#$ -m abe
#$ -M your_email
#$ -p -10

python build_csv_file.py /*.cov > my_sample_name.csv
wait
python build_phy_file.py my_sample_name.csv > my_sample_name.phy

