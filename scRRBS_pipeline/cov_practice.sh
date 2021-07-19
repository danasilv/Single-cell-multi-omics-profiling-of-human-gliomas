cov_file=${1}
sample_name=${2}
PE_report=${3}
PDR_file=${4}
Description=${5}

python cov_to_fullstats.py $cov_file $sample_name $PE_report $PDR_file SingleCell $Description

