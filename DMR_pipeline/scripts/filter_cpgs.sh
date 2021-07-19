for i in *cov;
do
count=$(wc -l $i | awk '{print $1}')
if [ $count -lt 50000 ];
then
echo $count
mv $i ./failed_QC
fi
done

