export LC_ALL=de_DE.utf-8
export LANG=de_DE.utf-8

ls /public/home/yyzou/sci-293/Cleandata/*/*R2.fq.gz |\
awk -F "/" '{print $7}' |\
while read i;\
do
echo "python /public/home/yyzou/sci-293/single_barcode.py /public/home/yyzou/sci-293/Cleandata/$i/${i}_R1.fq.gz /public/home/yyzou/sci-293/Cleandata/$i/${i}_R2.fq.gz /public/home/yyzou/sci-293/barcode 3 -O /public/home/yyzou/sci-293/Cleandata/$i/cell;
gzip /public/home/yyzou/sci-293/Cleandata/$i/cell/* ;"  >/public/home/yyzou/sci-293/Cleandata/$i/${i}_split.sh;
qsub <<EOF
#PBS -N $i
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o /public/home/yyzou/sci-293/log/$i.pbs
#PBS -V
bash /public/home/yyzou/sci-293/Cleandata/$i/${i}_split.sh
EOF
done

