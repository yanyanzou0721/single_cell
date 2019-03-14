#!/bin/bash

set -e

echo "You are running scritp : $0"
echo "input_data_dir : $1"
echo "ref : $2"
echo "gtf : $3"
echo "align_out_dir : $4"
echo "expression_out_dir : $5"

raw=$1;
genome=$2;
gtf=$3;
sc_dir=$4;
exp=$5;

if [ ! -d "${sc_dir}/" ];then
   mkdir ${sc_dir}/
fi

if [ ! -d "${exp}/" ];then
   mkdir ${exp}/
fi

if [ ! -d "${sc_dir}/map/" ];then
   mkdir ${sc_dir}/map/
fi

if [ ! -d "${sc_dir}/align_stat/" ];then
mkdir ${sc_dir}/align_stat/
fi

if [ ! -d "${sc_dir}/script/" ];then
mkdir ${sc_dir}/script/
fi

if [ ! -d "${sc_dir}/hisat_summary/" ];then
mkdir ${sc_dir}/hisat_summary/
fi

if [ ! -d "${sc_dir}/bam/" ];then
mkdir ${sc_dir}/bam/
fi

if [ ! -d "${exp}/raw_featureCount/" ];then
mkdir ${exp}/raw_featureCount/
fi

if [ ! -d "${exp}/exp/" ];then
mkdir ${exp}/exp/
fi

## input: path/group/cell , group_cell determine a single cell


ls ${raw}/*_R2.fq.gz|\
while read i;
do
sample=`echo ${raw} | awk -F "/" '{print $(NF-2)}'`;
cell=`echo $i | awk -F "/" '{print $NF'} | awk -F "_" '{print $1}'`;

zcat ${raw}/${cell}_R2.fq.gz | wc -l >>${sc_dir}/align_stat/reads_num &&

hisat2 -q -3 100  -x ${genome}  -U ${raw}/${cell}_R2.fq.gz -S ${sc_dir}/map/${sample}_${cell}.sam --summary-file ${sc_dir}/hisat_summary/${sample}_${cell}.summary &&

samtools view -bS ${sc_dir}/map/${sample}_${cell}.sam | samtools sort -o ${sc_dir}/bam/${sample}_${cell}.sorted.bam &&

samtools index ${sc_dir}/bam/${sample}_${cell}.sorted.bam &&

samtools flagstat ${sc_dir}/bam/${sample}_${cell}.sorted.bam >>${sc_dir}/align_stat/${sample}_${cell}.align &&

echo "${sample}_${cell}">> ${sc_dir}/align_stat/cell &&

cat ${sc_dir}/align_stat/${sample}_${cell}.align | grep 'mapped (' | awk -F "(" '{print $2}' | awk -F ":" '{print $1}'|cut -d "%" -f 1 >>${sc_dir}/align_stat/map_rate &&

featureCounts -a $gtf ${sc_dir}/bam/${sample}_${cell}.sorted.bam -t exon -o ${exp}/raw_featureCount/${sample}_${cell}_count.txt.full -T 8 -g gene_id -p -s 2 &&

cat ${exp}/raw_featureCount/${sample}_${cell}_count.txt.full | sed 1,2d | awk 'BEGIN{OFS="\t"}{print $1,$7}' | sort > ${exp}/exp/${sample}_${cell}.count.txt
done

bash /public/home/yyzou/sci-293/merge.sh ${exp}/exp/*.count.txt >${exp}/exp.txt &&

paste ${sc_dir}/align_stat/cell ${sc_dir}/align_stat/reads_num ${sc_dir}/align_stat/map_rate > ${sc_dir}/reads_align_stat.txt;

python /public/home/yyzou/sci-293/sc_plot.py -i ${sc_dir}/reads_align_stat.txt -e ${exp}/exp.txt -o ${exp}/test.pdf




