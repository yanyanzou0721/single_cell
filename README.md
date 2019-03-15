# **sci-293 single_cell data analysis**
## cell num 4896 of which 96 are not recognized

> **pre analysis steps**
```
cell split              ![method：single_barcode.py](https://github.com/yanyanzou0721/single_cell_split)
trim
hisat2                   align
samtools
featureCounts            expression statistic
```
----------
# example
## *Step one ： cell split for single barcode*
```
python single_barcode.py input_R1.fastq.gz input_R2.fastq.gz input.barcode mismatch_num -O outdir
```
> **input**
    input_R1.fastq.gz
    input_R2.fastq.gz
    input_barcode
    mismatch_num
    outdir

## *Step two ：trim*

## *Step three ：align & gene count & plot*
```
bash scrna-pipe.sh  data_dir  ref  gtf  align_ourdir  exp_outdir
```
      This .sh  script will automaticly calling the python script *sc_plot.py*
for sc_plot.py
```
python sc_plot.py -i reads_align_stat.txt -e exp.txt -o test.pdf
```
Result pre-display
![avatar](example/test.pdf）
