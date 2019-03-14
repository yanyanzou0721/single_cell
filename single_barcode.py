# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 13:36:27 2019

@author: yyzou

Initial file:
     batcode.txt   reads.fq.gz

"""

import io
import os
from os.path import exists
import gzip
import click
from Bio import SeqIO
from edlib import align


def fastq_iter(file_in, phred):
    """ return a fastq iterator """
    if str(phred) == '33':
        fastq_iter = SeqIO.parse(file_in, 'fastq')
    else:
        fastq_iter = SeqIO.parse(file_in, 'fastq-illumina')
    return fastq_iter

def R2_dict(file_in,phred):
    if str(phred) == '33':
        fastq_iter = SeqIO.to_dict(SeqIO.parse(file_in, 'fastq'))
    else:
        fastq_iter = SeqIO.to_dict(SeqIO.parse(file_in, 'fastq-illumina'))
    return fastq_iter

def readin_barcode(file_barcode):
    """ read in known barcode for split """
    with open(file_barcode,'rt') as b:
        barcodes = b.readlines()
        b_list = []
        for i in barcodes:
            code = str(i).split('\n')[0]
            b_list.append(code)
    return b_list

def extract_barcode(id_str):
    """ extract barcode from fastq id string """
    barcode_str = id_str[0:8]
    return barcode_str

def write_out(file_name,cont,input_fastq,outdir):
    file = file_name + '_'+(input_fastq.split('_')[-1]).split('.')[0]
    out = outdir+'/'+file+'.fq'
    output = open(out,'a')
    SeqIO.write(cont,output,'fastq')
    output.close()

def s_max(barcode_sample,barcode_list,mismatch):
    score = {}
    for i in barcode_list:
        align_ = align(str(barcode_sample),i)
        score[i] = align_["editDistance"]
    result = sorted(score,key=lambda x:score[x])[0]
    if score[result] <= int(mismatch):
        return result
    else:
        return "not_matched"

@click.command(name="reads_depth_normalization")
@click.argument("input1")
@click.argument("input2")
@click.argument("input_barcode")
@click.argument("mismatch")
@click.option("--outdir", "-O",
    default="./",
    help="path to output splited fastq files.")
def main_(input1,input2,input_barcode,mismatch,outdir):
    if not exists(outdir):
        os.mkdir(outdir)
    b_list = readin_barcode(input_barcode)
    with io.TextIOWrapper(gzip.open(input2,'rb')) as f_R2:
        raw_R2 = R2_dict(f_R2,33)
    with io.TextIOWrapper(gzip.open(input1,'rb')) as f_R1:
        raw_R1 = fastq_iter(f_R1,33)
        for reads in raw_R1:
            b_extract = extract_barcode(reads.seq)
            b_class = s_max(b_extract,b_list,mismatch)
            write_out(b_class,reads,input1,outdir)
            write_out(b_class,raw_R2[reads.id],input2,outdir)

if __name__ == "__main__":
    main_()
