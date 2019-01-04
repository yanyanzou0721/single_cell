# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 10:17:27 2018
rewrite https://github.com/Nanguage/double_index_spliter
@author: Administrator
"""

import re
from collections import namedtuple
import logging
import os
from os.path import join, exists
import gzip
import sys
from sys import stderr
from io import TextIOWrapper

import gzip
import click
from Bio import SeqIO

log = logging.getLogger(__name__)

class barcode(object):
    """ Denote a kind of index(a,b) combination """
    def __init__(self, name, phred=33):
        self.name = name
        self.phred=phred

    def fopen(self, path="./", suffix=".fq"):
        self.file = open_(join(path, self.name+suffix), 'w')
        self.writer = fastq_writer(self.file, self.phred)
        self.writer.write_header()


def open_(fname, mode):
    """ open file according to it's suffix """
    assert 'b' not in mode
    if fname.endswith(".gz"):
        if sys.version_info > (3, 0, 0):
            return TextIOWrapper(gzip.open(fname, mode+'b'))
        else:
            return gzip.open(fname, mode+'b')
    else:
        return open(fname, mode)
    
def fastq_iter(file_in, phred):
    """ return a fastq iterator """
    if str(phred) == '33':
        fastq_iter = SeqIO.parse(file_in, 'fastq')
    else:
        fastq_iter = SeqIO.parse(file_in, 'fastq-illumina')
    return fastq_iter

def fastq_writer(file_out, phred):
    """ return a fastq writer """
    if str(phred) == '33':
        writer = SeqIO.QualityIO.FastqPhredWriter(file_out)
    else:
        writer = SeqIO.QualityIO.FastqIlluminaWriter(file_out)
    return writer

def extract_barcode(id_str):
    """ extract barcode from fastq id string """
    barcode_str = (id_str.split("=")[1]).split(" ")[0]
    return barcode_str
    

def process_all(input_fq, phred):
    with open_(input_fq, 'r') as f:
        fq_iter = fastq_iter(f, phred)
        for rec in fq_iter:
            barcode = extract_barcode(rec.description)
            barcode.writer.write_record(rec)
                
            
@click.command(name="barcode_spliter")
@click.argument("input_fastq")
@click.option("--outdir", "-O",
    default="./",
    help="path to output splited fastq files.")
@click.option("--gzip/--no-gzip", "-z",
    default=False,
    help="compress output fastq files with gzip.")
@click.option("--phred",
    default=33,
    help="encode of fastq quality string.")

def main_(input_fastq, outdir, gzip, phred):

    if not exists(outdir):
        os.mkdir(outdir)
    # touch all sub fq file
    suffix = ".fq" if not gzip else ".fq.gz"
    
    with open_(input_fastq, 'r') as f:
        fq_iter = fastq_iter(f, phred)
        for rec in fq_iter:
            barcode = extract_barcode(rec)
            barcode.fopen(path=outdir, suffix=suffix)
            msg = "splited baocode: {} 's fastq will store to {}".format(
                str(barcode), join(outdir, barcode.name)+suffix )
        log.info(msg)
    
    process_all(input_fastq, phred)


if __name__ == "__main__":
    log.setLevel(logging.DEBUG)
    hdr = logging.StreamHandler(stream=stderr)
    fmt = logging.Formatter(
        fmt="[%(levelname)s][%(asctime)s] %(message)s",
        datefmt="%m/%d/%y %H:%M:%S")
    hdr.setFormatter(fmt)
    log.addHandler(hdr)
    main_()                