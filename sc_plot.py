#!/usr/bin python3

import pandas as pd
import numpy as np
import matplotlib.pyplot  as plt
import sys,getopt
from matplotlib.backends.backend_pdf import PdfPages

def count(series):
    count=0
    for i in series:
        if int(i) >0:
            count=count+1
    return count

def process(reads_align_stat,expression,outputfile):
    # cell reads & alignment data prepare
    cell_reads = pd.read_table(reads_align_stat,header=None,sep='\t',names=["cell","reads*4","align_rate"])
    cell_reads['reads'] = cell_reads.apply(lambda x: x['reads*4']/4 ,axis=1)
    reads = cell_reads.drop(['reads*4'],axis=1)
    read = np.array(reads["reads"])
    align_rate = np.array(reads["align_rate"])

    # gene num
    gene = pd.read_table(expression,index_col='gene',sep='\t')
    gene_num = np.array(gene.apply(lambda x: count(x) ,axis=0))

    # merge data 
    reads["gene_num"] = gene_num
    x = range(0,50)
    r_num = np.log(np.array(reads.sort_values('reads',ascending=True)['reads']))
    g_num = np.log(np.array(reads.sort_values('reads',ascending=True)['gene_num']))
    a_rate = np.log(np.array(reads.sort_values('reads',ascending=True)['align_rate']))
    
    plot(read,align_rate,gene_num,x,r_num,g_num,a_rate,outputfile)
    
def plot(read,align_rate,gene_num,x,r_num,g_num,a_rate,outputfile):    # plot
    plt.figure(figsize=(20,10))
    #subplot(nrows, ncols, plot_number)

    plt.subplot(221)
    plt.violinplot(read)
    plt.ylabel("reads_count")

    plt.subplot(222)
    plt.violinplot(align_rate)
    plt.ylabel("align_rate")

    plt.subplot(223)
    plt.violinplot(gene_num)
    plt.ylabel("gene_num")
    
    plt.subplot(224)
    plt.plot(x, r_num,label='reads_num')
    plt.plot(x, g_num,label='gene_num')
    plt.plot(x, a_rate,label='align_rate')
    plt.xlabel("sample")
    plt.ylabel("count")
    plt.legend()
    #plt.show()
    plt.savefig(outputfile)

def main(argv):
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:e:o:",["ifile=","efile=","ofile="])
   except getopt.GetoptError:
      print("sc_plot.py -i <reads_align_stat file> -e <expression file> -o <outputfile>")
      sys.exit(2)
   for opt,arg in opts:
      if opt == '-h':
         print("sc_plot.py -i <reads_align_stat file> -e <expression file> -o <outputfile>")
         sys.exit()
      elif opt in ("-i", "--ifile"):
         reads_align_stat = arg
      elif opt in ("-e", "--efile"):
         expression = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
   print('输入的文件为：', reads_align_stat,expression)
   print('输出的文件为：', outputfile)
   process(reads_align_stat,expression,outputfile)
      

if __name__ == "__main__":
   main(sys.argv[1:])
