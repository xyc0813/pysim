# -*- coding: cp936 -*-
import sys,os
import random
import string
from math import log,exp
from optparse import OptionParser
import getopt
import time
######生成模拟GC含量文件###########
def fun(x):
    if x>=0.4:
        return -2*x+1.8
    elif 0.2<=x<0.4:
        return 4*x-0.6
    else:
        return 0.1
        
    
def generate_GC_file(line=100,outfilename='GC_content.txt',start=0,seq=0.01):
    outfile=open(outfilename,'w')
    for i in range(line):
        outfile.write(str(start)+'\t'+str(fun(start))+'\n')
        start=start+seq
    outfile.close()   

#####计算序列GC含量###############
def Compute_GC (seq):
    if 'N' in seq:
        G=seq.count('G')
        C=seq.count('C')
        N=seq.count('N')
        if len(seq)-N!=0:
            return float(G+C)/(len(seq)-N)
        else:
            return 0
    else:
        G=seq.count('G')
        C=seq.count('C')
        if len(seq)>0:
            return float(G+C)/len(seq)
        else:
            return 0
   
        
#####读取GC文件####
def GC_file(GCfile='GC_content.txt'):
    #GC_function=[]
    GC_function={}
    for line in open(GCfile):
        if not line.startswith('#'):
            newline=line.rstrip().split('\t')
            #GC_function.append([float(newline[0]),float(newline[1])])
            GC_function[float(newline[0])]=float(newline[1])
        else:
            continue
    return GC_function
    #return sorted(GC_function,key=lambda d:d[0])
######计算f(GC)#######
def GC_value(GC,GC_function):
    tmp=0.0
    tag=0
    sort_GC=sorted(GC_function.keys())
    if GC in GC_function.keys():
        return GC_function[GC]
    else:
        
        for line in sort_GC:
            if GC>line:
                tmp=GC_function[line]
                continue
            else:
                result=tmp
                tag=1
                break
        if tag==0:
            result=tmp
        return result
######读取ref序列##############                   
def reference(ref_name):
    ref_dic={}
    chr_name=''
    for line in open(ref_name):
        newline=line.rstrip()
        if newline.startswith('>'):
            if chr_name!='':
                ref_dic[chr_name]=tmp_str
            chr_name=newline.split('>')[1]
            tmp_str=''
        else:
            tmp_str=tmp_str+newline
    ref_dic[chr_name]=tmp_str
    return ref_dic
######以概率生成0,1################
def weighted_choice_sub(weights):
    rnd = random.random() * sum(weights)
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return i
def reverse(str1):
    str2=''
    for line in str1:
        if line=='A':
            str2=str2+'T'
        elif line=='G':
            str2=str2+'C'
        elif line=='C':
            str2=str2+'G'
        elif line=='T':
            str2=str2+'A'
        else:
            str2=str2+line
    return str2
#######读取bam文件，确定1,0###########
def read_name_sort_bam(filename,ref_dic,GC_function,outfile_name='tmp',meta='n',keepsam='n'):
    if filename.split('.')[-1]=='bam':
        filename_out=filename+'.sam'
        cmd='samtools view -h '+filename +' > '+filename_out
        os.system(cmd)
    else:
        filename_out=filename
    
    outfile1=open(outfile_name+'_1.fq','w')
    outfile2=open(outfile_name+'_2.fq','w')
    i=0
    tmp=[]
    if keepsam=='y' or keepsam=='Y':
        outfile=open(outfile_name+'_sam','w')
        for line in open(filename_out):
            if line.startswith('@'):
                outfile.write(line)
            else:
                i=i+1
                if i%2==1:
                    tmp.append(line.rstrip().split('\t'))
                else:
                    tmp.append(line.rstrip().split('\t'))
                    pois=[]
                    #print tmp
                    if tmp[0][0]==tmp[1][0]:
                        #if tmp[0][2].startswith('chr') and tmp[1][2].startswith('chr') and int(tmp[0][4])>0 and int(tmp[1][4])>0:
                        if tmp[0][2].startswith('chr') and tmp[1][2].startswith('chr'):
                            #chrome=tmp[0][2]+'_hap1'
                            if meta=='y' or meta=='Y':
                                if int(tmp[0][3])<int(tmp[1][3]):
                                    pois=[int(tmp[0][3]),int(tmp[1][3]),len(tmp[1][9])]
                                    chrome=tmp[0][2]
                                else:
                                    a=tmp[0]
                                    tmp[0]=tmp[1]
                                    tmp[1]=a
                                    chrome=tmp[0][2]
                                    pois=[int(tmp[0][3]),int(tmp[1][3]),len(tmp[1][9])]
                            
                            else:
                                if int(tmp[0][3])<int(tmp[1][3]):
                                    pois=[int(tmp[0][3]),int(tmp[1][3]),len(tmp[1][9])]
                                    tmp[1][9]=reverse(tmp[1][9][::-1].upper())
                                    tmp[1][10]=tmp[1][10][::-1]
                                    chrome=tmp[0][2]
                                else:
                                    a=tmp[0]
                                    tmp[0]=tmp[1]
                                    tmp[1]=a
                                    pois=[int(tmp[0][3]),int(tmp[1][3]),len(tmp[1][9])]
                                    tmp[1][9]=reverse(tmp[1][9][::-1].upper())
                                    tmp[1][10]=tmp[1][10][::-1]
                                    chrome=tmp[0][2]
                            
                        else:
                            #print ['wrong',tmp]
                            i=0
                            tmp=[]
                            continue
                    else:
                        #print ['wrong',tmp]
                        i=0
                        tmp=[]
                        continue
                    seq=ref_dic[chrome][pois[0]-1:pois[1]+pois[2]]
                    #print [ref_dic.keys(),seq,pois[0]-1,pois[1]+pois[2]]
                    GC_content=Compute_GC(seq.upper())
                    f_GC=GC_value(GC_content,GC_function)
                    prob=exp(f_GC)/(1+exp(f_GC))
                    index=weighted_choice_sub([1-prob,prob])
                    #print([GC_content,f_GC,prob,index,pois])
                    if index==1:
                        outfile.write('\t'.join(tmp[0])+'\n')
                        outfile.write('\t'.join(tmp[1])+'\n')
                        outfile1.write('@'+tmp[0][0]+'\n')
                        outfile1.write(tmp[0][9]+'\n')
                        outfile1.write('+\n')
                        outfile1.write(tmp[0][10]+'\n')
                        outfile2.write('@'+tmp[1][0]+'\n')
                        outfile2.write(tmp[1][9]+'\n')
                        outfile2.write('+\n')
                        outfile2.write(tmp[1][10]+'\n')
                        tmp=[]
                        i=0
                    else:
                        tmp=[]
                        i=0
        outfile.close()
    else:
        for line in open(filename_out):
            if line.startswith('@'):
                continue
            else:
                i=i+1
                if i%2==1:
                    tmp.append(line.rstrip().split('\t'))
                else:
                    tmp.append(line.rstrip().split('\t'))
                    pois=[]
                    #print tmp
                    if tmp[0][0]==tmp[1][0]:
                        #if tmp[0][2].startswith('chr') and tmp[1][2].startswith('chr') and int(tmp[0][4])>0 and int(tmp[1][4])>0:
                        if tmp[0][2].startswith('chr') and tmp[1][2].startswith('chr'):
                            #chrome=tmp[0][2]+'_hap1'
                            if meta=='y' or meta=='Y':
                                if int(tmp[0][3])<int(tmp[1][3]):
                                    pois=[int(tmp[0][3]),int(tmp[1][3]),len(tmp[1][9])]
                                    chrome=tmp[0][2]
                                else:
                                    a=tmp[0]
                                    tmp[0]=tmp[1]
                                    tmp[1]=a
                                    chrome=tmp[0][2]
                                    pois=[int(tmp[0][3]),int(tmp[1][3]),len(tmp[1][9])]
                            
                            else:
                                if int(tmp[0][3])<int(tmp[1][3]):
                                    pois=[int(tmp[0][3]),int(tmp[1][3]),len(tmp[1][9])]
                                    tmp[1][9]=reverse(tmp[1][9][::-1].upper())
                                    tmp[1][10]=tmp[1][10][::-1]
                                    chrome=tmp[0][2]
                                else:
                                    a=tmp[0]
                                    tmp[0]=tmp[1]
                                    tmp[1]=a
                                    pois=[int(tmp[0][3]),int(tmp[1][3]),len(tmp[1][9])]
                                    tmp[1][9]=reverse(tmp[1][9][::-1].upper())
                                    tmp[1][10]=tmp[1][10][::-1]
                                    chrome=tmp[0][2]
                            
                        else:
                            #print ['wrong',tmp]
                            i=0
                            tmp=[]
                            continue
                    else:
                        #print ['wrong',tmp]
                        i=0
                        tmp=[]
                        continue
                    seq=ref_dic[chrome][pois[0]-1:pois[1]+pois[2]]
                    #print [ref_dic.keys(),seq,pois[0]-1,pois[1]+pois[2]]
                    GC_content=Compute_GC(seq.upper())
                    f_GC=GC_value(GC_content,GC_function)
                    prob=exp(f_GC)/(1+exp(f_GC))
                    index=weighted_choice_sub([1-prob,prob])
                    #print([GC_content,f_GC,prob,index,pois])
                    if index==1:
                        outfile1.write('@'+tmp[0][0]+'\n')
                        outfile1.write(tmp[0][9]+'\n')
                        outfile1.write('+\n')
                        outfile1.write(tmp[0][10]+'\n')
                        outfile2.write('@'+tmp[1][0]+'\n')
                        outfile2.write(tmp[1][9]+'\n')
                        outfile2.write('+\n')
                        outfile2.write(tmp[1][10]+'\n')
                        tmp=[]
                        i=0
                    else:
                        tmp=[]
                        i=0
    
    outfile1.close()
    outfile2.close()

def main():
    usage = """%prog -i <file>

GC v1.0.1
Author: Yuchao Xia	
Description: accord GC content to filter paired-end reads from name sorted simulate bam file.
python GC1.py -i <ref.fa> -g <GC_file> -b <name_sort_bam> -o <out_fastq> -sort_bam<y/n> -m <y/n> -k <y/n>
	"""
    parser = OptionParser(usage)
    parser.add_option("-i", "--inFile", dest="inFile", help="A reference fasta file.",metavar="FILE")
    parser.add_option("-g","--GC",dest='input_GC',help='GC_function file,splited by \t,the filst col is GC content',metavar="FILE")
    parser.add_option("-b","--bamfile",dest='input_bam',help='input bam file sorted by name',metavar="FILE")
    parser.add_option('-s','--sort_bam',dest='sort_bam',help='if convert sam to sorted bam, print y',metavar="y/n")
    parser.add_option('-o','--output',dest='output',help='output fastq file',metavar="file")
    parser.add_option('-m','--metasim',dest='meta',help='input metasim data',metavar='y/n')
    parser.add_option('-k','--keepsam',dest='keepsam',help='generate sam file',metavar='y/n')
    (opts, args) = parser.parse_args()
    if opts.inFile is None or opts.input_bam is None:
        parser.print_help()
    else:
        if opts.input_GC == None:
            print 'done'
            generate_GC_file()
            ref_dic=reference(opts.inFile)
            GC_function=GC_file()
            read_name_sort_bam(opts.input_bam,ref_dic,GC_function,opts.output,opts.meta,opts.keepsam)
        else:
            print 'done1'
            ref_dic=reference(opts.inFile)
            GC_function=GC_file(opts.input_GC)
            read_name_sort_bam(opts.input_bam,ref_dic,GC_function,opts.output,opts.meta,opts.keepsam)
        #if opts.sort_bam=='y' or opts.sort_bam=='Y' :
            #cmd1='samtools view -S -h -b '+ops.input_bam.split('.')[0]+'.sam -o '+ops.input_bam.split('.')[0]+'_new.bam'
            #cmd2='samtools sort ' +ops.input_bam.split('.')[0]+'_new.bam' +ops.input_bam.split('.')[0]+'_sorted'
            #os.system(cmd1)
            #os.system(cmd2)
if __name__ == "__main__":
    start = time.clock()
    main()
    end = time.clock()
    print("The function run time is : %.03f seconds" %(end-start))
    
    
                    
                
                
                
                
                
                
                
            
    
