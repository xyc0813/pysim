import sys
from optparse import OptionParser
import getopt
import time
import os
def read_config(config_filename):
    dic={}
    for line in open(config_filename):
        if not line.startswith('#'):
            newline=line.rstrip().split('\t')
            if 'chr' != newline[0][:2]: newline[0] = 'chr'+str(newline[0])
            dic[newline[0]]=int(newline[1])
        else:
            continue
    return dic
def ploidy(ref_dic,dic):
    tmp_dic=[]
    ref_dic_new={}
    for key in dic:
        for n in range(dic[key]):
            ref_dic_new[key+'-hap-'+str(n+1)]=ref_dic[key]
    return ref_dic_new
def read_fasta(filename):
    ref_dic={}
    chr_name=''
    for line in open(filename):
        newline=line.rstrip()
        if newline.startswith('>'):
            if chr_name!='':
                if not chr_name.startswith('chr'):
                    chr_name='chr'+chr_name
                ref_dic[chr_name]=tmp_str
            chr_name=newline.split('>')[1].split(' ').pop(0)
            if not chr_name.startswith('chr'):
                chr_name='chr'+chr_name
            tmp_str=''
        else:
            tmp_str=tmp_str+newline.upper()
    ref_dic[chr_name]=tmp_str
    return ref_dic
def output(ref,outfilename):
    outfile=open(outfilename,'w')
    tmp_key=sorted(ref.keys())
    for key in tmp_key:
        i=0
        str_len=len(ref[key])
        outfile.write('>'+key+'\n')
        while i+50<=str_len:
            outfile.write(ref[key][i:i+50]+'\n')
            i=i+50
    outfile.close()
def main():
    usage = """%prog -i <file> -c <config file>  -o <out_fasta> 

specify_ploidy
Author: Yuchao Xia	
Description: specify the number of ploidy for different chromesome
	"""
        
    parser = OptionParser(usage)
    parser.add_option("-i", "--inFile", dest="inFile", help="A reference fasta file.",metavar="FILE")
    parser.add_option("-c","--config",dest='config',help='the number of ploidy of different chromesome',metavar='FILE')
    parser.add_option('-o','--output',dest='output',help='output fasta file',metavar="file")
   
    (opts, args) = parser.parse_args()
    if opts.inFile is None or opts.config is None:
        parser.print_help()
    else:
        if opts.output is None:
            outfilename='output_ploidy.fa'
        else:
            outfilename=opts.output
        config_dic=read_config(opts.config)
        ref_dic=read_fasta(opts.inFile)
        ref_dic=ploidy(ref_dic,config_dic)
        output(ref_dic,outfilename)
if __name__ == "__main__":
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    start = time.clock()
    print('simulation begins at:'+str(start))
    main()
    end = time.clock()
    print('simulation ends at:'+str(end))
    print("The function run time is : %.03f seconds" %(end-start))
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) )
