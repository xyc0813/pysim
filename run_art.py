import random
import sys
from optparse import OptionParser
import getopt
import time
import os
def run_art(ref,read_length,coverage,insert,sd,output):
    cmd='art_illumina -sam -i '+ref+' -p -l '+str(read_length)+' -f '+str(coverage)+' -m '+str(insert)+' -s '+str(sd)+' -o '+output
    os.system(cmd)
def main():
    usage = """%prog -i <file> -l <read length>  -f <coverage> -m <insert size> -s <sd of insert size> -o <out_fastq> 

sim_SV v1.1
Author: Yuchao Xia	
Description: run art
	"""   
    parser = OptionParser(usage)
    parser.add_option("-i", "--inFile", dest="inFile", help="A reference fasta file.",metavar="FILE")
    parser.add_option("-l",dest='read_length',help='paired-end read length',metavar='100')
    parser.add_option("-f",dest='coverage',help='read coverage',metavar='30')
    parser.add_option("-m",dest='insert',help='library insert size',metavar='350')
    parser.add_option("-s",dest='sd',help='sd of insert size',metavar='50')
    parser.add_option('-o','--output',dest='output',help='output fastq file',metavar="file")
    (opts, args) = parser.parse_args()
    if opts.inFile is None:
        parser.print_help()
    else:
        if opts.read_length is None:
            read_length=100
        else:
            read_length=int(opts.read_length)
        if opts.coverage is None:
            coverage=30
        else:
            coverage=int(opts.coverage)
        if opts.insert is None:
            insert=350
        else:
            insert=int(opts.insert)
        if opts.sd is None:
            sd=50
        else:
            sd=float(opts.sd)
        if opts.output is None:
            output='output'
        else:
            output=opts.output
        run_art(opts.inFile,read_length,coverage,insert,sd,output)
if __name__ == "__main__":
    print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    start = time.clock()
    print('simulation begins at:'+str(start))
    main()
    end = time.clock()
    print('simulation ends at:'+str(end))
    print("The function run time is : %.03f seconds" %(end-start))
    print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            
    
