import random
from optparse import OptionParser
import getopt
import time
import os
import linecache
def read_line(filename):
    count = -1
    for count, line in enumerate(open(filename, 'rU')):
        pass
    count += 1
    return count
def readfq(filename):
    i=0
    tmp=[]
    fq_list=[]
    for line in open(filename):
        i=i+1
        if i%4==1 or i%4==2 or i%4==3:
            tmp.append(line.rstrip())
        else:
            tmp.append(line.rstrip())
            fq_list.append(tmp)
            tmp=[]
    return fq_list
        
def mixture(config,outfilename):
    outfile_1=open(outfilename+'_1.fq','w')
    outfile_1.truncate(0)
    
    outfile_2=open(outfilename+'_2.fq','w')
    outfile_2.truncate(0)
    outfiles=[outfile_1,outfile_2]
    
    filename=[]
    count=[]
    for line in open(config):
        if not line.startswith('#'):
            try:
                newline=line.rstrip().split()
                print( newline)
                filename_1=newline[0]+'1.fq'
                filename_2=newline[0]+'2.fq'
                filename.append([filename_1,filename_2,float(newline[1])])
            except:
                print( 'error parsing the config file') #lol what does this mean
                raise
        else:
            continue
    for line in filename:
        count.append([read_line(line[0])//4,line.pop()])
    # filename - [[germ_pair_read_1_filename, germ_pair_read_2][somatic_pair_read_1,somatic_pair_read_2]]
    print( 'count:',count) #[[number of reads, proportion]] for germline, somatic
    real_line=[int(line[0]*line[1]) for line in count] #number of reads from each file, w severe rounding error
    print( real_line)

    for i in range(len(real_line)):
        if real_line[i]>count[i][0]:
            raise ValueError('Proportions not less than 1')

    m=0
    for i in range(len(filename)): # for each sequencing run
        seqrun = filename[i]
        print('read pairs:', seqrun)
        l_list=list(range(1,count[i][0]*4,4))
        print('l_list',l_list)
        print( 'random_shuffle')
        random.shuffle(l_list)
        print( 'done')
        for random_line in l_list[:real_line[i]]:
            m+=1
            for k in range(len(seqrun)): # for each end of the paired end files
                paired_end = seqrun[k]            
                line2=linecache.getline(paired_end,random_line+1)
                line3=linecache.getline(paired_end,random_line+2)
                line4=linecache.getline(paired_end,random_line+3)
                line1='@r'+str(m)+'/'+str(k+1)+'\n'
                outfiles[k].write(line1)
                outfiles[k].write(line2)
                outfiles[k].write(line3)
                outfiles[k].write(line4)
        print( 'clearcache')
        linecache.clearcache()
        #m1=m1+real_line[i]
        #m2=m2+real_line[i]
    outfile_1.close()
    outfile_2.close()
def main():
     usage = """%prog -i <file>

simulate_copy_number v1.1
Author: Yuchao Xia	
Description: 
python """
     parser = OptionParser(usage)
     parser.add_option("-i", "--inFile", dest="inFile", help="config file.",metavar="FILE")
     parser.add_option('-o','--output',dest='output',help='output fastq file',metavar="file")
     (opts, args) = parser.parse_args()
     if opts.inFile is None:
         parser.print_help()
     else:
        if opts.output is None:
            out='output'
        else:
            out=opts.output
            mixture(opts.inFile,out)
if __name__ == "__main__":
    print( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    start = time.time()
    print('simulation begins at:'+str(start))
    #random.seed(42)
    main()
    end = time.time()
    print('simulation ends at:'+str(end))
    print("The function run time is : %.03f seconds" %(end-start))
    print( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            
                
