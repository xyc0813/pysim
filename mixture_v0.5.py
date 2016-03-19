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
    outfile_1.close()
    
    outfile_2=open(outfilename+'_2.fq','w')
    outfile_2.close()
    
    filename=[]
    count=[]
    for line in open(config):
        if not line.startswith('#'):
            try:
                newline=line.rstrip().split('\t')
                print newline
                filename_1=newline[0]+'_1.fq'
                filename_2=newline[0]+'_2.fq'
                filename.append([filename_1,filename_2,float(newline[1])])
            except:
                print 'error'
                continue
        else:
            continue
    for line in filename:
        count.append([read_line(line[0])/4,line[-1]])
    print count
    sum_count=0
    for line in count:
        sum_count=sum_count+line[0]
    real_line=[int(line[0]*line[1]) for line in count]
    #for line in count:
        #real_line.append(int(line[1]*sum_count))
    print real_line
    tag=0
    for i in range(len(real_line)):
        if real_line[i]<=count[i][0]:
            continue
        else:
            tag=1
            break
    if tag==1:
        print 'proportion error'
    else:
        i=0
        #linecache.clearcache()
        m1=0
        m2=0
        for line in filename:
            print line
            l_list=[j for j in range(count[i][0]*4) if j%4==1]
            print 'random_shuttle'
            random.shuffle(l_list)
            print 'done'
            k=0
            
            for tmp_line in line[:2]:
                k=k+1
                for random_line in l_list[:real_line[i]]:
                    #print line
                    j=j+1
                    if k==1:
                        m1=m1+1
                        line1=linecache.getline(tmp_line,random_line)
                        line2=linecache.getline(tmp_line,random_line+1)
                        line3=linecache.getline(tmp_line,random_line+2)
                        line4=linecache.getline(tmp_line,random_line+3)
                        #print line1
                        #print line2
                        #print line3
                        #print line4
                        line1='@r'+str(m1)+'/1\n'
                        with open(outfilename+'_1.fq', 'a+') as outfile_1:
                            outfile_1.write(line1)
                            outfile_1.write(line2)
                            outfile_1.write(line3)
                            outfile_1.write(line4)
                    else:
                        m2=m2+1
                        line1=linecache.getline(tmp_line,random_line)
                        line2=linecache.getline(tmp_line,random_line+1)
                        line3=linecache.getline(tmp_line,random_line+2)
                        line4=linecache.getline(tmp_line,random_line+3)
                        line1='@r'+str(m2)+'/2\n'
                        with open(outfilename+'_2.fq', 'a+') as outfile_2:
                            outfile_2.write(line1)
                            outfile_2.write(line2)
                            outfile_2.write(line3)
                            outfile_2.write(line4)
            
                print 'clearcache'
                linecache.clearcache()
            #m1=m1+real_line[i]
            #m2=m2+real_line[i]
            i=i+1
    #outfile_1.close()
    #outfile_2.close()
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
    print time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    start = time.clock()
    print('simulation begins at:'+str(start))
    main()
    end = time.clock()
    print('simulation ends at:'+str(end))
    print("The function run time is : %.03f seconds" %(end-start))
    print time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            
                
