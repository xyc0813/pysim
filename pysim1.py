import random
import sys
from optparse import OptionParser
import getopt
import time
import os
#import vcf
def read_vcf(dbsnp,chrome):
    snp_dic={}
    vcf_reader = vcf.Reader(open(dbsnp, 'r'))
    if chrome=='ALL':
        for record in vcf_reader:
            if record.CHROM not in snp_dic:
                snp_dic[record.CHROM]=[[record.POS,record.REF,record.ALT[0]]]
            else:
                snp_dic[record.CHROM].append([record.POS,record.REF,record.ALT])
    else:
        chrome_list=chrome.split(',')
        print( chrome_list)
        for record in  vcf_reader:
            if record.CHROM in chrome_list:
                if record.CHROM not in snp_dic:
                    snp_dic[record.CHROM]=[[record.POS,record.REF,record.ALT[0]]]
                else:
                    snp_dic[record.CHROM].append([record.POS,record.REF,record.ALT])
            else:
                continue
    return snp_dic
    
def read_dbsnp(dbsnp,chome):
    snp_dic={}
    if chome=='ALL':
        for line in open(dbsnp):
            if not line.startswith('#'):
                newline=line.rstrip().split()#whitespace friendly
                if len(newline[4])==1 and len(newline[3])==1:
                    if newline[0] not in snp_dic:
                        snp_dic[newline[0]]=[[newline[1],newline[3],newline[4]]]
                    else:
                        snp_dic[newline[0]].append([newline[1],newline[3],newline[4]])
    else:
        chrome_list=chome.split(',')
        print( chrome_list)
        for line in open(dbsnp):
            if not line.startswith('#'):
                newline=line.rstrip().split()#whitespace friendly
                if len(newline[4])==1 and len(newline[3])==1 and newline[0] in chrome_list:
                    if newline[0] not in snp_dic:
                        snp_dic[newline[0]]=[[newline[1],newline[3],newline[4]]]
                    else:
                        snp_dic[newline[0]].append([newline[1],newline[3],newline[4]])
    return snp_dic
def compute_range(ref):
    ref_range={}
    for key in ref:
        n=0
        k=0
        for str1 in ref[key][0]: #char in first string (out of one)
            #n=n+k
            if str1=='N':
                n=n+1
                k=0
            else:
                k=k+1
                if k>=20:
                    break
                else:
                    continue
        end_n=0
        end_k=0
        for str1 in ref[key][0][::-1]: #char in first string (out of one) from end
            #end_n=end_n+end_k
            if str1=='N':
                end_n=end_n+1
                end_k=0
            else:
                end_k=end_k+1
                if end_k>=20:
                    break
                else:
                    continue
        ref_range[key]=[n,len(ref[key][0])-end_n]
        print( [key,n,len(ref[key][0])-end_n,end_n])
    return ref_range
def read_indel(indel,chrome):
    snp_dic={}
    if chome=='ALL':
        for line in open(dbsnp):
            if not line.startswith('#'):
                newline=line.rstrip().split()
                if len(newline[4])==1 and len(newline[3])==1:
                    if newline[0] not in snp_dic:
                        snp_dic[newline[0]]=[[newline[1],newline[3],newline[4]]]
                    else:
                        snp_dic[newline[0]].append([newline[1],newline[3],newline[4]])
                else:
                    if len(newline[3])==1 and len(newline[4])!=1:
                        if newline[0] not in snp_dic:
                            snp_dic[newline[0]]=[[newline[1],newline[3],newline[4],'ins']]
                        else:
                            snp_dic[newline[0]].append([newline[1],newline[3],newline[4],'ins'])
                    else:
                        if len(newline[4])==1 and len(newline[3])!=1:
                            if newline[0] not in snp_dic:
                                snp_dic[newline[0]]=[[newline[1],newline[3],newline[4],'del']]
                            else:
                                snp_dic[newline[0]].append([newline[1],newline[3],newline[4],'del'])
                        
                        
    else:
        for line in open(dbsnp):
            if not line.startswith('#'):
                newline=line.rstrip().split()
                if len(newline[4])==1 and len(newline[3])==1 and newline[0]==chome:
                    if newline[0] not in snp_dic:
                        snp_dic[newline[0]]=[[newline[1],newline[3],newline[4]]]
                    else:
                        snp_dic[newline[0]].append([newline[1],newline[3],newline[4]])
                else:
                    if len(newline[3])==1 and len(newline[4])!=1:
                        if newline[0] not in snp_dic:
                            snp_dic[newline[0]]=[[newline[1],newline[3],newline[4],'ins']]
                        else:
                            snp_dic[newline[0]].append([newline[1],newline[3],newline[4],'ins'])
                    else:
                        if len(newline[4])==1 and len(newline[3])!=1:
                            if newline[0] not in snp_dic:
                                snp_dic[newline[0]]=[[newline[1],newline[3],newline[4],'del']]
                            else:
                                snp_dic[newline[0]].append([newline[1],newline[3],newline[4],'del'])
    return snp_dic

def reference(ref_name):
    ref_dic={}
    chr_name=None
    for line in open(ref_name):
        newline=line.rstrip()
        if newline.startswith('>'):
            if '_' in newline:
                chr_split=newline.split('_')
            elif '-' in newline:
                chr_split=newline.split('-')
            else:
                chr_split=[newline]
            if chr_name!=None:
                if chr_name not in ref_dic:
                    ref_dic[chr_name]=[sequence]
                else:
                    ref_dic[chr_name].append(sequence)

            chr_name=chr_split[0].split('>')[1]
            sequence=''
        else:
            sequence=sequence+newline
    if chr_name not in ref_dic:
        ref_dic[chr_name]=[sequence]
    else:
        ref_dic[chr_name].append(sequence)
    return ref_dic

def generate_normal(ref_dic,snp_dic,num,outfilename,hyp_rate=0.5):
    '''
    ref_dic: length 1. key: value -> chrom_num: [ref_seq] (in this case GRCh38_chr22.fasta). list of strings, len 1
    should be both haplotypes?
    snp_dic: length 1. key:value -> chrom_num: dbsnp, in the form [[pos, ref, alt]]
    ^ why are these dicts?
    num = number of SNPs we want to make
    return: ref_dic: dic containing chrom_num: [modified seq]
        snp_list: list of snps introduced, one-indexed to ref_dic
    
    '''
    outfile=open(outfilename,'w')
    outfile.write('#chr\tpois\tref\talt\thap\n')
    snp_list={}
    l1=[1 for i in range(int(num*hyp_rate))] # proportion homozygous
    l2=[0 for i in range(num-int(num*hyp_rate))] #proportion heterozygous
    total_list=l1+l2
    random.shuffle(total_list)
    num=int(num/len(ref_dic.keys()))
    i=0
    for key in ref_dic: #for each chrom
        hapl=list(range(len(ref_dic[key]))) # [0]
        alt_seq=ref_dic[key] # ref_seq as string. Will be subbing in SNPs

        #determine which positions we will make snps
        if len(snp_dic[key])>=num: # num of snps in db > num snps desired -> almost always true
            snp_positions=random.sample(snp_dic[key],num) # random set of snps from db
        else:
            snp_positioins=snp_dic[key]

        # make the snps
        for line in snp_positions:
            #raise Exception('this code is fucking awful')
            if total_list[i]==0:
                hap=random.sample(hapl,1)[0] # pick a 'N' in ref_seq
                index = int(line[0])-1
                old=alt_seq[hap][index]
                alt_seq[hap]=alt_seq[hap][:index]+line[-1]+alt_seq[hap][index+1:] # make an snp
                outfile.write(key+'\t'+line[0]+'\t'+old+'\t'+line[-1]+'\t'+str(hap+1)+'\n')
            else:
                for hap in hapl:
                    index = int(line[0])-1
                    old=alt_seq[hap][index]
                    alt_seq[hap]=alt_seq[hap][:index]+line[-1]+alt_seq[hap][index+1:] # make an snp
                outfile.write(key+'\t'+line[0]+'\t'+old+'\t'+line[-1]+'\thomozygous\n')
            i=i+1
            if key not in snp_list:
                snp_list[key]=[int(line[0])]
            else:
                snp_list[key].append(int(line[0]))
        ref_dic[key]=alt_seq
    outfile.close()
    return (ref_dic,snp_list)

def generate_somatic(ref_dic,snp_dic,num,outfilename,db,ref_range,hyp_rate=0.5):
    outfile=open(outfilename,'w')
    outfile.write('#chr\tpois\tref\talt\thap\n')
    l1=[1 for i in range(int(num*hyp_rate))]
    l2=[0 for i in range(num-int(num*hyp_rate))]
    total_list=l1+l2
    random.shuffle(total_list)
    ref_list=[]
    num=int(num/len(ref_dic.keys()))
    snp_list={}
    if db is 'N':
        i=0
        for key in ref_dic:
            snp_list[key]=[]
            hapl=[k for k in range(len(ref_dic[key]))]
            tmp_str_list=ref_dic[key]
            str_list=[]
            ref_range_list=ref_range[key]
            for tmp in tmp_str_list:
                str_list.append(list(tmp))
            for j in range(num):
                while True:
                    pois=random.randint(ref_range_list[0],ref_range_list[1])
                    if total_list[i]==0:
                        hap=random.sample(hapl,1)[0]
                        try:
                            if str_list[hap][pois-1]!='N':
                                old=str_list[hap][pois-1]
                                str_list[hap][pois-1]=SNP(str_list[hap][pois-1].upper())
                                outfile.write(key+'\t'+str(pois)+'\t'+old+'\t'+str_list[hap][pois-1]+'\t'+str(hap+1)+'\n')
                                i=i+1
                                snp_list[key].append(pois)
                                break
                        except:
                            print( [pois,hap])
                            continue
                    else:
                        try:
                            if str_list[0][pois-1]!='N':
                                old=str_list[0][pois-1]
                                new=SNP(str_list[0][pois-1].upper())
                                for t in hapl:
                                    str_list[t][pois-1]=new
                                outfile.write(key+'\t'+str(pois)+'\t'+old+'\t'+new+'\thomozygous\n')
                                i=i+1
                                snp_list[key].append(pois)
                                break
                        except:
                            print( [pois])
                            continue
            tmp_str=[]
            for tmp in str_list:
                tmp_str.append(''.join(tmp))
            ref_dic[key]=tmp_str                                
    else:
        i=0
        for key in ref_dic:
            hapl=[k for k in range(len(ref_dic[key]))]
            tmp_str_list=ref_dic[key]
            str_list=[]
            for tmp in tmp_str_list:
                str_list.append(list(tmp))
            if len(snp_dic[key])>=num:
                random_list=[]
                jj=0
                while jj<num:
                    while True:
                        tmp_snp_list=random.sample(snp_dic[key],1)[0]
                        if tmp_snp_list[0] not in snp_list[key]:
                            random_list.append(tmp_snp_list)
                            jj=jj+1
                            break
                        else:
                            continue
                for line in random_list:
                    if total_list[i]==0:
                        hap=random.sample(hapl,1)[0]
                        str_list[hap][int(line[0])-1]=line[-1]
                        outfile.write(key+'\t'+line[0]+'\t'+str_list[hap][int(line[0])-1]+'\t'+line[-1]+'\t'+str(hap+1)+'\n')
                    else:
                        old=str_list[0][int(line[0])-1]
                        new=SNP(str_list[t][int(line[0])-1].upper())
                        for t in hapl:
                            str_list[t][int(line[0])-1]=line[-1]
                        outfile.write(key+'\t'+line[0]+'\t'+old+'\t'+line[-1]+'\thomozygous\n')
                    i=i+1
                    if key not in snp_list:
                        snp_list[key]=[int(line[0])]
                    else:
                        snp_list[key].append(int(line[0]))
            else:
                for line in snp_dic[key]:
                    if line[0] not in snp_list[key]:
                        if total_list[i]==0:
                            hap=random.sample(hapl,1)[0]
                            str_list[hap][int(line[0])-1]=line[-1]
                            outfile.write(key+'\t'+line[0]+'\t'+str_list[hap][int(line[0])-1]+'\t'+line[-1]+'\t'+str(hap+1)+'\n')
                        else:
                            for t in hapl:
                                str_list[t][int(line[0])-1]=line[-1]
                            outfile.write(key+'\t'+line[0]+'\t'+str_list[t][int(line[0])-1]+'\t'+line[-1]+'\thomozygous\n')
                        i=i+1       
                        if key not in snp_list:
                            snp_list[key]=[int(line[0])]
                        else:
                            snp_list[key].append(int(line[0]))
            tmp_str=[]
            for tmp in str_list:
                tmp_str.append(''.join(tmp))
            ref_dic[key]=tmp_str
    outfile.close()
    return [ref_dic,snp_list]           
                    
            
        
def SNP(str1):
    l=['A','T','C','G','N']
    l.remove(str1)
    while True:
        tmp=random.sample(l,1)[0]
        if tmp!='N':
            break
        else:
            continue
    return tmp
def random_snp(l,n,up_down_stream):
    random_l=[i for i in range(up_down_stream)]
    random_int=random.sample(random_l,n)
    for line in random_int:
        try:
            l[line]=SNP(l[line].upper())
        except:
            continue
            #print ['error',line,l]
    return l
def intersection(l1,l2):
    tag=0
    for line in l2:
        if l1[0]<=line[0]<=l1[1]<=line[1] or line[0]<=l1[0]<=line[1]<=l1[1] or l1[0]<=line[0]<=line[1]<=l1[1] or line[0]<=l1[0]<=l1[1]<=line[1]:
            tag=1
            break
        else:
            continue
    return tag    
def deletion(ref,chrome,pois,length,snp_rate,outfile,specify,up_down_stream=50,indel_prob=1,min_indel_length=5,max_indel_length=15):
    list_str=[]
    for line in ref[chrome]:
        list_str.append(list(line))
    hapl=[i for i in range(len(ref[chrome]))]
    if specify>2:
        hapl_type=random.sample(hapl,1)[0]
        hap=random.sample(hapl,hapl_type+1)
        hap_str=[str(i+1) for i in hap]
        outfile.write(chrome+'\t'+str(pois)+'\t'+str(pois+length)+'\t'+str(length)+'\tdeletion\t'+','.join(hap_str)+'\n')
    elif specify==0:
        hap=random.sample(hapl,len(hapl))
        hap_str=[str(i+1) for i in hap]
        #outfile.write(chrome+'\t'+str(pois)+'\t'+str(pois+length)+'\t'+str(length)+'\tCNV\t'+','.join(hap_str)+'\n')
        #hap_str=['0' for i in range(len(hapl))]
        #outfile.write(chrome+'\t'+str(pois)+'\t'+str(pois+length)+'\t'+str(length)+'\tCNV\t'+','.join(hap_str)+'\n')
        outfile.write(chrome+'\t'+str(pois)+'\t'+str(pois+length)+'\t'+str(length)+'\tCNV\t')
        for t in hap_str:
            outfile.write(t+':0;')
        outfile.write('\n')
    elif specify==1:
        hap=random.sample(hapl,1)
        hap_str=[str(i+1) for i in hap]
        #outfile.write(chrome+'\t'+str(pois)+'\t'+str(pois+length)+'\t'+str(length)+'\tCNV\t'+','.join(hap_str)+'\n')
        outfile.write(chrome+'\t'+str(pois)+'\t'+str(pois+length)+'\t'+str(length)+'\tCNV\t')
        new_hapl=[str(j+1) for j in hapl]
        for t in hap_str:
            outfile.write(t+':0;')
            new_hapl.remove(t)
        for t in new_hapl:
            outfile.write(t+':1;')
        outfile.write('\n')
    else:
        hap=random.sample(hapl,1)
        hap_str=[str(i+1) for i in hap]
        #outfile.write(chrome+'\t'+str(pois)+'\t'+str(pois+length)+'\t'+str(length)+'\tCNV\t'+','.join(hap_str)+'\n')
        outfile.write(chrome+'\t'+str(pois)+'\t'+str(pois+length)+'\t'+str(length)+'\tCNV\t')
        for t in hap_str:
            outfile.write(t+':0;')
        for line in hap:
            l1=list_str[line][:pois-1]
            l2=list_str[line][pois-1:pois-1+length]
            l3=list_str[line][pois-1+length:]
            snp_number=int(snp_rate*up_down_stream)
            if snp_number>0:
                l1=random_snp(l1[::-1],snp_number,up_down_stream)[::-1]
                l3=random_snp(l3,snp_number,up_down_stream)
            if indel_prob==1:
                str1=[]
                indel_length=random.randint(min_indel_length,max_indel_length)
                for i in range(indel_length):
                    str1.append(random.sample(['A','T','C','G'],1)[0])
                str2=[]
                indel_length=random.randint(min_indel_length,max_indel_length)
                for i in range(indel_length):
                    str2.append(random.sample(['A','T','C','G'],1)[0])
                l1=l1+str1
                l3=str2+l3
                ref[chrome][line]=''.join(l1+l3)
            else:
                ref[chrome][line]=''.join(l1+l3)
    if specify!=2:
        for line in hap:
            l1=list_str[line][:pois-1]
            l2=list_str[line][pois-1:pois-1+length]
            l3=list_str[line][pois-1+length:]
            snp_number=int(snp_rate*up_down_stream)
            if snp_number>0:
                l1=random_snp(l1[::-1],snp_number,up_down_stream)[::-1]
                l3=random_snp(l3,snp_number,up_down_stream)
            if indel_prob==1:
                str1=[]
                indel_length=random.randint(min_indel_length,max_indel_length)
                for i in range(indel_length):
                    str1.append(random.sample(['A','T','C','G'],1)[0])
                str2=[]
                indel_length=random.randint(min_indel_length,max_indel_length)
                for i in range(indel_length):
                    str2.append(random.sample(['A','T','C','G'],1)[0])
                l1=l1+str1
                l3=str2+l3
                ref[chrome][line]=''.join(l1+l3)
            else:
                ref[chrome][line]=''.join(l1+l3)
    else:
        hapl.remove(hap[0])
        hap_t=random.sample(hapl,1)[0]
        list_str=[]
        for line in ref[chrome]:
            list_str.append(list(line))
        l1=list_str[hap_t][:pois-1]
        l3=list_str[hap_t][pois-1:]
        #l2=list_str[hap_t][pois-1:pois-1+length]
        outfile.write(str(hap_t+1)+':2\n')
        ref[chrome][hap_t]=''.join(l1+l2+l3)
    return ref
def insertion(ref,chrome,pois,length,snp_rate,outfile,up_down_stream=50,indel_prob=1,min_indel_length=5,max_indel_length=15):
    list_str=[]
    for line in ref[chrome]:
        list_str.append(list(line))
    hapl=[i for i in range(len(ref[chrome]))]
    hapl_type=random.sample(hapl,1)[0]
    hap=random.sample(hapl,hapl_type+1)
    hap_str=[str(i+1) for i in hap]
    outfile.write(chrome+'\t'+str(pois)+'\t'+str(pois+length)+'\t'+str(length)+'\tinversion\t'+','.join(hap_str)+'\n')
    l2=[]
    for i in range(length):
        l2.append(random.sample(['A','T','C','G'],1)[0])
    for line in hap:
        
        l1=list_str[line][:pois-1]
        l3=list_str[line][pois-1:]
    
        snp_number=int(snp_rate*up_down_stream)
        if snp_number>0:
            l1=random_snp(l1[::-1],snp_number,up_down_stream)[::-1]
            l3=random_snp(l3,snp_number,up_down_stream)
        ref[chrome][line]=''.join(l1+l2+l3)
    return ref
def translocation(ref,chrome,pois,length,str1_ins,snp_rate,outfile,up_down_stream=50,indel_prob=1,min_indel_length=5,max_indel_length=15):
    list_str=[]
    for line in ref[chrome]:
        list_str.append(list(line))
    hapl=[i for i in range(len(ref[chrome]))]
    hapl_type=random.sample(hapl,1)[0]
    hap=random.sample(hapl,hapl_type+1)
    hap_str=[str(i+1) for i in hap]
    outfile.write(chrome+'\t'+str(pois)+'\t'+str(pois+length)+'\t'+str(length)+'\tinversion\t'+','.join(hap_str)+'\n')
    for line in hap:
        l1=list_str[line][:pois-1]
        l2=list_str[line][pois-1:pois-1+length]
        l3=list_str[line][pois-1+length:]
        snp_number=int(snp_rate*up_down_stream)
        if snp_number>0:
            l1=random_snp(l1[::-1],snp_number,up_down_stream)[::-1]
            l3=random_snp(l3,snp_number,up_down_stream)
        if indel_prob==1:
            str1=[]
            indel_length=random.randint(min_indel_length,max_indel_length)
            for i in range(indel_length):
                str1.append(random.sample(['A','T','C','G'],1)[0])
            str2=[]
            indel_length=random.randint(min_indel_length,max_indel_length)
            for i in range(indel_length):
                str2.append(random.sample(['A','T','C','G'],1)[0])
            l1=l1+str1
            l3=str2+l3
            ref[chrome][line]=''.join(l1+str1_ins+l3)
        else:
            ref[chrome][line]=''.join(l1+str1_ins+l3)
    return ref
def inversion(ref,chrome,pois,length,snp_rate,outfile,up_down_stream=50,indel_prob=1,min_indel_length=5,max_indel_length=15):
    list_str=[]
    for line in ref[chrome]:
        list_str.append(list(line))
    hapl=[i for i in range(len(ref[chrome]))]
    hapl_type=random.sample(hapl,1)[0]
    hap=random.sample(hapl,hapl_type+1)
    hap_str=[str(i+1) for i in hap]
    outfile.write(chrome+'\t'+str(pois)+'\t'+str(pois+length)+'\t'+str(length)+'\tinversion\t'+','.join(hap_str)+'\n')
    for line in hap:
        l1=list_str[line][:pois-1]
        l2=list_str[line][pois-1:pois-1+length]
        l3=list_str[line][pois-1+length:]
        snp_number=int(snp_rate*up_down_stream)
        if snp_number>0:
            l1=random_snp(l1[::-1],snp_number,up_down_stream)[::-1]
            l3=random_snp(l3,snp_number,up_down_stream)
        if indel_prob==1:
            str1=[]
            indel_length=random.randint(min_indel_length,max_indel_length)
            for i in range(indel_length):
                str1.append(random.sample(['A','T','C','G'],1)[0])
            str2=[]
            indel_length=random.randint(min_indel_length,max_indel_length)
            for i in range(indel_length):
                str2.append(random.sample(['A','T','C','G'],1)[0])
            l1=l1+str1
            l3=str2+l3
            ref[chrome][line]=''.join(l1+l2[::-1]+l3)
        else:
            ref[chrome][line]=''.join(l1+l2[::-1]+l3)
    return ref
def tandem_duplication(ref,chrome,pois,length,snp_rate,outfile,up_down_stream=50,indel_prob=1,min_indel_length=5,max_indel_length=15,tandem_num=2):
    list_str=[]
    for line in ref[chrome]:
        list_str.append(list(line))
    hapl=[i for i in range(len(ref[chrome]))]
    hapl_type=random.sample(hapl,1)[0]
    hap=random.sample(hapl,hapl_type+1)
    hap_str=[str(i+1) for i in hap]
    outfile.write(chrome+'\t'+str(pois)+'\t'+str(pois+length)+'\t'+str(length)+'\ttandem_duplication\t'+','.join(hap_str)+'\n')
    for line in hap:
        l1=list_str[line][:pois-1]
        l2=list_str[line][pois-1:pois-1+length]
        l3=list_str[line][pois-1+length:]
        snp_number=int(snp_rate*up_down_stream)
        if snp_number>0:
            l1=random_snp(l1[::-1],snp_number,up_down_stream)[::-1]
            l3=random_snp(l3,snp_number,up_down_stream)
        if indel_prob==1:
            str1=[]
            indel_length=random.randint(min_indel_length,max_indel_length)
            for i in range(indel_length):
                str1.append(random.sample(['A','T','C','G'],1)[0])
            str2=[]
            indel_length=random.randint(min_indel_length,max_indel_length)
            for i in range(indel_length):
                str2.append(random.sample(['A','T','C','G'],1)[0])
            l1=l1+str1
            l3=str2+l3
            ref[chrome][line]=''.join(l1+l2*tandem_num+l3)
        else:
            ref[chrome][line]=''.join(l1+l2*tandem_num+l3)
    for line in ref[chrome]:
        print( len(line))
    return ref
def CNV(ref,chrome,pois,length,result,copy_number,snp_rate,up_down_stream=50,indel_prob=1,min_indel_length=5,max_indel_length=15):
    list_str=[]
    for line in ref[chrome]:
        list_str.append(list(line))
    hapl=[i for i in range(len(ref[chrome]))]
    hapl_type=random.sample(hapl,1)[0]
    CNV_str=list_str[hapl_type][pois-1:pois+length]
    if chrome not in result:
        result[chrome]=[[CNV_str,pois,length,hapl_type,copy_number]]
    else:
        result[chrome].append([CNV_str,pois,length,hapl_type,copy_number])
    return result
def read_config(filename):
    sv_list=[]
    for line in open(filename):
        if not line.startswith('#'):
            newline=line.rstrip().split()
            if newline[0].upper()!='CNV':
                for i in range(int(newline[-1])):
                    sv_list.append([int(newline[1]),newline[0]])
            else:
                for i in range(int(newline[-1])):
                    sv_list.append([int(newline[1]),int(newline[2]),newline[0]])
        else:
            continue
    return sv_list
def generate_pois(len_list,ref,ref_range):
    dic={}
    key=ref.keys()
    
    #pois_list=[]
    pois_dic={}
    for tmp_key in key:
        if tmp_key not in dic:
            dic[tmp_key]=[]
        else:
            continue
    #str_len=[len(ref[tmp][0]) for tmp in dic]
    for i in range(len(len_list)):
        k=i%len(key)
        count=1
        tmp_count=0
        while True:
            pois=random.randint(ref_range[key[k]][0],ref_range[key[k]][1])
            tmp=[pois,pois+len_list[i][0]]
            if ref_range[key[k]][0]<=pois+len_list[i][0]<=ref_range[key[k]][1]:
                #if pois_list==[]:
                #    pois_list.append(tmp)
                if key[k] not in pois_dic:
                    pois_dic[key[k]]=[tmp]
                else:
                    #if intersection(tmp,pois_list)==0:
                    if intersection(tmp,pois_dic[key[k]])==0:
                        if len_list[i][-1]!='CNV' and len_list[i][-1] !='trans':
                            dic[key[k]].append([pois,len_list[i][0],len_list[i][1]])
                            pois_dic[key[k]].append(tmp)
                            break
                        elif len_list[i][-1]=='CNV':
                            dic[key[k]].append([pois,len_list[i][0],len_list[i][1],len_list[i][-1]])
                            pois_dic[key[k]].append(tmp)
                            break
                        else:
                            #dic[key[k]].append([pois,len_list[i][0],len_list[i][1]])
                            str1=ref[key[k]][pois-1:pois+len_list[i][0]]
                            pois_dic[key[k]].append(tmp)
                            chrome_new=random.sample(key,1)[0]
                            while True:
                                pois_1=random.randint(ref_range[chrome_new][0],ref_range[chrome_new][1])
                                tmp=[pois_1,pois_1+len_list[i][0]]
                                if pois_1+len_list[i][0]<=ref_range[key[k]][1]:
                                    if chrome_new not in pois_dic:
                                        pois_dic[chrome_new]=[tmp]
                                    else:
                                        if intersection(tmp,pois_dic[chrome_new])==0:
                                            str2=ref[chrome_new][pois_1-1:pois_1+len_list[i][0]]
                                            
                                            pois_dic[chrome_new].append(tmp)
                                            break
                                        else:
                                            continue
                            dic[key[k]].append([pois,len_list[i][0],len_list[i][1],str2,len_list[i][1]])
                            dic[chrome_new].append([pois_1,len_list[i][0],str1,len_list[i][-1]])
                            break
                    else:
                        count=count+1
                        if count<30:
                            continue
                        else:
                            tmp_count=tmp_count+1
                            if tmp_count>=5:
                                break
                            else:
                                len_list[i][0]=int(len_list[i][0]/2)
            else:
                continue
    for key in dic:
        dic[key]=sorted(dic[key],key=lambda d:d[0],reverse=True)
    return dic
def generate_subclone_pois(len_list,ref,ref_range,germline_dic):
    dic={}
    key=ref.keys()
    #pois_list=[]
    pois_dic={}
    for t_key in germline_dic:
        for line in germline_dic[t_key]:
            if t_key not in pois_dic:
                pois_dic[t_key]=[[int(line[0]),int(line[0]+line[1])]]
            else:
                pois_dic[t_key].append([int(line[0]),int(line[0]+line[1])])
    for tmp_key in key:
        if tmp_key not in dic:
            dic[tmp_key]=[]
        else:
            continue
    #str_len=[len(ref[tmp][0]) for tmp in dic]
    for i in range(len(len_list)):
        k=i%len(key)
        #print [k,key,key[k]]
        while True:
            count=1
            pois=random.randint(ref_range[key[k]][0],ref_range[key[k]][1])
            tmp=[pois,pois+len_list[i][0]]
            if ref_range[key[k]][0]<=pois+len_list[i][0]<=ref_range[key[k]][1]:
                if key[k] not in pois_dic:
                    pois_dic[key[k]]=[tmp]
                else:
                    if intersection(tmp,pois_dic[key[k]])==0:
                        if len_list[i][-1]!='CNV' and len_list[i][-1] !='trans':
                            dic[key[k]].append([pois,len_list[i][0],len_list[i][1]])
                            pois_dic[key[k]].append(tmp)
                            break
                        elif len_list[i][-1]=='CNV':
                            dic[key[k]].append([pois,len_list[i][0],len_list[i][1],len_list[i][-1]])
                            pois_dic[key[k]].append(tmp)
                            break
                        else:
                            #dic[key[k]].append([pois,len_list[i][0],len_list[i][1]])
                            str1=ref[key[k]][pois-1:pois+len_list[i][0]]
                            pois_dic[key[k]].append(tmp)
                            chrome_new=random.sample(key,1)[0]
                            while True:
                                pois_1=random.randint(ref_range[chrome_new][0],ref_range[chrome_new][1])
                                tmp=[pois_1,pois_1+len_list[i][0]]
                                if pois_1+len_list[i][0]<=ref_range[chrome_new][1]:
                                    if chrome_new not in pois_dic:
                                        pois_dic[chrome_new]=[tmp]
                                    else:
                                        if intersection(tmp,pois_dic[chrome_new])==0:
                                            str2=ref[chrome_new][pois_1-1:pois_1+len_list[i][0]]
                                            
                                            pois_dic[chrome_new].append(tmp)
                                            break
                                        else:
                                            continue
                            dic[key[k]].append([pois,len_list[i][0],len_list[i][1],str2,len_list[i][1]])
                            dic[chrome_new].append([pois_1,len_list[i][0],str1,len_list[i][-1]])
                            break
                    else:
                        count=count+1
                        if count<20:
                            continue
                        else:
                            break
            else:
                continue
    for key in dic:
        dic[key]=sorted(dic[key],key=lambda d:d[0],reverse=True)
    return dic
def intersection1(pois,l):
    tag=0
    for line in l:
        if line[0]<=pois<=line[0]+line[1]:
            tag=1
            break
        else:
            continue
    return tag
def generate_CNV(ref,CNV_str,dic,ref_range,outfile):
    dic_chr={}
    dic_pois={}
    for key in CNV_str:
        dic_pois={}
        list_str=[]
        for str1 in ref[key]:
            list_str.append(list(str1))
        hapl=[i for i in range(len(ref[key]))]
        hapl_type_list=[]
        
        for line in CNV_str[key]:
            dic_chr={}
            for tmp_hapl in hapl:
                dic_chr[tmp_hapl]=0
            for i in range(int(line[-1])-2):
                hapl_type=random.sample(hapl,1)[0]
                while True:
                    pois=random.randint(ref_range[key][0],ref_range[key][1])
                    if intersection1(pois,dic[key])==0:
                        if hapl_type not in dic_pois:
                            dic_pois[hapl_type]=[pois]
                        else:
                            dic_pois[hapl_type].append(pois)
                        if hapl_type not in dic_chr:
                            dic_chr[hapl_type]=1
                            if hapl_type not in hapl_type_list:
                                hapl_type_list.append(hapl_type)
                        else:
                            dic_chr[hapl_type]=dic_chr[hapl_type]+1
                            if hapl_type not in hapl_type_list:
                                hapl_type_list.append(hapl_type)
                        list_str[hapl_type]=list_str[hapl_type][:pois]+line[0]+list_str[hapl_type][pois:]
                        break
                    else:
                        continue
            outfile.write(key+'\t'+str(line[1])+'\t'+str(line[1]+line[2])+'\t'+str(line[2])+'\tCNV\t')
            for tmp_hapl in hapl:
                outfile.write(str(tmp_hapl+1)+':'+str(dic_chr[tmp_hapl]+1)+';')
                
            outfile.write('\n')
        str_list=[]
        for list_line in list_str:
            str_list.append(''.join(list_line))
        ref[key]=str_list
    return ref
   
def generate_fasta(ref,dic,snp_rate,outfile,ref_range,up_down_stream=50,indel_prob=1,min_indel_length=5,max_indel_length=15):
    CNV_str={}
    if indel_prob==1 or indel_prob==0:
        for key in dic:
            for line in dic[key]:
                print( line)
                if line[-1]=='del':
                    ref=deletion(ref,key,line[0],line[1],snp_rate,outfile,3,up_down_stream,indel_prob,min_indel_length,max_indel_length)
                elif line[-1]=='inv':
                    ref=inversion(ref,key,line[0],line[1],snp_rate,outfile,up_down_stream,indel_prob,min_indel_length,max_indel_length)
                elif line[-1]=='ins':
                    ref=insertion(ref,key,line[0],line[1],snp_rate,outfile,up_down_stream,indel_prob,min_indel_length,max_indel_length)
                elif line[-1]=='tandem':
                    tandem_num=random.sample([2,3,4,5],1)[0]
                    ref=tandem_duplication(ref,key,line[0],line[1],snp_rate,outfile,up_down_stream,indel_prob,min_indel_length,max_indel_length,tandem_num)
                elif line[-1].upper()=='CNV':
                    if line[-2]==0:
                        ref= deletion(ref,key,line[0],line[1],snp_rate,outfile,0,up_down_stream=50,indel_prob=1,min_indel_length=5,max_indel_length=15)
                    elif line[-2]==1:
                        ref=deletion(ref,key,line[0],line[1],snp_rate,outfile,1,up_down_stream=50,indel_prob=1,min_indel_length=5,max_indel_length=15)
                    elif line[-2]==2:
                        ref=deletion(ref,key,line[0],line[1],snp_rate,outfile,2,up_down_stream=50,indel_prob=1,min_indel_length=5,max_indel_length=15)
                    else:
                        CNV_str=CNV(ref,key,line[0],line[1],CNV_str,line[-2],snp_rate,up_down_stream=50,indel_prob=1,min_indel_length=5,max_indel_length=15)
                elif line[-1]=='trans':
                    ref=translocation(ref,key,line[0],line[1],line[2],snp_rate,outfile,up_down_stream,indel_prob,min_indel_length,max_indel_length)
                else:
                    continue
                
    else:
        for key in dic:
            indel_num=int(indel_prob*len(dic[key]))
            k=0
            for line in dic[key]:
                k=k+1
                if k<=indel_num:
                    tmp=1
                else:
                    tmp=0
                if line[-1]=='del':
                    ref=deletion(ref,key,line[0],line[1],snp_rate,outfile,3,up_down_stream,tmp,min_indel_length,max_indel_length)
                elif line[-1]=='inv':
                    ref=inversion(ref,key,line[0],line[1],snp_rate,outfile,up_down_stream,tmp,min_indel_length,max_indel_length)
                elif line[-1]=='ins':
                    ref=insertion(ref,key,line[0],line[1],snp_rate,outfile,up_down_stream,tmp,min_indel_length,max_indel_length)
                elif line[-1]=='tandem':
                    tandem_num=random.sample([2,3,4,5],1)[0]
                    ref=tandem_duplication(ref,key,line[0],line[1],snp_rate,outfile,up_down_stream,tmp,min_indel_length,max_indel_length,tandem_num)
                elif line[-1].upper()=='CNV':
                    CNV_str=CNV(ref,key,line[0],line[1],CNV_str,line[2],snp_rate,up_down_stream=50,indel_prob=1,min_indel_length=5,max_indel_length=15)
                elif line[-1]=='trans':
                    ref=translocation(ref,key,line[0],line[1],line[2],snp_rate,outfile,up_down_stream,indel_prob,min_indel_length,max_indel_length)
                else:
                    print( line[-1])
                    continue
    if CNV_str!={}:
        ref=generate_CNV(ref,CNV_str,dic,ref_range,outfile)
        
    return ref
def output(ref,dic,outfilename):
    outfile=open(outfilename,'w')
    for key in ref:
        k=0
        for line in ref[key]:
            i=0
            k=k+1
            str_len=len(line)
            print( str_len)
            outfile.write('>'+key+'_'+str(k)+'\n')
            while i+50<=str_len:
                outfile.write(line[i:i+50]+'\n')
                i=i+50
    outfile.close()
    outfile=open(outfilename.split('.')[0]+'_SV.txt','w')
    outfile.write('#chr\tstart_pois\tend_pois\tlength\ttype\n')
    for key in dic:
        for line in dic[key]:

            outfile.write(key+'\t'+str(line[0])+'\t'+str(line[0]+line[1])+'\t'+str(line[1])+'\t'+line[-1]+'\n')
    outfile.close()
def add_dic(somatic_dic,germline_dic):
    new_dic={}
    tmp=[]
    if germline_dic !={}:
        for key in somatic_dic:
            if key in germline_dic:
                tmp=somatic_dic[key]+germline_dic[key]
                new_dic[key]=sorted(tmp,key=lambda d:d[0],reverse=True)
            else:
                continue
        for key in germline_dic:
            if key not in somatic_dic:
                new_dic[key]=germline_dic
    else:
        new_dic=somatic_dic
    return new_dic
            

    
    
    
