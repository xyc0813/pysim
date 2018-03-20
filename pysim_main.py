from pysim1 import *
import random
import sys
from optparse import OptionParser
import getopt
import time
import os
import configparser
def main():
    usage = """%prog -c ini.config
Author: Yuchao Xia	
Description: simulate germline and somatic SVs.
	"""
    parser = OptionParser(usage)
    parser.add_option("-c", "--config", dest="config", help="config file",metavar="FILE")
    (opts, args) = parser.parse_args()
    if opts.config is None:
        parser.print_help()
    else:
        cf = configparser.ConfigParser()
        cf.read(opts.config)
        SV_config_file=cf.get("pysim_settings", "SV_config_file")
        ref_fasta=cf.get("pysim_settings", "ref_fasta")
        dbsnp=cf.get("pysim_settings", "dbsnp")
        somatic=cf.get("pysim_settings", "somatic")
        germline_num=cf.getint("pysim_settings", "germline_num")
        somatic_num=cf.getint("pysim_settings", "somatic_num")
        db=cf.get("pysim_settings", "db")
        somatic_SNP_db=cf.get("pysim_settings", "somatic_SNP_db")
        germ_homozyg_rate=cf.getfloat("pysim_settings",'germline_homozygous_rate')
        soma_homozyg_rate=cf.getfloat("pysim_settings",'somatic_homozygous_rate')
        up_down_stream=cf.getint("pysim_settings", "up_down_stream")
        snp_rate=cf.getfloat("pysim_settings",'snp_rate')
        indel_prob=cf.getfloat("pysim_settings",'indel_prob')
        min_indel_length=cf.getint("pysim_settings", "min_indel_length")
        max_indel_length=cf.getint("pysim_settings", "max_indel_length")
        sub_clone=cf.getint("pysim_settings", 'sub_clone')
        out_prex=cf.get("pysim_settings", "out_prex") # prefix for output files
        chrome_can=cf.get("pysim_settings", "chrome") # TODO OSC not documented
        germ_ratio=cf.get("pysim_settings", "germ_ratio")
        sub_sub_clone=cf.getint("pysim_settings",'sub_sub_clone') # TODO OSC not documented
        len_list=read_config(SV_config_file)
        ref=reference(ref_fasta)
        ref_range=compute_range(ref)
        snp_dic=read_dbsnp(dbsnp,chrome_can)
        #snp_dic=read_vcf(dbsnp,chrome_can)
        if somatic=='Y':
            # generate germline SNPs
            ref_result=generate_normal(ref,snp_dic,germline_num,out_prex+'germline_SNP_pois.txt',germ_homozyg_rate) # creates SNP_pois file
            ref_germline=ref_result[0] # {chr : sequence}
            ref_germline_new=ref_germline
            snp_list=ref_result[1] #{chr : indices changed}
            germline_dic={}
            output(ref_germline_new,germline_dic,out_prex+'germline.fa')
            for chrom in snp_list:
                print( ['snp_list_germline',len(snp_list[chrom])])
            
            # generate somatic genomes
            print('ref_germline',ref_germline)
            ref_somatic_result=generate_somatic(ref_germline,snp_dic,
                                                somatic_num,out_prex+'somatic_SNP_pois.txt',db,ref_range,soma_homozyg_rate) # creates SNP_pois file
            ref_somatic=ref_somatic_result[0]
            snp_list_somatic=ref_somatic_result[1]
            for chrom in snp_list_somatic:
                print( ['snp_list',len(snp_list_somatic[chrom])])
            print('ref_somatic',ref_somatic)
            dic=generate_pois(len_list,ref,ref_range)
            germline_dic={}

            # what is germ_ratio?
            if float(germ_ratio)>0:
                for key in dic:
                    tmp=random.sample(dic[key],len(dic[key])/2)
                    tmp_sort=sorted(tmp,key=lambda d:d[0],reverse=True)
                    germline_dic[key]=tmp_sort
                dic_remove={}
                dic_tmp=[]
                for key in dic:
                    for line in dic[key]:
                        if line not in germline_dic[key]:
                            dic_tmp.append(line)
                        else:
                            continue
                    dic_remove[key]=dic_tmp
                    dic_tmp=[]
                outfile=open(out_prex+'SV_germline.txt','w')
                ref_germline=generate_fasta(ref_germline,germline_dic,snp_rate,
                                            outfile,ref_range,up_down_stream,indel_prob,min_indel_length,max_indel_length)
                outfile.close()
                output(ref_germline_new,germline_dic,out_prex+'germline.fa')
            #else:
            #    output(ref_germline_new,germline_dic,out_prex+'germline.fa')

            outfile=open(out_prex+'SV_somatic.txt','w')
            ref_somatic=generate_fasta(ref_somatic,dic,snp_rate,outfile,ref_range,
                                       up_down_stream,indel_prob,min_indel_length,max_indel_length)
            outfile.close()
            output(ref_somatic,dic,out_prex+'somatic.fa')

            # generate subclone genomes, if desired.
            if sub_clone>=2:
                for i in range(sub_clone-1):
                    print('run subclone '+str(sub_clone-1))
                    sub_len_list=random.sample(len_list,len(len_list))
                    sub_dic=generate_subclone_pois(sub_len_list,ref,ref_range,germline_dic)
                    new_sub_dic=add_dic(sub_dic,germline_dic)
                    ref_somatic_sub_result=generate_somatic(ref_germline_new,snp_dic,somatic_num,
                                                            out_prex+'somatic_SNP_pois_subclone_'+str(i+1)+'.txt',
                                                            db,ref_range,soma_homozyg_rate)
                    ref_sub_somatic=ref_somatic_sub_result[0]
                    ref_sub_somatic_new=ref_sub_somatic
                    snp_list_sub_somatic=ref_somatic_sub_result[1]
                    snp_list_sub_somatic_new=snp_list_sub_somatic
                    outfile=open(out_prex+'SV_somatic_subclone_'+str(i+1)+'.txt','w')
                    ref_sub_somatic=generate_fasta(ref_sub_somatic,new_sub_dic,snp_rate,outfile,ref_range,
                                                   up_down_stream,indel_prob,min_indel_length,max_indel_length)
                    outfile.close()
                    output(ref_sub_somatic,new_sub_dic,out_prex+'somatic_subclone_'+str(i+1)+'.fa')
                    if sub_sub_clone>0:
                        for j in range(sub_sub_clone):
                            sub_len_list=random.sample(len_list,len(len_list))
                            sub_dic=generate_subclone_pois(sub_len_list,ref,ref_range,germline_dic)
                            new_sub_dic_1=add_dic(sub_dic,new_sub_dic)
                            ref_somatic_sub_result=generate_somatic(ref_somatic_new,snp_dic,
                                                                    somatic_num,out_prex+'somatic_SNP_pois_sub_subclone_'+str(j+1)+'.txt',
                                                                    db,ref_range,soma_homozyg_rate)
                            ref_sub_somatic_new=ref_somatic_sub_result[0]
                            snp_list_sub_somatic_new=ref_somatic_sub_result[1]
                            outfile=open(out_prex+'SV_somatic_sub_subclone_'+str(i+1)+'.txt','w')
                            ref_sub_somatic=generate_fasta(ref_sub_somatic_new,new_sub_dic_1,snp_rate,outfile,ref_range,
                                                           up_down_stream,indel_prob,min_indel_length,max_indel_length)
                            outfile.close()
                            output(ref_sub_somatic,new_sub_dic,out_prex+'somatic_sub_subclone_'+str(i+1)+'.fa')
                        
                
                
        else:
            ref_result=generate_normal(ref,snp_dic,germline_num,out_prex+'germline_SNP_pois.txt',germ_homozyg_rate)
            ref_germline=ref_result[0]
            snp_list=ref_result[1]
            dic=generate_pois(len_list,ref,ref_range)
            germline_dic=dic
            outfile=open(out_prex+'SV_germline.txt','w')
            ref_germline=generate_fasta(ref_germline,germline_dic,snp_rate,outfile,
                                        ref_range,up_down_stream,indel_prob,min_indel_length,max_indel_length)
            print('ref_germline',ref_germline)
            outfile.close()
            output(ref_germline,germline_dic,out_prex+'germline.fa')
if __name__ == "__main__":
    start = time.clock()
    print('simulation begins at:'+str(start))
    random.seed(42)
    main()
    end = time.clock()
    print('simulation ends at:'+str(end))
    print("The function run time is : %.03f seconds" %(end-start))         
        
        
        
        
        
