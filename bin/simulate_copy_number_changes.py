#!/usr/bin/env python

from __future__ import print_function
from collections import defaultdict
import random, getopt, sys, pysam, copy, re, os, pprint, pdb, json
from math import *
from bx.intervals.intersection import IntervalTree
from subprocess import call

ScriptDir=os.path.abspath(os.path.dirname(__file__))



def main(argv):
    """reads in two haploid bam files so odd ploidy can be modelled"""

    (output_bam_file, bam_file_haplotype1, bam_file_haplotype2, query_chr,truth_file, mutate_somatic_variants, somatic_snv_file, somatic_indel_file, subsample_somatic_snvs, ploidy_depth,bam_depth) = read_in_arguments(argv)
    output_bam_file_haplotype1 = os.path.join(os.path.dirname(output_bam_file), 'haplotype1'+ os.path.basename(output_bam_file))
    output_bam_file_haplotype2 = os.path.join(os.path.dirname(output_bam_file), 'haplotype2'+ os.path.basename(output_bam_file))
  
    cl = re.search(".+" + query_chr + "_(.+).bam", output_bam_file).groups()[0]

    print(ploidy_depth)
    truth_set_cn_calls=read_in_truth_set(truth_file)

    snv_out_dir = os.path.dirname(output_bam_file) + "/forced_somatic_snv_frequencies"
    somatic_variants={}
    if mutate_somatic_variants ==True:
        somatic_variants=read_in_somatic_vcf_file(somatic_snv_file,somatic_indel_file,query_chr,cl,truth_set_cn_calls, snv_out_dir) #somatic_indel_file = '' if not defined and not used anyway at the mo


    create_synthetic_bam(bam_file_haplotype1, output_bam_file_haplotype1, query_chr, 'firsthaplotype', mutate_somatic_variants, somatic_variants, subsample_somatic_snvs, truth_set_cn_calls, ploidy_depth, bam_depth)
    print("finished first haplotype")

    create_synthetic_bam(bam_file_haplotype2, output_bam_file_haplotype2, query_chr, 'secondhaplotype', mutate_somatic_variants, somatic_variants, subsample_somatic_snvs, truth_set_cn_calls, ploidy_depth, bam_depth)
    print("finished second haplotype")

    print("concatenate bams ")
    samtools=os.path.join(ScriptDir,"../redist","samtools-1.2/samtools")
    headersam = os.path.join(ScriptDir , '../config/header.sam')
    cl=samtools + " cat -h " + headersam + " -o "+output_bam_file +" "+output_bam_file_haplotype1+ " "+ output_bam_file_haplotype2
    os.system(cl)

    print("deleting haplotype bams")
    rmh1='rm '+output_bam_file_haplotype1
    rmh2='rm '+output_bam_file_haplotype2
    os.system(rmh1)
    os.system(rmh2)



def read_in_arguments(argv):
    opts, args = getopt.getopt(argv,"a:b:o:c:t:m:s:i:d:e:",["i1file=","i2file=","ofile=","chr_no=","truth_file=","subsample_somatic_snvs=","depth=","bam_depth="])
    somatic_snv_file=''
    somatic_indel_file='' # not using for now will probably crash function below if used
    mutate_somatic_variants = False
    subsample_somatic_snvs=-1
    ploidy_depth = 20
    bam_depth = 174
    for opt, arg in opts:
        print("opt ",opt, " arg ",arg)
        if opt == '-h':
            print('simulate_copy_number_changes.py -a <first_haplotype_bam_file> -b <second_haplotype_bam_file> -o <output_bam_file> -c <chr_no> -t <truth_bed_file> -m <boolean mutate_somatic_variants> -d <depth> -e <bam_depth>')
            print("optional arguments -s (somatic snv file)")
            sys.exit()
        elif opt in ("-a", "--i1file"):
            bam_file_haplotype1 = arg
        elif opt in ("-b", "--i2file"):
            bam_file_haplotype2 = arg
        elif opt in ("-o", "--ofile"):
            output_bam_file = arg
        elif opt in ("-d", "--depth"):
            print("arg depth ",arg)
            ploidy_depth = float(arg)
        elif opt in ("-c", "--chr_no"):
            query_chr = arg
        elif opt in ("-t", "--truth_file"):
            truth_file = arg
        elif opt in ("-m", "--mutate_somatic_variants"):
            mutate_somatic_variants = bool(arg)
            print("mutate_somatic_variants ",mutate_somatic_variants)
        elif opt == '-s':
            somatic_snv_file=arg
        elif opt in ("-e", "--bam_depth"):
            print("bam depth ",arg)
            bam_depth=float(arg)
        elif opt == '-i':
            somatic_indel_file=arg
        elif opt == '--subsample_somatic_snvs':# sometimes we dont want to add in all somatic snvs eg to simulate low mutation samples
            subsample_somatic_snvs=float(arg)
            if subsample_somatic_snvs <0 or subsample_somatic_snvs>1:
                sys.exit('subsampling somatic snv value must be between 0 and 1')

    print("tfi "+truth_file)
    print("pc ss snvs",subsample_somatic_snvs)
    return(output_bam_file,bam_file_haplotype1,bam_file_haplotype2,query_chr,truth_file,mutate_somatic_variants,somatic_snv_file,somatic_indel_file,subsample_somatic_snvs,ploidy_depth,bam_depth)

def pileup(bamfile_name):
   samfile = pysam.Samfile(bamfile_name,"rb")
   output_pileup = open("test.cov", "w")
   sum_cov=0
   for pos_pileup in samfile.pileup():
       sum_cov+=pos_pileup.n
       if not pos_pileup.pos % 1e3:
          av_cov=int(sum_cov/1e3)
          sum_cov=0
          print(pos_pileup.tid,pos_pileup.pos,av_cov,file=output_pileup)


def read_in_truth_set(truth_file):
    truth_set_cn_calls={}
    truth_set=open(truth_file,'r')
    for line in truth_set:
        (chr,start,end,CN_hapA,CN_hapB)=line.strip().split("\t")
        (CN_hapA,CN_hapB)=(float(CN_hapA),float(CN_hapB)) 

        if not chr in truth_set_cn_calls:
            truth_set_cn_calls[chr]=[]
        truth_set_cn_calls[chr].append({'start':int(start),'end':int(end),'firsthaplotype_CN':CN_hapA,'secondhaplotype_CN':CN_hapB})

    print("ts ",truth_set_cn_calls)
    return truth_set_cn_calls

def create_synthetic_bam(bam_file,output_bam_file,query_chr, haplotype,mutate_somatic_variants,somatic_variants,subsample_somatic_snvs,truth_set_cn_calls,ploidy_depth,bam_depth):
    
    paired_read={}
    samfile = pysam.Samfile(bam_file,"rb")
    output_bam= pysam.Samfile(output_bam_file, "wb", template=samfile)


    for region in truth_set_cn_calls[query_chr]:
        print("qcr ",query_chr)
        h=haplotype+'_CN'
        subsample_fraction=float(ploidy_depth*region[h])/bam_depth
        if subsample_fraction > 1:
            sys.exit("Oversampling error: CN=%d, ploidy_depth=%d, bam_depth=%d, subsample_fraction=%f" % (region[h], ploidy_depth, bam_depth, subsample_fraction))

        print("region "+str(region),"pd ", str(ploidy_depth), " cn ", region[h], " ss ",subsample_fraction)
        c=cd=rl=xx=r2=0
        p=q={}
        chr='chr'+query_chr

        for read_line in samfile.fetch(query_chr,region['start'],region['end']):
            c+=1
            if not read_line.pos % 1e5:
                print("pos achieved ", read_line.pos)
            if read_line.pos<region['start']:
                continue # attempt to address subtle bug where reads are entered twice either side of breakpoint

            # IF FIRST READ ADD TO HASH TO BE WRITTEN TO FILE WITH SECOND READ
            if read_line.pos <read_line.pnext: #skipping reads with pair same start point for the mo
                rl+=1

                if random.random()<subsample_fraction:
                    if mutate_somatic_variants ==True:
                        if (subsample_somatic_snvs==-1 or random.random() <subsample_somatic_snvs):
                            somatic_variants_in_read=somatic_variants[query_chr].find(read_line.pos,read_line.pos+read_line.inferred_length)
                            read_line=add_in_somatic_variants(read_line,haplotype,somatic_variants_in_read) # CHANGE READS IN PLACE

                    if read_line.is_proper_pair:
                        if not read_line.pnext in paired_read:
                           paired_read[read_line.pnext]={}
                        paired_read[read_line.pnext][read_line.qname]=read_line

                    #WRITE OUT ORPHAN READS - need to check if gets secondary reads
                    else:
                        output_bam.write(read_line)
                        cd+=1

            # ELSE IF SECOND IN PAIR : WRITE TO FILE
            elif read_line.pos in paired_read:
                r2+=1
                if mutate_somatic_variants ==True:
                    somatic_variants_in_read=somatic_variants[query_chr].find(read_line.pos,read_line.pos+read_line.inferred_length)
                    read_line=add_in_somatic_variants(read_line,haplotype,somatic_variants_in_read) # CHANGE READS IN PLACE

                if read_line.qname in paired_read[read_line.pos]:
                    output_bam.write(read_line)
                    output_bam.write(paired_read[read_line.pos][read_line.qname])
                    cd+=1
                    d=paired_read[read_line.pos]
                    d.pop(read_line.qname,None)
                    if len(paired_read[read_line.pos])<1:
                        paired_read.pop(read_line.pos,None)
            else:
                xx+=1
            #    output_bam.write(read_line)

            #stop dict becoming too big by removing read ids (we know how far away the pair read is)
            if not read_line.pos % 1e4:
                tmpkeys=paired_read.keys()
                for pos in tmpkeys:
                   if pos<read_line.pos:
                       if len(paired_read[pos].keys())>0:
                           print("pnext ",pos, " pos ",read_line.pos)
                           print("missed read",paired_read[pos])
                           #sys.exit("missed read")
                           paired_read.pop(pos, None)

        print("rl ",rl,"c ",c, " cd ",cd, " xx",xx," r2 ",r2)
        print("end of region",region['end'])
        

    return()

def assign_freq_based_on_ploidy(copy_number):
    """
    copy number can now be non integer i.e. mixture of CN
    if CN=1.25 then == 3xhaploidCN and 1xdiploidCN
    model will always assume eg 2.5 is a mixture of 2&3 not say 1 and 4
    first is CN a fraction of two integers i.e. .5 .33 .25 .20 .1666
    """
    if copy_number <0.02:
        return(0)
    #if heterogenous round down to nearest copy
    CN_A=int(copy_number)
    
    freq=what_freq_for_haplotype_CN(CN_A)
    
    return(freq)

def what_freq_for_haplotype_CN(copy_number):
     freq=None
     if copy_number==0: freq=0
     if copy_number==1: freq=1
     elif copy_number==2:
            if random.random() >=0.5:
                freq=0.5
            else:
                freq=1
     elif copy_number==3:
        if random.random()<0.25:
            freq=0.33
        if random.random() >=0.25 and random.random<=0.75:
            freq=0.66
        if random.random>0.75:
            freq=1

     if freq :print("copyn ",copy_number, " freq ",freq)
     return(freq)

     
def read_in_somatic_vcf_file(somatic_snv_file, somatic_indel_file, query_chr,query_cl, truth_set_cn_calls, output_dir):
    print(os.path.join(output_dir,'forced_somatic_snv_frequencies_' + str(query_chr) + str(query_cl) + '.txt'))
    fsf=open(os.path.join(output_dir,'forced_somatic_snv_frequencies_' + str(query_chr) + str(query_cl) + '.txt'),'w')
    print("ri ",query_chr)
    h={}
    h2={}

    for file in [somatic_snv_file]:# for now just do SNVs as adding in indels involve increasing the size of reads which could cause issues; thinking about it probably wouldnt - quite faffy though
        FH=open(file,'r')
        for line in FH:
            if re.match('#',line):
                continue
            (chrom, pos, id, ref, alt, qual, filter, info, format, normal, tumor)=line.strip().split()
            pos=int(pos)
            if chrom != query_chr:
                continue

            if format != 'DP:FDP:SDP:SUBDP:AU:CU:GU:TU':
                sys.exit('vcf format not the usual'+'DP:FDP:SDP:SUBDP:AU:CU:GU:TU')
            print("tumor ",tumor)
            (DP,FDP,SDP,SUBDP,AU,CU,GU,TU)=tumor.strip().split(':')
            cov=float(DP)

            if ref =='A':
                l=[CU,GU,TU]
            if ref =='C':
                l=[AU,GU,TU]
            if ref =='G':
                l=[CU,AU,TU]
            if ref =='T':
                l=[CU,GU,AU] #should be a pithy python way to do this but this'll do for now


            (first,second,third)=sorted([int(cv.split(',')[0] ) for cv in l],reverse=True) #just using first tier reads for now

            if random.random()>0.5:
                somatic_haplotypeCN='firsthaplotype_CN'
                somatic_haplotype='firsthaplotype'
            else:
                somatic_haplotypeCN='secondhaplotype_CN'
                somatic_haplotype='secondhaplotype'

            print("pos ",pos, "shcn ", somatic_haplotypeCN , " r ")

            region_CN = 2
            for region in truth_set_cn_calls[query_chr]:
                if pos >=region['start'] and pos<=region['end']:
                    region_CN=region[somatic_haplotypeCN]
            somatic_mutation_freq=assign_freq_based_on_ploidy(region_CN)
            tree = IntervalTree()

            if not query_chr in h:
                h[query_chr] = tree

            h[query_chr].add(pos,pos,{'pos':pos,'ref':ref,'alt':alt,'line':line,'freq':somatic_mutation_freq,'somatic_haplotype':somatic_haplotype}) #theoretically bug: could have snv and indel at same pos #also bug if snp or indel last/first on read
            h2[pos]={'pos':pos,'ref':ref,'alt':alt,'line':line,'freq':somatic_mutation_freq,'somatic_haplotype':somatic_haplotype}

    pprint.pprint(h2)
    json.dump(h2,fsf,indent=4, sort_keys=True)
    return h

def add_in_somatic_variants(read_line,which_haplotype,somatic_variants_in_read): 
    """Add SNV variants
    NUmber of ploidy; at the mo max ploidy per haplotype is 3; but what are the prob of the somatic snvs being on 1,2 or 3 ploidy:
    if 1 ploidy 100%; if 2 ploidy 50/50; if 3 ploidy 0.33,0.66,0.166;
    """

    for d in somatic_variants_in_read:
        if d['somatic_haplotype'] =='both' or d['somatic_haplotype'] == which_haplotype:
            random_no=random.random()

            if d['somatic_haplotype'] == which_haplotype and random_no<d['freq']:

                pos_in_read=d['pos']-read_line.pos
                first_cigar_tuple=read_line.cigar[0]
                if first_cigar_tuple[0]==4 or first_cigar_tuple[0]==5: # hard or soft clipping
                    pos_in_read-=first_cigar_tuple[1]

                q=read_line.qual
                if re.search(',',d['alt']): continue # ignore tri-allelic positions for now

                #soft clipping not included with start of read
                upto=pos_in_read-1
                if upto ==-1: upto=0 #to fix bug
                read_line.seq=read_line.seq[:upto] +d['alt']  + read_line.seq[upto+1:]
                read_line.qual=q


    return(read_line) # should be mutated inline

if __name__ == "__main__":
    main(sys.argv[1:])
