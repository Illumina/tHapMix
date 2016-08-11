#!/usr/bin/env python

from __future__ import print_function
from collections import OrderedDict
from collections import deque
import random, getopt, sys, copy, re, os, pprint, pdb, json
ScriptDir=os.path.join(os.path.abspath(os.path.dirname(__file__)),"..")
sys.path.append(ScriptDir)
import pysam
from math import *
from bx.intervals.intersection import IntervalTree
from subprocess import call
from itertools import izip, izip_longest


def main(argv):
    """reads in two haploid bam files so odd ploidy can be modelled"""
    
    # get input params
    options = get_commandline_options()
    output_bam_file = options.output_bam_file
    bam_file_haplotype1 = options.bam_file_haplotype1
    bam_file_haplotype2 = options.bam_file_haplotype2
    query_chr = options.query_chr
    truth_file = options.truth_file
    mutate_somatic_variants = options.mutate_somatic_variants
    somatic_snv_files = options.somatic_snv_files
    clonal_percs = options.clonal_percs
    subsample_somatic_snvs = options.subsample_somatic_snvs
    ploidy_depth = options.ploidy_depth
    bam_depth = options.bam_depth
    truth_set_cn_calls = read_in_truth_set(truth_file)
    snv_out_dir = os.path.dirname(output_bam_file) + "/forced_somatic_snv_frequencies"
    somatic_variants = IntervalTree()
    if options.rand_seed:
        random.seed(options.rand_seed)
    if mutate_somatic_variants:
        somatic_variants = read_in_somatic_vcf_file(somatic_snv_files, clonal_percs, query_chr, truth_set_cn_calls, snv_out_dir, subsample_somatic_snvs) #somatic_indel_file = '' if not defined and not used anyway at the mo
    create_synthetic_bam(bam_file_haplotype1, bam_file_haplotype2, output_bam_file, query_chr, mutate_somatic_variants, somatic_variants,  truth_set_cn_calls, ploidy_depth, bam_depth)
    print("finished first haplotype")


def get_commandline_options():
    """Get user options from command line"""

    import argparse

    parser = argparse.ArgumentParser(prog='HapMix v0.6')
    parser.add_argument("-a", "--i1file", dest="bam_file_haplotype1", help="haplotype bam file 1")
    parser.add_argument("-b", "--i2file", dest="bam_file_haplotype2", help="haplotype bam file 2")
    parser.add_argument("-c", "--chr_no", dest="query_chr", help="chromosome number")
    parser.add_argument("-o", "--ofile", dest="output_bam_file", help="name of the output bam_file")
    parser.add_argument("-m", "--mutate_somatic_variants", action='store_true', help="use clonal somatic SNV files?")
    parser.add_argument("-t", "--truth_file", dest="truth_file", help="name of the CNV truth file")
    parser.add_argument("-e", "--bam_depth", dest="bam_depth", help="Depth of bam file", type=float)
    parser.add_argument("-d", "--ploidy_depth", dest="ploidy_depth", help="Simulated depth per haplotype", type=float)
    parser.add_argument("-s", "--somatic_snv_files", nargs="*", help="clonal somatic SNV truth file")
    parser.add_argument("--seed", dest="rand_seed", help="Seed for random number generation", default=None, type=int)
    parser.add_argument("-p", "--clonal_percs", nargs="*", help="clonal percentages")
    parser.add_argument("--subsample_somatic_snvs", dest="subsample_somatic_snvs", help="Subsample percentage for clonal somatic SNV truth file", type=float, default = 1)    
    options = parser.parse_args()

    if not options.bam_file_haplotype1 or not os.path.exists(options.bam_file_haplotype1):
        print("\nHaplotype bam file 1 is not specified or does not exist\n")
        parser.print_help()
        sys.exit(2)
    if not options.bam_file_haplotype2 or not os.path.exists(options.bam_file_haplotype2):
        print("\nHaplotype bam file 2 is not specified or does not exist\n")
        parser.print_help()
        sys.exit(2)
    if not options.query_chr:
        print("\nChromosome number is not specified\n")    
        parser.print_help()
        sys.exit(2)
    if not options.output_bam_file:
        print("\nOutput bam_file is not specified\n")    
        parser.print_help()
        sys.exit(2)   
    if not options.bam_depth:
        print("\nBam depth is not specified\n")    
        parser.print_help()
        sys.exit(2)   
    if not options.ploidy_depth:
        print("\nPloidy depth is not specified\n")    
        parser.print_help()
        sys.exit(2)           
    if  (not options.somatic_snv_files and options.clonal_percs) or (options.somatic_snv_files and not options.clonal_percs):
        sys.exit('Both somatic_snv_files and clonal_percs options should be specified')   
    if  options.somatic_snv_files and options.clonal_percs:
        if len(options.somatic_snv_files) != len(options.clonal_percs):
            sys.exit('somatic_snv_files and clonal_percs options should have the same length')     
    if options.subsample_somatic_snvs:
        subsample_somatic_snvs = options.subsample_somatic_snvs
        if subsample_somatic_snvs <0 or subsample_somatic_snvs>1:
            sys.exit('subsampling somatic snv value must be between 0 and 1')

    return options
   


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
        (chr, start, end, CN_hapA, CN_hapB) = line.strip().split("\t")
        (CN_hapA, CN_hapB) = (float(CN_hapA), float(CN_hapB)) 
        if not chr in truth_set_cn_calls:
            truth_set_cn_calls[chr] = []
        truth_set_cn_calls[chr].append({'start':int(start), 'end':int(end), 'firsthaplotype_CN':CN_hapA, 'secondhaplotype_CN':CN_hapB})
    print("ts ", truth_set_cn_calls)
    return truth_set_cn_calls
    
    
def write_eligible_read(read_line, region, paired_read, subsample_fraction, mutate_somatic_variants, output_bam, somatic_variants, haplotype):
    # debugging variables
    global cd, r1, xx, r2
    if read_line.pos < region['start']:
        xx+=1
        return # attempt to address subtle bug where reads are entered twice either side of breakpoint
     
    # IF FIRST READ ADD TO HASH TO BE WRITTEN TO FILE 
    if read_line.pos < read_line.pnext: #skipping reads with pair same start point for the mo
        if random.random() < subsample_fraction:
            if mutate_somatic_variants:
                somatic_variants_in_read = somatic_variants.find(read_line.pos,read_line.pos+read_line.inferred_length)
                read_line = add_in_somatic_variants(read_line,haplotype,somatic_variants_in_read) # CHANGE READS IN PLACE
            if read_line.is_proper_pair:
                r1 += 1
                if not read_line.qname in paired_read:
                    paired_read.add(read_line.qname)
                output_bam.write(read_line)
               
            # WRITE OUT ORPHAN READS - need to check if gets secondary reads
            else:
                output_bam.write(read_line)
                #print ("read pos "+ str(read_line.pos))
                cd += 1
        else:
            xx +=1        

    # ELSE IF SECOND IN PAIR : WRITE TO FILE
    elif read_line.qname in paired_read:
        if mutate_somatic_variants:
            somatic_variants_in_read = somatic_variants.find(read_line.pos, read_line.pos + read_line.inferred_length)
            read_line = add_in_somatic_variants(read_line, haplotype, somatic_variants_in_read) # CHANGE READS IN PLACE
        r2 += 1
        output_bam.write(read_line)
        #print ("read pos " + str(read_line.pos))
        paired_read.remove(read_line.qname)
    else:
        xx += 1
        


def write_eligible_reads(read_container, region, paired_read,  mutate_somatic_variants, output_bam, min_element, somatic_variants):
    first_occurence = True
    global cd, r1, xx, r2
    
    while read_container:
        read_line = read_container.popleft()
        if read_line.pos < region['start']:
            continue
            xx+=1           
        if read_line.pos > min_element:
            read_container.appendleft(read_line)
            break
         
        # IF FIRST READ ADD TO HASH TO BE WRITTEN TO FILE 
        if read_line.pos < read_line.pnext: #skipping reads with pair same start point for the mo
            r1 += 1

            if mutate_somatic_variants:
                if random.random() > 0.5:
                    haplotype = 'firsthaplotype_CN'
                else:
                    haplotype = 'secondhaplotype_CN'
                somatic_variants_in_read = somatic_variants.find(read_line.pos, read_line.pos + read_line.inferred_length)
                read_line = add_in_somatic_variants(read_line, haplotype, somatic_variants_in_read) # CHANGE READS IN PLACE

            if read_line.is_proper_pair:
                output_bam.write(read_line)

            # WRITE OUT ORPHAN READS - need to check if gets secondary reads
            else:
                output_bam.write(read_line)
                cd += 1

        # ELSE IF SECOND IN PAIR : WRITE TO FILE
        elif read_line.qname in paired_read:
            if mutate_somatic_variants:
                if random.random() > 0.5:
                    haplotype = 'firsthaplotype_CN'
                else:
                    haplotype = 'secondhaplotype_CN'
                somatic_variants_in_read = somatic_variants.find(read_line.pos, read_line.pos + read_line.inferred_length)
                read_line = add_in_somatic_variants(read_line, haplotype, somatic_variants_in_read) # CHANGE READS IN PLACE
            r2 += 1
            output_bam.write(read_line)
        #print ("read pos " + str(read_line.pos))
            paired_read.remove(read_line.qname)
        else:
            xx += 1
    

def create_synthetic_bam(bam_file_haplotype1, bam_file_haplotype2, output_bam_file, query_chr, mutate_somatic_variants, somatic_variants, truth_set_cn_calls, ploidy_depth, bam_depth):

    # paired_read holds read names for the first pair
    paired_read = set()
    samfile_haplotype1 = pysam.Samfile(bam_file_haplotype1,"rb")    
    samfile_haplotype2 = pysam.Samfile(bam_file_haplotype2,"rb")
    output_bam = pysam.Samfile(output_bam_file, "wb", template = samfile_haplotype1)

    # subsample each region from the bed file
    for region in truth_set_cn_calls[query_chr]:
        print("query_chr ", query_chr)
        hs = ['firsthaplotype_CN', 'secondhaplotype_CN']
        subsample_fractions = []
        
        for h in hs:
            subsample_fraction = float(ploidy_depth*region[h]) / bam_depth
            if subsample_fraction > 1:
                sys.exit("Oversampling error: CN=%d, ploidy_depth=%d, bam_depth=%d, subsample_fraction=%f" % (region[h], ploidy_depth, bam_depth, subsample_fraction))
            print("region " + str(region), "ploidy_depth ", str(ploidy_depth), " copy_number ", region[h], " subsample_fraction ", subsample_fraction)
            subsample_fractions.append(subsample_fraction)
            
        global cd, r1, xx, r2 # debug variables, r1/r2 : first and second read pair, cd - orphan reads, xx - the rest
        cd = r1 = xx = r2 = 0
        min_element = 0 # smaller read pos when reading in two bams simultaneously
        reserve_min_element = sys.maxint # smaller read pos of reads saved from previous readline queries
        reserve_read_container = deque()
        chr = 'chr' + query_chr

        # to cover shorter BAM with dummy reads
        dummy_read = pysam.AlignedSegment()
        dummy_read.pos = -1     
        # iterate over two bam files 
        for iterator in izip_longest(samfile_haplotype1.fetch(query_chr,region['start'],region['end']),samfile_haplotype2.fetch(query_chr,region['start'],region['end']), fillvalue = dummy_read):
            read_lines = list(iterator)
            # case1: first bam shorter   
            if read_lines[0].pos == -1:
                min_element = read_lines[1].pos
                write_eligible_read(read_lines[0], region, paired_read, subsample_fractions[0],  mutate_somatic_variants, output_bam, somatic_variants, hs[0])
            # case2: second bam shorter                    
            elif read_lines[1].pos == -1:
                min_element = read_lines[0].pos
                write_eligible_read(read_lines[1], region, paired_read, subsample_fractions[1],  mutate_somatic_variants, output_bam, somatic_variants, hs[1])
            # case3: reads from both bams
            else:
                min_element = min(read_lines[1].pos, read_lines[0].pos)
                # all reserve read positions are larger than currently processed read
                if reserve_min_element >= min_element:
                    if not reserve_read_container:
                        reserve_min_element = max(read_lines[1].pos, read_lines[0].pos)

                    if read_lines[1].pos == read_lines[0].pos:  
                        write_eligible_read(read_lines[0], region, paired_read, subsample_fractions[0],  mutate_somatic_variants, output_bam, somatic_variants, hs[0])
                        write_eligible_read(read_lines[1], region, paired_read, subsample_fractions[1],  mutate_somatic_variants, output_bam, somatic_variants, hs[1])  

                    elif read_lines[1].pos > read_lines[0].pos:
                        write_eligible_read(read_lines[0], region, paired_read, subsample_fractions[0],  mutate_somatic_variants, output_bam, somatic_variants, hs[0])
                        if random.random() < subsample_fractions[1]:
                            reserve_read_container.append(read_lines[1])
                            if read_lines[1].is_proper_pair:
                                if not read_lines[1].qname in paired_read:
                                    paired_read.add(read_lines[1].qname)
                        else:
                            xx += 1                        
                    else:
                        write_eligible_read(read_lines[1], region, paired_read, subsample_fractions[1],  mutate_somatic_variants, output_bam, somatic_variants, hs[1])
                        if random.random() < subsample_fractions[0]:
                            reserve_read_container.append(read_lines[0])
                            if read_lines[0].is_proper_pair:
                                if not read_lines[0].pnext in paired_read:
                                    paired_read.add(read_lines[0].qname)
                        
                        else:
                            xx += 1                                     
                # some reserve read positions are smaller than currently processed read                        
                else:
                    # write reserve reads with smaller read positions first
                    write_eligible_reads(reserve_read_container, region, paired_read,  mutate_somatic_variants, output_bam, min_element, somatic_variants)
                    # no reserve reads left
                    if not reserve_read_container:                    
                        reserve_min_element = max(read_lines[1].pos, read_lines[0].pos)
                    else:
                        reserve_min_element_tmp = reserve_read_container.popleft()
                        reserve_min_element = reserve_min_element_tmp.pos
                        reserve_read_container.appendleft(reserve_min_element_tmp)
                    
                    if read_lines[1].pos == read_lines[0].pos:  
                        write_eligible_read(read_lines[0], region, paired_read, subsample_fractions[0],  mutate_somatic_variants, output_bam, somatic_variants, hs[0])
                        write_eligible_read(read_lines[1], region, paired_read, subsample_fractions[1],  mutate_somatic_variants, output_bam, somatic_variants, hs[1])
                    elif read_lines[1].pos > read_lines[0].pos:
                        write_eligible_read(read_lines[0], region, paired_read, subsample_fractions[0],  mutate_somatic_variants, output_bam, somatic_variants, hs[0])
                        if random.random() < subsample_fractions[1]:
                            reserve_read_container.append(read_lines[1])
                            if read_lines[1].is_proper_pair:
                                if not read_lines[1].pnext in paired_read:
                                    paired_read.add(read_lines[1].qname)                        
                        else:
                            xx += 1                                     
                    else:
                        write_eligible_read(read_lines[1], region, paired_read, subsample_fractions[1],  mutate_somatic_variants, output_bam, somatic_variants, hs[1])
                        if random.random() < subsample_fractions[0]:
                            reserve_read_container.append(read_lines[0])
                            if read_lines[0].is_proper_pair:
                                if not read_lines[0].pnext in paired_read:
                                    paired_read.add(read_lines[0].qname)                        
                        else:
                            xx += 1                                     
                

            #if not read_lines[0].pos % 1e5:
            #    print("pos achieved ", read_lines[0].pos)
            #    print("itertime ", time.time())
            #    print("reads ", c)

        print("r1 ", r1, " r2 ", r2, " cd ", cd, " xx", xx)
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
     if copy_number== 0: freq=0
     if copy_number == 1: freq=1
     elif copy_number == 2:
            if random.random() >=0.5:
                freq=0.5
            else:
                freq=1
     elif copy_number >= 3:
        if random.random()<0.25:
            freq=0.33
        if random.random() >=0.25 and random.random <= 0.75:
            freq=0.66
        if random.random>0.75:
            freq=1

     if freq :print("copyn ",copy_number, " freq ",freq)
     return(freq)

     
def read_in_somatic_vcf_file(somatic_snv_files, clonal_percs, query_chr, truth_set_cn_calls, output_dir, subsample_somatic_snvs):
    """ read clonal somatic SNV vcf files"""
    fsf = open(os.path.join(output_dir,'forced_somatic_snv_frequencies_' + str(query_chr) + '.json'), 'w')
    print("ri ",query_chr)
    
    h = IntervalTree()
    h2={}

    for (somatic_snv_file, clonal_perc) in zip(somatic_snv_files, clonal_percs):
    # for now just do SNVs as adding in indels involve increasing the size of reads which could cause issues; 
    # thinking about it probably wouldnt - quite faffy though
        FH = open(somatic_snv_file,'r')
        for line in FH:
            if re.match('#',line):
                continue
            random_no = random.random()
            if random_no > subsample_somatic_snvs: 
                continue            
            (chrom, pos, id, ref, alt, qual, filter, info, format, normal, tumor)=line.strip().split()
            pos=int(pos)
            if chrom != query_chr:
                continue

            if format != 'DP:FDP:SDP:SUBDP:AU:CU:GU:TU':
                sys.exit('vcf format not the usual'+'DP:FDP:SDP:SUBDP:AU:CU:GU:TU')
            print("tumor ",tumor)
            (DP,FDP,SDP,SUBDP,AU,CU,GU,TU) = tumor.strip().split(':')
            cov=float(DP)

            if ref =='A':
                l=[CU,GU,TU]
            if ref =='C':
                l=[AU,GU,TU]
            if ref =='G':
                l=[CU,AU,TU]
            if ref =='T':
                l=[CU,GU,AU] #should be a pithy python way to do this but this'll do for now

            (first, second, third)=sorted([int(cv.split(',')[0] ) for cv in l], reverse=True) #just using first tier reads for now

            if random.random() > 0.5:
                somatic_haplotypeCN = 'firsthaplotype_CN'
            else:
                somatic_haplotypeCN = 'secondhaplotype_CN'

            print("pos ",pos, "shcn ", somatic_haplotypeCN , " r ")

            region_CN = 2
            for region in truth_set_cn_calls[query_chr]:
                if pos >= region['start'] and pos <= region['end']:
                    region_CN = region[somatic_haplotypeCN]
            somatic_mutation_freq = float(assign_freq_based_on_ploidy(region_CN))
            somatic_mutation_freq *= float(clonal_perc)
            h.add(pos, pos,{'pos':pos, 'ref':ref, 'alt':alt, 'line':line, 'freq':somatic_mutation_freq, 'somatic_haplotype':somatic_haplotypeCN}) #theoretically bug: could have snv and indel at same pos #also bug if snp or indel last/first on read
            h2[pos] = {'pos':pos, 'ref':ref, 'alt':alt, 'line':line, 'freq':somatic_mutation_freq, 'somatic_haplotype':somatic_haplotypeCN}

    pprint.pprint(h2)
    json.dump(h2, fsf, indent=4, sort_keys=True)
    return h

def add_in_somatic_variants(read_line,which_haplotype,somatic_variants_in_read): 
    """Add SNV variants
    Number of ploidy; at the mo max ploidy per haplotype is 3; but what are the prob of the somatic snvs being on 1,2 or 3 ploidy:
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
