import random,sys,pysam,re,subprocess,HTSeq,pdb,argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('-s','--sample_name',help="please give sample name e.g. NA12878",required=True)
parser.add_argument('-c','--chr_no',help="please give chr_no",required=True)
parser.add_argument('-b','--bam_file',help="please specify bam file",required=True)
parser.add_argument('-r','--ref_file',help="please specify reference directory",required=True)
parser.add_argument('-n','--snp_file',help="please specify tabixed haplotyped SNP file",required=True)
parser.add_argument('-i','--indel_file',help="please specify tabixed haplotyped indel file",required=True)
parser.add_argument('-x','--chr_prefix',help="does the chromsome need a prefix eg chr",required=False)

args = parser.parse_args()
sample_name=args.sample_name
chr_no=args.chr_no
bam_file=args.bam_file
ref_file=args.ref_file
snp_file=args.snp_file
indel_file=args.indel_file
if (args.chr_prefix):
    chr= args.chr_prefix+str(chr_no)
else:
    chr=str(chr_no)

sequence={}
for s in HTSeq.FastaReader(ref_file):
    sequence[s.name]=s
reference_seq=sequence["chr"+str(chr_no)]
pos_ref=0
samfile = pysam.Samfile(bam_file,"rb")

haplotyped_snp_file=subprocess.Popen(['tabix',snp_file,chr_no ],stdout=subprocess.PIPE)
haplotyped_indel_file=subprocess.Popen(['tabix',indel_file,chr_no ],stdout=subprocess.PIPE)


#d={'hc':0,'hd':0,'bt':0,'ot':0,'rf':0,'fr':0}

haplotypeC_bam= pysam.Samfile("haplotypeC_"+chr +".bam", "wb", template=samfile)
haplotypeD_bam= pysam.Samfile("haplotypeD_"+chr+".bam", "wb", template=samfile)
haplotype_count={}

def main():
    read_variant_dict={}
    paired_read={}

    (haplotype_dict_snvs,haplotype_dict_snvs_pos)=read_in_vcf(haplotyped_snp_file)
    (haplotype_dict_indels,haplotype_dict_indels_pos)=read_in_vcf(haplotyped_indel_file)

    chr_variant_dict={}
    chr_variant_dict['haplotypeC']=dict(haplotype_dict_snvs['haplotypeC'].items()+haplotype_dict_indels['haplotypeC'].items())
    chr_variant_dict['haplotypeD']=dict(haplotype_dict_snvs['haplotypeD'].items()+haplotype_dict_indels['haplotypeD'].items())

    haplotype_dict_pos=dict(haplotype_dict_snvs_pos.items()+haplotype_dict_indels_pos.items())

    for read_line in samfile.fetch(chr):

        if read_line.cigar == None:
            continue #SKIPPING READ AS UNMAPPED

        if not read_line.qname in read_variant_dict:
            read_variant_dict[read_line.qname]={}
    
        rvd=variant_count(read_line,haplotype_dict_pos)
        read_variant_dict[read_line.qname].update(rvd) #HYPOTHETICAL BUG IF INDEL AND SNP AT SAME POS

        if not read_line.qname in haplotype_count:
            haplotype_count[read_line.qname]={'other':{},'C':{},'D':{}}

        #COUNT NUMBER OF MUTATIONS FOR EACH READ WHICH CAN BE ASSIGNED TO EACH HAPLOTYPE
        for variant_pos in read_variant_dict[read_line.qname].keys():

            if variant_pos in chr_variant_dict['haplotypeC'] and  variant_pos in chr_variant_dict['haplotypeD'] and read_variant_dict[read_line.qname][variant_pos]['call']== chr_variant_dict['haplotypeC'][variant_pos]['alt'] and read_variant_dict[read_line.qname][variant_pos]['call']== chr_variant_dict['haplotypeD'][variant_pos]['alt']: #check hom/het and call:
                haplotype_count[read_line.qname]['C'][variant_pos]={}
                haplotype_count[read_line.qname]['D'][variant_pos]={}

            elif variant_pos in chr_variant_dict['haplotypeC'] and  read_variant_dict[read_line.qname][variant_pos]['call']== chr_variant_dict['haplotypeC'][variant_pos]['alt']:
                haplotype_count[read_line.qname]['C'][variant_pos]={'call':read_variant_dict[read_line.qname][variant_pos]['call']}

            elif variant_pos in chr_variant_dict['haplotypeD'] and read_variant_dict[read_line.qname][variant_pos]['call']== chr_variant_dict['haplotypeD'][variant_pos]['alt']:
                haplotype_count[read_line.qname]['D'][variant_pos]={}

            else:
                haplotype_count[read_line.qname]['other'][variant_pos]={}

        # IS IT THE SECOND/ORPHAN READ? CAN THE READ BE ASSIGNED UNIQUELY TO EITHER OF THE HAPLOTYPES?
        if not read_line.is_proper_pair or  (read_line.pnext in paired_read  and read_line.qname in paired_read[read_line.pnext])  :

            haplotype=assign_to_haplotype(haplotype_count,paired_read,read_line)
            write_to_bam_file(haplotype,paired_read,read_line)
            haplotype_count.pop(read_line.qname, None)
            read_variant_dict.pop(read_line.qname, None)

        # IS IT THE FIRST READ? ADD TO DICT
        if read_line.is_proper_pair and not read_line.pnext in paired_read:
            if not read_line.pos in paired_read:
                paired_read[read_line.pos]={}
            if  not read_line.qname in paired_read[read_line.pos]:
                paired_read[read_line.pos][read_line.qname]=read_line

        #FLUSH DICTIONARIES EVERY 10k bp
        if not read_line.pos % 1e4:
            tmpkeys=paired_read.keys()
            for pos in tmpkeys:
                if pos<read_line.pos:
                    paired_read.pop(pos, None)

        




def read_in_vcf(vcf_file):
    cd={'haplotypeC':{},'haplotypeD':{}}
    csdl={}

    for line in vcf_file.stdout:
        if re.match('#',line):
            continue
        if not re.search('bwa',line) and not re.search('isaac',line): # ONLY TRUST ISAAC & BWA BASED CALLS
            continue
        else:
            
            
            (chrom,pos,id,ref,alt,qual,filter,info,format,NA12877,NA12878,NA12879,NA12880,NA12881,NA12882,NA12883,NA12884,NA12885,NA12886,NA12887,NA12888,NA12889,NA12890,NA12891,NA12892,NA12893)=line.strip().split('\t')
            
            if re.match('chr',chr) and not re.match('chr',chrom):
                chrom='chr'+chrom
            if chrom != chr:
                continue
            pos=int(float(pos))
            format_columns=format.split(':') #JUST GENOTYPE AND EDIT DISTANCE
            format_columns_data=eval(sample_name).split(':')
            f_dict={}
            for i,k in enumerate(format_columns):
               f_dict[k]=format_columns_data[i]


            if 'GT' in f_dict:
                if re.search("/",f_dict['GT']):
                    continue
                (ploidyC,ploidyD)=f_dict['GT'].split('|')
                (ploidyC,ploidyD)=(int(ploidyC),int(ploidyD))
                ploidyC_base_call=''
                ploidyD_base_call=''
                if ploidyC ==0 and ploidyD ==0:
                    continue # not haplotyped so skip
                if ploidyC ==0:
                    ploidyC_base_call=ref
                elif ploidyC ==1:
                    ploidyC_base_call=alt
                if ploidyD ==0:
                    ploidyD_base_call=ref
                elif ploidyD ==1:
                    ploidyD_base_call=alt

                if len(ref)==len(alt)==1:
                    type='S'
                if len(ref)==len(alt)!=1:
                    type='SUB'
                if len(ref)>len(alt):
                    type='D'
                if len(ref)<len(alt):
                    type='I'

                cd['haplotypeC'][pos]={'pos':pos,'alt':ploidyC_base_call}
                cd['haplotypeD'][pos]={'pos':pos,'alt':ploidyD_base_call}
                csdl[pos]={'ref':ref,'alt':alt,'type':type}


            else:
                sys.exit("no genotyping on line")
    return(cd,csdl)




def variant_count(read_line,haplotype_dict_pos):
    pos_in_read=0
    pos_ref=read_line.pos
    read_variant_dict={}
    for cigar_operations in read_line.cigar:
        (type_cigar_op,length_cigar_op)=cigar_operations

        if type_cigar_op==0 or type_cigar_op==7: #MATCH
            ref_pos=pos_ref
            for ii in range(0,length_cigar_op):
                chr='chr'+str(read_line.tid)
                ref_base=reference_seq.seq[ref_pos].upper()
                pos_ref+=1
                if pos_ref in haplotype_dict_pos: # IF ITS A HAPLOTYPED READ
                    if  haplotype_dict_pos[pos_ref]['type']=='S':
                        read_variant_dict[pos_ref]={'type':haplotype_dict_pos[pos_ref]['type'],'call':read_line.seq[pos_in_read],'ref':ref_base}
                    if  haplotype_dict_pos[pos_ref]['type']=='D':
                        ref_del=reference_seq.seq[pos_ref-1:pos_ref+length_cigar_op].upper()
                        read_variant_dict[pos_ref]={'type':'D','alt':haplotype_dict_pos[pos_ref]['alt'],'call':haplotype_dict_pos[pos_ref]['ref'],'ln':len(haplotype_dict_pos[pos_ref]['alt'])} # deletions vcf ref will be longer than alt
                    if  haplotype_dict_pos[pos_ref]['type']=='I':
                        read_variant_dict[pos_ref]={'type':'I','alt':haplotype_dict_pos[pos_ref]['alt'],'call':haplotype_dict_pos[pos_ref]['ref']} # for indels this has to be base before as well
                pos_in_read+=1

        elif type_cigar_op==3 : #N
            pos_in_read+=length_cigar_op
            pos_ref+=length_cigar_op
        elif type_cigar_op==4: # SOFT CLIP
            pos_in_read+=length_cigar_op #BAM FILE START POS IS AFTER SOFT CLIPPING
        elif type_cigar_op==1 :# INSERTION
            if pos_ref in haplotype_dict_pos:
                read_variant_dict[pos_ref]={'type':'I','call':read_line.seq[pos_in_read-1:pos_in_read+length_cigar_op],'ref':read_line.seq[pos_in_read-1]} # for indels this has to be base before as well
            pos_in_read+=length_cigar_op
            pos_ref+=1
        elif type_cigar_op==2 :# DELETION
            if pos_ref in haplotype_dict_pos:
                ref_del=reference_seq.seq[pos_ref-1:pos_ref+length_cigar_op].upper()
                read_variant_dict[pos_ref]={'type':'D','call':read_line.seq[pos_in_read-1],'alt':read_line.seq[pos_in_read-1],'ref':ref_del,'ln':length_cigar_op} # deletions vcf ref will be longer than alt
            pos_ref+=length_cigar_op
    return read_variant_dict



def write_to_bam_file(haplotype,paired_read,read_line):
    if haplotype =='haplotypeC':
        haplotypeC_bam.write(read_line)
    elif haplotype =='haplotypeD':
        haplotypeD_bam.write(read_line)

    if read_line.is_proper_pair:
        other_read=paired_read[read_line.pnext][read_line.qname]
        if haplotype =='haplotypeC':
            haplotypeC_bam.write(other_read)
        elif haplotype =='haplotypeD':
            haplotypeD_bam.write(other_read)



def assign_to_haplotype(haplotype_count,paired_read,read_line):

    if len(haplotype_count[read_line.qname]['C']) != 0 and len(haplotype_count[read_line.qname]['D']) == 0 :
        haplotype='haplotypeC'

    if len(haplotype_count[read_line.qname]['C']) == 0 and len(haplotype_count[read_line.qname]['D']) != 0 :
        haplotype='haplotypeD'

    elif len(haplotype_count[read_line.qname]['C']) != 0 and len(haplotype_count[read_line.qname]['D']) != 0 :
        if random.random()<0.5:
            haplotype='haplotypeC'
        else:
            haplotype='haplotypeD'

    elif len(haplotype_count[read_line.qname]['C']) == 0 and len(haplotype_count[read_line.qname]['D']) == 0  and len(haplotype_count[read_line.qname]['other']) != 0:
        if random.random()<0.5:
            haplotype='haplotypeC'
        else:
            haplotype='haplotypeD'

    elif len(haplotype_count[read_line.qname]['C']) == 0 and len(haplotype_count[read_line.qname]['D']) == 0  and len(haplotype_count[read_line.qname]['other']) == 0:
        if random.random() <0.5:
            haplotype='haplotypeC'
        else:
            haplotype='haplotypeD'

    return haplotype

if __name__ == "__main__":
    main()

