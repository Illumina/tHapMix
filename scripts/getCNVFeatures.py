import subprocess
import os
import glob
import json
import argparse
import sys
import numpy

ScriptDir=os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(ScriptDir, "sskernel.py"))
from sskernel import sskernel

bcftools = "bcftools "

def parse_options():

    parser = argparse.ArgumentParser(description="Get CNV properties from somatic vcfs")
    parser.add_argument("-o", "--output_dir", help="Name of the output directory", required=True)
    parser.add_argument("-v", "--vcf_list", help="Comma-separated list of vcf files.\nDefault: look for *CNV.vcf* files in /illumina/build/CNA/DELTA/BioSample")
    parser.add_argument("-t", "--cn_thr", help="Copy number cutoff value", type=int)
    parser.add_argument("-s", "--save_sample_files", help="Save intermediate sample specific files", action="store_true")
    options = parser.parse_args()

    if options.vcf_list:
        options.sample_vcfs = options.vcf_list.split(",")
        options.sample_names = [os.path.basename(sample_vcf).split(".vcf")[0] for sample_vcf in options.sample_vcfs]
    else:
        options.sample_vcfs = glob.glob("/illumina/build/CNA/DELTA/BioSample/*/Analysis/*.somatic.CNV.vcf.gz")
        options.sample_names = [sample_vcf.split("BioSample/")[1].split("/Analysis")[0] for sample_vcf in options.sample_vcfs]
    if options.save_sample_files:
        duplicate_sample_names = set([x for x in options.sample_names if options.sample_names.count(x) > 1])
        print "[WARNING] Found duplicate sample names (sample specific file will be overwritten): " + ",".join([str(name) for name in duplicate_sample_names])

    return options


def get_sample_features(sample_vcf, save=False, sample_name=None, out_dir="."):

    cmd = "{0} filter {1} -i '(CN==2 & MCC!=1) | CN!=2' | " \
          "{0}  view --apply-filter PASS | " \
          "{0} query -f \"[%CHROM\t%POS\t%INFO/END\t%CN\t%MCC]\n\"".format(bcftools, sample_vcf)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=None, shell=True)
    output = process.communicate()

    if output[1]:
        print output[1]

    if save:
        if sample_name:
            out_file = os.path.join(out_dir, sample_name + ".tsv")
        else:
            out_file = os.path.join(out_dir, "vcf_features.tsv")
        with open(out_file, "w") as out:
            out.write("CHR\tSTART\tEND\tCN\tMCC\n")
            out.writelines(output[0])

    return output[0]


def get_density(feature_list, start=None, end=None, step=None):

    print feature_list
    if not start:
        start = min(feature_list)
    if not end:
        end = max(feature_list)
    if not step:
        step = 1

    # Calculate smoothed frequency
    # tin parameter should contain the data points at which density is estimated
    # (if not present, they are inferred from input data)
    # (tin should contain unique equi-distant data points)
    kernel = sskernel(numpy.array(sorted(feature_list))) #, tin=range(start, end, step))
    print kernel
    return dict([(x, density) for x, density in zip(kernel["t"], kernel["y"])])



def update_feature_summary(sample_info, CN_freq_dict, MCC_freq_dict, MCC_cond_freq_dict, CN_list, MCC_list, len_list, CN_thr=None):

    for row in sample_info.split("\n"):
        if row:
            chr, start, end, CN, MCC = row.split("\t")

            if not MCC == ".":
                if int(CN) < int(MCC) or int(MCC) < float(CN)/2:
                    continue
            if CN_thr and int(CN) > CN_thr:
                continue

            else:
                CN_list.append(int(CN))
                add_to_feature_dict(CN_freq_dict, int(CN))

                if not MCC == ".":
                    MCC_list.append(int(MCC))
                    add_to_feature_dict(MCC_freq_dict, int(MCC))
                    if not int(CN) in MCC_cond_freq_dict:
                        MCC_cond_freq_dict[int(CN)] = {int(MCC) : 1}
                    else:
                        add_to_feature_dict(MCC_cond_freq_dict[int(CN)], int(MCC))

                len_list.append(int(end) - int(start))


def add_to_feature_dict(dict, feature):

    if feature not in dict:
        dict[feature] = 1
    else:
        dict[feature] += 1


def print_list(list_feature, out_file):

    with open(out_file, "w") as out:
        out.writelines([str(feat) + "\n" for feat in list_feature])


def print_freq(dict_feature, out_file):

    dict_feature_freq = {}
    tot = sum(dict_feature.values())
    for k in dict_feature:
        dict_feature_freq[k] = float(dict_feature[k])/tot
    with open(out_file, "w") as out:
        json.dump(dict_feature_freq, out)


def print_conditional_freq(dict_feature, out_file):

    dict_feature_freq = {k:v for k,v in dict_feature.items()}
    for k1 in dict_feature:
        tot = sum(dict_feature[k1].values())
        for k2 in dict_feature[k1]:
            dict_feature_freq[k1][k2] = float(dict_feature[k1][k2])/tot
    with open(out_file, "w") as out:
        json.dump(dict_feature_freq, out, indent=4)


def main():

    options = parse_options()
    CN_list, MCC_list, len_list, num_list = [], [], [], []
    CN_freq_dict, MCC_freq_dict, len_freq_dict, num_freq_dict, MCC_cond_freq_dict = {},{},{},{},{}

    for sample_name, sample_vcf in zip(options.sample_names, options.sample_vcfs):
        #print sample_name, sample_vcf
        sample_info = get_sample_features(sample_vcf, options.save_sample_files, sample_name, options.output_dir)
        update_feature_summary(sample_info, CN_freq_dict, MCC_freq_dict, MCC_cond_freq_dict, CN_list, MCC_list, len_list, options.cn_thr)
        num_list.append(len(sample_info.split("\n")))

    len_freq_dict = get_density(len_list, step=1000)
    num_freq_dict = get_density(num_list)

    # print_list(CN_list, os.path.join(options.output_dir, "CN_list.txt"))
    # print_list(MCC_list, os.path.join(options.output_dir, "MCC_list.txt"))
    # print_list(len_list, os.path.join(options.output_dir, "len_list.txt"))
    # print_list(num_list, os.path.join(options.output_dir, "num_list.txt"))

    print_freq(CN_freq_dict, os.path.join(options.output_dir, "CN_freq.json"))
    print_freq(MCC_freq_dict, os.path.join(options.output_dir, "MCC_freq.json"))
    print_conditional_freq(MCC_cond_freq_dict, os.path.join(options.output_dir, "MCC_cond_freq.json"))
    print_freq(len_freq_dict, os.path.join(options.output_dir, "len_freq.json"))
    print_freq(num_freq_dict, os.path.join(options.output_dir, "num_freq.json"))



if __name__=="__main__":
    main()
