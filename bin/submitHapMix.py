#!/usr/bin/env python

import os, sys, json, copy
import multiprocessing
from HapMixWorkflow import *


class SimParams:
    """Object to store all the parameters of the simulation workflow"""
    
    pass


def get_HapMix_params():
    """Get simulation parameters from command line and configuration file into a SimParams object"""

    params = SimParams()

    # Parameters from command line
    options = get_commandline_options()
    params.output_dir = options.output_dir
    params.tmp_dir = options.tmp_dir
    params.mode = options.mode

    # Truth file configuration
    if options.truth_config_file:
        params.create_truth_bed = True
        parse_truth_config_file(options.truth_config_file, params)
    else:
        params.create_truth_bed = False
        params.tum_truth_file = options.tum_truth_file

    # Simulation configuration
    parse_sim_config_file(options.sim_config_file, params)

    return params


def parse_truth_config_file(truth_config_file, params):
    """Parse HapMix configuration file for the synthetic tumor truth bed file into params object"""

    config = json.load(open(truth_config_file))

    # Parameters defining the truth bed file to be created
    params.num_var_min = config["num_var_min"]
    params.num_var_max = config["num_var_max"]
    params.size_var_min = config["size_var_min"]
    params.size_var_max = config["size_var_max"]
    params.cn_min = config["cn_min"]
    params.cn_max = config["cn_max"]
    params.tum_truth_file = os.path.join(params.output_dir, config["truth_file_name"])



def parse_sim_config_file(sim_config_file, params):
    """Parse HapMix configuration file for the simulation parameters into params object"""

    config = json.load(open(sim_config_file))
    config_dir = os.path.dirname(os.path.realpath(sim_config_file))

    # Parameters defining the clonal evolution
    params.purity = float(config["purity"])
    params.num_clones = config["num_clones"]
    if params.num_clones == 1:
        params.perc_maj = -1
        params.var_het = -1
    else:
        params.perc_maj = config["perc_majority_clone"]
        params.var_het = config["perc_heterogeneous_variants"]
    params.tree_format = config["tree_format"]
    validate_tree_structure(params.tree_format, params.num_clones)
    if params.tree_format in ["random_binary", "random_single_level"]:
        if "tree_seed" in config:
            params.tree_seed = config["tree_seed"]
    elif isinstance(params.tree_format, (list, tuple)):
        params.tree_structure_file = os.path.join(params.tmp_dir, "tree_structure.json")
        json.dump(params.tree_format, open(params.tree_structure_file, "w"))
        if params.num_clones != 1:
            params.name_maj_clone = config["majority_clone"]
    params.chromosomes = config["chromosomes"]
    if isinstance(params.chromosomes, (list, tuple)):
        if len(params.chromosomes) == 0:
            print "\nchromosomes in config file should either be a comma separated integers or a keyword all\n"
            sys.exit(2)
    elif params.chromosomes.strip() != "all":
        print "\nchromosomes in config file should either be a comma separated integers or a keyword all\n"
        sys.exit(2)    

    # Parameters for the simulation workflow
    params.haplotyped_bam_dir = config["haplotyped_bam_dir"]
    params.output_bam_file = config["output_bam_file"]
    params.mutate_sv = config["mutate_sv"]
    params.somatic_vcf = os.path.join(config_dir, config["somatic_vcf"])
    params.hg_file = os.path.join(config_dir, config["hg_file"])
    params.bam_depth = config["bam_depth"]
    params.ploidy_depth = config["ploidy_depth"]

    if not os.path.exists(params.hg_file):
        print "\nGenome file does not exist\n"
        sys.exit(2)
    if not os.path.exists(params.haplotyped_bam_dir):
        print "\nDirectory for haplotyped bam files does not exist\n"
        sys.exit(2)
    if params.mutate_sv and not os.path.exists(params.somatic_vcf):
        print "\nSomatic mutation VCF file does not exist\n"
        sys.exit(2)


def get_commandline_options():
    """Get user options from command line"""

    import argparse

    parser = argparse.ArgumentParser(prog='HapMix v0.6')
    required_arguments = parser.add_argument_group('Required arguments (-t and -b mutually exclusive!)')
    required_arguments.add_argument("-c", "--sim_config_file", dest="sim_config_file", help="configuration file containing the HapMix simulation parameters")
    truth_file_args = required_arguments.add_mutually_exclusive_group(required=True)
    truth_file_args.add_argument("-b", "--truth_config_file", dest="truth_config_file", help="configuration file containing the parameters for creating a synthetic truth bed file\n(new truth file is created)")
    truth_file_args.add_argument("-t", "--tumor_truth_file", dest="tum_truth_file", help="name of the tumor truth file\n(existing truth file is used)")
    required_arguments.add_argument("-o", "--output_dir", dest="output_dir", help="name of the output directory")
    required_arguments.add_argument('-m', '--mode', dest='mode', default='local', choices=['sge', 'local'], help="select run mode (local|sge)")
    options = parser.parse_args()

    if not options.sim_config_file:
        print("\nConfiguration file for the HapMix simulation is not specified\n")
        parser.print_help()
        sys.exit(2)
    if not os.path.exists(options.sim_config_file):
        print("\nConfiguration file for the HapMix simulation does not exist\n")
        parser.print_help()
        sys.exit(2)
    if options.truth_config_file and (not os.path.exists(options.truth_config_file)):
        parser.print_help()
        print("\nConfiguration file for the tumor truth bed does not exist\n")
        sys.exit(2)
    if options.tum_truth_file and (not os.path.exists(options.tum_truth_file)):
        parser.print_help()
        print("\nTumor truth bed file does not exist\n")
        sys.exit(2)

    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)
    options.tmp_dir = os.path.join(options.output_dir, "Tmp")
    if not os.path.exists(options.tmp_dir):
        os.mkdir(options.tmp_dir)

    return options



def main():

    params = get_HapMix_params()
    print vars(params)

    if isinstance(params.num_clones, list) or isinstance(params.perc_maj, list) or isinstance(params.var_het, list):
        wflow = ProtocolFullWorkflow(params)
    else:
        wflow = SimFullWorkflow(params)

    exitpath=os.path.join(params.output_dir,"hapmix.workflow.exitcode.txt")

    if params.mode == "local":  
        nCores=multiprocessing.cpu_count()
    else:
        nCores=128

    try:
        retval=wflow.run(mode=params.mode,
                    nCores=nCores,
                    dataDirRoot=params.output_dir,
                    isContinue="Auto",
                    isForceContinue=True)
    finally:
        exitfp=open(exitpath,"w")
        exitfp.write("%i\n" % (retval))
        exitfp.close()



if __name__=="__main__":

    main()
