#!/usr/bin/env python
import os
import sys
import re
import copy
import json
from HapMixUtil import *

ScriptDir=os.path.abspath(os.path.dirname(__file__))
pyFlowPath=os.path.join(ScriptDir,"../redist","pyflow","src")
sys.path.append(pyFlowPath)
sys.path.append(ScriptDir)
from pyflow import WorkflowRunner


class TumorBEDWorkflow(WorkflowRunner):
    """Workflow for the creation of a tumor truth file"""

    def __init__(self, params):
        self.params = params

    def workflow(self):
        createTumorTruthFile = os.path.join(ScriptDir + "/createTumorTruthFile.py")
        ttTask = "python " + createTumorTruthFile
        ttTask += " -n %d" % self.params.num_var_min
        ttTask += " -m %d" % self.params.num_var_max
        ttTask += " -s %d" % self.params.size_var_min
        ttTask += " -t %d" % self.params.size_var_max
        ttTask += " -c %d" % self.params.cn_min
        ttTask += " -d %d" % self.params.cn_max
        ttTask += " -o %s" % self.params.tum_truth_file
        ttTask += " -g %s" % self.params.hg_file
        #print(ttTask)
        self.addTask("createTumorTruthFile", ttTask)



class ClonalFilesWorkflow(WorkflowRunner):
    """Workflow for the creation of clones truth files"""

    def __init__(self, params):
        self.params = params

    def workflow(self):
        createClonalFiles = os.path.join(ScriptDir + "/createClonalFiles.py")
        ctTask = "python " + createClonalFiles
        if self.params.tree_format in ["random_binary", "random_single_level"]:
            ctTask += " -f %s" % self.params.tree_format
        else:
             ctTask += " -f fixed"
        ctTask += " -t %s" % self.params.tum_truth_file
        ctTask += " -c %s" % " -c ".join(self.params.clones)
        if hasattr(self.params, "tree_structure_file"): ctTask += " -u %s" % self.params.tree_structure_file
        if hasattr(self.params, "tree_seed"): ctTask += " -e %d" % self.params.tree_seed
        ctTask += " -v %f" % self.params.var_het
        if self.params.mutate_sv: ctTask += " -m"
        ctTask += " -s %s" % self.params.somatic_vcf
        ctTask += " -o %s" % self.params.output_dir
        ctTask += " -g %s" % self.params.hg_file
        #print(ctTask)
        self.addTask("createClonalFiles", ctTask)


class RunHapMixWorkflow(WorkflowRunner):
    """Workflow for the creation of BAM files for simulated cancer evolution"""

    def __init__(self, params):
        self.params = params

    def workflow(self):

        samtools=os.path.join(ScriptDir,"../redist","samtools-1.2/samtools")


        # Truth files names
        clone_truth_file = {}
        for clID in self.params.clones_perc.keys():
             ctfile = os.path.join(self.params.output_dir, os.path.basename(self.params.tum_truth_file)[:-4]) + "_clone" + clID + ".bed"
             clone_truth_file[clID] = add_BED_complement(ctfile, self.params.hg_file, sort=True, out_dir=self.params.tmp_dir)
        norm_truth_file = "%s/%s_norm.bed" % (self.params.output_dir, os.path.basename(self.params.tum_truth_file)[:-4])

        # Compute clone ploidy depth for each clone and the normal sample
        overall_ploidy = 2*(1 - self.params.purity/100)
        for clID, perc in self.params.clones_perc.items():
            overall_ploidy += get_clone_ploidy(clone_truth_file[clID], self.params.hg_file, ["chrX", "chrY", "chrM"]) * perc * (self.params.purity/100)
        self.params.ploidy_depth = self.params.ploidy_depth/overall_ploidy*2
        print "Overall ploidy: %.2f\n" % overall_ploidy
        print "Ploidy depth: %.2f\n" % self.params.ploidy_depth
        clone_ploidy_depth = {}
        for clID, perc in self.params.clones_perc.items():
                clone_ploidy_depth[clID] = round(self.params.ploidy_depth * perc * (self.params.purity/100), 2)
        norm_ploidy_depth = round((self.params.ploidy_depth * (1 - self.params.purity/100)), 2)
        # # Write clone ploidy depths to files
        with open(os.path.join(self.params.tmp_dir, "ploidy_depths.txt"), "w") as pl_file:
                pl_file.writelines(["Overall ploidy: %.2f\n" % overall_ploidy, "Ploidy depth: %.2f\n" % self.params.ploidy_depth])
                pl_file.writelines(["Clone {}: {}\n".format(k, str(v)) for (k, v) in clone_ploidy_depth.items()] + ["Norm : " + str(norm_ploidy_depth) + "\n"])

        # If clone ploidy depths or normal ploidy depths are below threshold, do not launch the simulation
        if (self.params.purity!=0 and any([cpd < 4 for cpd in clone_ploidy_depth.values()])) or (norm_ploidy_depth != 0 and norm_ploidy_depth < 4):
            print("Ploidy depth is too low: ending simulation.\n")

        else:

            # Get header and chromosome X, Y, M bam files
            header = os.path.join(self.params.tmp_dir, "header.sam")
            self.addTask("getHeader", samtools + " view -H %s/sorted_haplotypeD_chr21.bam -o %s" % (self.params.haplotyped_bam_dir, header))
            self.addTask("sorted_chrX", samtools + " view -Sb %s -o %s/sorted_%s_chrX.bam" % (header, self.params.output_dir, self.params.output_bam_file), dependencies="getHeader")
            self.addTask("sorted_chrY", samtools + " view -Sb %s -o %s/sorted_%s_chrY.bam" % (header, self.params.output_dir, self.params.output_bam_file), dependencies="getHeader")
            self.addTask("sorted_chrM", samtools + " view -Sb %s -o %s/sorted_%s_chrM.bam" % (header, self.params.output_dir, self.params.output_bam_file), dependencies="getHeader")

            # Simulate each chromosome
            chrTaskIDs = []
            simulate_copy_number_changes = os.path.join(ScriptDir + "/simulate_copy_number_changes.py")

            chroms = []
            if not isinstance(self.params.chromosomes, (list, tuple)):
                chroms = range(1,23)
            else:
                temp = self.params.chromosomes
                for i in temp: 
                    if isinstance(i, int):
                        chroms.append(i)
                    else:
                        sys.exit("chromosomes %s in config file should either be a comma separated integers or a keyword all\n" % params.chromosomes)

            for chrID in chroms:
                chrTaskID = "simChr" + str(chrID)

                # Simulate clones
                clTaskIDs = []
                for clID, perc in self.params.clones_perc.items():
                    clTaskID = chrTaskID + "cl" + clID
                    clTaskIDs.append(clTaskID)
                    clTask = sys.executable
                    clTask += " " + simulate_copy_number_changes
                    clTask += " -a %s/sorted_haplotypeC_chr%d.bam" % (self.params.haplotyped_bam_dir, chrID)
                    clTask += " -b %s/sorted_haplotypeD_chr%d.bam" % (self.params.haplotyped_bam_dir, chrID)
                    clTask += " -o %s/%s_chr%d_cl%s.bam" % (self.params.output_dir, self.params.output_bam_file, chrID, clID)
                    clTask += " -c chr%d" % (chrID)
                    clTask += " -t %s" % (clone_truth_file[clID])
                    clTask += " -m %s" % str(self.params.mutate_sv)
                    clTask += " -s %s/%s_clone%s.vcf" % (self.params.output_dir, os.path.basename(self.params.somatic_vcf[:-4]), clID)
                    clTask += " -d %f" % (clone_ploidy_depth[clID])
                    clTask += " -e %d" % (self.params.bam_depth)
                    self.addTask(clTaskID, clTask, memMb=2024)
                    #print(clTask)

                # Simulate normal
                normTaskID = chrTaskID + "norm"
                clTaskIDs.append(normTaskID)
                normTask = sys.executable
                normTask += " " + simulate_copy_number_changes
                normTask += " -a %s/sorted_haplotypeC_chr%d.bam" % (self.params.haplotyped_bam_dir, chrID)
                normTask += " -b %s/sorted_haplotypeD_chr%d.bam" % (self.params.haplotyped_bam_dir, chrID)
                normTask += " -o %s/%s_chr%d_norm.bam" % (self.params.output_dir, self.params.output_bam_file, chrID)
                normTask += " -c chr%d" % (chrID)
                normTask += " -t %s" % norm_truth_file
                normTask += " -m False"
                normTask += " -s %s" % (self.params.somatic_vcf)
                normTask += " -d %f" % (norm_ploidy_depth)
                normTask += " -e %d" % (self.params.bam_depth)
                self.addTask(normTaskID, normTask, memMb=2024)
                #print(normTask)

                # Concatenate clones for the current chromosome
                catChrTask = samtools + " cat"
                catChrTask += " -h %s" % header
                catChrTask += " -o %s/%s_chr%d.bam" % (self.params.output_dir, self.params.output_bam_file, chrID)
                prefix = " %s/%s_chr%d_" % (self.params.output_dir, self.params.output_bam_file, chrID)
                catChrTask += " ".join([prefix + "cl" + cl + ".bam" for cl in self.params.clones_perc.keys()])
                catChrTask += " %s/%s_chr%d_norm.bam" % (self.params.output_dir, self.params.output_bam_file, chrID)
                self.addTask(chrTaskID + "cat", catChrTask, dependencies=clTaskIDs+["getHeader"])
                #print(catChrTask)

                # Remove clonal bam files
                rm_files = []
                for cl in self.params.clones_perc.keys():
                    rm_files.append("%s/%s_chr%d_cl%s.bam" % (self.params.output_dir, self.params.output_bam_file, chrID, cl))
                rm_files.append("%s/%s_chr%d_norm.bam" % (self.params.output_dir, self.params.output_bam_file, chrID))
                # print "rm " + " ".join(rm_files)
                self.addTask(chrTaskID + "ClRm", "rm " + " ".join(rm_files), dependencies=(chrTaskID + "cat"))

                # Sort chromosome bam files
                sortChrTask = samtools + " sort"
                sortChrTask += " %s/%s_chr%d.bam" % (self.params.output_dir, self.params.output_bam_file, chrID)
                sortChrTask += " -T %s/sorted_%s_chr%d" % (self.params.output_dir, self.params.output_bam_file, chrID)
                sortChrTask += " -o %s/sorted_%s_chr%d.bam" % (self.params.output_dir, self.params.output_bam_file, chrID)
                self.addTask(chrTaskID + "sort", sortChrTask, dependencies=(chrTaskID + "cat"))
                #print(sortChrTask)
                chrTaskIDs.append(chrTaskID + "sort")

                # Remove unsorted bam files
                # print "rm %s/%s_chr%d.bam" % (self.params.output_dir, self.params.output_bam_file, chrID)
                self.addTask(chrTaskID + "UnsortedRm", "rm %s/%s_chr%d.bam" % (self.params.output_dir, self.params.output_bam_file, chrID), dependencies=(chrTaskID + "sort"))

            # Concatenate chromosome bam files
            catTask = samtools + " cat"
            catTask += " -h %s" % header
            catTask += " -o %s/sorted_%s_som_var.bam " % (self.params.output_dir, self.params.output_bam_file)
            prefix = "%s/sorted_%s_chr" % (self.params.output_dir, self.params.output_bam_file)
            chr_names = ["M"] + [str(i) for i in chroms] + ["X", "Y"]
            catTask += " ".join([prefix + c + ".bam" for c in chr_names])
            self.addTask("catBamFiles", catTask, dependencies=chrTaskIDs)
            # self.addTask("catBamFiles", catTask)
            # print(catTask)

            # Remove chromosome bam files
            rm_files = []
            for chr_name in chr_names:
                rm_files.append(prefix + chr_name + ".bam")
            self.addTask("rmChrBamFiles", "rm " + " ".join(rm_files), dependencies=("catBamFiles"))

            # Index output bam file
            indTask = samtools + " index %s/sorted_%s_som_var.bam" % (self.params.output_dir, self.params.output_bam_file)
            self.addTask("indexBam", indTask, dependencies="catBamFiles")
            # print(indTask)



class SimFullWorkflow(WorkflowRunner):
    """Complete workflow for simulation of cancer evolution"""

    def __init__(self, params):
        self.params = params

    def workflow(self):

        self.params.perc_maj = float(self.params.perc_maj)
        self.params.var_het = float(self.params.var_het)
        self.params.purity = float(self.params.purity)

        self.params.clones = get_clones_names(self.params.tree_format, self.params.num_clones)
        print self.params.clones
        if self.params.tree_format in ["random_binary", "random_single_level"]:
            self.params.name_maj_clone = self.params.clones[0]
        if self.params.name_maj_clone not in self.params.clones:
            sys.exit("Majority clone name %s not in evolutionary tree %s\n" % self.params.name_maj_clone, str(self.params.clones))

        # Compute clonal percentages
        if self.params.num_clones == 1:
            self.params.clones_perc = {self.params.clones[0]: 1.0}
        else:
            if self.params.perc_maj > 0:
                perc_min = (1 - self.params.perc_maj / 100) / (self.params.num_clones - 1)
                clones_min = self.params.clones[:]
                clones_min.remove(self.params.name_maj_clone)
                self.params.clones_perc = dict(zip([self.params.name_maj_clone] + clones_min, [self.params.perc_maj / 100] + [perc_min] * (self.params.num_clones - 1)))
            else: # there is no major clone, split evenly all clones
                perc = round(1 / float((self.params.num_clones)), 2)
                self.params.clones_perc = dict(zip(self.params.clones, [perc] * self.params.num_clones))

        # Create output directories for the simulation
        if not os.path.exists(self.params.output_dir): # already exists of single simulation is run, to be created if simulation protocol is run
            print("Creating output directory: %s\n" % self.params.output_dir)
            os.mkdir(self.params.output_dir)
        self.params.tmp_dir = os.path.join(self.params.output_dir, "Tmp")
        if not os.path.exists(self.params.tmp_dir): # already exists of single simulation is run, to be created if simulation protocol is run
            print("Creating directory for tmp files: %s\n" % self.params.tmp_dir)
            os.mkdir(self.params.tmp_dir)
        snv_out_dir = os.path.join(self.params.output_dir, "forced_somatic_snv_frequencies")
        if not os.path.exists(snv_out_dir):
            os.mkdir(snv_out_dir)

        # Create filtered genome file (chr 1-23, X, Y, M)
        chr_names = ["chr" + str(i) for i in range(1,23)] + ["chrM", "chrX", "chrY"]
        with open(self.params.hg_file, "r") as curr_gf:
            self.params.hg_file = os.path.join(self.params.tmp_dir, "filtered_" + os.path.basename(self.params.hg_file))
            with open(self.params.hg_file, "w") as new_gf:
                new_gf.writelines([chr_line for chr_line in curr_gf if chr_line.strip().split("\t")[0] in chr_names])

        # Write parameters of current simulation to file
        with open(self.params.output_dir + "/sim_param.json", "w") as par_file:
            json.dump(vars(self.params), par_file, indent=4, separators=(',', ': '))

        # Run simulation
        if self.params.create_truth_bed:
            tumBED = self.addWorkflowTask("newTumorBED", TumorBEDWorkflow(self.params))
            clonFiles = self.addWorkflowTask("createClonesFiles", ClonalFilesWorkflow(self.params), dependencies=tumBED)
        else:
            clonFiles = self.addWorkflowTask("createClonesFiles", ClonalFilesWorkflow(self.params))
        runSim = self.addWorkflowTask("runSimulation", RunHapMixWorkflow(self.params), dependencies=clonFiles)


class ProtocolFullWorkflow(WorkflowRunner):
    """Complete workflow for simulation protocol (parameters variation: n, m, v)"""

    def __init__(self, params):
        self.all_params = params
        if type(self.all_params.num_clones) is not list: self.all_params.num_clones = [self.all_params.num_clones]
        if type(self.all_params.perc_maj) is not list: self.all_params.perc_maj = [self.all_params.perc_maj]
        if type(self.all_params.var_het) is not list: self.all_params.var_het = [self.all_params.var_het]

    def workflow(self):

        for num_clones in self.all_params.num_clones:
            perc_maj_range, var_het_range = self.all_params.perc_maj, self.all_params.var_het
            if num_clones == 1:
                perc_maj_range, var_het_range = [-1], [-1]
            for perc_maj in perc_maj_range:
                if num_clones==2 and perc_maj==50 and -1 in perc_maj_range:
                    continue
                for var_het in var_het_range:

                    sim_params = copy.copy(self.all_params)
                    sim_params.num_clones = num_clones
                    sim_params.perc_maj = perc_maj
                    sim_params.var_het = var_het

                    sim_output_dir = "sim%s_n%d_m%d_v%d" % (os.path.basename(sim_params.tum_truth_file)[:-4], num_clones, perc_maj, var_het)
                    sim_params.output_dir = os.path.join(self.all_params.output_dir, sim_output_dir)
                    if sim_params.tree_format == "random_single_level" or (isinstance(sim_params.tree_format, (list, tuple)) and len(sim_params.tree_format)>2): sim_params.output_dir += "_sl"

                    wfID = os.path.basename(sim_params.output_dir)
                    self.addWorkflowTask(wfID, SimFullWorkflow(sim_params))



