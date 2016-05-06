import os, re
from optparse import OptionParser

### prova


def run_FREEC(bam_dirs, out_gen_dir, bam_fname, config_template):

    with open(config_template, "r") as ct:
        template = ct.readlines()

    for bam_dir in bam_dirs:

        if not bam_dir.endswith("-1"):

            # not checked: simulation may still be running (indexBam step)
            if not os.path.exists(bam_dir):
                print("Simulation directory %s does not exist\n" % bam_dir)
            elif not os.path.exists("%s.bai" % os.path.join(bam_dir, bam_fname)):
                print("Simulated sample %s does not exist in %s\n" % (bam_fname, bam_dir))
            else:
                #
                # print "Sim dir: " + bam_dir
                sim_id = os.path.basename(bam_dir)
                out_dir = os.path.join(out_gen_dir, sim_id)
                config_file = os.path.join(out_dir, "config.txt")

                if not os.path.exists(out_dir):
                    # Create folder for results
                    os.mkdir(out_dir)

                # Write configuration file
                for ix, line in enumerate(template):
                    if line.startswith("mateFile"):
                        template[ix] = "mateFile = %s\n" % os.path.join(bam_dir, bam_fname)
                        break
                with (open(config_file, "w")) as cf:
                    cf.writelines(template)

                if not os.path.exists(os.path.join(out_dir, bam_fname + "_CNVs")):
                    # Run FREEC
                    print "Run FREEC for " + sim_id + "\n"
                    os.chdir(out_dir)
                    os.system("echo '/illumina/build/CNA/CNV/haplotype_simulation/FREEC/freec -conf config.txt' | qsub -cwd -N '%s_freec'" % sim_id)


def main():

    parser = OptionParser()
    parser.add_option("-s", "--sim_id", dest="sim_id", type="string", help="simulation id")
    parser.add_option("-d", "--sim_dirs", dest="sim_dirs", type="string", action="append", help="simulation directories")
    parser.add_option("-b", "--bam_fname", dest="bam_fname", type="string", help="name of the simulated bam file")
    parser.add_option("-o", "--out_dir", dest="out_dir", type="string", help="main output directory")
    parser.add_option("-c", "--config_template", dest="config_template", type="string", help="template file for config.txt")
    (options, args) = parser.parse_args()

    if options.sim_dirs is None and options.sim_id is None:
        parser.error("Specify simulation id or simulation directories")
    if options.sim_dirs is not None and options.sim_id is not None:
        parser.error("Simulation id and simulation directories are mutually exclusive options")
    if options.sim_dirs is not None:
        bam_dirs = options.sim_dirs
    if options.sim_id is not None:
        sim_dir = "/illumina/scratch/tmp/users/ccolombo/simulation/"
        print "Looking for simulation directories in: " + sim_dir
        bam_dirs = [sim_dir + dir for dir in os.listdir(sim_dir) if re.search("sim" + options.sim_id + ".*", dir)]
        if len(bam_dirs) == 0:
            print "No simulation directories found"
    if options.bam_fname is None:
        options.bam_fname = "sorted_Proteus_LP6007590_LP6007591_som_var.bam"
        print "Using default name for simulated bam: " + options.bam_fname
    if options.out_dir is None or not os.path.exists(options.out_dir):
        options.out_dir = "/illumina/build/CNA/CNV/haplotype_simulation/FREEC/"
        print "Using default output directory: " + options.out_dir
    if options.config_template is None or not os.path.exists(options.config_template):
        options.config_template = "/illumina/build/CNA/CNV/haplotype_simulation/FREEC/config_template.txt"
        print "Using default configuration template file: " + options.config_template

    run_FREEC(bam_dirs, options.out_dir, options.bam_fname, options.config_template)


if __name__ == "__main__":
    main()
