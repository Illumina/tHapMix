import os, sys, re, subprocess
from optparse import OptionParser


def runCanvas(bam_dirs):

    for bam_dir in bam_dirs:

        if not os.path.exists(bam_dir):

            print("Simulation directory %s does not exist\n" % bam_dir)

        else:

            if os.path.exists("%s/sorted_Proteus_LP6007590_LP6007591_som_var.bam.bai" % bam_dir):
            #         and not os.path.exists("%s/pyflow.data/active_pyflow_process.txt" % bam_dir)): # simulation is finished

                print "Sim dir: " + bam_dir
                out_dir = re.sub("simulation", "Canvas", bam_dir)
                #out_dir = "/illumina/scratch/tmp/users/ccolombo/Canvas/simNormOriginal"
                sh_file = out_dir + ".sh"

                # # Create folders
                if not os.path.exists(out_dir):
                    print "Prepare Canvas run for " + os.path.basename(bam_dir)
                    os.mkdir(out_dir)
                if not os.path.exists(out_dir + "/Analysis"):
                    os.mkdir(out_dir + "/Analysis")

                # Write Sample Sheet
                with (open(out_dir + "/SampleSheet.csv", "w")) as ss:
                    ss.writelines("[Header],,,,\n"
                                  "Workflow,TumorNormal,,,\n"
                                  ",,,,\n"
                                  "[Settings],,,,\n"
                                  "RetainTempFiles,1,,,\n"
                                  ",,,,\n"
                                  "[Data],,,,\n"
                                  "SampleID,SampleName,MatchedNormalID,AlignmentPath,GenomeFolder\n"
                                  "Normal,Normal,,/illumina/build/curium/builds/CNeval/PG_bams/NA12882-50x-TR2-p1_S1_30x.bam,Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta\n"
                                  "Tumor,Tumor,Normal,%s/sorted_Proteus_LP6007590_LP6007591_som_var.bam,"
                                  "Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta" % bam_dir)

                # Write bash script to launch Canvas
                with (open(sh_file, "w")) as sh:
                    sh.writelines("#!/bin/bash\n"
                                  "/illumina/development/Isis/2.6.26/Isis"
                                  " -r %s -a %s/Analysis -c 3" % (out_dir, out_dir))
                os.system("chmod +x " + sh_file)

                if os.path.exists("%s/Analysis/WorkflowLog.txt" % out_dir):
                    print "Canvas already running or Canvas finished for %s" % os.path.basename(bam_dir)

                else:
                    # Run Canvas
                    print "Run Canvas for " + os.path.basename(bam_dir)
                    subprocess.call(["/bin/bash", "-i", "-c", "qsubIt A " + sh_file])


def main():

    parser = OptionParser()
    parser.add_option("-s", "--sim_id", dest="sim_id", type="string", help="simulation id")
    parser.add_option("-d", "--sim_dirs", dest="sim_dirs", type="string", action="append", help="simulation directories")
    (options, args) = parser.parse_args()

    if options.sim_dirs is None and options.sim_id is None:
        parser.error("Specify simulation id or simulation directories")
    if options.sim_dirs is not None and options.sim_id is not None:
        parser.error("Simulation id and simulation directories are mutually exclusive options")
    if options.sim_dirs is not None:
        bam_dirs = options.sim_dirs
    if options.sim_id is not None:
        sim_dir = "/illumina/scratch/tmp/users/ccolombo/simulation/"
        print "Looking for simulation directories in %s, simulation id: %s" % (sim_dir, options.sim_id)
        bam_dirs = [sim_dir + dir for dir in os.listdir(sim_dir) if re.search(options.sim_id + ".*", dir)]
        print "Found: " + str(bam_dirs)

    runCanvas(bam_dirs)


if __name__ == "__main__":
    main()
