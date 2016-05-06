import os
import sys
import re
import json
from optparse import OptionParser
from getVarPerc import create_het_cn_file

sys.path.append("/illumina/development/pyflow/pyflow-1.1.6/src")
scriptDir = os.path.abspath(os.path.dirname(__file__))

from pyflow import WorkflowRunner


class CovBedWorkflow(WorkflowRunner):
    """Workflow for the creation of a tumor truth file"""

    def __init__(self, cnv_file, cnv_type, part_file, excl_file, perc, tum_pur=None):
        self.cnv_file = cnv_file
        self.cnv_type = cnv_type
        self.part_file = part_file
        self.excl_file = excl_file
        self.perc = perc
        self.tum_pur = tum_pur

    def workflow(self):

        cTask = "python /home/ccolombo/HapMix/getCanvasCovBed.py"
        cTask += " -i %s" % self.cnv_file
        cTask += " -c %s" % self.cnv_type
        cTask += " -d %s" % self.part_file
        if self.tum_pur:
            cTask += " -t %f" % self.tum_pur
        if self.excl_file:
            cTask += " -e %s" % self.excl_file
        if self.perc:
            cTask += " -p"
        print cTask
        self.addTask("createCovBed", cTask, memMb=5*1024)


class CovPlotWorkflow(WorkflowRunner):
    """Workflow for the creation of a tumor truth file"""

    def __init__(self, sim):
        self.sim = sim

    def workflow(self):

        if True: # and not os.path.exists("/illumina/scratch/tmp/users/ccolombo/Canvas/Het/plots/cov_" + self.sim + "/cov_" + sim + ".png"):
            pTask = "cat /home/ccolombo/HapMix/canvasHet.R | /illumina/thirdparty/R/R-3.1.2/bin/R --slave --quiet"
            pTask += " --args %s" % ("cov_" + self.sim)
            #self.addTask("plotCovBedCanvas", pTask, memMb=4*1024)
            self.addTask("plotCovBedCanvas", pTask,  memMb=4*1024)


class CovFullWorkflow(WorkflowRunner):
    """Workflow for the creation of a tumor truth file"""

    def __init__(self, canvas_sim_dirs, excl_file, plot):
        self.canvas_sim_dirs = canvas_sim_dirs
        self.excl_file = excl_file
        self.plot = plot

    def workflow(self):

        for canvas_sim_dir in self.canvas_sim_dirs:

            print canvas_sim_dir

            if not os.path.exists(canvas_sim_dir):
                sys.stderr.write("Canvas sim directory %s does not exist, coverage analysis interrupetd\n" % canvas_sim_dir)
                continue
            if not (os.path.exists(os.path.join(canvas_sim_dir, "Analysis/Normal_Tumor_G1_P1.somatic.SV.vcf.gz"))
                    or os.path.exists(os.path.join(canvas_sim_dir, "Analysis/Normal_Tumor_G1_P1.somatic.SV.vcf"))):
                sys.stderr.write("Canvas is not finished for %s\n" % os.path.basename(canvas_sim_dir))
                continue
            if os.path.basename(canvas_sim_dir).startswith("simNorm") or canvas_sim_dir.endswith("_purity100"):
                continue

            vcf_file = os.path.join(canvas_sim_dir, "Analysis/Normal_Tumor_G1_P1.somatic.SV.vcf")
            sim = os.path.basename(canvas_sim_dir)
            true_tum_purity = 0.8
            sim_dir = re.sub("Canvas", "simulation", canvas_sim_dir)

            #truth_file = create_het_cn_file(sim_dir, None, True)
            sim_par = json.load(open(os.path.join(sim_dir, "sim_param.json")))
            truth_file = sim_par["tum_truth_file"]
            print truth_file

            if os.path.basename(canvas_sim_dir) in ["simNorm", "simNormAll", "simNormOriginal"]:
                truth_file = "/home/ccolombo/HapMix/norm.bed"
            part_file = os.path.join(canvas_sim_dir, "Analysis/TempCNV_G1_P1/G1_P1.partitioned")

            #self.addWorkflowTask(sim + "TruthPlot", CovPlotWorkflow(sim + "_truth"))
            #self.addWorkflowTask(sim + "CanvasPlot", CovPlotWorkflow(sim))
            self.addWorkflowTask(sim + "TruthBed", CovBedWorkflow(truth_file, "truth", part_file, self.excl_file, True, true_tum_purity))
            self.addWorkflowTask(sim + "CanvasBed", CovBedWorkflow(vcf_file, "canvas", part_file, self.excl_file, True), dependencies=sim + "TruthBed")
            #self.addWorkflowTask(sim + "CanvasBed", CovBedWorkflow(vcf_file, "canvas", part_file, self.excl_file, True))
            if self.plot:
                self.addWorkflowTask(sim + "TruthPlot", CovPlotWorkflow(sim + "_truth"), dependencies=sim + "TruthBed")
                self.addWorkflowTask(sim + "CanvasPlot", CovPlotWorkflow(sim), dependencies=sim + "CanvasBed")

def main():

    parser = OptionParser()
    parser.add_option("-s", "--canvas_sim_id", dest="canvas_sim_id", type="string", help="simulation id (Canvas results directories: sim[sim_id])")
    parser.add_option("-d", "--canvas_sim_dir", dest="canvas_sim_dir", type="string", help="Canvas results directory")
    parser.add_option("-e", "--excl_file", dest="excl_file", type="string", help="optional bed file for excluded regions")
    parser.add_option("-p", "--plot", dest="plot", action="store_true", default=False, help="Plot results")
    #parser.add_option("-t", "--truth_file", dest="truth_file", type="string", help="truth file of the current simulation")
    (options, args) = parser.parse_args()

    if options.canvas_sim_dir is None and options.canvas_sim_id is None:
        parser.error("Specify Canvas results simulation id or directory")
    if options.canvas_sim_dir is not None and options.canvas_sim_id is not None:
        parser.error("Canvas results simulation id and directory are mutually exclusive options")
    if options.canvas_sim_dir is not None:
        canvas_sim_dirs = [options.canvas_sim_dir]
        options.canvas_sim_id = os.path.basename(canvas_sim_dirs[0])
    if options.canvas_sim_id is not None:
        canvas_dir = "/illumina/scratch/tmp/users/ccolombo/Canvas/"
        canvas_sim_dirs = [canvas_dir + cdir for cdir in os.listdir(canvas_dir) if os.path.isdir((canvas_dir + cdir)) and cdir.startswith(options.canvas_sim_id)]
        if len(canvas_sim_dirs) == 0:
            print "No Canvas results directories found"

    wflow = CovFullWorkflow(canvas_sim_dirs, options.excl_file, options.plot)
    py_dir = os.path.join("/illumina/scratch/tmp/users/ccolombo/evaluation/EvaluateCov", "pyflow." + options.canvas_sim_id + "plotcanvas")
    print "Pyflow directory: %s" % py_dir
    wflow.run(mode="sge", dataDirRoot=py_dir)





if __name__ == "__main__":
    main()




