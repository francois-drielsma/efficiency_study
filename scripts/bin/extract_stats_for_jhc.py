import glob
import json

class ExtractCoolingStats(object):
    def __init__(self):
        self.table = []
        self.rows = []
        self.headings = []
        self.job_name = ""
      
    def extract(self, source_glob, target_bin):
        self.target_bin = target_bin
        self.source_glob = source_glob
        self.load_json()

    def load_json(self):
        #print "Loading from", self.source_glob
        glob_list = sorted(glob.glob(self.source_glob))
        if len(glob_list) > 1:
            raise RuntimeError("Glob list too long")
        for fname in glob_list:
            job_name = fname.split("plots_")[1]
            job_name = job_name.split("/ampl")[0]
            job_name = job_name.replace("2017-2.7_", "")
            job_name = job_name.replace("_", " ")
            job_name = job_name.ljust(26)
            print job_name
            self.data = json.loads(open(fname).read())
            self.print_one_dataset()

    def print_one_dataset(self):
        for sample in "all_upstream", "all_downstream":
            print "   ", sample
            for key in "pdf", "corrected_pdf", "pdf_sys_errors":
                print "       ", key,
                for value in self.data['reco'][sample][key]:
                    print format(value, "4.4g")+",",
                print
        print


def main():
    prefix = "output/2017-02-7-v6/"
    stats = ExtractCoolingStats()
    stats.rows = ["None", "lH2_empty", "lH2_full", "LiH",]
    bin_headings = [(0, "4-140"), (0, "6-140"), (0, "10-140")]
    stats.headings = [head[1] for head in bin_headings]
    table = []
    for absorber in stats.rows:
        for bin_n, emittance in bin_headings:
            for real in ["_"]:
                suffix = "/amplitude/amplitude.json"
                a_glob = prefix+"plots"+real+"2017-2.7_"+emittance+"*"+absorber+"*"+suffix
                stats.extract(a_glob, bin_n)


if __name__ == "__main__":
    main()