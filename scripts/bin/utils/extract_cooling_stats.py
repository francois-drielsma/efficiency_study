import glob
import json

class ExtractCoolingStats(object):
    def __init__(self):
        self.table = []
        self.rows = []
        self.headings = []
        self.target_bin = None
      
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
            print job_name,
            data = json.loads(open(fname).read())
            self.find_significance(data)

    def find_significance(self, data):
        cdf_list = data['reco']['ratio']['corrected_cdf']
        cdf_stats_errors = data['reco']['ratio']['cdf_stats_errors']
        cdf_sys_errors = data['reco']['ratio']['cdf_sys_errors']
        cdf_sum_errors = [None]*len(cdf_list)
        max_sig = -1000.
        if self.target_bin != None:
            test_list = [self.target_bin]
        else:
            test_list = range(len(cdf_list))
        for i in test_list:
            cdf = cdf_list[i]
            cdf_stats = cdf_stats_errors[i]
            cdf_sys = cdf_sys_errors[i]
            cdf_err = (cdf_stats**2+cdf_sys**2)**0.5
            cdf_sum_errors[i] = cdf_err
            sig = (cdf-1.)/cdf_err
            if sig > max_sig:
                max_sig = sig
                max_i = i
        print "   ", format(max_i, "2d"), format(cdf_list[max_i], "8.4g"),\
              format(cdf_stats_errors[max_i], "8.4g"), \
              format(cdf_sys_errors[max_i], "8.4g"), \
              format(cdf_sum_errors[max_i], "8.4g")
        self.table[-1].append((cdf_list[max_i], cdf_stats_errors[max_i], cdf_sys_errors[max_i]))

    def new_row(self):
        self.table.append([])

    def print_table(self):
        print self.table_header+"{l|ccc}"
        print "          ",
        for heading in self.headings:
            print "&", heading.ljust(25),
        print "\\\\"
        for i, row in enumerate(self.table):
            print self.name_lookup[self.rows[i]].ljust(10),
            for cell in row:
                print "&", format(cell[0], "8.4g"), "$\pm$", format(cell[1], "8.1g"), "$\pm$", format(cell[2], "8.1g"),
            print "\\\\"
        print self.table_footer

    name_lookup = {
        "None":"No absorber",
        "lH2_empty":"Empty LH2",
        "lH2_full":"Full LH2",
        "LiH":"LiH"
    }

    table_header = """
\\begin{table*}
\\centering
\\caption{abc}
\\begin{tabular}[pos]"""
    table_footer =  """\\end{tabular}
\\end{table*}
"""



def main():
    prefix = "output/2017-02-7-v6/"
    stats = ExtractCoolingStats()
    stats.rows = ["None", "lH2_empty", "lH2_full", "LiH",]
    bin_headings = [(0, "4-140"), (0, "6-140"), (0, "10-140")]
    stats.headings = [head[1] for head in bin_headings]
    table = []
    for absorber in stats.rows:
        stats.new_row()
        for bin_n, emittance in bin_headings:
            for real in ["_"]:
                suffix = "/amplitude/amplitude.json"
                a_glob = prefix+"plots"+real+"2017-2.7_"+emittance+"*"+absorber+"*"+suffix
                stats.extract(a_glob, bin_n)
    stats.print_table()

if __name__ == "__main__":
    main()