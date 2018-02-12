import json

import numpy.random
from scripts.analysis_base import AnalysisBase

class Tof01Weighting(AnalysisBase):
    def __init__(self, config, config_anal, data_loader):
        super(Tof01Weighting, self).__init__(config, config_anal, data_loader)
        self.data_loader = data_loader
        self.mode = self.config_anal["weight_tof01_mode"]
        self.target_distribution = None
        self.weights = None
        self.n_bins = 101

    def get_bin_number(self, value):
        tof_min = self.config_anal["tof01_cut_low"]
        tof_step = (self.config_anal["tof01_cut_high"] - tof_min)/float(self.n_bins-1)
        bin_number = int((value-tof_min)/tof_step)
        if bin_number < 0:
            bin_number = 0
        if bin_number > self.n_bins-1:
            bin_number = self.n_bins-1
        return bin_number

    def birth_read_mode(self):
        source_distribution = open(self.config_anal["weight_tof01_source"])
        target_distribution = open(self.config_anal["weight_tof01_target"])

        weight_source = json.loads(source_distribution.read())
        weight_target = json.loads(target_distribution.read())
        print weight_source
        print weight_target
        #if weight_source["cuts"] != weight_target["cuts"]:
        #    raise RuntimeError("Reading data set with different cuts in TOF01 weighting")
        self.weights = [0. for src in weight_source["weights"]]
        for i, src in enumerate(weight_source["weights"]):
            tgt = weight_target["weights"][i]
            if tgt == 0.:
                self.weights[i] = 1.
            else:
                self.weights[i] = src/tgt
        max_weight = max(self.weights)
        self.weights = [wt/max_weight for wt in self.weights]
        self.target_distribution = None

    def birth_write_mode(self):
        self.target_distribution = open(self.config_anal["weight_tof01_target"], "w")
        self.weights = [0. for i in range(self.n_bins)]

    def process_read_mode(self):
        for event in self.data_loader.events:
            if not event["upstream_cut"]:
                a_bin = self.get_bin_number(event["tof01"])
                will_cut = numpy.random.uniform(0, 1) < self.weights[a_bin]
                event["will_cut"]["tof01_selection"] = will_cut

    def process_write_mode(self):
        for event in self.data_loader.events:
            if not event["upstream_cut"]:
                a_bin = self.get_bin_number(event["tof01"])
                self.weights[a_bin] += 1.

    def process(self):
        if self.mode == "sample_using_distribution":
            self.process_read_mode()
        elif self.mode == "build_distribution":
            self.process_write_mode()
        else:
            raise RuntimeError("Ddin't recognise mode "+str(self.mode))

    def birth(self):
        if self.mode == "sample_using_distribution":
            self.birth_read_mode()
        elif self.mode == "build_distribution":
            self.birth_write_mode()
        else:
            raise RuntimeError("Ddin't recognise mode "+str(self.mode))
        self.process()

    def death(self):
        if self.mode != "build_distribution":
            return
        target_distribution = open(self.config_anal["weight_tof01_target"], "w")
        data = {
            "cuts":self.config.upstream_cuts,
            "weights":self.weights
        }
        print >> target_distribution, json.dumps(data, indent = 2)
        target_distribution.close()
        print "Writing target distribution"
        print open(self.config_anal["weight_tof01_target"]).read()
