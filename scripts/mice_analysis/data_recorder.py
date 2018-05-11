import sys
import operator
import xboa.common
import json
import copy
import math
import numpy

import ROOT
from xboa.bunch import Bunch
import utilities.utilities
import utilities.cdb_tof_triggers_lookup

from analysis_base import AnalysisBase

class DataRecorder(AnalysisBase):
    def __init__(self, config, config_anal, data_loader):
        super(DataRecorder, self).__init__(config, config_anal, data_loader)
        self.data_loader = data_loader
        self.config = config
        self.config_anal = config_anal
        self.hit_vars = ("x", "px", "y", "py", "pz")
        self.upstream_detector = "tku_5"
        self.downstream_detector = "tkd_tp"
        self.hits_file_dict = {}

    def birth(self):
        self.set_plot_dir("data_recorder")
        
        self.hits_file_dict = {
          self.upstream_detector:{
            "vars":self.hit_vars,
            "cut":"upstream_cut",
            "file_name":self.plot_dir+"/"+self.upstream_detector+".json",
            "detector":self.upstream_detector,
            "file":None
          },
          self.downstream_detector:{
            "vars":self.hit_vars,
            "cut":"downstream_cut",
            "file_name":self.plot_dir+"/"+self.downstream_detector+".json",
            "detector":self.downstream_detector,
            "file":None
          }
        }
        for key, file_metadata in self.hits_file_dict.iteritems():
            file_name = file_metadata["file_name"]
            fout = open(file_name, "w")
            fout.write(json.dumps(file_metadata)+"\n")
            file_metadata["file"] = fout
        self.process()

    def process(self):
        for event in self.data_loader.events:
            for key, file_data in self.hits_file_dict.iteritems():
                cut = file_data["cut"]
                if event[cut]:
                    continue
                detector = file_data["detector"]
                fout = file_data["file"]
                for hit in event["data"]:
                    if hit["detector"] == detector:
                        self.write_hit(hit, fout)

    def write_hit(self, hit, fout):
        hit = hit["hit"]
        out = [hit[var] for var in self.hit_vars]
        fout.write(json.dumps(out)+"\n")

    def death(self):
        for file_data in self.hits_file_dict.values():
            file_data["file"].close()

