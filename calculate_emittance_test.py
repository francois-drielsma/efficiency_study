import sys
import os
import site
import json
import shutil
import importlib
import time

import numpy
import ROOT
ROOT.gROOT.SetBatch(True)

import xboa.common
import Configuration
import maus_cpp.globals
import maus_cpp.field
import maus_cpp.polynomial_map
import libxml2

try:
    import scripts.extrapolate_through_detectors
except ImportError:
    print "Note you are not set up for extrapolation... so this 'do_extrapolation' = True will throw an exception"
import scripts.amplitude_analysis
from scripts.tm_calculator import TMCalculator
from scripts.data_plotter import DataPlotter
from scripts.tm_calculator import TOF12Predicate
import scripts.data_loader
import scripts.mc_plotter_test
import scripts.data_plotter_test
import scripts.residual_fitter
import scripts.plot_residual_fitter
import scripts.utilities
config_file = None

class Analyser(object):
    def __init__(self):
        config_mod = sys.argv[1].replace(".py", "")
        config_mod = config_mod.replace("/", ".")
        print "Using configuration module", config_mod
        config_file = importlib.import_module(config_mod)
        scripts.utilities.set_palette()
        scripts.utilities.set_root_verbosity(5)
        self.config = config_file.Config()
        self.config_anal = None
        self.maus_globals(self.config)
        self.analysis_list = []

    def do_analysis(self):
        for i, self.config_anal in enumerate(self.config.analyses):
            self.init_phase()
            self.birth_phase()
            self.process_phase()
            self.print_phase()

    def maus_globals(self, config):
        try:
            os.makedirs("logs/tmp") # location for magnet cached maps
        except OSError:
            pass # probably the directory existed already

        str_conf = Configuration.Configuration().\
                                          getConfigJSON(command_line_args=False)
        json_conf = json.loads(str_conf)
        json_conf["simulation_geometry_filename"] = config.geometry
        json_conf["verbose_level"] = config.maus_verbose_level
        maus_conf = json.dumps(json_conf)
        maus_cpp.globals.birth(maus_conf)
        print maus_cpp.field.str(True)

    def file_mangle(self, config_file_name):
        print "Clearing old data"
        try:
            shutil.rmtree(self.config_anal["plot_dir"])
        except OSError:
            pass
        os.makedirs(self.config_anal["plot_dir"])
        shutil.copy(config_file_name, self.config_anal["plot_dir"])

    def init_phase(self):
        self.data_loader = scripts.data_loader.DataLoader(self.config, self.config_anal)
        self.data_loader.get_file_list()
        if self.config_anal["do_mc"]:
            print "Doing MC"
            self.analysis_list.append(scripts.mc_plotter_test.MCPlotter(self.config, self.config_anal, self.data_loader))
        if self.config_anal["do_plots"]:
            print "Doing plots"
            self.analysis_list.append(scripts.data_plotter_test.DataPlotter(self.config, self.config_anal, self.data_loader))

    def birth_phase(self):
        all_event_count = 0

        self.data_loader.load_spills(self.config.preanalysis_number_of_spills)
        self.config.maus_version = self.data_loader.maus_version

        self.file_mangle(sys.argv[1])

        for analysis in self.analysis_list:
            analysis.birth()
        self.data_loader.clear_data()

    def process_phase(self):
        now = time.time()
        while self.data_loader.load_spills(self.config.analysis_number_of_spills):
            for analysis in self.analysis_list:
                analysis.process()
            self.data_loader.clear_data()
            if now - time.time() > 1800: # 30 mins
                self.print_phase()

    def print_phase(self):
        for analysis in self.analysis_list:
            analysis.death()

    def finalise_phase(self):
        pass

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "Usage: python calculate_emittance.py </path/to/config/script>"
        sys.exit(1)
    analyser = Analyser()
    analyser.do_analysis()
    print "Done - press <CR> to finish"
    #raw_input()
