import sys
import os
import site
import json
import shutil
import importlib

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
import scripts.mc_plotter
import scripts.residual_fitter
import scripts.plot_residual_fitter
import scripts.utilities
config_file = None

def maus_globals(config):
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

def file_mangle(config_anal, config_file_name):
    print "Clearing old data"
    try:
        shutil.rmtree(config_anal["plot_dir"])
    except OSError:
        pass
    os.makedirs(config_anal["plot_dir"])
    shutil.copy(config_file_name, config_anal["plot_dir"])

def do_analysis(config, analysis_index):
    config_anal = config.analyses[analysis_index]
    print config_anal["name"]
    all_event_count = 0

    for max_spill in [config.number_of_spills]:
        if max_spill == 0:
            continue
        data_loader = scripts.data_loader.DataLoader(config, analysis_index)
        data_loader.load_data(0, 100000)
    config.maus_version = data_loader.maus_version
    tm_p_list = []

    file_mangle(config_anal, sys.argv[1])

    print "Using p bins", config_anal["p_bins"]
    if config_anal["do_magnet_alignment"]:
        print "Doing magnet alignment"
        fitter = scripts.residual_fitter.ResidualFitter(config, config_anal, data_loader)
        fitter.fit()
        #scripts.plot_residual_fitter.do_plots(config)
    if config_anal["do_extrapolation"]:
        # nb: also does the "aperture_ds", "aperture_us", "delta_tof01" cuts
        print "Doing track extrapolation"
        scripts.extrapolate_through_detectors.do_extrapolation(config, config_anal, data_loader)
    data_loader.update_cuts()
    if config_anal["do_amplitude"]:
        scripts.amplitude_analysis.do_amplitude_analysis(config, config_anal, data_loader)
    data_loader.update_cuts()
    if config_anal["do_mc"]:
        print "Doing MC"
        scripts.mc_plotter.MCPlotter.do_mc_plots(config, config_anal, data_loader)
    wiki_summary = "<plots disabled>"
    if config_anal["do_plots"]:
        print "Doing plots"
        wiki_summary = scripts.data_plotter.do_plots(config, config_anal, data_loader)
    return wiki_summary

def main():
    scripts.utilities.set_palette()
    config = config_file.Config()
    maus_globals(config)

    wiki_summary_list = []
    for i, anal in enumerate(config.analyses):
        wiki_summary_list.append(do_analysis(config, i))
    for summary in wiki_summary_list:
        print summary
    

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "Usage: python calculate_emittance.py </path/to/config/script>"
        sys.exit(1)
    config_mod = sys.argv[1].replace(".py", "")
    config_mod = config_mod.replace("/", ".")
    print "Using configuration module", config_mod
    config_file = importlib.import_module(config_mod)
    main()
    print "Done - press <CR> to finish"
    #raw_input()

