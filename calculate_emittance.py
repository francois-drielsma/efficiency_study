import sys
import os
import site
import json
import shutil

import numpy
import ROOT
ROOT.gROOT.SetBatch(True)

import xboa.common
import Configuration
import maus_cpp.globals
import maus_cpp.field
import maus_cpp.polynomial_map
import libxml2

import scripts.extrapolate_through_detectors
import scripts.amplitude_analysis
from scripts.tm_calculator import TMCalculator
from scripts.data_plotter import DataPlotter
from scripts.tm_calculator import TOF12Predicate
import scripts.data_loader
#import scripts.config_reco_2016_04_1_3_extrapolate as config_file
import scripts.config_reco_2016_04_1_2 as config_file


def maus_globals(config):
    str_conf = Configuration.Configuration().\
                                      getConfigJSON(command_line_args=False)
    json_conf = json.loads(str_conf)
    json_conf["simulation_geometry_filename"] = config.geometry
    json_conf["verbose_level"] = 1
    maus_conf = json.dumps(json_conf)
    maus_cpp.globals.birth(maus_conf)
    print maus_cpp.field.str(True)

def do_analysis(config, analysis_index):
    config_anal = config.analyses[analysis_index]
    print config_anal["name"]
    all_event_count = 0

    for max_spill in [config.number_of_spills]:
        if max_spill == 0:
            continue
        data_loader = scripts.data_loader.DataLoader(config, analysis_index)
        data_loader.load_data(0, 100000)

    tm_p_list = []
    print "Clearing old data"
    try:
        shutil.rmtree(config_anal["plot_dir"])
    except OSError:
        pass
    os.makedirs(config_anal["plot_dir"])
    print "Using p bins", config_anal["p_bins"]
    for p_low, p_high in config_anal["p_bins"]:
        if not config_anal["do_amplitude"]:
            continue
        p_bin = round((p_low+p_high)/2., 1)
        bunch_us, bunch_ds = scripts.amplitude_analysis.make_bunches(data_loader, p_low, p_high)
        canvas_pdf, canvas_ratio = scripts.amplitude_analysis.delta_amplitude_plot(bunch_us, bunch_ds, config.analyses[analysis_index]["plot_dir"], p_bin)
        bunch_us, bunch_ds = None, None
    if config_anal["do_extrapolation"]:
        scripts.extrapolate_through_detectors.do_extrapolation(config, config_anal, data_loader)

    xboa.common.clear_root()
    plotter = DataPlotter(config, analysis_index, data_loader.events, lambda event: event["any_cut"], data_loader.run_numbers)
    plotter.print_cuts_summary()
    cov_us, cov_ds = plotter.print_covariance_matrix()
    sys.stdout.flush()
    plotter.plot_tof01()
    plotter.plot_tof12()
    plotter.bunch_plots()
    plotter.plot_var_2d("x", "upstream", "px", "upstream")
    plotter.plot_var_2d("y", "upstream", "py", "upstream")
    plotter.plot_var_2d("x", "upstream", "y", "upstream", True)
    plotter.plot_var_2d("x", "upstream", "y", "upstream", False)
    plotter.plot_var_2d("x", "downstream", "px", "downstream")
    plotter.plot_var_2d("y", "downstream", "py", "downstream")
    plotter.plot_var_1d("x", "upstream")
    plotter.plot_var_1d("px", "upstream")
    plotter.plot_var_1d("y", "upstream")
    plotter.plot_var_1d("py", "upstream")
    plotter.plot_var_1d("x", "downstream")
    plotter.plot_var_1d("px", "downstream")
    plotter.plot_var_1d("y", "downstream")
    plotter.plot_var_1d("py", "downstream")
    plotter.plot_var_1d("p", "upstream")
    plotter.plot_var_1d("r", "upstream")
    plotter.plot_var_1d("p", "downstream")
    plotter.plot_var_1d("r", "downstream")
    plotter.plot_var_2d("r", "downstream", "pt", "downstream")
    plotter.plot_p_tot_res()
    plotter.plot_p_tot_vs_tof()
    accepted, accepted_data = plotter.get_cuts_summary()
    wiki_summary = plotter.print_wiki_summary()
    print wiki_summary, "\n\n"
    return wiki_summary

def main():
    config = config_file.Config()
    maus_globals(config)

    wiki_summary_list = []
    for i, anal in enumerate(config.analyses):
        wiki_summary_list.append(do_analysis(config, i))
    for summary in wiki_summary_list:
        print summary
    

if __name__ == "__main__":
    main()
    print "Done - press <CR> to finish"
    #raw_input()

