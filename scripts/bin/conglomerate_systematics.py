import copy
import sys
import os
import shutil

import ROOT

import utilities.utilities as utilities
import utilities.root_style as root_style
import conglomerate
from conglomerate.compare_config import CompareConfig
from conglomerate.conglomerate_merge import ConglomerateMerge
from conglomerate.conglomerate_one import ConglomerateOne
from conglomerate.conglomerate_one import ConglomerateContainer
from conglomerate.merge_cuts_summary_tex import MergeCutsSummaryTex

class CompareCutsSystematicConfig(CompareConfig):
    def __init__(self, beam, absorber, target_dir, top_labels, right_labels, sys_run, sys_abs, systematic, plot, hist, x_range, normalise):
        dir_list = [
            target_dir+"plots_"+beam+"_"+absorber+"/",
            target_dir+"plots_Simulated_"+beam+"_"+absorber+"/",
            sys_run+"/plots_Simulated_"+beam+"_"+sys_abs+"_Systematics_"+systematic+"/",
        ]
        output_dir = target_dir+"/compare_systematics/"+beam+"_"+absorber+"/"
        title = (beam+" "+absorber+"    "+systematic).replace("_", " ")
        modifiers = {
            "extra_labels":{
                "right":right_labels,
                "top":top_labels
            },
            "hist_title":title,
            "file_name":plot,
            "canvas_name":plot,
            "normalise_hist":normalise,
            "histogram_names":[hist+" us cut"],
            "mice_logo":False,
            "legend":False,
            "calculate_errors":[],
            "redraw":{
                "draw_option":["P", "", "P"], #E1 PLC
                "fill_color":[1, ROOT.kOrange-2, 1],
                "transparency":None,
                "line_color":[1, 1, 1],
                "marker_style":[20, 20, 25],
                "draw_order":[1, 2, 0],
                "x_range":x_range,
                "y_range":None,
                "graph_draw_option":None,
                "ignore_more_histograms":False,
            },
            "rescale_x":x_range,
            "rescale_y":True,
            "write_plots":{
                "dir":output_dir,
                "file_name":plot+"_vs_"+systematic
            },
            "axis_title":{"x":self.labels[plot]},
        }
        self.conglomerate_list = [
            self.get_conglomerate_0(modifiers = modifiers),
        ]
        self.data_caption = [[],]
        self.setup(beam, target_dir, "data_plots/", "compare_systematics/", dir_list)

    labels = {
        "tku_p":"p at TKU Reference Plane [MeV/c]",
        "chi2_tku":"#chi^{2}/D.o.F. in TKU",
        "chi2_tkd":"#chi^{2}/D.o.F. in TKD",
        "p_res":"p(TKU) - p(TKD) [MeV/c]",
        "tkd_p":"p at TKD Reference Plane [MeV/c]",
    }

def systematics():
    target_dir = "output/2017-02-7-Systematics-test-2/"
    dir_list = [
        "2017-2.7_3-140_lH2_empty_Systematics_tku_base",  "2017-2.7_4-140_lH2_empty_Systematics_tku_base",
        "2017-2.7_6-140_lH2_empty_Systematics_tku_base",  "2017-2.7_10-140_lH2_empty_Systematics_tku_base", 
    ]
    mc_cuts_summary = MergeCutsSummaryTex()
    mc_cuts_summary.table_ref_pre = "systematics_"
    for beam in dir_list:
        config = CompareConfig()
        dir_list = [
            #target_dir+"plots_"+beam+"/",
            target_dir+"plots_Simulated_"+beam
        ]
        config.setup(beam, target_dir, "cut_plots/", "compare_cuts/", dir_list)
        config.mc_caption = [["" for i in range(5)] for i in range(5)]
        mc_prefix = ["upstream reconstructed", "downstream reconstructed", "extrapolated reconstructed",
                     "upstream truth", "downstream truth"]
        for i, prefix in enumerate(mc_prefix):
            for j, item in enumerate(config.mc_caption[i]):
                config.mc_caption[i][j] = "The "+prefix+" simulated sample is listed. "+item
        mc_cuts_summary.append_summary(config, [0])

    mc_cuts_summary.caption = config.mc_caption
    mc_cuts_summary.merge_summaries(target_dir+"/cuts_summary/mc/", "mc_cuts_summary")
    #raw_input("Done systematics - press <CR> to continue to the rest")

def systematics_conglomerate(dir_lists,
                             target_dir,
                             top_labels,
                             right_labels):
    rows = len(dir_lists)
    cols = min([len(sub_list) for sub_list in dir_lists])
    dir_list = []
    for sub_list in dir_lists:
        dir_list += sub_list
    sys_run = "output/2017-02-7-Systematics-test/"
    sys_abs = "lH2_empty"
    systematic = "tku_base"
    sys_list =  []
    for systematic in ["tku_base", "tku_scale_SSUE2_plus", "tku_scale_SSUC_plus",
                       "tku_scale_SSUE1_plus", "tku_pos_plus", "tku_rot_plus"]:
        sys_list += [
            (systematic, "lH2_empty", "tku_p", "tku_p", [130., 150.], False),
            (systematic, "lH2_empty", "chi2_tku", "chi2", [0., 20.], True),
            (systematic, "lH2_empty", "p_res", "p_res", [-50., 50.], True),
        ]
    for systematic in ["tku_base", "tkd_scale_SSDE2_plus", "tkd_scale_SSDC_plus",
                       "tkd_scale_SSDE1_plus", "tkd_pos_plus", "tkd_rot_plus"]:
        sys_list += [
            (systematic, "lH2_empty", "chi2_tkd", "chi2", [0., 20.], True),
            (systematic, "lH2_empty", "p_res", "p_res", [-50., 50.], True),
            (systematic, "lH2_empty", "tkd_p", "tkd_p", [100., 160.], False),
        ]

    for systematic, sys_abs, fname, hist, x_range, normalise in sys_list:
        conglomerate_list = []
        for beam in dir_list:
            [beam, absorber] = beam.split("140_")
            beam += "140"
            try:
                config = CompareCutsSystematicConfig(beam,
                                                    absorber,
                                                    target_dir,
                                                    top_labels,
                                                    right_labels,
                                                    sys_run,
                                                    sys_abs,
                                                    systematic,
                                                    fname,
                                                    hist,
                                                    x_range,
                                                    normalise)
                cong = ConglomerateContainer(config)
                cong.conglomerate()
                conglomerate_list.append(cong)
            except Exception:
                sys.excepthook(*sys.exc_info())
        try:
            if len(dir_lists) > 1:
                merge = ConglomerateMerge(conglomerate_list)
                merge.merge_all(rows, cols)
        except Exception:
            sys.excepthook(*sys.exc_info())

def mkdirs(target_dir, dir_list):
    for beam_list in dir_list:
        for beam in beam_list:
            plot_dir = target_dir+"/compare_systematics/"+beam
            if os.path.exists(plot_dir):
                shutil.rmtree(plot_dir)
            os.makedirs(plot_dir)

def main():
    #systematics()
    root_style.setup_gstyle()
    ROOT.gROOT.SetBatch(True)
    target_dir = "output/2017-02-7-test/"
    my_dir_list = [
        #["2017-2.7_4-140_None",      "2017-2.7_6-140_None",      "2017-2.7_10-140_None",],
        ["2017-2.7_4-140_lH2_empty"]#, "2017-2.7_6-140_lH2_empty", "2017-2.7_10-140_lH2_empty",],
        #["2017-2.7_4-140_lH2_full",  "2017-2.7_6-140_lH2_full",  "2017-2.7_10-140_lH2_full",],
        #["2017-2.7_4-140_LiH",       "2017-2.7_6-140_LiH",       "2017-2.7_10-140_LiH"],
    ]
    top_labels = ["4-140", "6-140", "10-140"]
    right_labels = ["Empty\nLH2", "Full\nLH2", "LiH"] #"No\nabsorber", 
    mkdirs(target_dir, my_dir_list)
    systematics_conglomerate(my_dir_list,
                             target_dir,
                             top_labels,
                             right_labels)

if __name__ == "__main__":
    main()
    if not ROOT.gROOT.IsBatch():
        raw_input("Finished - press <CR> to end")
