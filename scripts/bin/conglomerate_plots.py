import sys

import ROOT

import utilities
import root_style
import conglomerate
from conglomerate.compare_config import CompareConfig
from conglomerate.conglomerate_merge import ConglomerateMerge
from conglomerate.conglomerate_one import ConglomerateOne
from conglomerate.conglomerate_one import ConglomerateContainer
from conglomerate.merge_cuts_summary_tex import MergeCutsSummaryTex


class MCConfig(CompareConfig):
    def __init__(self, mc_version):
        self.conglomerate_list = [
            self.get_conglomerate_1("tku_p", "tku_p", "all", [80., 200.], False),
            self.get_conglomerate_1("tof01", "tof01", "all", [24., 40.,], True),
            self.get_conglomerate_1("tof01", "tof01", "us cut", [28., 32.,], True),
            self.get_conglomerate_1("tku_p", "tku_p", "us cut", [130., 160.], False),
            self.get_conglomerate_1("tku_x", "tku_x", "us cut", None, False),
            self.get_conglomerate_1("tku_y", "tku_y", "us cut", None, False),
            self.get_conglomerate_1("tku_px", "tku_px", "us cut", None, False),
            self.get_conglomerate_1("tku_py", "tku_py", "us cut", None, False),
            self.get_conglomerate_1("chi2_tku", "chi2", "us cut", None, True),
            self.get_conglomerate_1("chi2_tkd", "chi2", "us cut", None, True),
            self.get_conglomerate_1("tof0_x", "tof0_x", "us cut", None, False),
            self.get_conglomerate_1("tof0_y", "tof0_y", "us cut", None, False),
            self.get_conglomerate_1("tof1_x", "tof1_x", "us cut", None, False),
            self.get_conglomerate_1("tof1_y", "tof1_y", "us cut", None, False),
        ]

class CompareCutsConfig(CompareConfig):
    def __init__(self, beam, target_dir):
        self.setup(beam, target_dir, "cut_plots/", "compare_cuts/", True, True)
        self.conglomerate_list = [
            self.get_conglomerate_2("global_through_virtual_diffuser_ds_r_9_0", None, "Radius at diffuser (upstream) [mm]", None, True, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("global_through_virtual_diffuser_us_r_9_0", None, "Radius at diffuser (downstream) [mm]", None, True, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("tkd_n_tracks_13_0", None, "Number of tracks in TKD", None, True, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("tkd_chi2_13_0", None, "#chi^{2}/D.o.F. in TKD", None, True, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("tkd_max_r_13_0", None, "Maximum radius in TKD [mm]", [0., 300.], True, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("tkd_p_13_0", None, "Momentum in TKD [MeV/c]", None, True, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("tku_chi2_9_0", None, "#chi^{2}/D.o.F. in TKU", None, True, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("tku_max_r_9_0", None, "Maximum radius in TKU [mm]", None, True, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("tku_n_tracks_8_0", None, "Number of tracks in TKU", None, True, [0.5, 0.5, 0.9, 0.9]), # disable two cuts
            self.get_conglomerate_2("tku_p_9_0", None, "Momentum in TKU [MeV/c]", [100., 220.], True, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("tof_tof0_n_sp_9_0", None, "Number of space points in TOF0", None, True, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("tof_tof1_n_sp_9_0", None, "Number of space points in TOF1", None, True, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("tof_tof01_9_0", None, "Time between TOF0 and TOF1 [ns]", [28., 33.], True, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("tof_delta_tof01_9_0", None, "t(TOF01) - extrapolated t(TOF01) [ns]", None, True, [0.1, 0.5, 0.5, 0.9]),

            self.get_conglomerate_2("tku_scifi_n_planes_with_clusters_10_0", None, "Number of planes with clusters in TKU", None, True, [0.1, 0.5, 0.5, 0.9]),
            self.get_conglomerate_2("tku_scifi_n_planes_with_clusters_8_0",  None, "Number of planes with clusters in TKU", None, True, [0.1, 0.5, 0.5, 0.9]),
            self.get_conglomerate_2("tku_scifi_n_planes_with_clusters_8_1",  None, "Number of planes with clusters in TKU", None, True, [0.1, 0.5, 0.5, 0.9]),
            self.get_conglomerate_2("tkd_scifi_n_planes_with_clusters_10_1", None, "Number of planes with clusters in TKD", None, True, [0.1, 0.5, 0.5, 0.9]),
            self.get_conglomerate_2("tkd_scifi_n_planes_with_clusters_11_1", None, "Number of planes with clusters in TKD", None, True, [0.1, 0.5, 0.5, 0.9]),
            self.get_conglomerate_2("tkd_scifi_n_planes_with_clusters_14_0", None, "Number of planes with clusters in TKD", None, True, [0.1, 0.5, 0.5, 0.9]),
        ]
        self.data_caption = """The sample selection criteria are listed sequentially, in the order 
that the selection was made. The number of events surviving the 
selection and all preceding selections is listed for each of the data runs studied
in this note."""
        self.mc_caption = "Simulated sample selection"

class CompareData1DConfig(CompareConfig):
    def __init__(self, beam, target_dir):
        self.setup(beam, target_dir, "data_plots/", "compare_data/", True, True)

        self.conglomerate_list = [
            self.get_conglomerate_1("tof0_slab_dt", "tof0_slab_dt", "us_cut", "Slab dt for TOF0 [ns]", [-2, 2], False, [0.55, 0.5, 0.9, 0.9]),
            self.get_conglomerate_1("tof1_slab_dt", "tof1_slab_dt", "us_cut", "Slab dt for TOF1 [ns]", [-2, 2], False, [0.55, 0.5, 0.9, 0.9]),
            self.get_conglomerate_1("tof2_slab_dt", "tof2_slab_dt", "ds_cut", "Slab dt for TOF2 [ns]", [-2, 2], False, [0.55, 0.5, 0.9, 0.9]),
            self.get_conglomerate_1("tku_x",  "tku_x",  "us_cut", "x in TKU [mm]",        [-149, 300], False, [0.55, 0.5, 0.9, 0.9]),
            self.get_conglomerate_1("tku_y",  "tku_y",  "us_cut", "y in TKU [mm]",        [-149, 300], False, [0.55, 0.5, 0.9, 0.9]),
            self.get_conglomerate_1("tku_px", "tku_px", "us_cut", "p_{x} in TKU [MeV/c]", [-95, 200], False, [0.55, 0.5, 0.9, 0.9]),
            self.get_conglomerate_1("tku_py", "tku_py", "us_cut", "p_{y} in TKU [MeV/c]", [-95, 200], False, [0.55, 0.5, 0.9, 0.9]),
            self.get_conglomerate_1("tku_p", "tku_p", "us_cut", "p in TKU [MeV/c]",       [132, 160], False, [0.55, 0.5, 0.9, 0.9]),
            self.get_conglomerate_1("tkd_x",  "tkd_x",  "ds_cut", "x in TKD [mm]",        [-149, 300], False, [0.55, 0.5, 0.9, 0.9]),
            self.get_conglomerate_1("tkd_y",  "tkd_y",  "ds_cut", "y in TKD [mm]",        [-149, 300], False, [0.55, 0.5, 0.9, 0.9]),
            self.get_conglomerate_1("tkd_px", "tkd_px", "ds_cut", "p_{x} in TKD [MeV/c]", [-95, 200], False, [0.55, 0.5, 0.9, 0.9]),
            self.get_conglomerate_1("tkd_py", "tkd_py", "ds_cut", "p_{y} in TKD [MeV/c]", [-95, 200], False, [0.55, 0.5, 0.9, 0.9]),
            self.get_conglomerate_1("tkd_p", "tkd_p", "ds_cut", "p in TKD [MeV/c]",       [89, 200], False, [0.55, 0.5, 0.9, 0.9]),
        ]


class CompareData2DConfig(CompareConfig):
    def __init__(self, beam, target_dir):
        self.setup(beam, target_dir, "data_plots/", "compare_data/", True, False)

        self.conglomerate_list = [
            self.get_conglomerate_3("p_res_vs_global_through_virtual_absorber_centre_r_ds_cut", "p_res_vs_global_through_virtual_absorber_centre_r_ds_cut", None, None),
            self.get_conglomerate_3("p_res_vs_global_through_virtual_absorber_centre_y_ds_cut", "p_res_vs_global_through_virtual_absorber_centre_y_ds_cut", None, None),
        ]


class CompareOpticsConfig(CompareConfig):
    def __init__(self, beam, target_dir):
        self.setup(beam, target_dir, "optics_plots/", "compare_optics/", True, False)

        self.conglomerate_list = [
            self.get_conglomerate_graph("beta_4d_ds", "z [m]", "#beta_{4D} [mm]"),
            self.get_conglomerate_graph("beta_x_ds", "z [m]", "#beta_{x} [mm]"),
            self.get_conglomerate_graph("beta_y_ds", "z [m]", "#beta_{y} [mm]"),
            self.get_conglomerate_graph("sigma_0_ds", "z [m]", "#sigma_{x} [mm]"),
            self.get_conglomerate_graph("sigma_2_ds", "z [m]", "#sigma_{y} [mm]"),
        ]



class CompareGlobalsConfig(CompareConfig):
    def __init__(self, beam, target_dir):
        self.setup(beam, target_dir, "global_plots/", "compare_globals/", True, False)

        self.conglomerate_list = [
            self.get_conglomerate_2("global_through_residual_tkd_tp_p", "res_ex_cut", "Residual Momentum in TKD [MeV/c]", [-20, 20], False, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("global_through_residual_tkd_tp_px", "res_ex_cut", "Residual P_{x} in TKD [MeV/c]", None, False, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("global_through_residual_tkd_tp_py", "res_ex_cut", "Residual P_{y} in TKD [MeV/c]", None, False, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("global_through_residual_tkd_tp_pz", "res_ex_cut", "Residual P_{z} in TKD [MeV/c]", [-20, 20], False, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("global_through_residual_tkd_tp_x", "res_ex_cut", "Residual x in TKD [mm]", None, False, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("global_through_residual_tkd_tp_y", "res_ex_cut", "Residual y in TKD [mm]", None, False, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("global_through_residual_tof0_t", "res_ex_cut", "Residual t in TOF0 [ns]", [-5, 5], False, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("global_through_residual_tof1_x", "res_ex_cut", "Residual x in TOF1 [mm]", None, False, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("global_through_residual_tof1_y", "res_ex_cut", "Residual y in TOF1 [mm]", None, False, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("global_ds_residual_tof2_x", "res_ex_cut", "Residual x in TOF2 [mm]", None, False, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("global_ds_residual_tof2_y", "res_ex_cut", "Residual y in TOF2 [mm]", None, False, [0.5, 0.5, 0.9, 0.9]),
            self.get_conglomerate_2("global_through_residual_tof2_t", "res_ex_cut", "Residual t in TOF2 [ns]", [-5, 5], False, [0.5, 0.5, 0.9, 0.9]),
        ]
        for item in self.conglomerate_list:
            if "tof0_t" in item["file_name"]:
                continue
            item["rebin"] = 4

class CompareAmplitudeConfigMC(CompareConfig): # MC corrections
    def __init__(self, beam, target_dir):
        self.setup(beam, target_dir, "amplitude/", "compare_amplitude_mc/", False, True)

        self.conglomerate_list = [
            self.get_conglomerate_3("crossing_probability_downstream", "crossing_probability_downstream", None, None),
            self.get_conglomerate_3("crossing_probability_upstream", "crossing_probability_upstream", None, None),
            self.get_conglomerate_graph("amplitude_inefficiency_all_upstream", "US True Amplitude [mm]",
                                        "#frac{Number in MC True (all)}{Number in MC True (reco)}",
                                        "inefficiency_all_upstream", ["inefficiency"], ["inefficiency"]), 
            self.get_conglomerate_graph("amplitude_inefficiency_all_downstream", "DS True Amplitude [mm]",
                                        "#frac{MC True (all) Number}{MC True (reco) Number}", 
                                        "inefficiency_all_downstream", ["inefficiency"], ["inefficiency"]), 
        ]

class CompareAmplitudeConfigData(CompareConfig): # data plots
    def __init__(self, beam, target_dir):
        self.setup(beam, target_dir, "amplitude/", "compare_amplitude_data/", True, False)

        self.conglomerate_list = [
            self.get_conglomerate_3("migration_matrix", "migration_matrix", None, None),
            self.get_conglomerate_3("amplitude_delta_reco_hist", "delta_amplitude_hist", None, None),
            #self.get_conglomerate_graph("amplitude_pdf_reco", "Reconstructed Amplitude [mm]",
            #                            "Number",
            #                            "amplitude_pdf_reco", ["Upstream_hist"],
            #                           ["Upstream", "Downstream", "Raw scraped", "Raw upstream", "Raw downstream"]),
            self.get_conglomerate_graph("amplitude_pdf_reco", "Reconstructed Amplitude [mm]",
                                        "Number",
                                        "amplitude_pdf_reco", ["Upstream_hist"],
                                        ["Upstream", "Downstream"], x_range=[1., 65.]),
            #self.get_conglomerate_graph("amplitude_pdf_reco", "Reconstructed Amplitude [mm]",
            #                            "Number",
            #                            "amplitude_pdf_reco", ["Upstream_hist"],
            #                            ["Upstream", "Downstream", "Raw scraped", "Raw upstream", "Raw downstream"]),
            self.get_conglomerate_graph("amplitude_cdf_reco", "Reconstructed Amplitude [mm]",
                                        "Cumulative density",
                                        "amplitude_cdf_reco", ["Upstream CDF_hist"],
                                        ["Upstream_CDF", "Downstream_CDF"], x_range=[1., 65.]),
        ]

class CompareAmplitudeConfigBoth(CompareConfig): # comparisons
    def __init__(self, beam, target_dir):
        self.setup(beam, target_dir, "amplitude/", "compare_amplitude_both/", True, True)

        self.conglomerate_list = [
            self.get_conglomerate_graph("pdf_ratio*", "Reconstructed Amplitude [mm]",
                                        "#frac{Number out}{Number in}",
                                        "pdf_ratio", ["PDF Ratio_hist", "PDF Ratio_hist"],
                                        ["PDF_Ratio"], x_range = [0.01, 110.], y_range = [0.601, 1.399], replace_hist = True,
                                        graph_marker_style=[20, 26], graph_marker_color=[1, ROOT.kRed]),

            self.get_conglomerate_graph("cdf_ratio*", "Reconstructed Amplitude [mm]",
                                        "#frac{Cumulative number out}{Cumulative number in}",
                                        "cdf_ratio", ["CDF Ratio_hist"],
                                        ["CDF_Ratio"], x_range = [0.01, 59.9], y_range = [0.601, 1.399], replace_hist = True,
                                        graph_marker_style=[20, 26], graph_marker_color=[1, ROOT.kRed]),
        ]

def cuts_summary(dir_list, target_dir):
    data_cuts_summary = MergeCutsSummaryTex()
    mc_cuts_summary = MergeCutsSummaryTex()
    for beam in dir_list:
        config = CompareCutsConfig(beam, target_dir)
        data_cuts_summary.append_summary(config, [0])
        mc_cuts_summary.append_summary(config, [1])
    data_cuts_summary.caption = config.data_caption
    mc_cuts_summary.caption = config.mc_caption
    data_cuts_summary.merge_summaries(target_dir+"/cuts_summary/data/", "data_cuts_summary.tex")
    mc_cuts_summary.merge_summaries(target_dir+"/cuts_summary/mc/", "mc_cuts_summary.tex")

def run_conglomerate(batch_level, config_list, dir_lists, do_cuts_summary, target_dir):
    rows = len(dir_lists)
    cols = min([len(sub_list) for sub_list in dir_lists])
    dir_list = []
    for sub_list in dir_lists:
        dir_list += sub_list
    if do_cuts_summary:
        cuts_summary(dir_list, target_dir)
    for ConfigClass in config_list:
        fail_list = []
        conglomerate_list = []
        print "Doing", ConfigClass.__name__
        for beam in dir_list:
            if batch_level > 5:
                ROOT.gROOT.SetBatch(False)
            try:
                config = ConfigClass(beam, target_dir)
                cong = ConglomerateContainer(config)
                cong.conglomerate()
                conglomerate_list.append(cong)
            except Exception:
                sys.excepthook(*sys.exc_info())
                fail_list.append(beam)
            ROOT.gROOT.SetBatch(True)
        if batch_level > 1:
            ROOT.gROOT.SetBatch(False)
        try:
            merge = ConglomerateMerge(conglomerate_list)
            merge.merge_all(rows, cols)
        except Exception:
            sys.excepthook(*sys.exc_info())
        print "Failed:"
        for fail in fail_list:
            print "    ", fail

def main(batch_level = 0):
    """
    Main program; 
    - batch_level tells how much output for ROOT: batch_level 0 is silent, 10 is most verbose
    """
    ROOT.gROOT.SetBatch(True)
    my_dir_list = [["2017-2.7_6-140_lH2_full"]]
    config_list = [CompareCutsConfig] #CompareData1DConfig]#, , CompareOpticsConfig, CompareData2DConfig]
    #config_list += [CompareAmplitudeConfigData, CompareAmplitudeConfigBoth]#, CompareAmplitudeConfigMC]
    target_dir = "output/2017-02/"
    batch_level = 10
    hide_root = False
    do_cuts_summary = False
    if batch_level < 10 and hide_root:
        ROOT.gErrorIgnoreLevel = 6000
    #run_conglomerate(batch_level, config_list, my_dir_list, do_cuts_summary, target_dir)
    my_dir_list = [
        ["2017-2.7_6-140_None", "2017-2.7_6-140_lH2_full", "2017-2.7_6-140_LiH",], #"2017-2.7_6-140_lH2_empty", 
        ["2017-2.7_10-140_None", "2017-2.7_10-140_lH2_full", "2017-2.7_10-140_LiH",], #"2017-2.7_10-140_lH2_empty", 
    ]
    config_list = [CompareAmplitudeConfigData]#, CompareAmplitudeConfigMC, CompareAmplitudeConfigBoth]
    run_conglomerate(batch_level, config_list, my_dir_list, do_cuts_summary, target_dir)

if __name__ == "__main__":
    main()
    if not ROOT.gROOT.IsBatch():
        raw_input("Finished - press <CR> to end")
