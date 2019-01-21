import xboa.common
import ROOT

import utilities.root_style

from mice_analysis.analysis_base import AnalysisBase


class EfficiencyPlotter(AnalysisBase):
    def __init__(self, config, config_anal, data_loader):
        super(EfficiencyPlotter, self).__init__(config, config_anal, data_loader)
        self.config = config
        self.config_anal = config_anal
        self.data_loader = data_loader
        self.reco_plot_list = []
        self.mc_plot_list = []
        self.data = {}
        self.get_reco_cut_list()

    def get_reco_cut_list(self):
        self.reco_cut_list = self.config.cut_report[1]
        self.reco_cut_list = set(self.reco_cut_list)
        self.reco_cut_list.discard("hline")
        self.reco_cut_list = list(self.reco_cut_list)

    def birth(self):
        self.set_plot_dir("efficiency_plots")
        self.get_data()
        self.setup_plots()
        self.process_plots()

    def process(self):
        self.get_data()
        self.process_plots()

    def death(self):
        for x_axis, y_axis, cut in self.reco_plot_list:
            self.finalise_reco_plot(x_axis, y_axis, cut)
        self.print_plots()

    def setup_plots(self):
        arg_list = [
            ("x",  20, [-200., 200.], "y",  20, [-200., 200.]),
            ("x",  20, [-200., 200.], "px", 20, [-100., 100.]),
            ("x",  20, [-200., 200.], "py", 20, [-100., 100.]),
            ("y",  20, [-200., 200.], "px", 20, [-100., 100.]),
            ("y",  20, [-200., 200.], "py", 20, [-100., 100.]),
            ("px", 20, [-100., 100.], "py", 20, [-100., 100.]),
        ]
        for arg in arg_list:
            self.setup_mc_plot(*arg)
        for cut in self.reco_cut_list:
            for arg in arg_list:
                self.setup_reco_plot(*arg, cut=cut)

    def process_plots(self):
        for x_axis, y_axis in self.mc_plot_list:
            self.process_mc_plot(x_axis, y_axis)
        for x_axis, y_axis, cut in self.reco_plot_list:
            self.process_reco_plot(x_axis, y_axis, cut)

    def death_plots(self):
        for x_axis, y_axis, cut in self.reco_plot_list:
            self.finalise_reco_plot(x_axis, y_axis, cut)

    def plot_name(self, x_axis, y_axis, sample, cut= ""):
        plot_name = "efficiency_"+sample+"_"+x_axis+"_vs_"+y_axis
        if cut != "":
            plot_name += "_"+cut
        return plot_name

    def get_data(self):
        self.data = {"mc":{}, "reco":{}}
        for var in self.var_list:
            self.data["mc"][var] = []
            self.data["reco"][var] = {}
            for cut in self.reco_cut_list:
                self.data["reco"][var][cut] = []
        station_ds = self.config.mc_plots["mc_stations"]["tkd_tp"][0]
        for event in self.data_loader.events:
            mc_hit = None
            # out of ds sample
            if event["mc_true_ds_cut"]:
                continue
            for detector_hit in event["data"]:
                if detector_hit["detector"] == station_ds:
                    mc_hit = detector_hit["hit"]
                    break
            for var in self.var_list:
                self.data["mc"][var].append(mc_hit[var])
            for cut in self.reco_cut_list:
                if cut in event["will_cut"]:
                    will_cut = event["will_cut"][cut]
                elif cut in event:
                    will_cut = event[cut]
                else:
                    raise KeyError("Did not recognise cut "+str(cut))
                if not will_cut:
                    for var in self.var_list:
                        self.data["reco"][var][cut].append(mc_hit[var])

    def setup_reco_plot(self, x_axis, nx_bins, x_min_max, y_axis, ny_bins, y_min_max, cut):
        self.reco_plot_list.append((x_axis, y_axis, cut))
        x_label = self.label_lookup[x_axis]
        y_label = self.label_lookup[y_axis]
        plot_name = self.plot_name(x_axis, y_axis, "reco", cut)
        hist = self.make_root_histogram(plot_name, plot_name+"_hist",
                                  [-1000], x_label, nx_bins,
                                  [-1000], y_label, ny_bins, [],
                                  x_min_max[0], x_min_max[1],
                                  y_min_max[0], y_min_max[1])
        ratio_plot_name = self.plot_name(x_axis, y_axis, "ratio", cut)
        hist.SetTitle("Sample: "+cut)
        hist.Draw("colz")
        hist = self.make_root_histogram(ratio_plot_name, ratio_plot_name+"_hist",
                                  [-1000], x_label, nx_bins,
                                  [-1000], y_label, ny_bins, [],
                                  x_min_max[0], x_min_max[1],
                                  y_min_max[0], y_min_max[1])
        hist.Draw("colz")

    def setup_mc_plot(self, x_axis, nx_bins, x_min_max, y_axis, ny_bins, y_min_max):
        self.mc_plot_list.append((x_axis, y_axis))
        x_label = self.label_lookup[x_axis]
        y_label = self.label_lookup[y_axis]
        plot_name = self.plot_name(x_axis, y_axis, "mc")
        hist = self.make_root_histogram(plot_name, plot_name+"_hist",
                                  [-1000], x_label, nx_bins,
                                  [-1000], y_label, ny_bins, [],
                                  x_min_max[0], x_min_max[1],
                                  y_min_max[0], y_min_max[1])
        hist.Draw("colz")

    def process_reco_plot(self, x_axis, y_axis, cut):
        plot_name = self.plot_name(x_axis, y_axis, "reco", cut)
        hist = self.get_plot(plot_name)["histograms"][plot_name+"_hist"]
        x_data = self.data["reco"][x_axis][cut]
        y_data = self.data["reco"][y_axis][cut]
        for i, x in enumerate(x_data):
            y = y_data[i]
            hist.Fill(x, y)

    def process_mc_plot(self, x_axis, y_axis):
        plot_name = self.plot_name(x_axis, y_axis, "mc")
        hist = self.get_plot(plot_name)["histograms"][plot_name+"_hist"]
        x_data = self.data["mc"][x_axis]
        y_data = self.data["mc"][y_axis]
        for i, x in enumerate(x_data):
            y = y_data[i]
            hist.Fill(x, y)

    def finalise_reco_plot(self, x_axis, y_axis, cut):
        truth_plot_name = self.plot_name(x_axis, y_axis, "mc")
        truth_hist = self.get_plot(truth_plot_name)["histograms"][truth_plot_name+"_hist"]
        reco_plot_name = self.plot_name(x_axis, y_axis, "reco", cut)
        reco_hist = self.get_plot(reco_plot_name)["histograms"][reco_plot_name+"_hist"]
        ratio_plot_name = self.plot_name(x_axis, y_axis, "ratio", cut)
        ratio_hist = self.get_plot(ratio_plot_name)["histograms"][ratio_plot_name+"_hist"]
        for i in range(truth_hist.GetNbinsX()):
            for j in range(truth_hist.GetNbinsY()):
                reco_count = reco_hist.GetBinContent(i, j)
                truth_count = truth_hist.GetBinContent(i, j)
                if truth_count > 0:
                    ratio = max(reco_count/truth_count, 1e-9)
                    ratio_hist.SetBinContent(i, j, ratio)
        ratio_hist.SetTitle(ratio_plot_name)
        ratio_hist.Draw("COLZ")

    var_list = ["x", "y", "px", "py", "pz"]
    label_lookup = {
        "x":"x [mm]",
        "y":"y [mm]",
        "px":"p_{x} [MeV/c]",
        "py":"p_{y} [MeV/c]",
        "pz":"p_{z} [MeV/c]",
    }
    root_objects = []
