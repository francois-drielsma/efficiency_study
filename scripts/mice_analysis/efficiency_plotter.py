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
        self.get_cut_list()

    def get_cut_list(self):
        self.us_cut_list = self.config.cut_report[0]
        self.us_cut_list = set(self.us_cut_list)
        self.us_cut_list.discard("hline")
        self.us_cut_list.discard("all events")
        self.us_cut_list = list(self.us_cut_list)
        self.ds_cut_list = self.config.cut_report[1]
        self.ds_cut_list = set(self.ds_cut_list)
        self.ds_cut_list.discard("hline")
        self.ds_cut_list = list(self.ds_cut_list)

    def birth(self):
        self.set_plot_dir("efficiency_plots")
        self.get_data()
        self.setup_plots()
        self.process_plots()

    def process(self):
        self.get_data()
        self.process_plots()

    def death(self):
        for x_axis, y_axis, cut, reco_data_type in self.reco_plot_list:
            mc_data_type = {"reco_us":"mc_us", "reco_ds":"mc_ds"}[reco_data_type]
            self.finalise_reco_plot(x_axis, y_axis, cut, reco_data_type)
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
            self.setup_mc_plot(*arg, data_type="mc_ds")
            self.setup_mc_plot(*arg, data_type="mc_us")
        for arg in arg_list:
            for cut in self.ds_cut_list:
                self.setup_reco_plot(*arg, data_type="reco_ds", cut=cut)
                self.setup_ratio_plot(*arg, data_type="ratio_ds", cut=cut)
            for cut in self.us_cut_list:
                self.setup_reco_plot(*arg, data_type="reco_us", cut=cut)
                self.setup_ratio_plot(*arg, data_type="ratio_us", cut=cut)


    def process_plots(self):
        for x_axis, y_axis, data_type in self.mc_plot_list:
            self.process_mc_plot(x_axis, y_axis, data_type)
        for x_axis, y_axis, cut, data_type in self.reco_plot_list:
            self.process_reco_plot(x_axis, y_axis, cut, data_type)

    def plot_name(self, x_axis, y_axis, sample, cut= ""):
        plot_name = "efficiency_"+sample+"_"+x_axis+"_vs_"+y_axis
        if cut != "":
            plot_name += "_"+cut
        return plot_name

    def get_data_upstream(self, event, us_data, mc_data):
        station_us = self.config.mc_plots["mc_stations"]["tku_tp"][0]
        mc_hit = None
        for detector_hit in event["data"]:
            if detector_hit["detector"] == station_us:
                mc_hit = detector_hit["hit"]
                break
        if mc_hit == None:
            return
        if not event["mc_true_us_cut"]:
            for var in self.var_list:
                mc_data[var].append(mc_hit[var])
        for cut1 in self.us_cut_list:
            will_cut = False
            for cut2 in self.us_cut_list:
                if cut2 == cut1:
                    continue
                if cut2 in event["will_cut"] and event["will_cut"][cut2]:
                    will_cut = True
                    break
            # passes all cuts except cut1
            if not will_cut:
                for var in self.var_list:
                    us_data[var][cut1].append(mc_hit[var])

    def get_data_downstream(self, event, reco_data, mc_data):
        mc_hit = None
        station_ds = self.config.mc_plots["mc_stations"]["tkd_tp"][0]
        # out of ds sample
        if event["mc_true_ds_cut"]:
            return
        for detector_hit in event["data"]:
            if detector_hit["detector"] == station_ds:
                mc_hit = detector_hit["hit"]
                break
        for var in self.var_list:
            mc_data[var].append(mc_hit[var])
        for cut in self.ds_cut_list:
            if cut in event["will_cut"]:
                will_cut = event["will_cut"][cut]
            elif cut in event:
                will_cut = event[cut]
            else:
                raise KeyError("Did not recognise cut "+str(cut))
            if not will_cut:
                for var in self.var_list:
                    reco_data[var][cut].append(mc_hit[var])

    def get_data(self):
        self.data = {"mc_ds":{}, "mc_us":{}, "reco_ds":{}, "reco_us":{}}
        for var in self.var_list:
            self.data["mc_ds"][var] = []
            self.data["mc_us"][var] = []
            self.data["reco_ds"][var] = {}
            self.data["reco_us"][var] = {}
            for cut in self.ds_cut_list:
                self.data["reco_ds"][var][cut] = []
            for cut in self.us_cut_list:
                self.data["reco_us"][var][cut] = []
        for event in self.data_loader.events:
            self.get_data_upstream(event, self.data["reco_us"], self.data["mc_us"])
            self.get_data_downstream(event, self.data["reco_ds"], self.data["mc_ds"])

    def setup_reco_plot(self, x_axis, nx_bins, x_min_max, y_axis, ny_bins, y_min_max, cut, data_type):
        ratio_data_type = {"reco_us":"ratio_us", "reco_ds":"ratio_ds"}[data_type]
        self.reco_plot_list.append((x_axis, y_axis, cut, data_type))
        x_label = self.label_lookup[x_axis]
        y_label = self.label_lookup[y_axis]
        plot_name = self.plot_name(x_axis, y_axis, data_type, cut)
        hist = self.make_root_histogram(plot_name, plot_name+"_hist",
                                  [-1000], x_label, nx_bins,
                                  [-1000], y_label, ny_bins, [],
                                  x_min_max[0], x_min_max[1],
                                  y_min_max[0], y_min_max[1])
        hist.SetTitle("Sample: "+cut)
        hist.Draw("colz")

    def setup_ratio_plot(self, x_axis, nx_bins, x_min_max, y_axis, ny_bins, y_min_max, cut, data_type):
        ratio_plot_name = self.plot_name(x_axis, y_axis, data_type, cut)
        x_label = self.label_lookup[x_axis]
        y_label = self.label_lookup[y_axis]
        hist = self.make_root_histogram(ratio_plot_name, ratio_plot_name+"_hist",
                                  [-1000], x_label, nx_bins,
                                  [-1000], y_label, ny_bins, [],
                                  x_min_max[0], x_min_max[1],
                                  y_min_max[0], y_min_max[1])
        hist.SetTitle("Sample: "+cut)
        hist.Draw("colz")

    def setup_mc_plot(self, x_axis, nx_bins, x_min_max, y_axis, ny_bins, y_min_max, data_type):
        self.mc_plot_list.append((x_axis, y_axis, data_type))
        x_label = self.label_lookup[x_axis]
        y_label = self.label_lookup[y_axis]
        plot_name = self.plot_name(x_axis, y_axis, data_type)
        hist = self.make_root_histogram(plot_name, plot_name+"_hist",
                                  [-1000], x_label, nx_bins,
                                  [-1000], y_label, ny_bins, [],
                                  x_min_max[0], x_min_max[1],
                                  y_min_max[0], y_min_max[1])
        hist.Draw("colz")

    def process_reco_plot(self, x_axis, y_axis, cut, data_type):
        plot_name = self.plot_name(x_axis, y_axis, data_type, cut)
        hist = self.get_plot(plot_name)["histograms"][plot_name+"_hist"]
        x_data = self.data[data_type][x_axis][cut]
        y_data = self.data[data_type][y_axis][cut]
        for i, x in enumerate(x_data):
            y = y_data[i]
            hist.Fill(x, y)

    def process_mc_plot(self, x_axis, y_axis, data_type):
        plot_name = self.plot_name(x_axis, y_axis, data_type)
        hist = self.get_plot(plot_name)["histograms"][plot_name+"_hist"]
        x_data = self.data[data_type][x_axis]
        y_data = self.data[data_type][y_axis]
        for i, x in enumerate(x_data):
            y = y_data[i]
            hist.Fill(x, y)

    def finalise_reco_plot(self, x_axis, y_axis, cut, reco_data_type):
        if reco_data_type == "reco_us":
            mc_data_type = "mc_us"
            ratio_data_type = "ratio_us"
        else:
            mc_data_type = "mc_ds"
            ratio_data_type = "ratio_ds"
        truth_plot_name = self.plot_name(x_axis, y_axis, mc_data_type)
        truth_hist = self.get_plot(truth_plot_name)["histograms"][truth_plot_name+"_hist"]
        reco_plot_name = self.plot_name(x_axis, y_axis, reco_data_type, cut)
        reco_hist = self.get_plot(reco_plot_name)["histograms"][reco_plot_name+"_hist"]
        ratio_plot_name = self.plot_name(x_axis, y_axis, ratio_data_type, cut)
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
