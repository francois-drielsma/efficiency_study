import sys
import operator
import xboa.common
import json
import copy
import math
import numpy

import ROOT
from xboa.bunch import Bunch
import utilities.cdb_tof_triggers_lookup as cdb_tof_triggers_lookup
import utilities.utilities
import utilities.cut_names

from analysis_base import AnalysisBase

class CutsPlotter(AnalysisBase):
    def __init__(self, config, config_anal, data_loader):
        super(CutsPlotter, self).__init__(config, config_anal, data_loader)
        self.data_loader = data_loader
        self.root_objects = []
        self.process_list = []
        self.cut_correlations = []
        self.cuts_report_full = {}


    def will_cut_except(self, wont_cut_list, will_cut_list, event):
        will_cut = False
        # require will_cut is FALSE if cut is on the wont_cut_list
        for a_cut in wont_cut_list:
            will_cut = will_cut or event["will_cut"][a_cut]
        # require will_cut is TRUE if cut is on the will_cut_list
        for a_cut in will_cut_list:
            will_cut = will_cut or not event["will_cut"][a_cut]
        return will_cut

    def birth(self):
        self.set_plot_dir("cut_plots")
        tof01_min = self.config_anal["tof01_cut_low"]
        tof01_max = self.config_anal["tof01_cut_high"]
        self.birth_cut_correlations()
        self.birth_var_1d("tof", "tof01", "us cut", ["tof01"], [1., 7.], 160, {}, [tof01_min, tof01_max])
        self.birth_var_1d("tof", "tof0_n_sp", "us cut", ["tof_0_sp"], [-0.5, 4.5], 5, {}, [0.5, 1.5])
        self.birth_var_1d("tof", "tof1_n_sp", "us cut", ["tof_1_sp"], [0.5, 4.5], 4, {}, [0.5, 1.5])
        self.birth_var_1d("tku", "n_tracks", "us cut", ["scifi_tracks_us", "global_through_tof0"], [-0.5, 4.5], 5, {}, [0.5, 1.5])
        # compare number of clusters in events with tracks to those without tracks; estimate impurity
        # number of clusters in all events in TKU that made a track
        self.birth_var_1d("tku", "scifi_n_planes_with_clusters", "us cut", [], [-0.5, 15.5], 16, {}, [])
        # number of clusters in all events in TKU with TOF (even those not making a track)
        self.birth_var_1d("tku", "scifi_n_planes_with_clusters", "us cut", ["scifi_tracks_us", "global_through_tof0"], [-0.5, 15.5], 16, {}, [])
        # number of clusters in TKU with TOF that DID NOT make a track
        # want to show that these events have few TKU clusters (hence impurity is unlikely)
        wont_cut_list = self.get_cut_list("us cut", ["scifi_tracks_us", "global_through_tof0"])
        will_cut_list = ["scifi_tracks_us"]
        self.birth_var_1d_cut_list("tku", "scifi_n_planes_with_clusters", wont_cut_list, will_cut_list, [-0.5, 15.5], 16, {}, [])
        # impurity end

        self.birth_var_1d("tku", "max_r", "us cut", ["scifi_fiducial_us"], [0., 300.], 100, {}, [150.])
        delta_tof01_min = self.config_anal["delta_tof01_lower"]
        delta_tof01_max = self.config_anal["delta_tof01_upper"]
        self.birth_var_1d("tof", "delta_tof01", "us cut", ["delta_tof01"], [-10., 5.], 100, {}, [delta_tof01_min, delta_tof01_max])

        p_min = min([min(a_bin) for a_bin in self.config_anal["p_bins"]])
        p_max = max([max(a_bin) for a_bin in self.config_anal["p_bins"]])
        self.birth_var_1d("tku", "p", "us cut", ["p_tot_us"], [p_min-35., p_max+35.], 100, {}, [p_min, p_max])
        chi2_max = self.config_anal["chi2_threshold"]
        self.birth_var_1d("tku", "chi2", "us cut", ["chi2_us"], [0., chi2_max*2], 100, {}, [chi2_max])
        diff = "global_through_virtual_diffuser"
        for aperture in self.config.upstream_aperture_cut:
            max_r = self.config.upstream_aperture_cut[aperture]
            if not self.config.upstream_cuts["upstream_aperture_cut"] or not self.config_anal["do_globals"]:
                continue
            self.birth_var_1d(aperture, "r", "us cut", ["upstream_aperture_cut"], [0., 200.], 100, {}, [max_r])
        p_min = self.config_anal["p_tot_ds_low"]
        p_max = self.config_anal["p_tot_ds_high"]
        delta_p = p_max-p_min
        self.birth_var_1d("tkd", "p", "ds cut", ["p_tot_ds"], [p_min-delta_p/2., p_max+delta_p/2], 100, {}, [p_min, p_max])
        self.birth_var_1d("tkd", "chi2", "ds cut", ["chi2_ds"], [0., chi2_max*2], 100, {}, [chi2_max])
        self.birth_var_1d("tkd", "n_tracks", "ds cut", ["scifi_tracks_ds"], [-0.5, 4.5], 5, {}, [0.5, 1.5])
        self.birth_var_1d("tkd", "max_r", "ds cut", ["scifi_fiducial_ds"], [0., 300.], 100, {}, [150.])
        # compare events with tracks to those without tracks; estimate noise impurity
        # number of clusters in all events in TKD that made a track
        self.birth_var_1d("tkd", "scifi_n_planes_with_clusters", "ds cut", [], [-0.5, 15.5], 16, {}, [])
        # number of TKD clusters in all events that did not make a TKD track
        # want to show that these events have few TKD clusters (hence impurity is unlikely)
        wont_cut_list = self.get_cut_list("us cut", [])
        will_cut_list = ["scifi_tracks_ds"]
        self.birth_var_1d_cut_list("tkd", "scifi_n_planes_with_clusters", wont_cut_list, will_cut_list, [-0.5, 15.5], 16, {}, [])

        # study efficiency by looking at events that make TOF2 space point but not corresponding TKD track
        # n clusters of all events passing upstream cuts, making tof2 sp but no tkd track
        wont_cut_list = self.get_cut_list("us cut", []) + ["tof_2_sp"]
        will_cut_list = ["scifi_tracks_ds"]
        self.birth_var_1d_cut_list("tkd", "scifi_n_planes_with_clusters", wont_cut_list, will_cut_list, [-0.5, 15.5], 16, {}, [])

        self.birth_cuts_summary()

    def process(self):
        for process in self.process_list:
            do_process = process["func"]
            args = process["process_args"]
            do_process(*args)
        self.process_cuts_summary()

    def death(self):
        self.normalise_cut_correlations()
        self.base_death()
        self.print_plots()
        self.death_cuts_summary()


    def get_n_hits(self, detector, event):
        return len([1 for hit in event["data"] if hit["detector"] == detector])

    def get_detector(self, detector, var, wont_cut_list, will_cut_list):
        data = []
        for event in self.data_loader.events:
            if self.will_cut_except(wont_cut_list, will_cut_list, event):
                continue
            for ev in event["data"]:
                if ev["detector"] == detector:
                    data.append(ev["hit"][var])
                    break
        return data

    def get_tof(self, detector, var, wont_cut_list, will_cut_list):
        data = []
        for event in self.data_loader.events:
            if self.will_cut_except(wont_cut_list, will_cut_list, event):
                continue
            if var in ("tof01", "tof12", "tof02", "delta_tof01", "delta_tof12"):
                data.append(event[var])
            elif "n_sp" in var:
                data.append(self.get_n_hits(var[0:4], event))
        return data

    def get_tracker(self, detector, var, wont_cut_list, will_cut_list):
        data = []
        for event in self.data_loader.events:
            if self.will_cut_except(wont_cut_list, will_cut_list, event):
                continue
            if var == "chi2":
                for hit in event["data"]:
                    if hit["detector"] == detector+"_tp":
                        data.append(hit["chi2"]/hit["ndf"])
                        break
            elif var == "n_tracks":
                tracker_key = {"tku":0, "tkd":1}[detector]
                data.append(event["scifi_n_tracks"][tracker_key])
            elif var == "scifi_n_planes_with_clusters":
                tracker_key = {"tku":0, "tkd":1}[detector]
                data.append(event["scifi_n_planes_with_clusters"][tracker_key])
            elif var == "max_r":
                max_r2 = 0.
                for hit in event["data"]:
                    if detector in hit["detector"] and "max_r2" in hit:
                        max_r2 = max(max_r2, hit["max_r2"])
                data.append(max_r2**0.5)
            else:
                data.append(event[detector][var])
        return data

    def get_cut_list(self, sample, cut_exception_list):
        cut_dict = {
            "us cut":self.config.upstream_cuts,
            "ds cut":self.config.downstream_cuts,
            "ex cut":self.config.extrapolation_cuts,
        }[sample]
        cut_list = [key for key in cut_dict if cut_dict[key]]
        cut_list = copy.deepcopy(cut_list)
        for cut_exception in cut_exception_list:
            try:
                cut_list.remove(cut_exception)
            except ValueError:
                print "WARNING:", cut_exception, "not found in sample", sample, "- IGNORING EXCEPTION!"
        return cut_list

    def get_data(self, detector, var, wont_cut_list, will_cut_list):
        try:
            getter = {
                "tof":self.get_tof,
                "tku":self.get_tracker,
                "tkd":self.get_tracker,
            }[detector]
        except KeyError:
            getter = self.get_detector

        data = getter(detector, var, wont_cut_list, will_cut_list)
        data = [value for value in data if value != None]
        return data

    def birth_cut_correlations(self):
        exclude = ("upstream_cut", "downstream_cut", "extrapolation_cut", "all events", "hline", "mc_true_us_cut", "mc_true_ds_cut")
        self.cut_correlations = []
        for a_report in self.config.cut_report:
            for item in a_report:
                if item not in exclude:
                    self.cut_correlations.append(item)
        n_bins = len(self.cut_correlations)
        name = "cut correlations"
        hist = ROOT.TH2D(name, "", n_bins, -0.5, n_bins-0.5, n_bins, -0.5, n_bins-0.5)
        for i, bin_name in enumerate(self.cut_correlations):
            hist.GetXaxis().SetBinLabel(i+1, bin_name)
            hist.GetYaxis().SetBinLabel(i+1, bin_name)
        self.get_plot(name)["histograms"][name] = hist
        hist.SetStats(False)
        hist.Draw("COLZ")
        self.get_plot(name)["config"]["background_fill"] = True
        self.process_cut_correlations()

    def process_cut_correlations(self):
        name = "cut correlations"
        hist = self.get_plot(name)["histograms"][name]
        for event in self.data_loader.events:
            for i1, item1 in enumerate(self.cut_correlations):
                for i2, item2 in enumerate(self.cut_correlations):
                    if event["will_cut"][item1] and event["will_cut"][item2]:
                        hist.Fill(i1, i2, 1)

    def normalise_cut_correlations(self):
        name = "cut correlations"
        hist = self.get_plot(name)["histograms"][name]
        bin_depth = []
        for i1 in range(len(self.cut_correlations)):
            bin_depth.append(hist.GetBinContent(i1, i1))
        for i1 in range(len(self.cut_correlations)):
            for i2 in range(len(self.cut_correlations)):
                a_bin = 1.*hist.GetBinContent(i1, i2)
                weight = 1.*min(bin_depth[i1], bin_depth[i2])
                new_content = 0.
                #if weight > 1e-9:
                #    new_content = a_bin/weight
                #hist.SetBinContent(i1, i2, new_content)
                print str(round(hist.GetBinContent(i1, i2))).rjust(5),
            print

    def get_cuts_box(self, cut_list):
        if len(cut_list) > 10:
            cut_list = ['...']+cut_list[-10:]
        y0 = 0.89 - 0.06*len(cut_list)
        text_box = ROOT.TPaveText(0.6, y0, 0.9, 0.89, "NDC")
        text_box.SetFillColor(0)
        text_box.SetBorderSize(0)
        text_box.SetTextSize(0.04)
        text_box.SetTextAlign(12)
        text_box.SetTextSize(0.03)

        text_box.AddText("Active cuts:")
        for cut in cut_list:
            text_box.AddText("  "+cut)
        text_box.Draw()
        return text_box

    def birth_var_1d(self, detector, var, sample, cut_exceptions, min_max, n_bins, options, verticals):
        wont_cut_list = self.get_cut_list(sample, cut_exceptions)
        will_cut_list = []
        self.birth_var_1d_cut_list(detector, var, wont_cut_list, will_cut_list, min_max, n_bins, options, verticals)

    def birth_var_1d_cut_list(self, detector, var, wont_cut_list, will_cut_list, min_max, n_bins, options, verticals):
        data = self.get_data(detector, var, wont_cut_list, will_cut_list)
        if len(data) == 0:
            if min_max[0] == None:
                for det in self.data_loader.detector_list():
                    print "   ", det
                raise RuntimeError("Failed to find var_1d data for "+detector+" "+var+" with no min set")
            else:
                data.append(min_max[0]-1.) # dummy underflow variable
        plot_name =  detector+"_"+var+"_"+str(len(wont_cut_list))+"_"+str(len(will_cut_list))
        units = utilities.utilities.default_units(var)
        label = var
        if units != '':
            label += ' ['+units+']'

        hist = self.make_root_histogram(plot_name, plot_name,
                                        data, label, n_bins,
                                        [], "", 50, [], min_max[0], min_max[1])
        hist.Draw("p e1")
        for x in verticals:
            hist, graph = self.make_root_graph(plot_name, plot_name+"_graph",
                      [x, x], "", [-1e9, 1e9], "", True,
                      None, None, None, None)
            graph.Draw("SAME")
        plot_config = self.get_plot(plot_name)["config"]
        plot_config["rescale"] = False
        plot_config["fit_1d_cuts"] = False
        plot_config["normalise"] = False
        plot_config["draw_1d_cuts"] = False
        #self.get_plot(plot_name)["misc"]["cuts_box"] = self.get_cuts_box(wont_cut_list)
        for key in options.keys():
            if key not in plot_config:
                raise KeyError("Did not recognise plot option "+str(key))
            plot_config[key] = options[key]
        self.process_list.append({
            "func":self.process_var_1d,
            "process_args":(plot_name, detector, var, wont_cut_list, will_cut_list),
        })

    def process_var_1d(self, plot_name, detector, var, wont_cut_list, will_cut_list):
        data = self.get_data(detector, var, wont_cut_list, will_cut_list)
        #if var == "max_r":
        #    print "cuts_plotter process_var_1d", data
        hist = self.get_plot(plot_name)["histograms"][plot_name]
        for item in data:
            hist.Fill(item)


    def birth_cuts_summary(self):
        self.global_cut = {"upstream_cut":0, "downstream_cut":0, "extrapolation_cut":0, "all_events":0}
        self.upstream_cut = {}
        self.downstream_cut = {}
        self.extrapolation_cut = {}
        self.cuts_report_full = []
        for key in self.data_loader.events[0]["will_cut"]:
            self.upstream_cut[key] = 0
            self.downstream_cut[key] = 0
            self.extrapolation_cut[key] = 0
        for a_report in self.config.cut_report:
            self.cuts_report_full.append({})
            for key in a_report:
                if key == "hline":
                    continue
                self.cuts_report_full[-1][key] = 0
        self.process_cuts_summary()

    def process_cuts_summary(self):
        for event in self.data_loader.events:
            for i, a_report in enumerate(self.config.cut_report):
                will_cut = False
                for key in a_report:
                    if key == "hline":
                        continue
                    elif key == "all events":
                        self.cuts_report_full[0][key] += 1
                        continue
                    if will_cut:
                        continue
                    if key in event["will_cut"]:
                        if event["will_cut"][key]:
                            will_cut = True
                    elif key in event:
                        # e.g. tof01 is a float, also a cut
                        if type(event[key]) == type(True) and event[key]:
                            will_cut = True
                    if not will_cut:
                        self.cuts_report_full[i][key] += 1


    def death_cuts_summary(self):
        for i, a_report in enumerate(self.config.cut_report):
            file_name = self.plot_dir+"/cuts_summary_"+str(i)+".tex"
            fout = open(file_name, "w")
            print >> fout, ("cut").ljust(20), "&", self.config_anal["name"], "//"
            for key in a_report:
                if key == "hline":
                    print >> fout, "\hline"
                    continue
                print >> fout, key.ljust(20), "&", str(self.cuts_report_full[i][key]), "//"
            fout.close()
            fout = open(file_name, "r")
            print fout.read()
            print

    cut_names = utilities.cut_names.cut_names

