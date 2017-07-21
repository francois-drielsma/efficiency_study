import sys
import operator
import xboa.common
import json
import copy
import math
import numpy

import cdb_tof_triggers_lookup
import ROOT
from xboa.bunch import Bunch
import scripts.utilities

from analysis_base import AnalysisBase

class DataPlotter(AnalysisBase):
    def __init__(self, config, config_anal, data_loader):
        super(DataPlotter, self).__init__(config, config_anal, data_loader)
        self.data_loader = data_loader
        self.max_tof01_bin = 0.
        self.max_tof12_bin = 0.
        self.max_p_bin = 0.
        self.run_numbers = set()
        self.will_cut_us = lambda event: event["upstream_cut"]
        self.will_cut_ds = lambda event: event["downstream_cut"]
        self.upstream_cut = {}
        self.downstream_cut = {}
        self.ellipse = {}

    def birth(self):
        self.run_numbers.update(self.data_loader.run_numbers)
        self.birth_tof("tof01")
        self.birth_tof("tof12")
        self.birth_pvalues("tku")
        self.birth_pvalues("tkd")
        self.birth_chi2("tku")
        self.birth_chi2("tkd")
        self.birth_ndf("tku")
        self.birth_ndf("tkd")
        self.birth_p_tot_vs_tof("tof01", "tku", "all")
        self.birth_p_tot_vs_tof("tof01", "tku", "us cut")
        self.birth_p_tot_vs_tof("tof12", "tkd", "all")
        self.birth_p_tot_vs_tof("tof12", "tkd", "ds cut")
        self.birth_p_tot_res()

        self.birth_var_2d("x", "tku", "px", "tku", [-150, 150], [-100, 100], "us cut")
        self.birth_var_2d("y", "tku", "py", "tku", [-150, 150], [-100, 100], "us cut")
        self.birth_var_2d("px", "tku", "py", "tku", [-100, 100], [-100, 100], "us cut")
        self.birth_var_2d("x", "tku", "y", "tku", [-150, 150], [-150, 150], "us cut")
        self.birth_var_2d("x", "tku", "y", "tku", [-150, 150], [-150, 150], "all")
        self.birth_var_2d("r", "tku", "pt", "tku", [0, 150*1.5], [0, 100*1.5], "all")

        self.birth_var_2d("x", "tkd", "px", "tkd", [-150, 150], [-100, 100], "ds cut")
        self.birth_var_2d("y", "tkd", "py", "tkd", [-150, 150], [-100, 100], "ds cut")
        self.birth_var_2d("px", "tkd", "py", "tkd", [-100, 100], [-100, 100], "ds cut")
        self.birth_var_2d("x", "tkd", "y", "tkd", [-150, 150], [-150, 150], "ds cut")
        self.birth_var_2d("x", "tkd", "y", "tkd", [-150, 150], [-150, 150], "all")
        self.birth_var_2d("r", "tkd", "pt", "tkd", [0, 150*1.5], [0, 100*1.5], "all")

        self.birth_var_1d("x", "tku")
        self.birth_var_1d("px", "tku")
        self.birth_var_1d("y", "tku")
        self.birth_var_1d("py", "tku")
        self.birth_var_1d("p", "tku") #, min_max=[80., 180.])
        self.birth_var_1d("r", "tku")
        self.birth_var_1d("SP Res(x)", "tku")
        self.birth_var_1d("SP Res(y)", "tku")

        self.birth_var_1d("r", "tku")
        self.birth_var_1d("x", "tkd")
        self.birth_var_1d("px", "tkd")
        self.birth_var_1d("y", "tkd")
        self.birth_var_1d("py", "tkd")
        self.birth_var_1d("p", "tkd") #, min_max=[80., 180.])
        self.birth_var_1d("r", "tkd")
        self.birth_var_1d("SP Res(x)", "tkd")
        self.birth_var_1d("SP Res(y)", "tkd")

        self.birth_ellipse("tku", "us cut")
        self.birth_ellipse("tku", "all")
        self.birth_ellipse("tkd", "ds cut")
        self.birth_ellipse("tkd", "all")

        self.birth_cuts_summary()

    def process(self):
        self.run_numbers.update(self.data_loader.run_numbers)
        self.process_tof("tof01")
        self.process_tof("tof12")
        self.process_tracker("tku", "pvalue")
        self.process_tracker("tkd", "pvalue")
        self.process_tracker("tku", "chi2")
        self.process_tracker("tkd", "chi2")
        self.process_tracker("tku", "ndf")
        self.process_tracker("tkd", "ndf")
        self.process_p_tot_vs_tof("tof01", "tku", "all")
        self.process_p_tot_vs_tof("tof01", "tku", "us cut")
        self.process_p_tot_vs_tof("tof12", "tkd", "all")
        self.process_p_tot_vs_tof("tof12", "tkd", "ds cut")
        self.process_p_tot_res()

        self.process_var_2d("x", "tku", "px", "tku", "us cut")
        self.process_var_2d("y", "tku", "py", "tku", "us cut")
        self.process_var_2d("px", "tku", "py", "tku", "us cut")
        self.process_var_2d("x", "tku", "y", "tku", "us cut")
        self.process_var_2d("x", "tku", "y", "tku", "all")
        self.process_var_2d("r", "tku", "pt", "tku", "all")

        self.process_var_2d("x", "tkd", "px", "tkd", "ds cut")
        self.process_var_2d("y", "tkd", "py", "tkd", "ds cut")
        self.process_var_2d("px", "tkd", "py", "tkd", "ds cut")
        self.process_var_2d("x", "tkd", "y", "tkd", "ds cut")
        self.process_var_2d("x", "tkd", "y", "tkd", "all")
        self.process_var_2d("r", "tkd", "pt", "tkd", "all")

        self.process_var_1d("x", "tku")
        self.process_var_1d("px", "tku")
        self.process_var_1d("y", "tku")
        self.process_var_1d("py", "tku")
        self.process_var_1d("p", "tku") #, min_max=[80., 180.])
        self.process_var_1d("r", "tku")
        self.process_var_1d("SP Res(x)", "tku")
        self.process_var_1d("SP Res(y)", "tku")


        self.process_var_1d("r", "tku")
        self.process_var_1d("x", "tkd")
        self.process_var_1d("px", "tkd")
        self.process_var_1d("y", "tkd")
        self.process_var_1d("py", "tkd")
        self.process_var_1d("p", "tkd") #, min_max=[80., 180.])
        self.process_var_1d("r", "tkd")
        self.process_var_1d("SP Res(x)", "tkd")
        self.process_var_1d("SP Res(y)", "tkd")

        self.process_ellipse("tku", "us cut")
        self.process_ellipse("tku", "all")
        self.process_ellipse("tkd", "ds cut")
        self.process_ellipse("tkd", "all")
        self.process_cuts_summary()

    def death(self):
        for plot_name in self.plots:
            min_value, max_value = None, None
            self.plots[plot_name]["canvas"].cd()
            for hist_name, hist in self.plots[plot_name]["histograms"].iteritems():
                if min_value == None:
                    min_value = hist.GetMinimum()/2.
                    max_value = hist.GetMaximum()*2.
                else:
                    min_value = min(hist.GetMinimum()/2., min_value)
                    max_value = max(hist.GetMaximum()*2., max_value)
            if min_value < 1.:
                min_value = 0.8
            hist_dict = self.get_plot(plot_name)["histograms"]
            for hist_name in sorted(hist_dict.keys()):
                hist.SetTitle(self.config_anal['name'])
                if len(hist_dict) != 3:
                    continue
                hist = hist_dict[hist_name]
                hist.GetYaxis().SetRangeUser(min_value, max_value)
                if "us cut" in hist_name:
                    hist.SetLineColor(2)
                    hist.Draw("SAME")
                elif "ds cut" in hist_name:
                    hist.SetLineColor(8)
                    hist.Draw("SAME")
                else:
                    hist.Draw()

        self.death_var_2d("x", "tku", "px", "tku", "us cut", True)
        self.death_var_2d("y", "tku", "py", "tku", "us cut", True)
        self.death_var_2d("px", "tku", "py", "tku", "us cut", True)
        self.death_var_2d("x", "tku", "y", "tku", "us cut", True)
        self.death_var_2d("x", "tku", "y", "tku", "all", True)

        self.death_var_2d("x", "tkd", "px", "tkd", "ds cut", True)
        self.death_var_2d("y", "tkd", "py", "tkd", "ds cut", True)
        self.death_var_2d("px", "tkd", "py", "tkd", "ds cut", True)
        self.death_var_2d("x", "tkd", "y", "tkd", "ds cut", True)
        self.death_var_2d("x", "tkd", "y", "tkd", "all", True)

        self.print_plots()
        self.death_wiki_summary()
        self.death_cuts_summary()


    def get_data_tof(self, tof):
        events = self.data_loader.events
        tof_us_cut = [event[tof] for event in events if not self.will_cut_us(event) and event[tof] != None]
        tof_ds_cut = [event[tof] for event in events if not self.will_cut_ds(event) and event[tof] != None]
        tof_all = [event[tof] for event in events if event[tof] != None]
        return tof_us_cut, tof_ds_cut, tof_all

    def birth_tof(self, tof):
        axis = {"tof01":"tof1 - tof0 [ns]", "tof12":"tof2 - tof1 [ns]"}[tof]
        tof_us_cut, tof_ds_cut, tof_all = self.get_data_tof(tof)
        xmin = min(25., self.config_anal[tof+"_cut_low"])
        xmax = max(45., self.config_anal[tof+"_cut_high"])
        hist = self.make_root_histogram(tof, tof+" all", tof_all, axis, 100, [], '', 0, [], xmin, xmax)
        hist = self.make_root_histogram(tof, tof+" us cut", tof_us_cut, axis, 100, [], '', 0, [], xmin, xmax)
        hist = self.make_root_histogram(tof, tof+" ds cut", tof_ds_cut, axis, 100, [], '', 0, [], xmin, xmax)
        self.get_plot(tof)["canvas"].SetLogy()

    def process_tof(self, tof):
        tof_us_cut, tof_ds_cut, tof_all = self.get_data_tof(tof)
        tof_hists = self.get_plot(tof)["histograms"]
        for data, hist_key in (tof_us_cut, "us cut"), (tof_ds_cut, "ds cut"), (tof_all, "all"):
            hist = tof_hists[tof+" "+hist_key]
            for item in data:
                hist.Fill(item)

    def get_tracker_data(self, tracker, var):
        data_all = []
        data_cut_us = []
        data_cut_ds = []
        for event in self.data_loader.events:
            if event[tracker] == None:
                continue
            for detector_hit in event["data"]:
                if detector_hit["detector"] != tracker+"_tp":
                    continue
                data = detector_hit[var]
                if var == "chi2":
                    data/float(detector_hit["ndf"])
                data_all.append(data)
                if not self.will_cut_us(event):
                    data_cut_us.append(data)
                if not self.will_cut_ds(event):
                    data_cut_ds.append(data)
        return data_cut_us, data_cut_ds, data_all

    def get_data_sp_residuals(self, tracker, axis):
        residuals_list = []
        residuals_cut_us_list = []
        residuals_cut_ds_list = []
        sp_detector = tracker+"_sp_"+str(self.config.tk_station)
        axis = axis.replace("SP Res(", "")
        axis = axis.replace(")", "")
        for event in self.data_loader.events:
            tp = event[tracker]
            if tp == None:
                continue # no track point recorded
            sp = None
            for point in event["data"]:
                if point["detector"] == sp_detector:
                    sp = point["hit"]
                    break
            if sp == None:
                continue # no space point recorded (e.g. due to dead channel)
            residuals_list.append(sp[axis] - tp[axis])        
            if not self.will_cut_us(event):
                residuals_cut_us_list.append(sp[axis] - tp[axis])        
            if not self.will_cut_ds(event):
                residuals_cut_ds_list.append(sp[axis] - tp[axis])        
        return residuals_cut_us_list, residuals_cut_ds_list, residuals_list

    def get_tracker_hit_data(self, tracker, var, require_both_trackers):
        if "SP Res" in var:
            return self.get_data_sp_residuals(tracker, var)
        data_all = []
        data_cut_us = []
        data_cut_ds = []
        other_tracker = {"tku":"tkd", "tkd":"tku"}[tracker]
        for event in self.data_loader.events:
            if event[tracker] == None:
                continue
            if require_both_trackers and event[other_tracker] == None:
                continue
            data = event[tracker][var]
            data_all.append(data)
            if not self.will_cut_us(event):
                data_cut_us.append(data)
            if not self.will_cut_ds(event):
                data_cut_ds.append(data)
        return data_cut_us, data_cut_ds, data_all

    def birth_pvalues(self, tracker):
        pvalues_cut_us, pvalues_cut_ds, pvalues_all = self.get_tracker_data(tracker, "pvalue")
        axis =  "P Value ("+tracker+")"
        hist = self.make_root_histogram("pvalue "+tracker, "pvalue all", pvalues_all, axis, 120, [], '', 0, [], -0.1, 1.1)
        hist_cut_us = self.make_root_histogram("pvalue "+tracker, "pvalue us cut", pvalues_cut_us, axis, 120, [], '', 0, [], -0.1, 1.1)
        hist_cut_ds = self.make_root_histogram("pvalue "+tracker, "pvalue ds cut", pvalues_cut_ds, axis, 120, [], '', 0, [], -0.1, 1.1)
        hist.Draw()
        hist_cut_us.Draw("SAME")
        hist_cut_ds.Draw("SAME")
        self.get_plot("pvalue "+tracker)["canvas"].SetLogy()

    def birth_chi2(self, tracker):
        chi2_cut_us, chi2_cut_ds, chi2_all = self.get_tracker_data(tracker, "chi2")
        axis =  "#chi^{2}/n_{df} ("+tracker+")"
        hist = self.make_root_histogram("chi2 "+tracker, "chi2 all", chi2_all, axis, 100, [], '', 0, [], 0., 20.)
        hist_cut_us = self.make_root_histogram("chi2 "+tracker, "chi2 us cut", chi2_cut_us, axis, 100, [], '', 0, [], 0., 20.)
        hist_cut_ds = self.make_root_histogram("chi2 "+tracker, "chi2 ds cut", chi2_cut_ds, axis, 100, [], '', 0, [], 0., 20.)
        hist.Draw()
        hist_cut_us.Draw("SAME")
        hist_cut_ds.Draw("SAME")
        self.get_plot("chi2 "+tracker)["canvas"].SetLogy()

    def birth_ndf(self, tracker):
        ndf_cut_us, ndf_cut_ds, ndf_all = self.get_tracker_data(tracker, "ndf")
        axis =  "n_{df} ("+tracker+")"
        hist = self.make_root_histogram("ndf "+tracker, "ndf all", ndf_all, axis, 100, [], '', 0, [], 0., 20.)
        hist_cut_us = self.make_root_histogram("ndf "+tracker, "ndf us cut", ndf_cut_us, axis, 100, [], '', 0, [], 0., 10.)
        hist_cut_ds = self.make_root_histogram("ndf "+tracker, "ndf ds cut", ndf_cut_ds, axis, 100, [], '', 0, [], 0., 10.)
        hist.Draw()
        hist_cut_us.Draw("SAME")
        hist_cut_ds.Draw("SAME")
        self.get_plot("ndf "+tracker)["canvas"].SetLogy()

    def process_tracker(self, tracker, key):
        data_cut_us, data_cut_ds, data_all = self.get_tracker_data(tracker, key)
        tracker_hists = self.get_plot(key+" "+tracker)["histograms"]
        for data, hist_key in (data_cut_us, "us cut"), (data_cut_ds, "ds cut"), (data_all, "all"):
            hist = tracker_hists[key+" "+hist_key]
            for item in data:
                hist.Fill(item)

    def get_p_tot_vs_tof_data(self, tof, tk, cuts):
        predicate_all = lambda event: event[tof] != None and event[tk] != None
        if cuts == "us cut":
            predicate = lambda event: predicate_all(event) and not self.will_cut_us(event)
        elif cuts == "ds cut":
            predicate = lambda event: predicate_all(event) and not self.will_cut_ds(event)
        elif cuts == "all":
            predicate = predicate_all

        tof_data = [event[tof] for event in self.data_loader.events if predicate(event)]
        p_data = [event[tk]["p"] for event in self.data_loader.events if predicate(event)]
        return tof_data, p_data

    def birth_p_tot_vs_tof(self, tof, tk, cuts):
        tof_data, p_data = self.get_p_tot_vs_tof_data(tof, tk, cuts)
        hist = self.make_root_histogram("p_"+tk+"_vs_"+tof+"_"+cuts,
                                        "p_tot_vs_tof", 
                                        tof_data, tof+" [ns]", 100, p_data,
                                        "p_{"+tk+"} [MeV/c]", 100, [],
                                        25., 45., 0., 300., )
        hist.Draw("COLZ")

    def process_p_tot_vs_tof(self, tof, tk, cuts):
        tof_data, p_data = self.get_p_tot_vs_tof_data(tof, tk, cuts)
        hist = self.get_plot("p_"+tk+"_vs_"+tof+"_"+cuts)["histograms"]
        hist = hist.values()[0]
        for i in range(len(tof_data)):
            hist.Fill(tof_data[i], p_data[i])

    def birth_p_tot_res(self):
        p_tku_cut_us, p_tku_cut_ds, p_tku_all = self.get_tracker_hit_data("tku", "p", True)
        p_tkd_cut_us, p_tkd_cut_ds, p_tkd_all = self.get_tracker_hit_data("tkd", "p", True)
        dp_cut_us = [p_tku - p_tkd_cut_us[i] for i, p_tku in enumerate(p_tku_cut_us)]
        dp_cut_ds = [p_tku - p_tkd_cut_ds[i] for i, p_tku in enumerate(p_tku_cut_ds)]
        dp_all = [p_tku - p_tkd_all[i] for i, p_tku in enumerate(p_tku_all)]

        axis = "p_{tku} - p_{tkd} [MeV/c]"
        hist = self.make_root_histogram("p_res", "p_res all", dp_all, axis, 100, [], '', 0, [], -50., 50.)
        hist_cut_us = self.make_root_histogram("p_res", "p_res us cut", dp_cut_us, axis, 100, [], '', 0, [], -50., 50.)
        hist_cut_ds = self.make_root_histogram("p_res", "p_res ds cut", dp_cut_ds, axis, 100, [], '', 0, [], -50., 50.)
        hist.Draw()
        hist_cut_us.Draw("SAME")
        hist_cut_ds.Draw("SAME")

        hist = self.make_root_histogram("p_res_vs_p_tku", "p_res_vs_p_tku", p_tku_cut_ds, "p_{tku} [MeV/c]", 50,
                                                          dp_cut_ds, "p_{tku} - p_{tkd} [MeV/c]", 50,
                                                          [], None, None, -50, 50)
        hist.Draw("COLZ")

    def process_p_tot_res(self):
        p_tku_cut_us, p_tku_cut_ds, p_tku_all = self.get_tracker_hit_data("tku", "p", True)
        p_tkd_cut_us, p_tkd_cut_ds, p_tkd_all = self.get_tracker_hit_data("tkd", "p", True)
        dp_cut_us = [p_tku - p_tkd_cut_us[i] for i, p_tku in enumerate(p_tku_cut_us)]
        dp_cut_ds = [p_tku - p_tkd_cut_ds[i] for i, p_tku in enumerate(p_tku_cut_ds)]
        dp_all = [p_tku - p_tkd_all[i] for i, p_tku in enumerate(p_tku_all)]

        p_res_hists = self.get_plot("p_res")["histograms"]
        for data, hist_key in (dp_cut_us, "us cut"), (dp_cut_ds, "ds cut"), (dp_all, "all"):
            hist = p_res_hists["p_res "+hist_key]
            for item in data:
                hist.Fill(item)

        p_res_vs_p_hist = self.get_plot("p_res_vs_p_tku")["histograms"]
        hist = p_res_vs_p_hist["p_res_vs_p_tku"]
        for i in range(len(p_tku_cut_ds)):
            hist.Fill(p_tku_cut_ds[i], dp_cut_ds[i])

    def birth_var_2d(self, var_1, us_ds_1, var_2, us_ds_2, min_max_1, min_max_2, cut):       
        data_cut_us, data_cut_ds, data_all = self.get_tracker_hit_data(us_ds_1, var_1, us_ds_1 == us_ds_2)
        data_1 = {"all":data_all, "us cut":data_cut_us, "ds cut":data_cut_ds}[cut]
        data_cut_us, data_cut_ds, data_all = self.get_tracker_hit_data(us_ds_2, var_2, us_ds_1 == us_ds_2)
        data_2 = {"all":data_all, "us cut":data_cut_us, "ds cut":data_cut_ds}[cut]
        name =  us_ds_1+"_"+var_1+"_"+ us_ds_2+"_"+var_2+"_"+cut
        lab_1 = us_ds_1+" "+var_1+" ["+scripts.utilities.default_units(var_1)+"]"
        lab_2 = us_ds_2+" "+var_2+" ["+scripts.utilities.default_units(var_2)+"]"
        hist = self.make_root_histogram(name, name,
                                        data_1, lab_1, 50,
                                        data_2, lab_2, 50, [],
                                        min_max_1[0], min_max_1[1],
                                        min_max_2[0], min_max_2[1])
        hist.Draw("COLZ")

    def process_var_2d(self, var_1, us_ds_1, var_2, us_ds_2, cut):
        data_cut_us, data_cut_ds, data_all = self.get_tracker_hit_data(us_ds_1, var_1, us_ds_1 == us_ds_2)
        data_1 = {"all":data_all, "us cut":data_cut_us, "ds cut":data_cut_ds}[cut]
        data_cut_us, data_cut_ds, data_all = self.get_tracker_hit_data(us_ds_2, var_2, us_ds_1 == us_ds_2)
        data_2 = {"all":data_all, "us cut":data_cut_us, "ds cut":data_cut_ds}[cut]
        name =  us_ds_1+"_"+var_1+"_"+ us_ds_2+"_"+var_2+"_"+cut
        hist = self.get_plot(name)["histograms"][name]
        for i, item_1 in enumerate(data_1):
            item_2 = data_2[i]
            hist.Fill(item_1, item_2)

    def death_var_2d(self, var_1, us_ds_1, var_2, us_ds_2, cut, build_ellipse):
        if not build_ellipse:
            return
        name =  us_ds_1+"_"+var_1+"_"+ us_ds_2+"_"+var_2+"_"+cut
        if us_ds_1 != us_ds_2:
            raise RuntimeError("We don't track correlations between upstream and downstream variables - ellipse not plotted")
        self.get_plot(name) # cd to the canvas
        index_1 = self.ellipse_variables.index(var_1)
        index_2 = self.ellipse_variables.index(var_2)
        mean_1 = self.ellipse["mean_"+us_ds_1+"_"+cut][index_1]
        mean_2 = self.ellipse["mean_"+us_ds_1+"_"+cut][index_2]
        cov_index = [[(index_1, index_1), (index_1, index_2)], [(index_1, index_2), (index_2, index_2)]]
        cov = self.ellipse["covariance_"+us_ds_1+"_"+cut]
        beam_cov = [[cov[i][j] for (i, j) in row] for row in cov_index]

        ellipse = self.make_ellipse_graph([mean_1, mean_2], beam_cov, 1)
        ellipse.Draw("SAMEL")

    @classmethod
    def make_ellipse_graph(cls, means, matrix, line_color):
        shell = xboa.common._common.make_shell(31, numpy.array(matrix))
        shell = [item.tolist()[0] for item in shell]
        shell = sorted(shell, key = lambda x: math.atan2(x[1], x[0]))
        x_list = [item[0]+means[0] for item in shell]
        y_list = [item[1]+means[1] for item in shell]
        hist, graph = xboa.common.make_root_graph("ellipse", x_list, "", y_list, "", sort = False)
        graph.SetLineColor(line_color)
        return graph

    def birth_ellipse(self, tracker, cut):
        n_var = len(self.ellipse_variables)
        self.ellipse["covariance_"+tracker+"_"+cut] = [[0. for j in range(n_var)] for i in range(n_var)]  
        self.ellipse["mean_"+tracker+"_"+cut] = [0. for i in range(n_var)]
        self.ellipse["nevents_"+tracker+"_"+cut] = 0
        self.process_ellipse(tracker, cut)

    def process_ellipse(self, tracker, cut):
        mean = self.ellipse["mean_"+tracker+"_"+cut]
        ellipse = self.ellipse["covariance_"+tracker+"_"+cut]
        n_events = self.ellipse["nevents_"+tracker+"_"+cut]
        n_var = len(self.ellipse_variables)
        # cut predicate the data
        if cut == "all":
            my_predicate = lambda event: event[tracker] == None
        elif cut == "us cut":
            my_predicate = lambda event: event[tracker] == None or self.will_cut_us(event)
        elif cut == "ds cut":
            my_predicate = lambda event: event[tracker] == None or self.will_cut_ds(event)
        # set up the matrix for this sample
        m_events = 0
        this_matrix = [[0. for j in range(n_var)] for i in range(n_var)]       
        this_mean = [0. for i in range(n_var)]
        for event in self.data_loader.events:
            if my_predicate(event):
                continue
            m_events += 1
            data = [event[tracker][var] for var in self.ellipse_variables]
            for i in range(n_var):
                this_mean[i] += data[i]
                for j in range(i, n_var):
                    this_matrix[i][j] += data[i]*data[j]

        # update the main ellipse
        for i in range(n_var):
            mean[i] = mean[i]*n_events/(n_events+m_events) + \
                      this_mean[i]/(n_events+m_events)
            for j in range(i, n_var):
                ellipse[i][j] = ellipse[i][j]*n_events/(n_events+m_events) + \
                                this_matrix[i][j]/(n_events+m_events)
                ellipse[j][i] = ellipse[i][j]
        self.ellipse["mean_"+tracker+"_"+cut] = mean
        self.ellipse["covariance_"+tracker+"_"+cut] = ellipse
        self.ellipse["nevents_"+tracker+"_"+cut] += m_events

    def birth_var_1d(self, var, us_ds):
        data_cut_us, data_cut_ds, data_all = self.get_tracker_hit_data(us_ds, var, False)
        name =  us_ds+"_"+var
        label = us_ds+" "+var+" ["+scripts.utilities.default_units(var)+"]"
        fit = scripts.utilities.fit_peak_data(data_all)
        mean = fit.GetParameter(1)
        sigma = fit.GetParameter(2)
        xmin, xmax = round(mean-sigma*5, 1), round(mean+sigma*5, 1)

        hist = self.make_root_histogram(name, name+" all",
                                        data_all, label, 50,
                                        [], "", 50, [], xmin, xmax)
        hist = self.make_root_histogram(name, name+" us cut",
                                        data_cut_us, label, 50,
                                        [], "", 50, [], xmin, xmax)
        hist = self.make_root_histogram(name, name+" ds cut",
                                        data_all, label, 50,
                                        [], "", 50, [], xmin, xmax)

    def process_var_1d(self, var, us_ds):       
        data_cut_us, data_cut_ds, data_all = self.get_tracker_hit_data(us_ds, var, False)
        name =  us_ds+"_"+var
        hist_dict = self.get_plot(name)["histograms"]
        for data, key in (data_all, "all"), (data_cut_ds, "ds cut"), (data_cut_us, "us cut"):
            hist_name = name+" "+key
            hist = hist_dict[hist_name]
            for i, item in enumerate(data):
                hist.Fill(item)

    def birth_cuts_summary(self):
        self.global_cut = {"upstream_cut":0, "downstream_cut":0, "all_events":0}
        self.upstream_cut = {}
        self.downstream_cut = {}
        for key in self.data_loader.events[0]["will_cut"]:
            self.upstream_cut[key] = 0
            self.downstream_cut[key] = 0

    def process_cuts_summary(self):
        for event in self.data_loader.events:
            self.global_cut["all_events"] += 1
            will_cut = event["will_cut"]
            for key in will_cut:
                if not will_cut[key]:
                    self.upstream_cut[key] += 1
                    if event["upstream_cut"]:
                        continue
                    self.downstream_cut[key] += 1
            if not event["upstream_cut"]:
                self.global_cut["upstream_cut"] += 1
            if not event["downstream_cut"]:
                self.global_cut["downstream_cut"] += 1

    def death_cuts_summary(self):
        print "========== cuts summary ============"
        for key in ["all_events", "upstream_cut", "downstream_cut"]:
            key_name = key.replace("_", " ")
            print "'"+key_name+":'", self.global_cut[key],
        print

        print "   ", "'cut name'".ljust(25), "us?".ljust(8), "ds?".ljust(8), "passed".ljust(8), "'upstream passed and passed'".ljust(8)
        for key in sorted(self.upstream_cut.keys()):
            key_name = "'"+key.replace("_", " ")+"'"
            is_active_us = self.config.upstream_cuts[key]
            is_active_ds = self.config.downstream_cuts[key]
            print "   ", key_name.ljust(25), str(is_active_us).ljust(8), str(is_active_ds).ljust(8), str(self.upstream_cut[key]).ljust(8), str(self.downstream_cut[key]).ljust(8)

    def death_wiki_summary(self):
        fout = open(self.plot_dir+"/wiki_summary.txt", "w")
        wiki_summary = "| "+self.config_anal['name']+" | "
        wiki_summary += " "+str(sorted(list(self.run_numbers)))+" |"
        try:
            cdb_dict = cdb_tof_triggers_lookup.parse_one_setting(cut_dict['runs'])
            cdb_dict["time"] = str(cdb_dict["time"][0])+" hrs "+str(cdb_dict["time"][1])+" mins"
            cdb_keys = ["lmc1234", "tof1_triggers", "tof2_triggers", "time"]
            for key in cdb_keys:
                wiki_summary += " "+str(cdb_dict[key]).ljust(8)+" |"
        except Exception:
            wiki_summary += " Failed to contact cdb |"
        for key in ["all_events", "upstream_cut", "downstream_cut"]:
            wiki_summary += " "+str(self.global_cut[key]).ljust(8)+" |"
        print >> fout, wiki_summary
        print wiki_summary

    ellipse_variables = ["x", "y", "px", "py"]


 
def do_plots(config, config_anal, data_loader):
    xboa.common.clear_root()
    plotter = DataPlotter(config, config_anal, data_loader.events, lambda event: event["upstream_cut"], lambda event: event["downstream_cut"], data_loader.run_numbers)
    sys.stdout.flush()

    plotter.plot_delta_tof("tof01")
    plotter.plot_delta_tof("tof12")

    plotter.bunch_plots("tku")

    plotter.bunch_plots("tkd")


