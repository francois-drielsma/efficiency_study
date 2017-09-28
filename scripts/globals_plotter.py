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

class GlobalsPlotter(AnalysisBase):
    def __init__(self, config, config_anal, data_loader):
        super(GlobalsPlotter, self).__init__(config, config_anal, data_loader)
        self.data_loader = data_loader
        self.config = config
        self.config_anal = config_anal
        self.root_objects = []

    def will_cut_us(self, event):
        return event["upstream_cut"]

    def will_cut_ds(self, event):
        return event["downstream_cut"]

    def birth(self):
        self.birth_residuals_1d("tkd_tp", "x")
        self.birth_residuals_1d("tkd_tp", "y")
        self.birth_residuals_1d("tkd_tp", "px")
        self.birth_residuals_1d("tkd_tp", "py")
        self.birth_residuals_1d("tkd_tp", "pz")
        self.birth_residuals_1d("tof1", "x")
        self.birth_residuals_1d("tof1", "y")
        self.birth_residuals_1d("tof2", "x")
        self.birth_residuals_1d("tof2", "y")

    def process(self):
        self.process_residuals_1d("tkd_tp", "x")
        self.process_residuals_1d("tkd_tp", "y")
        self.process_residuals_1d("tkd_tp", "px")
        self.process_residuals_1d("tkd_tp", "py")
        self.process_residuals_1d("tkd_tp", "pz")
        self.process_residuals_1d("tof1", "x")
        self.process_residuals_1d("tof1", "y")
        self.process_residuals_1d("tof2", "x")
        self.process_residuals_1d("tof2", "y")

    def get_text_box(self, fit_0, fit_1, fit_2, hist_0, hist_1, hist_2):
        text_box = ROOT.TPaveText(0.6, 0.4, 0.9, 0.9, "NDC")
        text_box.SetFillColor(0)
        text_box.SetBorderSize(0)
        text_box.SetTextSize(0.04)
        text_box.SetTextAlign(12)
        text_box.SetTextSize(0.03)

        text_box.AddText("All (Black)")
        text_box.AddText("  Number:    "+str(hist_0.GetEntries()))
        text_box.AddText("  Mean:        "+str(round(fit_0.GetParameter(1), 2)))
        text_box.AddText("  Std:           "+str(round(fit_0.GetParameter(2), 2)))
        text_box.AddText("US Cut (Red)")
        text_box.AddText("  Number:    "+str(hist_1.GetEntries()))
        text_box.AddText("  Mean:        "+str(round(fit_1.GetParameter(1), 2)))
        text_box.AddText("  Std:           "+str(round(fit_1.GetParameter(2), 2)))
        text_box.AddText("DS Cut (Green)")
        text_box.AddText("  Number:    "+str(hist_2.GetEntries()))
        text_box.AddText("  Mean:        "+str(round(fit_2.GetParameter(1), 2)))
        text_box.AddText("  Std:           "+str(round(fit_2.GetParameter(2), 2)))
        text_box.SetBorderSize(1)
        text_box.Draw()
        self.root_objects.append(text_box)
        return text_box


    def death_hist_generic(self):
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
                    fit_1 = scripts.utilities.fit_peak(hist, 2, "Q", "SAME")
                    hist.Draw("SAME")
                    hist_1 = hist
                elif "ds cut" in hist_name:
                    hist.SetLineColor(8)
                    fit_2 = scripts.utilities.fit_peak(hist, 2, "Q", "SAME")
                    hist.Draw("SAME")
                    hist_2 = hist
                else:
                    fit_0 = scripts.utilities.fit_peak(hist, 2, "Q", "SAME")
                    hist.Draw()
                    hist_0 = hist
            box = self.get_text_box(fit_0, fit_1, fit_2, hist_0, hist_1, hist_2)


    def death(self):
        self.death_hist_generic()
        self.print_plots()

    def detector_name_to_virtual_name(self, detector):
        #detector_dict = {"tof2":"tof2", "tkd_tp":"tkd_tp"}
        #detector = detector_dict[detector]
        for detector_z, dummy, det_name in self.config.detectors:
            if det_name == detector:
                break
        if det_name != detector:
            raise RuntimeError("Couldnt find "+str(detector)+" in detectors")
        for virtual_z, dummy, virtual_name in self.config.virtual_detectors:
            if abs(virtual_z - detector_z) < 3.:
                break
        if abs(virtual_z - detector_z) > 3.:
            raise RuntimeError("Couldnt find virtual detector with z near "+str(detector_z))
        return virtual_name

    def get_residual_data(self, detector, var):
        data_all = []
        data_cut_us = []
        data_cut_ds = []
        n_det_hits, n_virt_hits = 0, 0
        virtual_name = self.detector_name_to_virtual_name(detector)
        for event in self.data_loader.events:
            det_hit, virt_hit = None, None
            for hit in event["data"]:
                if hit["detector"] == detector:
                    det_hit = hit
                elif hit["detector"] == virtual_name:
                    virt_hit = hit
                if det_hit != None and virt_hit != None:
                    break
            if det_hit != None:
                n_det_hits += 1
            if virt_hit != None:
                n_virt_hits += 1
            if det_hit == None or virt_hit == None:
                continue
            data = det_hit["hit"][var] - virt_hit["hit"][var]
            data_all.append(data)
            if not self.will_cut_us(event):
                data_cut_us.append(data)
            if not self.will_cut_ds(event):
                data_cut_ds.append(data)
        if len(data_all) == 0:
            print self.data_loader.detector_list()
            print "number of det hits", n_det_hits, "virt hits", n_virt_hits
            raise RuntimeError("No data")
        return data_cut_us, data_cut_ds, data_all

    def birth_residuals_1d(self, detector, var):
        res_us_cut, res_ds_cut, res_all = self.get_residual_data(detector, var)
        units = xboa.hit.Hit.default_units()[var]
        if units != "":
            units = " ["+units+"]"

        fit = scripts.utilities.fit_peak_data(res_all)
        mean = fit.GetParameter(1)
        sigma = fit.GetParameter(2)
        xmin, xmax = round(mean-sigma*5, 1), round(mean+sigma*5, 1)
        axis =  detector+" "+var+units
        name = "residual "+detector+" "+var
        hist = self.make_root_histogram(name, "res all", res_all, axis, 120, [], '', 0, [], xmin, xmax)
        hist_us_cut = self.make_root_histogram(name, "res us cut", res_us_cut, axis, 120, [], '', 0, [], xmin, xmax)
        hist_ds_cut = self.make_root_histogram(name, "res ds cut", res_ds_cut, axis, 120, [], '', 0, [], xmin, xmax)
        hist.Draw()
        hist_us_cut.Draw("SAME")
        hist_ds_cut.Draw("SAME")
        self.get_plot(name)["canvas"].SetLogy()

    def process_residuals_1d(self, detector, var):
        res_us_cut, res_ds_cut, res_all = self.get_residual_data(detector, var)
        name = "residual "+detector+" "+var
        res_hists = self.get_plot(name)["histograms"]
        for data, hist_key in (res_us_cut, "us cut"), (res_ds_cut, "ds cut"), (res_all, "all"):
            hist = res_hists["res "+hist_key]
            for item in data:
                hist.Fill(item)

 
def do_plots(config, config_anal, data_loader):
    xboa.common.clear_root()
    plotter = DataPlotter(config, config_anal, data_loader.events, lambda event: event["upstream_cut"], lambda event: event["downstream_cut"], data_loader.run_numbers)
    sys.stdout.flush()

    plotter.plot_delta_tof("tof01")
    plotter.plot_delta_tof("tof12")

    plotter.bunch_plots("tku")

    plotter.bunch_plots("tkd")


