import sys
import operator
import xboa.common
import json
import copy
import math
import numpy

import utilities.cdb_tof_triggers_lookup as cdb_tof_triggers_lookup
import utilities.utilities as utilities
import ROOT
from xboa.bunch import Bunch

from analysis_base import AnalysisBase

class GlobalsPlotter(AnalysisBase):
    def __init__(self, config, config_anal, data_loader):
        super(GlobalsPlotter, self).__init__(config, config_anal, data_loader)
        self.data_loader = data_loader
        self.config = config
        self.config_anal = config_anal
        self.root_objects = []
        self.ds_counter = 0
        self.ex_counter = 0

    def will_cut_us(self, event):
        return event["upstream_cut"]

    def will_cut_ds(self, event):
        return event["downstream_cut"]

    def will_cut_ex(self, event):
        return event["extrapolation_cut"]

    def birth(self):
        self.set_plot_dir("global_plots")
        self.birth_residuals_1d("tkd_tp", "x", "global_through", -100., 200., 500)
        self.birth_residuals_1d("tkd_tp", "y", "global_through", -100., 200., 500)
        self.birth_residuals_1d("tkd_tp", "px", "global_through", -50., 100., 500)
        self.birth_residuals_1d("tkd_tp", "py", "global_through", -50., 100., 500)
        self.birth_residuals_1d("tkd_tp", "pz", "global_through", -50., 100., 500)
        self.birth_residuals_1d("tkd_tp", "p", "global_through", -50., 100., 500)
        self.birth_residuals_1d("tof0", "t", "global_through", -100., 100., 500)
        self.birth_residuals_1d("tof1", "x", "global_through", -400., 800., 500)
        self.birth_residuals_1d("tof1", "y", "global_through", -400., 800., 500)
        self.birth_residuals_1d("tof1", "t", "global_through", -1., 2.0, 500)
        self.birth_residuals_1d("tof2", "x", "global_through", -200., 400., 500)
        self.birth_residuals_1d("tof2", "y", "global_through", -200., 400., 500)
        self.birth_residuals_1d("tof2", "t", "global_through", -100., 100., 500)
        self.birth_residuals_1d("tof2", "x", "global_ds", -200., 400., 500)
        self.birth_residuals_1d("tof2", "y", "global_ds", -200., 400., 500)
        self.birth_residuals_1d("tof2", "t", "global_ds", -100., 100., 500)
        #self.birth_misses_2d("tkd_tp", "x", "y", "us cut", xmin=-300, xmax=300, nbinsx=40, ymin=-300, ymax=300., nbinsy=40)
        #self.birth_misses_2d("tof2", "x", "y", "ex cut", xmin=-400, xmax=400, nbinsx=40, ymin=-400, ymax=400., nbinsy=40)
        
        self.birth_event_display(3, False, ["global_through_tof1", "global_through_tkd_tp"], "ds cut")

    def process(self):
        self.process_residuals_1d("tkd_tp", "x", "global_through")
        self.process_residuals_1d("tkd_tp", "y", "global_through")
        self.process_residuals_1d("tkd_tp", "px", "global_through")
        self.process_residuals_1d("tkd_tp", "py", "global_through")
        self.process_residuals_1d("tkd_tp", "pz", "global_through")
        self.process_residuals_1d("tkd_tp", "p", "global_through")
        self.process_residuals_1d("tof0", "t", "global_through")
        self.process_residuals_1d("tof1", "x", "global_through")
        self.process_residuals_1d("tof1", "y", "global_through")
        self.process_residuals_1d("tof1", "t", "global_through")
        self.process_residuals_1d("tof2", "x", "global_through")
        self.process_residuals_1d("tof2", "y", "global_through")
        self.process_residuals_1d("tof2", "t", "global_through")
        self.process_residuals_1d("tof2", "x", "global_ds")
        self.process_residuals_1d("tof2", "y", "global_ds")
        self.process_residuals_1d("tof2", "t", "global_ds")
        #self.process_misses_2d("tkd_tp", "x", "y", "us cut")
        #self.process_misses_2d("tof2", "x", "y", "ex cut")

    def get_graph_min_max(self, plot_name):
        plot_dict = self.plots[plot_name]["graphs"]
        x_min_value, y_min_value = None, None
        for graph_name, graph in plot_dict.iteritems():
            if x_min_value == None:
                x_min_value = graph.GetXaxis().GetXmin()
                x_max_value = graph.GetXaxis().GetXmax()
            else:
                x_min_value = min(graph.GetXaxis().GetXmin(), x_min_value)
                x_max_value = max(graph.GetXaxis().GetXmax(), x_max_value)
        for graph_name, graph in plot_dict.iteritems():
            if y_min_value == None:
                y_min_value = graph.GetYaxis().GetXmin()
                y_max_value = graph.GetYaxis().GetXmax()
            else:
                y_min_value = min(graph.GetYaxis().GetXmin(), y_min_value)
                y_max_value = max(graph.GetYaxis().GetXmax(), y_max_value)
        return [x_min_value, x_max_value], [y_min_value, y_max_value]


    def death_graph_generic(self):
        for plot_name in self.plots:
            if len(self.get_plot(plot_name)["graphs"]) == 0:
                continue # it is a histogram
            x_min_max, y_min_max = self.get_graph_min_max(plot_name)
            self.plots[plot_name]["canvas"].cd()
            hist_dict = self.get_plot(plot_name)["histograms"]
            for hist_name in sorted(hist_dict.keys()):
                hist = hist_dict[hist_name]
                hist.SetTitle(self.config_anal['name'])
                hist.GetXaxis().SetRangeUser(x_min_max[0], x_min_max[1])
                hist.GetYaxis().SetRangeUser(y_min_max[0], y_min_max[1])
                hist.Draw()
                break
            graph_dict = self.get_plot(plot_name)["graphs"]
            for graph_name in sorted(graph_dict.keys()):
                graph = graph_dict[graph_name]
                if "global_" in graph_name:
                    graph.Draw("SAMEL")
                else:
                    graph.Draw("SAMEP")

    def death(self):
        self.base_death()
        self.death_graph_generic()
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

    def get_residual_data(self, detector, var, prefix):
        data_all = []
        data_cut_us = []
        data_cut_ds = []
        data_cut_ex = []
        n_det_hits, n_global_hits = 0, 0
        global_name = prefix+"_"+detector #self.detector_name_to_virtual_name(detector)
        global_offset = 0.
        if var == "t" and detector == "tof0":
            global_offset = utilities.electron_tof("tof0", "tof1", self.config) # note sign
        if var == "t" and detector == "tof2":
            global_offset = utilities.electron_tof("tof2", "tof1", self.config)
        for event in self.data_loader.events:
            det_hit, global_hit = None, None
            for hit in event["data"]:
                if hit["detector"] == detector:
                    det_hit = hit
                elif hit["detector"] == global_name:
                    global_hit = hit
                if det_hit != None and global_hit != None:
                    break
            if det_hit != None:
                n_det_hits += 1
            if global_hit != None:
                n_global_hits += 1
            if det_hit == None or global_hit == None:
                # either we missed the detector or we didn't get a global extrapolation
                continue
            data = det_hit["hit"][var] - global_hit["hit"][var] - global_offset

            data_all.append(data)
            if not self.will_cut_us(event):
                data_cut_us.append(data)
            if not self.will_cut_ds(event):
                data_cut_ds.append(data)
            if not self.will_cut_ex(event):
                data_cut_ex.append(data)
        if len(data_all) == 0:
            print "WARNING - No data in get_residual_data(", detector, var, prefix, ")"
            print "number of det hits", n_det_hits, "global hits", n_global_hits
        return data_cut_us, data_cut_ds, data_cut_ex, data_all

    def get_miss_data(self, detector, var):
        data_all = []
        data_cut_us = []
        data_cut_ds = []
        data_cut_ex = []

        n_det_hits, n_global_hits = 0, 0
        global_name = "global_through_virtual_"+detector
        for event in self.data_loader.events:
            det_hit, global_hit = None, None
            for hit in event["data"]:
                if hit["detector"] == detector:
                    det_hit = hit
                elif hit["detector"] == global_name:
                    global_hit = hit
                if det_hit != None and global_hit != None:
                    break
            if det_hit != None:
                n_det_hits += 1
            if global_hit != None:
                n_global_hits += 1
            if det_hit != None or global_hit == None:
                # either we hit the detector or we didn't get a global extrapolation
                continue
            data = global_hit["hit"][var]
            data_all.append(data)
            if not self.will_cut_us(event):
                data_cut_us.append(data)
            if not self.will_cut_ds(event):
                data_cut_ds.append(data)
            if not self.will_cut_ex(event):
                data_cut_ex.append(data)
        if len(data_all) == 0:
            print self.data_loader.detector_list()
            print "Get miss data failed: number of det hits", n_det_hits, "global hits", n_global_hits
            raise RuntimeError("No data")
        return data_cut_us, data_cut_ds, data_cut_ex, data_all

    def birth_event_display(self, n_events, aggregate, required_detectors, cut = None):
        data = []
        n_events_done = 1
        for event in self.data_loader.events:
            this_required_detectors = copy.deepcopy(required_detectors)
            for hit in event["data"]:
                for detector in this_required_detectors:
                    if detector == hit["detector"] and hit["hit"]["p"] > 1.:
                        this_required_detectors.remove(detector)
                        break
            if len(this_required_detectors) > 0:
                continue
            if cut == "us cut" and self.will_cut_us(event):
                continue
            elif cut == "ds cut" and self.will_cut_ds(event):
                continue
            elif cut == "ex cut" and self.will_cut_ex(event):
                continue
            canvas_name = "event_display_"+str(n_events_done)
            if aggregate:
                canvas_name = "event_display"
            #print "EVENT DISPLAY ##", n_events_done, "##"
            for hit in event["data"]:
                detector = hit["detector"]
                hit = hit["hit"]
                #print "       ", hit["x"], hit["y"], hit["z"], ";", hit["t"], "**", \
                #    hit["px"], hit["py"], hit["pz"], ";", hit["energy"], "**", \
                #    detector
            for hit_type, marker, color in [("global_through", 24, 2)]:#, ("global_ds", 25, 4,), ("reco", 26, 8)]:
                for var in ["x", "y", "energy", "t"]:
                    if var == "energy":
                        var_predicate = lambda hit: hit["hit"][var] > 106. and abs(hit["hit"]["z"]) > 1.
                    else:
                        var_predicate = lambda hit: abs(hit["hit"]["z"]) > 1.
                    if "global" in hit_type:
                        predicate = lambda hit: var_predicate(hit) and hit_type in hit["detector"]
                    else:
                        predicate = lambda hit: var_predicate(hit) and "global" not in hit["detector"]
                    z_list = [hit["hit"]["z"] for hit in event["data"] if predicate(hit)]
                    var_list = [hit["hit"][var] for hit in event["data"] if predicate(hit)]
                    if len(z_list) == 0:
                        continue
                    #print hit_type, var
                    #print "    z:  ", z_list
                    #print "    var:", var_list
                    #print
                    name = str(n_events_done)+"_"+var
                    units = {"energy":"MeV", "x":"mm", "y":"mm", "t":"ns"}[var]
                    hist, graph = self.make_root_graph(canvas_name+"_"+var, name+"_"+hit_type,
                                  z_list, "z [mm]", var_list, var+" ["+units+"]", True,
                                  None, None, None, None)
                    graph.SetMarkerStyle(marker)
                    graph.SetLineColor(color)
                    graph.SetMarkerColor(marker)
                    if aggregate:
                        graph.SetMarkerColor(2*n_events_done)
            n_events_done += 1
            if n_events_done > n_events:
                break

    def birth_residuals_1d(self, detector, var, prefix, xmin= None, xmax = None, nbins = None):
        res_us_cut, res_ds_cut, res_ex_cut, res_all = self.get_residual_data(detector, var, prefix)
        units = utilities.default_units(var)
        if units != "":
            units = " ["+units+"]"

        fit = utilities.fit_peak_data(res_all)
        mean = fit.GetParameter(1)
        sigma = fit.GetParameter(2)
        if xmin == None:
            xmin = round(mean-sigma*5, 1)
        if xmax == None:
            xmax = round(mean+sigma*5, 1)
        axis =  detector+" "+var+"(Reco - "+prefix+")"+units
        name = prefix+" residual "+detector+" "+var
        if nbins == None:
            nbins = 120
        hist = self.make_root_histogram(name, "res all", res_all, axis, nbins, [], '', 0, [], xmin, xmax)
        hist_us_cut = self.make_root_histogram(name, "res us cut", res_us_cut, axis, nbins, [], '', 0, [], xmin, xmax)
        hist_ds_cut = self.make_root_histogram(name, "res ds cut", res_ds_cut, axis, nbins, [], '', 0, [], xmin, xmax)
        hist_ex_cut = self.make_root_histogram(name, "res ex cut", res_ex_cut, axis, nbins, [], '', 0, [], xmin, xmax)
        hist.Draw()
        hist_us_cut.Draw("SAME")
        hist_ds_cut.Draw("SAME")
        hist_ex_cut.Draw("SAME")
        self.get_plot(name)["canvas"].SetLogy()
        self.get_plot(name)["config"]["rescale"] = True
        self.get_plot(name)["config"]["fit_1d_cuts"] = True
        self.get_plot(name)["config"]["normalise"] = True
        self.get_plot(name)["config"]["draw_1d_cuts"] = True

    def process_residuals_1d(self, detector, var, prefix):
        res_us_cut, res_ds_cut, res_ex_cut, res_all = self.get_residual_data(detector, var, prefix)
        name = prefix+" residual "+detector+" "+var
        res_hists = self.get_plot(name)["histograms"]
        for data, hist_key in (res_us_cut, "us cut"), (res_ds_cut, "ds cut"), (res_ex_cut, "ex cut"), (res_all, "all"):
            hist = res_hists["res "+hist_key]
            for item in data:
                hist.Fill(item)

    def birth_misses_2d(self, detector, var_1, var_2, cut, xmin= None, xmax = None, nbinsx = None, ymin = None, ymax = None, nbinsy = None):
        res_us_cut, res_ds_cut, res_ex_cut, res_all = self.get_miss_data(detector, var_1)
        res_1 = {"us cut":res_us_cut, "ds cut":res_ds_cut, "ex cut":res_ex_cut, "all":res_all}[cut]
        units_1 = xboa.hit.Hit.default_units()[var_1]
        if units_1 != "":
            units_1 = " ["+units_1+"]"

        res_us_cut, res_ds_cut, res_ex_cut, res_all = self.get_miss_data(detector, var_2)
        res_2 = {"us cut":res_us_cut, "ds cut":res_ds_cut, "ex cut":res_ex_cut, "all":res_all}[cut]
        units_2 = xboa.hit.Hit.default_units()[var_2]
        if units_2 != "":
            units_2 = " ["+units_2+"]"
        name = "misses "+detector+" "+var_1+" vs "+var_2+" "+cut

        if nbinsx == None:
            nbinsx = 50
        if nbinsy == None:
            nbinsy = 50
        axis_1 =  detector+" "+var_1+" (Extrap) "+units_1
        axis_2 =  detector+" "+var_2+" (Extrap) "+units_2
        if len(res_1) == 0:
            print "Failed to find any residuals for cut", cut
            print "    all:   ", len(res_all) 
            print "    us cut:", len(res_us_cut) 
            print "    ds cut:", len(res_ds_cut) 
            print "    ex cut:", len(res_ex_cut) 
            res_1, res_2 = [0.], [0.]
        hist = self.make_root_histogram(name, name,
                                        res_1, axis_1, nbinsx,
                                        res_2, axis_2, nbinsy, [],
                                        xmin, xmax, ymin, ymax)
        hist.Draw("COLZ")

    def process_misses_2d(self, detector, var_1, var_2, cut):
        res_us_cut, res_ds_cut, res_ex_cut, res_all = self.get_miss_data(detector, var_1)
        res_1 = {"us cut":res_us_cut, "ds cut":res_ds_cut, "ex cut":res_ex_cut, "all":res_all}[cut]
        units_1 = xboa.hit.Hit.default_units()[var_1]
        if units_1 != "":
            units_1 = " ["+units_1+"]"

        res_us_cut, res_ds_cut, res_ex_cut, res_all = self.get_miss_data(detector, var_2)
        res_2 = {"us cut":res_us_cut, "ds cut":res_ds_cut, "ex cut":res_ex_cut, "all":res_all}[cut]
        units_2 = xboa.hit.Hit.default_units()[var_2]
        if units_2 != "":
            units_2 = " ["+units_2+"]"
        name = "misses "+detector+" "+var_1+" vs "+var_2+" "+cut

        res_hist = self.get_plot(name)["histograms"][name]
        for i, r1 in enumerate(res_1):
            res_hist.Fill(r1, res_2[i])

