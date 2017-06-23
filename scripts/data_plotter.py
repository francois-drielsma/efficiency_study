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

class DataPlotter(object):
    def __init__(self, config, config_anal, events, will_cut_us, will_cut_ds, run_numbers = None):
        self.events = events
        self.will_cut_us = will_cut_us
        self.will_cut_ds = will_cut_ds
        self.config = config
        self.config_analysis = config_anal
        self.plot_dir = config_anal["plot_dir"]
        self.max_tof12_bin = None
        self.max_p_bin = None
        self.run_numbers = run_numbers
        for event in events:
            if self.will_cut_us(event) and event["tku"] != None:
                event["tku"]["local_weight"] = 0.
            if self.will_cut_ds(event) and event["tkd"] != None:
                event["tkd"]["local_weight"] = 0.

    def choose_data(self, key, cuts):
        if cuts == None:
            predicate = lambda event: event[key] != None
        else:
            predicate = lambda event: not cuts(event) and event[key] != None
        data = [event[key] for event in self.events if predicate(event)]
        return data

    def get_cuts_summary(self):
        # there are a set of cuts that are to do with basic "data quality" e.g.
        # exactly one hit in TOF, etc and there are a set of cuts that are to do
        # with "bias correction"/etc.
        # I report the number of events that pass each data quality cut;
        # and the number of events that pass each "bias correction" cut AND the data quality cuts
        accepted_dict = {"upstream_cut":0, "data_cut":None, "downstream_cut":0, "all_events":0}
        accepted_upstream_dict = {"upstream_cut":0, "data_cut":0, "downstream_cut":0, "all_events":0}
        for key in self.events[0]["will_cut"]:
            accepted_dict[key] = 0
            accepted_upstream_dict[key] = 0
        for event in self.events:
            accepted_dict["all_events"]  += 1 # total number of events considered
            accepted_upstream_dict["all_events"]  += 1 # total number of events considered
            will_cut = event["will_cut"]
            for key in will_cut:
                if not will_cut[key]:
                    accepted_dict[key] += 1
                    if event["upstream_cut"]:
                        continue
                    accepted_upstream_dict[key] += 1
            if not event["upstream_cut"]:
                accepted_dict["upstream_cut"]  += 1 # total number accepted by all cuts
            if not event["upstream_cut"]:
                accepted_upstream_dict["upstream_cut"]  += 1 # total number accepted by data cuts
                if not event["upstream_cut"]:
                    accepted_upstream_dict["upstream_cut"]  += 1 # total number accepted by data cuts
                if not event["downstream_cut"]:
                    accepted_upstream_dict["downstream_cut"]  += 1 # total number accepted by data cuts
            if not event["downstream_cut"]:
                accepted_dict["downstream_cut"]  += 1 # total number accepted by all cuts
        return accepted_dict, accepted_upstream_dict

    def print_cuts_summary(self):
        cut_dict, cut_data_dict = self.get_cuts_summary()
        print "========== cuts summary ============"
        print "   ", "cut name".ljust(25), "us?".ljust(5), "ds?".ljust(5), "pass".ljust(8), "pass and pass data cut".ljust(8)
        for key in sorted(cut_dict.keys()):
            if key in self.config.upstream_cuts:
                is_active_us = self.config.upstream_cuts[key]
            else:
                is_active_us = None
            if key in self.config.downstream_cuts:
                is_active_ds = self.config.downstream_cuts[key]
            else:
                is_active_ds = None
            print "   ", key.ljust(25), str(is_active_us).ljust(5), str(is_active_ds).ljust(5), str(cut_dict[key]).ljust(8), str(cut_data_dict[key]).ljust(8)

    def print_wiki_summary(self):
        cut_dict, cut_data_dict = self.get_cuts_summary()
        
        cut_dict['data_cut'] = cut_data_dict['data_cut']
        cut_dict['runs'] = sorted([run for run in self.run_numbers])
        cut_dict['max_p_bin'] = 0.
        cut_dict['max_tof01_bin'] = 0.
        try:
            cut_dict['max_p_bin'] = round(self.max_p_bin, 2)
            cut_dict['max_tof01_bin'] = round(self.max_tof01_bin, 2)
        except TypeError:
            pass
        cut_keys = ['tof01', 'tof_2_sp', 'data_cut', 'max_p_bin', 'max_tof01_bin']
        wiki_summary = "| "+self.config_analysis['name']+" | "
        try:
            cdb_dict = cdb_tof_triggers_lookup.parse_one_setting(cut_dict['runs'])
            cdb_dict["time"] = str(cdb_dict["time"][0])+" hrs "+str(cdb_dict["time"][1])+" mins"
            cdb_keys = ["runs", "lmc1234", "tof1_triggers", "tof2_triggers", "time"]
            for key in cdb_keys:
                wiki_summary += " "+str(cdb_dict[key]).ljust(8)+" |"
        except Exception:
            wiki_summary += "| Failed to contact cdb |"
        for key in cut_keys:
            wiki_summary += " "+str(cut_dict[key]).ljust(8)+" |"
        return wiki_summary

    def plot_pvalues(self):
        for tracker in "tku", "tkd":
            pvalue_canvas = xboa.common.make_root_canvas(tracker+" PValues")
            pvalues_all = []
            pvalues_cut_us = []
            pvalues_cut_ds = []
            for event in self.events:
                if event[tracker] == None:
                    continue
                for detector_hit in event["data"]:
                    if detector_hit["detector"] != tracker+"_tp":
                        continue
                    pvalues_all.append(detector_hit["pvalue"])
                    if not self.will_cut_us(event):
                        pvalues_cut_us.append(detector_hit["pvalue"])
                    if not self.will_cut_ds(event):
                        pvalues_cut_ds.append(detector_hit["pvalue"])
            hist = xboa.common.make_root_histogram("pvalues "+tracker, pvalues_all, "P Value ("+tracker+")", 120, xmin = -0.1, xmax = 1.1)
            hist.SetTitle(self.config_analysis['name'])
            hist_cut_us = xboa.common.make_root_histogram("pvalues "+tracker, pvalues_cut_us, "P Value ("+tracker+")", 120, xmin = -0.1, xmax = 1.1)
            hist_cut_us.SetLineColor(2)
            hist_cut_us.SetTitle(self.config_analysis['name'])
            hist_cut_ds = xboa.common.make_root_histogram("pvalues "+tracker, pvalues_cut_us, "P Value ("+tracker+")", 120, xmin = -0.1, xmax = 1.1)
            hist_cut_ds.SetLineColor(8)
            hist_cut_ds.SetTitle(self.config_analysis['name'])

            min_value = max(hist_cut_us.GetMinimum()/2., 0.8)
            max_value = hist.GetMaximum()*2.
            hist.GetYaxis().SetRangeUser(min_value, max_value)
            hist.Draw()
            hist_cut_us.Draw("SAME")
            hist_cut_ds.Draw("SAME")
            pvalue_canvas.SetLogy()
            pvalue_canvas.Update()
            for format in ["png", "root", "eps"]:
                pvalue_canvas.Print(self.plot_dir+"pvalue_"+tracker+"."+format)


    def plot_delta_tof(self, tofs):
        dtof01_us_cut = [event["delta_"+tofs] for event in self.events if not self.will_cut_us(event)]
        dtof01_us_cut = [dtof01 for dtof01 in dtof01_us_cut if dtof01 != None]
        dtof01_ds_cut = [event["delta_"+tofs] for event in self.events if not self.will_cut_ds(event)]
        dtof01_ds_cut = [dtof01 for dtof01 in dtof01_ds_cut if dtof01 != None]
        dtof01_all = [event["delta_"+tofs] for event in self.events]
        dtof01_all = [dtof01 for dtof01 in dtof01_all if dtof01 != None]
        if len(dtof01_all) == 0:
            print "No delta tof plot - perhaps extrapolation is disabled?"
            return
        dtof01_canvas = xboa.common.make_root_canvas("delta_tof01_canvas")
        dtof01_canvas.SetLogy()
        xmin = -15.
        xmax = 10.
        hist = xboa.common.make_root_histogram("delta "+tofs, dtof01_all, "Reco "+tofs+" - Extrapolated "+tofs+" [ns]", 100, xmin = xmin, xmax = xmax)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw()
        hist = xboa.common.make_root_histogram("delta "+tofs+" in cut", dtof01_us_cut, "Reco "+tofs+" - Extrapolated "+tofs+" [ns]", 100, xmin = xmin, xmax = xmax)
        hist.SetLineColor(2)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("SAME")
        hist = xboa.common.make_root_histogram("delta "+tofs+" in cut", dtof01_ds_cut, "Reco "+tofs+" - Extrapolated "+tofs+" [ns]", 100, xmin = xmin, xmax = xmax)
        hist.SetLineColor(8)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("SAME")
        dtof01_canvas.Update()
        for format in ["png", "root", "eps"]:
            dtof01_canvas.Print(self.plot_dir+"delta_"+tofs+"."+format)

    def plot_tof01(self):
        tof01_us_cut = [event["tof01"] for event in self.events if not self.will_cut_us(event)]
        tof01_ds_cut = [event["tof01"] for event in self.events if not self.will_cut_ds(event)]
        tof01_all = [event["tof01"] for event in self.events]
        tof01_all = [tof01 for tof01 in tof01_all if tof01 != None]
        tof01_canvas = xboa.common.make_root_canvas("tof01_canvas")
        tof01_canvas.SetLogy()
        xmin = min(25., self.config_analysis["tof01_cut_low"])
        xmax = max(45., self.config_analysis["tof01_cut_high"])
        hist = xboa.common.make_root_histogram("tof01", tof01_all, "tof1 - tof0 [ns]", 100, xmin = xmin, xmax = xmax)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw()
        hist = xboa.common.make_root_histogram("tof01 in cut", tof01_us_cut, "tof1 - tof0 [ns]", 100, xmin = xmin, xmax = xmax)
        hist.SetLineColor(2)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("SAME")
        hist = xboa.common.make_root_histogram("tof01 in cut", tof01_ds_cut, "tof1 - tof0 [ns]", 100, xmin = xmin, xmax = xmax)
        hist.SetLineColor(8)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("SAME")
        tof01_canvas.Update()
        self.max_tof01_bin = self.get_max(hist)
        for format in ["png", "root", "eps"]:
            tof01_canvas.Print(self.plot_dir+"tof01."+format)

    def plot_tof_position(self):
        for event in self.events:
            for hit in event:
                pass

    def plot_tof12(self):
        tof12_us_cut = [event["tof12"] for event in self.events if not self.will_cut_us(event) and event["tof12"] != None and event["tof12"] < 100.]
        tof12_ds_cut = [event["tof12"] for event in self.events if not self.will_cut_ds(event) and event["tof12"] != None and event["tof12"] < 100.]
        tof12_all = [event["tof12"] for event in self.events if event["tof12"] != None and event["tof12"] < 100.]
        tof12_canvas = xboa.common.make_root_canvas("tof12_canvas")
        tof12_canvas.SetLogy()
        xmin = min(25., self.config_analysis["tof12_cut_low"])
        xmax = max(45., self.config_analysis["tof12_cut_high"])
        hist = xboa.common.make_root_histogram("tof12", tof12_all, "tof2 - tof1 [ns]", 100, xmin = xmin, xmax = xmax)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw()
        hist = xboa.common.make_root_histogram("tof12 in cut", tof12_us_cut, "tof2 - tof1 [ns]", 100, xmin = xmin, xmax = xmax)
        hist.SetLineColor(2)
        hist.Draw("SAME")
        tof12_canvas.Update()
        hist = xboa.common.make_root_histogram("tof12 in cut", tof12_ds_cut, "tof2 - tof1 [ns]", 100, xmin = xmin, xmax = xmax)
        hist.SetLineColor(8)
        hist.Draw("SAME")
        tof12_canvas.Update()
        self.max_tof12_bin = max([hist.GetBinContent(i) for i in range(102)])
        for format in ["png", "root", "eps"]:
            tof12_canvas.Print(self.plot_dir+"tof12."+format)

    def plot_p_tot_vs_tof(self, tof, tk):
        canvas = xboa.common.make_root_canvas("p_tot vs tof no cuts")
        tof_no_cut = [event[tof] for event in self.events if event[tof] != None and event[tk] != None and event["tof12"] < 100.]
        p_tot_no_cut = [event[tk]["p"] for event in self.events if event[tof] != None and event[tk] != None and event["tof12"] < 100.]

        hist = xboa.common.make_root_histogram("Events in cut", tof_no_cut, tof+" [ns]", 100, p_tot_no_cut, "p_{"+tk+"} [MeV/c]", 100, ymin=0., ymax=300., xmin=25., xmax=45.)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("COLZ")
        canvas.Update()
        for format in ["png", "root", "eps"]:
            canvas.Print(self.plot_dir+tk+"_p_tot_vs_"+tof+"_all."+format)

        canvas = xboa.common.make_root_canvas("p_tot vs tof us cuts")
        tof_no_cut = [event[tof] for event in self.events if event[tof] != None and event[tk] != None and not self.will_cut_us(event) and event["tof12"] < 100.]
        p_tot_no_cut = [event[tk]["p"] for event in self.events if event[tof] != None and event[tk] != None and not self.will_cut_us(event) and event["tof12"] < 100.]

        hist = xboa.common.make_root_histogram("Events in cut", tof_no_cut, tof+" [ns]", 100, p_tot_no_cut, "p_{"+tk+"} [MeV/c]", 100)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("COLZ")
        canvas.Update()


        canvas = xboa.common.make_root_canvas("p_tot vs tof ds cuts")
        tof_no_cut = [event[tof] for event in self.events if event[tof] != None and event[tk] != None and not self.will_cut_ds(event) and event["tof12"] < 100.]
        p_tot_no_cut = [event[tk]["p"] for event in self.events if event[tof] != None and event[tk] != None and not self.will_cut_ds(event) and event["tof12"] < 100.]

        hist = xboa.common.make_root_histogram("Events in cut", tof_no_cut, tof+" [ns]", 100, p_tot_no_cut, "p_{"+tk+"} [MeV/c]", 100)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("COLZ")
        canvas.Update()

        for format in ["png", "root", "eps"]:
            canvas.Print(self.plot_dir+tk+"_p_tot_vs_"+tof+"_cuts."+format)

    def plot_p_tot_res(self):
        canvas = xboa.common.make_root_canvas("p_tot us vs ds")
        cut_pred = lambda event: event["tku"] != None and event["tkd"] != None
        p_tot_us_no_cut = [event["tku"]["p"] for event in self.events if cut_pred(event)]
        p_tot_ds_no_cut = [event["tkd"]["p"] for event in self.events if cut_pred(event)]
        p_tot_res_no_cut = [p_tot - p_tot_ds_no_cut[i] for i, p_tot in enumerate(p_tot_us_no_cut)]

        cut_pred = lambda event: event["tku"] != None and event["tkd"] != None and not self.will_cut_us(event)
        p_tot_us_cut_us = [event["tku"]["p"] for event in self.events if cut_pred(event)]
        p_tot_ds_cut_us = [event["tkd"]["p"] for event in self.events if cut_pred(event)]
        p_tot_res_cut_us = [p_tot - p_tot_ds_cut_us[i] for i, p_tot in enumerate(p_tot_us_cut_us)]
        #mean = numpy.mean(p_tot_res_cut)
        #std = numpy.std(p_tot_res_cut)

        cut_pred = lambda event: event["tku"] != None and event["tkd"] != None and not self.will_cut_ds(event)
        p_tot_us_cut_ds = [event["tku"]["p"] for event in self.events if cut_pred(event)]
        p_tot_ds_cut_ds = [event["tkd"]["p"] for event in self.events if cut_pred(event)]
        p_tot_res_cut_ds = [p_tot - p_tot_ds_cut_ds[i] for i, p_tot in enumerate(p_tot_us_cut_ds)]

        hist = xboa.common.make_root_histogram("All events", p_tot_res_no_cut, "p_{tku} - p_{tkd} [MeV/c]", 100, xmin=-50, xmax=50)
        hist.SetTitle(self.config_analysis['name'])
        #hist.SetTitle("mean: "+str(mean)+" std: "+str(std))
        hist.Draw()
        hist = xboa.common.make_root_histogram("Events in cut", p_tot_res_cut_us, "p_{tku} - p_{tkd} [MeV/c]", 100, xmin=-50, xmax=50)
        #hist.SetTitle("mean: "+str(mean)+" std: "+str(std))
        hist.SetLineColor(2)
        hist.Draw("SAME")
        hist = xboa.common.make_root_histogram("Events in cut", p_tot_res_cut_ds, "p_{tku} - p_{tkd} [MeV/c]", 100, xmin=-50, xmax=50)
        #hist.SetTitle("mean: "+str(mean)+" std: "+str(std))
        hist.SetLineColor(8)
        hist.Draw("SAME")
        canvas.Update()
        for format in ["png", "root", "eps"]:
            canvas.Print(self.plot_dir+"p_tot_res."+format)

        canvas = xboa.common.make_root_canvas("delta p vs p_tku")
        hist = xboa.common.make_root_histogram("Events in cut", p_tot_us_cut_ds, "p_{tku} [MeV/c]", 50,
                                                                p_tot_res_cut_ds, "p_{tku} - p_{tkd} [MeV/c]", 50,
                                                                ymin=-50, ymax=50)
        hist.Draw("COLZ")
        canvas.Update()
        for format in ["png", "root", "eps"]:
            canvas.Print(self.plot_dir+"delta_p_vs_p_tku."+format)

    def plot_analysis_b_field(self):
        pass

    # MAKE THIS A CONT PLOT AND SUPERIMPOSE EACH TOF BIN
    def plot_var_2d(self, var_1, us_ds_1, var_2, us_ds_2, min_max_1 = [None, None], min_max_2 = [None, None], cuts = -1):
        if cuts == -1:
            cuts =  self.will_cut_us
        cuts_string = {self.will_cut_us:"upstream", self.will_cut_ds:"downstream", None:"none"}
        name =  us_ds_1+"_"+var_1+"_"+ us_ds_2+"_"+var_2+"_cuts-"+cuts_string[cuts]
        canvas = xboa.common.make_root_canvas(name)
        data_1 = self.choose_data(us_ds_1, cuts)
        data_2 = self.choose_data(us_ds_2, cuts)
        data_1 = [u_in[var_1] for i, u_in in enumerate(data_1)]
        data_2 = [u_out[var_2] for i, u_out in enumerate(data_2)]
        lab_1 = us_ds_1+" "+var_1
        lab_2 = us_ds_2+" "+var_2
        print "Plot var", lab_1, "vs", lab_2, "for", len(data_1), "events"
        hist = xboa.common.make_root_histogram(name,
                                               data_1, lab_1, 50,
                                               data_2, lab_2, 50,
                                               xmin=min_max_1[0], xmax=min_max_1[1],
                                               ymin=min_max_2[0], ymax=min_max_2[1])
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("COL")
        canvas.Update()
        for format in ["png", "root", "eps"]:
            canvas.Print(self.plot_dir+name+"."+format)

    def get_max(self, hist):
        bins = [hist.GetBinContent(i) for i in range(1, 101)]
        max_index = bins.index(max(bins))
        bin_width = hist.GetBinWidth(max_index) # regular bins
        return bins[max_index]/bin_width

    def plot_kalman_pvalue(self, tracker):
        pass

    def plot_sp_residuals(self, tracker, axis):
        residuals_list = []
        residuals_cut_us_list = []
        residuals_cut_ds_list = []
        sp_detector = tracker+"_sp_"+str(self.config.tk_station)
        for event in self.events:
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

        xmin, xmax = -20., 20.
        name = "sp_residuals_"+tracker+"_"+axis
        canvas = xboa.common.make_root_canvas(name)
        for i in range(3):
            hist_raw = xboa.common.make_root_histogram(name, residuals_list, "Space Point "+axis+"-Track point "+axis+" [mm]", 500, xmin=xmin, xmax=xmax)
            hist_cut_us = xboa.common.make_root_histogram(name+"_cut", residuals_cut_us_list, "Space Point "+axis+"-Track point "+axis+" [mm]", 500, xmin=xmin, xmax=xmax)
            fit = scripts.utilities.fit_peak(hist_cut_us)
            mean = fit.GetParameter(1)
            sigma = fit.GetParameter(2)
            xmin, xmax = round(mean-sigma*5, 1), round(mean+sigma*5, 1)

        hist_raw.SetTitle(self.config_analysis['name']+" mean "+str(round(mean, 4))+" sigma "+str(round(sigma, 4)))
        hist_raw.Draw()
        hist_cut_us.SetLineColor(2)
        hist_cut_us.Draw("SAME")
        hist_cut_ds = xboa.common.make_root_histogram(name+"_cut", residuals_cut_ds_list, "Space Point "+axis+"-Track point "+axis+" [mm]", 500, xmin=xmin, xmax=xmax)
        hist_cut_us.SetLineColor(8)
        hist_cut_ds.Draw("SAME")

        canvas.Update()
        for format in ["png", "root", "eps"]:
            canvas.Print(self.plot_dir+name+"."+format)

    def plot_var_1d(self, var_1, us_ds_1, min_max = [None, None]):
        name =  us_ds_1+"_"+var_1
        canvas = xboa.common.make_root_canvas(name)
        data_1 = self.choose_data(us_ds_1, None)
        data_1 = [u_in[var_1] for i, u_in in enumerate(data_1)]
        lab_1 = us_ds_1+" "+var_1
        print "Plot var", lab_1, "for", len(data_1), "events"
        xmin = min_max[0]
        xmax = min_max[1]
        for i in range(3):
            hist = xboa.common.make_root_histogram(name, data_1, lab_1, 100, xmin=xmin, xmax=xmax)
            hist.SetTitle(self.config_analysis['name'])
            fit = ROOT.TF1("testfit", "gaus")
            hist.Fit(fit, "QN")
            print "Mean", fit.GetParameter(1), "rms", fit.GetParameter(2)
            xmin = round(fit.GetParameter(1)-3*fit.GetParameter(2))
            xmax = round(fit.GetParameter(1)+3*fit.GetParameter(2))
            hist.Draw()
            del fit
        if min_max[0] != None:
            xmin = min_max[0]
        if min_max[1] != None:
            xmax = min_max[1]

        hist = xboa.common.make_root_histogram(name, data_1, lab_1, 100, xmin=xmin, xmax=xmax)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw()
        data_1 = self.choose_data(us_ds_1, self.will_cut_us)
        data_1 = [u_in[var_1] for i, u_in in enumerate(data_1)]
        hist = xboa.common.make_root_histogram(name, data_1, lab_1, 100, xmin=xmin, xmax=xmax)
        hist.SetLineColor(2)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("SAME")
        data_1 = self.choose_data(us_ds_1, self.will_cut_ds)
        data_1 = [u_in[var_1] for i, u_in in enumerate(data_1)]
        hist = xboa.common.make_root_histogram(name, data_1, lab_1, 100, xmin=xmin, xmax=xmax)
        hist.SetLineColor(8)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("SAME")
        canvas.Update()
        for format in ["png", "root", "eps"]:
            canvas.Print(self.plot_dir+name+"."+format)
        if var_1 == "p" and us_ds_1 == "tku":
            self.max_p_bin = self.get_max(hist)

    def two_d_projection(self, matrix, i_1, i_2):
        my_cov = numpy.zeros([2, 2])
        for this_i, global_i in enumerate([i_1, i_2]):
          for this_j, global_j in enumerate([i_1, i_2]):
             my_cov[this_i, this_j] = matrix[global_i, global_j]
        return my_cov

    def cut_bunch(self, fout, cut_value, target_bunch):

        cut_weight = None
        cut_bunch = target_bunch.deepcopy()
        while cut_bunch.bunch_weight() != cut_weight:
            cut_weight = cut_bunch.bunch_weight()
            cut_bunch.cut({'amplitude x y':cut_value}, operator.gt)
            
        cut_bunch.clear_local_weights()
        mean_cut = cut_bunch.mean(self.mean_keys)
        beta_cut = cut_bunch.get_beta(['x', 'y'])
        alpha_cut = cut_bunch.get_alpha(['x', 'y'])
        beta_x_cut = cut_bunch.get_beta(['x'])
        alpha_x_cut = cut_bunch.get_alpha(['x'])
        beta_y_cut = cut_bunch.get_beta(['y'])
        alpha_y_cut = cut_bunch.get_alpha(['y'])
        cut_cov = cut_bunch.covariance_matrix(['x', 'px', 'y', 'py'])
        cut_weight = cut_bunch.bunch_weight()

        print >> fout, "Tails cut (cut_weight:", cut_weight, "):"
        print >> fout, "Cut means",
        for key in self.mean_keys:
            print >> fout, key+":", str(round(mean_cut[key], 2)),
        print >> fout, "beta", beta_cut, "alpha", alpha_cut
        print >> fout, "beta x", beta_x_cut, "alpha x", alpha_x_cut
        print >> fout, "beta y", beta_y_cut, "alpha y", alpha_y_cut
        print >> fout, cut_cov

        return cut_cov

    def bunch_plots(self, tracker):
        target_bunch = Bunch.new_from_hits([event[tracker] for event in self.events if event[tracker] != None])
        name = self.config_analysis['name']
        fname = name.replace(' ', '_') 
        fname = self.config_analysis["plot_dir"]+"/"+fname+"_summary_"+tracker
        fout = open(fname+".txt", "w")
        bz = 3e-3
        p = 140.
        if target_bunch.bunch_weight() == 0:
            print >> fout, "No data - check cuts"
            return
        beta = target_bunch.get_beta(['x', 'y'])
        alpha = target_bunch.get_alpha(['x', 'y'])
        eps = target_bunch.get_emittance(['x', 'y'])

        beta_x = target_bunch.get_beta(['x'])
        alpha_x = target_bunch.get_alpha(['x'])
        eps_x = target_bunch.get_emittance(['x'])

        beta_y = target_bunch.get_beta(['y'])
        alpha_y = target_bunch.get_alpha(['y'])
        eps_y = target_bunch.get_emittance(['y'])

        beam_mean = target_bunch.mean(self.mean_keys)
        beam_cov = target_bunch.covariance_matrix(['x', 'px', 'y', 'py'])
        beam_weight = target_bunch.bunch_weight()

        penn_beta = p/bz/0.15*1e-3
        penn_cov = Bunch.build_penn_ellipse(eps, xboa.common.pdg_pid_to_mass[13], penn_beta, 0., p, 0., bz, 1)
        print >> fout, "runs", self.run_numbers
        print >> fout, "name", name
        print >> fout, "Found", beam_weight, "events"
        print >> fout, "means",
        for key in self.mean_keys:
            print >> fout, key+":", str(round(beam_mean[key], 2)),
        print >> fout
        print >> fout, "beta", beta, "alpha", alpha, "emittance:", eps
        print >> fout, "beta x", beta_x, "alpha x", alpha_x, "emittance x:", eps_x
        print >> fout, "beta y", beta_y, "alpha y", alpha_y, "emittance y:", eps_y
        print >> fout, "Matrix (x, px, y, py):"
        print >> fout, beam_cov

        cut_cov = self.cut_bunch(fout, eps*4, target_bunch)

        print >> fout, "target beta", penn_beta, "target matrix:"
        print >> fout, penn_cov
        print >> fout, "theory matrix (based on measured beta, alpha, Lcan = 0):"
        print >> fout, Bunch.build_penn_ellipse(eps, xboa.common.pdg_pid_to_mass[13], beta, alpha, p, 0., bz, 1)
        fout.close()
        fin = open(fname+".txt")
        print fin.read()

        for i_1, i_2, var_1, var_2 in [(0, 1, "x", "px"), (2, 3, "y", "py"), (0, 2, "x", "y")]:
            units = {"x":"mm", "y":"mm", "px":"MeV/c", "py":"MeV/c", "pz":"MeV/c"}
            limits = {"x":150., "y":150., "px":100., "py":100., "pz":300.}
            canvas, hist = target_bunch.root_histogram(var_1, units[var_1],
                                                        var_2, units[var_2],
                                                        100, 100,
                                                        xmin=-limits[var_1], 
                                                        xmax=+limits[var_1], 
                                                        ymin=-limits[var_2], 
                                                        ymax=+limits[var_2])
            canvas.cd()
            hist.SetTitle(self.config_analysis['name'])
            hist.Draw("COLZ")

            ellipse = make_ellipse_graph([beam_mean[var_1], beam_mean[var_2]], self.two_d_projection(beam_cov, i_1, i_2), 1)
            ellipse.Draw("SAMEL")

            #ellipse = make_ellipse_graph([mean_cut[var_1], mean_cut[var_2]], self.two_d_projection(cut_cov, i_1, i_2), 14)
            #ellipse.Draw("SAMEL")

            ellipse = make_ellipse_graph([0., 0.], self.two_d_projection(penn_cov, i_1, i_2), 14)
            ellipse.Draw("SAMEL")
            canvas.Update()

            for format in ["png", "root", "eps"]:
                canvas.Print(self.plot_dir+"bunch_"+tracker+"_"+var_1+"_"+var_2+"."+format)

        canvas, hist = target_bunch.root_histogram("r", "mm", "pt", "MeV/c", nx_bins=50, ny_bins=50, xmin=0., xmax=200., ymin=0., ymax=100.)
        canvas.cd()
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("COLZ")
        canvas.Update()

        canvas.Update()
        for format in ["png", "root", "eps"]:
            canvas.Print(self.plot_dir+"bunch_r_pt."+format)

    # keys used for processing bunch data
    mean_keys = ['x', 'y', 'px', 'py', 'p', 'z', 'energy']

def make_ellipse_graph(means, matrix, line_color):
    shell = xboa.common._common.make_shell(31, matrix)
    shell = [item.tolist()[0] for item in shell]
    shell = sorted(shell, key = lambda x: math.atan2(x[1], x[0]))
    x_list = [item[0]+means[0] for item in shell]
    y_list = [item[1]+means[1] for item in shell]
    hist, graph = xboa.common.make_root_graph("ellipse", x_list, "", y_list, "", sort = False)
    graph.SetLineColor(line_color)
    return graph

def do_plots(config, config_anal, data_loader):
    xboa.common.clear_root()
    plotter = DataPlotter(config, config_anal, data_loader.events, lambda event: event["upstream_cut"], lambda event: event["downstream_cut"], data_loader.run_numbers)
    sys.stdout.flush()

    plotter.plot_tof01()
    plotter.plot_tof12()
    plotter.plot_delta_tof("tof01")
    plotter.plot_delta_tof("tof12")
    plotter.plot_pvalues()
    #plotter.plot_sp_residuals("tku", "x")
    #plotter.plot_sp_residuals("tku", "y")
    #plotter.plot_sp_residuals("tkd", "x")
    #plotter.plot_sp_residuals("tkd", "y")
    plotter.plot_var_2d("x", "tku", "px", "tku", [-150, 150], [-100, 100], plotter.will_cut_us)
    plotter.plot_var_2d("y", "tku", "py", "tku", [-150, 150], [-100, 100], plotter.will_cut_us)
    plotter.plot_var_2d("px", "tku", "py", "tku", [-100, 100], [-100, 100], plotter.will_cut_us)
    plotter.plot_var_2d("x", "tku", "y", "tku", [-150, 150], [-150, 150], plotter.will_cut_us)
    plotter.plot_var_2d("x", "tku", "y", "tku", [-150, 150], [-150, 150], None)
    plotter.plot_var_1d("x", "tku")
    plotter.plot_var_1d("px", "tku")
    plotter.plot_var_1d("y", "tku")
    plotter.plot_var_1d("py", "tku")
    plotter.plot_var_1d("p", "tku", min_max=[80., 180.])
    plotter.plot_var_1d("r", "tku")
    plotter.plot_sp_residuals("tku", "x")
    plotter.plot_p_tot_res()
    plotter.plot_p_tot_vs_tof("tof01", "tku")
    plotter.bunch_plots("tku")

    plotter.bunch_plots("tkd")
    plotter.plot_p_tot_vs_tof("tof12", "tkd")
    plotter.plot_var_2d("y", "tkd", "py", "tkd", [-150, 150], [-100, 100])
    plotter.plot_var_2d("px", "tkd", "py", "tkd", [-100, 100], [-100, 100])
    plotter.plot_var_2d("x", "tkd", "px", "tkd", [-150, 150], [-100, 100])
    plotter.plot_var_1d("x", "tkd")
    plotter.plot_var_1d("px", "tkd")
    plotter.plot_var_1d("y", "tkd")
    plotter.plot_var_1d("py", "tkd")
    plotter.plot_var_1d("p", "tkd", min_max=[80., 180.])
    plotter.plot_var_1d("r", "tkd")
    plotter.plot_var_2d("r", "tkd", "pt", "tkd")

    accepted, accepted_data = plotter.get_cuts_summary()
    wiki_summary = plotter.print_wiki_summary()
    print
    plotter.print_cuts_summary()

    print wiki_summary, "\n\n"
    return wiki_summary
