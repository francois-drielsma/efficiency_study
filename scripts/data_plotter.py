import xboa.common
import json
import copy
import numpy

import cdb_tof_triggers_lookup
import ROOT
from xboa.bunch import Bunch

class DataPlotter(object):
    def __init__(self, config, analysis_index, events, will_cut, run_numbers = None):
        self.events = events
        self.will_cut = will_cut
        self.config = config
        self.config_analysis = config.analyses[analysis_index]
        self.plot_dir = config.analyses[analysis_index]["plot_dir"]
        self.max_tof12_bin = None
        self.max_p_bin = None
        self.run_numbers = run_numbers
        for event in events:
            if self.will_cut(event) and event["tku"] != None:
                event["tku"]["local_weight"] = 0.
            if self.will_cut(event) and event["tkd"] != None:
                event["tkd"]["local_weight"] = 0.
        self.bunch_us = Bunch.new_from_hits([event["tku"] for event in events if event["tku"] != None])
        self.bunch_ds = Bunch.new_from_hits([event["tkd"] for event in events if event["tkd"] != None])

    def choose_data(self, choice, cuts):
        key_map = {
            "upstream":"tku",
            "downstream":"tkd",
        }
        key = key_map[choice]

        if cuts:
            predicate = lambda event: not self.will_cut(event) and event[key] != None
        else:
            predicate = lambda event: event[key] != None
        data = [event[key] for event in self.events if predicate(event)]
        return data

    def print_covariance_matrix(self):
        covariance_matrix_us = self.bunch_us.covariance_matrix(['x', 'px', 'y', 'py'])
        print "number upstream:"
        print self.bunch_us.bunch_weight()
        print "covariances upstream:"
        print covariance_matrix_us
        print
        data_ds = self.choose_data("downstream", cuts = True)
        covariance_matrix_ds = self.bunch_ds.covariance_matrix(['x', 'px', 'y', 'py'])
        print "number downstream:"
        print self.bunch_ds.bunch_weight()
        print "covariances downstream:"
        print covariance_matrix_ds
        print
        return covariance_matrix_us, covariance_matrix_ds

    def get_cuts_summary(self):
        # there are a set of cuts that are to do with basic "data quality" e.g.
        # exactly one hit in TOF, etc and there are a set of cuts that are to do
        # with "bias correction"/etc.
        # I report the number of events that pass each data quality cut;
        # and the number of events that pass each "bias correction" cut AND the data quality cuts
        accepted_dict = {"any_cut":0, "data_cut":None, "all_events":0}
        accepted_data_dict = {"any_cut":None, "data_cut":0, "all_events":0}
        for key in self.events[0]["will_cut"]:
            accepted_dict[key] = 0
            accepted_data_dict[key] = 0
        for event in self.events:
            accepted_dict["all_events"]  += 1 # total number of events considered
            accepted_data_dict["all_events"]  += 1 # total number of events considered
            will_cut = event["will_cut"]
            for key in will_cut:
                if not will_cut[key]:
                    accepted_dict[key] += 1
                    if event["data_cut"]:
                        continue
                    accepted_data_dict[key] += 1
            if not event["any_cut"]:
                accepted_dict["any_cut"]  += 1 # total number accepted by all cuts
            if not event["data_cut"]:
                accepted_data_dict["data_cut"]  += 1 # total number accepted by data cuts
        return accepted_dict, accepted_data_dict

    def print_cuts_summary(self):
        cut_dict, cut_data_dict = self.get_cuts_summary()
        for key in sorted(cut_dict.keys()):
            if key in self.config.cuts_active:
                is_active = self.config.cuts_active[key]
            else:
                is_active = None
            print "   ", key.ljust(25), str(is_active).ljust(5), str(cut_dict[key]).ljust(8), str(cut_data_dict[key]).ljust(8)

    def print_wiki_summary(self):
        cut_dict, cut_data_dict = self.get_cuts_summary()
        
        cut_dict['data_cut'] = cut_data_dict['data_cut']
        cut_dict['runs'] = sorted([run for run in self.run_numbers])
        cut_dict['max_p_bin'] = round(self.max_p_bin, 2)
        cut_dict['max_tof01_bin'] = round(self.max_tof01_bin, 2)
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


    def plot_tof01(self):
        tof01_no_cut = [event["tof01"] for event in self.events if not self.will_cut(event)]
        tof01_all = [event["tof01"] for event in self.events]
        tof01_all = [tof01 for tof01 in tof01_all if tof01 != None]
        tof01_canvas = xboa.common.make_root_canvas("tof01_canvas")
        tof01_canvas.SetLogy()
        xmin = min(25., self.config_analysis["tof01_cut_low"])
        xmax = max(45., self.config_analysis["tof01_cut_high"])
        hist = xboa.common.make_root_histogram("tof01", tof01_all, "tof1 - tof0 [ns]", 100, xmin = xmin, xmax = xmax)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw()
        hist = xboa.common.make_root_histogram("tof01 in cut", tof01_no_cut, "tof1 - tof0 [ns]", 100, xmin = xmin, xmax = xmax)
        hist.SetLineColor(4)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("SAME")
        tof01_canvas.Update()
        self.max_tof01_bin = self.get_max(hist)
        for format in ["png", "root", "pdf"]:
            tof01_canvas.Print(self.plot_dir+"tof01."+format)


    def plot_tof12(self):
        tof12_no_cut = [event["tof12"] for event in self.events if not self.will_cut(event) and event["tof12"] != None]
        tof12_all = [event["tof12"] for event in self.events if event["tof12"] != None]
        tof12_canvas = xboa.common.make_root_canvas("tof12_canvas")
        tof12_canvas.SetLogy()
        xmin = min(25., self.config_analysis["tof12_cut_low"])
        xmax = max(45., self.config_analysis["tof12_cut_high"])
        hist = xboa.common.make_root_histogram("tof12", tof12_all, "tof2 - tof1 [ns]", 100, xmin = xmin, xmax = xmax)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw()
        hist = xboa.common.make_root_histogram("tof12 in cut", tof12_no_cut, "tof2 - tof1 [ns]", 100, xmin = xmin, xmax = xmax)
        hist.SetLineColor(4)
        hist.Draw("SAME")
        tof12_canvas.Update()
        self.max_tof12_bin = max([hist.GetBinContent(i) for i in range(102)])
        for format in ["png", "root", "pdf"]:
            tof12_canvas.Print(self.plot_dir+"tof12."+format)

    def plot_p_tot_vs_tof(self):
        canvas = xboa.common.make_root_canvas("p_tot vs tof no cuts")
        tof_no_cut = [event["tof01"] for event in self.events if event["tof01"] != None and event["tku"] != None ]
        p_tot_no_cut = [event["tku"]["p"] for event in self.events if event["tof01"] != None and event["tku"] != None ]

        hist = xboa.common.make_root_histogram("Events in cut", tof_no_cut, "tof1 - tof0 [ns]", 100, p_tot_no_cut, "p_{tku} [MeV/c]", 100, ymin=0., ymax=300., xmin=25., xmax=45.)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("COLZ")
        canvas.Update()
        for format in ["png", "root", "pdf"]:
            canvas.Print(self.plot_dir+"p_tot_vs_tof01_all."+format)

        canvas = xboa.common.make_root_canvas("p_tot vs tof cuts")
        tof_no_cut = [event["tof01"] for event in self.events if not self.will_cut(event)]
        p_tot_no_cut = [event["tku"]["p"] for event in self.events if not self.will_cut(event)]

        hist = xboa.common.make_root_histogram("Events in cut", tof_no_cut, "tof1 - tof0 [ns]", 100, p_tot_no_cut, "p_{tku} [MeV/c]", 100)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("COLZ")
        canvas.Update()
        for format in ["png", "root", "pdf"]:
            canvas.Print(self.plot_dir+"p_tot_vs_tof01_cuts."+format)

    def plot_p_tot_res(self):
        canvas = xboa.common.make_root_canvas("p_tot us vs ds")
        cut_pred = lambda event: event["tku"] != None and event["tkd"] != None
        p_tot_us_no_cut = [event["tku"]["p"] for event in self.events if cut_pred(event)]
        p_tot_ds_no_cut = [event["tkd"]["p"] for event in self.events if cut_pred(event)]
        p_tot_res_no_cut = [p_tot - p_tot_ds_no_cut[i] for i, p_tot in enumerate(p_tot_us_no_cut)]

        cut_pred = lambda event: event["tku"] != None and event["tkd"] != None and not self.will_cut(event)
        p_tot_us_cut = [event["tku"]["p"] for event in self.events if cut_pred(event)]
        p_tot_ds_cut = [event["tkd"]["p"] for event in self.events if cut_pred(event)]
        p_tot_res_cut = [p_tot - p_tot_ds_cut[i] for i, p_tot in enumerate(p_tot_us_cut)]
        mean = numpy.mean(p_tot_res_cut)
        std = numpy.std(p_tot_res_cut)

        hist = xboa.common.make_root_histogram("All events", p_tot_res_no_cut, "p_{tku} - p_{tkd} [MeV/c]", 100, xmin=-50, xmax=50)
        hist.SetTitle(self.config_analysis['name'])
        hist.SetTitle("mean: "+str(mean)+" std: "+str(std))
        hist.Draw()
        hist = xboa.common.make_root_histogram("Events in cut", p_tot_res_cut, "p_{tku} - p_{tkd} [MeV/c]", 100, xmin=-50, xmax=50)
        hist.SetTitle("mean: "+str(mean)+" std: "+str(std))
        hist.SetLineColor(4)
        hist.Draw("SAME")
        canvas.Update()
        for format in ["png", "root", "pdf"]:
            canvas.Print(self.plot_dir+"p_tot_res."+format)

        canvas = xboa.common.make_root_canvas("delta p vs p_tku")
        hist = xboa.common.make_root_histogram("Events in cut", p_tot_us_cut, "p_{tku} [MeV/c]", 50,
                                                                p_tot_res_cut, "p_{tku} - p_{tkd} [MeV/c]", 50,
                                                                ymin=-50, ymax=50)
        hist.Draw("COLZ")
        canvas.Update()
        for format in ["png", "root", "pdf"]:
            canvas.Print(self.plot_dir+"delta_p_vs_p_tku."+format)


    def plot_var_2d(self, var_1, us_ds_1, var_2, us_ds_2, do_cuts = True): # MAKE THIS A CONT PLOT AND SUPERIMPOSE EACH TOF BIN
        name =  us_ds_1+"_"+var_1+"_"+ us_ds_2+"_"+var_2+"_cuts-"+str(do_cuts)
        canvas = xboa.common.make_root_canvas(name)
        data_1 = self.choose_data(us_ds_1, do_cuts)
        data_2 = self.choose_data(us_ds_2, do_cuts)
        data_1 = [u_in[var_1] for i, u_in in enumerate(data_1)]
        data_2 = [u_out[var_2] for i, u_out in enumerate(data_2)]
        lab_1 = us_ds_1+" "+var_1
        lab_2 = us_ds_2+" "+var_2
        print "Plot var", lab_1, "vs", lab_2, "for", len(data_1), "events"
        hist = xboa.common.make_root_histogram(name,
                                               data_1, lab_1, 50,
                                               data_2, lab_2, 50)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("COL")
        canvas.Update()
        for format in ["png", "root", "pdf"]:
            canvas.Print(self.plot_dir+name+"."+format)

    def get_max(self, hist):
        bins = [hist.GetBinContent(i) for i in range(1, 101)]
        max_index = bins.index(max(bins))
        bin_width = hist.GetBinWidth(max_index) # regular bins
        return bins[max_index]/bin_width

    def plot_var_1d(self, var_1, us_ds_1):
        name =  us_ds_1+"_"+var_1
        canvas = xboa.common.make_root_canvas(name)
        data_1 = self.choose_data(us_ds_1, False)
        data_1 = [u_in[var_1] for i, u_in in enumerate(data_1)]
        lab_1 = us_ds_1+" "+var_1
        print "Plot var", lab_1, "for", len(data_1), "events"
        xmin = None
        xmax = None
        for i in range(3):
            hist = xboa.common.make_root_histogram(name, data_1, lab_1, 100, xmin=xmin, xmax=xmax)
            hist.SetTitle(self.config_analysis['name'])
            fit = ROOT.TF1("testfit", "gaus")
            hist.Fit(fit, "QN")
            #print fit.GetParameter(0), fit.GetParameter(1), fit.GetParameter(2)
            xmin = round(fit.GetParameter(1)-5*fit.GetParameter(2))
            xmax = round(fit.GetParameter(1)+5*fit.GetParameter(2))
            hist.Draw()
            del fit

        hist = xboa.common.make_root_histogram(name, data_1, lab_1, 100, xmin=xmin, xmax=xmax)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw()
        data_1 = self.choose_data(us_ds_1, True)
        data_1 = [u_in[var_1] for i, u_in in enumerate(data_1)]
        hist = xboa.common.make_root_histogram(name, data_1, lab_1, 100, xmin=xmin, xmax=xmax)
        hist.SetLineColor(4)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("SAME")
        canvas.Update()
        for format in ["png", "root", "pdf"]:
            canvas.Print(self.plot_dir+name+"."+format)
        if var_1 == "p" and us_ds_1 == "upstream":
            self.max_p_bin = self.get_max(hist)

    def bunch_plots(self):
        name = self.config_analysis['name']
        fout = open(self.config.data_dir+"/"+name+".txt", "w")
        json_out = open(self.config.data_dir+"/"+name+".json", "w")
        bz = 3e-3
        p = 140.
        beta = self.bunch_us.get_beta(['x', 'y'])
        alpha = self.bunch_us.get_alpha(['x', 'y'])
        eps = self.bunch_us.get_emittance(['x', 'y'])
        beam_mean = self.bunch_us.mean(['x', 'px', 'y', 'py', 'p', 'z', 'energy'])
        beam_cov = self.bunch_us.covariance_matrix(['x', 'px', 'y', 'py'])
        penn_beta = p/bz/0.15*1e-3
        penn_cov = Bunch.build_penn_ellipse(eps, xboa.common.pdg_pid_to_mass[13], penn_beta, 0., p, 0., bz, 1)
        data = {
            "name":name,
            "weight":self.bunch_us.bunch_weight(),
            "means":beam_mean,
            "measured_cov":beam_cov.tolist(),
            "penn_cov":penn_cov.tolist(),
        }
        print >> json_out, json.dumps(data)
        print >> fout, "name", name
        print >> fout, "Found", self.bunch_us.bunch_weight(), "events"
        print >> fout, "means", beam_mean
        print >> fout, "beta", beta, "alpha", alpha, "emittance:", eps, "matrix:"
        print >> fout, beam_cov
        print >> fout, "target beta", penn_beta, "target matrix:"
        print >> fout, penn_cov
        print >> fout, "theory matrix (based on measured beta, alpha, Lcan = 0):"
        print >> fout, Bunch.build_penn_ellipse(eps, xboa.common.pdg_pid_to_mass[13], beta, alpha, p, 0., bz, 1)
        fout.close()
        fin = open(self.config.data_dir+"/"+name+".txt")
        print fin.read()

        for i_low, i_high, var_1, var_2 in [(0, 2, "x", "px"), (2, 4, "y", "py")]:
            canvas, hist = self.bunch_us.root_histogram(var_1, "mm", var_2, "MeV/c", xmin=-150., xmax=150., ymin=-100., ymax=100.)
            canvas.cd()
            hist.SetTitle(self.config_analysis['name'])
            hist.Draw("COLZ")
            canvas.Update()
            ellipse = xboa.common.make_root_ellipse_function([beam_mean[var_1], beam_mean[var_2]], beam_cov[i_low:i_high, i_low:i_high],
                                                            contours=[4])
            ellipse.SetNpx(1000)
            ellipse.SetNpy(1000)
            ellipse.SetLineColor(1)
            ellipse.Draw("SAME")
            canvas.Update()

            ellipse = xboa.common.make_root_ellipse_function([0., 0.], penn_cov[i_low:i_high, i_low:i_high], contours=[4.])
            ellipse.SetNpx(1000)
            ellipse.SetNpy(1000)
            ellipse.SetLineColor(8)
            ellipse.Draw("SAME")
            canvas.Update()

            canvas.Update()
            for format in ["png", "root", "pdf"]:
                canvas.Print(self.plot_dir+"bunch_"+var_1+"_"+var_2+"."+format)

        canvas, hist = self.bunch_us.root_histogram("r", "mm", "pt", "MeV/c", nx_bins=50, ny_bins=50, xmin=0., xmax=200., ymin=0., ymax=100.)
        canvas.cd()
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("COLZ")
        canvas.Update()

        canvas.Update()
        for format in ["png", "root", "pdf"]:
            canvas.Print(self.plot_dir+"bunch_r_pt."+format)


