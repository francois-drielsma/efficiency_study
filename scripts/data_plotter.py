import xboa.common
import json
import copy
import numpy

import cdb_tof_triggers_lookup
import ROOT
from xboa.bunch import Bunch
import scripts.utilities

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
        print "========== cuts summary ============"
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

    def plot_pvalues(self):
        for tracker in "tku", "tkd":
            pvalue_canvas = xboa.common.make_root_canvas(tracker+" PValues")
            pvalues_all = []
            pvalues_cut = []
            for event in self.events:
                if event[tracker] == None:
                    continue
                for detector_hit in event["data"]:
                    if detector_hit["detector"] != tracker+"_tp":
                        continue
                    pvalues_all.append(detector_hit["pvalue"])
                    if not self.will_cut(event):
                        pvalues_cut.append(detector_hit["pvalue"])
            hist = xboa.common.make_root_histogram("pvalues "+tracker, pvalues_all, "P Value ("+tracker+")", 120, xmin = -0.1, xmax = 1.1)
            hist.SetTitle(self.config_analysis['name'])
            hist.Draw()
            hist = xboa.common.make_root_histogram("pvalues "+tracker, pvalues_cut, "P Value ("+tracker+")", 120, xmin = -0.1, xmax = 1.1)
            hist.SetLineColor(4)
            hist.SetTitle(self.config_analysis['name'])
            hist.Draw("SAME")
            pvalue_canvas.SetLogy()
            pvalue_canvas.Update()
            for format in ["png", "root", "pdf"]:
                pvalue_canvas.Print(self.plot_dir+"pvalue_"+tracker+"."+format)


    def plot_delta_tof(self, tofs):
        dtof01_no_cut = [event["delta_"+tofs] for event in self.events if not self.will_cut(event)]
        dtof01_no_cut = [dtof01 for dtof01 in dtof01_no_cut if dtof01 != None]
        dtof01_all = [event["delta_"+tofs] for event in self.events]
        dtof01_all = [dtof01 for dtof01 in dtof01_all if dtof01 != None]
        if len(dtof01_all) == 0:
            print "No delta tof plot - perhaps extrapolation is disabled?"
            return
        dtof01_canvas = xboa.common.make_root_canvas("delta_tof01_canvas")
        dtof01_canvas.SetLogy()
        xmin = -10.
        xmax = 10.
        hist = xboa.common.make_root_histogram("delta "+tofs, dtof01_all, "Reco "+tofs+" - Extrapolated "+tofs+" [ns]", 100, xmin = xmin, xmax = xmax)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw()
        hist = xboa.common.make_root_histogram("delta "+tofs+" in cut", dtof01_no_cut, "Reco "+tofs+" - Extrapolated "+tofs+" [ns]", 100, xmin = xmin, xmax = xmax)
        hist.SetLineColor(4)
        hist.SetTitle(self.config_analysis['name'])
        hist.Draw("SAME")
        dtof01_canvas.Update()
        for format in ["png", "root", "pdf"]:
            dtof01_canvas.Print(self.plot_dir+"delta_"+tofs+"."+format)

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

    def plot_kalman_pvalue(self, tracker):
        pass

    def plot_sp_residuals(self, tracker, axis):
        residuals_list = []
        residuals_cut_list = []
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
            is_cut = event["any_cut"]
            if sp == None:
                continue # no space point recorded (e.g. due to dead channel)
            residuals_list.append(sp[axis] - tp[axis])        
            if not is_cut:
                residuals_cut_list.append(sp[axis] - tp[axis])        

        xmin, xmax = -20., 20.
        name = "sp_residuals_"+tracker+"_"+axis
        canvas = xboa.common.make_root_canvas(name)
        for i in range(3):
            hist_raw = xboa.common.make_root_histogram(name, residuals_list, "Space Point "+axis+"-Track point "+axis+" [mm]", 500, xmin=xmin, xmax=xmax)
            hist_cut = xboa.common.make_root_histogram(name+"_cut", residuals_cut_list, "Space Point "+axis+"-Track point "+axis+" [mm]", 500, xmin=xmin, xmax=xmax)
            fit = scripts.utilities.fit_peak(hist_cut)
            mean = fit.GetParameter(1)
            sigma = fit.GetParameter(2)
            xmin, xmax = round(mean-sigma*5, 1), round(mean+sigma*5, 1)

        hist_raw.SetTitle(self.config_analysis['name']+" mean "+str(round(mean, 4))+" sigma "+str(round(sigma, 4)))
        hist_raw.Draw()
        hist_cut.SetLineColor(4)
        hist_cut.Draw("SAME")

        canvas.Update()
        for format in ["png", "root", "pdf"]:
            canvas.Print(self.plot_dir+name+"."+format)

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
            print "Mean", fit.GetParameter(1), "rms", fit.GetParameter(2)
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

    def two_d_projection(self, matrix, i_1, i_2):
        my_cov = numpy.zeros([2, 2])
        for this_i, global_i in enumerate([i_1, i_2]):
          for this_j, global_j in enumerate([i_1, i_2]):
             my_cov[this_i, this_j] = matrix[global_i, global_j]
        return my_cov

    def bunch_plots(self):
        name = self.config_analysis['name']
        fname = name.replace(' ', '_') 
        fname = self.config_analysis["plot_dir"]+"/"+name+"_summary"
        fout = open(fname+".txt", "w")
        json_out = open(fname+".json", "w")
        bz = 3e-3
        p = 240.
        beta = self.bunch_us.get_beta(['x', 'y'])
        alpha = self.bunch_us.get_alpha(['x', 'y'])
        eps = self.bunch_us.get_emittance(['x', 'y'])

        beta_x = self.bunch_us.get_beta(['x'])
        alpha_x = self.bunch_us.get_alpha(['x'])
        eps_x = self.bunch_us.get_emittance(['x'])

        beta_y = self.bunch_us.get_beta(['y'])
        alpha_y = self.bunch_us.get_alpha(['y'])
        eps_y = self.bunch_us.get_emittance(['y'])

        mean_keys = ['x', 'y', 'px', 'py', 'p', 'z', 'energy']
        beam_mean = self.bunch_us.mean(mean_keys)
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
        print >> fout, "runs", self.run_numbers
        print >> fout, "name", name
        print >> fout, "Found", self.bunch_us.bunch_weight(), "events"
        print >> fout, "means",
        for key in mean_keys:
            print >> fout, key+":", str(round(beam_mean[key], 2)),
        print >> fout
        print >> fout, "beta", beta, "alpha", alpha, "emittance:", eps
        print >> fout, "beta x", beta_x, "alpha x", alpha_x, "emittance x:", eps_x
        print >> fout, "beta y", beta_y, "alpha y", alpha_y, "emittance y:", eps_y
        print >> fout, "Matrix (x, px, y, py):"
        print >> fout, beam_cov
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
            
            canvas, hist = self.bunch_us.root_histogram(var_1, units[var_1],
                                                        var_2, units[var_2],
                                                        xmin=-limits[var_1], 
                                                        xmax=+limits[var_1], 
                                                        ymin=-limits[var_2], 
                                                        ymax=+limits[var_2])
            canvas.cd()
            hist.SetTitle(self.config_analysis['name'])
            hist.Draw("COLZ")
            canvas.Update()
            # not tested
            ellipse = xboa.common.make_root_ellipse_function([beam_mean[var_1], beam_mean[var_2]], self.two_d_projection(beam_cov, i_1, i_2),
                                                            contours=[2.])
            ellipse.SetNpx(1000)
            ellipse.SetNpy(1000)
            ellipse.SetLineColor(1)
            ellipse.Draw("SAME")
            canvas.Update()

            ellipse = xboa.common.make_root_ellipse_function([0., 0.], self.two_d_projection(penn_cov, i_1, i_2), contours=[4.])
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


