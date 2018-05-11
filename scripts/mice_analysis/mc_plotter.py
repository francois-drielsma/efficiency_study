import xboa.common
import json
import copy
import numpy
import sys

import utilities.cdb_tof_triggers_lookup
import ROOT
from xboa.bunch import Bunch
import utilities.utilities

from analysis_base import AnalysisBase

class MCPlotter(AnalysisBase):
    def __init__(self, config, config_anal, data_loader):
        super(MCPlotter, self).__init__(config, config_anal, data_loader)
        self.process_args = {}
        self.mc_stations = {}
        self.failed_pids = {}

    def birth(self):
        self.set_plot_dir("mc_plots")
        self.mc_stations = self.config.mc_plots["mc_stations"]
        pid_colors = {211:ROOT.kGreen, -13:ROOT.kRed, -11:ROOT.kBlue-1, +11:ROOT.kBlue+1}
        for detector, virt_station_list in self.mc_stations.iteritems():
            virt_station = virt_station_list[0]
            self.birth_var_one_d("p_at_"+virt_station, virt_station, "pid", "p", pid_colors, "us cut", 0., 300.)
            self.birth_var_one_d("r_at_"+virt_station, virt_station, "pid", "r", pid_colors, "us cut", 0., 300.)
            self.birth_var_two_d_scatter("x_vs_px_at_"+virt_station, virt_station, "pid", "x", "px", pid_colors, True)
        self.birth_var_two_d_scatter("x_vs_z_of_mc_track_final", "mc_track_final", "pid", "z", "x", pid_colors)
        self.birth_var_two_d_scatter("y_vs_z_of_mc_track_final", "mc_track_final", "pid", "z", "y", pid_colors)
        self.birth_var_two_d_scatter("r_vs_z_of_mc_track_final", "mc_track_final", "pid", "z", "r", pid_colors)
        self.birth_var_two_d_scatter("x_vs_px_of_primary", "mc_primary", "pid", "x", "px", pid_colors)
        self.birth_var_two_d_scatter("y_vs_py_of_primary", "mc_primary", "pid", "y", "py", pid_colors)
        self.birth_var_one_d("p_at_mc_tof_0", "mc_tof_0", "pid", "p", pid_colors)
        self.birth_var_one_d("p_at_mc_tof_1", "mc_tof_1", "pid", "p", pid_colors)
        self.birth_var_one_d("e_dep_at_mc_tof_0", "mc_tof_0", "pid", "e_dep", pid_colors)
        self.birth_var_one_d("e_dep_at_mc_tof_0", "mc_tof_1", "pid", "e_dep", pid_colors)
        self.birth_var_one_d("p_at_mc_primary", "mc_primary", "pid", "p", pid_colors)
        self.birth_var_one_d("x_at_mc_primary", "mc_primary", "pid", "x", pid_colors)
        self.birth_var_one_d("y_at_mc_primary", "mc_primary", "pid", "y", pid_colors)
        self.birth_var_one_d("p_at_mc_track_final", "mc_track_final", "pid", "p", pid_colors)
        self.birth_var_one_d("x_at_mc_track_final", "mc_track_final", "pid", "x", pid_colors)
        self.birth_var_one_d("y_at_mc_track_final", "mc_track_final", "pid", "y", pid_colors)

        self.birth_var_one_d("tku_ke_in_st_1_plane_0", "mc_tk_111", "pid", "kinetic_energy", pid_colors, xmin = 0., xmax = 10.)
        self.birth_var_one_d("tku_e_dep_in_st_1_plane_0", "mc_tk_111", "pid", "e_dep", pid_colors, xmin = 0., xmax = 0.5)
        my_options = {
          "sub_dir":"accumulated_e_dep",
          "logy":True,
        }

        for tracker in range(1, 3):
            for station in range(1, 6):
                for plane in range(1, 4):
                    det_str = "mc_tk_"+str(100*tracker + 10*station + plane)
                    self.birth_var_one_d("tku_accumulated_e_dep_in_"+det_str, det_str, "pid", "accumulated_e_dep", pid_colors, cuts = "us cut", xmin = 0., xmax = 1.0)

        tku_hit_predicate = lambda hit: 13900 < hit["hit"]["z"] and hit["hit"]["z"] < 15100
        self.birth_var_one_d("tku_p_at_mc_track_initial", "mc_track_initial", "pid", "p", pid_colors, xmin = 0., xmax = 10., hit_predicate = tku_hit_predicate)
        self.birth_var_one_d("tku_z_at_mc_track_initial", "mc_track_initial", "pid", "z", pid_colors, hit_predicate = tku_hit_predicate)

        tkd_hit_predicate = lambda hit: 18800 < hit["hit"]["z"] and hit["hit"]["z"] < 20000
        self.birth_var_one_d("tkd_p_at_mc_track_initial", "mc_track_initial", "pid", "p", pid_colors, xmin = 0., xmax = 10., hit_predicate = tkd_hit_predicate)
        self.birth_var_one_d("tkd_z_at_mc_track_initial", "mc_track_initial", "pid", "z", pid_colors, hit_predicate = tkd_hit_predicate)
        self.birth_var_two_d_scatter("z_vs_p_of_track_initial", "mc_track_initial", "pid", "z", "p", pid_colors)
        self.birth_data_detector_residuals()

    def process(self):
        for name in sorted(self.process_args.keys()):
            process_function = self.process_args[name][0]
            process_args = self.process_args[name][1]
            process_function(name, *process_args)
        self.process_data_detector_residuals()

    def death(self):
        self.death_data_detector_residuals()
        print "Failed to put the following pids into scatter plots {pid:number}:"
        print "   ", self.failed_pids
        self.print_plots()

    def get_accumulated_e_dep(self, event, detector, hit_predicate, track_final):
        """
        Return a list of accumulated energy deposited
        
        Makes one entry per detector per pid per event
        """
        e_dep_dict = {} # pertains to this event/detector
        for detector_hit in event["data"]:
            if detector_hit["detector"] != detector:
                continue
            if hit_predicate != None and not hit_predicate(detector_hit):
                continue
            hit = detector_hit["hit"]
            pid = hit["pid"]
            if pid not in e_dep_dict:
                e_dep_dict[pid] = hit["e_dep"]
            else:
                e_dep_dict[pid] += hit["e_dep"]
        # now accumulate over the e_dep_dict
        for pid, e_dep in e_dep_dict.iteritems():
            if pid not in track_final:
                track_final[pid] = [e_dep]
            else:
                track_final[pid].append(e_dep)

 
    def get_data_var_one_d(self, detector, slice_variable, plot_variable, cuts, hit_predicate):
        track_final = {}
        cut_lambda = {
            "all":lambda event: False,
            "us cut":lambda event: event["upstream_cut"],
            "ds cut":lambda event: event["downstream_cut"],
        }[cuts]
        for event in self.data_loader.events:
            if cut_lambda(event):
                continue
            for detector_hit in event["data"]:
                if detector_hit["detector"] != detector:
                    continue
                if hit_predicate != None and not hit_predicate(detector_hit):
                    continue
                if plot_variable == "accumulated_e_dep":
                    if slice_variable != "pid":
                        raise RuntimeError("NO NO NO NO BADGER ATTACK")
                    self.get_accumulated_e_dep(event, detector, hit_predicate, track_final)
                    break
                else:
                    slice_var = detector_hit["hit"][slice_variable]
                    plot_var = detector_hit["hit"][plot_variable]
                    if slice_var not in track_final:
                        track_final[slice_var] = []
                    track_final[slice_var].append(plot_var)
                    break
        all_list = []
        for value in track_final.values():
            all_list += value
        return (all_list, track_final)

    def process_var_one_d(self, name, detector, slice_variable, plot_variable, color_dict, cuts, hit_predicate):
        all_list, track_final = self.get_data_var_one_d(detector, slice_variable, plot_variable, cuts, hit_predicate)
        for item in all_list:
            self.plots[name]["histograms"]["all"].Fill(item)
        for i, key in enumerate(sorted(track_final.keys())):
            hist_dict = self.plots[name]["histograms"]
            plot_name = str(slice_variable)+" = "+str(key)
            if plot_name not in hist_dict:
                if key not in self.failed_pids:
                    self.failed_pids[key] = 0
                self.failed_pids[key] += 1
                continue
            for item in track_final[key]:
                hist_dict[plot_name].Fill(item)

    def birth_var_one_d(self, canvas_name, detector, slice_variable, plot_variable, color_dict, cuts = "all", xmin = None, xmax = None, hit_predicate = None, options = {}):
        all_list, track_final = self.get_data_var_one_d(detector, slice_variable, plot_variable, cuts, hit_predicate)
        hist = self.make_root_histogram(canvas_name, "all", all_list, plot_variable+" at "+detector, 100, [], '', 0, [], xmin, xmax)
        hist.SetMarkerStyle(26)
        hist.Draw("P")
        hist.SetStats(True)
        for i, key in enumerate(sorted(track_final.keys())):
            var_list = track_final[key]
            name = slice_variable+" = "+str(key)
            label = plot_variable+" at "+detector+" ["+utilities.utilities.default_units(plot_variable)+"]"
            hist = self.make_root_histogram(canvas_name, name, var_list, label, 100, [], '', 0, [], xmin, xmax)
            hist.SetMarkerStyle(24)
            if key in color_dict:
                hist.SetMarkerColor(color_dict[key])
            hist.Draw("SAMEP")
        self.plots[canvas_name]["canvas"].SetLogy()
        #self.plots[canvas_name]["canvas"].BuildLegend()
        self.process_args[canvas_name] = [self.process_var_one_d, (detector, slice_variable, plot_variable, color_dict, cuts, hit_predicate)]
        for key in options.keys():
            if key not in self.get_plot(canvas_name)["config"]:
                raise KeyError("Did not recignise plot option "+str(key))
            self.get_plot(canvas_name)["config"][key] = options[key]


    def get_data_var_two_d_scatter(self, name, *args):
        track_final = {}
        for event in self.data_loader.events:
            if args[5] and event["downstream_cut"]:
                continue
            for detector_hit in event["data"]:
                if detector_hit["detector"] != args[0]:
                    continue
                pid = detector_hit["hit"][args[1]]
                x = detector_hit["hit"][args[2]]
                y = detector_hit["hit"][args[3]]
                if pid not in track_final:
                    track_final[pid] = ([], [])
                track_final[pid][0].append(x)
                track_final[pid][1].append(y)
        return track_final

    def process_var_two_d_scatter(self, name, *args):
        track_final = self.get_data_var_two_d_scatter(name, *args)
        for pid in sorted(track_final.keys()):
            graph_dict = self.plots[name]["graphs"]
            plot_name = str(args[1])+" = "+str(pid)
            if plot_name not in graph_dict:
                if pid not in self.failed_pids:
                    self.failed_pids[pid] = 0
                self.failed_pids[pid] += 1
                continue
            graph = graph_dict[plot_name]
            n_points = len(track_final[pid][0])
            n_orig = graph.GetN()
            graph.Set(n_orig+n_points)
            for i in range(n_points):
                graph_dict[plot_name].SetPoint(i + n_orig, track_final[pid][0][i], track_final[pid][1][i])

    def birth_var_two_d_scatter(self, canvas_name, detector, slice_variable, plot_variable_1, plot_variable_2, color_dict, cuts=False):
        track_final = self.get_data_var_two_d_scatter(canvas_name, detector, slice_variable, plot_variable_1, plot_variable_2, color_dict, cuts)
        if len(track_final) == 0:
            print "No tracks for", detector
            return
        xmin, xmax, ymin, ymax = 0., 0., 0., 0.
        for x_list, y_list in track_final.values():
            xmin = min([xmin]+x_list)
            xmax = max([xmax]+x_list)
            ymin = min([ymin]+y_list)
            ymax = max([ymax]+y_list)

        for i, pid in enumerate(sorted(track_final.keys())):
            x_list = track_final[pid][0]
            y_list = track_final[pid][1]
            name = slice_variable+" = "+str(pid)
            label_1 = plot_variable_1+" ["+utilities.utilities.default_units(plot_variable_1)+"]"
            label_2 = plot_variable_2+" ["+utilities.utilities.default_units(plot_variable_2)+"]"
            hist, graph = self.make_root_graph(canvas_name, name,
              x_list, label_1, y_list, label_2, True,
              xmin, xmax, ymin, ymax)
            if i == 0:
                hist.Draw()
            if pid in color_dict:
                graph.SetMarkerColor(color_dict[pid])
            graph.SetMarkerStyle(6)
            graph.Draw("PSAME")
        self.process_args[canvas_name] = [self.process_var_two_d_scatter, (detector, slice_variable, plot_variable_1, plot_variable_2, color_dict, cuts)]

    axis_labels = {"x":"x(meas) - x(true) [mm]", "y":"y(meas) - y(true) [mm]", "z":"z(meas) - z(true) [mm]",
                   "px":"p_{x}(meas) - p_{x}(true) [MeV/c]", "py":"p_{y}(meas) - p_{y}(true) [MeV/c]", "pz":"p_{z}(meas) - p_{z}(true) [MeV/c]"}
    def plot_detector_residuals(self, suffix):
        self.axis_min_max = {}
        for var in self.residual_dict.keys():
            if var == "amp_virt":
                continue  
            print "Plot detector residuals doing var", var
            data = self.residual_dict[var]
            canvas = xboa.common.make_root_canvas("MC residual "+var)
            canvas.Draw()
            xmin, xmax = utilities.utilities.fractional_axis_range(data, 0.95)
            hist = xboa.common.make_root_histogram(var, data, self.axis_labels[var], 100, xmin=xmin, xmax=xmax)
            hist.Draw()
            fit = utilities.utilities.fit_peak(hist, nsigma=8)
            mean = fit.GetParameter(1)
            sigma = fit.GetParameter(2)
            xmin, xmax = mean-5*sigma, mean+5*sigma
            self.axis_min_max[var] = (xmin, xmax)
            hist = xboa.common.make_root_histogram(var, data, self.axis_labels[var], 100, xmin=xmin, xmax=xmax)
            hist.Draw()
            fit = utilities.utilities.fit_peak(hist, nsigma=1)

            text_box = utilities.utilities.get_text_box(self.config, self.config_anal, data, fit)
            canvas.Update()
            for format in ["png", "eps", "root"]:
                canvas.Print(self.plot_dir+"/mc_residual_"+suffix+"_"+var+"."+format)

    def get_data_detector_residuals(self, tracker, virt_name):
        residual_dict = {"x":[], "y":[], "z":[], "px":[], "py":[], "pz":[]}
        for event in self.data_loader.events:
            if event["downstream_cut"]:
                continue
            if event[tracker] == None:
                continue
            for detector_hit in event["data"]:
                if detector_hit["detector"] != virt_name:
                    continue
                vhit = detector_hit["hit"]
                thit = event[tracker]
                for residual in residual_dict.keys():
                    residual_dict[residual].append(thit[residual] - vhit[residual])
        return residual_dict

    def process_data_detector_residuals(self):
        for detector, virtual_station_list in self.mc_stations.iteritems():
            virtual_station = virtual_station_list[0]
            residual_dict = self.get_data_detector_residuals(detector, virtual_station)
            for var in sorted(residual_dict.keys()):
                canvas_name = "mc_residual_"+detector+"_"+var
                data = residual_dict[var]
                hist = self.plots[canvas_name]["histograms"][var]
                for item in data:
                    hist.Fill(item)

    def birth_data_detector_residuals(self):
        dummy_canvas = xboa.common.make_root_canvas("dummy")
        for detector, virtual_station_list in self.mc_stations.iteritems():
            virtual_station = virtual_station_list[0]
            residual_dict = self.get_data_detector_residuals(detector, virtual_station)
            for var in sorted(residual_dict.keys()):
                canvas_name = "mc_residual_"+detector+"_"+var
                data = residual_dict[var]
                if len(data) == 0:
                    raise RuntimeError("Failed to find residual data for "+var+" "+detector+" "+virtual_station)
                xmin, xmax = utilities.utilities.fractional_axis_range(data, 0.95)
                dummy_canvas.cd() # make sure we don't accidentally overwrite "current" canvas
                hist = xboa.common.make_root_histogram(var, data, self.axis_labels[var], 100, [], '', 0, [], xmin, xmax)
                hist.Draw()
                fit = utilities.utilities.fit_peak(hist, nsigma=8)
                mean = fit.GetParameter(1)
                sigma = fit.GetParameter(2)
                xmin, xmax = mean-5*sigma, mean+5*sigma
                hist = self.make_root_histogram(canvas_name, var, data, self.axis_labels[var], 100, [], '', 0, [], xmin, xmax)
                hist.Draw()

    def death_data_detector_residuals(self):
        for name in self.plots:
            if "mc_residual_" not in name:
                continue
            self.plots[name]["canvas"].cd()
            for hist_name in self.plots[name]["histograms"]:
                hist = self.plots[name]["histograms"][hist_name]
                fit = utilities.utilities.fit_peak(hist, nsigma=1)
                text_box = utilities.utilities.get_text_box(self.config, self.config_anal, None, fit, hist)

    @staticmethod
    def do_mc_plots(config, config_anal, data_loader):
        plotter = MCPlotter(config, config_anal, data_loader)
        pid_colors = {211:2, -13:4, -11:8, +11:ROOT.kGray}
        for detector in ["tkd", "tku"]:
            plotter.plot_amplitude_residuals(detector)



