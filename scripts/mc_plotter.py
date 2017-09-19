import xboa.common
import json
import copy
import numpy
import sys

import cdb_tof_triggers_lookup
import ROOT
from xboa.bunch import Bunch
import scripts.utilities

from analysis_base import AnalysisBase

class MCPlotter(AnalysisBase):
    def __init__(self, config, config_anal, data_loader):
        super(MCPlotter, self).__init__(config, config_anal, data_loader)
        self.process_args = {}
        self.mc_stations = {}
        self.failed_pids = {}

    def birth(self):
        for key in self.config.mc_plots["mc_stations"].keys():
            self.mc_stations[key] = self.config.mc_plots["mc_stations"][key][0]
        pid_colors = {211:2, -13:4, -11:8, +11:ROOT.kGray}
        for detector in ["tkd", "tku"]:
            station = self.config.mc_plots["mc_stations"][detector][0]
            virt_station = "mc_virtual_"+str(station)
            self.birth_var_one_d("p_at_"+virt_station, virt_station, "pid", "p", pid_colors, True, 0., 300.)
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
        self.birth_data_detector_residuals()

    def process(self):
        for name in sorted(self.process_args.keys()):
            process_function = self.process_args[name][0]
            process_args = self.process_args[name][1]
            process_function(name, *process_args)
        self.process_data_detector_residuals()

    def death(self):
        self.death_data_detector_residuals()
        print "Failed to put the following pids into scatter plots:"
        print "   ", self.failed_pids
        self.print_plots()

    def get_data_var_one_d(self, name, *args):
        track_final = {}
        for event in self.data_loader.events:
            if args[5] and event["any_cut"]:
                continue
            for detector_hit in event["data"]:
                if detector_hit["detector"] != args[0]:
                    continue
                slice_var = detector_hit["hit"][args[1]]
                plot_var = detector_hit["hit"][args[2]]
                if slice_var not in track_final:
                    track_final[slice_var] = []
                track_final[slice_var].append(plot_var)
        all_list = []
        for value in track_final.values():
            all_list += value
        return (all_list, track_final)

    def process_var_one_d(self, name, *args):
        all_list, track_final = self.get_data_var_one_d(name, *args)
        for item in all_list:
            self.plots[name]["histograms"]["all"].Fill(item)
        for i, key in enumerate(sorted(track_final.keys())):
            hist_dict = self.plots[name]["histograms"]
            plot_name = str(args[1])+" = "+str(key)
            if plot_name not in hist_dict:
                if key not in self.failed_pids:
                    self.failed_pids[key] = 0
                self.failed_pids[key] += 1
                continue
            for item in track_final[key]:
                hist_dict[plot_name].Fill(item)

    def birth_var_one_d(self, canvas_name, detector, slice_variable, plot_variable, color_dict, cuts = False, xmin = None, xmax = None):
        all_list, track_final = self.get_data_var_one_d(canvas_name, detector, slice_variable, plot_variable, color_dict, cuts, xmin, xmax)
        hist = self.make_root_histogram(canvas_name, "all", all_list, plot_variable+" at "+detector, 100, [], '', 0, [], xmin, xmax)
        hist.Draw("P")
        for i, key in enumerate(sorted(track_final.keys())):
            var_list = track_final[key]
            name = slice_variable+" = "+str(key)
            label = plot_variable+" at "+detector+" ["+scripts.utilities.default_units(plot_variable)+"]"
            hist = self.make_root_histogram(canvas_name, name, var_list, label, 100, [], '', 0, [], xmin, xmax)
            hist.SetMarkerStyle(24)
            if key in color_dict:
                hist.SetMarkerColor(color_dict[key])
            hist.Draw("SAMEP")
        self.plots[canvas_name]["canvas"].SetLogy()
        self.plots[canvas_name]["canvas"].BuildLegend()
        self.process_args[canvas_name] = [self.process_var_one_d, (detector, slice_variable, plot_variable, color_dict, cuts, xmin, xmax)]

    def get_data_var_two_d_scatter(self, name, *args):
        track_final = {}
        for event in self.data_loader.events:
            if args[5] and event["any_cut"]:
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
            label_1 = plot_variable_1+" ["+scripts.utilities.default_units(plot_variable_1)+"]"
            label_2 = plot_variable_1+" ["+scripts.utilities.default_units(plot_variable_2)+"]"
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
            xmin, xmax = scripts.utilities.fractional_axis_range(data, 0.95)
            hist = xboa.common.make_root_histogram(var, data, self.axis_labels[var], 100, xmin=xmin, xmax=xmax)
            hist.Draw()
            fit = scripts.utilities.fit_peak(hist, nsigma=8)
            mean = fit.GetParameter(1)
            sigma = fit.GetParameter(2)
            xmin, xmax = mean-5*sigma, mean+5*sigma
            self.axis_min_max[var] = (xmin, xmax)
            hist = xboa.common.make_root_histogram(var, data, self.axis_labels[var], 100, xmin=xmin, xmax=xmax)
            hist.Draw()
            fit = scripts.utilities.fit_peak(hist, nsigma=1)

            text_box = scripts.utilities.get_text_box(self.config, self.config_anal, data, fit)
            canvas.Update()
            for format in ["png", "eps", "root"]:
                canvas.Print(self.plot_dir+"/mc_residual_"+suffix+"_"+var+"."+format)

    def get_data_detector_residuals(self, tracker, station):
        residual_dict = {"x":[], "y":[], "z":[], "px":[], "py":[], "pz":[]}
        virt_name = "mc_virtual_"+str(station)
        for event in self.data_loader.events:
            if event["any_cut"]:
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
        for detector, virtual_station in self.mc_stations.iteritems():
            residual_dict = self.get_data_detector_residuals(detector, virtual_station)
            for var in sorted(residual_dict.keys()):
                canvas_name = "mc_residual_"+detector+"_"+var
                data = residual_dict[var]
                hist = self.plots[canvas_name]["histograms"][var]
                for item in data:
                    hist.Fill(item)

    def birth_data_detector_residuals(self):
        for detector, virtual_station in self.mc_stations.iteritems():
            residual_dict = self.get_data_detector_residuals(detector, virtual_station)
            for var in sorted(residual_dict.keys()):
                canvas_name = "mc_residual_"+detector+"_"+var
                data = residual_dict[var]
                xmin, xmax = scripts.utilities.fractional_axis_range(data, 0.95)
                hist = xboa.common.make_root_histogram(var, data, self.axis_labels[var], 100, [], '', 0, [], xmin, xmax)
                hist.Draw()
                fit = scripts.utilities.fit_peak(hist, nsigma=8)
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
                fit = scripts.utilities.fit_peak(hist, nsigma=1)
                text_box = scripts.utilities.get_text_box(self.config, self.config_anal, None, fit, hist)

    @staticmethod
    def do_mc_plots(config, config_anal, data_loader):
        plotter = MCPlotter(config, config_anal, data_loader)
        pid_colors = {211:2, -13:4, -11:8, +11:ROOT.kGray}
        for detector in ["tkd", "tku"]:
            plotter.plot_amplitude_residuals(detector)



