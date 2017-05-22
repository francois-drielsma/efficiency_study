import xboa.common
import json
import copy
import numpy
import sys

import cdb_tof_triggers_lookup
import ROOT
from xboa.bunch import Bunch
import scripts.utilities

class MCPlotter(object):
    def __init__(self, config, config_anal, data_loader):
        self.data_loader = data_loader
        self.config = config
        self.config_anal = config_anal
        self.plot_dir = config_anal["plot_dir"]

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

    def plot_detector_var_one_d(self, detector, slice_variable, plot_variable, color_dict, cuts = False, xmin = None, xmax = None):
        track_final = {}
        for event in self.data_loader.events:
            if cuts and event["any_cut"]:
                continue
            for detector_hit in event["data"]:
                if detector_hit["detector"] != detector:
                    continue
                slice_var = detector_hit["hit"][slice_variable]
                plot_var = detector_hit["hit"][plot_variable]
                if slice_var not in track_final:
                    track_final[slice_var] = []
                track_final[slice_var].append(plot_var)
        all_list = []
        for value in track_final.values():
            all_list += value

        canvas = xboa.common.make_root_canvas(plot_variable+" at "+detector)
        canvas.SetLogy()
        hist = xboa.common.make_root_histogram("all", all_list, plot_variable+" at "+detector, 100, xmin=xmin, xmax=xmax)
        hist.SetTitle(self.config_anal['name'])
        hist.SetName("all")
        hist.Draw()
        for key in sorted(track_final.keys()):
            var_list = track_final[key]
            name = slice_variable+" = "+str(key)
            hist = xboa.common.make_root_histogram(name, var_list, plot_variable+" at "+detector, 100, xmin=xmin, xmax=xmax)
            hist.SetName(name)
            hist.SetMarkerStyle(24)
            if key in color_dict:
                hist.SetMarkerColor(color_dict[key])
            hist.Draw("SAMEP")
        canvas.BuildLegend()
        canvas.Update()
        for format in ["png", "eps", "root"]:
            canvas.Print(self.plot_dir+"/"+plot_variable+"_at_"+detector+"."+format)

    def plot_detector_var_two_d_scatter(self, detector, slice_variable, plot_variable_1, plot_variable_2, color_dict, cuts=False):
        track_final = {}
        for event in self.data_loader.events:
            if cuts and event["any_cut"]:
                continue
            for detector_hit in event["data"]:
                if detector_hit["detector"] != detector:
                    continue
                pid = detector_hit["hit"][slice_variable]
                x = detector_hit["hit"][plot_variable_1]
                y = detector_hit["hit"][plot_variable_2]
                if pid not in track_final:
                    track_final[pid] = ([], [])
                track_final[pid][0].append(x)
                track_final[pid][1].append(y)
        if len(track_final) == 0:
            print "No tracks for", detector
            return
        xmin, xmax, ymin, ymax = 0., 0., 0., 0.
        for x_list, y_list in track_final.values():
            xmin = min([xmin]+x_list)
            xmax = max([xmax]+x_list)
            ymin = min([ymin]+y_list)
            ymax = max([ymax]+y_list)

        title = plot_variable_1+" vs "+plot_variable_2+" at "+detector
        canvas = None
        for pid in sorted(track_final.keys()):
            x_list = track_final[pid][0]
            y_list = track_final[pid][1]
            hist, graph = xboa.common.make_root_graph(slice_variable+" = "+str(pid),
              x_list, plot_variable_1, y_list, plot_variable_2,
              xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
            if canvas == None:
                canvas = xboa.common.make_root_canvas(title)
                hist.SetTitle(title)
                hist.Draw()
            if pid in color_dict:
                graph.SetMarkerColor(color_dict[pid])
            graph.SetMarkerStyle(6)
            graph.Draw("PSAME")
        canvas.Update()
        for format in ["png", "eps", "root"]:
            title = title.replace(" ", "_")
            canvas.Print(self.plot_dir+"/"+title+"."+format)

    axis_labels = {"x":"x(meas) - x(true) [mm]", "y":"y(meas) - y(true) [mm]", "z":"z(meas) - z(true) [mm]",
                   "px":"p_{x}(meas) - p_{x}(true) [MeV/c]", "py":"p_{y}(meas) - p_{y}(true) [MeV/c]", "pz":"p_{z}(meas) - p_{z}(true) [MeV/c]",
                   "amp_res":"A_{#perp    }(meas) - A_{#perp    }(true) [mm]", "amp_virt":"A_{#perp    } [mm]"}
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

    root_objects = []
    def plot_amplitude_residuals(self, suffix):
        canvas = xboa.common.make_root_canvas("MC amplitude residual")
        canvas.Draw()
        amp_virt = self.residual_dict["amp_virt"]
        for var in self.residual_dict.keys():
            if var == "amp_virt":
                continue  
            res = self.residual_dict[var]
            ymin, ymax = self.axis_min_max[var]
            scatter_hist = xboa.common.make_root_histogram("amplitude residuals",
                                                      amp_virt, "A_{#perp} [mm]", 20,
                                                      res, self.axis_labels[var], 50,
                                                      xmin=0., xmax=100.,
                                                      ymin=ymin, ymax=ymax)
            hist = ROOT.TProfile("amplitude residuals", ";A_{#perp} [mm];#delta A_{#perp} [mm]", 20, 0., 100., ymin, ymax)
            hist.SetStats(False)
            hist.BuildOptions(-10., 10., "s")
            self.root_objects.append(hist)
            for i in range(len(amp_virt)):
                hist.Fill(amp_virt[i], res[i])
            scatter_hist.Draw("COLZ")
            hist.Draw("SAME")
            canvas.Update()
            for format in ["png", "eps", "root"]:
                canvas.Print(self.plot_dir+"/mc_residual_"+suffix+"_"+var+"_amplitude_2d."+format)

    # tku station is 53
    # tkd station is 63
    def get_detector_residuals(self, tracker, station):
        self.residual_dict = {"x":[], "y":[], "z":[], "px":[], "py":[], "pz":[], "amp_virt":[], "amp_res":[]}
        virtual_hits = []
        real_hits = []
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
                virtual_hits.append(vhit)
                real_hits.append(thit)
        real_bunch = Bunch.new_from_hits(real_hits)
        virtual_bunch = Bunch.new_from_hits(virtual_hits)
        real_bunch.clear_weights()
        virtual_bunch.clear_weights()
        if real_bunch.bunch_weight() < 1e-9:
            print "Bunch weight 0"
            return
        real_bunch.set_covariance_matrix()
        virtual_bunch.set_covariance_matrix()
        print "Doing amplitudes"
        for i in range(len(virtual_bunch)):
            virtual_amp = Bunch.get_amplitude(virtual_bunch, virtual_bunch[i], ["x", "y"])
            real_amp = Bunch.get_amplitude(real_bunch, real_bunch[i], ["x", "y"])
            self.residual_dict["amp_virt"].append(virtual_amp)
            self.residual_dict["amp_res"].append(real_amp - virtual_amp)             
            for residual in ["x", "y", "z", "px", "py", "pz"]:
                self.residual_dict[residual].append(real_bunch[i][residual] - virtual_bunch[i][residual])

    @staticmethod
    def do_mc_plots(config, config_anal, data_loader):
        plotter = MCPlotter(config, config_anal, data_loader)
        #plotter.plot_tracks_final()
        pid_colors = {211:2, -13:4, -11:8, +11:ROOT.kGray}
        for detector in ["tkd", "tku"]:
            station = config.mc_plots["mc_stations"][detector][0]
            virt_station = "mc_virtual_"+str(station)
            plotter.get_detector_residuals(detector, station)
            plotter.plot_detector_residuals(detector)
            plotter.plot_amplitude_residuals(detector)
            plotter.plot_detector_var_one_d(virt_station, "pid", "p", pid_colors, True, 0., 300.)
            plotter.plot_detector_var_two_d_scatter(virt_station, "pid", "x", "px", pid_colors, True)

        plotter.plot_detector_var_one_d("mc_tof_0", "pid", "p", pid_colors)
        plotter.plot_detector_var_one_d("mc_tof_1", "pid", "p", pid_colors)
        plotter.plot_detector_var_one_d("mc_tof_0", "pid", "e_dep", pid_colors)
        plotter.plot_detector_var_one_d("mc_tof_1", "pid", "e_dep", pid_colors)
        plotter.plot_detector_var_two_d_scatter("mc_track_final", "pid", "z", "x", pid_colors)
        plotter.plot_detector_var_two_d_scatter("mc_track_final", "pid", "z", "y", pid_colors)
        plotter.plot_detector_var_two_d_scatter("mc_track_final", "pid", "z", "r", pid_colors)
        plotter.plot_detector_var_two_d_scatter("mc_primary", "pid", "x", "px", pid_colors)
        plotter.plot_detector_var_two_d_scatter("mc_primary", "pid", "y", "py", pid_colors)
        plotter.plot_detector_var_one_d("mc_primary", "pid", "p", pid_colors)
        plotter.plot_detector_var_one_d("mc_primary", "pid", "x", pid_colors)
        plotter.plot_detector_var_one_d("mc_primary", "pid", "y", pid_colors)
        plotter.plot_detector_var_one_d("mc_track_final", "pid", "p", pid_colors)
        plotter.plot_detector_var_one_d("mc_track_final", "pid", "x", pid_colors)
        plotter.plot_detector_var_one_d("mc_track_final", "pid", "y", pid_colors)




