import xboa.common
import json
import copy
import numpy

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

    def plot_detector_var_one_d(self, detector, slice_variable, plot_variable, color_dict):
        track_final = {}
        for event in self.data_loader.events:
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
        hist = xboa.common.make_root_histogram("all", all_list, plot_variable, 100)
        hist.SetTitle(plot_variable+" at "+detector)
        hist.SetName("all")
        hist.Draw()
        for key in sorted(track_final.keys()):
            var_list = track_final[key]
            name = plot_variable+" = "+str(key)
            hist = xboa.common.make_root_histogram(name, var_list, plot_variable, 100)
            hist.SetName(name)
            hist.SetMarkerStyle(24)
            if key in color_dict:
                hist.SetMarkerColor(color_dict[key])
            hist.Draw("SAMEP")
        canvas.BuildLegend()
        canvas.Update()
        for format in ["png", "pdf", "root"]:
            canvas.Print(self.plot_dir+"/"+plot_variable+"_at_"+detector+"."+format)

    def plot_detector_var_two_d_scatter(self, detector, slice_variable, plot_variable_1, plot_variable_2, color_dict):
        track_final = {}
        for event in self.data_loader.events:
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
        for format in ["png", "pdf", "root"]:
            title = title.replace(" ", "_")
            canvas.Print(self.plot_dir+"/"+title+"."+format)

    def plot_tracks_final(self):
        track_final = {}
        for event in self.data_loader.events:
            for detector_hit in event["data"]:
                if detector_hit["detector"] != "mc_track_final":
                    continue
                pid = detector_hit["hit"]["pid"]
                x = detector_hit["hit"]["x"]
                y = detector_hit["hit"]["y"]
                y = detector_hit["hit"]["r"]
                z = detector_hit["hit"]["z"]
                if pid not in track_final:
                    track_final[pid] = ([], [], [], [])
                track_final[pid][0].append(z)
                track_final[pid][1].append(x)
                track_final[pid][2].append(y)
        pid_colors = {211:2, -13:4, -11:8}

      
        xmin, xmax, ymin, ymax = 0., 0., 0., 0.
        for z_list, x_list, y_list, p_list in track_final.values():
            xmin = min([xmin]+x_list)
            xmax = max([xmax]+x_list)
            ymin = min([ymin]+y_list)
            ymax = max([ymax]+y_list)

        xz_canvas, yz_canvas = None, None
        for pid in sorted(track_final.keys()):
            z_list = track_final[pid][0]
            x_list = track_final[pid][1]
            y_list = track_final[pid][2]
            x_hist, x_graph = xboa.common.make_root_graph(str(pid), z_list, "z [mm]", x_list, "x [mm]", xmin=0., xmax=20000., ymin=xmin, ymax=xmax)
            y_hist, y_graph = xboa.common.make_root_graph(str(pid), z_list, "z [mm]", y_list, "y [mm]", xmin=0., xmax=20000., ymin=ymin, ymax=ymax)
            if xz_canvas == None:
                xz_canvas = xboa.common.make_root_canvas("xz of tracks at death")
                x_hist.Draw()
                yz_canvas = xboa.common.make_root_canvas("yz of tracks at death")
                y_hist.Draw()
            if pid in pid_colors:
                x_graph.SetMarkerColor(pid_colors[pid])
                y_graph.SetMarkerColor(pid_colors[pid])
            x_graph.SetMarkerStyle(6)
            y_graph.SetMarkerStyle(6)
            xz_canvas.cd()
            x_graph.Draw("PSAME")
            yz_canvas.cd()
            y_graph.Draw("PSAME")
        xz_canvas.Update()
        yz_canvas.Update()
        for format in ["png", "pdf", "root"]:
            xz_canvas.Print(self.plot_dir+"/xz_track_position."+format)
            yz_canvas.Print(self.plot_dir+"/yz_track_position."+format)

    def plot_detector_residuals(self):
        pass

    @staticmethod
    def do_mc_plots(config, config_anal, data_loader):
        plotter = MCPlotter(config, config_anal, data_loader)
        #plotter.plot_tracks_final()
        pid_colors = {211:2, -13:4, -11:8}
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


