import shutil
import os

import ROOT

class CompareConfig(object):
    def __init__(self):
        pass
  
    def setup(self, beam, run_dir, src_plot_dir, target_plot_dir, will_do_data, will_do_mc):
        self.beam = beam
        self.name = beam.split("_")
        self.get_experiment_config()
        self.data_caption = ""
        self.mc_caption = ""
        self.cuts_tex = "cut_plots/cuts_summary_*.tex"
        self.conglomerate_dir = []
        if will_do_data:
            self.conglomerate_dir.append(run_dir+"plots_"+beam+"/")
        if will_do_mc:
            self.conglomerate_dir.append(run_dir+"plots_Simulated_"+beam+"/")
        for a_dir in self.conglomerate_dir:
            if not os.path.exists(a_dir):
                raise RuntimeError("Could not find "+str(a_dir))
        self.src_plot_dir = src_plot_dir
        self.target_plot_dir = run_dir+target_plot_dir
        plot_dir = self.target_plot_dir+beam+"/"
        if os.path.exists(plot_dir):
            shutil.rmtree(plot_dir)
        os.makedirs(plot_dir)
        self.beam_plot_dir = plot_dir

    def get_experiment_config(self):
        if len(self.name) > 2:
            self.channel = self.name[0]
            self.beamline = self.name[1]
            self.absorber = " ".join(self.name[2:])
            if self.absorber == "None":
                self.absorber = "No Absorber"
        else:
            self.channel = ""
            self.beamline = ""
            self.absorber = ""

    def get_conglomerate_1(self, canvas, hist, cut, axis_title, axis_range, normalise, legend_pos, rescale_y = True, modifiers = {}):
        hist_name = hist
        if cut != None:
            hist_name = hist+" "+cut
        y_range = None
        if type(rescale_y) == type([]):
            y_range = rescale_y
        my_config = {
                "file_name":canvas,
                "canvas_name":canvas,
                "histogram_names":[hist_name],
                "graph_names":[],
                "rescale_x":axis_range,
                "rescale_y":rescale_y,
                "normalise_hist":normalise,
                "normalise_graph":False,
                "calculate_errors":{
                    "histograms":[0],
                    "normalised":True,
                },
                "rebin":False,
                "mice_logo":True,
                "log_y":False,
                "hist_title":"",
                "replace_hist":False,
                "redraw":{
                    "draw_option":["P E1 PLC", ""],
                    "fill_color":[1, ROOT.kOrange-2],
                    "transparency":None,
                    "line_color":[1, 1],
                    "marker_style":[20, 20],
                    "draw_order":[1, 0],
                    "x_range":axis_range,
                    "y_range":y_range,
                    "graph_draw_option":["P", "P"],
                    "ignore_more_histograms":False,
                },
                "canvas_fill_color":None,
                "extra_lines":False,
                "extra_labels":False,
                "legend":{
                    "text":["data", "simulation"],
                    "draw_option":["p e1", "f l"],
                    "pos":legend_pos,
                },
                "defit":True,
                "write_plots":{
                    "dir":self.beam_plot_dir,
                    "file_name":None,
                    "formats":["png", "root", "eps", "pdf"],
                },
                "axis_title":{
                    "x":axis_title,
                    "y":None,
                },
            }
        self.recursive_modify_dict(my_config, modifiers)
        return my_config

    def get_conglomerate_2(self, canvas, hist_name, axis_title, axis_range, normalise, legend_pos, modifiers = {}):
        if hist_name == None:
            hist_name = canvas
        my_config = {
                "file_name":canvas,
                "canvas_name":canvas,
                "histogram_names":[hist_name],
                "graph_names":[],
                "rescale_x":axis_range,
                "rescale_y":True,
                "replace_hist":False,
                "normalise_hist":normalise,
                "normalise_graph":False,
                "calculate_errors":{
                    "histograms":[0],
                    "normalised":True,
                },
                "rebin":False,
                "mice_logo":True,
                "log_y":False,
                "hist_title":"",
                "redraw":{
                    "draw_option":["P E1 PLC", ""],
                    "fill_color":[1, ROOT.kOrange-2],
                    "transparency":None,
                    "line_color":[1, 1],
                    "marker_style":[20, 20],
                    "draw_order":[1, 0],
                    "x_range":axis_range,
                    "y_range":None,
                    "graph_draw_option":["P", "P"],
                    "ignore_more_histograms":False,
                },
                "canvas_fill_color":None,
                "extra_lines":False,
                "extra_labels":False,
                "legend":{
                    "text":["data", "simulation"],
                    "draw_option":["p e1", "f l"],
                    "pos":legend_pos,
                },
                "defit":True,
                "write_plots":{
                    "dir":self.beam_plot_dir,
                    "file_name":None,
                    "formats":["png", "root", "eps", "pdf"],
                },
                "axis_title":{
                    "x":axis_title,
                    "y":None,
                },
            }
        self.recursive_modify_dict(my_config, modifiers)
        return my_config

    def get_conglomerate_3(self, canvas, hist_name, x_axis_title, y_axis_title, modifiers = {}):
        my_config = {
                "file_name":canvas,
                "canvas_name":canvas,
                "histogram_names":[hist_name],
                "graph_names":[],
                "rescale_x":False,
                "rescale_y":False,
                "normalise_hist":False,
                "normalise_graph":False,
                "replace_hist":False,
                "calculate_errors":False,
                "rebin":False,
                "mice_logo":True,
                "log_y":False,
                "hist_title":"",
                "redraw":{
                    "draw_option":["COLZ"],
                    "fill_color":[1],
                    "transparency":None,
                    "line_color":[1],
                    "marker_style":[20],
                    "draw_order":[0],
                    "x_range":None,
                    "y_range":None,
                    "graph_draw_option":["P", "P"],
                    "ignore_more_histograms":False,
                },
                "canvas_fill_color":None,
                "legend":False,
                "extra_lines":False,
                "extra_labels":False,
                "defit":False,
                "write_plots":{
                    "dir":self.beam_plot_dir,
                    "file_name":None,
                    "formats":["png", "root", "eps", "pdf"],
                },
                "axis_title":{
                    "x":x_axis_title,
                    "y":y_axis_title,
                },
            }
        self.recursive_modify_dict(my_config, modifiers)
        return my_config

    def recursive_modify_dict(self, target_config, source_config):
        for key, value in source_config.iteritems():
            if key not in target_config:
                target_config[key] = value
            elif type(target_config[key]) != type(value):
                target_config[key] = value
            elif type(target_config[key]) == type({}):
                self.recursive_modify_dict(target_config[key], value)
            elif type(target_config[key]) == type([]):
                self.recursive_modify_list(target_config[key], value)
            else:
                target_config[key] = value

    def recursive_modify_list(self, target_config, source_config):
        for i, value in enumerate(source_config):
            if i >= len(target_config):
                target_config.append(value)
            elif type(target_config[i]) != type(value):
                target_config[i] = value
            elif type(target_config[i]) == type({}):
                self.recursive_modify_dict(target_config[i], value)
            elif type(target_config[i]) == type([]):
                self.recursive_modify_list(target_config[i], value)
            else:
                target_config[i] = value
                   
    def get_conglomerate_graph(self, file_name, x_axis_title, y_axis_title, canvas_name = None, 
                               hist_list = None, graph_list = None, x_range = None, y_range = None, replace_hist = False,
                               graph_draw_option = None, graph_marker_style = None, graph_marker_color = None, graph_draw_order = None,
                               modifiers = {}):
        if canvas_name == None:
            canvas_name = file_name
        if hist_list == None:
            hist_list = [canvas_name]
        if graph_list == None:
            graph_list = [canvas_name]
        if graph_draw_option == None:
            graph_draw_option = ["p"]*len(graph_list)*len(self.conglomerate_dir)
        
        my_config = {
                "file_name":file_name,
                "canvas_name":canvas_name,
                "histogram_names":hist_list,
                "graph_names":graph_list,
                "rescale_x":x_range,
                "rescale_y":y_range,
                "replace_hist":replace_hist,
                "normalise_hist":False,
                "normalise_graph":False,
                "calculate_errors":{
                    "histograms":[0],
                    "normalised":True,
                },
                "rebin":False,
                "mice_logo":True,
                "log_y":False,
                "hist_title":"",
                "extra_lines":False,
                "extra_labels":False,
                "redraw":{
                    "draw_option":[""],
                    "fill_color":[1],
                    "transparency":None,
                    "line_color":[1],
                    "marker_style":[20],
                    "draw_order":[0],
                    "x_range":None,
                    "y_range":y_range,
                    "ignore_more_histograms":True,
                    "graph":{
                        "draw_option":graph_draw_option,
                        "marker_style":graph_marker_style,
                        "marker_color":graph_marker_color,
                        "fill_color":None,
                        "fill_style":None,
                        "transparency":None,
                        "draw_order":graph_draw_order,
                    }
                },
                "canvas_fill_color":None,
                "legend":False,
                "extra_lines":False,
                "defit":True,
                "write_plots":{
                    "dir":self.beam_plot_dir,
                    "file_name":None,
                    "formats":["png", "root", "eps", "pdf"],
                },
                "axis_title":{
                    "x":x_axis_title,
                    "y":y_axis_title,
                },
            }
        self.recursive_modify_dict(my_config, modifiers)
        return my_config
