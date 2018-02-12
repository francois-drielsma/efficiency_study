import sys
import subprocess
import shutil
import os
import glob

import xboa.common
import ROOT


class ResidualConfig(object):
    residual_conglomerate_dir = ["output/2017-12-07_2017-02_reco/plots_*"]
    residual_conglomerate_list = [{
            "file_name":"global_through_residual_tof2_t",
            "canvas_name":"global_through residual tof2 t",
            "histogram_names":["res ex cut"],
            "rescale_x":[-2., 2.],
            "rescale_y":True,
            "log_y":True,
            "hist_title":"",
        }, {
            "file_name":"global_through_residual_tkd_tp_p",
            "canvas_name":"global_through residual tkd_tp p",
            "histogram_names":["res ex cut"],
            "rescale_x":[-5., 5.],
            "rescale_y":True,
            "log_y":True,
            "hist_title":"",
        }, {
            "file_name":"global_through_residual_tof0_t",
            "canvas_name":"global_through residual tof0 t",
            "histogram_names":["res ex cut"],
            "rescale_x":[-2., 2.],
            "rescale_y":True,
            "log_y":True,
            "hist_title":"",
        },]

def get_conglomerate(canvas, hist, cut, axis_range, normalise, plot_dir):
    hist_name = hist
    if cut != None:
        hist_name = hist+" "+cut
    return {
            "file_name":canvas,
            "canvas_name":canvas,
            "histogram_names":[hist_name],
            "rescale_x":axis_range,
            "rescale_y":True,
            "normalise":normalise,
            "calculate_errors":{
                "histograms":[0],
                "normalised":True,
            },
            "rebin":False,
            "mice_logo":False,
            "log_y":False,
            "hist_title":"",
            "redraw":{
                "draw_option":["P E1 PLC", ""],
                "fill_color":[1, ROOT.kOrange-2],
                "line_color":[1, 1],
                "marker_style":[20, 20],
                "draw_order":[1, 0],
            },
            "legend":{
                "text":["data", "simulation"],
                "draw_option":["p e1", "f l"],
            },
            "defit":True,
            "write_plots":{
                "dir":plot_dir,
                "formats":["png", "root", "eps", "pdf"],
            }
        }

class MCConfig(object):
    def __init__(self, mc_version):
        self.conglomerate_dir = ["output/2017-02_reco_weighting/plots_10-140-10052/", "output/2017-02_mc_10052_"+mc_version+"/plots_10-140-empty/"]
        plot_dir = "output/data_mc_comparison_"+mc_version+"/"
        if os.path.exists(plot_dir):
            shutil.rmtree(plot_dir)
        os.makedirs(plot_dir)
        self.conglomerate_list = [
            get_conglomerate("tku_p", "tku_p", "all", [80., 200.], False, plot_dir),
            get_conglomerate("tof01", "tof01", "all", [24., 40.,], True, plot_dir),
            get_conglomerate("tof01", "tof01", "us cut", [28., 32.,], True, plot_dir),
            get_conglomerate("tku_p", "tku_p", "us cut", [130., 160.], False, plot_dir),
            get_conglomerate("tku_x", "tku_x", "us cut", None, False, plot_dir),
            get_conglomerate("tku_y", "tku_y", "us cut", None, False, plot_dir),
            get_conglomerate("tku_px", "tku_px", "us cut", None, False, plot_dir),
            get_conglomerate("tku_py", "tku_py", "us cut", None, False, plot_dir),
            get_conglomerate("chi2_tku", "chi2", "us cut", None, True, plot_dir),
            get_conglomerate("chi2_tkd", "chi2", "us cut", None, True, plot_dir),
            get_conglomerate("tof0_x", "tof0_x", "us cut", None, False, plot_dir),
            get_conglomerate("tof0_y", "tof0_y", "us cut", None, False, plot_dir),
            get_conglomerate("tof1_x", "tof1_x", "us cut", None, False, plot_dir),
            get_conglomerate("tof1_y", "tof1_y", "us cut", None, False, plot_dir),
        ]


def get_conglomerate_2(canvas, hist, cut, axis_range, normalise, plot_dir):
    hist_name = hist
    if cut != None:
        hist_name = hist+" "+cut
    return {
            "file_name":canvas,
            "canvas_name":canvas,
            "histogram_names":[hist_name],
            "graph_names":[hist_name+"_graph"],
            "rescale_x":axis_range,
            "rescale_y":True,
            "normalise":normalise,
            "calculate_errors":{
                "histograms":[0],
                "normalised":True,
            },
            "rebin":False,
            "mice_logo":False,
            "log_y":False,
            "hist_title":"",
            "redraw":{
                "draw_option":["P E1 PLC", ""],
                "fill_color":[1, ROOT.kOrange-2],
                "line_color":[1, 1],
                "marker_style":[20, 20],
                "draw_order":[1, 0],
            },
            "legend":{
                "text":["Recon data", "Recon sim"],
                "draw_option":["p e1", "f l"],
            },
            "defit":True,
            "write_plots":{
                "dir":plot_dir,
                "formats":["png", "root", "eps", "pdf"],
            }
        }


class CompareCutsConfig(object):
    def __init__(self, beam):
        self.beam = beam
        self.conglomerate_dir = ["output/2017-02/plots_"+beam+"/",
                                 "output/2017-02/plots_Simulated_"+beam+"/"]
        for a_dir in self.conglomerate_dir:
            if not os.path.exists(a_dir):
                raise RuntimeError("Could not find "+str(a_dir))
        self.plot_dir = "/cut_plots/"
        self.cuts_tex = "/data_plots/cuts_summary.tex"
        plot_dir = "output/2017-02/compare_cuts/"+beam+"/"
        if os.path.exists(plot_dir):
            shutil.rmtree(plot_dir)
        os.makedirs(plot_dir)
        self.conglomerate_list = [
            get_conglomerate_2("global_through_virtual_diffuser_ds_r_us_cut", "global_through_virtual_diffuser_ds_r_us_cut", None, None, True, plot_dir),
            get_conglomerate_2("global_through_virtual_diffuser_us_r_us_cut", "global_through_virtual_diffuser_us_r_us_cut", None, None, True, plot_dir),
            get_conglomerate_2("tkd_n_tracks_ds_cut", "tkd_n_tracks_ds_cut", None, None, True, plot_dir),
            get_conglomerate_2("tkd_chi2_ds_cut", "tkd_chi2_ds_cut", None, None, True, plot_dir),
            get_conglomerate_2("tkd_max_r_ds_cut", "tkd_max_r_ds_cut", None, None, True, plot_dir),
            get_conglomerate_2("tkd_p_ds_cut", "tkd_p_ds_cut", None, None, True, plot_dir),
            get_conglomerate_2("tku_chi2_us_cut", "tku_chi2_us_cut", None, None, True, plot_dir),
            get_conglomerate_2("tku_max_r_us_cut", "tku_max_r_us_cut", None, None, True, plot_dir),
            get_conglomerate_2("tku_n_clusters_us_cut", "tku_n_clusters_us_cut", None, None, True, plot_dir),
            get_conglomerate_2("tku_n_tracks_us_cut", "tku_n_tracks_us_cut", None, None, True, plot_dir),
            get_conglomerate_2("tku_p_us_cut", "tku_p_us_cut", None, None, True, plot_dir),
            get_conglomerate_2("tof_tof0_n_sp_us_cut", "tof_tof0_n_sp_us_cut", None, None, True, plot_dir),
            get_conglomerate_2("tof_tof1_n_sp_us_cut", "tof_tof1_n_sp_us_cut", None, None, True, plot_dir),
            get_conglomerate_2("tof_tof01_us_cut", "tof_tof01_us_cut", None, None, True, plot_dir),
            get_conglomerate_2("tof_delta_tof01_us_cut", "tof_delta_tof01_us_cut", None, None, True, plot_dir),
        ]
        self.data_caption = """The sample selection criteria are listed sequentially, in the order 
that the selection was made. The number of events surviving the 
selection and all preceding selections is listed for each of the data runs studied
in this note."""
        self.mc_caption = "Simulated sample selection"

class CompareGlobalsConfig(object):
    def __init__(self, beam):
        self.beam = beam
        self.conglomerate_dir = ["output/2017-02/plots_"+beam+"/",
                                 "output/2017-02/plots_Simulated_"+beam+"/"]
        for a_dir in self.conglomerate_dir:
            if not os.path.exists(a_dir):
                raise RuntimeError("Could not find "+str(a_dir))
                
        self.plot_dir = "/global_plots/"
        plot_dir = "output/2017-02/compare_globals/"+beam+"/"
        if os.path.exists(plot_dir):
            shutil.rmtree(plot_dir)
        os.makedirs(plot_dir)
        self.conglomerate_list = [
            get_conglomerate_2("global_through_residual_tkd_tp_p", "res_ex_cut", None, [-20, 20], False, plot_dir),
            get_conglomerate_2("global_through_residual_tkd_tp_px", "res_ex_cut", None, None, False, plot_dir),
            get_conglomerate_2("global_through_residual_tkd_tp_py", "res_ex_cut", None, None, False, plot_dir),
            get_conglomerate_2("global_through_residual_tkd_tp_pz", "res_ex_cut", None, [-20, 20], False, plot_dir),
            get_conglomerate_2("global_through_residual_tkd_tp_x", "res_ex_cut", None, None, False, plot_dir),
            get_conglomerate_2("global_through_residual_tkd_tp_y", "res_ex_cut", None, None, False, plot_dir),
            get_conglomerate_2("global_through_residual_tof0_t", "res_ex_cut", None, [-5, 5], False, plot_dir),
            get_conglomerate_2("global_through_residual_tof1_x", "res_ex_cut", None, None, False, plot_dir),
            get_conglomerate_2("global_through_residual_tof1_y", "res_ex_cut", None, None, False, plot_dir),
            get_conglomerate_2("global_through_residual_tof2_x", "res_ex_cut", None, None, False, plot_dir),
            get_conglomerate_2("global_through_residual_tof2_y", "res_ex_cut", None, None, False, plot_dir),
        ]
        for item in self.conglomerate_list:
            if "tof0_t" in item["file_name"]:
                continue
            item["rebin"] = 4
        print self.conglomerate_list[0]["rebin"]

class ConglomerateBase(object):
    def __init__(self, config):
        self.config = config
    
    def conglomerate(self):
        for item in self.config.conglomerate_list:
            cong = ConglomerateOne(item, self.config)
            cong.dir_path = [a_dir+self.config.plot_dir for a_dir in self.config.conglomerate_dir]
            cong.conglomerate()

class ConglomerateOne(object):
    def __init__(self, options, config):
        self.options = options
        self.config = config

    def get_hist_list(self, canvas):
        hist_list = []
        for an_object in canvas.GetListOfPrimitives():
            name = str(an_object.GetName()).replace(" ", "_")
            for hist_name in self.options["histogram_names"]:
                if hist_name in name:
                    if type(an_object) == type(ROOT.TGraph()):
                        continue
                    hist_list.append(an_object)
                    break # we can only add each object once
        if len(hist_list) == 0:
            print "Did not find a histogram from list", self.options["histogram_names"]
            for an_object in canvas.GetListOfPrimitives():
                print an_object.GetName()
            raise KeyError("Failed to find histogram")
        return hist_list

    def get_graph_list(self, canvas):
        graph_list = []
        for an_object in canvas.GetListOfPrimitives():
            name = str(an_object.GetName()).replace(" ", "_")
            for graph_name in self.options["graph_names"]:
                if graph_name in name:
                    if type(an_object) != type(ROOT.TGraph()):
                        continue
                    graph_list.append(an_object)
                    break # we can only add each object once
        return graph_list

    def get_canvas(self, file_name):
        fin = ROOT.TFile(file_name, "READ")
        a_list = fin.GetListOfKeys()
        for key in a_list:
            old_canvas = key.ReadObj()
            target_name = self.options["canvas_name"].replace(" ", "_")
            canvas_name = (old_canvas.GetName()).replace(" ", "_")
            if target_name in canvas_name:
                return old_canvas
        print "Keys found in file", file_name, "while searching for", self.options["canvas_name"]
        for key in a_list:
            print key.ReadObj().GetName()
        raise KeyError("Failed to find canvas")

    def defit(self, canvas, hist_list, graph_list):
        if not self.options["defit"]:
            return
        for hist in hist_list:
            for func in hist.GetListOfFunctions():
                func.SetRange(-1e9, -0.9e9)
            #function.SetRange(-1e9, -0.9e9)
            #hist.Fit("gaus", "", "", -1e9, -0.9e9)

    def rescale_x(self, canvas, hist_list, graph_list):
        if not self.options["rescale_x"]:
            return
        x_range = self.options["rescale_x"]
        for hist in hist_list:
            hist.GetXaxis().SetRangeUser(x_range[0], x_range[1])

    def rescale_y(self, canvas, hist_list, graph_list):
        if not self.options["rescale_y"]:
            return
        hist = hist_list[0]
        min_bin, max_bin = hist.GetMinimum(), hist.GetMaximum()
        for hist in hist_list[1:]:
            min_bin = max(min_bin, hist.GetMinimum())
            max_bin = max(max_bin, hist.GetMaximum())
        if self.options["log_y"]:
            min_bin = max(1e-4, min_bin)
        for hist in hist_list:
            hist.GetYaxis().SetRangeUser(min_bin, max_bin*1.1)

    def calculate_errors(self, canvas, hist_list, graph_list):
        if not self.options["calculate_errors"]:
            return
        for i in self.options["calculate_errors"]["histograms"]:
            hist = hist_list[i]
            for index in range(hist.GetNbinsX()+2):
                if hist.GetBinContent(index) < 1e-6:
                    continue
                err = hist.GetBinContent(index)**0.5
                if self.options["calculate_errors"]["normalised"]:
                    err = err/hist.GetEntries()**0.5
                hist.SetBinError(index, err)
            print

    def normalise(self, canvas, hist_list, graph_list):
        if not self.options["normalise"]:
            return
        for hist in hist_list:
            try:
                n_entries = hist.GetEntries()
            except AttributeError: # not a histogram (maybe a graph?)
                continue
            if n_entries == 0:
                continue
            hist.Scale(1./n_entries)

    def rebin(self, canvas, hist_list, graph_list):
        if not self.options["rebin"]:
            return
        for hist in hist_list:
            hist.Rebin(self.options["rebin"])


    def redraw(self, canvas, hist_list, graph_list):
        redraw = self.options["redraw"]
        draw_option = ["" for hist in hist_list]
        draw_order = [i for i, hist in enumerate(hist_list)]
        if redraw:
            if len(hist_list) != len(redraw["line_color"]):
                print "Failed to find all the histograms for redraw(...); found", len(hist_list), "expected", len(redraw["line_color"])
                raise RuntimeError("Failed to find all histograms for redraw")
                return
            for i, hist in enumerate(hist_list):
                if type(hist) == type(ROOT.TGraph()):
                    continue
                hist.SetLineColor(redraw["line_color"][i])
                hist.SetFillColor(redraw["fill_color"][i])
                hist.SetMarkerStyle(redraw["marker_style"][i])
            draw_option = redraw["draw_option"]
            draw_order = redraw["draw_order"]
        same = ""
        for i in draw_order:
            hist_list[i].Draw(same+draw_option[i])
            same = "SAME "
        for graph in graph_list:
            graph.Draw("SAME L")

    def legend(self, canvas, hist_list, graph_list):
        if not self.options["legend"]:
            return
        legend = ROOT.TLegend(0.6,0.6,0.89,0.89)
        legend.SetBorderSize(1)
        legend.SetHeader(self.config.beam)
        if len(hist_list) != len(self.options["legend"]["text"]):
            print "Failed to find all the histograms for legend(...)"
            return
        for i, text in enumerate(self.options["legend"]["text"]):
            if text == None:
                continue
            hist = hist_list[i]
            draw_option = self.options["legend"]["draw_option"][i]
            legend.AddEntry(hist, text, draw_option)
        legend.Draw()
        self.root_objects.append(legend)

    def mice_logo(self, canvas, hist_list, graph_list):
        if not self.options["mice_logo"]:
            return
        text_box = ROOT.TPaveText(0.7, 0.2, 0.9, 0.3, "NDC")
        text_box.SetFillColor(0)
        text_box.SetBorderSize(0)
        text_box.SetTextSize(0.04)
        text_box.SetTextAlign(12)
        text_box.AddText("MICE")
        text_box.Draw()
        self.root_objects.append(text_box)
        text_box = ROOT.TPaveText(0.7, 0.12, 0.9, 0.2, "NDC")
        text_box.SetFillColor(0)
        text_box.SetBorderSize(0)
        text_box.SetTextSize(0.03)
        text_box.SetTextAlign(12)
        text_box.AddText("ISIS Cycle 2017/02")
        text_box.AddText("MAUS v3.0")
        text_box.Draw()
        self.root_objects.append(text_box)

    def hist_title(self, canvas, hist_list, graph_list):
        if self.options["hist_title"] == None:
            return
        for hist in hist_list:
            hist.SetTitle(self.options["hist_title"])

    def log_y(self, canvas, hist_list, graph_list):
        if self.options["log_y"]:
            canvas.SetLogy()

    def write_plots(self, canvas, hist_list, graph_list):
        if not self.options["write_plots"]:
            return
        formats = self.options["write_plots"]["formats"]
        plot_dir = self.options["write_plots"]["dir"]
        canvas.Update()
        name = str(canvas.GetName()).replace(" ", "_")
        name = name.replace(".",  "_")
        name = name.split("-")
        name = "-".join(name[:-1])
        name = os.path.join(plot_dir, name)
        for fmt in formats:
            canvas.Print(name+"."+fmt)

    def murgle_histograms(self, canvas, hist_list, graph_list):
        self.rebin(canvas, hist_list, graph_list)
        self.normalise(canvas, hist_list, graph_list)
        self.defit(canvas, hist_list, graph_list)
        self.rescale_x(canvas, hist_list, graph_list)
        self.rescale_y(canvas, hist_list, graph_list)
        self.calculate_errors(canvas, hist_list, graph_list)
        self.hist_title(canvas, hist_list, graph_list)
        self.log_y(canvas, hist_list, graph_list)
        self.redraw(canvas, hist_list, graph_list)
        self.mice_logo(canvas, hist_list, graph_list)
        self.legend(canvas, hist_list, graph_list)
        self.write_plots(canvas, hist_list, graph_list)

    def conglomerate(self):
        file_list = []
        for a_path in self.dir_path:
            path = a_path+"/"+self.options["file_name"]+".root"
            file_list += sorted(glob.glob(path))
        print file_list

        hist_list = []
        graph_list = []
        for file_name in file_list:
            old_canvas = self.get_canvas(file_name)
            hist_list += self.get_hist_list(old_canvas)
            graph_list += self.get_graph_list(old_canvas)
        if len(hist_list) == 0:
            print "Error - failed to find plots for", self.options["file_name"]
        print "Found", len(hist_list), "histograms"
        for hist in hist_list:
            print "   ", hist.GetName()

        new_canvas = xboa.common.make_root_canvas(self.options["canvas_name"])
        new_canvas.cd()
        try:
            self.murgle_histograms(new_canvas, hist_list, graph_list)
        except Exception:
            sys.excepthook(*sys.exc_info())
            print "Failed with dir path", self.dir_path
            print "       and file list", file_list

    dir_path = []
    root_objects = []

class MergeCutsSummaryTex(object):
    def __init__(self):
        self.summary_list = []
        self.headings = []
        self.data = []
        self.caption = ""

    def append_summary(self, config, dir_selection):
        for i, a_dir in enumerate(config.conglomerate_dir):
            if i not in dir_selection:
                continue
            self.summary_list.append(a_dir+config.cuts_tex)
    
    def parse_one_line(self, line):
        words_tmp = line.split("&")
        words = []
        for word in words_tmp:
            #word = word.replace(" ", "")
            word = word.replace("//", "")
            word = word.replace("\n", "")
            words.append(word)
        return words

    def update_headings_list(self, line_number, words):
        if line_number >= len(self.headings):
            self.headings.append(words[0])
        elif self.headings[line_number] != words[0]:
            raise ValueError("could not match heading "+str(self.headings[line_number])+" to input "+str(words[0]))

    def update_data(self, file_number, line_number, words):
        row = line_number
        column = file_number

        if row == len(self.data):
            self.data.append([])
        if len(words) < 2:
            self.data[row].append(None)
        else:
            self.data[row].append(words[1])

    def print_data(self, folder, file_name):
        try:
            shutil.rmtree(folder)
        except OSError:
            pass
        os.makedirs(folder)
        fout = open(os.path.join(folder, file_name), "w")
        n_rows = len(self.data[0])
        head_matter = """
\\documentclass{letter}
\\usepackage{pdflscape}
\\begin{document}
"""
        tail_matter = """
\\end{document}
"""

        print >> fout, """
\\newcommand{\splitcell}[2][c]{%
    \\begin{tabular}[#1]{@{}c@{}}#2\end{tabular}}

\\begin{landscape}
\\begin{table}
\\centering
\\caption{"""+self.caption+"""}
\\begin{tabular}[pos]{l|""",
        for i in range(len(self.data[0])):
            fout.write("c")
        print >> fout, "}"
        for row in range(len(self.data)):
            #if row == 0:
            #    continue
            a_head = self.headings[row].replace("_", " ")
            print >> fout, a_head,
            for column in range(len(self.data[row])):
                if self.data[row][column] == None:
                    continue
                item = self.data[row][column]
                if row == 0:
                    item = "\\splitcell{"+item.replace(" ", "\\\\")+"}"
                print >> fout, "&", item,
            if self.headings[row] != "\hline":
                print >> fout, "\\\\",
            print >> fout
        print >> fout, """
\\end{tabular}
\\end{table}
\\end{landscape}
"""


    def latex(self, folder, file_name):
        return
        here = os.getcwd()
        os.chdir(folder)
        subprocess.check_output(["pdflatex", file_name])
        os.chdir(here)

    def merge_summaries(self, folder, file_name):
        for file_number, summary in enumerate(self.summary_list):
            try:
                fin = open(summary)
            except IOError:
                print "Failed to open", summary
                continue
            for line_number, line in enumerate(fin.readlines()):
                words = self.parse_one_line(line)
                if len(words) == 0:
                    continue
                self.update_headings_list(line_number, words)
                self.update_data(file_number, line_number, words)
        self.print_data(folder, file_name)
        self.latex(folder, file_name)

def main_cuts():
    dir_list = [
      "2017-2.7_3-140_lH2_empty", "2017-2.7_3-140_lH2_full",
      "2017-2.7_3-140_LiH", "2017-2.7_3-140_None",
      "2017-2.7_6-140_lH2_empty", "2017-2.7_6-140_lH2_full",
      "2017-2.7_6-140_LiH", "2017-2.7_6-140_None",
      "2017-2.7_10-140_lH2_empty", "2017-2.7_10-140_lH2_full",
      "2017-2.7_10-140_LiH", "2017-2.7_10-140_None",
    ]
    dir_list = ["2017-2.7_10-140_LiH"]
    data_cuts_summary = MergeCutsSummaryTex()
    mc_cuts_summary = MergeCutsSummaryTex()
    for beam in dir_list:
        try:
            config = CompareCutsConfig(beam)
            cong = ConglomerateBase(config)
            cong.conglomerate()
            data_cuts_summary.append_summary(config, [0])
            mc_cuts_summary.append_summary(config, [1])
        except Exception:
            sys.excepthook(*sys.exc_info())
    data_cuts_summary.caption = config.data_caption
    mc_cuts_summary.caption = config.mc_caption
    data_cuts_summary.merge_summaries("output/2017-02/cuts_summary/data/", "data_cuts_summary.tex")
    mc_cuts_summary.merge_summaries("output/2017-02/cuts_summary/mc/", "mc_cuts_summary.tex")

def main_globals():
    dummy = [
        "2017-2.7_3-140_lH2_empty", "2017-2.7_3-140_lH2_full",
        "2017-2.7_3-140_LiH", "2017-2.7_3-140_None",
        "2017-2.7_6-140_lH2_empty", "2017-2.7_6-140_lH2_full",
        "2017-2.7_6-140_LiH", "2017-2.7_6-140_None",
        "2017-2.7_10-140_lH2_empty", "2017-2.7_10-140_lH2_full",
        "2017-2.7_10-140_LiH", "2017-2.7_10-140_None",
        ]
    dir_list = ["2017-2.7_10-140_LiH"]
    for beam in dir_list:
        try:
            config = CompareGlobalsConfig(beam)
            cong = ConglomerateBase(config)
            cong.conglomerate()
        except Exception:
            sys.excepthook(*sys.exc_info())

if __name__ == "__main__":
    main_cuts()
    #main_globals()
    raw_input("Finished - press <CR> to end")
