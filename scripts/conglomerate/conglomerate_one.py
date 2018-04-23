import glob
import os
import copy
import ROOT


class ConglomerateContainer(object):
    def __init__(self, config):
        self.config = config
        self.conglomerations = []
    
    def conglomerate(self):
        for item in self.config.conglomerate_list:
            cong = ConglomerateOne(item, self.config)
            cong.dir_path = [a_dir+self.config.src_plot_dir for a_dir in self.config.conglomerate_dir]
            cong.conglomerate()
            self.conglomerations.append(cong)

class ConglomerateOne(object):
    def __init__(self, options, config):
        self.options = options
        self.config = config
        self.canvas = None
        self.pad = None
        self.hist_list = []
        self.graph_list = []
        self.legend = None
        self.x_axis = None
        self.y_axis = None
        self.label_size = 0.05

    def get_pad(self, canvas):
        pad = canvas # if we cant find pad, use canvas instead
        for an_object in canvas.GetListOfPrimitives():
            name = an_object.GetName()
            if "pad" in name:
                pad = an_object
                break
        return pad

    def get_hist_list(self, canvas):
        hist_list = []
        pad = self.get_pad(canvas)
        for an_object in pad.GetListOfPrimitives():
            name = str(an_object.GetName()).replace(" ", "_")
            for hist_name in self.options["histogram_names"]:
                hist_name = hist_name.replace(" ", "_")
                if hist_name in name:
                    if type(an_object) == type(ROOT.TGraph()):
                        continue
                    hist_list.append(an_object)
                    an_object.SetName(an_object.GetName()+self.uid())
                    break # we can only add each object once
        if len(hist_list) == 0:
            print "Did not find a histogram from list", self.options["histogram_names"]
            for an_object in pad.GetListOfPrimitives():
                print an_object.GetName()
            raise KeyError("Failed to find histogram")
        return hist_list

    def get_graph_list(self, canvas):
        graph_list = []
        pad = self.get_pad(canvas)
        graph_types = [type(ROOT.TGraph()), type(ROOT.TGraphAsymmErrors())]
        graph_name_list = copy.deepcopy(self.options["graph_names"])
        graph_name_list = [name.replace(" ", "_").lower() for name in graph_name_list]
        for an_object in pad.GetListOfPrimitives():
            name = str(an_object.GetName()).replace(" ", "_").lower()
            for graph_name in graph_name_list:
                if graph_name in name:
                    if type(an_object) not in graph_types:
                        print name, "matches", graph_name, "but type", type(an_object), "not a graph"
                        continue
                    graph_list.append(an_object)
                    graph_name_list.remove(graph_name)
                    break # we can only add each object once
        if len(graph_name_list) != 0:
            print "Did not find graphs", graph_name_list, "from:"
            for an_object in pad.GetListOfPrimitives():
                print an_object.GetName()
        return graph_list

    def get_canvas(self, file_name):
        file_name = glob.glob(file_name)[0]
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
            print key.ReadObj().GetName(), "of type", type(key.ReadObj())
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
        print "RESCALED", x_range

    def replace_hist(self, canvas, hist_list, graph_list):
        if not self.options["replace_hist"]:
            return
        x_range = self.options["rescale_x"]
        y_range = self.options["rescale_y"]
        if not x_range or not y_range:
            print "Error - need to define rescale_x and rescale_y if replace_hist"
        hist = ROOT.TH2D(canvas.GetName()+" hist replacement", "", 1000, x_range[0], x_range[1], 1000, y_range[0], y_range[1])
        hist.SetStats(False)
        for axis in hist.GetXaxis(), hist.GetYaxis():
            axis.SetNdivisions(5, 5, 0)
            axis.SetLabelSize(self.label_size)
        hist_list.insert(0, hist)
        
    def rescale_y(self, canvas, hist_list, graph_list):
        if not self.options["rescale_y"]:
            return
        hist = hist_list[0]
        if self.options["rescale_y"] == True or self.options["rescale_y"] == "auto":
            min_bin, max_bin = hist.GetMinimum(), hist.GetMaximum()
            for hist in hist_list[1:]:
                min_bin = max(min_bin, hist.GetMinimum())
                max_bin = max(max_bin, hist.GetMaximum())
            if self.options["log_y"]:
                min_bin = max(1e-4, min_bin)
            for hist in hist_list:
                hist.GetYaxis().SetRangeUser(min_bin, max_bin*1.1)
        else:
            hist.GetYaxis().SetRangeUser(self.options["rescale_y"][0], self.options["rescale_y"][1])

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

    def axis_title(self, canvas, hist_list, graph_list):
        if not self.options["axis_title"]:
            return
        for hist in hist_list:
            hist.GetXaxis().SetTitle("")
            hist.GetYaxis().SetTitle("")
        canvas.cd()
        if self.options["axis_title"]["x"] != None:
            x_text_box = ROOT.TPaveText(0.10, 0.02, 0.97, 0.04, "NDC")
            x_text_box.SetFillColor(0)
            x_text_box.SetBorderSize(0)
            x_text_box.SetTextSize(0.05)
            x_text_box.SetTextAlign(22)
            x_text_box.AddText(self.options["axis_title"]["x"])
            x_text_box.Draw()
            self.x_axis = x_text_box
        if self.options["axis_title"]["y"] != None:
            y_text_box = ROOT.TPaveText(0.015, 0.10, 0.04, 0.97, "NDC")
            y_text_box = ROOT.TLatex(0.05, 0.5, self.options["axis_title"]["y"])
            y_text_box.SetTextSize(0.05)
            y_text_box.SetTextAlign(22)
            y_text_box.SetTextAngle(90)
            y_text_box.Draw()
            self.y_axis = y_text_box

    unique_id = 0
    @classmethod
    def uid(cls):
        cls.unique_id += 1
        return str(cls.unique_id)

    def redraw(self, canvas, hist_list, graph_list):
        redraw = self.options["redraw"]
        draw_option = ["" for hist in hist_list]
        draw_order = [i for i, hist in enumerate(hist_list)]
        if redraw:
            if len(hist_list) != len(redraw["line_color"]) and not redraw["ignore_more_histograms"]:
                print "Failed to find all the histograms for redraw(...); found", len(hist_list), "expected", len(redraw["line_color"])
                raise RuntimeError("Failed to find all histograms for redraw")
            for i, hist in enumerate(hist_list):
                if type(hist) == type(ROOT.TGraph()):
                    continue
                if i >= len(redraw["line_color"]):
                    continue # ignore
                hist.SetLineColor(redraw["line_color"][i])
                hist.SetFillColor(redraw["fill_color"][i])
                hist.SetMarkerStyle(redraw["marker_style"][i])
                for axis in hist.GetXaxis(), hist.GetYaxis():
                    axis.SetNdivisions(5, 5, 0)
                    axis.SetLabelSize(self.label_size)
                hist.SetName(hist.GetName()+"_"+self.uid())
            draw_option = redraw["draw_option"]
            draw_order = redraw["draw_order"]
            if redraw["x_range"] != None:
                hist.GetXaxis().SetRangeUser(redraw["x_range"][0], redraw["x_range"][1])
            if redraw["y_range"] != None:
                hist.GetYaxis().SetRangeUser(redraw["y_range"][0], redraw["y_range"][1])
        same = ""
        for i in draw_order:
            hist_list[i].Draw(same+draw_option[i])
            same = "SAME "
        if "graph" in redraw:
            graph_draw = redraw["graph"]
        else:
            graph_draw = {
                "marker_style":None,
                "marker_color":None,
                "draw_option":["p"]*len(graph_list)
            }
        for i, graph in enumerate(graph_list):
            graph.SetName(graph.GetName()+"_"+self.uid())
            if graph_draw["marker_style"] != None:
                graph.SetMarkerStyle(graph_draw["marker_style"][i])
            if graph_draw["marker_color"] != None:
                graph.SetMarkerColor(graph_draw["marker_color"][i])
            graph.Draw("SAME "+graph_draw["draw_option"][i])

    def do_legend(self, canvas, hist_list, graph_list):
        if not self.options["legend"]:
            return
        pos = tuple(self.options["legend"]["pos"])
        legend = ROOT.TLegend(*pos)
        legend.SetBorderSize(1)
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
        self.legend = legend

    def mice_logo(self, canvas, hist_list, graph_list):
        if not self.options["mice_logo"]:
            return
        text_box = ROOT.TPaveText(0.55, 0.3, 0.9, 0.4, "NDC")
        text_box.SetFillColor(0)
        text_box.SetBorderSize(0)
        text_box.SetTextSize(0.04)
        text_box.SetTextAlign(12)
        text_box.AddText("Internal Only")
        text_box.Draw()
        self.root_objects.append(text_box)
        text_box = ROOT.TPaveText(0.55, 0.15, 0.9, 0.3, "NDC")
        text_box.SetFillColor(0)
        text_box.SetBorderSize(0)
        text_box.SetTextSize(0.03)
        text_box.SetTextAlign(12)
        text_box.AddText("ISIS Cycle 2017/02 and 2017/03")
        text_box.AddText(str(self.config.channel))
        text_box.AddText(str(self.config.beamline))
        text_box.AddText(str(self.config.absorber))
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

    def write_plots(self, canvas, pad, hist_list, graph_list):
        if not self.options["write_plots"]:
            return
        formats = self.options["write_plots"]["formats"]
        plot_dir = self.options["write_plots"]["dir"]
        canvas.Update()
        canvas.cd()
        canvas.Draw()
        canvas.Update()
        name = str(canvas.GetName()).replace(" ", "_")
        name = name.replace(".",  "_")
        name = os.path.join(plot_dir, name)
        for fmt in formats:
            canvas.Print(name+"."+fmt)

    def murgle_many(self, canvas, hist_list, graph_list):
        self.label_size = 0.1
        #self.rebin(canvas, hist_list, graph_list)
        #self.normalise(canvas, hist_list, graph_list)
        #self.defit(canvas, hist_list, graph_list)
        #self.rescale_x(canvas, hist_list, graph_list)
        #self.rescale_y(canvas, hist_list, graph_list)
        #self.calculate_errors(canvas, hist_list, graph_list)
        #self.hist_title(canvas, hist_list, graph_list)
        #self.log_y(canvas, hist_list, graph_list)
        self.redraw(canvas, hist_list, graph_list)
        #self.mice_logo(canvas, hist_list, graph_list)
        self.do_legend(canvas, hist_list, graph_list)
        #self.axis_title(canvas, hist_list, graph_list)

    def murgle_histograms(self, canvas, hist_list, graph_list):
        self.label_size = 0.05
        self.rebin(canvas, hist_list, graph_list)
        self.normalise(canvas, hist_list, graph_list)
        self.defit(canvas, hist_list, graph_list)
        self.replace_hist(canvas, hist_list, graph_list)
        self.rescale_x(canvas, hist_list, graph_list)
        self.rescale_y(canvas, hist_list, graph_list)
        self.calculate_errors(canvas, hist_list, graph_list)
        self.hist_title(canvas, hist_list, graph_list)
        self.log_y(canvas, hist_list, graph_list)
        self.redraw(canvas, hist_list, graph_list)
        self.mice_logo(canvas, hist_list, graph_list)
        self.do_legend(canvas, hist_list, graph_list)
        self.axis_title(canvas, hist_list, graph_list)

    def conglomerate(self):
        file_list = []
        for a_path in self.dir_path:
            path = a_path+"/"+self.options["file_name"]+".root"
            file_list += sorted(glob.glob(path))
        print "Searching", a_path, "for", self.options["file_name"], "found", len(file_list), "files"

        for file_name in file_list:
            old_canvas = self.get_canvas(file_name)
            self.hist_list += self.get_hist_list(old_canvas)
            self.graph_list += self.get_graph_list(old_canvas)
        if len(self.hist_list) == 0:
            print "Error - failed to find plots for", self.options["file_name"]
        #print "Found", len(hist_list), "histograms"
        #for hist in hist_list:
        #    print "   ", hist.GetName()
        self.canvas = ROOT.TCanvas(self.options["canvas_name"]+"_"+self.uid(), self.options["canvas_name"], 1400, 1000)
        self.canvas.Draw()
        self.root_objects.append(self.canvas)
        self.pad = ROOT.TPad(self.options["file_name"]+"-pad", "pad", 0.10, 0.05, 0.97, 1.0)
        self.pad.Draw()
        self.pad.cd()
        try:
            self.murgle_histograms(self.canvas, self.hist_list, self.graph_list)
        except Exception:
            print "Failed with dir path", self.dir_path
            print "       and file list", file_list
            raise
        self.write_plots(self.canvas, self.pad, self.hist_list, self.graph_list)

    dir_path = []
    root_objects = []
