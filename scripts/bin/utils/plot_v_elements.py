import json
import os

import ROOT
import xboa.common


class PlotVElements(object):
    def __init__(self, file_name):
        self.file_name = file_name
        self.root_objects = []
        self.load_file()
        self.do_plot(True)
        self.do_plot(False)

    def load_file(self):
        fin = open(self.file_name)
        self.data = json.loads(fin.read())
        self.plot_dir = self.file_name
        for i in range(2):
            self.plot_dir = os.path.split(self.plot_dir)[0]
        print "Plotting to", self.plot_dir

    def do_plot(self, normalised):
        self.canvas = ROOT.TCanvas("v_elements", "v_elements", 1000, 1000)
        self.canvas.Divide(4, 4, 0., 0.)
        for i in range(4):
            for j in range(i, 4):
                canvas_index = i*4+j+1
                self.canvas.cd(canvas_index)
                self.plot_one(i, j, True)
        self.canvas.Update()
        norm = "raw"
        if normalised:
            norm = "norm"
        self.canvas.cd(13)
        self.plot_n()
        self.canvas.Print(self.plot_dir+"/covariance_canvas_"+norm+".png")

    def plot_one(self, i, j, normalised):
        multigraph = ROOT.TMultiGraph()
        for series_i, series in enumerate(self.data):
            max_amp_list = []
            var_list = []
            for amp_bin in series:
                if normalised and i == j:
                    norm = amp_bin["emittance"]
                else:
                    norm = 1.
                max_amp_list.append(amp_bin["max_amplitude"])
                var = amp_bin["cov"][i][j]/norm
                if i != j:
                    var /= amp_bin["cov"][i][i]**0.5
                    var /= amp_bin["cov"][j][j]**0.5
                var_list.append(var)
            nbins = len(max_amp_list)
            graph = ROOT.TGraph(nbins)
            graph.SetMarkerStyle(series_i+24)
            for bin_i in range(nbins):
                graph.SetPoint(bin_i, max_amp_list[bin_i], var_list[bin_i])
            multigraph.Add(graph)
        multigraph.Draw("A P")
        multigraph.GetXaxis().SetLimits(0., 105.)
        xmin = multigraph.GetYaxis().GetXmin()
        xmax = multigraph.GetYaxis().GetXmax()
        if xmin > 0:
            multigraph.SetMinimum(0.)
        elif xmax < 0:
            multigraph.SetMaximum(0.)
        self.make_box(i, j, normalised)
        self.root_objects.append(multigraph)

    def plot_n(self):
        multigraph = ROOT.TMultiGraph()
        for series_i, series in enumerate(self.data):
            max_amp_list = []
            var_list = []
            for amp_bin in series:
                max_amp_list.append(amp_bin["max_amplitude"])
                var_list.append(amp_bin["n_events"])
            nbins = len(max_amp_list)
            graph = ROOT.TGraph(nbins)
            graph.SetMarkerStyle(series_i+24)
            for bin_i in range(nbins):
                graph.SetPoint(bin_i, max_amp_list[bin_i], var_list[bin_i])
            multigraph.Add(graph)
        multigraph.Draw("A P")
        multigraph.GetXaxis().SetLimits(0., 105.)
        multigraph.SetMinimum(-0.2*multigraph.GetYaxis().GetXmax())
        self.make_box(-1, 0, True)
        self.root_objects.append(multigraph)

    def make_box(self, i, j, normalised):
        index = ["x", "px", "y", "py"]
        text_box = ROOT.TPaveText(0.6, 0.01, 0.99, 0.2, "NDC")
        text_box.SetFillColor(0)
        text_box.SetBorderSize(0)
        text = "Var("+index[i]+", "+index[j]+")"
        if normalised:
            text += "/#varepsilon_{rms}"
        if i != j:
            text = "Corr("+index[i]+", "+index[j]+")"
        if i == -1:
            text = "n_events"
        text_box.AddText(text)
        text_box.Draw()
        self.root_objects.append(text_box)


if __name__ == "__main__":
    PlotVElements("output/2017-02-7-test/plots_2017-2.7_10-140_LiH/amplitude/data/amp_data_recon_us.json")