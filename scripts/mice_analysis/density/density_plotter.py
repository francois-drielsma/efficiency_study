import math
import numpy
import ROOT

import xboa.common
import utilities.root_style
from mice_analysis.density.knn_density_estimator import kNNDensityEstimator

class DensityPlotter(object):

    def __init__(self, plot_dir, name):
        """
        Initialise the DensityPlotter class for drawing density estimates
        * plot_dir is the target directory
        * name is the name of the subset being represented
        """
        self.plot_dir = plot_dir
        self.name = name

    def plot_phase_space(self, de):
        """
        Calls the plotting subroutines
        * de is the density estimator
        """
        for i in range(4):
            for j in range(i+1, 4):
                self.plot_section(de, i, j, 21)

    def plot_section(self, de, i, j, n, plot_option="CONTZ"):
        """
        Plots a single 2D Poincare section of phase space
        * de is the density estimator
        * i is the index of the abscissae axis variable
        * j is the index of the ordinate axis variable
        * n is the number of points to evaluate the density in
        * plot_option is the drawing option
        """
        # Get the histogram from the density estimator
        hist = de.section([i, j], [0., 0., 0., 0.], [-99.99, -99.99], [99.99, 99.99])

        # Draw it
        x_label = self.psv_labels[i]
        y_label = self.psv_labels[j]
        title = "density_phase_space_"+self.name+"_"+self.psv_names[i]+"_"+self.psv_names[j]
        canvas = xboa.common.make_root_canvas(title)
        canvas.SetFrameFillColor(utilities.root_style.get_frame_fill())
        hist.SetTitle(";%s;%s"%(x_label, y_label))
        hist.Draw(plot_option)

        for fmt in ["root", "png", "pdf"]:
            canvas.Print(self.plot_dir+"/phase_space/"+title+"."+fmt)

    def plot_corrections(self, inefficiency, response):
        """
        Plots the correction factors
        * Array of corrections
        """
        # Initialize a graph for each correction category
        npoints = len(inefficiency)
        graph_ineff = ROOT.TGraph(npoints)
        graph_ineff.SetLineColor(2)
        graph_resp = ROOT.TGraph(npoints)
        graph_resp.SetLineColor(4)
        graph_corr = ROOT.TGraph(npoints)
        for i in range(npoints):
            x = float(i+1.)/(npoints+1)
            graph_ineff.SetPoint(i, x, inefficiency[i])
            graph_resp.SetPoint(i, x, response[i])
            graph_corr.SetPoint(i, x, inefficiency[i]*response[i])

        # Initialize a legend
        leg = ROOT.TLegend(.6, .65, .8, .85)
        leg.AddEntry(graph_ineff, "Inefficiency", "l")
        leg.AddEntry(graph_resp, "Response", "l")
        leg.AddEntry(graph_corr, "Combined", "l")

        # Draw the graphs, the legend, print
        title = "density_corrections_"+self.name
        canvas = xboa.common.make_root_canvas(title)
        graph_ineff.SetTitle(";Fraction #alpha;Correction factor")
        graph_ineff.Draw("AL")
        graph_ineff.SetMinimum(.9)
        graph_ineff.SetMaximum(1.3)
        graph_resp.Draw("LSAME")
        graph_corr.Draw("LSAME")
        leg.Draw("SAME")

        for fmt in ["root", "png", "pdf"]:
            canvas.Print(self.plot_dir+"/corrections/"+title+"."+fmt)


    def plot_systematics(self, graphs):
        """
        Plots the correction factors
        * Array of corrections
        """
        # Initialize a graph for each correction category, add them
        # to the multigraph and the legend
        mg = ROOT.TMultiGraph(self.name+"_systematics", ";Fraction #alpha;Relative residual");
        leg = ROOT.TLegend(.6, .65, .8, .85)
        gid = 0
        for key, graph in graphs.iteritems():
            graph.SetLineColor(2+gid)
            leg.AddEntry(graph, key, "l")
            mg.Add(graph, "l")
            gid += 1

        # Initialize envelope graphs and add them to the mutligraph
        graph_ref = graphs.itervalues().next()
        npoints = graph_ref.GetN()
        graph_upper = ROOT.TGraph(npoints)
        graph_lower = ROOT.TGraph(npoints)
        graph_upper.SetLineWidth(2)
        graph_lower.SetLineWidth(2)
        for i in range(npoints):
            quad_sum = 0.
            for graph in graphs.itervalues():
                quad_sum = (quad_sum**2 + graph.GetY()[i]**2)**0.5

            alpha = graph_ref.GetX()[i]
            graph_upper.SetPoint(i, alpha, quad_sum)
            graph_lower.SetPoint(i, alpha, -quad_sum)

        mg.Add(graph_upper, "l")
        mg.Add(graph_lower, "l")

        # Draw the graphs, the legend, print
        title = "density_systematics_"+self.name
        canvas = xboa.common.make_root_canvas(title)
        mg.Draw("A")
        mg.SetMinimum(-.2)
        mg.SetMaximum(.2)
        leg.Draw("SAME")

        for fmt in ["root", "png", "pdf"]:
            canvas.Print(self.plot_dir+"/systematics/"+title+"."+fmt)

    psv_names = ["x", "px", "y", "py"]
    psv_labels = [
        "x [mm]",
        "p_{x} [MeV/c]",
        "y [mm]",
        "p_{y} [MeV/c]",
    ]
