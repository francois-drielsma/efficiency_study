import math
import numpy
import os

import xboa.common
import ROOT

from array import array

import utilities.root_style

class PlotAmplitudeData(object):
    def __init__(self, amplitude_data, plot_dir, key):
        self.data = amplitude_data
        self.plot_dir = plot_dir
        self.key = key

    def plot(self):
        self.plot_data_1d("emittance_vs_n_events_"+self.key, self.emittance_4d_lambda, "#varepsilon_{4D} [mm]", self.n_events_lambda, "Number of Events")
        self.plot_data_1d("emittance_vs_beta_x_"+self.key, self.emittance_4d_lambda, "#varepsilon_{4D} [mm]", self.beta_x_lambda, "#beta_{x} [mm]")
        self.plot_data_1d("emittance_vs_beta_y_"+self.key, self.emittance_4d_lambda, "#varepsilon_{4D} [mm]", self.beta_y_lambda, "#beta_{y} [mm]")
        for i in range(4):
            for j in range(i+1, 4):
                self.plot_phase_space(i, j)
		self.plot_phase_space_scatter(i, j)

    def plot_data_1d(self, plot_name, plot_lambda_x, x_label, plot_lambda_y, y_label):
        if len(self.data.state_list) == 0:
            print "Warning - no data for emittance vs beta plots/etc"
            return
        marker = 24
        x_axis = []
        y_axis = []
        for sample_states in self.data.state_list:
            for state in sample_states:
                x_axis.append(plot_lambda_x(state))
                y_axis.append(plot_lambda_y(state))
            hist, graph = xboa.common.make_root_graph(plot_name, x_axis, x_label, y_axis, y_label)
            canvas = xboa.common.make_root_canvas("plot")
            canvas.Draw()
            hist.Draw()
        for sample_states in self.data.state_list:
            x_axis = []
            y_axis = []
            for state in sample_states:
                x_axis.append(plot_lambda_x(state))
                y_axis.append(plot_lambda_y(state))
            hist, graph = xboa.common.make_root_graph(plot_name, x_axis, x_label, y_axis, y_label)
            graph.SetMarkerStyle(marker)
            marker += 1
            graph.Draw("SAME P")
        canvas.Update()
        for fmt in ["png", "pdf", "root"]:
            canvas.Print(self.plot_dir+"/phase_space/"+plot_name+"."+fmt)

    def plot_phase_space(self, x_var, y_var):
        x_label = self.psv_labels[x_var]
        y_label = self.psv_labels[y_var]
        title = "amplitude_phase_space_"+self.key+"_"+self.psv_names[x_var]+"_"+self.psv_names[y_var]
        canvas = xboa.common.make_root_canvas(title)
        x_data, y_data = [], []
        for sample in range(2):
            for a_bin in range(21):
                for run, spill, evt, psv, amp in self.data.retrieve(a_bin, sample):
                    x_data.append(psv[x_var])
                    y_data.append(psv[y_var])

        hist = xboa.common.make_root_histogram(title, x_data, x_label, 50, y_data, y_label, 50)
        canvas.SetFrameFillColor(utilities.root_style.get_frame_fill())
        hist.Draw("COLZ")
        delta_list = [-10, -7, -4, 1, 3]
        for color, sample in [(ROOT.kGreen, 0), (ROOT.kRed, 1)]: 
            step = len(self.data.state_list[sample])/len(delta_list)+1
            for i, ellipse in enumerate(self.data.state_list[sample][::step]):
                my_color = color+4
                if i < len(delta_list):
                    my_color = color+delta_list[i]
                graph = self.plot_ellipse(ellipse, x_var, y_var)
                graph.SetLineColor(my_color)
        canvas.Update()
        for fmt in ["root", "png", "pdf"]:
            canvas.Print(self.plot_dir+"/phase_space/"+title+"."+fmt)

    def plot_phase_space_scatter(self, x_var, y_var):
	# Initialize the canvas
        x_label = self.psv_labels[x_var]
        y_label = self.psv_labels[y_var]
        title = "amplitude_phase_space_scatter_"+\
		self.key+"_"+self.psv_names[x_var]+"_"+self.psv_names[y_var]
        canvas = xboa.common.make_root_canvas(title)

	# Fill the data
        data = []
        for sample in range(2):
            for a_bin in range(21):
                for run, spill, evt, psv, amp in self.data.retrieve(a_bin, sample):
                    data.append([psv[x_var], psv[y_var], amp])

	# Sort by descending order of amplitude (low amplitude points drawn last)
	data.sort(key=lambda el: el[2], reverse=True)

	# Initalize and draw the graph
	x_data = [el[0] for el in data]
	y_data = [el[1] for el in data]
	amps = [el[2] for el in data]
        graph = ROOT.TGraph2D(len(x_data), array('d', x_data), array('d', y_data), array('d', amps))
	graph.SetTitle(";;;A_{#perp}  [mm]")
        graph.Draw("PCOLZ")

	# Set the view right (from above), set the style
	ROOT.gPad.SetTheta(90);
	ROOT.gPad.SetPhi(0);

	graph.SetMarkerStyle(20);
	graph.SetMarkerSize(1);

	# Remove default x and y axes, replace them by correctly placed TGaxis
	graph.GetXaxis().SetTitleOffset(999);
	graph.GetXaxis().SetLabelOffset(999);
	graph.GetXaxis().SetTickSize(0);

	graph.GetYaxis().SetTitleOffset(999);
	graph.GetYaxis().SetLabelOffset(999);
	graph.GetYaxis().SetTickSize(0);

	wmin = graph.GetXaxis().GetXmin();
	wmax = graph.GetXaxis().GetXmax();
	xaxis = ROOT.TGaxis(-.5775, -.5775, .5775, -.5775, wmin, wmax, 505);
        xaxis.SetTitle(x_label)
	xaxis.SetLabelFont(42);
	xaxis.SetTextFont(42);
	xaxis.Draw("SAME")

	wmin = graph.GetYaxis().GetXmin();
	wmax = graph.GetYaxis().GetXmax();
	yaxis = ROOT.TGaxis(-.5775, -.5775, -.5775, .5775, wmin, wmax, 505);
        yaxis.SetTitle(y_label)
	yaxis.SetLabelFont(42);
	yaxis.SetTextFont(42);
	yaxis.Draw("SAME")

	# Output
        for fmt in ["root", "png", "pdf"]:
            canvas.Print(self.plot_dir+"/phase_space/"+title+"."+fmt)

    @classmethod
    def plot_ellipse(cls, ellipse, var_1, var_2):
        mean = [ellipse["mean"][var_1], ellipse["mean"][var_2]]
        cov = [[ellipse["cov"][i][j] for i in [var_1, var_2]] for j in [var_1, var_2]]
        try:
            points = xboa.common.make_shell(41, numpy.array(cov))
        except Exception:
            graph = ROOT.TGraph()
            return graph
        graph = ROOT.TGraph(len(points)+1)
        points = [(a_point[0, 0], a_point[0, 1]) for a_point in points]
        points = sorted(points, key = lambda points: math.atan2(points[1], points[0]))
        points.append(points[0])
        for i, a_point in enumerate(points):
            graph.SetPoint(i, a_point[0]+mean[0], a_point[1]+mean[1])
        graph.SetLineWidth(2)
        graph.Draw("SAME L")
        cls.root_objects.append(graph)
        return graph


    @classmethod
    def emittance_4d_lambda(cls, state):
        return state["emittance"]

    @classmethod
    def n_events_lambda(cls, state):
        return state["n_events"]

    @classmethod
    def beta_x_lambda(cls, state):
        return cls.beta_2d(state, 0)

    @classmethod
    def beta_y_lambda(cls, state):
        return cls.beta_2d(state, 2)

    @classmethod
    def beta_2d(cls, state, axis):
        twod_matrix = [item[axis:axis+2] for item in state["cov"][axis:axis+2]]
        emit = numpy.linalg.det(twod_matrix)**0.5/cls.mu_mass
        beta = twod_matrix[0][0]/emit
        return beta

    psv_names = ["x", "px", "y", "py"]
    psv_labels = [
        "x [mm]",
        "p_{x} [MeV/c]",
        "y [mm]",
        "p_{y} [MeV/c]",
    ]

    mu_mass = xboa.common.pdg_pid_to_mass[13]
    root_objects = []
