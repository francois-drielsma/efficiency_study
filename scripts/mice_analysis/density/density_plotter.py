import math
import numpy
import ROOT

import xboa.common
import utilities.root_style
from mice_analysis.density.knn_density_estimator import kNNDensityEstimator

class DensityPlotter(object):

    def __init__(self, plot_dir, typ):
        """
        Initialise the DensityPlotter class for drawing density estimates
	* plot_dir is the target directory
	* typ is the type of data being represented (all_mc, reco_mc, reco)
        """
        self.plot_dir = plot_dir
	self.typ = typ

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
        """
	# Get the histogram from the density estimator
	hist = de.section([i, j], [0., 0., 0., 0.], [-99.99, -99.99], [99.99, 99.99])

	# Draw it
        x_label = self.psv_labels[i]
        y_label = self.psv_labels[j]
        title = "density_phase_space_"+self.typ+"_"+self.psv_names[i]+"_"+self.psv_names[j]
        canvas = xboa.common.make_root_canvas(title)
        canvas.SetFrameFillColor(utilities.root_style.get_frame_fill())
	hist.SetTitle(";%s;%s"%(x_label, y_label))
        hist.Draw(plot_option)

        for fmt in ["root", "png", "pdf"]:
            canvas.Print(self.plot_dir+"/phase_space/"+title+"."+fmt)

    psv_names = ["x", "px", "y", "py"]
    psv_labels = [
        "x [mm]",
        "p_{x} [MeV/c]",
        "y [mm]",
        "p_{y} [MeV/c]",
    ]
