import tempfile
import json
import ROOT
import copy
import os
import numpy as np

from mice_analysis.density.density_data import DensityData
from mice_analysis.density.density_plotter import DensityPlotter
from mice_analysis.density.knn_density_estimator import kNNDensityEstimator
from mice_analysis.analysis_base import AnalysisBase

class DensityAnalysis(AnalysisBase):

    def __init__(self, config, config_anal, data_loader):
        """
        Initialise the DensityAnalysis class for the nonparametric density estimation analysis
        * config and config_anal are the configuration files
	* data_loader extracts the data from the MAUS ROOT files
        """
        super(DensityAnalysis, self).__init__(config, config_anal, data_loader)

	# Initialize the configuration
        self.config = config
        self.config_anal = config_anal
        self.data_loader = data_loader

	# Initialize the density estimator parameters
        self.nthreads = self.config.density_nthreads
        self.rotate = self.config.density_knn_rotate
	self.uncertainty = self.config.density_uncertainty
	self.npoints = self.config.density_npoints

	# Initialize the data containers
        self.a_dir = tempfile.mkdtemp()
        file_name = self.a_dir+"/density_data_"
        self.data_types = ("all_mc", "reco_mc", "reco")
        self.locations = ("us", "ds")
	self.data = {}
        for typ in self.data_types:
	    self.data[typ] = {}
	    for loc in self.locations:
		self.data[typ][loc] = DensityData(file_name+"_"+typ+"_"+loc)

	# Calculate corrections if required
        self.calculate_corrections = self.config_anal["density_corrections"] == None

    def birth(self):
        """
        Sets the output directory for the density plots
	Resets all the data containers
        """
        self.set_plot_dir("density")
        for typ in self.data_types:
	    for loc in self.locations:
		self.data[typ][loc].clear()

	self.clear_density_data()
	self.load_errors()
        self.append_data()

        try:
            os.mkdir(self.plot_dir+"/phase_space")
        except OSError:
            pass
        try:
            os.mkdir(self.plot_dir+"/corrections")
        except OSError:
            pass

    def process(self):
        """
        Called on each spill
        """
        self.append_data()

    def append_data(self):
        """
        Add data to the density calculation (done at death time)
        
        If density_mc is false, we take 'tku' data from upstream_cut sample 
        and 'tkd' data from downstream_cut sample
        
        If density_mc is true, then we build a couple of additional samples
        1. all_mc is for MC truth of all events that should have been included 
           in the sample (some of which might have been missed by recon; some in
           recon sample maybe should have been excluded)
        2. reco_mc is for MC truth of all events that should have been included
           in the sample and were reconstructed
        We use "mc_true_us_cut"/"mc_true_ds_cut" for the all_mc sample and we
        use ("mc_true_us_cut" AND "upstream_cut") for the reco_mc sample
        """
	hits = {}
        for typ in self.data_types:
	    hits[typ] = {}
	    for loc in self.locations:
		hits[typ][loc] = []

        if self.config_anal["density_mc"]:
            station_us = self.config.mc_plots["mc_stations"]["tku_tp"][0]
            station_ds = self.config.mc_plots["mc_stations"]["tkd_tp"][0]

        for event in self.data_loader.events:
            if event['upstream_cut']:
                continue
            hits["reco"]["us"].append(event['tku'])
            if not event['downstream_cut']:
                hits["reco"]["ds"].append(event['tkd'])

            if self.config_anal["density_mc"]:
                hit_mc_us, hit_mc_ds = None, None
                for detector_hit in event["data"]:
                    if detector_hit["detector"] == station_us:
                        hit_mc_us = detector_hit["hit"]
                    if detector_hit["detector"] == station_ds:
                        hit_mc_ds = detector_hit["hit"]
                if not event['mc_true_us_cut']:
                    hits["all_mc"]["us"].append(hit_mc_us)# no inefficiency upstream
                    hits["reco_mc"]["us"].append(hit_mc_us)
                if not event['mc_true_ds_cut']:
                    hits["all_mc"]["ds"].append(hit_mc_ds)
                    if not event['downstream_cut']:
                        hits["reco_mc"]["ds"].append(hit_mc_ds)

        for typ in self.data_types:
	    for loc in self.locations:
		self.data[typ][loc].append_hits(hits[typ][loc])

	for loc in self.locations:
            print "Loaded %s (density):" % loc
	    for typ in self.data_types:
		print "    %s:   " % typ, len(hits[typ][loc])

        return

    def death(self):
        """
        Called when all the data has been loaded
	Feeds the data to the density estimator, extracts density profiles
	Saves the profiles as a json file (TODO)
        """
	# Initialize the plotter, use it to produce the requested elements
        self.make_profiles()

	# Save the density profiles to a json file
        self.save()

    def make_profiles(self):
        """
        Produces density profile comparisons between upstream and downstream
	for each of the categories of data
        """
	# For each type of data, extract a numpy ndarray that contains all of the
        # phase space vectors, initialize the kNN density estimator and extract the
        # density profile with its uncertainty
	graphs = {}
        for typ in self.data_types:
	    # Skip if no data
	    if not self.data[typ]["us"].get_n_events():
		continue

	    graphs[typ] = {}
	    transmission = float(self.data[typ]["ds"].get_n_events())/self.data[typ]["us"].get_n_events()
	    for loc in self.locations:
		# Extract the data as a numpy array
        	ps_data = np.ndarray((self.data[typ][loc].get_n_events(), 4))
		for i, event in enumerate(self.data[typ][loc].retrieve()):
	    	    ps_data[i] = event[-1].tolist()

		# Use the kNN module to get a density profile
		norm = 1.
		if loc == "ds":
		    norm = transmission
		density_estimator = kNNDensityEstimator(ps_data, self.rotate, self.nthreads, norm)
		graphs[typ][loc] = density_estimator.profile(self.npoints, self.uncertainty)

		# Plot the Poincare sections if requested
		if self.config_anal["density_sections"]:
		    plotter = DensityPlotter(self.plot_dir, typ+"_"+loc)
		    plotter.plot_phase_space(density_estimator)

	# Save the density profiles to a text file
	self.set_levels(graphs)

	# Compute the corrections if necessary
	self.set_corrections()

	# Produce plots that compare the density profiles upstream and downstream
        # of the absorber. Produce one for each category of data
        for typ in self.data_types:
	    # Skip if no data
	    if not self.data[typ]["us"].get_n_events():
		continue

	    # Print up/down comparison
            canvas_name = 'density_profile_%s' % typ
            canvas = self.get_plot(canvas_name)["pad"]

	    mg = self.make_multigraph(typ, graphs[typ], {"us":1,"ds":4}, "LE3", "density_profile")
	    leg = self.make_multigraph_legend(graphs[typ])
            for fmt in ["pdf", "png", "root"]:
                canvas.Print(self.plot_dir+"/"+canvas_name+"."+fmt)

	    # Print ratios
            canvas_name = 'density_ratio_%s' % typ
            canvas = self.get_plot(canvas_name)["pad"]

	    gratio = self.make_ratio(typ, graphs[typ], 1, "LE3", "density_ratio")
            for fmt in ["pdf", "png", "root"]:
                canvas.Print(self.plot_dir+"/"+canvas_name+"."+fmt)

    def make_multigraph(self, typ, graphs, colors, plot_option, name):
        """
        Initializes a multigraph, draws it
        * typ speficies the type of data (all_mc, reco_mc, reco)
	* graphs is a dictionary containing the up and downstream graphs
	* colors is a dictionary containing the up and downstream graph colors
	* plot_option is the graph drawing option
	* name of the type of graph being output
        """

	mg = ROOT.TMultiGraph(name+"_"+typ, ";Fraction #alpha;#rho_{#alpha} [mm^{-2}(MeV/c)^{-2}]");
	for loc in self.locations:
	    graphs[loc].SetLineColor(colors[loc])
            graphs[loc].SetFillColorAlpha(colors[loc], .25)
	    mg.Add(graphs[loc], plot_option)
            self.plots[name+"_"+typ]["graphs"][loc] = graphs[loc]

        mg.Draw("A")
	return mg

    def make_multigraph_legend(self, graphs):
        """
        Initializes a multigraph legend, draws it
	* graphs is a dictionary containing the up and downstream graphs
        """

	leg = ROOT.TLegend(.6, .6, .8, .8)
	leg.AddEntry(graphs["us"], "Upstream", "LF")
	leg.AddEntry(graphs["ds"], "Downstream", "LF")
	leg.Draw("SAME")
	return leg

    def make_ratio(self, typ, graphs, color, plot_option, name):
        """
        Initializes a graph ratio, draws it
        * typ speficies the type of data (all_mc, reco_mc, reco)
	* graphs is a dictionary containing the up and downstream graphs
	* color is the graph color
	* plot_option is the graph drawing option
	* name of the type of graph being output
        """
	gratio = ROOT.TGraphErrors(self.npoints)
	gratio.SetTitle(";Fraction #alpha;#rho_{#alpha}^{d} /#rho_{#alpha}^{u}")
	for i in range(self.npoints):
	    ratio = graphs["ds"].GetY()[i]/graphs["us"].GetY()[i]
	    gratio.GetX()[i] = graphs["us"].GetX()[i]
	    gratio.GetEX()[i] = graphs["us"].GetEX()[i]
	    gratio.GetY()[i] = ratio
	    gratio.GetEY()[i] = 0.
	    if graphs["ds"].GetY()[i] > 0.:
		gratio.GetEY()[i] = ratio*graphs["ds"].GetEY()[i]/graphs["ds"].GetY()[i]

        self.plots[name+"_"+typ]["graphs"]["ratio"] = gratio

	gratio.SetLineColor(color)
	gratio.SetFillColorAlpha(color, .25)

	gratio.Draw("A"+plot_option)
	return gratio

    def performance_sys_error(self):
        pass

    def reco_sys_error(self):
        pass

    def save(self):
        fout = open(self.plot_dir+"/density.json", "w")
        print >> fout, json.dumps(self.density_data, sort_keys=True, indent=4)

    def set_levels(self, graphs):
        """
        Sets the levels in the dictionary to be output
        """
	self.density_data["npoints"] = self.npoints
        for typ in self.data_types:
	    for loc in self.locations:
		self.density_data[typ][loc]["levels"] = [1. for i in range(self.npoints)]
		for i in range(self.npoints):
		    self.density_data[typ][loc]["levels"][i] = graphs[typ][loc].GetY()[i]

    def set_corrections(self):
        """
        Calculate the profile corrections
        
        Uses the MC to generate corrections
        """
	for loc in self.locations:
            all_mc_levels = self.density_data["all_mc"][loc]["levels"]
            reco_levels = self.density_data["reco"][loc]["levels"]
            corrections = []
            for i in range(len(all_mc_levels)):
                if reco_levels[i] == 0:
                    corrections.append(1.)
                else:
                    corrections.append(float(all_mc_levels[i])/reco_levels[i])

            cutoff = self.config_anal["density_corrections_cutoff"]
            cutoff_index = int(cutoff*(self.npoints+1.))
            all_mc_sum = sum(all_mc_levels[cutoff_index:])
            reco_sum = sum(reco_levels[cutoff_index:])
            corrections_averaged = copy.deepcopy(corrections)
            for i in range(cutoff_index, len(corrections_averaged)):
                corrections_averaged[i] = all_mc_sum/reco_sum

            self.density_data["correction"][loc] = {
                "level_ratio":corrections,
                "level_ratio_averaged":corrections_averaged,
            }

	    if self.config_anal["density_corrections_draw"]:
		plotter = DensityPlotter(self.plot_dir, loc)
		plotter.plot_corrections(corrections)
		    

    def clear_density_data(self):
        """
        Initializes the dictionary that contains all the information
	used to make corrections down the line
        """
        self.density_data = {
          "correction":{
              "us":{
                  "level_ratio":[1. for i in range(self.npoints)],
            	  "level_ratio_averaged":1.
              },
              "ds":{
                  "level_ratio":[1. for i in range(self.npoints)],
            	  "level_ratio_averaged":1.
              },
          },
          "source":"",
          "scale":1.,
	  "npoints":0
        }
        level_data = {
            "performance_reference":None,
            "detector_reference":None,
            "us":{
            	"levels":[],
                "detector_systematics":[],
                "performance_systematics":[],
                "detector_systematics_output":{},
                "performance_systematics_output":{}
            },
            "ds":{
            	"levels":[],
                "detector_systematics":[],
                "performance_systematics":[],
                "detector_systematics_output":{},
                "performance_systematics_output":{}
            }
        }

        for typ in self.data_types:
            self.density_data[typ] = copy.deepcopy(level_data)

    def load_errors(self):
        """
        Two "classes" of systematic errors;
        * systematic errors on the reconstruction are contained in the
          correction factors. For these we store the correction factors and 
          compare to the reference correction factors
        * systematic errors on the performance are contained in the actual
          density profile. For these we store the point-by-point fractional
          difference between the density profile and reference.
        """
        if self.calculate_corrections:
            return

        # set base correction factors
        self.load_corrections(self.config_anal["density_corrections"])
        systematics = self.config_anal["density_systematics"]
        for suffix in systematics:
            print "Loading", suffix
            if suffix not in self.density_data:
                self.density_data[suffix] = {}
            for ref_key in ["detector_reference", "performance_reference"]:
                ref_src = systematics[suffix][ref_key]
                if ref_src == None:
                    self.density_data[suffix][ref_key] = None
                else:
                    self.density_data[suffix][ref_key] = \
                                              self.load_one_error(ref_src, None)
                print "  Loaded reference", suffix, ref_key, ref_src, \
                                          type(self.density_data[suffix][ref_key])
            for loc in ["us", "ds"]:
                if loc not in self.density_data[suffix]:
                    self.density_data[suffix][loc] = {}
                for key in ["detector_systematics", "performance_systematics"]:
                    err_src_dict = systematics[suffix][loc][key]
                    self.density_data[suffix][loc][key] = [
                        self.load_one_error(err_src, scale) \
                              for err_src, scale in err_src_dict.iteritems()
                    ]
                    print "  Loaded", len(self.density_data[suffix][loc][key]), loc, key, "systematics"

    def load_one_error(self, file_name, scale):
        """
        Load the density analysis output for a given uncertainty source
        """
        fin = open(file_name)
        density_str = fin.read()
        density = json.loads(density_str)
        density["source"] = file_name
        density["scale"] = scale
        return density

    def load_corrections(self, file_name):
        """
        Load the density corrections to be applied during this density
        analysis. Loads the correction factors
        """
        fin = open(file_name)
        density_str = fin.read()
        src_density = json.loads(density_str)
        src_density["source"] = file_name
        self.density_data["correction"] = src_density["correction"]
