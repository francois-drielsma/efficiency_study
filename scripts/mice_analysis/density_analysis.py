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

    def death(self):
        """
        Called when all the data has been loaded
        """
	# First, build the profiles, store the corresponding arrays
	# If requested, draw Poincare sections of phase space in density/phase_space
	self.make_profiles("reco")
        if self.config_anal["density_mc"]:
            self.make_profiles("all_mc")
            self.make_profiles("reco_mc")

	# Evaluate corrections if necessary
	# If requested, draw the corrections in density/corrections
        if self.calculate_corrections:
	    self.set_corrections()

	# Apply the corrections, evaluate the systematic uncertainties (TODO)
	self.corrections_and_uncertainties("reco")
        if self.config_anal["density_mc"]:
            self.corrections_and_uncertainties("all_mc")
            self.corrections_and_uncertainties("reco_mc")

	# Draw the density profiles
	self.draw_profiles()

	# Save everything to a json file
        self.save()

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

    def make_profiles(self, typ):
        """
        Produces density profiles. Extract a numpy ndarray that contains all of 
	the phase space vectors, initialize the kNN density estimator and extract
	the density profile with its uncertainty
        * typ specifies the type of data (all_mc, reco_mc, reco)
        """
	# Skip if no data
	if not self.data[typ]["us"].get_n_events():
	    return

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
	    levels, errors = density_estimator.profile(self.npoints, self.uncertainty)

	    # Store the profiles and their statistical uncertainties
	    self.density_data[typ][loc]["levels"] = levels
	    self.density_data[typ][loc]["levels_stat_errors"] = errors

	    # Plot the Poincare sections if requested
	    if self.config_anal["density_sections"]:
		plotter = DensityPlotter(self.plot_dir, typ+"_"+loc)
		plotter.plot_phase_space(density_estimator)

    def corrections_and_uncertainties(self, typ):
        """
        Calculate corrected profiles and uncertainties
        * typ specifies the type of data (all_mc, reco_mc, reco)
        * corrected profiles is given by rho(alpha) = c(alpha)*rho(alpha) where c(alpha)
            is the inverse reconstruction response function
        * statistical errors are provided by the density estimator
        * total errors are given by sum in quadrature of statistical errors and
          systematic errors
        """
	# Set the upstream statistical uncertainty to 0 as it is the given
	# profile to which to compare the downstream profile
        data = self.density_data[typ]
        data["us"]["levels_stat_errors"] = [0. for bin in range(self.npoints)]

	# Do the correction for each of the tracker locations
        for loc in self.locations:
            print "Doing density level correction for", typ, loc
            levels = self.do_corrections(typ, loc)
            data[loc]["corrected_levels"] = levels

	# Evaluate the systematic uncertainties for each of the tracker locations
        for loc in self.locations:
	    continue # TODO
            print "Finding systematic errors for", typ, loc
            reco_sys_list = self.calculate_detector_systematics(typ, loc)
            perf_sys_list = self.calculate_performance_systematics(typ, loc)
            sys_error_list = [(reco_sys_list[i]**2+perf_sys_list[i]**2)**0.5 \
                                                  for i in range(len(pdf_list))]
            data[loc]["pdf_sys_errors"] = sys_error_list
            print "    sys errors:   ", data[loc]["pdf_sys_errors"] 
            print "    stats errors: ", data[loc]["pdf_stats_errors"]
            print "    pdf:          ", pdf_list

        self.density_data[loc] = data

    def do_corrections(self, typ, loc, use_capped = True):
        """
        Applies the corrections to the requested density profile
	Only applies response correction to the reconstructed sample
        * typ specifies the type of data (all_mc, reco_mc, reco)
	* loc specifies the location of the tracker (us, ds)
        """
	levels = np.array(self.density_data[typ][loc]["levels"])
	corr_key = "level_ratio"
        if use_capped:
            corr_key = "level_ratio_capped"
        if typ == "reco":
            response = np.array(self.density_data["response"][loc][corr_key])
            levels = levels*response
        if typ == "reco" or typ == "reco_mc":
            inefficiency = np.array(self.density_data["inefficiency"][loc][corr_key])
            levels = levels*inefficiency

        return levels.tolist()

    def draw_profiles(self):
        """
        Produce plots that compare the density profiles upstream and downstream
	of the absorber. Produce one for each category of data
        """
        for typ in self.data_types:
	    # Skip if no data
	    if not self.data[typ]["us"].get_n_events():
		continue

	    # Initialize the graphs
	    graphs = {}
	    for loc in self.locations:
	        graphs[loc] = self.make_graph(typ, loc, True)

	    # Print up/down comparison
            canvas_name = 'density_profile_%s' % typ
            canvas = self.get_plot(canvas_name)["pad"]

	    mg = self.make_multigraph(typ, graphs, {"us":1,"ds":4}, "LE3", "density_profile")
	    leg = self.make_multigraph_legend(graphs)
            for fmt in ["pdf", "png", "root"]:
                canvas.Print(self.plot_dir+"/"+canvas_name+"."+fmt)

	    # Print ratios
            canvas_name = 'density_ratio_%s' % typ
            canvas = self.get_plot(canvas_name)["pad"]

	    gratio = self.make_ratio(typ, graphs, 1, "LE3", "density_ratio")
            for fmt in ["pdf", "png", "root"]:
                canvas.Print(self.plot_dir+"/"+canvas_name+"."+fmt)

    def make_graph(self, typ, loc, corr):
        """
        Builds a TGraphErrors for the requested data type and location
        * typ specifies the type of data (all_mc, reco_mc, reco)
	* loc specifies the location of the tracker (us, ds)
        * corr is True if the corrected levels are to be represented
        """
	level_type = "levels"
	if corr:
	    level_type = "corrected_levels"

	graph = ROOT.TGraphErrors(self.npoints)
	for i in range(self.npoints):
	    alpha = (float(i+1.)/(self.npoints+1.))
	    graph.SetPoint(i, alpha, self.density_data[typ][loc][level_type][i])
	    graph.SetPointError(i, 0., self.density_data[typ][loc]["levels_stat_errors"][i])

	return graph

    def make_multigraph(self, typ, graphs, colors, plot_option, name):
        """
        Initializes a multigraph, draws it
        * typ specifies the type of data (all_mc, reco_mc, reco)
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

	leg = ROOT.TLegend(.6, .65, .8, .85)
	leg.AddEntry(graphs["us"], "Upstream", "LF")
	leg.AddEntry(graphs["ds"], "Downstream", "LF")
	leg.Draw("SAME")
	return leg

    def make_ratio(self, typ, graphs, color, plot_option, name):
        """
        Initializes a graph ratio, draws it
        * typ specifies the type of data (all_mc, reco_mc, reco)
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

    def save(self):
        """
        Saves the data dictionary to a json file
        """
        fout = open(self.plot_dir+"/density.json", "w")
        print >> fout, json.dumps(self.density_data, sort_keys=True, indent=4)

    def set_corrections(self):
        """
        Calculate the profile corrections, i.e. the inefficiency and response function
        * uses the Monte Carlo to generate corrections
        """
	for loc in self.locations:
	    # Initialize the data	
            all_mc_levels = self.density_data["all_mc"][loc]["levels"]
	    reco_mc_levels = self.density_data["reco_mc"][loc]["levels"]
            reco_levels = self.density_data["reco"][loc]["levels"]

	    # Calculate the corrections
            inefficiency, response = [], []
            for i in range(len(all_mc_levels)):
		# Inherent detector inefficiency
		if reco_mc_levels[i] == 0:
                    inefficiency.append(1.)
	        else:
                    inefficiency.append(float(all_mc_levels[i])/reco_mc_levels[i])

		# Detector response function
                if reco_levels[i] == 0:
                    response.append(1.)
                else:
                    response.append(float(reco_mc_levels[i])/reco_levels[i])

	    # Produce a capped version of the corrections
            cutoff = self.config_anal["density_corrections_cutoff"]
            cutoff_index = int(cutoff*(self.npoints+1.))
	    ineff_cap = all_mc_levels[cutoff_index]/reco_mc_levels[cutoff_index]
	    resp_cap = reco_mc_levels[cutoff_index]/reco_levels[cutoff_index]

            inefficiency_capped = copy.deepcopy(inefficiency)
            response_capped = copy.deepcopy(response)
            for i in range(cutoff_index, len(inefficiency_capped)):
                inefficiency_capped[i] = ineff_cap
                response_capped[i] = resp_cap

	    # Store the correction factors
            self.density_data["inefficiency"][loc] = {
                "level_ratio":inefficiency,
                "level_ratio_capped":inefficiency_capped,
            }
            self.density_data["response"][loc] = {
                "level_ratio":response,
                "level_ratio_capped":response_capped,
            }

	    # Draw the corrections, if requested
	    if self.config_anal["density_corrections_draw"]:
		plotter = DensityPlotter(self.plot_dir, loc)
		plotter.plot_corrections(inefficiency, response)

		plotter = DensityPlotter(self.plot_dir, "capped_"+loc)
		plotter.plot_corrections(inefficiency_capped, response_capped)

    def clear_density_data(self):
        """
        Initializes the dictionary that contains all the information
	used to make corrections down the line
        """
        self.density_data = {
          "inefficiency":{
              "us":{
                  "level_ratio":[1. for i in range(self.npoints)],
            	  "level_ratio_capped":[1. for i in range(self.npoints)]
              },
              "ds":{
                  "level_ratio":[1. for i in range(self.npoints)],
            	  "level_ratio_capped":[1. for i in range(self.npoints)]
              },
          },
          "response":{
              "us":{
                  "level_ratio":[1. for i in range(self.npoints)],
            	  "level_ratio_capped":[1. for i in range(self.npoints)]
              },
              "ds":{
                  "level_ratio":[1. for i in range(self.npoints)],
            	  "level_ratio_capped":[1. for i in range(self.npoints)]
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
        for typ in systematics:
            print "Loading", typ
            if typ not in self.density_data:
                self.density_data[typ] = {}
            for ref_key in ["detector_reference", "performance_reference"]:
                ref_src = systematics[typ][ref_key]
                if ref_src == None:
                    self.density_data[typ][ref_key] = None
                else:
                    self.density_data[typ][ref_key] = \
                                              self.load_one_error(ref_src, None)
                print "  Loaded reference", typ, ref_key, ref_src, \
                                          type(self.density_data[typ][ref_key])
            for loc in ["us", "ds"]:
                if loc not in self.density_data[typ]:
                    self.density_data[typ][loc] = {}
                for key in ["detector_systematics", "performance_systematics"]:
                    err_src_dict = systematics[typ][loc][key]
                    self.density_data[typ][loc][key] = [
                        self.load_one_error(err_src, scale) \
                              for err_src, scale in err_src_dict.iteritems()
                    ]
                    print "  Loaded", len(self.density_data[typ][loc][key]), loc, key, "systematics"

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
