import tempfile
import json
import ROOT
import numpy as np

from mice_analysis.density_estimator.density_data import DensityData
from mice_analysis.density_estimator.knn_density_estimator import kNNDensityEstimator
from mice_analysis.analysis_base import AnalysisBase

class DensityAnalysis(AnalysisBase):

    def __init__(self, config, config_anal, data_loader):
        """
        Initialise the DensityAnalysis class for the density estimation
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
	self.graph_npoints = self.config.density_graph_npoints

	# Initialize the data containers
	self.json_data = [] # data for plotting/storing on disk
        self.a_dir = tempfile.mkdtemp()
        file_name = self.a_dir+"/density_data_"
        self.data_types = ("mc", "reco_mc", "recon")
        self.locations = ("us", "ds")
	self.data = {}
        for typ in self.data_types:
	    self.data[typ] = {}
	    for loc in self.locations:
		self.data[typ][loc] = DensityData(file_name+"_"+typ+"_"+loc)

    def birth(self):
        """
        Sets the output directory for the density plots
	Resets all the data containers
        """
        self.set_plot_dir("density")
        for typ in self.data_types:
	    for loc in self.locations:
		self.data[typ][loc].clear()
        self.append_data()

    def process(self):
        """
        Called on each spill
        """
        self.append_data()

    def append_data(self):
        """
        Add data to the density calculation (done at death time)
        
        If amplitude_mc is false, we take 'tku' data from upstream_cut sample 
        and 'tkd' data from downstream_cut sample
        
        If amplitude_mc is true, then we build a couple of additional samples
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

        if self.config_anal["amplitude_mc"]:
            station_us = self.config.mc_plots["mc_stations"]["tku_tp"][0]
            station_ds = self.config.mc_plots["mc_stations"]["tkd_tp"][0]

        for event in self.data_loader.events:
            if event['upstream_cut']:
                continue
            hits["recon"]["us"].append(event['tku'])
            if not event['downstream_cut']:
                hits["recon"]["ds"].append(event['tkd'])

            if self.config_anal["amplitude_mc"]:
                hit_mc_us, hit_mc_ds = None, None
                for detector_hit in event["data"]:
                    if detector_hit["detector"] == station_us:
                        hit_mc_us = detector_hit["hit"]
                    if detector_hit["detector"] == station_ds:
                        hit_mc_ds = detector_hit["hit"]
                if not event['mc_true_us_cut']:
                    hits["mc"]["us"].append(hit_mc_us)# no inefficiency upstream
                    hits["reco_mc"]["us"].append(hit_mc_us)
                if not event['mc_true_ds_cut']:
                    hits["mc"]["ds"].append(hit_mc_ds)
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

    def make_multigraph(self, typ, graphs, colors, plot_option, name):
        """
        Initializes a multigraph, draws it
        * typ speficies the type of data (mc, reco_mc, recon)
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
        * typ speficies the type of data (mc, reco_mc, recon)
	* graphs is a dictionary containing the up and downstream graphs
	* color is the graph color
	* plot_option is the graph drawing option
	* name of the type of graph being output
        """
	gratio = ROOT.TGraphErrors(self.graph_npoints)
	gratio.SetTitle(";Fraction #alpha;#rho_{#alpha}^{d} /#rho_{#alpha}^{u}")
	for i in range(self.graph_npoints):
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
		graphs[typ][loc] = density_estimator.profile(self.graph_npoints, self.uncertainty)

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
		

    def performance_sys_error(self):
        pass

    def reco_sys_error(self):
        pass

    def save(self):
        for item in self.json_data:
            del item['values']
        fout = open(self.plot_dir+"/density.json", "w")
        print >> fout, json.dumps(self.json_data)

    def death(self):
        """
        Called when all the data has been loaded
	Feeds the data to the density estimator, extracts density profiles
	Saves the profiles as a json file (TODO)
        """
        self.make_profiles()
        self.save()

