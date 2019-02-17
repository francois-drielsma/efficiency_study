import tempfile
import json
import ROOT
import copy
import os
import xboa.common

from mice_analysis.amplitude.amplitude_data_binned import AmplitudeDataBinned
from mice_analysis.fractional_emittance.amplitude_quantile import AmplitudeQuantile
from mice_analysis.analysis_base import AnalysisBase

class FractionalAnalysis(AnalysisBase):

    def __init__(self, config, config_anal, data_loader):
        """
        Initialise the FractionalAnalysis class for the fractional emittance analysis
        * config and config_anal are the configuration files
	* data_loader extracts the data from the MAUS ROOT files
        """
        super(FractionalAnalysis, self).__init__(config, config_anal, data_loader)

	# Initialize the configuration and the data loader
        self.config = config
        self.config_anal = config_anal
        self.data_loader = data_loader

	# Initialize the fractional emittance parameters
        self.bin_edges = self.config.fractional_emittance_bins
        self.fraction = self.config.fractional_emittance_fraction
        self.uncertainty = self.config.fractional_emittance_uncertainty

	# Initialize the data containers
        self.a_dir = tempfile.mkdtemp()
        self.data_types = ("all_mc", "reco")
	self.data = {}
        self.amp_dict = {}
	for typ in self.data_types:
	    self.data[typ] = {}
	    self.amp_dict[typ] = {}

	# Calculate corrections if required
        self.calculate_corrections = self.config_anal["fractional_emittance_corrections"] == None

	# Initialize the systematics graphs if necessary
	if self.config_anal["fractional_emittance_systematics_draw"]:
	    self.syst_graphs = {}
	    for typ in self.data_types:
	    	self.syst_graphs[typ] = {}

    def birth(self):
        """
        Sets the output directory for the fractional emittance plots
	Resets all the data containers and loads the uncertainties
        """
	# Create directory, clear the data containers
        self.set_plot_dir("fractional_emittance")
	for typ in self.data_types:
	    self.data[typ].clear()

        self.get_planes()
        self.clear_feps_data()

	# Load the systematic corrections/uncertainties from the file if they are provided
	self.load_errors()

	# Load the data
        self.append_data()

        try:
            os.mkdir(self.plot_dir+"/corrections")
        except OSError:
            pass
        try:
            os.mkdir(self.plot_dir+"/systematics")
        except OSError:
            pass

    def process(self):
        """
        Imports the data to the fractional emittance format
        """
        self.append_data()

    def death(self):
        """
        Called when all the data has been loaded
        """
	# First, evaluate the amplitude quantiles and their statistical uncertainties
        self.make_amplitude_quantiles()

	# Evaluate corrections if necessary. Import them if they are provided
	# Draw the corrections in fractional_emittance/corrections
        if self.calculate_corrections:
	    self.set_corrections()
	if self.config_anal["fractional_emittance_corrections_draw"]:
	    self.draw_corrections()

	# Apply the corrections, evaluate the systematic uncertainties
	# If requested, draw the systematics in fractional_emittance/systematics
	self.corrections_and_uncertainties("reco")
        if self.config_anal["amplitude_mc"]:
            self.corrections_and_uncertainties("all_mc")
	if self.config_anal["fractional_emittance_systematics_draw"]:
	    self.draw_systematics()

	# Draw the fractional emittance evolution graphs
        self.draw_evolution()

	# Save everything to a json file
        self.save()

    def get_planes(self):
        """
        Initialize a data container for each of the data type.
	* all_mc has a container for each virtual plane between the two ref. planes included
	* reco and reco has a container for each ref. plane
        """
	# Initialize z position container, muon mass
        self.z_pos = {}
        mass = xboa.common.pdg_pid_to_mass[13]

	# Initialize data container for each reco plane
        for name in ['tku', 'tkd']:
            for z, dummy, z_name in self.config.detectors:
                if z_name == name+"_tp":
                    self.z_pos[name] = z
                    break

            file_name = self.a_dir+"/amp_data_"+name+"_"
            self.data["reco"][name] = AmplitudeDataBinned(file_name, self.bin_edges, mass, False)
            self.data["reco"][name].clear()

        print "Set up", len(self.data["reco"]), "reco planes"

	# If MC is present, initialize their data containers as well
	# Skip if the virtual plane is before tku or after tkd
        if self.config_anal["amplitude_mc"]:
            found_tku = False
            found_tkd = False
            for z, dummy, name in self.config.virtual_detectors:
                if "virtual_tku_tp" in name:
                    found_tku = True
                if not found_tku or found_tkd:
                    continue
                if "virtual_tkd_tp" in name:
                    found_tkd = True

                name = "mc_"+name
                file_name = self.a_dir+"/amp_data_"+name+"_"
                self.data["all_mc"][name] = AmplitudeDataBinned(file_name, self.bin_edges, mass, False)
                self.data["all_mc"][name].clear()
                self.z_pos[name] = z

            print "Set up", len(self.data["all_mc"]), "mc planes."


    def append_data(self):
        """
        Add data to the amplitude calculation (done at death time)
        
        If amplitude_mc is false, we take 'tku' data from upstream_cut sample 
        and 'tkd' data from downstream_cut sample
        
        If amplitude_mc is true, then we build an additional samples
        1. all_mc is for MC truth of all events that should have been included 
           in the sample (some of which might have been missed by recon; some in
           recon sample maybe should have been excluded)
        We use "mc_true_us_cut"/"mc_true_ds_cut" for the all_mc sample
        """
	hits = {}
        for typ in self.data_types:
	    hits[typ] = {}
	    for loc in self.data[typ].keys():
		hits[typ][loc] = []

        for event in self.data_loader.events:
            if event['upstream_cut']:
                continue
            hits["reco"]["tku"].append(event['tku'])
            if not event['downstream_cut']:
                hits["reco"]["tkd"].append(event['tkd'])

            if self.config_anal["amplitude_mc"]:
                for detector_hit in event["data"]:
                    det = detector_hit["detector"]
                    if det in self.data["all_mc"].keys():
                        hits["all_mc"][det].append(detector_hit["hit"])

        for typ in self.data_types:
	    for loc in self.data[typ].keys():
		self.data[typ][loc].append_hits(hits[typ][loc])

    def make_amplitude_quantiles(self):
        """
        Calculates an amplitude array for each of the datum point.
	Feeds the amplitude arrays to the AmplitudeQuantile object
	which computes the quantiles and their stat. uncertainties
        """
	# Build the amplitude dictionaries
	self.build_amplitude_dictionaries()

	# Get the number of particles in the upstream tracker
        us_events = len(self.amp_dict["reco"]["tku"])

	# Loop over the planes, get the quantiles
	for typ in self.data_types:
	    plane_id = 0
            for name, amps in self.amp_dict[typ].iteritems():
	        # Initialize the quantile object
                amp_values = amps.values()
                n_events = len(amp_values)		
	        norm = float(n_events)/us_events
	        quantile = AmplitudeQuantile(4, norm, self.uncertainty)

	        # Get the quantile
                feps = -1.
	        feps_stat = 0.
	        if self.fraction < norm:
		    quantile.set_quantile(amp_values, self.fraction)
                    feps = quantile.value()
                    feps_stat = quantile.error()

	        # Add a datum point for the current plane
                datum = {
                  'n_events':n_events,
                  'fraction':self.fraction,
                  'value':feps,
                  'stat_error':feps_stat,
	          'syst_error':0.,
	          'corrected':feps_stat,
                  'z_pos':self.z_pos[name],
		  'plane_id':plane_id
                }
		plane_id += 1
                json.dumps(datum)
                self.feps_data[typ]["quantiles"][name] = datum

    def build_amplitude_dictionaries(self):
        """
        Temporarily build amplitude dictionaries
        """
	for typ in self.data_types:
            for name, amp_data in self.data[typ].iteritems():
                try:
                    self.amp_dict[typ][name] = amp_data.fractional_amplitude()
                except Exception:
                    print "Failed to get amplitudes for", name

    def set_corrections(self):
        """
        Calculate the corrections, i.e. the inefficiency and response function
	The corrections are only relevant to the tku and tkd ref. planes
        * uses the Monte Carlo to generate corrections
        """
	for loc in ("tku", "tkd"):
	    # Evaluate the correction
	    mc_quantile = self.feps_data["all_mc"]["quantiles"]["mc_virtual_"+loc+"_tp"]["value"]
	    rec_quantile = self.feps_data["reco"]["quantiles"][loc]["value"]
	    corr = mc_quantile/rec_quantile

	    # Store it in the dictionary
	    self.feps_data["correction"][loc] = corr

    def corrections_and_uncertainties(self, typ):
        """
        Calculate corrected fractional emittance and uncertainties
        * typ specifies the type of data (all_mc, reco)
        * corrected quantile is given by A_alpha = c_alpha*A_alpha where c_alpha
            is the correction factor for this specific quantile
        * statistical errors are provided by the amplitude quantile calculator
        * total errors are given by the sum in quadrature of statistical errors and
          systematic errors
        """
	# Set the upstream statistical uncertainty to 0, if it is not an
	# MC data set, as it corresponds to the measurement to which to compare 
	# the downstream measurement.
        data = self.feps_data[typ]
        if not self.config_anal["amplitude_mc"] and typ == "reco":
	    data["reco"]["quantiles"]["tku"]["stat_error"] = 0.

	# Do the reco correction for each of the tracker locations
	if typ == "reco":
            for loc in ("tku", "tkd"):
                print "Doing fractional emittance correction for", typ, loc
	        source = self.feps_data
                quantile_value = self.do_corrections(typ, loc, source)
                data["quantiles"][loc]["corrected"] = quantile_value

	# If the systematics graphs are requested, initialize a graph per systematic sample
	if self.config_anal["fractional_emittance_systematics_draw"]:
	    for loc in ("us", "ds"):
	        for key in ["detector_systematics", "performance_systematics"]:
		    for source in data[loc][key]:
		        name = self.get_syst_name(source["source"])
			npoints = len(data["quantiles"])
		        self.syst_graphs[typ][name] = ROOT.TGraph(npoints)
		        self.syst_graphs[typ][name].SetName(name)
			for quantile in data["quantiles"].itervalues():
			    self.syst_graphs[typ][name].SetPoint(\
				quantile["plane_id"], quantile["z_pos"]/1e3, 0.)

	# Evaluate the systematic uncertainties for each of the measurement planes.
	# For the all_mc sample, consider the influence of the performance
	# systematics on all of the virtual planes (MC model uncertainty).
        for plane, quantile in data["quantiles"].iteritems():
            print "Finding fractional emittance systematic errors for", typ, plane
            reco_syst = self.calculate_detector_systematics(typ, plane)
            perf_syst = self.calculate_performance_systematics(typ, plane)
            syst_error = (reco_syst**2+perf_syst**2)**0.5
            quantile["syst_error"] = syst_error

        self.feps_data[typ] = data

    def do_corrections(self, typ, loc, source):
        """
        Applies the corrections to the requested quantile
        * typ specifies the type of data (all_mc, reco_mc, reco)
	* loc specifies the location of the tracker (tku, tkd)
	* source specifies the source of the corrections to be used
        """
	quantile = source[typ]["quantiles"][loc]["value"]
        corr = source["correction"][loc]
	quantile *= corr
	return quantile

    def calculate_detector_systematics(self, typ, plane):
        """
        Calculate the systematic errors in the reconstruction of tracks. The uncertainty
	on fractional corresponds to the sum in quadrature of the residuals between the reference
	reconstruction set and the data sets that are shifted from the reference set
        * typ specifies the type of data (all_mc, reco)
	* plane specifies the name of the measurement plane
        """
	# If it is all_mc, nothing to do here as there is no reconstruction uncertainty
	# If there is no reference specified, skip as there is no systematic uncertainty
        syst_error = 0.
        data = self.feps_data[typ]
        if typ == "all_mc" or data["detector_reference"] == None:
            return syst_error

        print "\nEvaluating fractional emittance reconstruction systematic errors in", plane

	# Correct the fractional emittance with the reference correction
        source = data["detector_reference"]
        ref_feps = self.do_corrections(typ, plane, source)

	# Loop over the detector systematics list for each measurement plane
	loc = {"tku":"us", "tkd":"ds"}[plane]
        systematics_list = data[loc]["detector_systematics"]
        for source in systematics_list:
	    # Evaluate the levels with the corresponding systematic shift
            syst_feps = self.do_corrections(typ, plane, source)

	    # Add in quadrature an uncertainty that corresponds to the level shift due
	    # to the use of a different set of corrections
            scale = source["scale"]
	    err = (syst_feps - ref_feps)*scale
            syst_error = (syst_error**2+err**2)**0.5

	    if self.config_anal["fractional_emittance_systematics_draw"]:
		name = self.get_syst_name(source["source"])
		quantile = data["quantiles"][plane]
		val = err/ref_feps
		self.syst_graphs[typ][name].GetY()[quantile["plane_id"]] = val

        return syst_error

    def calculate_performance_systematics(self, typ, plane):
        """
        Calculate the systematic errors in the channel performance. The experiment measures
	ratios of downstream fractional emittance over upstream fractional emittance.
	The uncertainty is evaluated as the shifts in the fractional emittance for
	variable deviations from the expected cooling channel performance.
        * typ specifies the type of data (all_mc, reco)
	* plane specifies the name of the measurement plane
        """
	# If the plane is in TKU, nothing to do here as there is no performance uncertainty
	# If there is no reference specified, skip as there is no systematic uncertainty
	syst_error = 0.
        data = self.feps_data[typ]
	if "tku" in plane or data["performance_reference"] == None:
            return syst_error

        print "\nEvaluating fractional emittance performance systematic errors in", plane

	# Get the reference ratio
        source = data["performance_reference"]
	quantiles = source[typ]["quantiles"]
	ref_key = {"all_mc":"mc_virtual_tku_tp", "reco":"tku"}[typ]
        ref_ratio = quantiles[plane]["value"]/quantiles[ref_key]["value"]

	# Loop over the performance systematics list
        systematics_list = data["ds"]["performance_systematics"]
	ratio_error = 0.
        for source in systematics_list:
	    # Evaluate the ratio with the corresponding systematic shift
	    quantiles = source[typ]["quantiles"]
            syst_ratio = quantiles[plane]["value"]/quantiles[ref_key]["value"]

	    # Add in quadrature an uncertainty that corresponds to the ratio shift due
	    # to the use of a different cooling channel
            scale = source["scale"]
            err = (syst_ratio - ref_ratio)*scale
            ratio_error = (ratio_error**2+err**2)**0.5

	    if self.config_anal["fractional_emittance_systematics_draw"]:
		name = self.get_syst_name(source["source"])
		quantile = data["quantiles"][plane]
		self.syst_graphs[typ][name].GetY()[quantile["plane_id"]] = err

	# Convert the uncertainties in terms of fractional emittance
        ref_feps = data["quantiles"][plane]["value"]
        syst_error = ratio_error * ref_feps

        return syst_error

    def get_syst_name(self, path):
	"""
	Convert systematic path to a systematic name
	"""
	suffix = path.split("Systematics_",1)[1]
	name = suffix.split("/")[0]
	return name

    def draw_systematics(self):
        """
        Draws the systematic errors. The uncertainty on each level corresponds to the 
	residuals between the reference reconstruction set and the data sets that 
	are shifted from the reference set
        """
	# Draw graphs that represent the deviations from the baseline
	# One multigraph per data category
        for typ in self.data_types:
	    # Initialize the multigraph
            title = "fractional_emittance_systematics_"+typ
            canvas = xboa.common.make_root_canvas(title)
	    mg = ROOT.TMultiGraph("mg_syst", ";z [m];Systematic residual")

	    # Loop over the systematic graphs, add them to the multigraph and the legend
	    leg = ROOT.TLegend(.6, .65, .8, .85)

	    gid = 0
	    for key, graph in self.syst_graphs[typ].iteritems():
		graph.Sort()
	        graph.SetLineColor(2+gid)
	        graph.SetMarkerStyle(21+gid)
		mg.Add(graph, "lp")
		leg.AddEntry(graph, key, "lp")
		gid += 1

	    mg.Draw("A")
	    leg.Draw("SAME")

	    # Add an envelope (quadratic sum of all sources of uncertainty)
	    graph_ref = self.syst_graphs[typ].itervalues().next()
            npoints = graph_ref.GetN()
            graph_upper = ROOT.TGraph(npoints)
            graph_lower = ROOT.TGraph(npoints)
	    graph_upper.SetLineWidth(2)
	    graph_lower.SetLineWidth(2)
	    graph_upper.SetMarkerStyle(20)
	    graph_lower.SetMarkerStyle(20)
	    names = []
	    for i in range(npoints):
	        quad_sum = 0.
	        for graph in self.syst_graphs[typ].itervalues():
	             quad_sum = (quad_sum**2 + graph.GetY()[i]**2)**0.5

	        z = graph_ref.GetX()[i]
	        graph_upper.SetPoint(i, z, quad_sum)
	        graph_lower.SetPoint(i, z, -quad_sum)
		names.append(graph.GetName())

	    mg.Add(graph_upper, "lp")
	    mg.Add(graph_lower, "lp")

	    # Add labels to the top of the envelope (measurement plane tags)
	    labels = {}
	    for plane, quantile in self.feps_data[typ]["quantiles"].iteritems():
		i = quantile["plane_id"]
	        labels[plane] = ROOT.TLatex(graph_upper.GetX()[i], graph_upper.GetY()[i], "  "+plane)
	        labels[plane].SetTextSize(0.025)
	        labels[plane].SetTextAngle(90)
		labels[plane].Draw("SAME")

	    # Print the canvas
            for fmt in ["root", "png", "pdf"]:
                canvas.Print(self.plot_dir+"/systematics/"+title+"."+fmt)

    def draw_evolution(self):
        """
        Produce a plot that shows the evolution of the fractional emittance
	across the MICE cooling channel. Displays the different types of data on
	shared canvas using the ROOT TMultiGraph class
        """
	# Initialize the canvas and the TMultiGraph
        canvas_name = 'fractional_emittance'
        canvas = self.get_plot(canvas_name)["pad"]
	mg = ROOT.TMultiGraph("mg_feps", ";z [m];#epsilon_{%d}  [mm]" % int(1e2*self.fraction))

	# Initialize the graphs
	graphs, graphs_tot = {}, {}
	markers = {"all_mc":24, "reco":20}
	colors = {"all_mc":4, "reco":1}
	draw_options= {"all_mc":"ple3", "reco":"p"}
	for typ in self.data_types:
	    # Get the quantiles
	    quantiles = self.feps_data[typ]["quantiles"]

	    # Add a graph that contains only statistical uncertainties
            point_list = [(datum['z_pos'], datum['value'], datum['stat_error']) \
                                		for datum in quantiles.itervalues()]
            graphs[typ] = self.make_graph(point_list, colors[typ], markers[typ], typ)
	    graphs[typ].SetLineWidth(2)
	    mg.Add(graphs[typ], draw_options[typ])

	    # Add a graph that contains the full uncertainties (quadratic sum of stat. and syst.)
            point_list = [(datum['z_pos'], datum['value'],\
				(datum['stat_error']**2+datum['syst_error']**2)**0.5)\
                                for datum in quantiles.itervalues()]
            graphs_tot[typ] = self.make_graph(point_list, colors[typ], markers[typ], typ)
	    mg.Add(graphs_tot[typ], draw_options[typ])

	# Initialize a legend
	leg = ROOT.TLegend(.6, .65, .8, .85)
	leg.AddEntry(graphs["all_mc"], graphs["all_mc"].GetName(), "plf")
	leg.AddEntry(graphs["reco"], graphs["reco"].GetName(), "p")

	# Draw
	mg.Draw("A")
	leg.Draw("SAME")
        for fmt in ["pdf", "png", "root"]:
            canvas.Print(self.plot_dir+"/fractional_emittance."+fmt)

    def make_graph(self, point_list, color, style, name):
        """
        Produce a single graph associated with a list of points
	* point_list is a list of points to be respresented, array of (z, val, err)
	* color specifies the color of the graph and markers
	* style specifies the marker style
	* name specifies the name of the graph
        """
        n_points = len(point_list)
        graph = ROOT.TGraphErrors(n_points)
	graph.SetTitle(";z [m];#epsilon_{%d}  [mm]" % int(1e2*self.fraction))
        graph.SetLineColor(color)
        graph.SetFillColorAlpha(color, .25)
        graph.SetMarkerColor(color)
        graph.SetMarkerStyle(style)
        graph.SetName(name)

        point_list = sorted(point_list)
        for i, point in enumerate(point_list):
            z = point[0]/1e3
            val = point[1]
            err = point[2]
            graph.SetPoint(i, z, val)
            graph.SetPointError(i, 0, err)

        self.plots["fractional_emittance"]["graphs"][name] = graph

        return graph

    def save(self):
        fout = open(self.plot_dir+"/fractional_emittance.json", "w")
        print >> fout, json.dumps(self.feps_data, sort_keys=True, indent=4)

    def clear_feps_data(self):
        """
        Initializes the dictionary that is used to store
	data and make corrections down the line
        """
	self.feps_data = {
          "correction":{
              "tku":1.,
              "tkd":1.
          },
          "source":"",
          "scale":1.
        }

        quantile_data = {
            "performance_reference":None,
            "detector_reference":None,
            "quantiles":{},
            "us":{
                "detector_systematics":[],
                "performance_systematics":[],
            },
            "ds":{
                "detector_systematics":[],
                "performance_systematics":[],
            }
        }
        for typ in self.data_types:
            self.feps_data[typ] = copy.deepcopy(quantile_data)

    def draw_corrections(self):
	"""
	Draw the correction factors used
	"""
	quantiles = self.feps_data["reco"]["quantiles"]
	corr_list = self.feps_data["correction"]

        title = "fractional_emittance_corrections"
        canvas = xboa.common.make_root_canvas(title)
	graph = ROOT.TGraph(2)
	graph.SetMarkerStyle(20)
	graph.SetTitle(";z [m];Correction factor")
	labels = {}
	for i, loc in enumerate(["tku", "tkd"]):
	    graph.SetPoint(i, quantiles[loc]["z_pos"]/1e3, corr_list[loc])
	    labels[loc] = ROOT.TLatex(graph.GetX()[i], graph.GetY()[i], loc)
	    labels[loc].SetTextSize(0.025)
	    labels[loc].SetTextAngle(90)
	
	graph.Draw("APL")
	for loc in quantiles.keys():
	    labels[loc].Draw("SAME")

        for fmt in ["root", "png", "pdf"]:
            canvas.Print(self.plot_dir+"/corrections/"+title+"."+fmt)
	

    def load_errors(self):
        """
        Two "classes" of systematic errors;
        * systematic errors on the reconstruction
        * systematic errors on the performance
        """
	# If the corrections are to calculated in this analysis, skip this step
        if self.calculate_corrections:
            return

        # Set base correction factors
        self.load_corrections(self.config_anal["fractional_emittance_corrections"])

	# Load systematic uncertainties
        systematics = self.config_anal["fractional_emittance_systematics"]
        for typ in systematics:
            print "Loading fractional emitttance systematic errors for", typ
            if typ not in self.feps_data:
                self.feps_data[typ] = {}
            for ref_key in ["detector_reference", "performance_reference"]:
                ref_src = systematics[typ][ref_key]
                if ref_src == None:
                    self.feps_data[typ][ref_key] = None
                else:
                    self.feps_data[typ][ref_key] = \
                                              self.load_one_error(ref_src, None)
                print "  Loaded reference", typ, ref_key, ref_src, \
                                          type(self.feps_data[typ][ref_key])
            for loc in ["us", "ds"]:
                if loc not in self.feps_data[typ]:
                    self.feps_data[typ][loc] = {}
                for key in ["detector_systematics", "performance_systematics"]:
                    err_src_dict = systematics[typ][loc][key]
                    self.feps_data[typ][loc][key] = [
                        self.load_one_error(err_src, scale) \
                              for err_src, scale in err_src_dict.iteritems()
                    ]
                    print "  Loaded", len(self.feps_data[typ][loc][key]), loc, key

    def load_corrections(self, file_name):
        """
        Load the frational emittance corrections to be applied during this 
        fractional emittance analysis. Loads the correction factors
        """
        fin = open(file_name)
        feps_str = fin.read()
        src_feps = json.loads(feps_str)
        src_feps["source"] = file_name
        self.feps_data["correction"] = src_feps["correction"]

    def load_one_error(self, file_name, scale):
        """
        Load the fractional emittance analysis output for a given uncertainty source
        """
        fin = open(file_name)
        feps_str = fin.read()
        feps = json.loads(feps_str)
        feps["source"] = file_name
        feps["scale"] = scale
        return feps
