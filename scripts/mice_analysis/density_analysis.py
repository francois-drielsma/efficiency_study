import math
import tempfile
import xboa.common
import json

import ROOT
import scipy.interpolate

from mice_analysis.density_estimator.density_data import DensityData
from mice_analysis.density_estimator.knn_density_estimator import kNNDensityEstimator
from mice_analysis.analysis_base import AnalysisBase


class DensityAnalysis(AnalysisBase):

    def __init__(self, config, config_anal, data_loader):
        """
        Initialise the DensityData class for the basic amplitude calculation
        * File name is the name of the (temporary) file to which the density data is written

        Internally, we store using numpy.memmap (which is a buffered file-based 
        array thingy). We have several memmaps:
        * run_array stores the run number of each event
        * spill_array stores the spill number of each event
        * event_array stores the event number of each event
        * ps_matrix stores the 4D phase space vector of each event
        Each one is stored in a memmap file, based on "file_name"+suffix
        During the amplitude calculation, for events whose amplitude we have
        not yet finished counting:
        * number of events in the covariance matrix
        """
        super(DensityAnalysis, self).__init__(config, config_anal, data_loader)
        self.config = config
        self.config_anal = config_anal
        self.data_loader = data_loader
        self.nthreads = self.config.density_nthreads
        self.knn_rotate = self.config.density_knn_rotate
        self.a_dir = tempfile.mkdtemp()
        file_name = self.a_dir+"/density_data_"
        self.all_mc_data_us = DensityData(file_name+"mc_us")
        self.reco_mc_data_us = DensityData(file_name+"reco_mc_us")
        self.reco_data_us = DensityData(file_name+"recon_us")
        self.all_mc_data_ds = DensityData(file_name+"mc_ds")
        self.reco_mc_data_ds = DensityData(file_name+"reco_mc_ds")
        self.reco_data_ds = DensityData(file_name+"recon_ds")

    def birth(self):
        self.set_plot_dir("density")
        self.all_mc_data_us.clear()
        self.reco_mc_data_us.clear()
        self.reco_data_us.clear()
        self.all_mc_data_ds.clear()
        self.reco_mc_data_ds.clear()
        self.reco_data_ds.clear()
        self.append_data()

    def process(self):
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
        hits_reco_us = []
        hits_all_mc_us = []
        hits_reco_mc_us = []
        hits_reco_ds = []
        hits_all_mc_ds = []
        hits_reco_mc_ds= []

        if self.config_anal["amplitude_mc"]:
            station_us = self.config.mc_plots["mc_stations"]["tku_tp"][0]
            station_ds = self.config.mc_plots["mc_stations"]["tkd_tp"][0]

        for event in self.data_loader.events:
            if event['upstream_cut']:
                continue
            hits_reco_us.append(event['tku'])
            if not event['downstream_cut']:
                hits_reco_ds.append(event['tkd'])

            if self.config_anal["amplitude_mc"]:
                hit_mc_us, hit_mc_ds = None, None
                for detector_hit in event["data"]:
                    if detector_hit["detector"] == station_us:
                        hit_mc_us = detector_hit["hit"]
                    if detector_hit["detector"] == station_ds:
                        hit_mc_ds = detector_hit["hit"]
                if not event['mc_true_us_cut']:
                    hits_all_mc_us.append(hit_mc_us)# no inefficiency upstream
                    hits_reco_mc_us.append(hit_mc_us)
                if not event['mc_true_ds_cut']:
                    hits_all_mc_ds.append(hit_mc_ds)
                    if not event['downstream_cut']:
                        hits_reco_mc_ds.append(hit_mc_ds)

        self.all_mc_data_us.append_hits(hits_all_mc_us)
        self.reco_mc_data_us.append_hits(hits_reco_mc_us)
        self.reco_data_us.append_hits(hits_reco_us)
        self.all_mc_data_ds.append_hits(hits_all_mc_ds)
        self.reco_mc_data_ds.append_hits(hits_reco_mc_ds)
        self.reco_data_ds.append_hits(hits_reco_ds)
        print "Loaded upstream (density):"
        print "    reco:   ", len(hits_reco_us)
        print "    all mc: ", len(hits_all_mc_us)
        print "    reco mc:", len(hits_reco_mc_us)
        print "Loaded downstream (density):"
        print "    reco:   ", len(hits_reco_ds)
        print "    all mc: ", len(hits_all_mc_ds)
        print "    reco mc:", len(hits_reco_mc_ds)
        return

    def make_graph(self, point_list, color, style, plot_option, name):
        n_points = len(point_list)
        point_list = sorted(point_list)
        graph = ROOT.TGraphAsymmErrors(n_points)
        for i, point in enumerate(point_list):
            z = point[0]
            amp = point[1]
            err_low = amp-point[2][0]
            err_high = point[2][1]-amp            
            graph.SetPoint(i, z, amp)
            graph.SetPointError(i, 0, 0, err_low, err_high)
        graph.SetLineColor(color)
        graph.SetMarkerColor(color)
        graph.SetMarkerStyle(style)
        graph.SetName(name)
        graph.Draw(plot_option)
        self.plots["fractional_amplitude"]["graphs"][name] = graph
        return graph

    def calculate_stats_errors(self):
	# Estimates the statistical uncertainty at each point of the graph (TODO)
        return

    def make_axes(self, canvas_name):
        hist = self.make_root_histogram(canvas_name, canvas_name,
                                        [-1], "z [mm]", 1000,
                                        [-1], "Amplitude [mm]", 1000, [],
                                        15000., 19000.,
                                        0., 100.)
        hist.Draw()

    def make_plot(self):
        n_fractions = len(self.fractions)
        mc_amp_list = []
        reco_amp_list = []
        canvas_name = 'fractional_amplitude'
        canvas = self.get_plot(canvas_name)["pad"]
        self.make_axes(canvas_name)
        reco_predicate = lambda datum: datum['name'] == "tku" or datum['name'] == "tkd"
        for i in range(n_fractions):
            name = "Reco "+str(self.fractions[i]*100)+" %"
            point_list = [(datum['z_pos'], datum['amplitude_bounds'][i], datum['stats_error'][i]) \
                                for datum in self.data if reco_predicate(datum)]
            self.make_graph(point_list, 1, 20, "P SAME", name)

        for i in range(n_fractions):
            name = "MC "+str(self.fractions[i]*100)+" %"
            point_list = [(datum['z_pos'], datum['amplitude_bounds'][i], datum['stats_error'][i]) \
                                for datum in self.data if not reco_predicate(datum)]
            self.make_graph(point_list, 4, 24, "P L SAME", name)

        for fmt in ["pdf", "png", "root"]:
            canvas.Print(self.plot_dir+"/density."+fmt)

    def performance_sys_error(self):
        pass

    def reco_sys_error(self):
        pass

    def save(self):
        for item in self.data:
            del item['values']
        fout = open(self.plot_dir+"/amplitude.json", "w")
        print >> fout, json.dumps(self.data)

    def death(self):
	print "TODO"
#        self.calculate_stats_errors()
#        self.make_plot()
#        self.save()

