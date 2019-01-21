import math
import tempfile
import xboa.common
import json

import ROOT
import scipy.interpolate

from utilities.binomial_confidence_interval import BinomialConfidenceInterval

from mice_analysis.amplitude.amplitude_data_binned import AmplitudeDataBinned
from mice_analysis.analysis_base import AnalysisBase


class FractionalEmittance(AnalysisBase):
    def __init__(self, config, config_anal, data_loader):
        super(FractionalEmittance, self).__init__(config, config_anal, data_loader)
        self.config = config
        self.config_anal = config_anal
        self.data_loader = data_loader
        self.bin_edges = self.config.fractional_emittance_bins
        self.fractions = sorted(self.config.fractional_emittance_fractions)
        self.a_dir = tempfile.mkdtemp()
        self.recon_planes = {} # reco raw data for amplitude calculation
        self.mc_planes = {} # mc raw data for amplitude calculation
        self.dict_list = [] # the calculated amplitudes
        self.data = [] # data for plotting/storing on disk
        self.z_pos = {}
        self.hist_drawn = False
        self.get_planes()

    def birth(self):
        self.set_plot_dir("fractional_emittance")
        self.data = []
        self.append_data()

    def get_planes(self):
        self.recon_planes = {}
        self.mc_planes = {}
        self.z_pos = {}
        mass = xboa.common.pdg_pid_to_mass[13]
        found_tku = False
        found_tkd = False
        for name in ['tku', 'tkd']:
            file_name = self.a_dir+"/amp_data_"+name+"_"
            self.recon_planes[name] = AmplitudeDataBinned(file_name, self.bin_edges, mass, False)
            self.recon_planes[name].clear()
            for z, dummy, z_name in self.config.detectors:
                if z_name == name+"_tp":
                    self.z_pos[name] = z
                    break
        print "fractional_emittance... Got recon planes...",
        if not self.config_anal["amplitude_mc"]: # not doing amplitude mc
            return
        for z, dummy, name in self.config.virtual_detectors:
            if "virtual_tku_tp" in name:
                found_tku = True
            if not found_tku or found_tkd:
                continue
            if "virtual_tkd_tp" in name:
                found_tkd = True
            name = "mc_"+name
            file_name = self.a_dir+"/amp_data_"+name+"_"
            self.mc_planes[name] = AmplitudeDataBinned(file_name, self.bin_edges, mass, False)
            self.mc_planes[name].clear()
            self.z_pos[name] = z
        print "Set up", len(self.mc_planes), "mc planes."

    def process(self):
        self.append_data()

    def append_data(self):
        recon_hits_tku = []
        recon_hits_tkd = []
        mc_dets = self.mc_planes.keys()
        mc_hits = dict([(key, []) for key in mc_dets])

        for event in self.data_loader.events:
            if event['upstream_cut']:
                continue
            recon_hits_tku.append(event['tku'])
            if not event['downstream_cut']:
                recon_hits_tkd.append(event['tkd'])

            if self.config_anal["amplitude_mc"]:
                for detector_hit in event["data"]:
                    det = detector_hit["detector"]
                    if det in mc_dets:
                        mc_hits[det].append(detector_hit["hit"])
        self.recon_planes['tku'].append_hits(recon_hits_tku)
        self.recon_planes['tkd'].append_hits(recon_hits_tkd)
        for det in mc_dets:
            self.mc_planes[det].append_hits(mc_hits[det])

    def build_amplitude_dictionaries(self):
        dict_tku = self.recon_planes['tku'].fractional_amplitude()
        self.dict_list = [
            ('tku', dict_tku),
            ('tkd', self.recon_planes['tkd'].fractional_amplitude())
        ]
        for name, amp_data in self.mc_planes.iteritems():
            try:
                self.dict_list.append((name, amp_data.fractional_amplitude()))
            except Exception:
                print "Failed to get amplitudes for", name

    def get_amplitude_bounds(self):
        if self.dict_list[0][0] != 'tku':
            raise RuntimeError("Found event with no tku data while doing fractional emittance")
        n_events = len(self.dict_list[0][1])
        indices = [int(math.floor(n_events*frac)) for frac in self.fractions]
        for name, amp_values in self.dict_list:
            amp_values = sorted(amp_values.values())
            n_events = len(amp_values)
            fractions = []
            for index in indices:
                if index < n_events:
                    fractions.append(float(amp_values[index]))
                else:
                    fractions.append(-1)
            datum = {
              'name':name,
              'n_events':n_events,
              'n_bounds':indices,
              'amplitude_bounds':fractions,
              'z_pos':self.z_pos[name],
              'values':[float(value) for value in amp_values],
            }
            json.dumps(datum)
            self.data.append(datum)

    def calculate_stats_errors(self):
        interval = 0.68 # 1 sigma
        # number of particles in each bin
        n_bounds = self.data[0]['n_bounds']
        # binomial distn uncertainty in number of particles in the bin
        n_errors = [BinomialConfidenceInterval.binomial_confidence_interval(
                                  n,
                                  n_bounds[-1],
                                  interval) for n in n_bounds]
        print "fractional_emittance ... calculated stats errors in number", n_errors
        # corresponding uncertainty in amplitude
        for datum in self.data:
            datum['stats_error'] = []
            y = datum['values']
            x = [float(i) for i in range(len(y))]
            interp = scipy.interpolate.interp1d(x, y, bounds_error=False, fill_value=(0., 1000.))
            for an_n_error in n_errors:
                amp_error = [float(interp(an_n_error[0])),
                             float(interp(an_n_error[1]))]
                datum['stats_error'].append(amp_error)
            print "fractional_emittance ... calculated stats errors in amp",
            print datum['stats_error']
        for datum in self.data:
            print datum.keys()
        return

    def make_graph(self, z_list, amp_list, err_list, color, style, plot_option, name):
        n_points = len(z_list)
        graph = ROOT.TGraphAsymmErrors(len(z_list))
        for i in range(n_points):
            graph.SetPoint(i, z_list[i], amp_list[i])
            y_low = amp_list[i]-err_list[i][0]
            y_high = err_list[i][1]-amp_list[i]
            graph.SetPointError(i, 0, 0, y_low, y_high)
        graph.SetLineColor(color)
        graph.SetMarkerColor(color)
        graph.SetMarkerStyle(style)
        graph.SetName(name)
        graph.Draw(plot_option)
        self.plots["fractional_amplitude"]["graphs"][name] = graph
        return graph

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
            z_list = [datum['z_pos'] for datum in self.data if reco_predicate(datum)]
            amp_list = [datum['amplitude_bounds'][i] for datum in self.data \
                             if reco_predicate(datum)]
            err_list = [datum['stats_error'][i] for datum in self.data \
                             if reco_predicate(datum)]
            self.make_graph(z_list, amp_list, err_list, 1, 20, "P SAME", name)

        for i in range(n_fractions):
            name = "MC "+str(self.fractions[i]*100)+" %"
            z_list = [datum['z_pos'] for datum in self.data if not reco_predicate(datum)]
            amp_list = [datum['amplitude_bounds'][i] for datum in self.data \
                             if not reco_predicate(datum)]
            err_list = [datum['stats_error'][i] for datum in self.data \
                             if not reco_predicate(datum)]
            self.make_graph(z_list, amp_list, err_list, 4, 24, "P L SAME", name)

        for fmt in ["pdf", "png", "root"]:
            canvas.Print(self.plot_dir+"/fractional_emittance."+fmt)

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
        self.build_amplitude_dictionaries()
        self.get_amplitude_bounds()
        self.calculate_stats_errors()
        self.make_plot()
        self.save()

