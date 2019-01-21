import math
import tempfile
import xboa.common
import json

from mice_analysis.amplitude.amplitude_data_binned import AmplitudeDataBinned
from mice_analysis.analysis_base import AnalysisBase


class FractionalEmittance(AnalysisBase):
    def __init__(self, config, config_anal, data_loader):
        super(FractionalEmittance, self).__init__(config, config_anal, data_loader)
        self.config = config
        self.config_anal = config_anal
        self.data_loader = data_loader
        self.bin_edges = self.config.fractional_emittance_bins
        self.fractions = self.config.fractional_emittance_fractions
        self.a_dir = tempfile.mkdtemp()
        self.recon_planes = {}
        self.mc_planes = {}
        self.data = []
        self.z_pos = {}
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

    def performance_sys_error(self):
        pass

    def reco_sys_error(self):
        pass

    def stats_error(self, n_bin, n_tot):
        interval = 0.68 # 1 sigma
        bounds = BinomialConfidenceInterval.binomial_confidence_interval(
                                  int(n_bin),
                                  int(n_tot),
                                  interval)
        error = (bounds[0] - bounds[1])/2.
        return error

    def make_plot(self):
        z_list = [datum['z_pos'] for datum in self.data]
        n_fractions = len(self.fractions)
        amp_list = []
        frac_range = [0, 0]
        canvas_name = 'fractional_amplitude'
        canvas = self.get_plot(canvas_name)["pad"]
        for i in range(n_fractions):
            amp_list.append([datum['fractional_amplitude'][i] for datum in self.data])
            frac_range[1] = max(amp_list+[frac_range[1]])
        z_min = math.floor(min(z_list)/1000.)*1000.-1.
        z_max = math.floor(max(z_list)/1000.+1.)*1000.+1.
        for i in range(n_fractions):
            graph_name = self.config.fractional_emittance_fractions[i]*100
            graph_name = str(graph_name)+" %"
            hist, graph = self.make_root_graph(canvas_name, graph_name,
                z_list, "z [mm]",
                amp_list[i], "Amplitude [mm]", True,
                z_min, z_max, 0., 50.)
            if i == 0:
                hist.Draw()
            graph.SetMarkerStyle(24+i)
            graph.SetMarkerColor(41+i)
            graph.Draw("SAME P")
        for fmt in ["pdf", "png", "root"]:
            canvas.Print(self.plot_dir+"/fractional_emittance."+fmt)

    def do_analysis(self, dict_list):
        if dict_list[0][0] != 'tku':
            raise RuntimeError("Ooblicheck")
        n_events = len(dict_list[0][1])
        indices = [int(math.floor(n_events*frac)) for frac in self.fractions]
        for name, amp_values in dict_list:
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
              'fractional_amplitude':fractions,
              'stats_error':stats_error,
              'z_pos':self.z_pos[name],
              #'values':[float(value) for value in amp_values],
            }
            json.dumps(datum)
            self.data.append(datum)
        fout = open(self.plot_dir+"/amplitude.json", "w")
        print >> fout, json.dumps(self.data)

    def death(self):
        dict_tku = self.recon_planes['tku'].fractional_amplitude()
        dict_list = [
            ('tku', dict_tku),
            ('tkd', self.recon_planes['tkd'].fractional_amplitude())
        ]
        for name, amp_data in self.mc_planes.iteritems():
            try:
                dict_list.append((name, amp_data.fractional_amplitude()))
            except Exception:
                print "Failed to get amplitudes for", name
        self.do_analysis(dict_list)
        self.make_plot()

