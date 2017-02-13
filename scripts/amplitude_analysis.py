import sys

import operator
import xboa.common as common
from xboa.hit import Hit
from xboa.bunch import Bunch

class AmplitudeAnalysis(object):
    def __init__(self, config, config_anal):
        self.config = config
        self.config_anal = config_anal

    def delta_weight(self, bunch, delta):
        weight_sum_accepted = 0
        weight_sum_rejected = 0
        for hit in bunch:
            if (hit['spill'], hit['event_number']) in delta:
                weight_sum_accepted += hit['weight']
            else:
                weight_sum_rejected += hit['weight']
        return weight_sum_accepted, weight_sum_rejected

    def fractional_amplitude(self, bunch, amplitude_list, bunch_delta):
        amplitude_list = sorted(amplitude_list)[::-1]
        weight_accepted_list, weight_rejected_list = [], []
        bunch = bunch.deepcopy()
        delta = [(hit['spill'], hit['event_number']) for hit in bunch_delta]
        for amplitude in amplitude_list:
            old_weight = -1.
            new_weight = bunch.bunch_weight()
            while old_weight != new_weight:
                try:
                    old_weight = new_weight
                    bunch.cut({'amplitude x y':amplitude}, operator.gt)
                    new_weight = bunch.bunch_weight()
                except ZeroDivisionError:
                    sys.excepthook(*sys.exc_info())
                    new_weight = 0.
                    old_weight = 0.
            weight_sum_accepted, weight_sum_rejected = self.delta_weight(bunch, delta)
            weight_accepted_list.append(weight_sum_accepted)
            weight_rejected_list.append(weight_sum_rejected)
        weight_accepted_list = weight_accepted_list[::-1]
        weight_accepted_list = [0.]+weight_accepted_list
        weight_rejected_list = weight_rejected_list[::-1]
        weight_rejected_list = [0.]+weight_rejected_list
        return weight_accepted_list, weight_rejected_list

    #Add some histograms for sanity checks...
    def delta_amplitude_plot(self, bunch_0, bunch_1, plot_dir, p_bin):
        for hit in bunch_0[0:10]:
            print (hit['spill'], hit['event_number'], hit['particle_number'])
        for hit in bunch_1[0:10]:
            print (hit['spill'], hit['event_number'], hit['particle_number'])
        xmax = 100.
        nbins = 50
        bunch_good = bunch_0.deepcopy()
        bunch_good.transmission_cut(bunch_1)
        print "Momentum bin", p_bin


        print "bunch_0"
        print "    z:        ", bunch_0[0]["z"]
        print "    beta:     ", bunch_0.get_beta(["x", "y"])
        print "    emittance:", bunch_0.get_emittance(["x", "y"])
        print "    weight:   ", bunch_0.bunch_weight()
        print bunch_0.covariance_matrix(["x", "px", "y", "py"])
        print "bunch_good"
        print "    z:        ", bunch_good[0]["z"]
        print "    beta:     ", bunch_good.get_beta(["x", "y"])
        print "    emittance:", bunch_good.get_emittance(["x", "y"])
        print "    weight:   ", bunch_good.bunch_weight()
        print bunch_good.covariance_matrix(["x", "px", "y", "py"])
        print "bunch_1"
        print "    z:        ", bunch_1[0]["z"]
        print "    beta:     ", bunch_1.get_beta(["x", "y"])
        print "    emittance:", bunch_1.get_emittance(["x", "y"])
        print "    weight:   ", bunch_1.bunch_weight()
        print bunch_1.covariance_matrix(["x", "px", "y", "py"])
        axes = common.make_root_histogram("cdf 0", [-1.], "A_{#perp}   [mm]", nbins, xmin = 0., xmax = xmax)
        bin_edge_list = [0.]+[axes.GetBinCenter(bin)+axes.GetBinWidth(bin)/2. for bin in range(1, nbins+2)]

        print "Upstream amplitude calculation..."
        cdf_list_0, cdf_list_2 = self.fractional_amplitude(bunch_0, bin_edge_list, bunch_1)
        print "Downstream amplitude calculation..."
        cdf_list_1, dummy = self.fractional_amplitude(bunch_1, bin_edge_list, bunch_1)

        cdf_list_3 = []
        for bin_0, bin_2 in zip(cdf_list_0, cdf_list_2):
            cdf_list_3.append(bin_0+bin_2)
        print "CDF 0", cdf_list_0
        print "CDF 1", cdf_list_1
        print "CDF 2", cdf_list_2
        print "CDF 3", cdf_list_3
        print "DUMMY", dummy
        canvas_pdf = common.make_root_canvas("amplitudes")

        # Set up the histograms
        hist_3 = common.make_root_histogram("cdf 3", [-1.], "A_{#perp}   [mm]", nbins, xmin = 0., xmax = xmax)
        for bin, weight in enumerate(cdf_list_3[:-1]):
            hist_3.SetBinContent(bin+1, cdf_list_3[bin+1]-cdf_list_3[bin])
        hist_3.SetName("z = "+str(round(bunch_0[0]['z'], 1))+" mm (all)")

        hist_0 = common.make_root_histogram("cdf 0", [-1.], "A_{#perp}   [mm]", nbins, xmin = 0., xmax = xmax)
        for bin, weight in enumerate(cdf_list_0[:-1]):
            hist_0.SetBinContent(bin+1, cdf_list_0[bin+1]-cdf_list_0[bin])
        hist_0.SetName("z = "+str(round(bunch_0[0]['z'], 1))+" mm (observed downstream)")
        hist_0.SetLineColor(8)

        hist_1 = common.make_root_histogram("cdf 1", [-1.], "A_{#perp}   [mm]", nbins, xmin = 0., xmax = xmax)
        for bin, weight in enumerate(cdf_list_1[:-1]):
            hist_1.SetBinContent(bin+1, cdf_list_1[bin+1]-cdf_list_1[bin])
        hist_1.SetName("z = "+str(round(bunch_1[0]['z'], 1))+" mm (all)")
        hist_1.SetLineColor(4)

        hist_2 = common.make_root_histogram("cdf 2", [-1.], "A_{#perp}   [mm]", nbins, xmin = 0., xmax = xmax)
        for bin, weight in enumerate(cdf_list_2[:-1]):
            hist_2.SetBinContent(bin+1, cdf_list_2[bin+1]-cdf_list_2[bin])
        hist_2.SetName("z = "+str(round(bunch_0[0]['z'], 1))+" mm (not observed downstream)")
        hist_2.SetLineColor(2)

        # Draw the histograms (highest max bin first)
        hist_list = sorted([hist_3, hist_0, hist_1, hist_2], key = lambda hist: -hist.GetMaximum())
        hist_list[0].SetTitle(self.config_anal["name"])
        hist_list[0].Draw()
        for hist in hist_list[1:]:
            hist.Draw("SAME")
        hist_list[0].Draw("SAME") # draw again, on top of other histograms

        legend = common.make_root_legend(canvas_pdf, [hist_3, hist_0, hist_2, hist_1])
        legend.SetX1NDC(0.6)
        legend.SetX2NDC(0.89)
        legend.SetY1NDC(0.7)
        legend.SetY2NDC(0.89)
        legend.Draw()
        canvas_pdf.Update()

        canvas_ratio = common.make_root_canvas("amplitude ratio")
        hist_ratio = common.make_root_histogram("amplitude 1", [-1.], "A_{#perp}   [mm]", nbins, xmin = 0., xmax = xmax)
        sum_0, sum_1 = 0., 0.
        for bin in range(1, nbins+2):
            if hist_0.GetBinContent(bin) == 0:
                ratio = 0.
            else:
                ratio = hist_1.GetBinContent(bin)/float(hist_0.GetBinContent(bin))
            hist_ratio.SetBinContent(bin, ratio)
        cdf_ratio_list = []
        for i, cdf_3 in enumerate(cdf_list_3):
            if i == 0:
                continue
            if cdf_3 > 0.:
                cdf_ratio_list.append(cdf_list_1[i]/cdf_3)
            else:
                cdf_ratio_list.append(0.)
        print "cdf ratio list", cdf_ratio_list
        hist_cdf, graph = common.make_root_graph("ratio", bin_edge_list, "A_{#perp}   [mm]", cdf_ratio_list, "Ratio of cdf", xmin=0., xmax=xmax, ymin=0.8, ymax=1.4)
        hist_cdf.Draw()
        hist_cdf.SetTitle(self.config_anal["name"])
        graph.Draw("SAMEL")
        hist_ratio.SetMarkerStyle(24)
        hist_ratio.Draw("SAMEP")
        canvas_ratio.Update()
        for format in ["pdf", "png", "root"]:
            canvas_pdf.Print(plot_dir+"amplitude_pdf_p-"+str(p_bin)+"."+format)
            canvas_ratio.Print(plot_dir+"amplitude_ratio_p-"+str(p_bin)+"."+format)
        return canvas_pdf, canvas_ratio

    #['eventNumber', 'event_number', 'particleNumber', 'particle_number', 'pid', 'spill', 'station', 'status', '', "x'", "y'", 'p', "t'", 'kinetic_energy', 'ct']

    def make_bunches(self, data_loader, p_min, p_max):
        hit_list_us = []
        hit_list_ds = []
        hit_list_scraped = []

        for event in data_loader.events:
            if event['any_cut'] or event['tku']['p'] < p_min or event['tku']['p'] > p_max:
                continue

            if event['tku'] == None:
                continue
            hit_list_us.append(event['tku'])

            if event['tkd'] == None:
                continue
            hit_list_ds.append(event['tkd'])

        return Bunch.new_from_hits(hit_list_us),  Bunch.new_from_hits(hit_list_ds)


