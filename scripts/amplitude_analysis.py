import sys
import datetime
import bisect
import operator
import json
import numpy

import ROOT
import xboa.common as common
from xboa.hit import Hit
from xboa.bunch import Bunch
import scripts.chi2_distribution
from scripts.binomial_confidence_interval import BinomialConfidenceInterval

class AmplitudeAnalysis(object):
    def __init__(self, config, config_anal):
        self.config = config
        self.config_anal = config_anal
        self.amplitudes = {}
        self.bin_edge_list = [float(i) for i in range(0, 101, config.amplitude_bin_width)]

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
        """
        bunch: all particles, on which amplitude is calculated
        amplitude_list: the set of bin edges over which amplitude is calculated
        bunch_delta: the set of particles which make it downstream
        """
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
                except Exception:
                    sys.excepthook(*sys.exc_info())
                    new_weight = 0.
                    old_weight = 0.
            weight_sum_accepted, weight_sum_rejected = self.delta_weight(bunch, delta)
            weight_accepted_list.append(weight_sum_accepted)
            weight_rejected_list.append(weight_sum_rejected)
        weight_accepted_list = weight_accepted_list[::-1]
        #weight_accepted_list = [0.]+weight_accepted_list
        weight_rejected_list = weight_rejected_list[::-1]
        #weight_rejected_list = [0.]+weight_rejected_list
        return weight_accepted_list, weight_rejected_list

    def fractional_amplitude_2(self, bunch):
        """
        bunch: all particles, on which amplitude is calculated
        amplitude_list: the set of bin edges over which amplitude is calculated
        bunch_delta: the set of particles which make it downstream
        """
        amplitude_bin_edges = sorted(self.bin_edge_list)[::-1] # bins
        bunch = bunch.deepcopy()
        bunch.clear_weights()
        amplitude_dict = dict([((hit['spill'], hit['event_number']), None) for hit in bunch]) # mapping of spill/event to amplitude

        for bin_upper_edge in amplitude_bin_edges[:-1]:
            new_weight = -1
            old_weight = 1
            print "Calc bin", bin_upper_edge,
            while old_weight != new_weight:
                old_weight = new_weight
                try:
                    amplitude_list = bunch.list_get_hit_variable(['amplitude x y'])[0]
                    for i, amplitude in enumerate(amplitude_list):
                        #print amplitude,
                        hit = bunch[i]
                        if amplitude > bin_upper_edge:
                            hit['local_weight'] = 0
                        else:
                            amplitude_dict[(hit['spill'], hit['event_number'])] = amplitude
                    new_weight = bunch.bunch_weight()
                except Exception:
                    sys.excepthook(*sys.exc_info())
                    break
            print new_weight
        #for key, value in amplitude_dict.iteritems():
        #    print key, value
        return amplitude_dict

    def delta_amplitude_calc_2(self, bunch_0, bunch_1, suffix):
        print "Amplitude", suffix
        print "  starting delta_amplitude_calc_2    ", datetime.datetime.now()
        data = {}
        print "  Upstream amplitude calculation...  ", datetime.datetime.now()
        fractional_amplitude_dict_0 = self.fractional_amplitude_2(bunch_0)
        print "  Downstream amplitude calculation...", datetime.datetime.now()
        fractional_amplitude_dict_1 = self.fractional_amplitude_2(bunch_1)
        print "  ... done                           ", datetime.datetime.now()


        data["all_upstream"] = self.get_pdfs(fractional_amplitude_dict_0)
        data["all_downstream"] = self.get_pdfs(fractional_amplitude_dict_1)
        for name, bunch in [("all_upstream", bunch_0),
                           ("all_downstream", bunch_1)]:
            data[name]["z"] = bunch[0]["z"]
            data[name]["beta_4D"] = bunch.get_beta(["x", "y"])
            data[name]["emittance"] = bunch.get_emittance(["x", "y"])
            data[name]["weight"] = bunch.bunch_weight()
            data[name]["covariance_matrix"] = bunch.covariance_matrix(["x", "px", "y", "py"]).tolist()
            data[name]["bin_edge_list"] = self.bin_edge_list

        data["migration_matrix"] = self.migration_matrix(fractional_amplitude_dict_0, fractional_amplitude_dict_1)
        data["upstream_scraped"], data["upstream_not_scraped"] = self.get_delta_pdfs(fractional_amplitude_dict_0, fractional_amplitude_dict_1)
        self.amplitudes[suffix] = data
        print "  done amplitude calc                ", datetime.datetime.now()

    def get_bin_centre_list(self):
        bin_width = self.config.amplitude_bin_width
        bin_centre_list = [bin_edge+bin_width/2. for bin_edge in self.bin_edge_list[:-1]]
        return bin_centre_list

    def get_delta_pdfs(self, a_dict_us, a_dict_ds):
        scraped_dict = {}
        not_scraped_dict = {}
        for key, value in a_dict_us.iteritems():
            if key in a_dict_ds:
                not_scraped_dict[key] = value
            else:
                scraped_dict[key] = value
        upstream_scraped = self.get_pdfs(scraped_dict)
        upstream_not_scraped = self.get_pdfs(not_scraped_dict)
        return upstream_scraped, upstream_not_scraped


    def get_pdfs(self, fractional_amplitude_dict):
        bin_centre_list = self.get_bin_centre_list()
        data = {}
        data["bin_centre_list"] = bin_centre_list
        data["cdf"] = [0]*(len(self.bin_edge_list))
        data["pdf"] = [0]*(len(self.bin_edge_list))
        amplitude_index = 0
        for amplitude in sorted(fractional_amplitude_dict.values()):
            while amplitude > self.bin_edge_list[amplitude_index]:
                amplitude_index += 1
            data["cdf"][amplitude_index-1] += 1
            data["pdf"][amplitude_index-1] += 1
        return data

    def migration_matrix(self, fractional_amplitude_dict_0, fractional_amplitude_dict_1):
        migration_matrix = [[0. for i in range(len(self.bin_edge_list)+1)]  for j in range(len(self.bin_edge_list))]
        # migration dict maps event number onto amplitude tuple (us, ds)
        migration_dict = {}
        for key, amp_0 in fractional_amplitude_dict_0.iteritems():
            amp_1 = None
            if key in fractional_amplitude_dict_1:
                amp_1 = fractional_amplitude_dict_1[key]
            migration_dict[key] = (amp_0, amp_1)
        for key, amp_pair in migration_dict.iteritems():
            amp_0_bin = bisect.bisect_left(self.bin_edge_list[1:], amp_pair[0])
            if amp_pair[1] == None:
                amp_1_bin = len(self.bin_edge_list)
            else:
                amp_1_bin = bisect.bisect_left(self.bin_edge_list[1:], amp_pair[1])
            migration_matrix[amp_0_bin][amp_1_bin] += 1
        return migration_matrix

    #Add some histograms for sanity checks...
    def delta_amplitude_calc(self, bunch_0, bunch_1, suffix, bin_width):
        data = {}
        bunch_good = bunch_0.deepcopy()
        bunch_good.transmission_cut(bunch_1)
        print "Amplitude", suffix
        for name, bunch in [("all_upstream", bunch_0),
                           ("all_downstream", bunch_1)]:
            data[name] = {
              "z":bunch[0]["z"],
              "beta_4D":bunch.get_beta(["x", "y"]),
              "emittance":bunch.get_emittance(["x", "y"]),
              "weight":bunch.bunch_weight(),
              "covariance_matrix":bunch.covariance_matrix(["x", "px", "y", "py"]),
              "bin_edge_list":self.bin_edge_list,
            }

        print "Upstream amplitude calculation..."
        cdf_list_0, cdf_list_2 = self.fractional_amplitude(bunch_0, self.bin_edge_list, bunch_1)
        print "Downstream amplitude calculation..."
        cdf_list_1, dummy = self.fractional_amplitude(bunch_1, self.bin_edge_list, bunch_1)

        cdf_list_3 = []
        for bin_0, bin_2 in zip(cdf_list_0, cdf_list_2):
            cdf_list_3.append(bin_0+bin_2)

        data["all_upstream"]["cdf"] = cdf_list_3
        data["all_upstream"]["bin_centre_list"], data["all_upstream"]["pdf"] = self.pdf(cdf_list_3)

        data["all_downstream"]["cdf"] = cdf_list_1
        data["all_downstream"]["bin_centre_list"], data["all_downstream"]["pdf"] = self.pdf(cdf_list_1)
        self.amplitudes[suffix] = data

    def pdf(self, cdf_list):
        y_value = [cdf_list[i+1]-cdf_list[i] for i, lower in enumerate(self.bin_edge_list[:-1])]
        x_value = [(lower+bin_edges_list[i+1])/2. for i, lower in enumerate(self.bin_edge_list[:-1])]
        return x_value, y_value

    def amplitude_pdf_graph(self, x_value, y_value, normalisation, name, interval):
        print x_value, y_value
        x_err = (x_value[1]-x_value[0])/2.
        n_total = sum(y_value)
        y_errs = []
        for n_in_bin in y_value:
            print "calculating confidence interval for", n_in_bin, n_total, interval
            y_errs.append(BinomialConfidenceInterval.binomial_confidence_interval(int(n_in_bin), int(n_total), interval))
            print "   ...", y_errs[-1]

        graph = ROOT.TGraphAsymmErrors(len(x_value))
        for i, x in enumerate(x_value):
            graph.SetPoint(i, x, y_value[i]/normalisation)
            print i, y_value[i], y_errs[i][0], y_errs[i][1], y_value[i]/normalisation, y_errs[i][0]/normalisation, y_errs[i][1]/normalisation
            graph.SetPointError(i, x_err, x_err, y_errs[i][0]/normalisation, y_errs[i][1]/normalisation)
            graph.SetPointError(i, x_err, x_err, y_errs[i][0]/normalisation, y_errs[i][1]/normalisation)
            graph.SetPointEYlow(i, y_errs[i][0]/normalisation)
            graph.SetPointEYhigh(i, y_errs[i][1]/normalisation)
        self.root_objects.append(graph)
        graph.SetName(name)
        graph.Draw("SAMEP")
        return graph

    def field_uncertainty_migration(self, i, pdf_1, pdf_2):
        """Fraction migrating from bin i_1 to bin i_1 + 1"""
        # fraction of the bin that overlaps i_1+1
        # (b2*(1+alpha)-b2*alpha)*(b2 - b1) = b2*alpha/(b2-b1)
        alpha = self.config_anal["field_uncertainty"]
        b1 = self.bin_edge_list[i]
        b2 = self.bin_edge_list[i+1]
        number_moving_up = pdf_1*b2*alpha/(b2-b1)
        return number_moving_up

    def systematics_calc(self, target):
        all_mc_pdf = self.amplitudes["all_mc"][target]["pdf"]
        reco_mc_pdf = self.amplitudes["reco_mc"][target]["pdf"]
        reco_pdf = self.amplitudes["reco"][target]["pdf"]

        efficiency_list = []
        efficiency_err = []
        for i, mc in enumerate(all_mc_pdf):
            if mc != 0:
                delta = mc-reco_mc_pdf[i]
                err = BinomialConfidenceInterval.binomial_confidence_interval(int(delta), int(mc), 0.68)
                efficiency_list.append((delta)/float(mc))
                efficiency_err.append([(delta-err[0])/float(mc), (err[1]-delta)/float(mc)]) # sqrt(n)/n
            else:
                efficiency_list.append(0.)
                efficiency_err.append([0., 0.])
        if "inefficiency" not in self.amplitudes:
            self.amplitudes["inefficiency"] = {}
        self.amplitudes["inefficiency"][target] = {}
        self.amplitudes["inefficiency"][target]["pdf"] = efficiency_list
        self.amplitudes["inefficiency"][target]["err"] = efficiency_err

        crossing_probability_list = []
        crossing_probability_err = []
        for i, reco_mc in enumerate(reco_mc_pdf):
            if reco_mc != 0:
                delta = reco_mc-reco_pdf[i]
                crossing_probability_list.append(delta/float(reco_mc))
                crossing_probability_err.append(float(abs(delta))**0.5/reco_mc) # sqrt(n)/n_total
            else:
                crossing_probability_list.append(0.)
                crossing_probability_err.append(0.)
        if "crossing_probability" not in self.amplitudes:
            self.amplitudes["crossing_probability"] = {}
        self.amplitudes["crossing_probability"][target] = {}
        self.amplitudes["crossing_probability"][target]["pdf"] = crossing_probability_list
        self.amplitudes["crossing_probability"][target]["err"] = crossing_probability_err

    def field_uncertainty_calc(self, target):
        reco_pdf = self.amplitudes["reco"][target]["pdf"]
        n_bins = len(reco_pdf)
        migration_list = [0.]+[self.field_uncertainty_migration(i, reco_pdf[i], reco_pdf[i+1]) for i in range(n_bins-1)]+[0.]
        momentum_uncertainty_list = []
        momentum_uncertainty_err = [0. for i in range(n_bins)]
        momentum_uncertainty = self.config_anal["field_uncertainty"]
        for i in range(n_bins):
            if reco_pdf[i] != 0.:
                momentum_uncertainty_list.append(abs(migration_list[i] - migration_list[i+1])/reco_pdf[i])
            else:
                momentum_uncertainty_list.append(0.)
        if "momentum_uncertainty" not in self.amplitudes:
            self.amplitudes["momentum_uncertainty"] = {}
        self.amplitudes["momentum_uncertainty"][target] = {}
        self.amplitudes["momentum_uncertainty"][target]["pdf"] = momentum_uncertainty_list
        self.amplitudes["momentum_uncertainty"][target]["err"] = momentum_uncertainty_err
        

    def systematics_plot(self, target):
        canvas_efficiency = common.make_root_canvas("inefficiency")
        bin_centre_list = self.amplitudes["all_mc"][target]["bin_centre_list"]
        efficiency_list = self.amplitudes["inefficiency"][target]["pdf"]
        err_list = self.amplitudes["inefficiency"][target]["err"]
        graph = ROOT.TGraphAsymmErrors(len(bin_centre_list))
        self.root_objects.append(graph)
        ymin = -1.
        ymax = +1.
        bin_width = bin_centre_list[1]-bin_centre_list[0]
        for i, bin_centre in enumerate(bin_centre_list):
            graph.SetPoint(i, bin_centre, efficiency_list[i])
            graph.SetPointError(i, bin_width/2., bin_width/2., err_list[i][0], err_list[i][1])
            ymin = min(ymin, efficiency_list[i]-err_list[i][0])
            ymax = max(ymax, efficiency_list[i]+err_list[i][1])

        # Set up the histograms
        hist = common.make_root_histogram("inefficiency", [0.], "A_{#perp}   [mm]", 100,
                                          [0.], "Inefficiency", 100,
                                          xmin=min(self.bin_edge_list), xmax=max(self.bin_edge_list),
                                          ymin=ymin*1.1, ymax=ymax*1.1) # last entry is the overflow bin
        hist.Draw()
        graph.SetMarkerStyle(24)
        graph.Draw("SAMEP")
        canvas_efficiency.Update()

        canvas_crossing = common.make_root_canvas("crossing probability")
        bin_centre_list = self.amplitudes["all_mc"][target]["bin_centre_list"]
        crossing_probability_list = self.amplitudes["crossing_probability"][target]["pdf"]
        err_list = self.amplitudes["crossing_probability"][target]["err"]
        graph = ROOT.TGraphErrors(len(bin_centre_list))
        self.root_objects.append(graph)
        ymin = -1.
        ymax = +1.
        bin_width = bin_centre_list[1]-bin_centre_list[0]
        for i, bin_centre in enumerate(bin_centre_list):
            graph.SetPoint(i, bin_centre, crossing_probability_list[i])
            graph.SetPointError(i, bin_width/2., err_list[i])
            ymin = min(ymin, crossing_probability_list[i]-err_list[i])
            ymax = max(ymax, crossing_probability_list[i]+err_list[i])

        # Set up the histograms
        hist = common.make_root_histogram("inefficiency", [0.], "A_{#perp}   [mm]", 100,
                                          [0.], "Tracker Resolution Crossing Probability", 100,
                                          xmin=min(self.bin_edge_list), xmax=max(self.bin_edge_list),
                                          ymin=ymin*1.1, ymax=ymax*1.1)
        hist.Draw()
        graph.SetMarkerStyle(24)
        graph.Draw("SAMEP")
        canvas_crossing.Update()
        plot_dir = self.config_anal["plot_dir"]
        for format in ["eps", "png", "root"]:
            canvas_efficiency.Print(plot_dir+"amplitude_efficiency_"+target+"."+format)
            canvas_crossing.Print(plot_dir+"amplitude_crossing_"+target+"."+format)

    def field_uncertainty_plot(self, target):
        canvas_field_uncertainty = common.make_root_canvas("field uncertainty")
        bin_centre_list = self.amplitudes["reco"][target]["bin_centre_list"]
        field_uncertainty_list = self.amplitudes["momentum_uncertainty"][target]["pdf"]
        err_list = self.amplitudes["momentum_uncertainty"][target]["err"]
        graph = ROOT.TGraphErrors(len(bin_centre_list))
        self.root_objects.append(graph)
        ymin = -1.
        ymax = +1.
        bin_width = bin_centre_list[1]-bin_centre_list[0]
        for i, bin_centre in enumerate(bin_centre_list):
            graph.SetPoint(i, bin_centre, field_uncertainty_list[i])
            graph.SetPointError(i, bin_width/2., err_list[i])
            ymin = min(ymin, field_uncertainty_list[i]-err_list[i])
            ymax = max(ymax, field_uncertainty_list[i]+err_list[i])

        # Set up the histograms
        hist = common.make_root_histogram("field uncertainty", [0.], "A_{#perp}   [mm]", 100,
                                          [0.], "Field Uncertainty Crossing Probability", 100,
                                          xmin=min(self.bin_edge_list), xmax=max(self.bin_edge_list),
                                          ymin=ymin*1.1, ymax=ymax*1.1)
        hist.Draw()
        graph.SetMarkerStyle(24)
        graph.Draw("SAMEP")
        canvas_field_uncertainty.Update()


        plot_dir = self.config_anal["plot_dir"]
        for format in ["eps", "png", "root"]:
            canvas_field_uncertainty.Print(plot_dir+"amplitude_field_uncertainty_"+target+"."+format)

    def matrix_str(self, matrix, rounding=2, width=10):
        matrix_str = ""
        for row in matrix:
            try:
                for element in row:
                    if type(element) == type(1.):
                        element = round(element, rounding)
                    matrix_str += str(element).rjust(width)
                matrix_str += "\n"
            except TypeError:
                if type(row) == type(1.):
                    row = round(row, rounding)
                    matrix_str += str(row).rjust(width) # maybe it was a vector
        return matrix_str

    def sys_errors(self, target, pdf_list):
        sys_errors = [0. for bin in self.bin_edge_list]
        if "crossing_probability" not in self.amplitudes:
            return sys_errors
        crossing_probability = self.amplitudes["crossing_probability"][target]["pdf"]
        efficiency = self.amplitudes["inefficiency"][target]["pdf"]
        momentum_uncertainty = self.amplitudes["momentum_uncertainty"][target]["pdf"]
        for i, pdf in enumerate(pdf_list):
            sys_errors[i] = (efficiency[i]**2+
                             crossing_probability[i]**2+
                             momentum_uncertainty[i]**2)**0.5
            sys_errors[i] *= pdf
        return sys_errors

    def stats_errors(self, migration_matrix):
        # migration_matrix M_ij is number of events having 
        # upstream bin i AND downstream bin j
        # statistical error on number in bin M_ij is binomially distributed
        # statistical error on number in column j is sum of errors in each M_ij
        # (where we sum errors in quadrature)
        row_sum = [sum(row) for row in migration_matrix]
        n = len(migration_matrix)
        interval = 0.68 # "1 sigma confidence interval"
        error_matrix = [[0. for i in range(n)] for j in range(n)]
        for i in range(n):
            if int(row_sum[i]) == 0:
                for j in range(n):
                    error_matrix[j][i] = 0.
                continue
            for j in range(n):
                bounds = BinomialConfidenceInterval.binomial_confidence_interval(
                                          int(migration_matrix[i][j]),
                                          int(row_sum[i]),
                                          interval)
                error = (bounds[0] - bounds[1])/2.
                # Note reversed indices because we want the error on downstream
                # bins (columns of migration_matrix which become rows of
                # error_matrix)
                # Square because we add errors in quadrature
                error_matrix[j][i] = error**2
        # Add errors in quadrature
        downstream_errors = [sum(row)**0.5 for row in error_matrix]
        print "upstream pdf\n", self.matrix_str(row_sum)
        print "migration matrix\n", self.matrix_str(migration_matrix)
        print "error matrix\n", self.matrix_str(error_matrix)
        print "downstream errors\n", self.matrix_str(downstream_errors)
        return downstream_errors

    def get_delta_graphs(self, suffix):
        data = self.amplitudes[suffix]
        pdf_list_tku = data["all_upstream"]["pdf"]
        scraped_list = data["upstream_scraped"]["pdf"]
        not_scraped_list = data["upstream_not_scraped"]["pdf"]
        bin_centre_list = data["all_upstream"]["bin_centre_list"]
        bin_width = bin_centre_list[1] - bin_centre_list[0]
        norm = 1. #float(sum(data["all_upstream"]["pdf"]))
        upstream_scraped_graph = ROOT.TGraphAsymmErrors(len(bin_centre_list))
        upstream_not_scraped_graph = ROOT.TGraphAsymmErrors(len(bin_centre_list))
        sys_errors_us = self.sys_errors("all_upstream", pdf_list_tku)
        for i, bin_centre in enumerate(bin_centre_list):
            upstream_scraped_graph.SetPoint(i, bin_centre, scraped_list[i]/norm)
            upstream_not_scraped_graph.SetPoint(i, bin_centre, not_scraped_list[i]/norm)
            for graph in upstream_scraped_graph, upstream_not_scraped_graph:
                graph.SetPointError(i, bin_width/2., bin_width/2., sys_errors_us[i]/norm, sys_errors_us[i]/norm)
        upstream_scraped_graph.SetName("upstream scraped")
        upstream_not_scraped_graph.SetName("upstream not scraped")
        self.root_objects.append(upstream_scraped_graph)
        self.root_objects.append(upstream_not_scraped_graph)
        return upstream_scraped_graph, upstream_not_scraped_graph

    def get_graphs(self, suffix):
        data = self.amplitudes[suffix]
        pdf_list_tku = data["all_upstream"]["pdf"]
        pdf_list_tkd = data["all_downstream"]["pdf"]
        migration_matrix = data["migration_matrix"]
        bin_centre_list = self.amplitudes[suffix]["all_upstream"]["bin_centre_list"]
        norm = 1. #float(sum(pdf_list_tku))
        upstream_graph = ROOT.TGraphAsymmErrors(len(bin_centre_list))
        downstream_graph = ROOT.TGraphAsymmErrors(len(bin_centre_list))
        bin_width = bin_centre_list[1] - bin_centre_list[0]
        stats_errors = self.stats_errors(migration_matrix)
        sys_errors_us = self.sys_errors("all_upstream", pdf_list_tku)
        sys_errors_ds = self.sys_errors("all_downstream", pdf_list_tkd)
        sys.stdout.flush()
        y_max_list = []
        for i, bin_centre in enumerate(bin_centre_list):
            upstream_graph.SetPoint(i, bin_centre, pdf_list_tku[i]/norm)
            downstream_graph.SetPoint(i, bin_centre, pdf_list_tkd[i]/norm)
            upstream_graph.SetPointError(i, bin_width/2., bin_width/2., sys_errors_us[i]/norm, sys_errors_us[i]/norm)
            downstream_graph.SetPointError(i, bin_width/2., bin_width/2.,
                                              (stats_errors[i]**2+sys_errors_ds[i]**2)**0.5/norm,
                                              (stats_errors[i]**2+sys_errors_ds[i]**2)**0.5/norm)

            y_max_list.append(pdf_list_tku[i]/norm)
            y_max_list.append((pdf_list_tkd[i] + stats_errors[i])/norm)
        xmax = max(self.bin_edge_list)
        ymax = max(y_max_list)*1.1
        upstream_graph.SetName("upstream")
        downstream_graph.SetName("downstream")
        hist = common.make_root_histogram("amplitude pdf", [0.], "Amplitude [mm]", 100, [0.], "", 100,
                                          xmin=0., xmax=xmax,
                                          ymin=0., ymax=ymax)
        self.root_objects.append(upstream_graph)
        self.root_objects.append(downstream_graph)
        print "Histogram with x/y max:", xmax, ymax
        return hist, upstream_graph, downstream_graph

    def chi2_graph(self, suffix, distribution):
        data = self.amplitudes[suffix]
        pdf_list_tku = data[distribution]["pdf"]
        bin_centre_list = data[distribution]["bin_centre_list"]
        n_bins = len(bin_centre_list)
        weighted_integral = [bin_centre_list[i]*pdf_list_tku[i] for i in range(n_bins)]
        integral = sum(pdf_list_tku)#*(bin_centre_list[1]-bin_centre_list[0])
        max_bin = max(pdf_list_tku)
        emittance = sum(weighted_integral)/integral/4.
        dummy, graph = scripts.chi2_distribution.chi2_graph(emittance, max_bin, 100, 0., 100.)
        return graph

    def text_box(self, upstream_graph, downstream_graph, scraped_graph):
        legend = ROOT.TLegend(0.6, 0.65, 0.89, 0.89)
        for graph in [upstream_graph, downstream_graph, scraped_graph]:
            legend.AddEntry(graph, graph.GetName(), "lep")
        legend.SetBorderSize(0)
        legend.Draw()
        self.root_objects.append(legend)
        text_box = ROOT.TPaveText(0.6, 0.55, 0.89, 0.65, "NDC")
        text_box.SetFillColor(0)
        text_box.SetBorderSize(0)
        text_box.SetTextSize(0.04)
        text_box.SetTextAlign(12)
        text_box.AddText("MICE Preliminary")
        text_box.Draw()
        self.root_objects.append(text_box)
        text_box = ROOT.TPaveText(0.6, 0.45, 0.89, 0.55, "NDC")
        text_box.SetFillColor(0)
        text_box.SetBorderSize(0)
        text_box.SetTextSize(0.03)
        text_box.SetTextAlign(12)
        text_box.AddText("ISIS Cycle 2016/04")
        text_box.AddText("MAUS v2.8.5")
        text_box.Draw()
        self.root_objects.append(text_box)

    def delta_amplitude_plot_2(self, suffix):
        canvas = common.make_root_canvas("amplitude pdf")
        hist, upstream_graph, downstream_graph = self.get_graphs(suffix)
        upstream_graph.SetMarkerColor(ROOT.kBlue)
        downstream_graph.SetMarkerColor(ROOT.kRed)
        upstream_graph.SetMarkerStyle(24)
        downstream_graph.SetMarkerStyle(26)
        hist.Draw()
        hist.GetXaxis().SetRangeUser(0, self.config.amplitude_max)
        upstream_graph.Draw("SAMEP")
        downstream_graph.Draw("SAMEP")


        upstream_scraped_graph, upstream_not_scraped_graph = self.get_delta_graphs(suffix)

        upstream_scraped_graph.SetMarkerColor(ROOT.kViolet)
        upstream_scraped_graph.SetMarkerStyle(25)
        upstream_scraped_graph.Draw("SAMEP")

        self.text_box(upstream_graph, downstream_graph, upstream_scraped_graph)
        canvas.Update()
        for format in ["eps", "png", "root"]:
            canvas.Print(self.config_anal["plot_dir"]+"amplitude_pdf_"+suffix+"."+format)
        hist.GetXaxis().SetRangeUser(0, 100.)
        upstream_not_scraped_graph.SetMarkerColor(ROOT.kGray+2)
        upstream_not_scraped_graph.SetMarkerStyle(25)
        upstream_not_scraped_graph.Draw("SAMEP")

        upstream_chi2 = self.chi2_graph(suffix, "all_upstream")
        upstream_chi2.SetLineColor(ROOT.kBlue+2) 
        upstream_chi2.Draw("SAMEL")
        downstream_chi2 = self.chi2_graph(suffix, "all_downstream")
        downstream_chi2.SetLineColor(ROOT.kRed+2) 
        downstream_chi2.Draw("SAMEL")

        canvas.Update()
        for format in ["eps", "png", "root"]:
            canvas.Print(self.config_anal["plot_dir"]+"amplitude_pdf_"+suffix+"_extras."+format)

        return canvas

    def delta_amplitude_plot(self, suffix):
        data = self.amplitudes[suffix]
        pdf_list_tku = data["all_upstream"]["pdf"]
        pdf_list_tkd = data["all_downstream"]["pdf"]
        normalisation = float(sum(pdf_list_tku))
        xmax = data["all_upstream"]["bin_edge_list"][-1]
        bin_centre_list = data["all_upstream"]["bin_centre_list"]

        canvas_pdf = common.make_root_canvas("amplitudes")

        # Set up the histograms
        print xmax
        hist = common.make_root_histogram("amplitudes", [-1.], "A_{#perp}   [mm]", 100, xmin = 0., xmax = xmax, ymin = -0.1, ymax=1.1)
        hist.Draw()

        graph_3 = self.amplitude_pdf_graph(bin_centre_list, pdf_list_tku, normalisation, "Upstream", 0.682)
        graph_3.SetMarkerColor(4)
        graph_3.SetMarkerStyle(24)
        graph_1 = self.amplitude_pdf_graph(bin_centre_list, pdf_list_tkd, normalisation, "Downstream", 0.682)
        graph_1.SetMarkerStyle(25)
        graph_1.SetMarkerColor(8)

        legend = ROOT.TLegend(0.6, 0.7, 0.89, 0.89)
        for graph in [graph_3, graph_1]:
            legend.AddEntry(graph, graph.GetName(), "lep")
        legend.Draw()
        self.root_objects.append(legend)
        canvas_pdf.Update()

        plot_dir = self.config_anal["plot_dir"]
        for format in ["eps", "png", "root"]:
            canvas_pdf.Print(plot_dir+"amplitude_pdf_"+suffix+"."+format)
        return canvas_pdf

    #['eventNumber', 'event_number', 'particleNumber', 'particle_number', 'pid', 'spill', 'station', 'status', '', "x'", "y'", 'p', "t'", 'kinetic_energy', 'ct']

    def tk_truth_is_okay(self, tk_mc_dict):
        radius_okay = [key for key, value in tk_mc_dict.iteritems() if value['r'] < 150.]     
        pid_okay = [key for key, value in tk_mc_dict.iteritems() if value['pid'] == -13]
        return [len(tk_mc_dict), len(radius_okay), len(pid_okay)]

    def make_bunches(self, data_loader):
        hits_reco_us = []
        hits_reco_ds = []
        hits_all_mc_us = []
        hits_all_mc_ds = []
        hits_reco_mc_us = []
        hits_reco_mc_ds = []
        if self.config_anal["do_mc"]:
            station_us = ["mc_virtual_"+str(station) for station in self.config.mc_plots["mc_stations"]["tku"]]
            station_ds = ["mc_virtual_"+str(station) for station in self.config.mc_plots["mc_stations"]["tkd"]]

        for event in data_loader.events:
            if 'any_cut' not in event:
                print "Didnt find any_cut in", event.keys()
                continue

            if event['any_cut']:
                continue
            spill = event['tku']['spill']
            ev = event['tku']['event_number']

            mc_hit_us = None
            mc_hit_ds = None
            mc_hit_dict_us = {}
            mc_hit_dict_ds = {}
            if self.config_anal["do_mc"]:
                for detector_hit in event["data"]:
                    if detector_hit["detector"] in station_us:
                        mc_hit_dict_us[detector_hit["detector"]] = detector_hit["hit"]
                    elif detector_hit["detector"] in station_ds:
                        mc_hit_dict_ds[detector_hit["detector"]] = detector_hit["hit"]
                    else:
                        continue
                tku_truth = self.tk_truth_is_okay(mc_hit_dict_us)
                tkd_truth = self.tk_truth_is_okay(mc_hit_dict_ds)
                if tku_truth == [15, 15, 15]:
                    mc_hit_us = mc_hit_dict_us[station_us[0]]
                    mc_hit_us['spill'] = spill
                    hits_all_mc_us.append(mc_hit_us)
                if tkd_truth == [15, 15, 15]:
                    mc_hit_ds = mc_hit_dict_ds[station_ds[0]]
                    mc_hit_ds['spill'] = spill
                    mc_hit_ds['event_number'] = ev
                    hits_all_mc_ds.append(mc_hit_ds)
            if event['tku'] == None:
                continue
            hits_reco_us.append(event['tku'])
            if mc_hit_us != None:
                hits_reco_mc_us.append(mc_hit_us)
            else:
                print "tku impurity", tku_truth

            if event['downstream_cut'] or event['tkd'] == None:
                continue
            hits_reco_ds.append(event['tkd'])
            if mc_hit_ds != None:
                hits_reco_mc_ds.append(mc_hit_ds)
            else:
                print "tkd impurity", tkd_truth
                if tkd_truth[2] != 15:
                    print mc_hit_dict_ds[station_ds[0]]
        print "Loaded upstream:"
        print "    reco:", len(hits_reco_us)
        print "    all mc:", len(hits_all_mc_us)
        print "    reco mc:", len(hits_reco_mc_us)
        print "Loaded downstream:"
        print "    reco:", len(hits_reco_ds)
        print "    all mc:", len(hits_all_mc_ds)
        print "    reco mc:", len(hits_reco_mc_ds)
        reco_us = Bunch.new_from_hits(hits_reco_us)
        reco_ds = Bunch.new_from_hits(hits_reco_ds)
        all_mc_us = Bunch.new_from_hits(hits_all_mc_us)
        all_mc_ds = Bunch.new_from_hits(hits_all_mc_ds)
        reco_mc_us = Bunch.new_from_hits(hits_reco_mc_us)
        reco_mc_ds = Bunch.new_from_hits(hits_reco_mc_ds)

        return reco_us, reco_ds, all_mc_us, all_mc_ds, reco_mc_us, reco_mc_ds

    def print_data(self):
        fout = open(self.config_anal["plot_dir"]+"/amplitude.json", "w")
        out_str = json.dumps(self.amplitudes, indent=2)
        fout.write(out_str)
        fout.close()

    def load_errors(self):
        if self.config_anal["amplitude_source"] != None and self.config_anal["do_mc"]:
            raise RuntimeError("Conflicted - do I get errors from mc or from amplitude source?")
        if self.config_anal["amplitude_source"] == None:
            return 
        fin = open(self.config_anal["amplitude_source"])
        amp_str = fin.read()
        self.amplitudes = json.loads(amp_str)
        # we only want to keep the errors; not the actual pdfs
        del self.amplitudes["reco"]
        del self.amplitudes["all_mc"]
        del self.amplitudes["reco_mc"]
        try:
            del self.amplitudes["field_uncertainty"]
        except KeyError:
            pass
        print self.amplitudes.keys()

    root_objects = []

def do_amplitude_analysis(config, config_anal, data_loader):
    print "Doing amplitude"
    amp_anal = AmplitudeAnalysis(config, config_anal)
    reco_us, reco_ds, all_mc_us, all_mc_ds, reco_mc_us, reco_mc_ds = amp_anal.make_bunches(data_loader)
    amp_anal.load_errors()
    amp_anal.delta_amplitude_calc_2(reco_us, reco_ds, "reco")
    amp_anal.field_uncertainty_calc("all_upstream")
    amp_anal.field_uncertainty_calc("all_downstream")
    if config_anal["do_mc"]:
        amp_anal.delta_amplitude_calc_2(all_mc_us, all_mc_ds, "all_mc")
        amp_anal.delta_amplitude_calc_2(reco_mc_us, reco_mc_ds, "reco_mc")
        amp_anal.systematics_calc("all_upstream")
        amp_anal.systematics_calc("all_downstream")
    for suffix in ["reco"]:
        amp_anal.delta_amplitude_plot_2(suffix)
    amp_anal.field_uncertainty_plot("all_upstream")
    amp_anal.field_uncertainty_plot("all_downstream")
    if config_anal["do_mc"]:
        for suffix in ["reco_mc", "all_mc"]:
            amp_anal.delta_amplitude_plot_2(suffix)
        for target in "all_upstream", "all_downstream":
            amp_anal.systematics_plot(target)
    amp_anal.print_data()

