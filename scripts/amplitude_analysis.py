import sys
import datetime
import bisect
import operator
import json
import numpy
import tempfile

import ROOT
import xboa.common as common
from xboa.hit import Hit
from xboa.bunch import Bunch
import scripts.chi2_distribution
from scripts.amplitude_data import AmplitudeData
from scripts.binomial_confidence_interval import BinomialConfidenceInterval
from analysis_base import AnalysisBase

#REPHRASE SYSTEMATICS;
#* efficiency and purity is the migration matrix from mc to reco mc
#** Better to split efficiency and purity
#* crossing probability is the migration matrix from recon_mc to recon
#* field uncertainty is the migration matrix from recon to "recon with bfield" 

class AmplitudeAnalysis(AnalysisBase):
    def __init__(self, config, config_anal, data_loader):
        super(AmplitudeAnalysis, self).__init__(config, config_anal, data_loader)
        self.config = config
        self.config_anal = config_anal
        self.amplitudes = {"momentum_uncertainty":{}, "inefficiency":{}, "crossing_probability":{}}
        self.bin_edge_list = [float(i) for i in range(0, 101, config.amplitude_bin_width)]
        self.a_dir = tempfile.mkdtemp()
        file_name = self.a_dir+"/amp_data_"
        mu_mass = common.pdg_pid_to_mass[13]
        self.all_mc_data_us = AmplitudeData(file_name+"mc_us", self.bin_edge_list, mu_mass)
        self.reco_mc_data_us = AmplitudeData(file_name+"reco_mc_us", self.bin_edge_list, mu_mass)
        self.reco_data_us = AmplitudeData(file_name+"recon_us", self.bin_edge_list, mu_mass)
        self.all_mc_data_ds = AmplitudeData(file_name+"mc_ds", self.bin_edge_list, mu_mass)
        self.reco_mc_data_ds = AmplitudeData(file_name+"reco_mc_ds", self.bin_edge_list, mu_mass)
        self.reco_data_ds = AmplitudeData(file_name+"recon_ds", self.bin_edge_list, mu_mass)

    def delta_amplitude_calc(self, suffix):
        data_0, data_1 = {
            "all_mc":(self.all_mc_data_us, self.all_mc_data_ds),
            "reco":(self.reco_data_us, self.reco_data_ds),
            "reco_mc":(self.reco_mc_data_us, self.reco_mc_data_ds),
        }[suffix]
        print "Amplitude", suffix
        print "  starting delta_amplitude_calc    ", datetime.datetime.now()
        data = {}
        print "  Upstream amplitude calculation...  ", datetime.datetime.now()
        fractional_amplitude_dict_0 = data_0.fractional_amplitude()
        print "  n_events", len(fractional_amplitude_dict_0)
        print "  Downstream amplitude calculation...", datetime.datetime.now()
        fractional_amplitude_dict_1 = data_1.fractional_amplitude()
        print "  n_events", len(fractional_amplitude_dict_1)
        print "  ... done                           ", datetime.datetime.now()

        data["amplitude_dict_upstream"] = fractional_amplitude_dict_0
        data["amplitude_dict_downstream"] = fractional_amplitude_dict_1
        data["all_upstream"] = self.get_pdfs(fractional_amplitude_dict_0)
        data["all_downstream"] = self.get_pdfs(fractional_amplitude_dict_1)
        for name, amp_data in [("all_upstream", data_0),
                           ("all_downstream", data_1)]:
            data[name]["emittance"] = amp_data.get_emittance()
            data[name]["weight"] = sum(amp_data.n_events)
            data[name]["covariance_matrix"] = amp_data.cov.tolist()
            data[name]["bin_edge_list"] = self.bin_edge_list
        print "  upstream  ", data["all_upstream"]["pdf"], "sum", sum(data["all_upstream"]["pdf"])
        print "  downstream", data["all_downstream"]["pdf"], "sum", sum(data["all_downstream"]["pdf"])

        data["migration_matrix"] = self.migration_matrix(fractional_amplitude_dict_0, fractional_amplitude_dict_1, False)
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
        data["pdf"] = [0]*(len(self.bin_edge_list))
        amplitude_index = 0
        for amplitude in fractional_amplitude_dict.values():
            amp_0_bin = bisect.bisect_left(self.bin_edge_list, amplitude)-1
            data["pdf"][amp_0_bin] += 1
        return data

    def migration_matrix(self, fractional_amplitude_dict_0, fractional_amplitude_dict_1, normalise):
        n_bins = len(self.bin_edge_list)
        migration_matrix = [[0. for i in range(n_bins)]  for j in range(n_bins)]
        # migration dict maps event number onto amplitude tuple (us, ds)
        migration_dict = {}
        for key, amp_0 in fractional_amplitude_dict_0.iteritems():
            amp_1 = None
            if key in fractional_amplitude_dict_1:
                amp_1 = fractional_amplitude_dict_1[key]
            migration_dict[key] = (amp_0, amp_1)
        for key, amp_pair in migration_dict.iteritems():
            if amp_pair[0] == None or amp_pair[1] == None:
                continue # could do an overflow bin (e.g. events that scraped between upstream and downstream)
            amp_0_bin = bisect.bisect_left(self.bin_edge_list, amp_pair[0])-1
            amp_1_bin = bisect.bisect_left(self.bin_edge_list, amp_pair[1])-1
            migration_matrix[amp_0_bin][amp_1_bin] += 1
        if normalise:
            for i, row in enumerate(migration_matrix):
                row_sum = float(sum(row))
                for j, element in enumerate(row):
                    if abs(row_sum) < 1e-3:
                        if i == j:
                            migration_matrix[i][j] = 1. # identity
                        else:
                            migration_matrix[i][j] = 0.
                    else:
                        migration_matrix[i][j] /= row_sum
        return migration_matrix

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

    def systematics_calc(self, target):
        amp_key = {"all_upstream":"amplitude_dict_upstream",
                   "all_downstream":"amplitude_dict_downstream"}[target]
        all_mc_pdf = self.amplitudes["all_mc"][target]["pdf"]
        reco_mc_pdf = self.amplitudes["reco_mc"][target]["pdf"]
        reco_pdf = self.amplitudes["reco"][target]["pdf"]
        print "reco pdf", target, reco_pdf

        fractional_mc = self.amplitudes["all_mc"][amp_key]
        fractional_reco_mc = self.amplitudes["reco_mc"][amp_key]
        fractional_reco = self.amplitudes["reco"][amp_key]
        print "mc pdf", target, all_mc_pdf
        print "mc reco pdf", target, reco_mc_pdf
        inefficiency = []
        for i in range(len(all_mc_pdf)):
            if reco_mc_pdf[i] == 0:
                inefficiency.append(1.)
            else:
                inefficiency.append(float(all_mc_pdf[i])/reco_mc_pdf[i])
        self.amplitudes["inefficiency"][target] = {"pdf_ratio":inefficiency}
        self.amplitudes["crossing_probability"][target] = {}
        self.amplitudes["crossing_probability"][target]["migration_matrix"] = \
                      self.migration_matrix(fractional_reco, fractional_reco_mc, True)

    def field_uncertainty_calc(self, target):
        amp_key = {"all_upstream":"amplitude_dict_upstream",
                   "all_downstream":"amplitude_dict_downstream"}[target]
        fractional_reco = self.amplitudes["reco"][amp_key]
        fractional_reco_uncertainty = {}
        for key, value in fractional_reco.iteritems():
            #print "field uncertainty calc", key, value
            if value == None:
                fractional_reco_uncertainty[key] = None
            else:
                fractional_reco_uncertainty[key] = value*(1.+self.config_anal["field_uncertainty"])
        self.amplitudes["momentum_uncertainty"][target] = {}
        self.amplitudes["momentum_uncertainty"][target]["migration_matrix"] = \
             self.migration_matrix(fractional_reco, fractional_reco_uncertainty, True)

    def systematics_plot(self, target):
        for key in ["momentum_uncertainty", "crossing_probability"]:
            matrix = self.amplitudes[key][target]["migration_matrix"]
            title = (key+"_"+target).replace("_all", "").replace("_", " ")
            self.matrix_plot(matrix, title)
        inefficiency = [1-x for x in self.amplitudes["inefficiency"][target]["pdf_ratio"]]
        canvas = common.make_root_canvas("inefficiency_"+target)
        print "Plotting", len(self.get_bin_centre_list()), len(inefficiency)
        hist, graph = common.make_root_graph("inefficiency",
                                             self.get_bin_centre_list(), "Amplitude [mm]",
                                             inefficiency[:-1], "Inefficiency [fractional]")
        hist.SetTitle(self.config_anal['name'])
        hist.Draw()
        graph.SetMarkerStyle(24)
        graph.Draw("p")
        for format in "eps", "root", "png":
            canvas.Print(self.config_anal["plot_dir"]+"amplitude_inefficiency_"+target+"."+format)
        hist.GetYaxis().SetRangeUser(-0.1, 0.1)
        hist.Draw()
        graph.SetMarkerStyle(24)
        graph.Draw("p")
        for format in "eps", "root", "png":
            canvas.Print(self.config_anal["plot_dir"]+"amplitude_inefficiency_zoom_"+target+"."+format)

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
        self.matrix_plot(migration_matrix, "migration matrix")
        downstream_errors = [sum(row)**0.5 for row in error_matrix]
        error_matrix_sqrt = [[element**0.5 for element in row] for row in error_matrix]
        self.matrix_plot(error_matrix_sqrt, "stats error")
        print "upstream pdf\n", self.matrix_str(row_sum)
        print "migration matrix\n", self.matrix_str(migration_matrix)
        print "error matrix\n", self.matrix_str(error_matrix)
        print "downstream errors\n", self.matrix_str(downstream_errors)
        return downstream_errors

    def matrix_plot(self, matrix, title, normalise=False):
        bin_centre_list = self.amplitudes["reco"]["all_upstream"]["bin_centre_list"]
        xmax = max(self.bin_edge_list)
        nbins = len(bin_centre_list)
        matrix_hist = ROOT.TH2D(title, title+";Amplitude Downstream [mm];Amplitude Upstream [mm]", nbins, 0., xmax, nbins, 0., xmax)
        for i in range(nbins):
            for j in range(nbins):
                matrix_hist.Fill(bin_centre_list[i], bin_centre_list[j], matrix[i][j])
        self.root_objects.append(matrix_hist)
        canvas = common.make_root_canvas(title)
        matrix_hist.SetTitle(self.config_anal['name'])
        matrix_hist.SetStats(False)
        matrix_hist.Draw("COLZ")
        title = title.replace(" ", "_")
        for format in ["eps", "png", "root"]:
            canvas.Print(self.config_anal["plot_dir"]+"amplitude_"+title+"."+format)

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
        sys_errors_us = [0. for i in range(21)] #self.sys_errors("all_upstream", pdf_list_tku)
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

    def do_migrations(self, target):
        print "Doing migrations for", target
        reco_pdfs = numpy.array(self.amplitudes["reco"][target]["pdf"])
        print "Reco             ", reco_pdfs
        reco_mc_matrix = self.amplitudes["crossing_probability"][target]["migration_matrix"]
        reco_mc_matrix = numpy.transpose(numpy.array(reco_mc_matrix))
        mc_reco_pdfs = numpy.dot(reco_mc_matrix, reco_pdfs)
        print "Recovered MC Reco", mc_reco_pdfs
        try:
            print "Actual MC Reco   ", numpy.array(self.amplitudes["reco_mc"][target]["pdf"])
        except KeyError:
            print "<Nothing recorded>"
        mc_pdfs = mc_reco_pdfs*numpy.array(self.amplitudes["inefficiency"][target]["pdf_ratio"])
        print "Recovered All MC ", mc_pdfs
        try:
            print "Actual all MC    ", numpy.array(self.amplitudes["all_mc"][target]["pdf"])
        except KeyError:
            print "<Nothing recorded>"
        sys_errors = [0. for i in reco_pdfs]
        return mc_pdfs.tolist(), sys_errors

    def get_graphs(self, suffix):
        data = self.amplitudes[suffix]
        raw_pdf_list_tku = data["all_upstream"]["pdf"]
        raw_pdf_list_tkd = data["all_downstream"]["pdf"]
        migration_matrix = data["migration_matrix"]
        bin_centre_list = self.amplitudes[suffix]["all_upstream"]["bin_centre_list"]
        norm = 1. #float(sum(pdf_list_tku))
        upstream_graph = ROOT.TGraphAsymmErrors(len(bin_centre_list))
        downstream_graph = ROOT.TGraphAsymmErrors(len(bin_centre_list))
        raw_upstream_graph = ROOT.TGraphAsymmErrors(len(bin_centre_list))
        raw_downstream_graph = ROOT.TGraphAsymmErrors(len(bin_centre_list))
        bin_width = bin_centre_list[1] - bin_centre_list[0]
        stats_errors = self.stats_errors(migration_matrix)
        if suffix == "reco":
            pdf_list_tku, sys_errors_us = self.do_migrations("all_upstream")
            pdf_list_tkd, sys_errors_ds = self.do_migrations("all_downstream")
        else:
            pdf_list_tku, sys_errors_us = raw_pdf_list_tku, [0. for i in raw_pdf_list_tku]
            pdf_list_tkd, sys_errors_ds = raw_pdf_list_tkd, [0. for i in raw_pdf_list_tkd]
        sys.stdout.flush()
        y_max_list = []
        for i, bin_centre in enumerate(bin_centre_list):
            raw_upstream_graph.SetPoint(i, bin_centre, raw_pdf_list_tku[i]/norm)
            raw_downstream_graph.SetPoint(i, bin_centre, raw_pdf_list_tkd[i]/norm)
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
        raw_upstream_graph.SetName("raw upstream")
        raw_downstream_graph.SetName("raw downstream")
        hist = common.make_root_histogram("amplitude pdf", [0.], "Amplitude [mm]", 100, [0.], "", 100,
                                          xmin=0., xmax=xmax,
                                          ymin=0., ymax=ymax)
        self.root_objects.append(upstream_graph)
        self.root_objects.append(downstream_graph)
        self.root_objects.append(raw_upstream_graph)
        self.root_objects.append(raw_downstream_graph)
        print "Histogram with x/y max:", xmax, ymax
        return hist, upstream_graph, downstream_graph, raw_upstream_graph, raw_downstream_graph

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

    def text_box(self, graph_list):
        legend = ROOT.TLegend(0.6, 0.65, 0.89, 0.89)
        for graph in graph_list:
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
        text_box.AddText(self.config_anal['name'])
        text_box.Draw()
        self.root_objects.append(text_box)

    def delta_amplitude_plot(self, suffix):
        canvas = common.make_root_canvas("amplitude pdf")
        hist, upstream_graph, downstream_graph, raw_upstream_graph, raw_downstream_graph = self.get_graphs(suffix)
        upstream_graph.SetMarkerColor(ROOT.kBlue)
        downstream_graph.SetMarkerColor(ROOT.kRed)
        upstream_graph.SetMarkerStyle(24)
        downstream_graph.SetMarkerStyle(26)
        canvas.cd()
        hist.Draw()
        hist.GetXaxis().SetRangeUser(0, self.config.amplitude_max)
        upstream_graph.Draw("SAMEP")
        downstream_graph.Draw("SAMEP")

        upstream_scraped_graph, upstream_not_scraped_graph = self.get_delta_graphs(suffix)

        upstream_scraped_graph.SetMarkerColor(ROOT.kViolet)
        upstream_scraped_graph.SetMarkerStyle(25)
        upstream_scraped_graph.Draw("SAMEP")

        self.text_box([upstream_graph, downstream_graph, upstream_scraped_graph])
        canvas.Update()
        for format in ["eps", "png", "root"]:
            canvas.Print(self.config_anal["plot_dir"]+"amplitude_pdf_"+suffix+"."+format)

        hist.GetXaxis().SetRangeUser(0, 100.)
        upstream_not_scraped_graph.SetMarkerColor(ROOT.kGray+2)
        upstream_not_scraped_graph.SetMarkerStyle(25)
        upstream_not_scraped_graph.Draw("SAMEP")

        raw_upstream_graph.SetMarkerStyle(20)
        raw_upstream_graph.SetMarkerColor(ROOT.kBlue)
        raw_upstream_graph.Draw("SAMEP")

        raw_downstream_graph.SetMarkerStyle(22)
        raw_downstream_graph.SetMarkerColor(ROOT.kRed)
        raw_downstream_graph.Draw("SAMEP")

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

    def amplitude_scatter_plot(self, suffix):
        amp_dict_ds = self.amplitudes[suffix]["amplitude_dict_downstream"]
        amp_dict_us = self.amplitudes[suffix]["amplitude_dict_upstream"]
        amp_list_us = []
        amp_list_ds = []
        amp_list_delta = []
        for key, amp_ds in amp_dict_ds.iteritems():
            try:
                amp_us = amp_dict_us[key]
            except KeyError:
                sys.excepthook(*sys.exc_info())
                continue
            amp_delta = amp_us-amp_ds
            amp_list_us.append(amp_us)
            amp_list_ds.append(amp_ds)
            amp_list_delta.append(amp_delta)
        canvas = common.make_root_canvas("delta_amplitude_scatter")
        canvas.Draw()
        n_points = min(len(amp_list_us), 10000) # no more than 10k points in the scatter
        hist, graph = common.make_root_graph("delta amplitude scatter",
                                             amp_list_us, "Amplitude (us) [mm]",
                                             amp_list_delta, "Amplitude (us) - Amplitude (ds) [mm]", xmin=0., xmax=100., ymin=-50., ymax=50.)
        hist.SetTitle(self.config_anal['name'])
        hist.Draw()
        graph.Draw("P")
        canvas.Update()
        for format in ["eps", "png", "root"]:
            canvas.Print(self.config_anal["plot_dir"]+"amplitude_delta_"+suffix+"_scatter."+format)

        canvas = common.make_root_canvas("delta_amplitude_hist")
        canvas.Draw()
        hist = common.make_root_histogram("delta amplitude hist",
                                          amp_list_us, "Amplitude (us) [mm]", 100,
                                          amp_list_delta, "Amplitude (us) - Amplitude (ds) [mm]", 100,
                                          xmin=0., xmax=100., ymin=-50., ymax=50.)
        hist.SetTitle(self.config_anal['name'])
        hist.Draw("COLZ")
        canvas.Update()
        for format in ["eps", "png", "root"]:
            canvas.Print(self.config_anal["plot_dir"]+"amplitude_delta_"+suffix+"_hist."+format)
      

    def tk_truth_is_okay(self, tk_mc_dict):
        radius_okay = [key for key, value in tk_mc_dict.iteritems() if value['r'] < 150.]     
        pid_okay = [key for key, value in tk_mc_dict.iteritems() if value['pid'] == -13]
        return [len(tk_mc_dict), len(radius_okay), len(pid_okay)]

    def append_data(self):
        hits_reco_us = []
        hits_all_mc_us = []
        hits_reco_mc_us = []
        hits_reco_ds = []
        hits_all_mc_ds = []
        hits_reco_mc_ds= []

        if self.config_anal["do_mc"]:
            station_us = ["mc_virtual_"+str(station) for station in self.config.mc_plots["mc_stations"]["tku"]]
            station_ds = ["mc_virtual_"+str(station) for station in self.config.mc_plots["mc_stations"]["tkd"]]

        for event in self.data_loader.events:
            if 'any_cut' not in event:
                print "Didnt find any_cut in", event.keys()
                continue

            if event['upstream_cut']:
                continue

            mc_hit_us = None
            mc_hit_ds = None
            mc_hit_dict_us = {}
            mc_hit_dict_ds = {}
            tku_truth = [None]*3
            tkd_truth = [None]*3
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
                    hits_all_mc_us.append(mc_hit_us)
                if tkd_truth == [15, 15, 15]:
                    mc_hit_ds = mc_hit_dict_ds[station_ds[0]]
                    hits_all_mc_ds.append(mc_hit_ds)
            hits_reco_us.append(event['tku'])
            if mc_hit_us != None:
                hits_reco_mc_us.append(mc_hit_us)
            elif tku_truth[0] != None:
                pass #print "tku impurity", tku_truth

            if event['downstream_cut']:
                continue
            #if event['tkd'] == None:
            #    print event['tkd']
            #    print [hit for hit in event['data'] if 'tkd' in hit['detector']]
            #    print 'cuts', event['will_cut']
            hits_reco_ds.append(event['tkd'])
            if mc_hit_ds != None:
                hits_reco_mc_ds.append(mc_hit_ds)
            elif tkd_truth[0] != None:
                #print "tkd impurity", tkd_truth
                if tkd_truth[2] != 15 and tkd_truth[2] != None:
                    print
                    try:
                        print mc_hit_dict_ds[station_ds[0]]
                    except KeyError:
                        print mc_hit_dict_ds
        self.all_mc_data_us.append_hits(hits_all_mc_us)
        self.reco_mc_data_us.append_hits(hits_reco_mc_us)
        self.reco_data_us.append_hits(hits_reco_us)
        self.all_mc_data_ds.append_hits(hits_all_mc_ds)
        self.reco_mc_data_ds.append_hits(hits_reco_mc_ds)
        self.reco_data_ds.append_hits(hits_reco_ds)
        return
        print "Loaded upstream:"
        print "    reco:   ", len(hits_reco_us)
        print "    all mc: ", len(hits_all_mc_us)
        print "    reco mc:", len(hits_reco_mc_us)
        print "Loaded downstream:"
        print "    reco:   ", len(hits_reco_ds)
        print "    all mc: ", len(hits_all_mc_ds)
        print "    reco mc:", len(hits_reco_mc_ds)

    def print_data(self):
        fout = open(self.config_anal["plot_dir"]+"/amplitude.json", "w")
        for suffix in self.amplitudes:
            print self.amplitudes[suffix].keys()
            try:
                del self.amplitudes[suffix]["amplitude_dict_upstream"]
                del self.amplitudes[suffix]["amplitude_dict_downstream"]
            except KeyError:
                pass
        out_str = json.dumps(self.amplitudes, indent=2)
        fout.write(out_str)
        fout.close()

    def load_errors(self):
        if self.config_anal["amplitude_source"] != None and self.config_anal["do_mc"]:
            raise RuntimeError("Conflicted - do I get errors from mc or from amplitude source?")
        if self.config_anal["amplitude_source"] == None and not self.config_anal["do_mc"]:
            raise RuntimeError("No source for errors was selected - abort")
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

    def birth(self):
        self.load_errors()
        self.append_data()

    def process(self):
        self.append_data()

    def death(self):
        self.delta_amplitude_calc("reco")
        self.field_uncertainty_calc("all_upstream")
        self.field_uncertainty_calc("all_downstream")
        if self.config_anal["do_mc"]:
            self.delta_amplitude_calc("all_mc")
            self.delta_amplitude_calc("reco_mc")
            self.systematics_calc("all_upstream")
            self.systematics_calc("all_downstream")
        if self.config_anal["do_mc"]:
            for target in "all_upstream", "all_downstream":
                self.systematics_plot(target)
            for suffix in ["reco_mc", "all_mc"]:
                self.delta_amplitude_plot(suffix)
                self.amplitude_scatter_plot(suffix)
        for suffix in ["reco"]:
            self.delta_amplitude_plot(suffix)
            self.amplitude_scatter_plot(suffix)
        self.print_data()

