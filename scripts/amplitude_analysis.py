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
from scripts.amplitude_data_binned import AmplitudeDataBinned
from scripts.amplitude_data_binless import AmplitudeDataBinless
from scripts.binomial_confidence_interval import BinomialConfidenceInterval
from analysis_base import AnalysisBase
from scripts.analysis.amplitude.plot_amplitude_data import PlotAmplitudeData

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
        AmplitudeData = self.get_amplitude_algorithm()
        self.all_mc_data_us = AmplitudeData(file_name+"mc_us", self.bin_edge_list, mu_mass, self.cov_fixed_us)
        self.reco_mc_data_us = AmplitudeData(file_name+"reco_mc_us", self.bin_edge_list, mu_mass, self.cov_fixed_us)
        self.reco_data_us = AmplitudeData(file_name+"recon_us", self.bin_edge_list, mu_mass, self.cov_fixed_us)
        self.all_mc_data_ds = AmplitudeData(file_name+"mc_ds", self.bin_edge_list, mu_mass, self.cov_fixed_ds)
        self.reco_mc_data_ds = AmplitudeData(file_name+"reco_mc_ds", self.bin_edge_list, mu_mass, self.cov_fixed_ds)
        self.reco_data_ds = AmplitudeData(file_name+"recon_ds", self.bin_edge_list, mu_mass, self.cov_fixed_ds)

    def get_amplitude_algorithm(self):
        if self.config_anal['amplitude_algorithm'] == 'binned':
            AmplitudeData = AmplitudeDataBinned
            self.cov_fixed_us = None
            self.cov_fixed_ds = None
        elif self.config_anal['amplitude_algorithm'] == 'binless':
            AmplitudeData = AmplitudeDataBinless
            self.cov_fixed_us = None
            self.cov_fixed_ds = None
        elif self.config_anal['amplitude_algorithm'] == 'fixed':
            AmplitudeData = AmplitudeDataBinless
            self.cov_fixed_us = self.config_anal['cov_fixed_us']
            self.cov_fixed_ds = self.config_anal['cov_fixed_ds']
        else:
            raise ValueError("Did not understand amplitude algorithm "+str(config_anal['amplitude_algorithm']))
        return AmplitudeData

    def birth(self):
        self.set_plot_dir("amplitude")
        self.load_errors()
        self.append_data()

    def process(self):
        self.append_data()

    def death(self):
        self.delta_amplitude_calc("reco")
        self.field_uncertainty_calc("all_upstream")
        self.field_uncertainty_calc("all_downstream")
        if self.config_anal["amplitude_mc"]:
            self.delta_amplitude_calc("all_mc")
            self.delta_amplitude_calc("reco_mc")
            self.systematics_calc("all_upstream")
            self.systematics_calc("all_downstream")
        if self.config_anal["amplitude_mc"]:
            for target in "all_upstream", "all_downstream":
                self.systematics_plot(target)
            for suffix in ["reco_mc", "all_mc"]:
                self.pdf_plot(suffix)
                self.amplitude_scatter_plot(suffix)
        for suffix in ["reco"]:
            self.corrections_and_uncertainties(suffix)
            self.cdf_data(suffix)
            self.ratio_data(suffix)
            self.pdf_plot(suffix)
            self.cdf_plot(suffix)
            self.pdf_ratio_plot(suffix)
            self.cdf_ratio_plot(suffix)
            self.amplitude_scatter_plot(suffix)
        self.print_data()

    def delta_amplitude_calc(self, suffix):
        """
        Calculate amplitudes for upstream and downstream sample
        * suffix: "reco", "all_mc" or "reco_mc" for different sampling criteria,
          as defined during the append_data step. reco is supposed to be
          recontructed data; all_mc is all monte carlo data; reco_mc is all
          monte carlo data that is also in the reco sample
        Calculates amplitudes using the AmplitudeDataBinless algorithms; fills
        dicts; then fills pdfs and migration matrices
        """
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
        PlotAmplitudeData(data_0, self.plot_dir, "us").plot()
        PlotAmplitudeData(data_1, self.plot_dir, "ds").plot()
        print "  done amplitude calc                ", datetime.datetime.now()

    def get_bin_centre_list(self):
        """
        Return the list of bin centres, calculated using self.bin_edge_list
        """
        bin_width = self.config.amplitude_bin_width
        bin_centre_list = [bin_edge+bin_width/2. for bin_edge in self.bin_edge_list[:-1]]
        return bin_centre_list

    def get_delta_pdfs(self, a_dict_us, a_dict_ds):
        """
        Return binned data corresponding to scraped events
        
        * a_dict_us: dictionary mapping <event_id>:amplitude
        * a_dict_ds: dictionary mapping <event_id>:amplitude

        Returns a tuple of two lists (upstream_scraped, upstream_not_scraped)
        * upstream_scraped corresponds to the events that are in a_dict_us but
                            are not in a_dict_ds
        * upstream_not_scraped corresponds to the events that are in a_dict_us
                            and are in a_dict_ds
        """
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
        """
        Bin the data in the amplitude dict and return as a list
        """
        bin_centre_list = self.get_bin_centre_list()
        data = {}
        data["bin_centre_list"] = bin_centre_list
        data["pdf"] = [0]*(len(self.bin_edge_list))
        amplitude_index = 0
        for amplitude in fractional_amplitude_dict.values():
            amp_0_bin = bisect.bisect_left(self.bin_edge_list, amplitude)-1
            data["pdf"][amp_0_bin] += 1
        return data

    def plot_amplitude_data(self, amplitude_data, plot_lambda):
        """
        Plot 1D histogram of the ellipse together with the ellipses used for the
        amplitude calculation
        """
        pass


    def migration_matrix(self, fractional_amplitude_dict_0, fractional_amplitude_dict_1, normalise):
        """
        Generate a migration matrix for the bins
        * fractional_amplitude_dict_0: dictionary of <event_id>:amplitude
        * fractional_amplitude_dict_1: dictionary of <event_id>:amplitude
        <event_id> can be any hashable id; it has to be the same for dict_0 and dict_1
        * normalise: if True, rows are set so that sum(row) = 1. If sum(row) = 0
          before normalisation then the diagonal element is set to 1.

        If normalise is set to False, generate a matrix m_{ij} with elements 
        given by the number of events in sample 0 bin i AND sample 1 bin j.
        
        If normalise is set to true, generate a matrix m_{ij} with elements given
        by the (number of events in sample 0 bin i AND sample 1 bin j)/(number of
        events in sample 0 bin i). Event is not in both sample 0 and sample 1
        are ignored.

        One use is that the probability of crossing from bin i to bin j is
        given by (Number in bin i AND bin j)/(Number in bin i)

        Returns matrix M represented as a list of lists length n_bins by n_bins
        """
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
             self.migration_matrix(fractional_reco_uncertainty, fractional_reco, True)

    def systematics_plot(self, target):
        """
        Plot systematic corrections:
        * Momentum correction matrix
        * Crossing probability due to reconstruction
        * Inefficiency/impurity
        """
        if "upstream" in target:
            label = "US"
        elif "downstream" in target:
            label = "DS"
        matrix = self.amplitudes["momentum_uncertainty"][target]["migration_matrix"]
        title = ("momentum_uncertainty_"+target).replace("_all", "").replace("_", " ")
        self.matrix_plot(matrix, title,
                         label+" "+self.get_suffix_label("reco")+" amplitude before momentum correction [mm]",
                         label+" "+self.get_suffix_label("reco")+" amplitude after momentum correction [mm]")

        matrix = self.amplitudes["crossing_probability"][target]["migration_matrix"]
        title = ("crossing_probability_"+target).replace("_all", "").replace("_", " ")
        self.matrix_plot(matrix, title,
                         label+" "+self.get_suffix_label("reco")+" amplitude [mm]",
                         label+" "+self.get_suffix_label("reco_mc")+" amplitude [mm]")

        inefficiency = self.amplitudes["inefficiency"][target]["pdf_ratio"]
        canvas = common.make_root_canvas("inefficiency_"+target)

        hist, graph = common.make_root_graph("inefficiency",
                                             self.get_bin_centre_list(), label+" MC Truth Amplitude [mm]",
                                             inefficiency[:-1], label+" "+self.get_suffix_label("all_mc")+" number/"+self.get_suffix_label("reco_mc")+" number")
        hist.GetYaxis().SetRangeUser(0., 2.0)
        hist.SetTitle(self.config_anal['name'])
        hist.Draw()
        graph.SetMarkerStyle(24)
        graph.Draw("p")
        for format in "eps", "root", "png":
            canvas.Print(self.plot_dir+"amplitude_inefficiency_"+target+"."+format)
        hist.GetYaxis().SetRangeUser(0.9, 1.1)
        hist.Draw()
        graph.SetMarkerStyle(24)
        graph.Draw("p")
        for format in "eps", "root", "png":
            canvas.Print(self.plot_dir+"amplitude_inefficiency_zoom_"+target+"."+format)

    def matrix_str(self, matrix, rounding=2, width=10):
        """Represent the matrix as a string (space separated)"""
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

    def stats_errors(self, migration_matrix, suffix):
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
        suffix_label = self.get_suffix_label(suffix)
        self.matrix_plot(migration_matrix, "migration matrix",
                         "DS "+suffix_label+" Amplitude [mm]",
                         "US "+suffix_label+" Amplitude [mm]")
        downstream_errors = [sum(row)**0.5 for row in error_matrix]
        error_matrix_sqrt = [[element**0.5 for element in row] for row in error_matrix]
        self.matrix_plot(error_matrix_sqrt, "stats error",
                         "DS "+suffix_label+" Amplitude [mm]",
                         "US "+suffix_label+" Amplitude [mm]")
        print "upstream pdf\n", self.matrix_str(row_sum)
        print "migration matrix\n", self.matrix_str(migration_matrix)
        print "error matrix\n", self.matrix_str(error_matrix)
        print "downstream errors\n", self.matrix_str(downstream_errors)
        return downstream_errors

    def matrix_plot(self, matrix, title, axis_1, axis_2):
        """
        Plot the matrix entries as a TH2D
        
        * matrix: matrix data to plot; list n_bins by n_bins
        * title: the title of the hist (used in the file name also)
        * axis_1: x axis title
        * axis_2: y axis title
        """
        bin_centre_list = self.amplitudes["reco"]["all_upstream"]["bin_centre_list"]
        xmax = max(self.bin_edge_list)
        nbins = len(bin_centre_list)
        matrix_hist = ROOT.TH2D(title, title+";"+axis_1+";"+axis_2, nbins, 0., xmax, nbins, 0., xmax)
        for i in range(nbins):
            for j in range(nbins):
                matrix_hist.Fill(bin_centre_list[i], bin_centre_list[j], matrix[i][j])
        self.root_objects.append(matrix_hist)
        canvas = common.make_root_canvas(title)
        canvas.SetFrameFillColor(scripts.utilities.get_frame_fill())
        matrix_hist.SetTitle(self.config_anal['name'])
        matrix_hist.SetStats(False)
        matrix_hist.Draw("COLZ")
        title = title.replace(" ", "_")
        for format in ["eps", "png", "root"]:
            canvas.Print(self.plot_dir+title+"."+format)

    def do_migrations(self, target, systematic):
        #print "Doing migrations for", target
        if systematic == None:
            migrations = self.amplitudes
        else:
            migrations = self.amplitudes["systematics"][systematic]
        reco_pdfs = numpy.array(self.amplitudes["reco"][target]["pdf"])
        #print "Reco             ", reco_pdfs
        reco_mc_matrix = migrations["crossing_probability"][target]["migration_matrix"]
        reco_mc_matrix = numpy.transpose(numpy.array(reco_mc_matrix))
        mc_reco_pdfs = numpy.dot(reco_mc_matrix, reco_pdfs)
        #print "Recovered MC Reco", mc_reco_pdfs
        #try:
        #    print "Actual MC Reco   ", numpy.array(self.amplitudes["reco_mc"][target]["pdf"])
        #except KeyError:
        #    print "<Nothing recorded>"
        mc_pdfs = mc_reco_pdfs*numpy.array(migrations["inefficiency"][target]["pdf_ratio"])
        #print "Recovered All MC ", mc_pdfs
        #try:
        #    print "Actual all MC    ", numpy.array(self.amplitudes["all_mc"][target]["pdf"])
        #except KeyError:
        #    print "<Nothing recorded>"
        return mc_pdfs.tolist()

    def corrections_and_uncertainties(self, suffix):
        """
        Calculate corrected pdfs and uncertainties
        * corrected pdf is given by bin_i = Sum_j(Eff_i*M_ij*bin_j) where M_ij
            is the migration matrix and Eff_i is the efficiency
        * statistical errors are given by sum in quadrature of binomial errors
            on each bin migration
        * systematic errors are given by estimated momentum uncertainty (CHECK,
            NEEDS WORK)
        * total errors are given by sum in quadrature of statistical errors and
          systematic errors
        """
        data = self.amplitudes[suffix]
        raw_pdf_list_tku = data["all_upstream"]["pdf"]
        raw_pdf_list_tkd = data["all_downstream"]["pdf"]
        bins = range(len(raw_pdf_list_tku))
        migration_matrix = data["migration_matrix"]
        data["all_upstream"]["pdf_stats_errors"] = [0. for bin in raw_pdf_list_tku]
        data["all_downstream"]["pdf_stats_errors"] = self.stats_errors(migration_matrix, suffix)
        for key in ["all_upstream", "all_downstream"]:
            if suffix == "reco":
                pdf_list = self.do_migrations(key, None)
            else:
                pdf_list = raw_pdf_list_tku
            sys_error_list = [0 for i in bins]
            if len(self.amplitudes["systematics"]) > 0:
                print "Systematic errors"
                ref_pdf_list = self.do_migrations(key, 0)
                for i in range(len(self.amplitudes["systematics"][1:])):
                    i += 1
                    scale = self.amplitudes["systematics"][i]["scale"]
                    print " ", i, self.amplitudes["systematics"][i]["source"], "with scale", scale
                    sys_pdf_list = self.do_migrations(key, i)
                    print "    err: ",
                    for j in bins:
                        err = (sys_pdf_list[j] - ref_pdf_list[j])*scale
                        print round(err),
                        sys_error_list[j] = (sys_error_list[j]**2+err**2)**0.5
                    print "\n    sum: ", [round(err) for err in sys_error_list]
            data[key]["corrected_pdf"] = pdf_list
            data[key]["pdf_sys_errors"] = sys_error_list
            data[key]["pdf_tot_errors"] = [0. for i in bins]
            for i in bins:
                sys = data[key]["pdf_sys_errors"][i]
                stats = data[key]["pdf_stats_errors"][i]
                data[key]["pdf_tot_errors"][i] = (sys**2+stats**2)**0.5
            print "    sys errors:   ", data[key]["pdf_sys_errors"] 
            print "    stats errors: ", data[key]["pdf_stats_errors"] 
            print "    total errors: ", data[key]["pdf_tot_errors"] 
        self.amplitudes[suffix] = data

    def cdf_data(self, suffix):
        """
        Calculated Cumulative Density Function and associated uncertainties
        * cdf is done by adding bins
        * cdf errors are done by adding in quadrature errors of preceding bins
          so for pdf errors e^p_i the cdf error
              e^c_i = sum_{j=0}^{j=i}(e^p_j**2)**0.5
        """
        data = self.amplitudes[suffix]
        for key in ["all_upstream", "all_downstream"]:
            pdf_list = data[key]["corrected_pdf"]
            err_list = data[key]["pdf_tot_errors"]
            bins = range(len(pdf_list))

            cdf_list = [sum(pdf_list[:i+1]) for i in bins]
            cdf_list_errs = [0. for i in bins]
            for i in bins:
                err_sq = [err**2 for err in err_list[:i+1]]
                cdf_list_errs[i] = sum(err_sq)**0.5

            data[key]["corrected_cdf"] = cdf_list
            data[key]["cdf_tot_errors"] = cdf_list_errs


    def ratio_data(self, suffix):
        """
        Calculate ratio of PDF and CDF
        * ratio is just downstream/upstream (so > 1 means increased density)
        * errors are done by adding fractional errors in quadrature, so for a
          given bin with content c, error e, 
              e = c*((e_us/c_us)**2+(e_ds/c_ds)**2)**0.5
        """
        data = self.amplitudes[suffix]
        data["ratio"] = {}
        for key, err_key in ("corrected_pdf", "pdf_tot_errors"), ("corrected_cdf", "cdf_tot_errors"):
            pdf_list_tku = data["all_upstream"][key]
            pdf_list_tkd = data["all_downstream"][key]
            bins = range(len(pdf_list_tku))

            err_list_tku = [0. for i in bins]
            err_list_tkd = [0. for i in bins]
            ratio_pdf = [1. for i in bins]
            for i in bins:
                if pdf_list_tku[i] > 0.5:
                    err_list_tku[i] = data["all_upstream"][err_key][i]/pdf_list_tku[i]
                    ratio_pdf[i] = pdf_list_tkd[i]/pdf_list_tku[i]
                if pdf_list_tkd[i] > 0.5:
                    err_list_tkd[i] = data["all_downstream"][err_key][i]/pdf_list_tkd[i]

            data["ratio"][key] = ratio_pdf
            data["ratio"][err_key] = [(err_list_tku[i]**2+err_list_tkd[i]**2)*ratio_pdf[i] for i in bins]
        self.amplitudes[suffix] = data

    def get_asymm_error_graph(self, points, errors=None, norm=1., style=None, color=None, name="Graph"):
        graph = ROOT.TGraphAsymmErrors(len(points)-1)
        if errors == None:
            errors = [0. for i in points]
        for i, low_edge in enumerate(self.bin_edge_list[:-1]):
            high_edge = self.bin_edge_list[i+1]
            centre = (low_edge+high_edge)/2.
            graph.SetPoint(i, centre, points[i]/norm)
            graph.SetPointError(i, centre-low_edge, high_edge-centre, errors[i]/norm, errors[i]/norm)
        if style != None:
            graph.SetMarkerStyle(style)
        if color != None:
            graph.SetMarkerColor(color)
        graph.SetName(name)
        self.root_objects.append(graph)
        return graph

    def get_hist(self, graph_list, x_axis, y_axis, min_x = None, max_x = None, min_y = None, max_y = None):
        if min_x == None:
            min_x = min([graph.GetXaxis().GetXmin() for graph in graph_list])
        if max_x == None:
            max_x = max([graph.GetXaxis().GetXmax() for graph in graph_list])
        if min_y == None:
            min_y = min([graph.GetYaxis().GetXmin() for graph in graph_list])
        if max_y == None:
            max_y = max([graph.GetYaxis().GetXmax() for graph in graph_list])
        print "HIST", min_x, max_x, min_y, max_y
        hist = ROOT.TH2D(graph_list[0].GetName()+"_hist", ";"+x_axis+";"+y_axis,
                          1000, min_x, max_x, 1000, min_y, max_y)
        hist.SetStats(False)
        self.root_objects.append(hist)
        return hist

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
        text_box.AddText("MICE INTERNAL")
        text_box.Draw()
        self.root_objects.append(text_box)
        text_box = ROOT.TPaveText(0.6, 0.45, 0.89, 0.55, "NDC")
        text_box.SetFillColor(0)
        text_box.SetBorderSize(0)
        text_box.SetTextSize(0.03)
        text_box.SetTextAlign(12)
        text_box.AddText(self.config_anal['name'])
        text_box.Draw()
        self.root_objects.append(text_box)

    def pdf_plot(self, suffix):
        data = self.amplitudes[suffix]
        us_data = self.amplitudes[suffix]["all_upstream"]
        ds_data = self.amplitudes[suffix]["all_downstream"]
        name = "amplitude_pdf_"+suffix
        canvas = common.make_root_canvas(name)
        canvas.SetName(name)
        raw_upstream_graph = self.get_asymm_error_graph(us_data["pdf"],
                                                    style=20, color=ROOT.kBlue, name="Raw upstream")
        raw_downstream_graph = self.get_asymm_error_graph(ds_data["pdf"],
                                                    style=22, color=ROOT.kRed, name="Raw downstream")
        scraped_graph = self.get_asymm_error_graph(data["upstream_scraped"]["pdf"],
                                                    style=21, color=ROOT.kViolet+2, name="Raw scraped") # 25 for not raw
        
        if suffix == "reco":
            upstream_graph = self.get_asymm_error_graph(us_data["corrected_pdf"],
                                                        us_data["pdf_tot_errors"],
                                                        style=24, color=ROOT.kBlue, name = "Upstream")
            downstream_graph = self.get_asymm_error_graph(ds_data["corrected_pdf"],
                                                          ds_data["pdf_tot_errors"],
                                                        style=26, color=ROOT.kRed, name = "Downstream")
            graph_list = [upstream_graph, downstream_graph, scraped_graph, raw_upstream_graph, raw_downstream_graph]
        else:
            graph_list = [scraped_graph, raw_upstream_graph, raw_downstream_graph]
        hist = self.get_hist(graph_list, self.get_suffix_label(suffix)+" Amplitude [mm]", "Number")
        hist.Draw()
        for graph in graph_list:
            graph.Draw("SAMEP")
        self.text_box(graph_list)

        if self.config_anal["amplitude_chi2"]:
            upstream_chi2 = self.chi2_graph(suffix, "all_upstream")
            upstream_chi2.SetLineColor(ROOT.kBlue+2) 
            upstream_chi2.Draw("SAMEL")
            downstream_chi2 = self.chi2_graph(suffix, "all_downstream")
            downstream_chi2.SetLineColor(ROOT.kRed+2)
            downstream_chi2.Draw("SAMEL")

        canvas.Update()
        for a_format in ["eps", "pdf", "root", "png"]:
            canvas.Print(self.plot_dir+"amplitude_pdf_"+suffix+"."+a_format)
        return canvas

    def cdf_plot(self, suffix):
        us_data = self.amplitudes[suffix]["all_upstream"]
        ds_data = self.amplitudes[suffix]["all_downstream"]
        name = "amplitude_cdf_"+suffix
        canvas = common.make_root_canvas(name)
        canvas.SetName(name)

        upstream_graph = self.get_asymm_error_graph(us_data["corrected_cdf"],
                                                    us_data["cdf_tot_errors"],
                                                    style=24, color=ROOT.kBlue, name = "Upstream CDF")
        downstream_graph = self.get_asymm_error_graph(ds_data["corrected_cdf"],
                                                      ds_data["cdf_tot_errors"],
                                                    style=26, color=ROOT.kRed, name = "Downstream CDF")
        graph_list = [upstream_graph, downstream_graph]

        hist = self.get_hist(graph_list, self.get_suffix_label(suffix)+" Amplitude [mm]", "Cumulative Number")
        hist.Draw()
        for graph in graph_list:
            graph.Draw("SAMEP")
        self.text_box(graph_list)

        canvas.Update()
        for a_format in ["eps", "pdf", "root", "png"]:
            canvas.Print(self.plot_dir+"amplitude_cdf_"+suffix+"."+a_format)
        return canvas

    def cdf_ratio_plot(self, suffix):
        data = self.amplitudes[suffix]["ratio"]
        name = "cdf_ratio_"+suffix
        canvas = common.make_root_canvas(name)
        canvas.SetName(name)

        cdf_graph = self.get_asymm_error_graph(data["corrected_cdf"],
                                               data["cdf_tot_errors"],
                                               style=24, color=ROOT.kBlue, name = "CDF Ratio")
        graph_list = [cdf_graph]

        hist = self.get_hist(graph_list, self.get_suffix_label(suffix)+" Amplitude [mm]", "Cumulative Number Ratio")
        hist.Draw()
        for graph in graph_list:
            graph.Draw("SAMEP")
        self.text_box(graph_list)

        canvas.Update()
        for a_format in ["eps", "pdf", "root", "png"]:
            canvas.Print(self.plot_dir+canvas.GetName()+"."+a_format)
        return canvas

    def pdf_ratio_plot(self, suffix):
        data = self.amplitudes[suffix]["ratio"]
        name = "pdf_ratio_"+suffix
        canvas = common.make_root_canvas(name)
        canvas.SetName(name)

        pdf_graph = self.get_asymm_error_graph(data["corrected_pdf"],
                                               data["pdf_tot_errors"],
                                               style=24, color=ROOT.kBlue, name = "PDF Ratio")
        graph_list = [pdf_graph]

        hist = self.get_hist(graph_list, self.get_suffix_label(suffix)+" Amplitude [mm]", "Number Ratio")
        hist.Draw()
        for graph in graph_list:
            graph.Draw("SAMEP")
        self.text_box(graph_list)

        canvas.Update()
        for a_format in ["eps", "pdf", "root", "png"]:
            canvas.Print(self.plot_dir+canvas.GetName()+"."+a_format)
        return canvas


    def amplitude_scatter_plot(self, suffix):
        """
        Make a scatter plot of upstream versus downstream data
        """
        amp_dict_ds = self.amplitudes[suffix]["amplitude_dict_downstream"]
        amp_dict_us = self.amplitudes[suffix]["amplitude_dict_upstream"]
        amp_list_us = []
        amp_list_ds = []
        amp_list_delta = []

        suffix_label = self.get_suffix_label(suffix)
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
                                             amp_list_us, "US "+suffix_label+" Amplitude [mm]",
                                             amp_list_delta, suffix_label+" US Amplitude - DS Amplitude [mm]", xmin=0., xmax=100., ymin=-50., ymax=50.)
        hist.SetTitle(self.config_anal['name'])
        hist.Draw()
        graph.Draw("P")
        canvas.Update()
        for format in ["eps", "png", "root"]:
            canvas.Print(self.plot_dir+"amplitude_delta_"+suffix+"_scatter."+format)

        canvas = common.make_root_canvas("delta_amplitude_hist")
        canvas.Draw()
        canvas.SetFrameFillColor(scripts.utilities.get_frame_fill())
        hist = common.make_root_histogram("delta amplitude hist",
                                          amp_list_us, suffix_label+" US Amplitude [mm]", 100,
                                          amp_list_delta, suffix_label+" US Amplitude - DS Amplitude [mm]", 100,
                                          xmin=0., xmax=100., ymin=-50., ymax=50.)
        hist.SetTitle(self.config_anal['name'])
        hist.Draw("COLZ")
        canvas.Update()
        for format in ["eps", "png", "root"]:
            canvas.Print(self.plot_dir+"amplitude_delta_"+suffix+"_hist."+format)

    def get_suffix_label(self, suffix):
        return {"all_mc":"MC truth (all)",
                "reco_mc":"MC truth (recon)",
                "reco":"Reconstructed"}[suffix]

    def append_data(self):
        """
        Add data to the amplitude calculation (done at death time)
        
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
            station_us = self.config.mc_plots["mc_stations"]["tku"][0]
            station_ds = self.config.mc_plots["mc_stations"]["tkd"][0]

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
        fout = open(self.plot_dir+"/amplitude.json", "w")
        for suffix in self.amplitudes:
            if type(self.amplitudes[suffix]) != type({}):
                continue
            print self.amplitudes[suffix].keys()
            try:
                del self.amplitudes[suffix]["amplitude_dict_upstream"]
                del self.amplitudes[suffix]["amplitude_dict_downstream"]
            except KeyError:
                pass
        out_str = json.dumps(self.amplitudes, indent=2)
        fout.write(out_str)
        fout.close()

    def empty_data(self):
        n_bins = 21
        diagonal = [[0. for i in range(n_bins)] for j in range(n_bins)]
        for i in range(n_bins):
            diagonal[i][i] = 1.
        self.amplitudes = {
          "momentum_uncertainty":{
              "all_upstream":{
                  "migration_matrix":diagonal
              },
              "all_downstream":{
                  "migration_matrix":diagonal
              },
          }, 
          "inefficiency":{
              "all_upstream":{
                  "pdf_ratio":[1. for i in range(n_bins)]
              },
              "all_downstream":{
                  "pdf_ratio":[1. for i in range(n_bins)]
              },
          },
          "crossing_probability":{
              "all_upstream":{
                  "migration_matrix":diagonal
              },
              "all_downstream":{
                  "migration_matrix":diagonal
              },
          }, 
          "systematics":[],
          "source":"",
          "scale":1.,
        }

    def load_one_error(self, file_name, scale):
        fin = open(file_name)
        amp_str = fin.read()
        amplitudes = json.loads(amp_str)
        amplitudes["source"] = file_name
        amplitudes["scale"] = scale
        # we only want to keep the errors; not the actual pdfs
        del amplitudes["reco"]
        del amplitudes["all_mc"]
        del amplitudes["reco_mc"]
        try:
            del amplitudes["field_uncertainty"]
        except KeyError:
            pass
        return amplitudes

    def load_errors(self):
        if self.config_anal["amplitude_source"] != None and self.config_anal["amplitude_mc"]:
            raise RuntimeError("Conflicted - do I get errors from mc or from amplitude source?")
        if self.config_anal["amplitude_source"] == None:
            self.empty_data()
            return
        # base correction factors
        self.amplitudes = self.load_one_error(self.config_anal["amplitude_source"], 1.)
        # systematic uncertainties
        self.amplitudes["systematics"] = []
        err_src = self.config_anal["amplitude_systematic_reference"]
        error = self.load_one_error(err_src, None)
        self.amplitudes["systematics"].append(error)
        for err_src, scale in self.config_anal["amplitude_systematic_sources"].iteritems():
            error = self.load_one_error(err_src, scale)
            self.amplitudes["systematics"].append(error)
            
    root_objects = []
