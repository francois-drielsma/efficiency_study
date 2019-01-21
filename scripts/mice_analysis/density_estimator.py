import sys
import datetime
import bisect
import operator
import json
import numpy
import tempfile
import glob
import shutil
import os
import copy


import ROOT
import xboa.common as common
from xboa.hit import Hit
from xboa.bunch import Bunch

import utilities.chi2_distribution
from utilities.binomial_confidence_interval import BinomialConfidenceInterval

from mice_analysis.amplitude.amplitude_data_binned import AmplitudeDataBinned
from mice_analysis.amplitude.amplitude_data_binless import AmplitudeDataBinless
from mice_analysis.analysis_base import AnalysisBase
from mice_analysis.amplitude.plot_amplitude_data import PlotAmplitudeData
from mice_analysis.amplitude.plot_amplitude_analysis import PlotAmplitudeAnalysis


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
        self.clear_amplitude_data()
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
        self.plotter = PlotAmplitudeAnalysis(self)
        self.calculate_corrections = self.config_anal["amplitude_corrections"] == None

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
        self.plotter.plot_dir = self.plot_dir
        self.all_mc_data_us.clear()
        self.reco_mc_data_us.clear()
        self.reco_data_us.clear()
        self.all_mc_data_ds.clear()
        self.reco_mc_data_ds.clear()
        self.reco_data_ds.clear()
        self.clear_amplitude_data()
        self.load_errors()
        try:
            os.mkdir(self.plot_dir+"/phase_space")
        except OSError:
            pass
        try:
            os.mkdir(self.plot_dir+"/weighting")
        except OSError:
            pass
        self.append_data()

    def process(self):
        self.append_data()

    def death(self):
        self.delta_amplitude_calc("reco")
        if self.config_anal["amplitude_mc"]:
            self.delta_amplitude_calc("all_mc")
            self.delta_amplitude_calc("reco_mc")
        if self.calculate_corrections:
            self.corrections_calc("all_upstream")
            self.corrections_calc("all_downstream")
            for target in "all_upstream", "all_downstream":
                self.plotter.systematics_plot(target)
        if self.config_anal["amplitude_mc"]:
            for suffix in ["reco_mc", "all_mc"]:
                self.corrections_and_uncertainties(suffix)
                self.cdf_data(suffix)
                self.ratio_data(suffix)
                self.plotter.pdf_plot(suffix)
                self.plotter.cdf_plot(suffix)
                self.plotter.pdf_ratio_plot(suffix)
                self.plotter.cdf_ratio_plot(suffix)
                self.plotter.amplitude_scatter_plot(suffix)
        if self.calculate_corrections:
            self.plotter.amplitude_residuals_plot("upstream")
            self.plotter.amplitude_residuals_plot("downstream")
        for suffix in ["reco"]:
            self.corrections_and_uncertainties(suffix)
            self.cdf_data(suffix)
            self.ratio_data(suffix)
            self.plotter.pdf_plot(suffix)
            self.plotter.cdf_plot(suffix)
            self.plotter.pdf_ratio_plot(suffix)
            self.plotter.cdf_ratio_plot(suffix)
            self.plotter.amplitude_scatter_plot(suffix)
        self.print_data()
        self.move_data()

    def move_data(self):
        target_dir = self.plot_dir+"/data"
        os.mkdir(target_dir)
        file_list = glob.glob(self.a_dir+"/*")
        for a_file in file_list:
            file_name = os.path.split(a_file)[1]
            shutil.copy(a_file, target_dir+"/"+file_name)

    def delta_amplitude_calc(self, suffix):
        """
        Calculate amplitudes for upstream and downstream sample
        * suffix: "reco", "all_mc" or "reco_mc" for different sampling criteria,
          as defined during the append_data step. reco is supposed to be
          recontructed data; all_mc is all monte carlo data; reco_mc is all
          monte carlo data that is also in the reco sample
        Calculates amplitudes using the AmplitudeData algorithms; fills
        dicts; then fills pdfs and migration matrices
        """
        data_0, data_1 = {
            "all_mc":(self.all_mc_data_us, self.all_mc_data_ds),
            "reco":(self.reco_data_us, self.reco_data_ds),
            "reco_mc":(self.reco_mc_data_us, self.reco_mc_data_ds),
        }[suffix]
        print "Amplitude", suffix
        print "  starting delta_amplitude_calc    ", datetime.datetime.now()
        data = self.amplitudes[suffix]
        print "  Upstream amplitude calculation...  ", datetime.datetime.now()
        fractional_amplitude_dict_0 = data_0.fractional_amplitude()
        print "  n_events", len(fractional_amplitude_dict_0)
        print "  Downstream amplitude calculation...", datetime.datetime.now()
        fractional_amplitude_dict_1 = data_1.fractional_amplitude()
        print "  n_events", len(fractional_amplitude_dict_1)
        print "  ... done                           ", datetime.datetime.now()

        data["amplitude_dict_upstream"] = fractional_amplitude_dict_0
        data["amplitude_dict_downstream"] = fractional_amplitude_dict_1
        self.set_pdfs(data, "all_upstream", fractional_amplitude_dict_0)
        self.set_pdfs(data, "all_downstream", fractional_amplitude_dict_1)
        for name, amp_data in [("all_upstream", data_0),
                           ("all_downstream", data_1)]:
            data[name]["emittance"] = amp_data.get_emittance()
            data[name]["weight"] = amp_data.get_n_events()
            data[name]["covariance_matrix"] = amp_data.cov.tolist()
            data[name]["bin_edge_list"] = self.bin_edge_list
        print "  upstream  ", data["all_upstream"]["pdf"], "sum", sum(data["all_upstream"]["pdf"])
        print "  downstream", data["all_downstream"]["pdf"], "sum", sum(data["all_downstream"]["pdf"])

        data["migration_matrix"] = self.migration_matrix(fractional_amplitude_dict_0, fractional_amplitude_dict_1, False)
        self.set_delta_pdfs(data, fractional_amplitude_dict_0, fractional_amplitude_dict_1)
        self.amplitudes[suffix] = data

        us_plotter = PlotAmplitudeData(data_0, self.plot_dir, "us")
        us_plotter.plot()

        ds_plotter = PlotAmplitudeData(data_1, self.plot_dir, "ds")
        ds_plotter.plot()

        print "  done amplitude calc                ", datetime.datetime.now()

    def get_bin_centre_list(self):
        """
        Return the list of bin centres, calculated using self.bin_edge_list
        """
        bin_width = self.config.amplitude_bin_width
        bin_centre_list = [bin_edge+bin_width/2. for bin_edge in self.bin_edge_list[:-1]]
        return bin_centre_list

    def set_delta_pdfs(self, data, a_dict_us, a_dict_ds):
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
        self.set_pdfs(data, "upstream_scraped", scraped_dict)
        self.set_pdfs(data, "upstream_not_scraped", not_scraped_dict)

    def set_pdfs(self, data, us_ds, fractional_amplitude_dict):
        """
        Bin the data in the amplitude dict and return as a list
        """
        bin_centre_list = self.get_bin_centre_list()
        if us_ds not in data:
            data[us_ds] = {}
        data[us_ds]["bin_centre_list"] = bin_centre_list
        data[us_ds]["pdf"] = [0.]*(len(self.bin_edge_list))
        amplitude_index = 0
        for amplitude in fractional_amplitude_dict.values():
            amp_0_bin = bisect.bisect_left(self.bin_edge_list, amplitude)-1
            data[us_ds]["pdf"][amp_0_bin] += 1.

    def migration_matrix_alt(self, fractional_amplitude_dict_0, fractional_amplitude_dict_1):
        """
        Generate a migration matrix to migrate from distribution 0 to distribution 1
        * fractional_amplitude_dict_0: dictionary of <event_id>:amplitude
        * fractional_amplitude_dict_1: dictionary of <event_id>:amplitude
        <event_id> can be any hashable id; it has to be the same for dict_0 and dict_1

        Builds the number of events N_0 and N_1 in each distribution; set matrix 
        so that N_1 = M N_0; overflow bin is allowed
        """
        n_bins = len(self.bin_edge_list)
        # bin the data
        n_0 = [0. for i in range(n_bins)]
        n_1 = [0. for i in range(n_bins)]
        for amp in fractional_amplitude_dict_0.values():
            amp_bin = bisect.bisect_left(self.bin_edge_list, amp)-1
            n_0[amp_bin] += 1.
        for amp in fractional_amplitude_dict_1.values():
            amp_bin = bisect.bisect_left(self.bin_edge_list, amp)-1
            n_1[amp_bin] += 1.

        migration_matrix = [[0. for i in range(n_bins)]  for j in range(n_bins)]
        for i in range(n_bins):
            for j in range(n_bins):
                if n_0[j] < 0.5:
                    if i == j:
                        migration_matrix[i][j] = 1.
                    else:
                        migration_matrix[i][j] = 0.
                else:
                    migration_matrix[i][j] = n_1[i]/n_0[j]

        return migration_matrix


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
            migration_matrix = self.row_normalise_matrix(migration_matrix, 1.)
        return migration_matrix

    def row_normalise_matrix(self, matrix, default):
        """
        Return a copy of matrix that is normalised to the number in row so that
        sum_i(m_ij) = 1.
        - matrix: matrix to be normalised
        - default: if row_sum is 0., the diagonal terms go to "default"
        """
        migration_matrix = copy.deepcopy(matrix)
        for i, row in enumerate(migration_matrix):
            row_sum = float(sum(row))
            for j, element in enumerate(row):
                if abs(row_sum) < 1e-3:
                    if i == j:
                        migration_matrix[i][j] = default # identity
                    else:
                        migration_matrix[i][j] = 0.
                else:
                    migration_matrix[i][j] /= row_sum
        return migration_matrix


    def corrections_calc(self, target):
        """
        Calculate the pdf corrections (migration matrix and efficiency matrix)
        
        Uses the MC to generate corrections
        """
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

        cutoff = self.config_anal["amplitude_inefficiency_cutoff"]
        cutoff_index = bisect.bisect_left(self.bin_edge_list, cutoff)
        all_mc_sum = sum(all_mc_pdf[cutoff_index:])
        reco_mc_sum = sum(reco_mc_pdf[cutoff_index:])
        inefficiency_averaged = copy.deepcopy(inefficiency)
        for i in range(cutoff_index, len(inefficiency_averaged)):
            inefficiency_averaged[i] = all_mc_sum/reco_mc_sum

        ineff_migration = self.migration_matrix_alt(fractional_reco_mc, fractional_reco)
        self.amplitudes["inefficiency"][target] = {
            "pdf_ratio":inefficiency,
            "pdf_ratio_averaged":inefficiency_averaged,
            "migration_matrix":ineff_migration
        }
        self.amplitudes["crossing_probability"][target] = {}
        self.amplitudes["crossing_probability"][target]["migration_matrix"] = \
                      self.migration_matrix(fractional_reco, fractional_reco_mc, True)

    def matrix_str(self, matrix, rounding=2, width=10):
        """Represent a matrix as a string (space separated)"""
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
        suffix_label = self.plotter.get_suffix_label(suffix)
        self.plotter.matrix_plot(migration_matrix, "migration matrix",
                         "DS "+suffix_label+" Amplitude [mm]",
                         "US "+suffix_label+" Amplitude [mm]")
        downstream_errors = [sum(row)**0.5 for row in error_matrix]
        error_matrix_sqrt = [[element**0.5 for element in row] for row in error_matrix]
        self.plotter.matrix_plot(error_matrix_sqrt, "stats error",
                         "DS "+suffix_label+" Amplitude [mm]",
                         "US "+suffix_label+" Amplitude [mm]")
        print "upstream pdf\n", self.matrix_str(row_sum)
        print "migration matrix\n", self.matrix_str(migration_matrix)
        print "error matrix\n", self.matrix_str(error_matrix)
        print "downstream errors\n", self.matrix_str(downstream_errors)
        return downstream_errors

    def do_migrations(self, suffix, us_ds, migrations):
        """
        Calculate the pdf that results from applying the specified migrations.
        - suffix: Controls the source data and how the migration is applied.
        If suffix is "reco", this will apply the correction matrix 
        (for detector resolution) and then the efficiency correction; if suffix 
        is "reco_mc" this will apply only the efficiency correction; else no 
        correction is applied at all and the pdf is returned with no change
        - us_ds: control whether to apply the data on upstream or downstream
        - migrations: source for the migrations
        """
        pdf = numpy.array(self.amplitudes[suffix][us_ds]["pdf"])
        if suffix == "reco":
            reco_mc_matrix = migrations["crossing_probability"][us_ds]["migration_matrix"]
            reco_mc_matrix = numpy.transpose(numpy.array(reco_mc_matrix))
            pdf = numpy.dot(reco_mc_matrix, pdf)
        if suffix == "reco" or suffix == "reco_mc":
            pdf = pdf*numpy.array(migrations["inefficiency"][us_ds]["pdf_ratio"])
        return pdf.tolist()

    def calculate_detector_systematics(self, suffix, us_ds):
        """
        Calculate the systematic errors in the reconstruction of tracks
        """
        data = self.amplitudes[suffix]
        bins = range(len(data[us_ds]["pdf"]))
        sys_error_list = [0 for i in bins]
        if data["detector_reference"] == None:
            return [0. for i in bins]

        print "\nReconstruction systematic errors", us_ds
        migrations = data["detector_reference"]
        ref_pdf_list = self.do_migrations(suffix, us_ds, migrations)

        print 'ref pdf: ', [format(p, '6.4g') for p in ref_pdf_list]

        systematics_list = data[us_ds]["detector_systematics"]
        for i, migrations in enumerate(systematics_list):
            scale = migrations["scale"]
            source = migrations["source"]
            sys_pdf_list = self.do_migrations(suffix, us_ds, migrations)
            print "sys pdf: ", [format(p, '6.4g') for p in sys_pdf_list]
            print "  ", i, str(scale).ljust(6), "  err:",
            for j in bins:
                err = (sys_pdf_list[j] - ref_pdf_list[j])*scale
                print format(err, '6.2e'),
                sys_error_list[j] = (sys_error_list[j]**2+err**2)**0.5
            print
        print "\nsum: ", [format(err, '6.2e') for err in sys_error_list]
        return sys_error_list

    def calculate_performance_systematics(self, suffix, us_ds):
        """
        Calculate the systematic errors in the channel performance. We study the
        migration from tku bin i to tkd bin j. We characterise the uncertainty 
        in the downstream distributions from those uncertainties.
        
        The uncertainty in the migration from ith bin to jth bin m_ij is 
        calculated by taking
            m(err_k)_ij = m(err)_ij - m(ref)_ij
            m(err)_ij = (sum_k m**2(err_k)_ij)**0.5
            n(err)_i = sum_j m(err)_ij n(tku)_j
        Stored migration matrix m_ij is number in (TKU bin i) AND (TKD bin j).
        """
        data = self.amplitudes[suffix]
        tku_ref = data["all_upstream"]["pdf"]
        bins = range(len(tku_ref))
        if data["performance_reference"] == None:
            return [0. for i in bins]
        print "\nPerformance systematic errors", suffix
        migrations = data["performance_reference"]
        ref_matrix = migrations[suffix]["migration_matrix"]
        ref_matrix = self.row_normalise_matrix(ref_matrix, 0.)
        err_matrix = [[0. for i in bins] for j in bins]

        systematics_list = data[us_ds]["performance_systematics"]
        for i, systematic in enumerate(systematics_list):
            sys_matrix = systematic[suffix]["migration_matrix"]
            sys_matrix = self.row_normalise_matrix(sys_matrix, 0.)
            scale = systematic["scale"]
            source = systematic["source"]
            print "  ", i, scale, source, "  error matrix:"
            for i in bins:
                print "   ",
                for j in bins:
                    err = sys_matrix[i][j] - ref_matrix[i][j]
                    print format(err, '10.4g'),
                    err_matrix[i][j] += (err*scale)**2
                print
            print
        err_matrix = [[cell**0.5 for cell in row] for row in err_matrix]
        ref = data["performance_reference"][suffix]
        print 'Finally TKU pdf:\n   ',
        for i in bins:
              pdf = ref["all_upstream"]["pdf"]
              print format(pdf[i], '6.4g'),
        print '\nFinally TKD pdf:\n   ',
        for i in bins:
              pdf = ref["all_downstream"]["pdf"]
              print format(pdf[i], '6.4g'),
        tku_list = ref["all_upstream"]["pdf"]
        sys_error_list = [0. for i in bins]
        for i in bins:
            for j in bins:
                sys_error_list[j] += err_matrix[i][j] * tku_list[i]
        print 'Finally err pdf:\n   ',
        for i in bins:
              pdf = sys_error_list
              print format(pdf[i], '6.4g'),
        print "\nFinally error matrix\n   ",
        print self.matrix_str(err_matrix)
        print "\nFinally reference migration matrix\n   ",
        print self.matrix_str(ref_matrix)

        sys_error_list = [0. for i in bins]
        for i in bins:
            for j in bins:
                sys_error_list[j] += err_matrix[i][j] * tku_ref[i]
        return sys_error_list

    def corrections_and_uncertainties(self, suffix):
        """
        Calculate corrected pdfs and uncertainties
        * corrected pdf is given by bin_i = Sum_j(Eff_i*M_ij*bin_j) where M_ij
            is the migration matrix and Eff_i is the efficiency
        * statistical errors are given by sum in quadrature of binomial errors
            on each bin migration
        * total errors are given by sum in quadrature of statistical errors and
          systematic errors
        """
        data = self.amplitudes[suffix]
        raw_pdf_list_tku = data["all_upstream"]["pdf"]
        raw_pdf_list_tkd = data["all_downstream"]["pdf"]
        migration_matrix = data["migration_matrix"]
        data["all_upstream"]["pdf_stats_errors"] = [0. for bin in raw_pdf_list_tku]
        data["all_downstream"]["pdf_stats_errors"] = \
                                     self.stats_errors(migration_matrix, suffix)
        for us_ds in ["all_upstream", "all_downstream"]:
            # apply basic migrations
            print "Doing error correction for", suffix, us_ds
            migrations = self.amplitudes
            pdf_list = self.do_migrations(suffix, us_ds, migrations)
            data[us_ds]["corrected_pdf"] = pdf_list
            print "Finding systematic errors for", suffix, us_ds
            reco_sys_list = self.calculate_detector_systematics(suffix, us_ds)
            perf_sys_list = self.calculate_performance_systematics(suffix, us_ds)
            sys_error_list = [(reco_sys_list[i]**2+perf_sys_list[i]**2)**0.5 \
                                                  for i in range(len(pdf_list))]
            data[us_ds]["pdf_sys_errors"] = sys_error_list
            print "    sys errors:   ", data[us_ds]["pdf_sys_errors"] 
            print "    stats errors: ", data[us_ds]["pdf_stats_errors"]
            print "    pdf:          ", pdf_list
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
            stats_err_list = data[key]["pdf_stats_errors"]
            sys_err_list = data[key]["pdf_sys_errors"]
            bins = range(len(pdf_list))

            cdf_list = [sum(pdf_list[:i+1]) for i in bins]
            cdf_list_stats_errs = [0. for i in bins]
            cdf_list_sys_errs = [0. for i in bins]
            for i in bins:
                stats_err_sq = [err**2 for err in stats_err_list[:i+1]]
                sys_err_sq = [err**2 for err in sys_err_list[:i+1]]
                cdf_list_stats_errs[i] = sum(stats_err_sq)**0.5
                cdf_list_sys_errs[i] = sum(sys_err_sq)**0.5

            data[key]["corrected_cdf"] = cdf_list
            data[key]["cdf_stats_errors"] = cdf_list_stats_errs
            data[key]["cdf_sys_errors"] = cdf_list_sys_errs

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
        for key, stats_err_key, sys_err_key in [("corrected_pdf", "pdf_stats_errors", "pdf_sys_errors",),
                                                ("corrected_cdf", "cdf_stats_errors", "cdf_sys_errors",)]:
            pdf_list_tku = data["all_upstream"][key]
            pdf_list_tkd = data["all_downstream"][key]
            bins = range(len(pdf_list_tku))

            ratio_pdf = [1. for i in bins]
            for i in bins:
                if pdf_list_tku[i] > 0.5:
                    ratio_pdf[i] = pdf_list_tkd[i]/pdf_list_tku[i]
            data["ratio"][key] = ratio_pdf

            for err_key in stats_err_key, sys_err_key:
                err_list_tku = [0. for i in bins]
                err_list_tkd = [0. for i in bins]
                data["ratio"][err_key] = [0. for i in bins]
                for i in bins:
                    if pdf_list_tku[i] > 0.5:
                        err_list_tku[i] = data["all_upstream"][err_key][i]/pdf_list_tku[i]
                    if pdf_list_tkd[i] > 0.5:
                        err_list_tkd[i] = data["all_downstream"][err_key][i]/pdf_list_tkd[i]
                    data["ratio"][err_key][i] = (err_list_tku[i]**2+err_list_tkd[i]**2)**0.5*ratio_pdf[i]
                print "amplitude_analysis.ratio_data", suffix, key, err_key+":", data["ratio"][err_key]

        self.amplitudes[suffix] = data

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
        print "Loaded upstream:"
        print "    reco:   ", len(hits_reco_us)
        print "    all mc: ", len(hits_all_mc_us)
        print "    reco mc:", len(hits_reco_mc_us)
        print "Loaded downstream:"
        print "    reco:   ", len(hits_reco_ds)
        print "    all mc: ", len(hits_all_mc_ds)
        print "    reco mc:", len(hits_reco_mc_ds)
        return

    def print_data(self):
        """
        Write the amplitude distributions to disk.
        """
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

    def clear_amplitude_data(self):
        n_bins = 21
        diagonal = [[0. for i in range(n_bins)] for j in range(n_bins)]
        for i in range(n_bins):
            diagonal[i][i] = 1.
        self.amplitudes = {
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
          "source":"",
          "scale":1.,
        }
        data = {
            "performance_reference":None,
            "detector_reference":None,
            "all_upstream":{
                "detector_systematics":[],
                "performance_systematics":[],
                "detector_systematics_output":{},
                "performance_systematics_output":{},
            },
            "all_downstream":{
                "detector_systematics":[],
                "performance_systematics":[],
                "detector_systematics_output":{},
                "performance_systematics_output":{},
            },
            "pdf":[],
            "cdf":[],
            
        }
        for suffix in "reco", "all_mc", "reco_mc":
            self.amplitudes[suffix] = copy.deepcopy(data)


    def load_corrections(self, file_name):
        """
        Load the amplitude corrections to be applied during this amplitude
        analysis. Loads the inefficiency and crossing_probability
        """
        fin = open(file_name)
        amp_str = fin.read()
        src_amplitudes = json.loads(amp_str)
        src_amplitudes["source"] = file_name
        for key in "inefficiency", "crossing_probability":
            self.amplitudes[key] = src_amplitudes[key]
        

    def load_one_error(self, file_name, scale):
        """
        Load the amplitude analysis output for a given uncertainty source
        """
        fin = open(file_name)
        amp_str = fin.read()
        amplitudes = json.loads(amp_str)
        amplitudes["source"] = file_name
        amplitudes["scale"] = scale
        # for performance error we want the pdfs
        # for recon error we want the corrections
        # for storage concerns we don't store the amplitude dictionaries
        for tgt_1 in "reco", "all_mc", "reco_mc":
            for tgt_2 in "amplitude_dict_upstream", "amplitude_dict_downstream":
                try:
                    del amplitudes[tgt_1][tgt_2]
                except KeyError:
                    pass
                    #print "Could not find", tgt_1, tgt_2, "while loading errors"
        for key in "field_uncertainty", "momentum_uncertainty":
            # legacy systematic corrections
            try:
                del amplitudes[key]
            except KeyError:
                pass
        return amplitudes

    def load_errors(self):
        """
        Two "classes" of systematic errors;
        * systematic errors on the reconstruction are contained in the
          correction factors. For these we store the correction factors and 
          compare to the reference correction factors
        * systematic errors on the performance are contained in the actual
          amplitude pdfs. For these we store the bin-by-bin fractional
          difference between the amplitude pdf and reference.
        """
        if self.calculate_corrections:
            self.clear_amplitude_data()
            return

        # set base correction factors
        self.load_corrections(self.config_anal["amplitude_corrections"])

        systematics = self.config_anal["amplitude_systematics"]
        for suffix in systematics:
            print "Loading", suffix
            if suffix not in self.amplitudes:
                self.amplitudes[suffix] = {}
            for ref_key in ["detector_reference", "performance_reference"]:
                ref_src = systematics[suffix][ref_key]
                if ref_src == None:
                    self.amplitudes[suffix][ref_key] = None
                else:
                    self.amplitudes[suffix][ref_key] = \
                                              self.load_one_error(ref_src, None)
                print "  Loaded reference", suffix, ref_key, ref_src, \
                                          type(self.amplitudes[suffix][ref_key])
            for us_ds in ["all_upstream", "all_downstream"]:
                if us_ds not in self.amplitudes[suffix]:
                    self.amplitudes[suffix][us_ds] = {}
                for key in ["detector_systematics", "performance_systematics"]:
                    err_src_dict = systematics[suffix][us_ds][key]
                    self.amplitudes[suffix][us_ds][key] = [
                        self.load_one_error(err_src, scale) \
                              for err_src, scale in err_src_dict.iteritems()
                    ]
                    print "  Loaded", len(self.amplitudes[suffix][us_ds][key]), us_ds, key, "systematics"

