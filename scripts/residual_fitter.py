import sys
import json
import copy
import numpy

import ROOT
import xboa.common as common
from xboa.hit import Hit
from xboa.bunch import Bunch
import Configuration
import maus_cpp.field
import maus_cpp.globals

from scripts.extrapolate_through_detectors import ExtrapolateTrackPoints
import scripts.field_murgler

class ResidualFitter(object):
    """
    Fit a transfer matrix to a transfer matrix
    """
    def __init__(self, config, config_anal, data_loader):
        """
        Initialises the fitter
        """
        self.config = config
        self.config_anal = config_anal
        self.data_loader = data_loader
        self.event_list = []
        self.extrapolate_list = []
        # best guess i.e. output from latest fit iteration
        self.best_guess = {}
        # maximum number of iterations for each step of the fitter
        self.max_iterations = self.config.residual_fitter["max_iterations"]
        # fitter will attempt to minimise sum chi2 within +- 0.1
        self.resolution = self.config.residual_fitter["resolution"]
        # gets filled with the calculated residuals on each iteration of the
        # fitter
        self.residuals = {}
        # iteration counter
        self._iteration = 0
        self.track_source = ""
        self.residuals_source = ""
        self.direction = "both"
        self.residuals_variables = []
        self.optimisation_axes = {}
        self.optimisation_currents = {}

        self.murgler = scripts.field_murgler.Murgler(False)
        self.fout = None

    def fit(self):
        """
        Attempt to fit the transfer matrix, given statistical errors on each
        term in the transfer matrix.

        The fit proceeds in a number of steps
        1. Attempt to fit scale factor; looks at the upper quadrant of the
           focussing part of the matrix i.e. terms (0,1; 0,2; 1,2, 1,2,) in the
           transfer matrix.
        2. Attempt to fit x and x'; looks at the horizontal kick terms i.e.
           terms (0,0; 1,0) in the transfer matrix.
        3. Attempt to fit y and y'; looks at the vertical kick terms i.e.
           terms (2,0; 3,0) in the transfer matrix.
        """
        self.pick_events(self.data_loader)
        self.best_guess = {}
        self._iteration = 0
        ResidualFitter.fitter = self
        print "Running minuit to find scale_factor"
        print "Test run:"
        self.run_once(True, True)
        for fit in self.config.residual_fitter["fits"]:
            self.optimisation_axes = {fit["name"]:{"z_pos":fit["z_pos"]}}
            self.track_source = fit["track_source"]
            self.direction = fit["direction"]
            self.residuals_source = fit["detector"]
            self.residuals_variables = fit["variables"]
            self.optimisation_currents = {}
                      #"SSD_M2":[{"name":"MatchCoil2_1", "j":None}],
                      #"SSD_ECE":[{"name":"EndCoil1_1", "j":None},
                      #           {"name":"CenterCoil_1", "j":None},
                      #           {"name":"EndCoil2_1", "j":None}]
            print "Run optimisation"
            self._do_one_fit()
            print "Finished optimisation - rounding up"
            self.run_once(True, True)

    def _do_one_fit(self):
        """
        Run a single iteration of the fitting code
        - fit_vars: list of vars from self.seed to fit for
        - score_function_terms: terms in the transfer matrix that will be used
          for calculating sum(chi2).
        """
        ierr = ROOT.Long()
        self.minuit = ROOT.TMinuit(len(self.optimisation_axes)*4 + len(self.optimisation_currents))
        # tell minuit the parameters that we vary
        parameter_index = 0
        for key in self.optimisation_axes:
            x_min, x_ref, x_max, x_err = -20., 0., 20., 1. # mm
            y_min, y_ref, y_max, y_err = -20., 0., 20., 1. # mm
            xp_min, xp_ref, xp_max, xp_err = -0.02, 0., 0.02, 0.001 # rad
            yp_min, yp_ref, yp_max, yp_err = -0.02, 0., 0.02, 0.001 # rad
            self.minuit.mnparm(parameter_index, key+" x", x_ref, x_err, x_min, x_max, ierr)
            parameter_index +=1
            self.minuit.mnparm(parameter_index, key+" y", y_ref, y_err, y_min, y_max, ierr)
            parameter_index +=1
            self.minuit.mnparm(parameter_index, key+" xp", xp_ref, xp_err, xp_min, xp_max, ierr)
            #self.minuit.FixParameter(parameter_index)
            parameter_index +=1
            self.minuit.mnparm(parameter_index, key+" yp", yp_ref, yp_err, yp_min, yp_max, ierr)
            #self.minuit.FixParameter(parameter_index)
            parameter_index +=1
        for key in self.optimisation_currents:
            j_min, j_ref, j_max, j_err = 0.95, 1., 1.05, 0.01 # fractional
            self.minuit.mnparm(parameter_index, key+" j", j_ref, j_err, j_min, j_max, ierr)
            parameter_index +=1
        self.minuit.SetFCN(self.score_function)
        #self.minuit.SetPrintLevel(2) # -1 to silence
        args = numpy.array([self.max_iterations, self.resolution],
                            dtype=numpy.float64)
        # run the minimiser
        self.do_tracking()
        self.minuit.mnexcm("SIMPLEX", args, 2, ierr)
        # update the best guess from minuit and errors and print them        
        print "Final fits and estimated errors"
        self.update_from_minuit(True)

    def pick_events(self, data_loader):
        print "Picking events..."
        index = 0
        n_selected = 0
        target_count = self.config.residual_fitter["n_events"]
        p_bins =  copy.deepcopy(self.config_anal["p_bins"])
        self.event_list = [[] for bin in p_bins]
        extrapolate = ExtrapolateTrackPoints(self.config, self.config_anal, self.data_loader, False)
        # while we still have events to consider and not all bins are full       
        while index < len(data_loader.events) and \
              n_selected < len(self.event_list)*target_count:
            event = data_loader.events[index]
            index += 1
            if ExtrapolateTrackPoints.is_cut(event, self.config):
                continue
            for bin_index, bounds in enumerate(p_bins):
                if len(self.event_list[bin_index]) == target_count:
                    continue
                tku_p = event["tku"]["p"]
                # check if momentum is in the bin
                if tku_p < bounds[0] or tku_p > bounds[1]:
                    continue
                # check for extrapolation cuts (i.e. PID)
                event_copy = copy.deepcopy(event) # clean copy with no extrapolation
                extrapolate.extrapolate_event(event)
                extrapolate.extrapolation_cuts(event)
                if ExtrapolateTrackPoints.is_cut(event, self.config):
                    continue
                # keep the clean copy; we don't want the extrapolation
                self.event_list[bin_index].append(event_copy)
                n_selected += 1
                if n_selected % 1000 == 0:
                    print n_selected, "of", len(self.event_list)*target_count
                # if the bin is full, remove it from consideration
        print "Picked events after", index, "of", len(data_loader.events), "events" 
        for i in range(len(self.event_list)):
            print "   ", p_bins[i], len(self.event_list[i])

    def update_from_minuit(self, print_parameters):
        """
        update the best guess and estimated error based on minuit parameters
        """
        current, error = ROOT.Double(), ROOT.Double()
        parameter_index = 0
        self.best_guess = {"axes":{}, "currents":{}}
        for key in self.optimisation_axes:
            self.best_guess["axes"][key] = {}
            if print_parameters:
                print key,
            for var in ["x", "y", "xp", "yp"]:
                self.minuit.GetParameter(parameter_index, current, error)
                self.best_guess["axes"][key][var] = {"value":float(current), "error":float(error)}
                parameter_index += 1
                if print_parameters:
                    print var, current, "+-", error,
        for key in self.optimisation_currents:
            self.best_guess["currents"][key] = {}
            if print_parameters:
                print key,
            for var in ["j"]:
                self.minuit.GetParameter(parameter_index, current, error)
                self.best_guess["currents"][key][var] = {"value":current, "error":error}
                parameter_index += 1
                if print_parameters:
                    print var, current, "+-", error,
        if print_parameters:
            print

    def set_lattice(self):
        for key in self.optimisation_currents:
            # fractional change in magnet current
            new_factor = self.best_guess["currents"][key]["j"]["value"]
            for module in self.optimisation_currents[key]:
                module_name = module["name"]
                if module["j"] == None:
                    module["j"] = self.murgler.get_scale_factor(module_name)
                if module["j"] == None:
                    raise RuntimeError("Failed to get current for module "+module_name)
                new_scale_factor = module["j"]*new_factor
                self.murgler.set_scale_factor(module_name, new_scale_factor)
        for module_name in self.optimisation_axes:
            z = self.optimisation_axes[module_name]["z_pos"]
            alignment = self.best_guess["axes"][module_name]
            x = alignment["x"]["value"]
            y = alignment["y"]["value"]
            xp = alignment["xp"]["value"]
            yp = alignment["yp"]["value"]
            self.murgler.misalign_mice_module(x, xp, y, yp, z, module_name)
        self.murgler.update()

    def do_tracking(self):
        self.extrapolate_list = []
        print "Doing tracking for", self.track_source, "tracking", self.direction, "to", self.residuals_source
        for event_bin in self.event_list: # momentum bins
            extrapolate = ExtrapolateTrackPoints(self.config, self.config_anal, self.data_loader, False)
            self.extrapolate_list.append(extrapolate)
            for event in event_bin:
                event = copy.deepcopy(event)
                try:
                    extrapolate.extrapolate_event(event, self.track_source, self.direction)
                except ValueError:
                    pass #sys.excepthook(*sys.exc_info())
                for detector in [self.residuals_source]:
                    for axis in self.residuals_variables:
                        extrapolate.append_residual(event, detector, axis)

    def calculate_score_1(self):
        """
        Calculate score based on sum of square of deviation of all track positions
        """
        chi2 = 0.
        residual_dicts = self.extrapolate.residual_dicts
        for detector in [self.residuals_source]:
            for axis in self.residuals_variables:
                residual = residual_dicts[detector][axis]
                for value in residual:
                    chi2 += value**2/len(residual)
        self.best_guess["score"] = {"score":chi2}
        return chi2

    def calculate_score_2(self):
        """
        Calculate score based on sum of square of mean deviation of track positions
        """
        chi2 = 0.
        results = []
        for i, extrapolate in enumerate(self.extrapolate_list):
            residual_dicts = extrapolate.residual_dicts
            means = {"p_bin":self.config_anal["p_bins"][i]}
            for detector in [self.residuals_source]:
                print "   ", detector, means["p_bin"],
                if detector not in residual_dicts:
                    print "Failed to find event residuals for detector", detector, "for momenta", means["p_bin"]
                    continue
                for axis in self.residuals_variables:
                    means[axis] = numpy.mean(residual_dicts[detector][axis])
                    means[axis+"_err"] = numpy.std(residual_dicts[detector][axis])
                    means["n"] = len(residual_dicts[detector][axis])
                    chi2 += means[axis]**2
                    print axis+":", round(means[axis], 2), "std:", round(means[axis+"_err"], 2), "n:", means["n"],
                print
            results.append(means)
        score = {"score":chi2, "results":results}
        self.best_guess["score"] = score
        return chi2

    def print_output(self):
        if self.fout == None:
            self.fout = open(self.config_anal["plot_dir"]+"/alignment.json", "w")
        print >> self.fout, json.dumps(self.best_guess)
        self.fout.flush()

    @staticmethod
    def run_once(print_field = False, print_parameters = True, print_score = True):
        self = ResidualFitter.fitter
        # update the best_guess
        self.update_from_minuit(print_parameters)
        self.set_lattice()
        self.do_tracking()
        if print_field:
            print maus_cpp.field.str(True)
        score = self.calculate_score_2()
        self.print_output()
        return score


    @staticmethod
    def score_function(n_pars, pars, score, buff, err):
        """
        Calculate the transfer matrix and compare with the measured tranfer 
        matrix; use to calculate a chi2
        """
        self = ResidualFitter.fitter
        score[0] = ResidualFitter.fitter.run_once()
        self._iteration += 1
        print "Iteration", self._iteration, "Score", score[0]
        print

    fitter = None

