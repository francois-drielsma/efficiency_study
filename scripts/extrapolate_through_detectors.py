import copy
import sys
import os
import site
import json
import datetime
import bisect
import libxml2

import numpy
import ROOT
import xboa.common
from xboa.hit import Hit
from xboa.bunch import Bunch
import Configuration
import maus_cpp.global_error_tracking
import maus_cpp.field
import xboa.common

import scripts.utilities

class ExtrapolateTrackPoints(object):
    def __init__(self, config, config_anal, data_loader, verbose = True):
        self.config = config
        self.config_anal = config_anal
        self.current_space_points = [] # list of space points in current event
        self.fitted_track_points = [] # list of fitted space points in current event
        self.best_guess = None
        self.iteration = 0
        self.mu_mass = xboa.common.pdg_pid_to_mass[13]
        self.setup_maus_config()
        self.residual_dicts = {}
        self.miss_lists = {}
        self.weighted_residual_dicts = {}
        self.data_loader = data_loader
        self.csv_out = None
        self.csv_headers = ["spill", "event_number", "t", "x", "y", "z", "px", "py", "pz"]
        self.apertures = sorted(config.apertures+config.detectors)
        if verbose:
            print "Found apertures like"
            print "z".rjust(8), "radius".rjust(8), "name"
            for z, radius, name in self.apertures:
                print str(z).rjust(10), str(radius).rjust(8), name.rjust(12),
                if radius != None:
                    print "o"+"|".rjust(int(radius/10))
                else:
                    print "o"

    def setup_maus_config(self):
        self.tracking = maus_cpp.global_error_tracking.GlobalErrorTracking()
        self.tracking.set_min_step_size(self.config.global_min_step_size)
        self.tracking.set_max_step_size(self.config.global_max_step_size)
        self.tracking.set_energy_loss_model("bethe_bloch_forwards")
        self.tracking.set_scattering_model("moliere_backwards")
        self.tracking.set_geometry_model("axial_lookup")
        for material in ["AIR", "POLYSTYRENE", "POLYCARBONATE", "POLYVINYL_TOLUENE",
                        "POLYVINYL_ACETATE", "AEROGEL_112a", "AEROGEL_107a",
                        "Al", "CELLULOSE_CELLOPHANE", "POLYVINYL_CHLORIDE_LOWD",
                        "POLYVINYL_CHLORIDE", "He", "LITHIUM_HYDRIDE"]:
            # disable apparent ALUMINUM between 17317.1 and 17230.1
            maus_cpp.global_error_tracking.enable_material(material)

    def get_point(self, event, detector):
        if len(event) == 0:
            raise ValueError("Event was empty")
        # t can be None if time_offset fails to find the reference time detector
        # enforced here so that we always plot the same sample
        for rec_point in event:
            if rec_point["detector"] == detector and \
               rec_point["hit"]["t"] != None:
                break
        if rec_point["detector"] != detector:
            detectors = str([rec["detector"] for rec in event])
            raise ValueError("Failed to find rec point for "+str(detector)+" from "+ detectors)
        return rec_point


    def make_rec(self, value, ellipse, detector, reference):
        hit = Hit.new_from_dict( {
            "t":value[0],
            "x":value[1],
            "y":value[2],
            "z":value[3],
            "px":value[5],
            "py":value[6],
            "pz":value[7],
            "pid":reference["hit"]["pid"],
            "charge":reference["hit"]["charge"],
            "mass":xboa.common.pdg_pid_to_mass[13],
            "spill":reference["hit"]["spill"],
            "event_number":reference["hit"]["event_number"],
            "particle_number":reference["hit"]["particle_number"],
          }, "energy"
        )
        rec_point = {
            "detector":self.global_key(detector),
            "covariance":ellipse,
            "hit":hit
        }
        return rec_point

    def is_global(self, detector):
        return len(detector) > 6 and detector[:7] == "global_"

    def global_key(self, detector):
        return "global_"+detector

    def get_delta_tof(self, event, det_1, det_2):
        """Get the difference between tof for extrapolation and recon
           (recon_tof_2 - recon_tof_1) - (extrap_tof_2 - extrap_tof_1)
        """
        times = [None, None, None, None]
        dets = (det_1, det_2, self.global_key(det_1), self.global_key(det_2))
        # find times for each detector type in dets and fill into times
        for hit in event["data"]:
            for i, det in enumerate(dets):
                if hit["detector"] == det:
                    times[i] = hit["hit"]["t"]
        # return None or the difference between tofs
        if None in times:
            return None
        else:
            delta_tof = (times[1] - times[0]) - (times[3] - times[2])
            return delta_tof


    def time_offset(self, detector_hits):
        """subtract time from specified detector or set to None if not found"""
        reco_delta_t = None
        global_delta_t = None
        for hit in detector_hits:
            if hit["detector"] == self.config.time_from:
                reco_delta_t = hit["hit"]["t"]
            elif hit["detector"] == self.global_key(self.config.time_from):
                global_delta_t = hit["hit"]["t"]
        if reco_delta_t == None or global_delta_t == None: # time_from detector was not found; default to None
            for hit in detector_hits:
                hit["hit"]["t"] = 1e9 # time must be a float - so this is an error code (urk)!
            return
        delta_t = -global_delta_t + reco_delta_t
        for hit in detector_hits:
            if self.is_global(hit["detector"]):
                hit["hit"]["t"] += delta_t

    def extrapolate_event(self, event, reference_detector = "tku_tp", direction = "both"):
        detector_hits = event["data"]
        detector_hits = sorted(detector_hits, key = lambda tp: tp["hit"]["z"]) # check that it is sorted
        tku_rec = self.get_point(detector_hits, reference_detector)
        energy = (tku_rec["hit"]["px"]**2 + tku_rec["hit"]["py"]**2 + tku_rec["hit"]["pz"]**2 + self.mu_mass**2)**0.5
        seed = [0.,     tku_rec["hit"]["x"],  tku_rec["hit"]["y"],  tku_rec["hit"]["z"],
                tku_rec["hit"]["energy"], tku_rec["hit"]["px"], tku_rec["hit"]["py"], tku_rec["hit"]["pz"],]
        ellipse = tku_rec["covariance"]
        void = [
            [0., 0.,  0., 0.,  0., 0.],
            [0., 0.1, 0., 0.,  0., 0.],
            [0., 0., 0.1, 0.,  0., 0.],
            [0., 0.,  0., 16., 0., 0.],
            [0., 0.,  0., 0.,  4., 0.],
            [0., 0.,  0., 0.,  0., 4.],
        ]
        first = bisect.bisect_left(self.apertures, (seed[3], None, None))
        # walking upstream from reference_detector
        point, error = seed, ellipse
        #self.tracking.set_tracking_model("em_forwards_dynamic")
        self.tracking.set_energy_loss_model("bethe_bloch_forwards")

        if direction == "both" or direction == "upstream":
            self.tracking.set_scattering_model("moliere_backwards") 
            try:
                for z_pos, radius, name in reversed(self.apertures[:first]):
                    point, error = self.tracking.propagate_errors(point, error, z_pos)
                    hit = self.make_rec(point, error, name, tku_rec)
                    if radius != None and hit["hit"]["r"] > radius:
                        event["apertures_us"].append(name)
                    detector_hits.append(hit)
            except RuntimeError:
                pass #sys.excepthook(*sys.exc_info())

        if direction == "both" or direction == "downstream":
            self.tracking.set_scattering_model("moliere_forwards")
            # walking downstream from reference_detector
            point, error = seed, ellipse
            try:
                for z_pos, radius, name in self.apertures[first:]:
                    point, error = self.tracking.propagate_errors(point, error, z_pos)
                    hit = self.make_rec(point, error, name, tku_rec)
                    if radius != None and hit["hit"]["r"] > radius:
                        event["apertures_ds"].append(name)
                    detector_hits.append(hit)
            except RuntimeError:
                pass # sys.excepthook(*sys.exc_info())

        # use time from detector defined in datacards; move all times relatively
        self.time_offset(detector_hits)
        event["data"] = detector_hits

    def extrapolation_cuts(self, event):
        delta_tof01 = self.get_delta_tof(event, "tof0", "tof1")
        lower = self.config_anal["delta_tof01_lower"]
        upper = self.config_anal["delta_tof01_upper"]
        event["delta_tof01"] = delta_tof01
        event["will_cut"]["delta_tof01"] = delta_tof01 < lower or \
                                           delta_tof01 > upper
        delta_tof12 = self.get_delta_tof(event, "tof1", "tof2")
        lower = self.config_anal["delta_tof12_lower"]
        upper = self.config_anal["delta_tof12_upper"]
        event["delta_tof12"] = delta_tof12
        event["will_cut"]["delta_tof12"] = delta_tof12 < lower or \
                                           delta_tof12 > upper
        event["will_cut"]["aperture_us"] = len(event["apertures_us"]) > 0
        event["will_cut"]["aperture_ds"] = len(event["apertures_ds"]) > 0
        return event

    def alt_min_max(self, float_list, n_sigma):
        length = len(float_list)+1
        while len(float_list) < length:
            mean = numpy.mean(float_list)
            sigma = numpy.std(float_list)
            # e.g. detector has sigma t exactly == 0 due to the time_offset algorithm
            if sigma == 0.:
                return [-1., 1.]
            min_max = [
                mean-n_sigma*sigma,
                mean+n_sigma*sigma,
            ]
            length = len(float_list)
            float_list = [x for x in float_list if x < min_max[1] and x > min_max[0]]
        return min_max

    def append_misses(self, event, detector):
        if detector not in self.miss_lists:
            self.miss_lists[detector] = []
        event = event["data"]
        try:
            det_rec = self.get_point(event, detector)
            return
        except ValueError:
            # we successfully extrapolated the track, but it is not observed in
            # the detector... I wonder why?
            pass

        try:
            glob_rec = self.get_point(event, self.global_key(detector))["hit"]
        except ValueError:
            # we failed to extrapolate a track, nothing more to do
            return
        self.miss_lists[detector].append(glob_rec)

    def clear_residuals(self):
        self.residual_dicts = {}

    def append_residual(self, event, detector, axis_list):
        event = event["data"]
        if detector not in self.residual_dicts:
            self.residual_dicts[detector] = {}
            self.weighted_residual_dicts[detector] = {}
        for axis in axis_list:
            if axis not in self.residual_dicts[detector]:
                self.residual_dicts[detector][axis] = []
                self.weighted_residual_dicts[detector][axis] = []
        
        glob_rec, det_rec = None, None
        try:
            glob_rec = self.get_point(event, self.global_key(detector))
        except ValueError:
            # we failed to extrapolate a track, nothing more to do
            # sys.excepthook(*sys.exc_info())
            return

        try:
            det_rec = self.get_point(event, detector)
        except ValueError:
            # we successfully extrapolated the track, but it is not observed in
            # the detector... I wonder why?
            # sys.excepthook(*sys.exc_info())
            return

        # we have an extrapolated track and a detector hit. How close is the
        # extrapolated track to the measured particle?
        for axis in axis_list:
            residual_list = self.residual_dicts[detector][axis]
            residual = det_rec["hit"][axis]-glob_rec["hit"][axis]
            residual_list.append(residual)

        # gymnastics to ensure that weighted residuals plots have same number of events
        # at a given detector
        weighted_residual_calc = []
        for axis in axis_list:
            if axis in self.cov_key_list and "covariance" in det_rec and "covariance"  in glob_rec:
                key = self.cov_key_list.index(axis)
                try:
                    sigma = (glob_rec["covariance"][key][key]+det_rec["covariance"][key][key])**0.5 # -ve sqrt exception?
                    if sigma < 1e-12:
                        raise ValueError("Zero RMS")
                    weighted_residual_calc.append(residual/sigma)
                except ValueError:
                    #print axis, "glob:", glob_rec["covariance"][key][key], detector, det_rec["covariance"][key][key]
                    return
            else:
                return

        for i, axis in enumerate(axis_list):
            weighted_residual_list = self.weighted_residual_dicts[detector][axis]
            weighted_residual_list.append(weighted_residual_calc[i])

    def print_canvas(self, canvas, name):
        canvas.Update()
        plot_name = os.path.join(self.config_anal["plot_dir"], name.replace(' ', '_'))
        for format in "png", "pdf", "root":
            canvas.Print(plot_name+"."+format)

    def get_geometry(self):
        path_to_info_file = self.config.info_file
        info = libxml2.parseFile(path_to_info_file)
        path = "gdml/MICE_Information/Configuration_Information/GeometryID"
        geo_id = info.xpathEval(path)[0].prop("value")
        info.freeDoc()
        return geo_id

    text_boxes = []
    def get_text_box(self, residual_list, fit = None):
        text_box = ROOT.TPaveText(0.6, 0.4, 0.9, 0.9, "NDC")
        text_box.SetFillColor(0)
        text_box.SetBorderSize(0)
        text_box.SetTextSize(0.04)
        text_box.SetTextAlign(12)
        text_box.SetTextSize(0.03)
        if fit == None:
            mean = round(numpy.mean(residual_list), 2)
            sigma = round(numpy.std(residual_list), 2)
        else:
            mean = round(fit.GetParameter(1), 2)
            sigma = round(fit.GetParameter(2), 2)
        text_box.AddText("MICE Internal")
        text_box.AddText(self.config_anal["name"])
        text_box.AddText("Recon: "+self.data_loader.maus_version)
        text_box.AddText("Geometry: "+str(self.get_geometry()))
        text_box.AddText("All events (black)")
        text_box.AddText("Number: "+str(len(residual_list)))
        text_box.AddText("Mean:   "+str(mean))
        text_box.AddText("Std:    "+str(sigma))
        text_box.SetBorderSize(1)
        text_box.Draw()
        self.text_boxes.append(text_box)
        return text_box


    def plot_residuals(self, detector, axis):
        if detector not in self.residual_dicts:
            print "No residuals for", detector
            return
        residual_list = self.residual_dicts[detector][axis]
        weighted_residual_list = self.weighted_residual_dicts[detector][axis]
        if False: #axis == "t" and detector == "tof1": #axis == "pz":
            print "Residuals", residual_list
            print "Weighted Residuals", weighted_residual_list
            print "Det", detector, "axis", axis
        if len(residual_list) == 0:
            print "No residuals were found for "+detector+" "+axis+"? Something screwy?"
            return

        nbins = self.config.residuals_plots_nbins
        name = "residuals - "+detector+" "+axis
        label = "Res("+axis+") ["+self.units[axis]+"]"
        min_max = self.alt_min_max(residual_list, 10)
        canvas = xboa.common.make_root_canvas(name)
        hist = xboa.common.make_root_histogram(name, residual_list, label, nbins, xmin = min_max[0], xmax = min_max[1])
        hist.SetTitle(detector+": "+axis)
        fit = scripts.utilities.fit_peak(hist, nsigma=1)
        hist.Draw()
        self.get_text_box(residual_list, fit)
        canvas.Update()
        self.print_canvas(canvas, name)

        name = "normalised "+name
        canvas = xboa.common.make_root_canvas(name)
        min_max = self.alt_min_max(weighted_residual_list, 10)
        if len(weighted_residual_list) > 0:
            hist = xboa.common.make_root_histogram(name, weighted_residual_list, "Res("+axis+")/#sigma("+axis+")", nbins, xmin = min_max[0], xmax = min_max[1])
            hist.SetTitle(detector+": "+axis)
            fit = scripts.utilities.fit_peak(hist, nsigma=1)
            hist.Draw()
            self.get_text_box(weighted_residual_list, fit)
            self.print_canvas(canvas, name)

    def get_miss_text_box(self, bunch, detector):
        text_box = ROOT.TPaveText(0.6, 0.4, 0.9, 0.9, "NDC")
        text_box.SetFillColor(0)
        text_box.SetBorderSize(0)
        text_box.SetTextSize(0.04)
        text_box.SetTextAlign(12)
        text_box.SetTextSize(0.03)
        text_box.AddText("MICE Internal")
        text_box.AddText("Recon: "+self.data_loader.maus_version)
        text_box.AddText("Geometry: "+str(self.get_geometry()))
        text_box.AddText("All events (black)")
        text_box.AddText("Number missed: "+str(len(bunch)))
        n_rec = len(self.residual_dicts[detector].values()[0])
        text_box.AddText("Number recorded: "+str(n_rec))
        text_box.Draw()
        self.text_boxes.append(text_box)
        return text_box

    def draw_box(self, x0, y0, x1, y1):
        x = [x0, x0, x1, x1, x0]
        y = [y0, y1, y1, y0, y0]
        hist, graph = xboa.common.make_root_graph("", x, y, "", "")
        graph.Draw()
        return graph

    def plot_misses(self):
        for detector in "tof0", "tof1", "tof2":
            if detector not in self.miss_lists or len(self.miss_lists[detector]) == 0:
                print "No misses for", detector
                continue
            print detector
            bunch = Bunch.new_from_hits(self.miss_lists[detector])
            canvas, hist = bunch.root_histogram("x", "mm", nx_bins=20, xmin=-500., xmax=+500.)
            hist.SetTitle("Misses - "+detector)
            self.get_miss_text_box(bunch, detector)
            self.print_canvas(canvas, "misses_"+detector+"_x")
            canvas, hist = bunch.root_histogram("y", "mm", nx_bins=20, xmin=-500., xmax=+500.)
            hist.SetTitle("Misses - "+detector)
            self.get_miss_text_box(bunch, detector)
            self.print_canvas(canvas, "misses_"+detector+"_y")
            canvas, hist = bunch.root_histogram("x", "mm", "y", "mm", nx_bins=20, ny_bins=20, xmin=-500., xmax=+500., ymin=-500., ymax=+500.)
            hist.SetTitle("Misses - "+detector)
            self.get_miss_text_box(bunch, detector)
            self.print_canvas(canvas, "misses_"+detector+"_x-y")
        for detector in "tkd_tp",:
            print self.miss_lists.keys()
            if detector not in self.miss_lists or len(self.miss_lists[detector]) == 0:
                print "No misses for", detector
                continue
            print detector
            bunch = Bunch.new_from_hits(self.miss_lists[detector])
            canvas, hist = bunch.root_histogram("r", "mm", nx_bins=10, xmin=0., xmax=+300.)
            hist.SetTitle("Misses - "+detector)
            self.get_miss_text_box(bunch, detector)
            self.print_canvas(canvas, "misses_"+detector+"_r")
            canvas, hist = bunch.root_histogram("pt", "MeV/c", nx_bins=10, xmin=0., xmax=+100.)
            hist.SetTitle("Misses - "+detector)
            self.get_miss_text_box(bunch, detector)
            self.print_canvas(canvas, "misses_"+detector+"_pt")
            canvas, hist = bunch.root_histogram("r", "mm", "pt", "MeV/c", nx_bins=50, ny_bins=50, xmin=0., xmax=+200., ymin=0., ymax=+100.)
            hist.SetTitle("Misses - "+detector)
            self.get_miss_text_box(bunch, detector)
            self.print_canvas(canvas, "misses_"+detector+"_r-pt")
            canvas, hist = bunch.root_histogram("x", "mm", "y", "mm", nx_bins=10, ny_bins=10, xmin=-300., xmax=+300., ymin=-300., ymax=+300.)
            hist.SetTitle("Misses - "+detector)
            self.get_miss_text_box(bunch, detector)
            self.print_canvas(canvas, "misses_"+detector+"_x-y")

    units = {"t":"ns", "x":"mm", "y":"mm", "z":"mm", "energy":"MeV", "px":"MeV/c", "py":"MeV/c", "pz":"MeV/c"}
    cov_key_list = ["t", "x", "y", "energy", "px", "py"]
    mu_mass = xboa.common.pdg_pid_to_mass[13]

    def print_hits_csv(self):
        if "csv_output_filename" not in self.config_anal or self.config_anal["csv_output_filename"] == None:
            return
        if self.csv_out == None:
            self.csv_out = open(self.config_anal["csv_output_filename"], "w")
            for head in self.csv_headers:
                if head in self.units:
                    units = " ["+self.units[head]+"]"
                else:
                    units = ""
                print >> self.csv_out, str(head+units).ljust(15),
            print >> self.csv_out
        print "    writing data...",
        sys.stdout.flush()
        for event in self.data_loader.events:
            for detector in self.config_anal["csv_output_detectors"]:
                for hit in event["data"]:
                    if hit["detector"] != self.global_key(detector):
                        continue
                    for head in self.csv_headers:
                        if head in hit["hit"].get_variables():
                            print >> self.csv_out, str(hit["hit"][head]).ljust(15),
                        elif head in event.keys():
                            print >> self.csv_out, str(event[head]).ljust(15),
                        else:
                            raise KeyError("Couldn't find key "+head+" from event keys: "+str(event.keys())+" hit keys: "+str(hit["hit"].get_variables()))
                    print >> self.csv_out
                    break
        self.csv_out.flush()
        print "written"

    root_objects = []

    @classmethod
    def is_cut(cls, event, config):
        for cut in config.extrapolation_cuts:
            if config.extrapolation_cuts[cut] and event["will_cut"][cut]:
                return True
        return False

def do_plots(extrapolate):
        xboa.common.clear_root()
        extrapolate.plot_misses()
        for detector in "tof0", "tof1", "tof2":
            for axis in "x", "y", "t":
                extrapolate.plot_residuals(detector, axis)      
        for detector in "tkd_tp",:
            for axis in "x", "y", "px", "py", "pz":
                extrapolate.plot_residuals(detector, axis)      

def do_extrapolation(config, config_anal, data_loader):
    print "Doing extrapolation"
    index = 0
    pass_first_cuts = 0
    pass_extrapolate = 0
    pass_second_cuts = 0
    step = 500
    will_do_plots = False
    extrapolate = ExtrapolateTrackPoints(config, config_anal, data_loader)
    print "using model:"
    print extrapolate.tracking.geometry_model_string()
    print "starting extrapolation..."
    sys.stdout.flush()
    while index < len(data_loader.events):
        try:
            while will_do_plots == False:
                index += 1
                event = data_loader.events[index]
                pass_first_cuts += 1
                try:
                    extrapolate.extrapolate_event(event, config_anal["extrapolation_source"]) #, "tkd_tp")  
                    extrapolate.extrapolation_cuts(event)
                    pass_extrapolate += 1
                except ValueError:
                    pass #sys.excepthook(*sys.exc_info())
                if extrapolate.is_cut(event, config):
                    continue
                pass_second_cuts += 1
                for detector in "tof1", "tof0", "tof2":
                    extrapolate.append_misses(event, detector)
                    extrapolate.append_residual(event, detector, ["x", "y", "t"])
                for detector in "tkd_tp",:
                    extrapolate.append_misses(event, detector)
                    extrapolate.append_residual(event, detector, [ "x", "y", "px", "py", "pz"])
                will_do_plots = pass_second_cuts % step == 0
        except IndexError:
            index = len(data_loader.events) # force to finish
            sys.excepthook(*sys.exc_info())
        print "event", index, "of", len(data_loader.events), "with", pass_first_cuts, \
              "passing first cuts,", pass_extrapolate, "extrapolating okay and", \
              pass_second_cuts, "passing second cuts"
        do_plots(extrapolate)
        will_do_plots = False
    extrapolate.print_hits_csv()

