import sys
import glob
import itertools
import math
import time

import ROOT
import libMausCpp

import xboa.common
from xboa.hit import Hit

# Fix tracker cuts!

class DataLoader(object):
    def __init__(self, config, analysis_index):
        self.config = config
        self.config_anal = config.analyses[analysis_index]
        self.detector_count = {}
        self.spill_count = 0
        self.event_count = 0
        self.accepted_count = 0
        self.events = []
        self.maus_version = ""
        self.run_numbers = set([])
        self.this_run = 0
        self.this_spill = 0
        self.this_event = 0
        self.cuts = {}
        self.event_id = 0
        self.det_pos = self.get_det_pos()

    def get_det_pos(self):
        det_pos = {}
        for det in self.config.detectors:
            det_pos[det[2]] = det[0]
        return det_pos

    def load_data(self, min_spill, max_spill):
        if min_spill >= max_spill:
            return
        self.spill_count = 0
        file_name_list = []
        for fname in self.config_anal["reco_files"]:
            file_name_list += glob.glob(fname)
        file_name_list = sorted(file_name_list)
        print "Found", file_name_list, "files"
        if len(file_name_list) == 0:
            raise IOError("Couldn't find files for reco_files "+str(self.config_anal["reco_files"]))
        for file_name in file_name_list:
            print "Loading ROOT file", file_name
            
            sys.stdout.flush()
            root_file = ROOT.TFile(file_name, "READ") # pylint: disable = E1101

            data = ROOT.MAUS.Data() # pylint: disable = E1101
            tree = root_file.Get("Spill")
            try:
                tree.SetBranchAddress("data", data)
            except AttributeError:
                continue
            print tree.GetEntries(), "spills..."
            events_per_file = [] # worried that list.append is consuming cpu
            new_t = time.time()
            for i in range(tree.GetEntries()):
                events_per_spill = [] # worried that list.append is consuming cpu
                if i % 100 == 0:
                    old_t = new_t
                    new_t = time.time()
                    print i, "("+str(round(new_t - old_t, 2))+")", 
                    sys.stdout.flush()
                tree.GetEntry(i)
                spill = data.GetSpill()
                self.this_run = spill.GetRunNumber()
                self.this_spill = spill.GetSpillNumber()
                if i == 0:
                    self.run_numbers.add(spill.GetRunNumber())
                if spill.GetDaqEventType() == "physics_event":
                    if self.config.number_of_spills != None and \
                       self.spill_count > self.config.number_of_spills:
                        #file_name_list = []
                        break
                    self.spill_count += 1
                    for ev_number, reco_event in enumerate(spill.GetReconEvents()):
                        self.this_event = reco_event.GetPartEventNumber()
                        try:
                            event = self.load_reco_event(reco_event)
                        except ValueError:
                            #sys.excepthook(*sys.exc_info())
                            print "spill", spill.GetSpillNumber(), "particle_number", reco_event.GetPartEventNumber()
                        except ZeroDivisionError:
                            pass
                        if event == None: # missing TOF1 - not considered further
                            continue 
                        self.event_count += 1
                        event["spill"] = spill.GetSpillNumber()
                        events_per_spill.append(event)
                        event["event_number"] = ev_number
                        if self.config_anal["do_mc"] and len(spill.GetMCEvents()) > ev_number:
                            self.load_mc_event(event, spill.GetMCEvents()[ev_number])
                events_per_file += events_per_spill
            self.events += events_per_file
            print "\n  ...loaded", self.spill_count, "'physics_event' spills and", self.event_count, "events"
            sys.stdout.flush()
        print "Finalising cuts"
        self.update_cuts()

    def load_mc_event(self, loaded_event, mc_event):
        # mc load functions return a list of hits
        primary = mc_event.GetPrimary()
        primary_loaded = self.load_primary(primary)
        track_vector = mc_event.GetTracks()
        track_loaded = self.load_tracks(track_vector)
        virtual_vector = mc_event.GetVirtualHits()
        virtual_loaded = self.load_virtuals(virtual_vector)
        tof_hit_vector = mc_event.GetTOFHits()
        tof_hit_loaded = self.load_tof_mc_hits(tof_hit_vector)
        temp_event = sorted(primary_loaded+track_loaded+virtual_loaded+tof_hit_loaded, key = lambda tp: tp['hit']['z'])
        loaded_event["data"] += temp_event
        
    def load_tracks(self, track_vector):
        loaded_track_vector = []
        for track in track_vector:
            if abs(track.GetParticleId()) == 211 and track.GetKillReason() != "":
                print "Pion killed because", track.GetKillReason()
            hit = {}
            hit["x"] = track.GetInitialPosition().x()
            hit["y"] = track.GetInitialPosition().y()
            hit["z"] = track.GetInitialPosition().z()
            hit["px"] = track.GetInitialMomentum().x()
            hit["py"] = track.GetInitialMomentum().y()
            hit["pz"] = track.GetInitialMomentum().z()
            hit["pid"] = track.GetParticleId()
            try:
                hit["mass"] = xboa.common.pdg_pid_to_mass[abs(hit["pid"])]
            except KeyError:
                hit["mass"] = 0.
            try:
                hit["charge"] = xboa.common.pdg_pid_to_charge[hit["charge"]]
            except KeyError:
                hit["charge"] = 0.

            hit = {}
            hit["x"] = track.GetFinalPosition().x()
            hit["y"] = track.GetFinalPosition().y()
            hit["z"] = track.GetFinalPosition().z()
            hit["px"] = track.GetFinalMomentum().x()
            hit["py"] = track.GetFinalMomentum().y()
            hit["pz"] = track.GetFinalMomentum().z()
            hit["pid"] = track.GetParticleId()
            try:
                hit["mass"] = xboa.common.pdg_pid_to_mass[abs(hit["pid"])]
            except KeyError:
                hit["mass"] = 0.
            try:
                hit["charge"] = xboa.common.pdg_pid_to_charge[hit["charge"]]
            except KeyError:
                hit["charge"] = 0.
       
            loaded_track_initial = {
                "hit":Hit.new_from_dict(hit, "energy"),
                "detector":"mc_track_initial"
            }
            loaded_track_final = {
                "hit":Hit.new_from_dict(hit, "energy"),
                "detector":"mc_track_final",
            }
            loaded_track_vector += [loaded_track_initial, loaded_track_final]
        return loaded_track_vector

    def load_primary(self, primary):
        hit = {}
        hit["x"] = primary.GetPosition().x()
        hit["y"] = primary.GetPosition().y()
        hit["z"] = primary.GetPosition().z()
        hit["px"] = primary.GetMomentum().x()
        hit["py"] = primary.GetMomentum().y()
        hit["pz"] = primary.GetMomentum().z()
        hit["pid"] = primary.GetParticleId()
        hit["energy"] = primary.GetEnergy()
        hit["t"] = primary.GetTime()
        try:
            hit["mass"] = xboa.common.pdg_pid_to_mass[abs(hit["pid"])]
        except KeyError:
            hit["mass"] = 0.
        try:
            hit["charge"] = xboa.common.pdg_pid_to_charge[hit["charge"]]
        except KeyError:
            hit["charge"] = 0.
        loaded_primary = {
            "hit":Hit.new_from_dict(hit),
            "detector":"mc_primary"
        }
        return [loaded_primary]

    def load_virtuals(self, virtual_vector):
        loaded_virtual_vector = [None]*len(virtual_vector)
        for i, virtual_hit in enumerate(virtual_vector):
            hit = xboa.hit.Hit()
            hit["x"] = virtual_hit.GetPosition().x()
            hit["y"] = virtual_hit.GetPosition().y()
            hit["z"] = virtual_hit.GetPosition().z()
            hit["px"] = virtual_hit.GetMomentum().x()
            hit["py"] = virtual_hit.GetMomentum().y()
            hit["pz"] = virtual_hit.GetMomentum().z()
            hit["bx"] = virtual_hit.GetBField().x()
            hit["by"] = virtual_hit.GetBField().y()
            hit["bz"] = virtual_hit.GetBField().z()
            hit["ex"] = virtual_hit.GetEField().x()
            hit["ey"] = virtual_hit.GetEField().y()
            hit["ez"] = virtual_hit.GetEField().z()
            hit["pid"] = virtual_hit.GetParticleId()
            hit["station"] = virtual_hit.GetStationId()
            hit["particle_number"] = virtual_hit.GetTrackId()
            hit["t"] = virtual_hit.GetTime()
            hit["mass"] = virtual_hit.GetMass()
            hit["charge"] = virtual_hit.GetCharge()
            hit.mass_shell_condition('energy')
            loaded_virtual = {
                "hit":hit,
                "detector":"mc_virtual_"+str(hit["station"])
            }
            loaded_virtual_vector[i] = loaded_virtual
        return loaded_virtual_vector

    def load_tof_mc_hits(self, tof_mc_vector):
        loaded_tof_mc_vector = [None]*len(tof_mc_vector)
        for i, tof_mc_hit in enumerate(tof_mc_vector):
            hit = xboa.hit.Hit()
            hit["x"] = tof_mc_hit.GetPosition().x()
            hit["y"] = tof_mc_hit.GetPosition().y()
            hit["z"] = tof_mc_hit.GetPosition().z()
            hit["px"] = tof_mc_hit.GetMomentum().x()
            hit["py"] = tof_mc_hit.GetMomentum().y()
            hit["pz"] = tof_mc_hit.GetMomentum().z()
            hit["pid"] = tof_mc_hit.GetParticleId()
            hit["station"] = tof_mc_hit.GetChannelId().GetStation()
            hit["particle_number"] = tof_mc_hit.GetTrackId()
            hit["t"] = tof_mc_hit.GetTime()
            hit["energy"] = tof_mc_hit.GetEnergy()
            hit["e_dep"] = tof_mc_hit.GetEnergyDeposited()
            try:
                hit["mass"] = xboa.common.pdg_pid_to_mass[abs(tof_mc_hit.GetParticleId())]
            except KeyError:
                hit["mass"] = (hit["energy"]**2-hit["p"]**2)**0.5
            loaded_tof_mc = {
                "hit":hit,
                "detector":"mc_tof_"+str(hit["station"]),
            }
            loaded_tof_mc_vector[i] = loaded_tof_mc
        return loaded_tof_mc_vector

    def load_tof_sp(self, tof_sp, station):
        xerr = tof_sp.GetGlobalPosXErr()
        yerr = tof_sp.GetGlobalPosYErr()
        cov = [
            [0.0025, 0.,     0., 0., 0., 0.],
            [0.,     xerr,   0., 0., 0., 0.],
            [0.,     0.,   yerr, 0., 0., 0.],
            [0.,     0., 0., 0., 0., 0.],
            [0.,     0., 0., 0., 0., 0.],
            [0.,     0., 0., 0., 0., 0.],
        ]
        sp_dict = {
            "x":tof_sp.GetGlobalPosX(),
            "y":tof_sp.GetGlobalPosY(),
            "z":tof_sp.GetGlobalPosZ(),
            "t":tof_sp.GetTime(),
        }
        loaded_sp = {
            "hit":Hit.new_from_dict(sp_dict),
            "detector":station,
            "covariance":cov,
        }
        return loaded_sp


    def load_tof_event(self, tof_event):
        space_points = tof_event.GetTOFEventSpacePoint()
        tof_sp_list = []
        #if len(space_points.GetTOF2SpacePointArray()) > 0:
        #    print len(space_points.GetTOF0SpacePointArray()), len(space_points.GetTOF1SpacePointArray()), len(space_points.GetTOF2SpacePointArray())
        for tof_sp in space_points.GetTOF0SpacePointArray():
            tof_sp_list.append(self.load_tof_sp(tof_sp, "tof0"))
        for tof_sp in space_points.GetTOF1SpacePointArray():
            tof_sp_list.append(self.load_tof_sp(tof_sp, "tof1"))
        for tof_sp in space_points.GetTOF2SpacePointArray():
            tof_sp_list.append(self.load_tof_sp(tof_sp, "tof2"))
        detectors = [x["detector"] for x in tof_sp_list]

        tof01_cut = False
        if "tof1" in detectors and "tof0" in detectors:
            tof01 = tof_sp_list[detectors.index("tof1")]["hit"]["t"] - tof_sp_list[detectors.index("tof0")]["hit"]["t"]
            if tof01 > self.config_anal["tof01_cut_high"] or \
               tof01 < self.config_anal["tof01_cut_low"]:
                tof01_cut = True
        else:
            tof01_cut = True

        tof12_cut = False
        if "tof2" in detectors and "tof1" in detectors:
            tof12 = tof_sp_list[detectors.index("tof2")]["hit"]["t"] - tof_sp_list[detectors.index("tof1")]["hit"]["t"]
            if tof12 > self.config_anal["tof12_cut_high"] or \
               tof12 < self.config_anal["tof12_cut_low"]:
                tof12_cut = True
        else:
            tof12_cut = True

        return (tof_sp_list, {
            "tof_0_sp":space_points.GetTOF0SpacePointArray().size() != 1,
            "tof_1_sp":space_points.GetTOF1SpacePointArray().size() != 1,
            "tof_2_sp":space_points.GetTOF2SpacePointArray().size() != 1,
            "tof01":tof01_cut,
            "tof12":tof12_cut,
          })

    def load_scifi_space_points(self, scifi_event):
        space_points_out = []
        for space_point in scifi_event.spacepoints():
            # if we require triplets and space point is not a triplet, skip it
            if self.config.will_require_triplets and \
                space_point.get_channels().GetEntries() != 3:
                continue
            position = space_point.get_global_position()
            sp_dict = {
                  "x":position.x(),
                  "y":position.y(),
                  "z":position.z(),
            }
            space_points_out.append({
              "hit":Hit.new_from_dict(sp_dict),
              "detector":["tku_sp_", "tkd_sp_"][space_point.get_tracker()]
            })
            space_points_out[-1]["detector"] += str(space_point.get_station())
        return space_points_out

    def load_track_points(self, scifi_event):
        track_points_out = []
        for track in scifi_event.scifitracks():
            detector = ["tku_tp", "tkd_tp"][track.tracker()]
            for track_point in track.scifitrackpoints():
                #print track.tracker(), track_point.station(), track_point.pos().z()
                if track_point.station() != self.config.tk_station or track_point.plane() != self.config.tk_plane:
                    continue
                position = track_point.pos() # maybe global coordinates?
                momentum = track_point.mom() # maybe global coordinates?
                if abs(position.z() - self.det_pos[detector]) > 5.: # detector station is a couple mm thick
                    det_str = " - detector "+str(detector)+" track pos "+str(position.z())+" det pos "+str(self.det_pos[detector])
                    raise RuntimeError("Track point z position not consistent with config"+det_str)
                index = 0
                cov = [[0. for i in range(6)] for j in range(6)]
                index_mapping = [1, 4, 2, 5, 3] # pz, y, x, py, px
                if track_point.covariance().size() == 25:
                    try:
                        self.nan_check(track_point.covariance())
                    except Exception:
                        #print "\nCov nan in run", self.this_run, "spill", self.this_spill, "event", self.this_event
                        #print "    ", [x for x in track_point.covariance()]
                        break
                    for i in range(5):
                        k = index_mapping[i]
                        for j in range(5):
                            l = index_mapping[j]
                            cov[k][l] = track_point.covariance()[index] # cov[4][1]
                            index += 1
                tp_dict = {
                  "x":position.x(),
                  "y":position.y(),
                  "z":position.z(),
                  "px":momentum.x(),
                  "py":momentum.y(),
                  "pz":momentum.z(),
                  "pid":-13,
                  "mass":xboa.common.pdg_pid_to_mass[13],
                }
                try:
                    self.nan_check(tp_dict.values())
                except:
                    #print "\nTP nan in run", self.this_run, "spill", self.this_spill, "event", self.this_event, ":", tp_dict
                    break
                track_points_out.append({
                  "hit":Hit.new_from_dict(tp_dict, "energy"),
                  "detector":detector,
                  "covariance":cov,
                  "pvalue":track.P_value(),
                })
        track_points_out = sorted(track_points_out, key = lambda tp: tp["hit"]["z"])
        return track_points_out

    def load_scifi_event(self, scifi_event):
        points = []
        will_cut_on_scifi_cluster = self.will_cut_on_scifi_clusters(scifi_event)
        will_cut_on_scifi_space_points = self.will_cut_on_scifi_space_points(scifi_event)
        will_cut_on_scifi_tracks = self.will_cut_on_scifi_tracks(scifi_event)
        will_cut_on_scifi_track_points = self.will_cut_on_scifi_track_points(scifi_event)
        will_cut_on_pvalue = self.will_cut_on_pvalue(scifi_event)

        if self.config.will_load_tk_space_points:
            points += self.load_scifi_space_points(scifi_event)
        if self.config.will_load_tk_track_points:
            points += self.load_track_points(scifi_event)
        return [
            points, {
                "scifi_space_points":will_cut_on_scifi_space_points,
                "scifi_space_clusters":will_cut_on_scifi_cluster, 
                "scifi_tracks_us":will_cut_on_scifi_tracks[0],
                "scifi_tracks_ds":will_cut_on_scifi_tracks[1],
                "scifi_track_points_us":will_cut_on_scifi_track_points[0],
                "scifi_track_points_ds":will_cut_on_scifi_track_points[1],
                "pvalue_us":will_cut_on_pvalue[0],
                "pvalue_ds":will_cut_on_pvalue[1],
            }
        ]

    def will_cut_on_scifi_clusters(self, scifi_event):
        """Require exactly one cluster in each plane"""
        n_clusters = [[[0 for k in range(3)] for j in range(5)] for i in range(2)]
        for cluster in scifi_event.clusters():
            tracker = cluster.get_tracker()
            station = cluster.get_station()-1
            plane = cluster.get_plane()
            n_clusters[tracker][station][plane] += 1
        for i, tracker_cluster in enumerate(n_clusters):
            for station_cluster in tracker_cluster:
                if station_cluster != [1, 1, 1]:
                    return True
        return False

    def will_cut_on_scifi_space_points(self, scifi_event):
        """Require exactly one space point in each station"""
        n_space_points = [[0 for j in range(5)] for i in range(2)]
        for space_point in scifi_event.spacepoints():
            tracker = space_point.get_tracker()
            station = space_point.get_station()-1
            if not self.config.will_require_triplets:
                n_space_points[tracker][station] += 1
            # require triplet space point
            elif space_point.get_channels().GetEntries() == 3:
                n_space_points[tracker][station] += 1
        #print n_space_points
        for i, tracker_space_point in enumerate(n_space_points):
            if i in self.config.required_trackers:
                if tracker_space_point != [1, 1, 1, 1, 1]:
                    return True
        return False

    def will_cut_on_scifi_track_points(self, scifi_event):
        """Require minimum number of track points in at least one track in each tracker; require track point in station 1"""
        okay = [True, True]
        for track in scifi_event.scifitracks():
            points = track.scifitrackpoints()
            if len(points) >= self.config.required_number_of_track_points:
                for track_point in points:
                    if track_point.station() == 1:
                        okay[track.tracker()] = False
        return okay

        #if okay == [False, False]:
        #    return False # will NOT cut
        #else:
        #    return True

    def will_cut_on_scifi_tracks(self, scifi_event):
        """Require exactly one track in each tracker"""
        n_tracks = [0, 0]
        for track in scifi_event.scifitracks():
            n_tracks[track.tracker()] += 1
        n_tracks = [i != 1 for i in n_tracks]
        return n_tracks

    def will_cut_on_pvalue(self, scifi_event):
        """Require any track in a tracker to have pvalue greater than threshold"""
        will_cut = [True, True]
        for track in scifi_event.scifitracks():
            will_cut[track.tracker()] = will_cut[track.tracker()] and \
                                track.P_value() < self.config_anal["pvalue_threshold"]
        return will_cut


    def load_reco_event(self, reco_event):
        tof_event = reco_event.GetTOFEvent()
        global_event = reco_event.GetGlobalEvent()
        scifi_event = reco_event.GetSciFiEvent()
        tof_loaded = self.load_tof_event(tof_event)
        if self.config.will_require_tof1 or self.config.will_require_tof2:
            tofs = [ev["detector"] for ev in tof_loaded[0]]
            if "tof1" not in tofs and self.config.will_require_tof1:
                return None
            if "tof2" not in tofs and self.config.will_require_tof2:
                return None
        scifi_loaded = self.load_scifi_event(scifi_event)
        event = {"data":sorted(scifi_loaded[0]+tof_loaded[0], key = lambda tp: tp['hit']['z'])}
        event["particle_number"] = reco_event.GetPartEventNumber()
        
        cuts_chain = itertools.chain(tof_loaded[1].iteritems(), scifi_loaded[1].iteritems())
        cuts = (elem for elem in cuts_chain)
        event["will_cut"] = dict(cuts) # True means cut the thing
        event = self.tm_data(event)
        return event

    def update_cuts(self):
        for event in self.events:
            event["any_cut"] = False
            event["downstream_cut"] = False
            event["data_cut"] = False
            for key, value in event["will_cut"].iteritems():
                if value and self.config.cuts_active[key]:
                    event["any_cut"] = True
                    event["data_cut"] = True
                if value and self.config.downstream_cuts[key]:
                    event["downstream_cut"] = True

    def will_do_p_cut_us(self, event):
        p_bins = self.config_anal["p_bins"]
        p_low = min(min(p_bins))
        p_high = max(max(p_bins))
        return event["tku"] == None or event["tku"]["p"] < p_low or event["tku"]["p"] > p_high

    def will_do_p_cut_ds(self, event):
        return event["tkd"] == None or \
               event["tkd"]["p"] < self.config_anal["p_tot_ds_low"] or \
               event["tkd"]["p"] > self.config_anal["p_tot_ds_high"]

    def nan_check(self, float_list):
        if float_list == None:
            return
        for x in float_list:
            if math.isnan(x) or math.isinf(x) or x != x:
                raise ZeroDivisionError("Nan found in tracker: "+str(float_list))
            #if x > 150.:
            #    raise ValueError("Tracker recon out of range: "+str(float_list))

    def tm_data(self, event):
        tku = None
        tkd = None
        tof0 = None
        tof1 = None
        tof2 = None
        for point in event["data"]:
            point["hit"]["particle_number"] = self.this_run
            point["hit"]["event_number"] = self.this_event
            point["hit"]["spill"]  = self.this_spill
            if point["detector"] == "tku_tp":
                tku = point["hit"]
            elif point["detector"] == "tkd_tp":
                tkd = point["hit"]
            elif point["detector"] == "tof0":
                tof0 = point["hit"]["t"]
            elif point["detector"] == "tof1":
                tof1 = point["hit"]["t"]
            elif point["detector"] == "tof2":
                tof2 = point["hit"]["t"]
        try:
            event["tof12"] = tof2 - tof1
        except TypeError:
            event["tof12"] = None
        try:
            event["tof01"] = tof1 - tof0
        except TypeError:
            event["tof01"] = None
        event["apertures_us"] = [] # list of apertures that the event hit upstream of tku
        event["apertures_ds"] = [] # list of apertures that the event hit downstream of tku
        event["tku"] = tku
        event["tkd"] = tkd
        event["delta_tof01"] = None # loaded during track extrapolation
        event["delta_tof12"] = None # loaded during track extrapolation
        event["will_cut"]["p_tot_us"] = self.will_do_p_cut_us(event)
        event["will_cut"]["p_tot_ds"] = self.will_do_p_cut_ds(event)
        event["will_cut"]["aperture_us"] = False # loaded during track extrapolation
        event["will_cut"]["aperture_ds"] = False # loaded during track extrapolation
        event["will_cut"]["delta_tof01"] = False # loaded during track extrapolation
        event["will_cut"]["delta_tof12"] = False # loaded during track extrapolation
        self.event_id += 1
        return event

    cuts = None


