import sys
import math
import itertools
import bisect

import xboa.common
from xboa.hit import Hit

import utilities


class LoadReco(object):
    def __init__(self, config, config_anal):
        self.config = config
        self.config_anal = config_anal

        self.det_pos = self.get_det_pos()
        self.time_offsets = {"tof0":self.config.tof0_offset,
                            "tof1":self.config.tof1_offset,
                            "tof2":self.config.tof2_offset}


    def load(self, event, spill, ev_number):
        reco_event = spill.GetReconEvents()[ev_number]
        spill_number = spill.GetSpillNumber()
        self.load_reco_event(event, reco_event, spill_number)


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
            "t":tof_sp.GetTime()+self.time_offsets[station],
        }
        loaded_sp = {
            "hit":Hit.new_from_dict(sp_dict),
            "detector":station,
            "covariance":cov,
            "dt":tof_sp.GetDt(),
        }
        return loaded_sp

    def load_n_slabs(self, tof_event):
        tof_event_slabs = tof_event.GetTOFEventSlabHit()
        tof0_slabs = [0, 0]
        tof1_slabs = [0, 0]
        tof2_slabs = [0, 0]
        for slab in tof_event_slabs.GetTOF0SlabHitArray():
            tof0_slabs[slab.GetPlane()] += 1
        for slab in tof_event_slabs.GetTOF1SlabHitArray():
            tof1_slabs[slab.GetPlane()] += 1
        for slab in tof_event_slabs.GetTOF2SlabHitArray():
            tof2_slabs[slab.GetPlane()] += 1
        return {"tof0":tof0_slabs, "tof1":tof1_slabs, "tof2":tof2_slabs}

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

        tof12_cut = False
        if "tof2" in detectors and "tof1" in detectors:
            tof12 = tof_sp_list[detectors.index("tof2")]["hit"]["t"] - tof_sp_list[detectors.index("tof1")]["hit"]["t"]
            if tof12 > self.config_anal["tof12_cut_high"] or \
               tof12 < self.config_anal["tof12_cut_low"]:
                tof12_cut = True

        return (tof_sp_list, {
            "tof_0_sp":space_points.GetTOF0SpacePointArray().size() != 1,
            "tof_1_sp":space_points.GetTOF1SpacePointArray().size() != 1,
            "tof_2_sp":len([det for det in detectors if det == "tof2"]) != 1,
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
              "detector":["tku_sp_", "tkd_sp_"][space_point.get_tracker()],
              "n_channels":space_point.get_channels().GetEntries(),
              "is_used":space_point.is_used(),
              "channels":[chan.get_channel() for chan in space_point.get_channels()],
              "npe":[chan.get_npe() for chan in space_point.get_channels()],
            })
            space_points_out[-1]["detector"] += str(space_point.get_station())
        return space_points_out

    def load_rejects(self, pr_track):
        reject_keys = ["sz_chi2", "xy_chi2", "sz_r", "n_sp", "tracker", ]
        cut_values = [65., 10., 150., 99, 99]
        #print "Load rejects\n   ",
        #for key in reject_keys:
        #    print key.rjust(8),
        #print
        reject_values = [
            pr_track.get_rejects_sz_chi2(),
            pr_track.get_rejects_xy_chi2(),
            pr_track.get_rejects_sz_r(),
            pr_track.get_rejects_n_sp(),
            pr_track.get_rejects_tracker(),
        ]
        # Complicated bit:
        # we sort by (sz chi2 cut, xy chi2 cut, sz r cut, sz_chi2, xy_chi2, sz_r)
        # so that passed events are always first

        # convert -1 flag (for failed to even try to calc chi2/etc) to "infinity"
        for i, item in enumerate(reject_values):
            item = [x for x in item]
            for j, x in enumerate(item):
                if item[j] < -1e-3:
                    item[j] = 999
            reject_values[i] = tuple(item)

        # make list of "was rejected" flags followed by list of actual chi2
        was_rejected = []
        for i, item in enumerate(reject_values):
            # False < True so (not "was rejected" events) < ("was rejected" events) => not rejected events evaluate to lower
            was_rejected.append(tuple([cut_values[i] < x for x in item]))

        sort_list = was_rejected+reject_values
        sort_list = zip(*sort_list)
        n_sp = 5 #pr_track.get_spacepoints_pointers().size()
        sort_list = [item for item in sort_list if item[8] == n_sp]
        sort_list = sorted(sort_list)
        if len(sort_list) > 10000:
            print reject_keys
            print cut_values
            for i, item in enumerate(sort_list):
                print "REJECTS ", i, item
            print
        for i, item in enumerate(sort_list[0:0]):
            if i > 2:
                break
            print "   ",
            for value in item:
                if value == False or value == True:
                    continue
                print str(round(value, 2)).rjust(8),
            print
        sort_list = zip(*sort_list)
        if len(sort_list) > 0:
            sort_list = dict([(key, sort_list[i+len(reject_keys)]) for i, key in enumerate(reject_keys)])
        else:
            sort_list = dict([(key, []) for i, key in enumerate(reject_keys)])
        return sort_list

    def load_track_points(self, scifi_event):
        track_points_out = []
        for track in scifi_event.scifitracks():
            detector = ["tku_tp", "tkd_tp"][track.tracker()]
            this_track_points = []
            for track_point in track.scifitrackpoints():
                if track_point.station() != self.config.tk_station:
                    if track_point.plane() != 2:
                        continue
                    detector = detector[:3]+"_"+str(track_point.station())
                elif track_point.plane() != self.config.tk_plane:
                    continue
                position = track_point.pos() # maybe global coordinates?
                momentum = track_point.mom() # maybe global coordinates?
                if abs(position.z() - self.det_pos[detector]) > 10.: # detector station is a couple mm thick
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
                        this_track_points.append({
                          "detector":detector+"_nan",
                          "hit":{"z":999},
                        })
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
                    this_track_points.append({
                      "detector":detector+"_nan",
                      "hit":{"z":999},
                    })
                    break
                pr_tracks = [prtk for prtk in scifi_event.helicalprtracks() if prtk.get_tracker() == track.tracker()]
                #print detector, "with", len(pr_tracks), "helical pr tracks"
                this_track_points.append({
                  "hit":Hit.new_from_dict(tp_dict, "energy"),
                  "detector":detector,
                  "covariance":cov,
                  "pvalue":track.P_value(),
                  "chi2":track.chi2(),
                  "ndf":track.ndf(),
                  "max_r2":-99., # maximum radius squared between this track point and the next one in the same tracker
                  "rejects":self.load_rejects(track.pr_track_pointer()),
                })
            this_track_points = sorted(this_track_points)
            this_track_points = self.set_max_r(this_track_points)
            track_points_out += this_track_points
        track_points_out = sorted(track_points_out, key = lambda tp: tp["hit"]["z"])
        return track_points_out

    def get_xyphi(self, a_track_point):
        # x = x0 + r0 cos(\phi) ## BUG sine vs cosine
        # y = y0 + r0 sin(\phi)
        # r^2 = x^2 + y^2 = x0^2+y0^2+r0^2 + 2 r0 [x0 cos(\phi) + y0 sin(\phi)]
        # dr/dphi = 0 => d(r^2)/dphi = 0 when \phi = \phi_{max}
        # d(r^2)/dphi = 2 r0 [- x0 sin(\phi_{max}) + y0 cos(\phi_{max})] = 0
        # => tan(\phi_{max}) = y0/x0
        return 0., 0., 0.

    def set_max_r(self, track_point_list):
        for i, tp in enumerate(track_point_list):
            track_point_list[i]["max_r2"] = tp["hit"]["x"]**2 + tp["hit"]["y"]**2
        return track_point_list
        
    def will_cut_on_scifi_fiducial(self, track_point_list):
        max_r2 = [-111., -111.]
        for tp in track_point_list:
            if "tku" in tp["detector"]:
                max_r2[0] = max(max_r2[0], tp["max_r2"])
            else:
                max_r2[1] = max(max_r2[1], tp["max_r2"])
        fiducial_r2 = self.config_anal["tracker_fiducial_radius"]**2
        return [max_r2[0] > fiducial_r2, max_r2[0] > fiducial_r2]
                
    def load_scifi_event(self, scifi_event):
        will_cut_on_scifi_cluster = self.will_cut_on_scifi_clusters(scifi_event)
        will_cut_on_scifi_space_points = self.will_cut_on_scifi_space_points(scifi_event)
        will_cut_on_scifi_tracks = self.will_cut_on_scifi_tracks(scifi_event)
        will_cut_on_scifi_track_points = self.will_cut_on_scifi_track_points(scifi_event)
        will_cut_on_pvalue = self.will_cut_on_pvalue(scifi_event)
        will_cut_on_chi2 = self.will_cut_on_chi2(scifi_event)

        if self.config.will_load_tk_space_points:
            space_points = self.load_scifi_space_points(scifi_event)
        if self.config.will_load_tk_track_points:
            track_points = self.load_track_points(scifi_event)
        # a bit hacky - we flag the nans and exclude from recon
        nan_cut_upstream = [point for point in track_points if point["detector"] == "tku_tp_nan"]
        nan_cut_downstream = [point for point in track_points if point["detector"] == "tkd_tp_nan"]
        track_points = [point for point in track_points if point["detector"] != "tku_tp_nan" and \
                                                           point["detector"] != "tkd_tp_nan"]
        n_tkd = len([point for point in track_points if "tkd" in point["detector"] ])
        if n_tkd == 0 and will_cut_on_scifi_tracks[1] == False:
            print "NAN?"
            print nan_cut_downstream
            raw_input()
        will_cut_on_scifi_fiducial = self.will_cut_on_scifi_fiducial(track_points)
        points = space_points + track_points
        return [
            points, {
                "scifi_space_points":will_cut_on_scifi_space_points,
                "scifi_space_clusters":will_cut_on_scifi_cluster, 
                "scifi_tracks_us":will_cut_on_scifi_tracks[0],
                "scifi_tracks_ds":will_cut_on_scifi_tracks[1],
                "scifi_fiducial_us":will_cut_on_scifi_fiducial[0], # BUG - not finished yet
                "scifi_fiducial_ds":will_cut_on_scifi_fiducial[1],
                "scifi_track_points_us":will_cut_on_scifi_track_points[0],
                "scifi_track_points_ds":will_cut_on_scifi_track_points[1],
                "scifi_nan_us":len(nan_cut_upstream) != 0,
                "scifi_nan_ds":len(nan_cut_downstream) != 0,
                "pvalue_us":will_cut_on_pvalue[0],
                "pvalue_ds":will_cut_on_pvalue[1],
                "chi2_us":will_cut_on_chi2[0],
                "chi2_ds":will_cut_on_chi2[1],
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

    def get_n_tracks(self, scifi_event):
        n_tracks = [0, 0]
        for track in scifi_event.scifitracks():
            n_tracks[track.tracker()] += 1
        return n_tracks

    def get_n_planes_with_clusters(self, scifi_event):
        """
        Number of planes with at least one cluster
        """
        n_clusters = [[0 for i in range(15)], [0 for i in range(15)]]
        for cluster in scifi_event.clusters():
            tracker = cluster.get_tracker()
            plane = (cluster.get_station()-1)*3+cluster.get_plane()
            n_clusters[tracker][plane] += 1
        n_clusters[0] = len([i for i in n_clusters[0] if i > 0])
        n_clusters[1] = len([i for i in n_clusters[1] if i > 0])
        return n_clusters

    def will_cut_on_scifi_tracks(self, scifi_event):
        """Require exactly one track in each tracker"""
        n_tracks = self.get_n_tracks(scifi_event)
        n_tracks = [i != 1 for i in n_tracks]
        return n_tracks

    def will_cut_on_pvalue(self, scifi_event):
        """Require any track in a tracker to have pvalue greater than threshold"""
        will_cut = [True, True]
        for track in scifi_event.scifitracks():
            will_cut[track.tracker()] = will_cut[track.tracker()] and \
                                track.P_value() < self.config_anal["pvalue_threshold"]
        return will_cut

    def will_cut_on_chi2(self, scifi_event):
        """Cut if all tracks in a tracker to have chi2 greater than threshold"""
        will_cut = [0, 0] #-1 passes; 0 no result; 1 fails
        chi2_cut = self.config_anal["chi2_threshold"]
        for track in scifi_event.scifitracks():
            tracker = track.tracker()
            chi2_df = track.chi2()/track.ndf()
            if chi2_df > chi2_cut: # this track fails cut
                if will_cut[tracker] != -1: # no previous track passed cut; so chuck out
                    will_cut[tracker] = 1
                else: # a previous track passed cut; so passes
                    will_cut[tracker] = -1
            else: # this track passed cut so passes
                will_cut[tracker] = -1
        for i in range(2): # will_cut is true if all tracks failed
            will_cut[i] = will_cut[i] == 1
        return will_cut

    detectors = {
        32:"_tof0",
        33:"_tof0",
        34:"_tof0",
        35:"_ckov",
        36:"_ckov",
        37:"_tof1",
        38:"_tof1",
        39:"_tof1",
        40:"_tku_tp",
        41:"_tku_tp",
        42:"_tku_2",
        43:"_tku_3",
        44:"_tku_4",
        45:"_tku_5",
        46:"_tkd_tp",
        47:"_tkd_tp",
        48:"_tkd_2",
        49:"_tkd_3",
        50:"_tkd_4",
        51:"_tkd_5",
        52:"_tof2",
        53:"_tof2",
        54:"_tof2",
        55:"_kl",
        56:"_emr",
    }

    global_set = set()
    def global_detector_lookup(self, z_pos, detector_type, prefix):
        if detector_type == 0:
            raise KeyError("Failed to recognise global track point with undefined type.")
        elif detector_type != 1: # 1 is virtual
            if detector_type not in self.detectors.keys():
                #print "no id for detector", detector_type
                return "global_real"
            return prefix+self.detectors[detector_type]
        # virtual detectors
        det_list = self.config.virtual_detectors
        detector = bisect.bisect_left(det_list, (z_pos, None, ""))
        if detector == 0:
            det_pos, dummy, det_name = det_list[detector]
        elif detector == len(det_list):
            det_pos, dummy, det_name = det_list[-1]
        elif det_list[detector][0] - z_pos < z_pos - det_list[detector-1][0]:
            det_pos, dummy, det_name = det_list[detector]
        else:
            det_pos, dummy, det_name = det_list[detector-1]
        #if abs(z_pos - det_pos) > 1.:
            #self.virtual_fail_list.append(z_pos)
            #print "Failed to find virtual at", z_pos
        my_det = prefix+"_"+det_name
        #self.global_set.add(my_det)
        return my_det

    def get_global_prefix(self, track):
        tku_pos = utilities.detector_position("tku_tp", self.config)
        tkd_pos = utilities.detector_position("tkd_tp", self.config)
        has_us, has_ds = False, False
        for track_point in track.GetTrackPoints():
            pos = track_point.get_position()
            has_us = has_us or pos.Z() < tku_pos # needs to clear tku by 100 mm...
            has_ds = has_ds or pos.Z() > tkd_pos # needs to clear tku by 100 mm...
        if has_us and has_ds:
            return "global_through"
        elif has_us:
            return "global_us"
        elif has_ds:
            return "global_ds"
        return None # not useful for this analysis

    def get_global_prefix_alt(self, track):
        tku_pos = utilities.detector_position("tku_tp", self.config)
        has_us = False
        has_ds = False
        print "    global prefix alt",
        has_us = track.GetTrackPoints()[0].get_position().Z() < tku_pos - 100.
        has_ds = track.GetTrackPoints().back().get_position().Z() > tku_pos + 100.
        
        print

        has_tof1 = track.GetTrackPoints(37).size() > 0
        has_tku = track.GetTrackPoints(41).size() > 0
        has_tkd = track.GetTrackPoints(47).size() > 0

        has_us = has_us and has_tof1 and has_tku
        has_ds = has_ds and has_tkd
        prefix = None
        if has_us and has_ds:
            prefix = "global_through"
        elif has_us:
            prefix = "global_us"
        elif has_ds:
            prefix = "global_ds"
        print "    ... prefix allocated:", prefix
        return prefix # not useful for this analysis

    def load_global_event(self, global_event, verbose = False):
        # verbose = True
        cuts = {"upstream_aperture_cut":False, "downstream_aperture_cut":False}
        match_cuts = ["global_through_tof0", "global_through_tof1", "global_through_tof2",
                      "global_through_tku_tp", "global_through_tkd_tp",]
        for a_cut in match_cuts: # require global match of given detector type
            cuts[a_cut] = True
        if not verbose:
            sys.stdout = open("/dev/null", "w")
        print "Load global event"
        pid = self.config_anal["pid"]
        track_points_out = []
        mass = xboa.common.pdg_pid_to_mass[abs(pid)]
        for track in global_event.get_tracks():
            # we require a TKU hit on the track
            this_track_points = []
            prefix = self.get_global_prefix(track)
            if prefix == None or prefix == "global_us":
                continue
            print prefix
            for track_point in track.GetTrackPoints():
                pos = track_point.get_position()
                mom = track_point.get_momentum()
                tp_dict = {
                  "x":pos.X(),
                  "y":pos.Y(),
                  "z":pos.Z(),
                  "t":pos.T(),
                  "pid":pid,
                  "mass":mass,
                }
                if mom.Mag2() > 0.1: # TOF does not have any momentum info stored
                    if False:# and abs(mom.Mag2()/mass - mass) > 1.0: # mass-shell check DISABLED
                        print "\nFailed to load track point; failed mass shell check, with", mom.Mag2()**0.5,
                        print "detector", track_point.get_detector()
                        continue # next track point please
                    tp_dict["px"] = mom.X()
                    tp_dict["py"] = mom.Y()
                    tp_dict["pz"] = mom.Z()
                    tp_dict["energy"] = mom.T()
                detector = self.global_detector_lookup(pos.Z(), track_point.get_detector(), prefix)
                if detector == "global_real":
                    continue
                if "tof" in detector:
                    # through-going tracks get time from tof1
                    if "global_through" in detector or "global_us":
                        tp_dict["t"] -= self.config.tof1_offset
                    # ds tracks get time from tof2
                    elif "global_ds" in detector:
                        tp_dict["t"] -= self.config.tof2_offset
                try:
                    self.nan_check(tp_dict.values())
                except:
                    continue
                hit = Hit.new_from_dict(tp_dict, "energy")
                print "       ", hit["x"], hit["y"], hit["z"], ";", hit["t"], "**", \
                                hit["px"], hit["py"], hit["pz"], ";", hit["energy"], "**", \
                                detector
                this_track_points.append({
                  "hit":hit,
                  "detector":detector,
                })
                if detector in self.config.upstream_aperture_cut:
                    if hit["r"] > self.config.upstream_aperture_cut[detector]:
                       cuts["upstream_aperture_cut"] = True
                if detector in self.config.downstream_aperture_cut:
                    if hit["r"] > self.config.downstream_aperture_cut[detector]:
                       cuts["downstream_aperture_cut"] = True
                for ref_detector in match_cuts:
                    if detector == ref_detector:
                        cuts[ref_detector] = False
            print "   ", track.GetTrackPoints().size(), "track points including", track.GetTrackPoints(1).size(), "virtual hits"
            track_points_out += this_track_points
        sys.stdout = sys.__stdout__
        return track_points_out, cuts

    global_passes = []
    def load_reco_event(self, event, reco_event, spill_number):
        #print "Spill", spill_number, "Event", reco_event.GetPartEventNumber()
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
        verbose = False #spill_number < 5
        global_loaded = self.load_global_event(global_event, verbose)
        event["data"] += scifi_loaded[0]+tof_loaded[0]+global_loaded[0]
        event["tof_slabs"] = self.load_n_slabs(tof_event)
        event["scifi_n_tracks"] = self.get_n_tracks(scifi_event)
        event["scifi_n_planes_with_clusters"] = self.get_n_planes_with_clusters(scifi_event)
        event["particle_number"] = reco_event.GetPartEventNumber()
        
        cuts_chain = itertools.chain(tof_loaded[1].iteritems(), scifi_loaded[1].iteritems(), global_loaded[1].iteritems())
        cuts = (elem for elem in cuts_chain)
        event["will_cut"] = dict(cuts) # True means cut the thing
        event = self.tm_data(event)
        return event

    def will_do_p_cut_us(self, event):
        p_bins = self.config_anal["p_bins"]
        p_low = min(min(p_bins))
        p_high = max(max(p_bins))
        if event["tku"] == None:
            return False
        return event["tku"]["p"] < p_low or event["tku"]["p"] > p_high

    def will_do_p_cut_ds(self, event):
        if event["tkd"] == None:
            return False
        return event["tkd"]["p"] < self.config_anal["p_tot_ds_low"] or \
               event["tkd"]["p"] > self.config_anal["p_tot_ds_high"]

    def will_do_delta_tof01_cut(self, event):
        delta_tof01 = event["delta_tof01"]
        if delta_tof01 == None:
            return False
        return delta_tof01 < self.config_anal["delta_tof01_lower"] or \
               delta_tof01 > self.config_anal["delta_tof01_upper"]
        
    def will_do_delta_tof12_cut(self, event):
        delta_tof12 = event["delta_tof12"]
        if delta_tof12 == None:
            return False
        return delta_tof12 < self.config_anal["delta_tof12_lower"] or \
               delta_tof12 > self.config_anal["delta_tof12_upper"]

    def nan_check(self, float_list):
        if float_list == None:
            return
        for x in float_list:
            if math.isnan(x) or math.isinf(x) or x != x:
                raise ZeroDivisionError("Nan found in tracker: "+str(float_list))
            #if x > 150.:load_file
            #    raise ValueError("Tracker recon out of range: "+str(float_list))

    def tm_data(self, event):
        tku = None
        tkd = None
        tof0 = None
        tof1 = None
        tof2 = None
        global_tof0 = None
        global_tof1 = None
        global_tof2 = None
        for point in event["data"]:
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
            elif point["detector"] == "global_through_tof0":
                global_tof0 = point["hit"]["t"]
            elif point["detector"] == "global_through_tof1":
                global_tof1 = point["hit"]["t"]
            elif point["detector"] == "global_through_tof2":
                global_tof2 = point["hit"]["t"]
        try:
            event["tof12"] = tof2 - tof1
        except TypeError:
            event["tof12"] = None
        try:
            event["tof01"] = tof1 - tof0
        except TypeError:
            event["tof01"] = None
        try:
            event["delta_tof12"] = (global_tof2 - global_tof1) - (tof2 - tof1)
        except TypeError:
            event["delta_tof12"] = None
        try:
            event["delta_tof01"] = (global_tof1 - global_tof0) - (tof1 - tof0)
        except TypeError:
            event["delta_tof01"] = None
        event["apertures_us"] = [] # list of apertures that the event hit upstream of tku
        event["apertures_ds"] = [] # list of apertures that the event hit downstream of tku
        event["tku"] = tku
        event["tkd"] = tkd
        if tku == None: # BUG - some events are missing the tracks cut
            event["will_cut"]["scifi_tracks_us"] = True
        if tkd == None:
            event["will_cut"]["scifi_tracks_ds"] = True
        event["will_cut"]["delta_tof01"] = self.will_do_delta_tof01_cut(event)
        event["will_cut"]["delta_tof12"] = self.will_do_delta_tof12_cut(event)
        event["will_cut"]["p_tot_us"] = self.will_do_p_cut_us(event)
        event["will_cut"]["p_tot_ds"] = self.will_do_p_cut_ds(event)
        event["will_cut"]["tof01_selection"] = False
        return event

    def get_det_pos(self):
        det_pos = {}
        for det in self.config.detectors:
            det_pos[det[2]] = det[0]
        return det_pos

