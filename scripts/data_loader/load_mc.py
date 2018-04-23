import bisect

import xboa.common
from xboa.hit import Hit


def get_tk_station(hit):
    chan = hit.GetChannelId()
    tk = chan.GetTrackerNumber()+1
    plane = chan.GetPlaneNumber()+1
    st = chan.GetStationNumber()
    station = tk*100+st*10+plane
    return station


class LoadMC(object):
    def __init__(self, config, config_anal):
        self.config = config
        self.config_anal = config_anal

    def load(self, event, spill, ev_number):
        mc_event_vector = spill.GetMCEvents()
        do_mc = self.config_anal["do_mc"] or self.config_anal["amplitude_mc"] 
        if not do_mc or len(mc_event_vector) < ev_number:
            for cut in self.virtual_cut_list:
                event["will_cut"][cut] = False
            return
        # mc load functions return a list of hits
        mc_event = mc_event_vector[ev_number]
        primary = mc_event.GetPrimary()
        primary_cuts, primary_loaded = self.load_primary(primary)
        track_vector = mc_event.GetTracks()
        track_cuts, track_loaded = self.load_tracks(track_vector)
        virtual_vector = mc_event.GetVirtualHits()
        virtual_cuts, virtual_loaded = self.load_virtuals(virtual_vector)
        tof_hit_vector = mc_event.GetTOFHits()
        tof_hit_lambda = lambda hit: hit.GetChannelId().GetStation()
        tof_cuts, tof_hit_loaded = self.load_mc_hits(tof_hit_vector, "mc_tof", tof_hit_lambda)
        tk_hit_vector = mc_event.GetSciFiHits()
        tk_cuts, tk_hit_loaded = self.load_mc_hits(tk_hit_vector, "mc_tk", get_tk_station)
        temp_event = sorted(primary_loaded+track_loaded+virtual_loaded+tof_hit_loaded+tk_hit_loaded, key = lambda tp: tp['hit']['z'])
        event["data"] += temp_event
        for item in primary_cuts, track_cuts, virtual_cuts, tof_cuts:
            for key, value in item.iteritems():
                event["will_cut"][key] = value

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
        return {}, loaded_track_vector

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
        return {}, [loaded_primary]

    global_set = set()
    def virtual_detector_lookup(self, z_pos):
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
        my_det = "mc_"+det_name
        #self.global_set.add(my_det)
        return my_det

    def virtual_cuts(self, detector, hit):
        # will_cut if not a muon
        if detector == "mc_virtual_tku_tp":
            self.virtual_cuts_tmp[0] = ("mc_muon_us", abs(hit["pid"]) != 13)
        elif detector == "mc_virtual_tkd_tp":
            self.virtual_cuts_tmp[3] = ("mc_muon_ds", abs(hit["pid"]) != 13)
         # will cut if > 150 mm
        if "mc_virtual_tku" in detector:
            self.virtual_cuts_tmp[2] = ("mc_scifi_fiducial_us", hit["r"] > 150. or self.virtual_cuts_tmp[1][1])
        elif "mc_virtual_tku" in detector:
            self.virtual_cuts_tmp[5] = ("mc_scifi_fiducial_ds", hit["r"] > 150. or self.virtual_cuts_tmp[1][1])
        mc_stations = self.config.mc_plots["mc_stations"]
        for station in mc_stations["tku"]:
            if station == detector:
                self.station_count["tku"].add(station)
        self.virtual_cuts_tmp[1] = ("mc_stations_us", len(self.station_count["tku"]) != len(mc_stations["tku"]) )
        for station in mc_stations["tkd"]:
            if station == detector:
                self.station_count["tkd"].add(station)
        self.virtual_cuts_tmp[4] = ("mc_stations_ds", len(self.station_count["tkd"]) != len(mc_stations["tkd"]) )

    def load_virtuals(self, virtual_vector):
        loaded_virtual_vector = [None]*len(virtual_vector)
        self.virtual_cuts_tmp = [(cut, False) for cut in self.virtual_cut_list]
        self.station_count = {"tku":set(), "tkd":set()}
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
            detector = self.virtual_detector_lookup(hit["z"])
            loaded_virtual = {
                "hit":hit,
                "detector":detector,
            }
            loaded_virtual_vector[i] = loaded_virtual
            self.virtual_cuts(detector, hit)
        return dict(self.virtual_cuts_tmp), loaded_virtual_vector

    def load_mc_hits(self, mc_vector, detector, station_lambda):
        loaded_mc_vector = [None]*len(mc_vector)
        for i, mc_hit in enumerate(mc_vector):
            hit = xboa.hit.Hit()
            hit["x"] = mc_hit.GetPosition().x()
            hit["y"] = mc_hit.GetPosition().y()
            hit["z"] = mc_hit.GetPosition().z()
            hit["px"] = mc_hit.GetMomentum().x()
            hit["py"] = mc_hit.GetMomentum().y()
            hit["pz"] = mc_hit.GetMomentum().z()
            hit["pid"] = mc_hit.GetParticleId()
            hit["station"] = station_lambda(mc_hit)
            hit["particle_number"] = mc_hit.GetTrackId()
            hit["t"] = mc_hit.GetTime()
            hit["energy"] = mc_hit.GetEnergy()
            hit["e_dep"] = mc_hit.GetEnergyDeposited()
            try:
                hit["mass"] = xboa.common.pdg_pid_to_mass[abs(mc_hit.GetParticleId())]
            except KeyError:
                hit["mass"] = (hit["energy"]**2-hit["p"]**2)**0.5
            loaded_mc = {
                "hit":hit,
                "detector":detector+"_"+str(hit["station"]),
            }
            loaded_mc_vector[i] = loaded_mc
        #print "Loaded", [mc["detector"] for mc in loaded_mc_vector]
        return {}, loaded_mc_vector

    virtual_cut_list = ["mc_muon_us", "mc_stations_us", "mc_scifi_fiducial_us",
                        "mc_muon_ds", "mc_stations_ds", "mc_scifi_fiducial_ds"]
