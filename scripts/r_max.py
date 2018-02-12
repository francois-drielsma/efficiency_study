import json
import math

import xboa.hit

import Configuration
import maus_cpp.globals
import maus_cpp.global_error_tracking


def get_r_max(track_point_list, bz):
    #x1, y1, phi1 = self.get_xyphi(track_point_list[0])
    # x0, y0, phi0 = x1, y1, phi1
    # x1, y1, phi1 = self.get_xyphi(tp)
    # max_phi = -math.atan2(y0, x0) # BUG CHECK signs
    # subtle logic; need to be sure that we get the pi -> -pi boundary right/etc
    
    for i, tp_0 in enumerate(track_point_list[:-1]):
        pt = tp_0["hit"]["pt"]
        phi_0 = math.atan2(tp_0["hit"]["px"], tp_0["hit"]["py"])
        bz = 0.
        q = 1.
        # coordinates of helix centre
        x0 = q*bz*tp_0["hit"]["x"] - pt*math.sin(phi_0)
        y0 = q*bz*tp_0["hit"]["y"] - pt*math.cos(phi_0)
        # azimuthal angle of maximum extent
        phi_max = math.atan2(x0, y0)
        tp_1 = track_point_list[i+1]
        phi_1 = math.atan2(tp_1["hit"]["px"], tp_1["hit"]["py"])
        while phi_0 < 0.: # phi_0 is always in domain 0 < 2 pi
            phi_0 += 2*math.pi
        while phi_1 < phi_0: # phi_1 is always > phi_0 (BUG do muons go in +ve phi direction?)
            phi_1 += 2*math.pi
        while phi_max < phi_0:
            phi_max += 2*math.pi
        max_r2 = tp_0["hit"]["x"]**2+tp_0["hit"]["y"]**2
        if phi_max < phi_1:
            max_r2 = x0**2 + y0**2 + r0**2 + 2*x0*r0*math.sin(phi_max) + 2*y0*x0*math.cos(phi_max)
        tp_0["max_r2"] = max_r2

def hit_from_psv(psv):
    hit_dict = {
        "t":psv[0],
        "x":psv[1],
        "y":psv[2],
        "z":psv[3],
        "energy":psv[4],
        "px":psv[5],
        "py":psv[6],
        "pz":psv[7]
    }
    hit = xboa.hit.Hit.new_from_dict(hit_dict, "")
    detector_hit = {
        "hit":hit,
        "max_r2":0.,
    }
    return detector_hit

def test_get_r_max(bz, x0, y0, px0, py0, pz0, z_list):
    energy = (px0**2+py0**2+pz0**2+105.658**2)**0.5
    var = [0., x0, y0, 0., energy, px0, py0, pz0]
    ellipse = [[0. for i in range(6)] for j in range(6)]
    tracking = maus_cpp.global_error_tracking.GlobalErrorTracking()
    max_r2 = 0.
    track_list = [hit_from_psv(var)]
    for z in z_list:
        for i in range(z):
            var, ellipse = tracking.propagate_errors(var, ellipse, i)
            r2 = var[1]**2+var[2]**2
            max_r2 = max(r2, max_r2)
            #print i, var[1], var[2], max_r2**0.5
        track_list.append(hit_from_psv(var))
    print [(tp["hit"]["x"], tp["hit"]["y"], tp["hit"]["z"]) for tp in track_list]
    get_r_max(track_list, 3e-3)
    print "calculated", [tp["max_r2"]**0.5 for tp in track_list]
    print "solution  ", max_r2**0.5

def setup():
    config = Configuration.Configuration().getConfigJSON()
    config = json.loads(config)
    config["simulation_geometry_filename"] = "scripts/dedx/test_track.dat"
    config = json.dumps(config)
    maus_cpp.globals.birth(config)

def main():
    setup()
    test_get_r_max(3e-3, 50., 100., 20., 30., 140., [5000])

if __name__ == "__main__":
    main()
  