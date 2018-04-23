import time
import importlib
import sys

import ROOT

import Configuration
import maus_cpp.globals
import maus_cpp.material

import xboa.common

def initialise_maus():
    configuration = Configuration.Configuration().\
                                          getConfigJSON(command_line_args=True)
    maus_cpp.globals.birth(configuration)

MATERIAL_LIST = []
def material_to_colour(material):
    global MATERIAL_LIST
    if material[0:3] == "G4_":
        material = material[3:]
    if material not in MATERIAL_LIST:
        MATERIAL_LIST.append(material)
    if material in ("Galactic"):
        return None
    elif material in ("AIR", "He"):
        return ROOT.kYellow
    if material in ("Fe"): # "kill volumes"
        return 1
    if material in ("MYLAR", "POLYSTYRENE", "NYLON-6-6", "POLYCARBONATE", "POLYVINYL_TOLUENE", "POLYURETHANE", "G10", "TUFNOL"):
        return 8
    if material in ("Zn", "Cu", "W", "Al", "ALUMINUM", "TUNGSTEN", "BRASS", "STEEL", "IRON"):
        return 2
    if material in ("lH2", "MICE_LITHIUM_HYDRIDE", "LITHIUM_HYDRIDE", "RenCast6400"):
        return 4
    print "UNRECOGNISED MATERIAL", material
    return 1

def get_materials(radius, z_start, z_end, z_step):
    x = radius
    material = None
    material_start = []
    n_steps = int((z_end-z_start)/z_step)
    for i in range(n_steps):
        z = z_step*i+z_start
        maus_cpp.material.set_position(x, 0., z)
        material_data = maus_cpp.material.get_material_data()
        new_material = material_data['name']
        if new_material != material:
            material = new_material
            material_start.append({"x":x, "z":z, "material":material})
    return material_start

ROOT_GRAPHS = []
def plot_materials(r_start, r_end, r_step, z_start, z_end, z_step, name):
    global ROOT_GRAPHS
    if name == "" or name == None:
        name = "materials"
    canvas = xboa.common.make_root_canvas(name)
    canvas.SetWindowSize(1900, 1000)
    n_steps = int((r_end-r_start)/r_step)
    title = name+" z_step: "+str(z_step)+" r_step: "+str(r_step)
    hist = ROOT.TH2D(name, title+";z [mm]; x [mm]", 1000, z_start, z_end, 1000, r_start, r_end)
    hist.SetStats(False)
    hist.Draw()
    ROOT_GRAPHS.append(hist)
    for i in range(n_steps):
        r = r_step*i+r_start
        materials = get_materials(r, z_start,z_end, z_step)
        print "At radius", r, "found", len(materials), "materials using", len(ROOT_GRAPHS), "root objects"
        for i, material in enumerate(materials):
            colour = material_to_colour(material["material"])
            if colour == None:
                continue
            z_min = material["z"]
            radius = material["x"]
            if i+1 >= len(materials):
                z_max = z_end+1
            else:
                z_max = materials[i+1]["z"]
            if i == 0:
                z_min -= 1
            graph = ROOT.TGraph(2)
            graph.SetPoint(0, z_min, radius)
            graph.SetPoint(1, z_max, radius)
            graph.SetLineColor(colour)
            graph.SetMarkerColor(colour)
            graph.SetMarkerStyle(6)
            graph.SetLineWidth(2)
            graph.Draw("plsame")
            ROOT_GRAPHS.append(graph)
            if i % 10 == 0:
                canvas.Update()

    canvas.Update()
    for format in "png", "eps", "root":
        canvas.Print("plots/"+name+"."+format)

def get_z_tk(config_mod):
    config = importlib.import_module("config_reco").Config
    z_list = [(det[0], det[2]) for det in config.detectors if "tku" in det[2] or "tkd" in det[2]]
    print "Found tracker detectors at", z_list
    return z_list

def plot_trackers():
    z_tk_list = get_z_tk("scripts/config_reco.py")
    initialise_maus()
    old_time = time.time()
    for z_tk, name in z_tk_list[4:5]:
        plot_materials(-1.5, 1.5, 0.01, z_tk-3., z_tk+3., 0.01, name = name)
    print "Plotting took", time.time() - old_time, "seconds"
    print "Found the following materials", MATERIAL_LIST 

def main():
    plot_trackers()
    raw_input()

if __name__ == "__main__":
    main()
