import time
import importlib
import sys
import math

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
    if material in ("MYLAR", "POLYSTYRENE", "NYLON-6-6", "POLYCARBONATE", "POLYVINYL_TOLUENE", "POLYURETHANE", "G10"):
        return 8
    if material in ("Zn", "Cu", "W", "Al", "ALUMINUM", "TUNGSTEN", "BRASS", "STEEL", "IRON", "TAM1000"):
        return 2
    if material in ("lH2", "MICE_LITHIUM_HYDRIDE", "LITHIUM_HYDRIDE", "RenCast6400", "TUFNOL"):
        return 4
    print "UNRECOGNISED MATERIAL", material
    return 1

PRINT_VOLUMES = True
THETA = math.pi/4.
def get_materials(radius, z_start, z_end, z_step):
    global PRINT_VOLUMES, THETA
    x = radius*math.cos(THETA)
    y = radius*math.sin(THETA)
    material = None
    volume = None
    material_start = []
    n_steps = int((z_end-z_start)/z_step)
    for i in range(n_steps):
        z = z_step*i+z_start
        maus_cpp.material.set_position(x, y, z)
        material_data = maus_cpp.material.get_material_data()
        new_material = material_data['name']
        new_volume = material_data['volume']
        if new_material != material:
            material = new_material
            material_start.append({"x":x, "y":y, "z":z, "material":material})
        if PRINT_VOLUMES and new_volume != volume:
            volume = new_volume
            print (volume+" "+material+" "+str(round(z, 1))+"  "),
    if PRINT_VOLUMES:
        print
    return material_start

ROOT_GRAPHS = []
def plot_materials(r_start, r_end, r_step, z_start, z_end, z_step, name):
    global ROOT_GRAPHS, PRINT_VOLUMES
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
    for i in range(n_steps): # radial steps
        r = r_step*i+r_start
        materials = get_materials(r, z_start,z_end, z_step)
        print "At radius", r, "found", len(materials), "materials using", len(ROOT_GRAPHS), "root objects"
        for i, material in enumerate(materials):
            colour = material_to_colour(material["material"])
            if colour == None:
                continue
            z_min = material["z"]
            if abs(r) > 1e-9:
                radius = (material["x"]**2+material["y"]**2)**0.5*r/abs(r)
            else:
                radius = 0.
            if i+1 >= len(materials):
                z_max = z_end+1
            else:
                z_max = materials[i+1]["z"]
            if i == 0:
                z_min -= 1
            print material["material"], round(z_min), round(z_max), colour, "   ",
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
    print

    canvas.Update()
    for format in "png", "eps", "root":
        canvas.Print("plots/"+name+"."+format)

def get_z_tk():
    config = importlib.import_module("config.config_reco").Config
    z_list = [(det[0], det[2]) for det in config.detectors if "tku_tp" in det[2] or "tkd_tp" in det[2]]
    print "Found tracker detectors at", z_list
    return z_list

def plot_trackers():
    initialise_maus()
    old_time = time.time()
    plot_materials(-200.1, 200.1, 1., 13000, 20000., 1., name = "geometry")
    tk_list = get_z_tk()
    print tk_list
    for name, z_tk in tk_list:
        plot_materials(-2.1, 2.1, 0.01, z_tk-2., z_tk+2, 0.1, name = name)
    print "Plotting took", time.time() - old_time, "seconds"
    print "Found the following materials", MATERIAL_LIST 

def main():
    plot_trackers()
    raw_input()

if __name__ == "__main__":
    main()
