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
    elif material in ("AIR"):
        return ROOT.kYellow
    elif material in ("He"):
        return ROOT.kYellow+1
    if material in ("Fe"): # "kill volumes"
        return 1 # black
    if material in ("MYLAR", "POLYSTYRENE", "NYLON-6-6", "POLYCARBONATE", "POLYVINYL_TOLUENE", "POLYURETHANE", "G10"):
        return 8 # dark green
    if material in ("Al", "ALUMINUM"):
        return ROOT.kOrange # orange/mustard
    if material in ("Zn", "Cu", "W", "TUNGSTEN", "BRASS", "STEEL", "IRON", "TAM1000"):
        return 2 # red
    if material in ("lH2", "MICE_LITHIUM_HYDRIDE", "LITHIUM_HYDRIDE", "TrackerGlue", "TUFNOL"):
        return 4 # blue
    print "UNRECOGNISED MATERIAL", material
    return 1

PRINT_VOLUMES = True
THETA = math.pi/4.
def get_material_start_recursive(x, y, z0, z1, mat0, mat1, z_tolerance):
    if abs(z0-z1) < z_tolerance/2.:
        #print "\n\n"
        return z0
    new_z = (z0+z1)/2.
    maus_cpp.material.set_position(x, y, new_z)
    new_material = maus_cpp.material.get_material_data()['name']
    #print "    0:", z0, mat0, "** 1:", z1, mat1, "** new:", new_z, new_material,
    if new_material == mat0:
        z0 = new_z
    else:
        z1 = new_z
    #print " ... ", z0, z1
    return get_material_start_recursive(x, y, z0, z1, mat0, mat1, z_tolerance)

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
            if material == None:
                z_boundary = z
            else:
                z_boundary = get_material_start_recursive(x, y, z-z_step, z,
                                                       material, new_material,
                                                       1e-6)
            material = new_material
            material_start.append({"x":x, "y":y, "z":z_boundary,
                                   "material":material,
                                   "volume":new_volume})
        if PRINT_VOLUMES and new_volume != volume:
            volume = new_volume
            print (volume+" "+material+" "+str(round(z_boundary, 4))+"  "),
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
    print "\n\n"+"*"*20, "Plotting materials", name, "*"*20
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
            #print material["material"], round(z_min), round(z_max), colour, "   ",
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

def plot_thickness(volume_list, r_start, r_end, r_step, z_start, z_end, z_step, name):
    print "\n\n"+"*"*20, "Plotting thickness", name, "*"*20
    z_list = []
    r_list = []
    n_steps = int((r_end-r_start)/r_step)
    for i in range(n_steps): # radial steps
        r = r_step*i+r_start
        materials = get_materials(r, z_start,z_end, z_step)
        z_min, z_max = None, None
        for i, material in enumerate(materials):
            if material['volume'] in volume_list:
                z_min = material['z']
                break
        if z_min == None:
            continue
        i += 1
        if i >= len(material):
            continue
        z_max = materials[i]['z']
        r_list.append(r)
        z_list.append(z_max-z_min)
        print r_list[-1], z_list[-1]
    canvas = xboa.common.make_root_canvas(name)
    hist, graph = xboa.common.make_root_graph(name,
                                              z_list, "thickness [mm]",
                                              r_list, "radius [mm]",
                                              xmin=0.)
    hist.SetTitle(name)
    canvas.Draw()
    hist.Draw()
    graph.Draw("L SAME")
    for format in "png", "eps", "root":
        canvas.Print("plots/"+name+"."+format)
    return canvas, hist, graph

def analytical_tracker_window_thickness():
    y_list = [float(i) for i in range(106)]
    thickness_list = []
    for y in y_list:
        z_outer = (173.26**2-y**2)**0.5-2.76
        z_inner = (170.33**2-y**2)**0.5
        thickness_list.append(z_outer-z_inner)
    hist, graph = xboa.common.make_root_graph("", thickness_list, "", y_list, "")
    graph.SetLineColor(ROOT.kBlue)
    return graph

def get_z_tk():
    config = importlib.import_module("config.config_reco").Config
    z_list = [(det[0], det[2]) for det in config.detectors if "tku_tp" in det[2] or "tkd_tp" in det[2]]
    print "Found tracker detectors at", z_list
    return z_list

def get_z_diffuser():
    config = importlib.import_module("config.config_reco").Config
    z_list = [z for z, dummy, name in config.virtual_detectors if "diffuser" in name]
    z_list = sorted(z_list)
    print "Found diffusers at", z_list
    return z_list



def plot_trackers():
    old_time = time.time()
    #plot_materials(-250.1, 250.1, 1., 12000, 23000., 1., name = "materials")
    try:
        tk_list = get_z_tk()
        print tk_list
        #for z_tk, name in tk_list:
        #    plot_materials(-2.1, 2.1, 0.01, z_tk-2., z_tk+2, 0.01, name = name)
        diff_z = get_z_diffuser()
        plot_materials(-500.1, 500.1, 1., min(diff_z)-800., max(diff_z)+500., 0.1, "diffuser")
    except ImportError:
        print "Failed to do tracker zoom"
        sys.excepthook(*sys.exc_info())
    print "Plotting took", time.time() - old_time, "seconds"
    print "Found the following materials", MATERIAL_LIST 

def plot_tracker_windows():
    #plot_materials(0., 200.1, 1., 15180., 15240., 0.1, name = "tku_window")
    volumes = ["HeWindowOuterMaterial_phys", "HeWindowInnerMaterial_phys"]
    canvas, hist, graph = plot_thickness(volumes, 0., 200.1, 1., 15180., 15240., 0.1, name = "tku_window_thickness")
    graph = analytical_tracker_window_thickness()
    graph.Draw("L SAME")
    canvas.Update()
    for format in "png", "eps", "root":
        canvas.Print("plots/tku_window_thickness."+format)
    #plot_materials(0., 200.1, 1., 18670., 18750., 0.1, name = "tkd_window")
    canvas, hist, graph = plot_thickness(volumes, 0., 200.1, 1., 18670., 18750., 0.1, name = "tkd_window_thickness")
    graph = analytical_tracker_window_thickness()
    graph.Draw("L SAME")
    canvas.Update()
    for format in "png", "eps", "root":
        canvas.Print("plots/tkd_window_thickness."+format)

def plot_apertures():
    plot_materials( 0., 300.1, 1., 16500., 17500., 0.1, name = "absorber")
    plot_materials( 0., 300.1, 1., 17500., 18500., 0.1, name = "ssd_aperture")

def main():
    initialise_maus()
    plot_apertures()
    #raw_input()

if __name__ == "__main__":
    main()
