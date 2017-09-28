import Configuration
import maus_cpp.globals
import maus_cpp.material
import xboa.common
import ROOT
import time

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
    if material in ("Galactic", "AIR", "He"):
        return None
    if material in ("Fe"): # "kill volumes"
        return 1
    if material in ("MYLAR", "POLYSTYRENE", "RenCast6400", "NYLON-6-6", "POLYCARBONATE", "POLYVINYL_TOLUENE", "POLYURETHANE", "G10", "TUFNOL"):
        return 8
    if material in ("Zn", "Cu", "W", "Al", "ALUMINUM", "TUNGSTEN", "BRASS", "STEEL"):
        return 2
    if material in ("lH2", "MICE_LITHIUM_HYDRIDE", "LITHIUM_HYDRIDE"):
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
def plot_materials(r_start, r_end, r_step, z_start, z_end, z_step):
    global ROOT_GRAPHS
    canvas = xboa.common.make_root_canvas("materials")
    n_steps = int((r_end-r_start)/r_step)
    hist = ROOT.TH2D("materials", ";z [mm]; r [mm]", 1000, z_start, z_end, 1000, r_start, r_end)
    hist.SetStats(False)
    hist.Draw()
    ROOT_GRAPHS.append(hist)
    for i in range(n_steps):
        r = r_step*i+r_start
        materials = get_materials(r, z_start,z_end, z_step)
        print "At radius", r, "found", len(materials), "materials"
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
        canvas.Print("plots/materials."+format)

def main():
    initialise_maus()
    old_time = time.time()
    plot_materials(0, 200.5, 1., 12000., 21000., 0.1)
    print "Plotting took", time.time() - old_time, "seconds"
    print "Found the following materials", MATERIAL_LIST

if __name__ == "__main__":
    main()
    raw_input()
