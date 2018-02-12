import sys
import copy
import numpy
import json

import xboa.common
import maus_cpp.field
import lattice
try:
    import run_sim
except ImportError:
    print "cant do import from csv"

"""
Try to fit B_0 - B(m_1 m_2) = j_{e1} e_1 + j2 f2 + j3 f3
"""


class TrimCoils(object):
    def __init__(self, a_lattice, z_range):
        self.lattice = a_lattice
        self.z_range = z_range

    @classmethod
    def new_from_coils(cls, momentum, coils):
        coils["MatchCoil1_DS"] = 0.
        a_lattice = lattice.Lattice(momentum, 0.)
        a_lattice.set_magnet_scale_factors(coils)
        z_range = range(-3000, -1800, 20)
        return TrimCoils(a_lattice, z_range)

    @classmethod
    def new_from_csv(cls, index, filename):
        print "Trimming coils", index, filename
        settings = run_sim.get_run_settings_from_csv(index, filename)
        settings["MatchCoil1_DS"] = 0.
        if settings == None:
            print "Failed to load csv filename", filename, "index", index
            return None
        a_lattice = lattice.Lattice(settings["momentum"], 0.)
        a_lattice.set_magnet_scale_factors(settings)
        z_range = range(-3000, -1800, 20)
        return TrimCoils(a_lattice, z_range)

    def trim_coils(self, coils_to_trim, target_uniform_field):
        
        coil_currents = self.lattice.coil_currents # this is fixed
        trimmed = copy.deepcopy(coil_currents) # this we fiddle with
        print "Target field", target_uniform_field, "kT"
        print "Initial coil currents\n", json.dumps(trimmed, indent=2)
        for coil in coils_to_trim:
            if coil not in trimmed:
                raise KeyError("Coil "+str(coil)+" does not exist")
            trimmed[coil] = 0.
        self.lattice.set_magnet_scale_factors(trimmed)
        # target field - field from match coils
        y_list = [target_uniform_field - maus_cpp.field.get_field_value(0, 0, z, 0)[2] for z in self.z_range]
        print "Y_LIST", y_list
        for coil in trimmed:
              trimmed[coil] = 0.
        f_list = {}
        for coil in coils_to_trim:
              trimmed[coil] = 1.
              self.lattice.set_magnet_scale_factors(trimmed)
              f_list[coil] = [maus_cpp.field.get_field_value(0, 0, z, 0)[2] for z in self.z_range]
              trimmed[coil] = 0.
              print "F_LIST", coil, f_list[coil]

        y = numpy.transpose(numpy.array([y_list]))
        x = numpy.transpose(numpy.array([f_list[coil] for coil in coils_to_trim]))
        fitted_currents = numpy.linalg.lstsq(x, y)
        trimmed = copy.deepcopy(coil_currents)
        for i, coil in enumerate(coils_to_trim):
            trimmed[coil] = fitted_currents[0][i][0]
        print "Fitted coil currents\n", json.dumps(trimmed, indent=2)
        self.lattice.set_magnet_scale_factors(trimmed)
        return trimmed

def test_lls():
    y_list = [[3.*z**2 + 4.*z for z in range(10)]]        
    x_list = []
    x_list.append([z**2 for z in range(10)])
    x_list.append([z for z in range(10)])

    y = numpy.transpose(numpy.array(y_list))
    x = numpy.transpose(numpy.array(x_list))

    j = numpy.linalg.lstsq(x, y)
    print "Return", j[0].shape
    print j[0]
    print j[1]

def plot_each_coil(lattice, field_canvas, z_list):
    coil_currents = lattice.coil_currents
    trimmed = copy.deepcopy(coil_currents)

    lattice.plot_field(canvas = field_canvas, line_color = 1, z_list = z_list)

    for coil in trimmed:
        trimmed[coil] = 0.
    for coil in trimmed:
        trimmed[coil] = coil_currents[coil]
        lattice.set_magnet_scale_factors(trimmed)
        lattice.plot_field(canvas = field_canvas, line_color = 8, z_list = z_list)
        trimmed[coil] = 0.

    lattice.set_magnet_scale_factors(coil_currents)
    return field_canvas

def test_field_fit():
    lattice = Lattice(200., 0.)
    lattice.coil_currents["MatchCoil2_US"] = 172.94
    lattice.coil_currents["MatchCoil1_US"] = 100.
    lattice.coil_currents["MatchCoil2_DS"] = 0.
    lattice.coil_currents["MatchCoil1_DS"] = 0.
    lattice.coil_currents["FocusCoil_US"] = 40.77
    lattice.coil_currents["FocusCoil_DS"] = 40.77

    # full plot before optimisation
    full_z_list = range(-3500, +3500, 10)
    full_canvas = lattice.plot_field(line_color = 4, z_list = full_z_list)

    # upstream fit
    us_z_list = sorted(range(-2990, -1890, 10))
    trimmer = TrimCoils(lattice, us_z_list)
    trimmer.trim_coils(["EndCoil2_US", "CenterCoil_US", "EndCoil1_US"], 3e-3)
    field_canvas = lattice.plot_field(line_color = 1, z_list = trimmer.z_range)

    # downstream fit
    ds_z_list = sorted(range(1890, 2990, 10))
    trimmer = TrimCoils(lattice, ds_z_list)
    trimmer.trim_coils(["EndCoil2_DS", "CenterCoil_DS", "EndCoil1_DS"], 3e-3)
    field_canvas = lattice.plot_field(line_color = 1, z_list = trimmer.z_range)

    # full plot after optimisation including coils
    plot_each_coil(lattice, full_canvas, full_z_list)
    raw_input()

def main():
    file_name = sys.argv[1].split('/')[-1][:-4]+"_"+sys.argv[2]
    print "Writing trimmed setup to:"
    print "run_plans/trim_coils/"+file_name+"_trims.csv" # runner.py uses this line
    sys.stdout.flush()
    try:
        trimmer = TrimCoils.new_from_csv(int(sys.argv[2]), sys.argv[1])
        fout = open("run_plans/trim_coils/"+file_name+"_trims.csv", "w") # empty file
    except Exception:
        fout = open("run_plans/trim_coils/"+file_name+"_trims.csv", "w") # empty file
        return
    target_field = maus_cpp.field.get_field_value(0., 0., -2450., 0.)[2] # fix bz at centre of CC
    output = trimmer.trim_coils(["EndCoil2_US", "CenterCoil_US", "EndCoil1_US"], target_field)

    field_canvas = xboa.common.make_root_canvas("field")
    #field_canvas = plot_each_coil(trimmer.lattice, field_canvas, trimmer.z_range)
    field_canvas.Print("plots/trim_coils_"+file_name+".png")
    field_canvas.Print("plots/trim_coils_"+file_name+".pdf")

    print output
    print >> fout, " ", trimmer.lattice.p, "0 0 0 0 0",
    for coil in run_sim.CSV_KEYS[6:-1]:
        print >> fout, output[coil],
    print >> fout, 0


if __name__ == "__main__":
    main()


