#!/usr/bin/env python

#  This file is part of MAUS: http://micewww.pp.rl.ac.uk/projects/maus
#
#  MAUS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  MAUS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with MAUS.  If not, see <http://www.gnu.org/licenses/>.

"""
Example to demonstrate how to make a field map.

In this example, MAUS will build a geometry specified in the usual way by the
"simulation_geometry_filename" datacard and then plot the magnetic (b) field
parallel to the beam axis (bz) in the region z=-8 m to z=8 m.
"""

import Configuration  # MAUS configuration (datacards) module
import maus_cpp.globals as maus_globals # MAUS C++ calls
import maus_cpp.field as field # MAUS field map calls
import xboa.common  # xboa post-processor library

def initialise_maus():
    configuration = Configuration.Configuration().\
                                          getConfigJSON(command_line_args=True)
    maus_globals.birth(configuration)

def plot_hps(z_list):
    hp_z_positions = {
        "hp_77":14104.,
        "hp_79":14429.,
        #"hp_67":15286.,
        "hp_65":14429.,
        "hp_66":18625.,
        "hp_72":19482.,
        "hp_73":19807.,
        "hp_80":19482.,
    }
    hp_btot = {
        "hp_77":3.0751005,
        "hp_79":3.0408025,
        #"hp_67":0.,
        "hp_65":3.0430975,
        "hp_66":3.100803,
        "hp_72":3.011501,
        "hp_73":3.0571025,
        "hp_80":3.012287,
    }
    graph_list = []
    hp_color_list = [1, 2, 4, 6, 8]
    i = 0
    for hp in sorted(hp_z_positions.keys()):
        z_pos = hp_z_positions[hp]
        btot = hp_btot[hp]
        name = "Hall probe "+hp[-2:]
        if z_pos > z_list[0] and z_pos < z_list[-1]:
            hist, graph = xboa.common.make_root_graph(name, [z_pos], "", [btot], "")
            graph.SetMarkerStyle(20+len(graph_list))
            hp_color = hp_color_list[len(graph_list)%5]
            graph.SetLineColor(10)
            graph.SetMarkerColor(hp_color)
            graph.Draw("p")
            graph_list.append(graph)
    return graph_list

def plot_tracker_stations(z_list, btot_list):
    tracker_z_positions = [
        (13968.0, "tku_5"),
        (14319.0, "tku_4"),
        (14618.0, "tku_3"),
        (14867.0, "tku_2"),
        (15068.0, "tku_tp"),
        (18756.0, "tkd_tp"),    
        (18955.0, "tkd_2"),
        (19206.0, "tkd_3"),
        (19505.0, "tkd_4"),
        (19855.0, "tkd_5"),
    ]
    graph_list = []
    for z_pos, station in tracker_z_positions:
        if z_pos > z_list[0] and z_pos < z_list[-1]:
            hist, graph = xboa.common.make_root_graph("Tracker stations", [z_pos, z_pos], "", [-1e6, 1e6], "")
            graph.SetLineColor(1)
            graph.SetLineStyle(2)
            graph.Draw("l")
            graph_list.append(graph)
    col = 10
    i = 1
    print "z_u".rjust(col), "z_d".rjust(col), "<b>".rjust(col), "dphi".rjust(col), "dphi [3T]".rjust(col), "dphi/dl".rjust(col)
    for z_pos, station in sorted(tracker_z_positions):
        if i == len(z_list):
            break
        print str(z_list[i]).rjust(col),
        bdl = 0.
        dl = 0.
        while i < len(z_list) and z_list[i] < z_pos:
            bdl += btot_list[i]*(z_list[i]-z_list[i-1])
            dl += (z_list[i]-z_list[i-1])
            i += 1
        if abs(dl) < 1e-9:
            print  "and", z_list[i-1], "is undef"
        else:
            mean_b = bdl/dl
            dphi = bdl/130.*1e-3*300. #kT
            dphi_3 = 3.0*dl/130.*1e-3*300. #kT
            print str(z_list[i-1]).rjust(col), str(round(mean_b, 4)).rjust(col), str(round(dphi, 4)).rjust(col), str(round(dphi_3, 4)).rjust(col), str(round(dphi/dl, 6)).rjust(col)
    return graph_list


def plot_z_range(z_list, name):
    canvas = None
    graph_list = []
    for r_pos, line_color in [(140., 6), (160., 4), (180., 8), (0., 1)]:
        #z_list = [float(z_pos) for z_pos in range(19000, 20001, 10)]
        btot_list = []
        for z_pos in z_list:
            (bx_field, by_field, bz_field, ex_field, ey_field, ez_field) = \
                                     field.get_field_value(r_pos, 0., z_pos, 0.)
            btot = (bx_field**2+by_field**2+bz_field**2)**0.5
            btot_list.append(btot*1e3)  # btot in T
            #print 'z:', z_pos, ' ** b:', bx_field, by_field, bz_field, \
            #                       'e:', ex_field, ey_field, ez_field

        [ymin, ymax] = xboa.common.min_max(btot_list+[3.1])
        [xmin, xmax] = xboa.common.min_max(z_list)
        xmax += (xmax-xmin)*0.3
        # now make a ROOT graph of bz against z
        hist, graph = xboa.common.make_root_graph("x="+str(r_pos)+" mm", z_list, "z [m]",
                                                   btot_list, "B_{tot} [T]",
                                                  xmin=xmin, xmax=xmax,
                                                  ymin=ymin, ymax=ymax)
        graph.SetLineColor(line_color)
        if canvas == None:
            canvas = xboa.common.make_root_canvas("bz vs z")
            hist.Draw()
        graph_list.append(graph)
        graph.Draw('l')
        canvas.Update()
    graph_list += plot_hps(z_list)
    plot_tracker_stations(z_list, btot_list)
    legend = xboa.common.make_root_legend(canvas, graph_list)
    print legend.GetX1NDC(), legend.GetX2NDC()
    legend.SetX1NDC(0.7)
    legend.SetX2NDC(0.9)
    legend.SetBorderSize(1)
    legend.Draw()
    canvas.Print(name) #'plots/bfield_vs_z_ssd.png'

def main():
    """
    Make a plot of z, bz
    """
    # set up datacards (geometry specified e.g. on command line using --simulation_geometry_filename)
    initialise_maus()
    # make plots
    plot_z_range(range(13900, 15101, 1), "bfield_vs_z_ssu.png")
    plot_z_range(range(18600, 20001, 1), "bfield_vs_z_ssd.png")
    plot_z_range(range(13000, 21001, 1), "bfield_vs_z.png")
    # Clean up
    maus_globals.death()
    # Finished

if __name__ == "__main__":
    main()
    print "Finished - press <CR> to end"
    raw_input()
    
