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

import os
import shutil
import ROOT

import Configuration  # MAUS configuration (datacards) module
import maus_cpp.globals as maus_globals # MAUS C++ calls
import maus_cpp.field as field # MAUS field map calls
import xboa.common  # xboa post-processor library

Z_OFFSET = 0. #(15068.0+18756.0)/2.

def initialise_maus():
    configuration = Configuration.Configuration().\
                                          getConfigJSON(command_line_args=True)
    maus_globals.birth(configuration)

def plot_hps(z_list, different_markers):
    ssu_offset = 6 # mm, compare MICE Note 518 with MICE Note 501
    ssd_offset = 10 # mm, compare MICE Note 518 with MICE Note 501
    hp_z_positions = {
        "hp_77":14104.+ssu_offset,
        "hp_79":14429.+ssu_offset,
        "hp_65":14429.+ssu_offset,
        #"hp_67":15286.+ssu_offset,
        #"hp_66":18625.+ssd_offset,
        "hp_72":19482.+ssd_offset,
        "hp_73":19807.+ssd_offset,
        "hp_80":19482.+ssd_offset,
    }
    hp_btot = {
        "hp_77":3.0717,
        "hp_79":3.0361,
        #"hp_67":3.2280,
        "hp_65":3.0379,
        #"hp_66":2.3315,
        "hp_72":2.1196,
        "hp_73":-1e9, #3.0571025,
        "hp_80":-1e9, #3.012287,
    }
    graph_list = []
    hp_color_list = [1, 2, 4, 6, 8]
    i = 0
    for hp in sorted(hp_z_positions.keys()):
        z_pos = hp_z_positions[hp]
        btot = hp_btot[hp]
        name = "Hall probe "+hp[-2:]
        if z_pos > z_list[0] and z_pos < z_list[-1]:
            print "Plotting", name, "at", z_pos-Z_OFFSET, "with", btot, "T"
            hist, graph = xboa.common.make_root_graph(name, [z_pos-Z_OFFSET], "", [btot], "")
            graph.SetMarkerStyle(20+len(graph_list))
            hp_color = hp_color_list[len(graph_list)%5]
            graph.SetLineColor(10)
            graph.SetMarkerColor(hp_color)
            graph.Draw("p")
            if different_markers:
                graph_list.append(graph)
    if different_markers:
        return graph_list
    else:
        graph.SetName("Hall Probes")
        return [graph]

def plot_tracker_stations(z_list, btot_list):
    tku_trp = 15068
    tkd_trp = 18836.8+8.
    tracker_z_positions = [ # Ref email from Melissa 12th Jan 2017
        (tku_trp-1100., "tku_5"),
        (tku_trp-750., "tku_4"),
        (tku_trp-450., "tku_3"),
        (tku_trp-200., "tku_2"),
        (tku_trp, "tku_tp"),
        ((tku_trp+tkd_trp)/2., "absorber"),
        (tkd_trp, "tkd_tp"),    
        (tkd_trp+200., "tkd_2"),
        (tkd_trp+450., "tkd_3"),
        (tkd_trp+750., "tkd_4"),
        (tkd_trp+1100., "tkd_5"),
    ]
    graph_list = []
    for z_pos, station in tracker_z_positions:
        if z_pos > z_list[0] and z_pos < z_list[-1]:
            hist, graph = xboa.common.make_root_graph("Tracker stations", [z_pos-Z_OFFSET, z_pos-Z_OFFSET], "", [-1e6, 1e6], "")
            graph.SetLineColor(4)
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

TEXT_BOXES = []
def text_box(graph_list):
    legend = ROOT.TLegend(0.65, 0.7, 0.89, 0.89)
    for graph in graph_list:
        legend.AddEntry(graph, graph.GetName(), "lp")
    legend.SetBorderSize(0)
    legend.Draw()
    TEXT_BOXES.append(legend)
    text_box = ROOT.TPaveText(0.65, 0.6, 0.89, 0.7, "NDC")
    text_box.SetFillColor(0)
    text_box.SetBorderSize(0)
    text_box.SetTextSize(0.03)
    text_box.SetTextAlign(12)
    text_box.AddText("MICE Preliminary")
    text_box.Draw()
    TEXT_BOXES.append(text_box)
    text_box = ROOT.TPaveText(0.65, 0.5, 0.89, 0.6, "NDC")
    text_box.SetFillColor(0)
    text_box.SetBorderSize(0)
    text_box.SetTextSize(0.024)
    text_box.SetTextAlign(12)
    text_box.AddText("Setting 2017/02-7")
    text_box.Draw()
    TEXT_BOXES.append(text_box)

def plot_z_range(z_list, b_min_max, name):
    canvas = None
    graph_list = []
    for r_pos, line_color in [(160., 8), (0., 1)]:
        #z_list = [float(z_pos) for z_pos in range(19000, 20001, 10)]
        btot_list = []
        for z_pos in z_list:
            (bx_field, by_field, bz_field, ex_field, ey_field, ez_field) = \
                                     field.get_field_value(r_pos, 0., z_pos, 0.)
            btot = (bx_field**2+by_field**2+bz_field**2)**0.5
            btot_list.append(btot*1e3)  # btot in T
            #print 'z:', z_pos, ' ** b:', bx_field, by_field, bz_field, \
            #                       'e:', ex_field, ey_field, ez_field
        gz_list = [z - Z_OFFSET for z in z_list]
        [ymin, ymax] = [b_min_max[0], b_min_max[1]] # xboa.common.min_max(btot_list+[3.1])
        [xmin, xmax] = xboa.common.min_max(gz_list)
        xmax += (xmax-xmin)*0.3
        # now make a ROOT graph of bz against z
        hist, graph = xboa.common.make_root_graph("x="+str(r_pos)+" mm", gz_list, "z [mm]",
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
    graph_list += plot_hps(z_list, False)
    plot_tracker_stations(z_list, btot_list)
    text_box(graph_list)
    for format in ["root", "eps", "png", "pdf"]:
        canvas.Print(name+"."+format)

def main():
    """
    Make a plot of z, bz
    """
    # set up datacards (geometry specified e.g. on command line using --simulation_geometry_filename)
    my_dir = "field/"
    if os.path.exists(my_dir):
        shutil.rmtree(my_dir)
    os.makedirs(my_dir)
    initialise_maus()
        
    # make plots
    plot_z_range(range(13700, 15401, 1), [2.95, 3.3], my_dir+"bfield_vs_z_ssu")
    plot_z_range(range(18600, 20201, 1), [2.0, 2.4], my_dir+"bfield_vs_z_ssd")
    plot_z_range(range(13000, 21001, 1), [None, None], my_dir+"bfield_vs_z")
    # Clean up
    maus_globals.death()
    # Finished

if __name__ == "__main__":
    main()
    #raw_input()
    print "Finished - press <CR> to end"
    
