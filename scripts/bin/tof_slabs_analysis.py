#!/usr/bin/env python


"""
Example to load a ROOT file and make a histogram showing the beam profile at
TOF1
"""

import os
import subprocess

# basic PyROOT definitions
import ROOT 
import xboa.common
import utilities

# definitions of MAUS data structure for PyROOT
import libMausCpp #pylint: disable = W0611

def get_slab_delta_t(slab_hits):
    delta_t, raw_delta_t = "no slabs", "no slabs"
    slab_hit_list_0 = [hit for hit in slab_hits if hit.GetPlane() == 0]
    slab_hit_list_1 = [hit for hit in slab_hits if hit.GetPlane() == 1]

    # loop over all combinations of slab hits; return the smallest possible
    # value of delta_t or None if no delta_t could be found
    for hit_0 in slab_hit_list_0:
        for hit_1 in slab_hit_list_1:
            this_delta_t = hit_0.GetTime() - hit_1.GetTime()
            this_raw_delta_t = hit_0.GetRawTime() - hit_1.GetRawTime()
            if abs(this_delta_t) < 1e-12: # calibration failed
                if delta_t == "no slabs":
                    delta_t = "no calibration"
            elif type(delta_t) == type("") or abs(this_delta_t) < abs(delta_t):
                delta_t = this_delta_t
            if type(raw_delta_t) == type("") or abs(this_raw_delta_t) < abs(raw_delta_t):
                raw_delta_t = this_raw_delta_t

    return delta_t, raw_delta_t


def main():
    """
    Generates some data and then attempts to load it and make a simple histogram
    """
    # first off, we try to generate some data based on some default data file
    # let's generate some data by running the reconstruction...
    print "Generating some data"
    #my_file_name = "/home/cr67/MAUS/work/reco/mc/rogers/data_0.62_0.30_v2.8.0/0/maus_output.root"
    #my_file_name = "/home/cr67/MAUS/work/reco/other/08448/08448_recon.root"
    my_file_name = "/home/cr67/work/reco/MAUS-v3.0.0/09805/09805_recon.root"
    n_spills = None
    #target_pid = 211
    target_station = 2

    x0, x1 = -2., 2.
    tof0_counters = [0, 0, 0]
    tof1_counters = [0, 0, 0]
    tof2_counters = [0, 0, 0]
    dt_counter = 0
    dt_canvas = xboa.common.make_root_canvas("TOF dt")
    dt_canvas.SetLogy()
    dt_hist_0 = xboa.common.make_root_histogram("TOF0 dt", [0.], "TOF dt [ns]", 100, xmin=x0, xmax=x1)
    dt_hist_0.SetLineColor(2)
    dt_hist_1 = xboa.common.make_root_histogram("TOF1 dt", [0.], "TOF dt [ns]", 100, xmin=x0, xmax=x1)
    dt_hist_1.SetLineColor(4)
    dt_hist_2 = xboa.common.make_root_histogram("TOF2 dt", [0.], "TOF dt [ns]", 100, xmin=x0, xmax=x1)
    dt_hist_2.SetLineColor(8)

    dt_hist_0_slab = xboa.common.make_root_histogram("TOF0 slabs dt", [0.], "TOF dt", 100, xmin=x0, xmax=x1)
    dt_hist_0_slab.SetLineColor(2)
    dt_hist_0_slab.SetLineStyle(2)
    dt_hist_1_slab = xboa.common.make_root_histogram("TOF1 slabs dt", [0.], "TOF dt", 100, xmin=x0, xmax=x1)
    dt_hist_1_slab.SetLineColor(4)
    dt_hist_1_slab.SetLineStyle(2)
    dt_hist_2_slab = xboa.common.make_root_histogram("TOF2 slabs dt", [0.], "TOF dt", 100, xmin=x0, xmax=x1)
    dt_hist_2_slab.SetLineColor(8)
    dt_hist_2_slab.SetLineStyle(2)

    dt_raw_canvas = xboa.common.make_root_canvas("TOF raw dt")
    dt_raw_canvas.SetLogy()

    dt_hist_0_slab_raw = xboa.common.make_root_histogram("TOF0 slabs raw dt", [0.], "TOF dt", 100, xmin=x0*2, xmax=x1*2)
    dt_hist_0_slab_raw.SetLineColor(2)
    dt_hist_0_slab_raw.SetLineStyle(2)
    dt_hist_1_slab_raw = xboa.common.make_root_histogram("TOF1 slabs raw dt", [0.], "TOF dt", 100, xmin=x0*2, xmax=x1*2)
    dt_hist_1_slab_raw.SetLineColor(4)
    dt_hist_1_slab_raw.SetLineStyle(2)
    dt_hist_2_slab_raw = xboa.common.make_root_histogram("TOF2 slabs raw dt", [0.], "TOF dt", 100, xmin=x0*2, xmax=x1*2)
    dt_hist_2_slab_raw.SetLineColor(8)
    dt_hist_2_slab_raw.SetLineStyle(2)

    dt_hist_0_slab_raw_missed = xboa.common.make_root_histogram("TOF0 slabs raw dt missed", [0.], "TOF dt", 100, xmin=x0*2, xmax=x1*2)
    dt_hist_0_slab_raw_missed.SetLineColor(2)
    dt_hist_1_slab_raw_missed = xboa.common.make_root_histogram("TOF1 slabs raw dt missed", [0.], "TOF dt", 100, xmin=x0*2, xmax=x1*2)
    dt_hist_1_slab_raw_missed.SetLineColor(4)
    dt_hist_2_slab_raw_missed = xboa.common.make_root_histogram("TOF2 slabs raw dt missed", [0.], "TOF dt", 100, xmin=x0*2, xmax=x1*2)
    dt_hist_2_slab_raw_missed.SetLineColor(8)



    # now load the ROOT file
    print "Loading ROOT file", my_file_name
    root_file = ROOT.TFile(my_file_name, "READ") # pylint: disable = E1101

    # and set up the data tree to be filled by ROOT IO
    print "Setting up data tree"
    data = ROOT.MAUS.Data() # pylint: disable = E1101
    tree = root_file.Get("Spill")
    tree.SetBranchAddress("data", data)

    print "Getting some data"
    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        spill = data.GetSpill()
        if spill.GetDaqEventType() == "physics_event":
            for reco_event in spill.GetReconEvents():
                tof_event = reco_event.GetTOFEvent()
                if not tof_event:
                    continue
                if len(tof_event.GetTOFEventSpacePoint().GetTOF1SpacePointArray()) == 0:
                    continue

                tof0_slab_dt, tof0_slab_raw_dt = get_slab_delta_t(tof_event.GetTOFEventSlabHit().GetTOF0SlabHitArray())
                if tof0_slab_dt != "no slabs":
                    dt_hist_0_slab_raw.Fill(tof0_slab_raw_dt)
                    if tof0_slab_dt == "no calibration":
                        tof0_counters[0] += 1
                        dt_hist_0_slab_raw_missed.Fill(tof0_slab_raw_dt)
                    elif abs(tof0_slab_dt) > 0.5:
                        tof0_counters[1] += 1
                        dt_hist_0_slab.Fill(tof0_slab_dt)
                    else:
                        tof0_counters[2] += 1
                        dt_hist_0_slab.Fill(tof0_slab_dt)

                tof1_slab_dt, tof1_slab_raw_dt = get_slab_delta_t(tof_event.GetTOFEventSlabHit().GetTOF1SlabHitArray())
                if tof1_slab_dt != "no slabs":
                    dt_hist_1_slab_raw.Fill(tof1_slab_raw_dt)
                    if tof1_slab_dt == "no calibration":
                        tof1_counters[0] += 1
                        dt_hist_1_slab_raw_missed.Fill(tof1_slab_raw_dt)
                    elif abs(tof1_slab_dt) > 0.5:
                        tof1_counters[1] += 1
                        dt_hist_1_slab.Fill(tof1_slab_dt)
                    else:
                        tof1_counters[2] += 1
                        dt_hist_1_slab.Fill(tof1_slab_dt)

                tof2_slab_dt, tof2_slab_raw_dt = get_slab_delta_t(tof_event.GetTOFEventSlabHit().GetTOF2SlabHitArray())
                if tof2_slab_dt != "no slabs":
                    dt_hist_2_slab_raw.Fill(tof2_slab_raw_dt)
                    if tof2_slab_dt == "no calibration":
                        dt_hist_2_slab_raw_missed.Fill(tof2_slab_raw_dt)
                        tof2_counters[0] += 1
                    elif abs(tof2_slab_dt) > 0.5:
                        tof2_counters[1] += 1
                        dt_hist_2_slab.Fill(tof2_slab_dt)
                    else:
                        tof2_counters[2] += 1
                        dt_hist_2_slab.Fill(tof2_slab_dt)
                
                for tof_sp in tof_event.GetTOFEventSpacePoint().GetTOF0SpacePointArray()[:1]:
                    dt_hist_0.Fill(tof_sp.GetDt())
                for tof_sp in tof_event.GetTOFEventSpacePoint().GetTOF1SpacePointArray()[:1]:
                    dt_hist_1.Fill(tof_sp.GetDt())
                for tof_sp in tof_event.GetTOFEventSpacePoint().GetTOF2SpacePointArray()[:1]:
                    dt_hist_2.Fill(tof_sp.GetDt())

            if i % 100 == 0:
                dt_canvas.cd()
                dt_counter = 0
                #dt_hist_1_slab.Draw()
                #dt_hist_0_slab.Draw("SAME")
                #dt_hist_2_slab.Draw("SAME")
                dt_hist_1.Draw("SAME")
                dt_hist_0.Draw("SAME")
                dt_hist_2.Draw("SAME")
                dt_canvas.Update()

                dt_raw_canvas.cd()
                dt_hist_1_slab_raw.Draw()
                dt_hist_0_slab_raw.Draw("SAME")
                dt_hist_2_slab_raw.Draw("SAME")
                dt_hist_1_slab_raw_missed.Draw("SAME")
                dt_hist_0_slab_raw_missed.Draw("SAME")
                dt_hist_2_slab_raw_missed.Draw("SAME")
                dt_raw_canvas.Update()

                print "DAQ Event", i, "of", tree.GetEntries()
                print "TOF0 failed calibration", round(100.*tof0_counters[0]/sum(tof0_counters), 2),  \
                      "% out of range", round(100.*tof0_counters[1]/sum(tof0_counters), 2), \
                      "% Okay", round(100.*tof0_counters[2]/sum(tof0_counters), 2), "%"
                print "TOF1 failed calibration", round(100.*tof1_counters[0]/sum(tof1_counters), 2),  \
                      "% out of range", round(100.*tof1_counters[1]/sum(tof1_counters), 2), \
                      "% Okay", round(100.*tof1_counters[2]/sum(tof1_counters), 2), "%"
                print "TOF2 failed calibration", round(100.*tof2_counters[0]/sum(tof2_counters), 2),  \
                      "% out of range", round(100.*tof2_counters[1]/sum(tof2_counters), 2), \
                      "% Okay", round(100.*tof2_counters[2]/sum(tof2_counters), 2), "%"
                if n_spills != None and i > n_spills:
                    break

    dt_canvas.cd()
    dt_fit_1 = utilities.fit_peak(dt_hist_1)
    dt_fit_0 = utilities.fit_peak(dt_hist_0)
    dt_fit_2 = utilities.fit_peak(dt_hist_2)
    get_text_box(dt_fit_0, dt_fit_1, dt_fit_2, dt_hist_0, dt_hist_1, dt_hist_2)
    dt_canvas.Update()
    for format in "png", "pdf", "root":
        dt_canvas.Print("tof_validation."+format)
        dt_raw_canvas.Print("tof_validation_raw."+format)
    root_file.Close()

TEXT_BOXES = []
def get_text_box(fit_0, fit_1, fit_2, hist_0, hist_1, hist_2):
    text_box = ROOT.TPaveText(0.6, 0.4, 0.9, 0.9, "NDC")
    text_box.SetFillColor(0)
    text_box.SetBorderSize(0)
    text_box.SetTextSize(0.04)
    text_box.SetTextAlign(12)
    text_box.SetTextSize(0.03)

    text_box.AddText("TOF0 (Red)")
    #text_box.AddText("  Number:    "+str(hist_0.GetEntries()))
    text_box.AddText("  Mean:        "+str(round(fit_0.GetParameter(1), 2))+" ns")
    text_box.AddText("  Std:           "+str(round(fit_0.GetParameter(2), 2))+" ns")
    text_box.AddText("TOF1 (Blue)")
    #text_box.AddText("  Number:    "+str(hist_1.GetEntries()))
    text_box.AddText("  Mean:        "+str(round(fit_1.GetParameter(1), 2))+" ns")
    text_box.AddText("  Std:           "+str(round(fit_1.GetParameter(2), 2))+" ns")
    text_box.AddText("TOF2 (Green)")
    #text_box.AddText("  Number:    "+str(hist_2.GetEntries()))
    text_box.AddText("  Mean:        "+str(round(fit_2.GetParameter(1), 2))+" ns")
    text_box.AddText("  Std:           "+str(round(fit_2.GetParameter(2), 2))+" ns")
    text_box.SetBorderSize(1)
    text_box.Draw()
    TEXT_BOXES.append(text_box)
    return text_box


if __name__ == "__main__":
    main()
    print "Press <CR> to finish"
    #raw_input()

