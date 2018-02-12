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

def main():
    """
    Generates some data and then attempts to load it and make a simple histogram
    """
    # first off, we try to generate some data based on some default data file
    # let's generate some data by running the reconstruction...
    print "Generating some data"
    my_file_name = "/home/cr67/work/reco/MAUS-v2.9.1/08645/08645_recon_global_test.root"
    n_spills = None
    target_station = 2
    # now load the ROOT file
    print "Loading ROOT file", my_file_name
    root_file = ROOT.TFile(my_file_name, "READ") # pylint: disable = E1101

    # and set up the data tree to be filled by ROOT IO
    print "Setting up data tree"
    data = ROOT.MAUS.Data() # pylint: disable = E1101
    tree = root_file.Get("Spill")
    tree.SetBranchAddress("data", data)

    print "Getting some data"
    target_dict = {3:[0, 4, 27]}
    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        spill = data.GetSpill()
        if spill.GetDaqEventType() != "physics_event" or \
           spill.GetSpillNumber() not in target_dict.keys():
            if len(target_dict) > 0:
                continue
        part_numbers = []
        if len(target_dict) > 0:
            part_numbers = target_dict[spill.GetSpillNumber()]
        is_dodgy = False
        for reco_event in spill.GetReconEvents():
            global_event = reco_event.GetGlobalEvent()
            for track in global_event.get_tracks():
                for tp in track.GetTrackPoints():
                    if tp.get_detector() == 45 and abs(tp.get_position().T()) > 1000:
                        is_dodgy = True
            if is_dodgy and len(target_dict) > 0:
                print "DODGY", spill.GetSpillNumber(), reco_event.GetPartEventNumber()
                is_dodgy = False
            if reco_event.GetPartEventNumber() not in part_numbers:
                continue
            #if not is_dodgy:
            #    continue
            print "\n\nFound global event for spill", spill.GetSpillNumber(), "and particle event", reco_event.GetPartEventNumber()
            through_chains = global_event.GetThroughPrimaryChains()
            orphan_chains = global_event.GetPrimaryChainOrphans()
            for chain in through_chains:
                print "  Through Chain:"
                read_chain(chain)
            for chain in orphan_chains:
                print "  Orphan Chain:"
                read_chain(chain)

def read_chain(chain):
              for track in chain.GetMatchedTracks():
                print "    Track with", track.GetTrackPoints().size(), "track points"
                for tp in track.GetTrackPoints():
                    pos = tp.get_position()
                    mom = tp.get_momentum()
                    print_data = [
                        ("x", pos.X()),
                        ("y", pos.Y()),
                        ("z", pos.Z()),
                        ("pz", mom.Z()),
                        ("energy", mom.T()),
                        ("t", pos.T()),
                    ]
                    print "   tp",
                    for key, value in print_data:
                        print key+":", str(round(value, 3)).ljust(12),
                    print "det:", str(tp.get_detector()).ljust(4)
                    continue
                    pos = tp.get_space_point().get_position()
                    #mom = tp.get_space_point().get_momentum()
                    print_data = [
                        ("x", pos.X()),
                        ("y", pos.Y()),
                        ("z", pos.Z()),
                        ("pz", 0),
                        ("energy", 0),
                        ("t", pos.T()),
                    ]
                    print "   sp",
                    for key, value in print_data:
                        print key+":", str(round(value, 3)).ljust(12),
                    print "det:", str(tp.get_detector()).ljust(4)
              if chain.GetUSDaughter():
                  print "  US Daughter:"
                  read_chain(chain.GetUSDaughter())
              if chain.GetDSDaughter():
                  print "  DS Daughter:"
                  read_chain(chain.GetDSDaughter())

if __name__ == "__main__":
    main()
    print "Press <CR> to finish"
    #raw_input()

