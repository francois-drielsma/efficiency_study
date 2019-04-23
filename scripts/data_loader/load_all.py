import sys
import glob
import time
import bisect

import ROOT
import libMausCpp
import numpy

import xboa.common
from xboa.hit import Hit

import utilities.utilities
# Fix tracker cuts!
from load_mc import LoadMC
from load_reco import LoadReco

class LoadAll(object):
    """
    Top level data loader object. The aim is to:
    * Load data into a source-independent format, so we don't have to have all
       sorts of hacky functions dependent on what detector we want to study
    * Load data from multiple files "invisibly". We don't want to have to
       worry about file handling etc during the analysis.
    The main method in this class is the "load_spills" function, which will load
    a bunch of data.

    The resultant data is stored in the "events" list. The events list is 
    considered to be ephemeral, i.e. purged at every call to load_spills (to
    avoid big memory footprint). The events list format goes like:
      event = {
        "tku":None or <Hit type> with hit from tku station 1
        "tkd":None or <Hit type> with hit from tkd station 1
        "data":list of all detector hits
        "tof01":None or float
        "tof12":None or float
        "p_tof01":float estimate of momentum through tof01
        "will_cut":{"<cut name>":True (will cut) or False (will not be cut)
        # The following sample flags are found by logical OR of will_cut entries
        # True means they should be cut from the sample (excluded)
        "upstream_cut":<bool> cut from the upstream reconstructed sample
        "data_recorder_cut":<bool> cut from the sample for hybrid mc
        "downstream_cut":<bool> cut from the downstream reconstructed sample
        "extrapolation_cut":<bool> cut from the sample for global extrapolation plots
        "mc_true_us_cut":<bool> cut from the upstream truth sample
        "mc_true_ds_cut":<bool> cut from the downstream truth sample
      }
    
    The detector hits are stored in the "data" list. This has a format like
      data = [
        {
          "detector":<string> name of the detector that made the hit
          "hit":<Hit type>
          <some detector specific data>
        },
        {
          ...
        }
      ]
    The Hit type is a (mainly) C object implemented in xboa.Hit/hitcore

    load_all handles the overall file and spill loops. Parsing of the recon
    event is delegated to reco_loader (of LoadReco type) and parsing of the mc
    event is delegated to mc_loader (of LoadMC type).
    """
    def __init__(self, config, config_anal):
        """
        Initialise empty data
        """
        self.config = config
        self.config_anal = config_anal
        self.maus_version = ""
        self.run_numbers = set([])
        self.file_name_list = []

        self.spill_count = 0
        self.suspect_spill_count = 0
        self.event_count = 0
        self.accepted_count = 0
        self.start_time = time.time()

        self.this_file_name = "a"
        self.this_file_number = -1
        self.this_root_file = None
        self.this_run = 0
        self.this_spill = 0
        self.this_event = 0
        self.this_daq_event = 0
        self.this_tree = None
        self.all_root_files = [None] # f@@@ing root deletes histograms randomly if I let files go out of scope

        self.mc_loader = LoadMC(config, config_anal)
        self.reco_loader = LoadReco(config, config_anal)

        self.events = []

    def load_spills(self, number_of_daq_events):
        """
        Load a number of spills from the files
        - number_of_daq_events: number of daq events to load (daq events, not
                                physics_events)
        If we run out of spills from this file, try the next file. If the file
        won't load, keep on with the next one. Print status every 60 seconds or
        every file, whichever is shorter.

        Call "load_one_spill" subfunction to do the spill parsing.
        """
        load_spills_daq_event = 0 # number of daq events loaded during this call
                                  # to load_spills
        self.load_new_file()
        while load_spills_daq_event < number_of_daq_events and \
              self.this_file_name != "" and \
              self.this_tree != None and \
              self.check_spill_count():
            sys.stdout.flush()
            if self.this_file_name == "" or self.this_tree == None:
                break # ran out of files
            old_t = time.time()
            while self.this_daq_event < self.this_tree.GetEntries() and \
                  load_spills_daq_event < number_of_daq_events:
                new_t = time.time()
                if new_t - old_t > 60.:
                    print "Spill", self.this_daq_event, "Time", round(new_t - self.start_time, 2)
                    old_t = new_t
                    sys.stdout.flush()
                try:
                    self.this_tree.GetEntry(self.this_daq_event)
                except SystemError: # abort the file
                    sys.excepthook(*sys.exc_info())
                    print "Aborting file", self.this_file_name
                    self.this_daq_event = self.this_tree.GetEntries()
                    break
                spill = self.this_data.GetSpill()
                self.load_one_spill(spill)
                load_spills_daq_event += 1
                self.this_daq_event += 1
            if self.this_daq_event >= self.this_tree.GetEntries():
                self.next_file()
                self.load_new_file()
            print "  ...loaded", load_spills_daq_event, "'daq events'", \
                  self.spill_count, "'physics_event' spills, ", \
                  self.event_count,"events and", \
                  self.reco_loader.nan_count, "tracker nans"
            #self.mc_loader.print_virtual_detectors()
            if self.this_tree != None:
                print " at", self.this_daq_event, "/", self.this_tree.GetEntries(), "spills from file", self.this_file_name, self.this_run
            else:
                print
            sys.stdout.flush()
        self.this_root_file.Close()
        self.this_tree = None
        self.update_cuts()
        # return True if there are more events
        if True:
            virt_data = [(numpy.mean(self.mc_loader.virtual_dict[key]),
                          numpy.std(self.mc_loader.virtual_dict[key]), key) 
                          for key in self.mc_loader.virtual_dict.keys()]
            for data in sorted(virt_data):
                print data[2].ljust(40), format(data[0], "10.6g"), format(data[1], "10.6g")
        return self.this_file_name != ""

    def clear_data(self):
        """Clear any ephemeral data"""
        self.events = []

    def load_data(self, min_daq_event, max_daq_event):
        """
        Alternative to load_spills, enabling start at a particular daq event
        """
        self.get_file_list()
        self.root_event = min_daq_event
        self.load_spills(max_daq_event - min_daq_event)

    def get_file_list(self):
        """
        Store the list of files, based on glob of config_anal["reco_files"] and 
        do some pre-loading setup.
        """
        self.file_name_list = []
        for fname in self.config_anal["reco_files"]:
            self.file_name_list += glob.glob(fname)
        self.file_name_list = sorted(self.file_name_list)
        if len(self.file_name_list) == 0:
            raise RuntimeError("No files from "+str(self.config_anal["reco_files"]))
        print "Found", len(self.file_name_list), "files"
        print "    ", self.file_name_list[0:3], "...", self.file_name_list[-3:]
        self.next_file()
        self.this_daq_event = 0
        self.spill_count = 0

    def next_file(self):
        """
        Move on to the next file
        """
        try:
            self.this_file_name = self.file_name_list.pop(0)
            self.this_file_number += 1
            print "Loading ROOT file", self.this_file_name, self.this_file_number
        except IndexError:
            self.this_file_name = ""
            print "No more files to load"
        self.this_tree = None
        self.this_daq_event = 0

    def check_spill_count(self):
        """
        Helper function; check whether we have loaded the number of files 
        specified in config
        """
        return self.spill_count < self.config.number_of_spills or \
               self.config.number_of_spills == None

    def load_new_file(self):
        """
        Open a new file for reading
        """
        while self.this_tree == None and self.this_file_name != "":
            self.all_root_files[0] = self.this_root_file
            self.this_root_file = ROOT.TFile(self.this_file_name, "READ") # pylint: disable = E1101
            self.this_data = ROOT.MAUS.Data() # pylint: disable = E1101
            self.this_tree = self.this_root_file.Get("Spill")
            self.this_run = None
            try:
                self.this_tree.SetBranchAddress("data", self.this_data)
            except AttributeError:
                print "Failed to load 'Spill' tree for file", self.this_file_name, "maybe it isnt a MAUS output file?"
                self.this_tree = None
                self.next_file()
                continue

    def load_one_spill(self, spill):
        """
        Load the contents of one spill. If physics_event, loop over reco_events
        and mc_events; get reco_loader mc_loader to load the respective event
        type. 
        """
        old_this_run = self.this_run
        try:
            self.this_run = max(spill.GetRunNumber(), self.this_file_number) # mc runs all have run number 0
        except ReferenceError:
            print "WARNING: Spill was NULL"
            self.suspect_spill_count += 1
            return
        self.run_numbers.add(self.this_run)
        self.this_spill = spill.GetSpillNumber()
        if old_this_run != None and old_this_run != self.this_run:
            # Nb: Durga figured out this issue was related to DAQ saturating
            # and failing to fill the "run number" int for some spills
            print "WARNING: run number changed from", old_this_run, "to", self.this_run,
            print "in file", self.this_file_name, "daq event", self.this_daq_event,
            print "spill", spill.GetSpillNumber(), "n recon events", spill.GetReconEvents().size(), "<------------WARNING"
            self.suspect_spill_count += 1
        if spill.GetDaqEventType() == "physics_event":
            self.spill_count += 1
            for ev_number, reco_event in enumerate(spill.GetReconEvents()):
                self.this_event = reco_event.GetPartEventNumber()
                event = {"data":[]}
                try:
                    self.reco_loader.load(event, spill, ev_number)
                    if len(event["data"]) == 0: # missing TOF1 - not considered further
                        continue 
                    self.mc_loader.load(event, spill, ev_number)
                except ValueError:
                    print "spill", spill.GetSpillNumber(), "particle_number", reco_event.GetPartEventNumber()
                    sys.excepthook(*sys.exc_info())
                except ZeroDivisionError:
                    pass
                event["run"] = self.this_run
                self.event_count += 1
                event["spill"] = spill.GetSpillNumber()
                self.events.append(event)
                event["event_number"] = ev_number
                for hit in event["data"]:
                    hit["hit"]["event_number"] = ev_number
                    hit["hit"]["spill"] = spill.GetSpillNumber()
                    hit["hit"]["particle_number"] = self.this_run
                event["data"] = sorted(event["data"], key = lambda hit: hit["hit"]["z"])


    def update_cuts(self):
        """For a given event, update the cuts based on MC and reco data"""
        for event in self.events:
            event["upstream_cut"] = False
            event["data_recorder_cut"] = False
            event["downstream_cut"] = False
            event["extrapolation_cut"] = False
            event["mc_true_us_cut"] = False
            event["mc_true_ds_cut"] = False
            for key, value in event["will_cut"].iteritems():
                if value and self.config.upstream_cuts[key]:
                    event["upstream_cut"] = True
                if value and self.config.data_recorder_cuts[key]:
                    event["data_recorder_cut"] = True
                if value and self.config.downstream_cuts[key]:
                    event["downstream_cut"] = True
                if value and self.config.extrapolation_cuts[key]:
                    event["extrapolation_cut"] = True
                if value and self.config.mc_true_us_cuts[key]:
                    event["mc_true_us_cut"] = True
                if value and self.config.mc_true_ds_cuts[key]:
                    event["mc_true_ds_cut"] = True

