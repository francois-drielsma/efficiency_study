import copy

def mc_file_names(datasets):
    file_list = ["/home/cr67/work/reco/mc/"+datasets+"/*.root"]
    return file_list

def reco_file_names(run_number_list, maus):
    file_list = []
    for run in run_number_list:
        run = str(run).rjust(5, '0')
        a_file = "/home/cr67/MAUS/work/reco/"+maus+"/"+run+"/"+run+"_recon.root"
        file_list.append(a_file)
    return file_list

def get_analysis(job_name, name, tof01_max, data_dir):
    plot_dir = data_dir+"plots_"+name+"/"
    plot_dir = plot_dir.replace(" ", "_")
    return {
            "plot_dir":plot_dir, # makedirs and then put plots in this directory. Removes any old plots there!!!
            "tof12_cut_low":33., # TOF12 cut lower bound
            "tof12_cut_high":45., # TOF12 cut upper bound
            "delta_tof01_lower":-5., # Delta TOF01 cut lower bound 
            "delta_tof01_upper":5., # Delta TOF01 cut upper bound 
            "delta_tof12_lower":-5., # Delta TOF01 cut lower bound 
            "delta_tof12_upper":5., # Delta TOF01 cut upper bound 
            "tof01_cut_low":28., # TOF01 cut lower bound
            "tof01_cut_high":tof01_max, # TOF01 cut upper bound
            "p_bins":[[135, 145]], # set of momentum bins; for now really it is just a lower and upper bound
            "p_tot_ds_low":80., # downstream momentum cut lower bound
            "p_tot_ds_high":200., # downstream momentum cut upper bound
            "reco_files":mc_file_names(job_name), # list of strings to be handed to glob
            "name":"2016-04 1.2 "+name, # appears on plots
            "color":4, # not used
            "pid":-13, # assume pid of tracks following TOF cut
            "pvalue_threshold":0.02, # minimum allowed pvalue for pvalue cut
            "chi2_threshold":5.0, # maximum allowed chi2/dof for chi2 cut
            "amplitude_source":None,
            "field_uncertainty":0.02,
            "do_magnet_alignment":False,
            "do_globals":False,
            "do_amplitude":True,
            "do_extrapolation":False, #name == "3-140+M3-Test2",
            "do_mc":True, #True,
            "do_plots":True,
            "csv_output_detectors":["tof1", "diffuser_us", "diffuser_mid", "diffuser_ds"], # write data at listed detector locations
            "csv_output_filename":None, #"8590_mc_extrapolated_tracks.csv", # write a summary output of data in flat text format to listed filename; set to None to do nothing
            "extrapolation_source":"tku_tp",
        }


class Config(object):
    geometry = "Test.dat" #"geometry_08681/ParentGeometryFile.dat" 
    # location to which data and plots will be dumped following analysis
    data_dir = "output/2017-02_1_mc/"
    info_file = "geometry_08681/Maus_Information.gdml"
    will_require_tof1 = True # require at least one TOF1 Space point to even load the data
    will_require_tof2 = False # require at least one TOF2 Space point to even load the data
    tk_station = 1 # load data from a particular tracker station
    tk_plane = 0
    # prerequisite for space point cut
    will_require_triplets = False #True # require triplet space points
    upstream_cuts = { # Set to true to make data_plotter and amplitude_analysis use these cuts; False to ignore the cut
          "downstream_cut":None,
          "any_cut":None,
          "scifi_space_clusters":False,
          "scifi_space_points":False,
          "scifi_tracks_us":True,
          "scifi_nan_us":True,
          "scifi_track_points_us":False,
          "aperture_us":False,
          "pvalue_us":False,
          "chi2_us":True,
          "aperture_ds":False,
          "scifi_tracks_ds":False,
          "scifi_track_points_ds":False,
          "scifi_nan_ds":False,
          "pvalue_ds":False,
          "chi2_ds":False,
          "tof01":True,
          "tof12":False,
          "delta_tof01":False, #extrapolatedtof01 compared to recon tof01
          "delta_tof12":False, #extrapolatedtof12 compared to recon tof12
          "p_tot_us":True,
          "p_tot_ds":False,
          "tof_0_sp":True,
          "tof_1_sp":True,
          "tof_2_sp":False,
    }
    downstream_cuts = copy.deepcopy(upstream_cuts)
    downstream_cuts["p_tot_ds"] = True
    downstream_cuts["pvalue_ds"] = False
    downstream_cuts["chi2_ds"] = True
    downstream_cuts["tof12"] = True # if TOF12 is out of range chuck it (but ignore "no TOF2" events)
    downstream_cuts["scifi_nan_ds"] = True
    downstream_cuts["scifi_tracks_ds"] = True
    extrapolation_cuts = upstream_cuts

    analyses = []
    #analyses.append(get_analysis("rogers/data_8685_franchini_scale-d1=1.02_d2=1.02_ds=1.0_Br12.87_W8.4_hi-stats/*", "10-140 MC", 30, "MAUS-v2.8.2", data_dir))
    #analyses.append(get_analysis("rogers/data_8699_franchini_scale-d1=1.02_d2=1.02_ds=1.0_Br2.0_W1.6_hi-stats//*", "6-140 MC", 31, "MAUS-v2.8.2", data_dir))
    #analyses.append(get_analysis("franchini/", "3-140 MC Franchini", 32, "MAUS-v2.8.2", data_dir))
    analyses.append(get_analysis("000067/00???", "10-140 MC", 30, data_dir)) #
    analyses.append(get_analysis("000068/00???", "6-140 MC", 31, data_dir)) #
    analyses.append(get_analysis("000071/00???", "3-140 MC", 32, data_dir)) #
    #analyses.append(get_analysis("franchini/", "3-140 MC TestAbstraction2", 32, "MAUS-v2.8.2", data_dir)) #
    amplitude_bin_width = 5

    required_trackers = [0, 1] # for space points
    required_number_of_track_points = 12 # doesnt do anything
    global_min_step_size = 1. # for extrapolation, set the extrapolation step size
    global_max_step_size = 100. # for extrapolation, set the extrapolation step size
    will_load_tk_space_points = True # determines whether data loader will attempt to load tracker space points
    will_load_tk_track_points = True # determines whether data loader will attempt to load tracker track points
    preanalysis_number_of_spills = 10 # number of spills to analyse during "pre-analysis"
    analysis_number_of_spills = 10 # number of spills to analyse during each "analysis" step
    number_of_spills = None # maximum number of spills to analyse for each sub-analysis
    momentum_from_tracker = True # i.e. not from TOFs
    time_from = "tof1"
    maus_version = ""

    residuals_plots_nbins = 100 # used for track extrapolation plots
    extrapolation_does_apertures = True # set to True in order to include apertures in track extrapolation
    maus_verbose_level = 5
    amplitude_max = 25

    magnet_alignment = {
        "n_events":10,
        "max_iterations":100,
        "resolution":1.,
    }

    mc_plots = {
        "mc_stations" : { # one virtual plane for each tracker view
            "tku":[53, 52, 51, 49, 48, 47, 46, 45, 44, 42, 41, 40, 38, 37, 36],
            "tkd":[63, 64, 65, 66, 67, 68, 70, 71, 72, 74, 75, 76, 77, 78, 79],
        }
    }

    # z position of central absorber (used for offsetting
    z_apertures = 0.
    # z position of detectors (used for track extrapolation) (z, name)
    tkd_offset = 81.#81.
    detectors = [
        (5287.2, None, "tof0"),
        (12929.6, None, "tof1"),
        (13968.0, None, "tku_5"),
        (14318.0, None, "tku_4"),
        (14618.0, None, "tku_3"),
        (14867.0, None, "tku_2"),
        (15068.0, None, "tku_tp"),
        (18756.+tkd_offset, None, "tkd_tp"),
        (18855.+tkd_offset, None, "tkd_2"),
        (19205.+tkd_offset, None, "tkd_3"),
        (19505.+tkd_offset, None, "tkd_4"),
        (19855.+tkd_offset, None, "tkd_5"),
        (21139.4, None, "tof2"),
    ]


    z_afc = 16955.74
    # z position of apertures (z, maximum radius, name)
    # Notes from Jason: 209.6 to fixed flange
    # 282.5 to movable flange
    # what about tracker butterfly us and ds?
    apertures = sorted(
      [(-3325.74+z_afc, 100., "diffuser_us"),
       (-3276.74+z_afc, 100., "diffuser_mid"),
       (-3227.74+z_afc, 100., "diffuser_ds"),]+\
      [(float(z)+z_afc, 150., "tku_"+str(z)) for z in range(-2918, -1817, 100)]+\
      [(float(z)+z_afc, 150., "tkd_"+str(z)) for z in range(1818, 2919, 100)]+\
      [(float(z)+z_afc, 200., "pipe_"+str(z)) for z in range(-1800, 1801, 100)]+\
      [(+209.6+z_afc, 160., "afc_209.6")],
    )


