import copy

def mc_file_names(job_name, datasets):
    file_list = ["/home/cr67/MAUS/work/reco/mc/"+job_name+"/"+datasets+"/*_sim.root"]
    return file_list

def reco_file_names(run_number_list, maus):
    file_list = []
    for run in run_number_list:
        run = str(run).rjust(5, '0')
        a_file = "/home/cr67/MAUS/work/reco/"+maus+"/"+run+"/"+run+"_recon.root"
        file_list.append(a_file)
    return file_list

def get_analysis(run_list, name, tof01_max, maus_version, data_dir, extrap, amplitude_source):
    plot_dir = data_dir+"plots_"+name+"/"
    plot_dir = plot_dir.replace(" ", "_")
    return {
            "plot_dir":plot_dir, # makedirs and then put plots in this directory. Removes any old plots there!!!
            "tof12_cut_low":32., # TOF12 cut lower bound
            "tof12_cut_high":39., # TOF12 cut upper bound
            "delta_tof01_lower":-2., # Delta TOF01 cut lower bound 
            "delta_tof01_upper":1.5, # Delta TOF01 cut upper bound 
            "delta_tof12_lower":-5., # Delta TOF01 cut lower bound 
            "delta_tof12_upper":5., # Delta TOF01 cut upper bound 
            "tof01_cut_low":28., # TOF01 cut lower bound
            "tof01_cut_high":tof01_max, # TOF01 cut upper bound
            "p_bins":[[135, 145]], # set of momentum bins; for now really it is just a lower and upper bound
            "p_tot_ds_low":80., # downstream momentum cut lower bound
            "p_tot_ds_high":200., # downstream momentum cut upper bound
            "reco_files":reco_file_names(run_list, maus_version), # list of strings to be handed to glob
            "name":"2016-04 1.2 "+name, # appears on plots
            "color":4, # not used
            "pid":-13, # assume pid of tracks following TOF cut
            "pvalue_threshold":0.02, # minimum allowed pvalue for pvalue cut
            "amplitude_source":"output/2016-04_1.2_mc_ds-cuts/plots_"+amplitude_source+"/amplitude.json",
            "field_uncertainty":0.02,
            "do_magnet_alignment":False,
            "do_amplitude":True, #True,
            "do_extrapolation":False, #name == "3-140+M3-Test2",
            "do_mc":False,
            "do_plots":True,
            "csv_output_detectors":["tof1", "diffuser_us", "diffuser_mid", "diffuser_ds"], # write data at listed detector locations
            "csv_output_filename":None, #"8590_mc_extrapolated_tracks.csv", # write a summary output of data in flat text format to listed filename; set to None to do nothing
            "extrapolation_source":extrap
        }


class Config(object):
    geometry = "geometry_08681/ParentGeometryFile.dat" #"Test.dat" # "geometry_08681/ParentGeometryFile.dat" #
    # location to which data and plots will be dumped following analysis
    data_dir = "output/2016-04_1.2_reco/"
    info_file = "geometry_08681/Maus_Information.gdml"
    will_require_tof1 = True # require at least one TOF1 Space point to even load the data
    will_require_tof2 = False # require at least one TOF2 Space point to even load the data
    tk_station = 1 # load data from a particular tracker station
    tk_plane = 0
    # prerequisite for space point cut
    will_require_triplets = False #True # require triplet space points
    upstream_cuts = { # Set to true to make data_plotter and amplitude_analysis use these cuts; False to ignore the cut
          "any_cut":None,
          "scifi_space_clusters":False,
          "scifi_space_points":False,
          "scifi_tracks_us":True,
          "scifi_nan_us":True,
          "scifi_track_points_us":False,
          "aperture_us":False,
          "pvalue_us":True,
          "aperture_ds":False,
          "scifi_tracks_ds":False,
          "scifi_nan_ds":False,
          "scifi_track_points_ds":False,
          "pvalue_ds":False,
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
    downstream_cuts["tof12"] = True # if TOF12 is out of range chuck it (but ignore "no TOF2" events)
    downstream_cuts["scifi_tracks_ds"] = True
    downstream_cuts["scifi_nan_ds"] = True
    extrapolation_cuts = upstream_cuts

    analyses = []
    analyses.append(get_analysis([8681], "3-140", 32, "MAUS-v2.8.5", data_dir, "tku_tp", "3-140_MC_Prod"))
    #analyses.append(get_analysis([8699], "6-140", 31., "MAUS-v2.8.5", data_dir, "tku_tp", "6-140_MC"))
    #analyses.append(get_analysis([8685], "10-140", 30., "MAUS-v2.8.5", data_dir, "tku_tp", "10-140_MC"))
    amplitude_bin_width = 5
    amplitude_max = 25

    required_trackers = [0, 1] # for space points
    required_number_of_track_points = 12 # doesnt do anything
    global_min_step_size = 1. # for extrapolation, set the extrapolation step size
    global_max_step_size = 100. # for extrapolation, set the extrapolation step size
    will_load_tk_space_points = True # determines whether data loader will attempt to load tracker space points
    will_load_tk_track_points = True # determines whether data loader will attempt to load tracker track points
    number_of_spills = None #100 # if set to an integer, limits the number of spills loaded for each sub-analysis
    momentum_from_tracker = True # i.e. not from TOFs
    time_from = "tof1"

    residuals_plots_nbins = 100 # used for track extrapolation plots
    extrapolation_does_apertures = True # set to True in order to include apertures in track extrapolation
    maus_verbose_level = 5

    magnet_alignment = {
        "n_events":10,
        "max_iterations":100,
        "resolution":1.,
    }

    # z position of central absorber (used for offsetting
    z_apertures = 0.
    # z position of detectors (used for track extrapolation) (z, name)
    detectors = [
        (5287.2, None, "tof0"),
        (12929.6, None, "tof1"),
        (13968.0, None, "tku_5"),
        (14318.0, None, "tku_4"),
        (14618.0, None, "tku_3"),
        (14867.0, None, "tku_2"),
        (15068.0, None, "tku_tp"),
        (18836.8, None, "tkd_tp"),
        (18855., None, "tkd_2"),
        (19205., None, "tkd_3"),
        (19505., None, "tkd_4"),
        (19855., None, "tkd_5"),
        (21139.4, None, "tof2"),
    ]

    mc_plots = {
        "mc_stations" : {
            "tku":53,
            "tkd":63,
        }
    }

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


