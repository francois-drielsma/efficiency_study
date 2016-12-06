def mc_file_names(job_name, datasets):
    file_list = ["/home/cr67/MAUS/work/reco/mc/"+job_name+"/"+datasets+"/*_sim.root"]
    return file_list

def reco_file_names(run_number_list, maus="MAUS-v2.6.5"):
    file_list = []
    for run in run_number_list:
        run = str(run).rjust(5, '0')
        a_file = "/home/cr67/MAUS/work/reco/"+maus+"/"+run+"/"+run+"_recon.root"
        file_list.append(a_file)
    return file_list

class Config(object):
    geometry = "geometry_08590/ParentGeometryFile.dat" # "Step_IV_Coils_2016_04_1.2.dat" # "Test.dat" # "Test.dat" #
    # location to which data and plots will be dumped following analysis
    data_dir = "output/8681/"
    info_file = "geometry_08681/Maus_Information.gdml"
    will_require_tof1 = True # require at least one TOF1 Space point to even load the data
    will_require_tof2 = False # require at least one TOF2 Space point to even load the data
    tk_station = 1 # load data from a particular tracker station
    # prerequisite for space point cut
    will_require_triplets = False #True # require triplet space points
    cuts_active = { # Set to true to make data_plotter and amplitude_analysis use these cuts; False to ignore the cut
          "any_cut":None,
          "scifi_space_clusters":False,
          "scifi_space_points":False,
          "scifi_tracks_us":True,
          "scifi_track_points_us":True,
          "aperture_us":False,
          "aperture_ds":True,
          "scifi_tracks_ds":False,
          "scifi_track_points_ds":False,
          "tof01":True,
          "tof12":False,
          "p_tot_us":True,
          "p_tot_ds":False,
          "tof_0_sp":True,
          "tof_1_sp":True,
          "tof_2_sp":False,
    }
    extrapolation_cuts = { # Set to true to make extrapolate_tracks use the cuts; False to ignore the cut
          "any_cut":False,
          "scifi_space_clusters":False,
          "scifi_space_points":False,
          "scifi_tracks_us":True,
          "scifi_track_points_us":False,
          "aperture_us":False,
          "aperture_ds":True,
          "scifi_tracks_ds":False,
          "scifi_track_points_ds":False,
          "tof01":False,
          "tof12":False,
          "p_tot_us":False,
          "p_tot_ds":False,
          "tof_0_sp":False,
          "tof_1_sp":True,
          "tof_2_sp":True,
    }

    analyses =[{
            "plot_dir":data_dir+"plots_data/", # makedirs and then put plots in this directory. Removes any old plots there!!!
            "tof12_cut_low":32., # TOF12 cut lower bound
            "tof12_cut_high":39., # TOF12 cut upper bound
            "tof01_cut_low":28., # TOF01 cut lower bound
            "tof01_cut_high":32., # TOF01 cut upper bound
            "p_bins":[[135, 145]], # set of momentum bins; for now really it is just a lower and upper bound
            "p_tot_ds_low":80., # downstream momentum cut lower bound
            "p_tot_ds_high":200., # downstream momentum cut upper bound
            "reco_files":reco_file_names([8681], "MAUS-v2.6.5"), # mc_file_names("000025", "0000?"), # list of strings to be handed to glob
            "name":"Run 8681", # appears on plots
            "color":4, # not used
            "pid":-13, # assume pid of tracks following TOF cut
            "do_extrapolation":True,
            "do_amplitude":False,
            "csv_output_detectors":["tof1", "diffuser_us", "diffuser_mid", "diffuser_ds"], # write data at listed detector locations
            "csv_output_filename":None #"8590_mc_extrapolated_tracks.csv", # write a summary output of data in flat text format to listed filename; set to None to do nothing
        },
    ]
    void = [{
            "plot_dir":data_dir+"plots_mc/", # makedirs and then put plots in this directory. Removes any old plots there!!!
            "tof12_cut_low":32., # TOF12 cut lower bound
            "tof12_cut_high":39., # TOF12 cut upper bound
            "tof01_cut_low":27., # TOF01 cut lower bound
            "tof01_cut_high":30., # TOF01 cut upper bound
            "p_bins":[[190, 210]], # set of momentum bins; for now really it is just a lower and upper bound
            "p_tot_ds_low":80., # downstream momentum cut lower bound
            "p_tot_ds_high":200., # downstream momentum cut upper bound
            "reco_files":mc_file_names("000025", "0000?"), # list of strings to be handed to glob
            "name":"8590-mc", # appears on plots
            "color":4, # not used
            "pid":-13, # assume pid of tracks following TOF cut
            "do_amplitude":False, # set to True to do amplitude calculation
            "do_extrapolation":True, # set to True to do extrapolation to different detectors
            "csv_output_detectors":["tof1", "diffuser_us", "diffuser_mid", "diffuser_ds"], # write data at listed detector locations
            "csv_output_filename":"8590_mc_extrapolated_tracks.csv", # write a summary output of data in flat text format to listed filename; set to None to do nothing
        },
    ]
    required_trackers = [0, 1] # for space points
    required_number_of_track_points = 12 # doesnt do anything
    global_min_step_size = 1. # for extrapolation, set the extrapolation step size
    global_max_step_size = 100. # for extrapolation, set the extrapolation step size
    will_load_tk_space_points = False # determines whether data loader will attempt to load tracker space points
    will_load_tk_track_points = True # determines whether data loader will attempt to load tracker track points
    number_of_spills = None # if set to an integer, limits the number of spills loaded for each sub-analysis
    momentum_from_tracker = True # i.e. not from TOFs

    residuals_plots_nbins = 100 # used for track extrapolation plots
    extrapolation_does_apertures = True # set to True in order to include apertures in track extrapolation

    # z position of central absorber (used for offsetting
    z_apertures = 0.
    # z position of detectors (used for track extrapolation) (z, name)
    detectors = [
        (5287.2, None, "tof0"),
        (12929.6, None, "tof1"),
        (15062.0, None, "tku_tp"),
        (18849.5, None, "tkd_tp"),
        (19049.1, None, "tkd_2"),
        (19299.0, None, "tkd_3"),
        (19599.0, None, "tkd_4"),
        (19948.9, None, "tkd_5"),
        (21139.4, None, "tof2"),
    ]

    z_afc = 16955.74
    # z position of apertures (z, maximum radius, name)
    apertures = sorted(
      [(-3325.74+z_afc, 100., "diffuser_us"),
       (-3276.74+z_afc, 100., "diffuser_mid"),
       (-3227.74+z_afc, 100., "diffuser_ds"),]+\
      [(float(z)+z_afc, 150., "tku_"+str(z)) for z in range(-2918, -1817, 100)]+\
      [(float(z)+z_afc, 150., "tkd_"+str(z)) for z in range(1818, 2919, 100)]+\
      [(float(z)+z_afc, 200., "pipe_"+str(z)) for z in range(-1800, 1801, 100)]+\
      [(+225.0+z_afc, 160., "afc_225")],
    )


