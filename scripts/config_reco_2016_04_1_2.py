def file_names(run_number_list, maus="MAUS-v2.6.5"):
    file_list = []
    for run in run_number_list:
        run = str(run).rjust(5, '0')
        a_file = "/home/cr67/MAUS/work/reco/"+maus+"/"+run+"/"+run+"_recon.root"
        file_list.append(a_file)
    return file_list

class Config(object):
    geometry = "Step_IV_Coils_2016_04_1.2.dat" # "Test.dat" # 
    # location to which data and plots will be dumped following analysis
    data_dir = "output/2016-04_1.2/"
    info_file = "geometry_08445/Maus_Information.gdml"
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
          "scifi_track_points_us":True,
          "scifi_tracks_ds":False,
          "scifi_track_points_ds":False,
          "tof01":True,
          "tof12":False,
          "p_tot_us":True,
          "p_tot_ds":True,
          "tof_0_sp":True,
          "tof_1_sp":True,
          "tof_2_sp":True,
    }

    analyses =[{
            "plot_dir":data_dir+"plots_140MeVc-3/", # makedirs and then put plots in this directory. Removes any old plots there!!!
            "tof12_cut_low":32., # TOF12 cut lower bound
            "tof12_cut_high":39., # TOF12 cut upper bound
            "tof01_cut_low":28., # TOF01 cut lower bound
            "tof01_cut_high":32., # TOF01 cut upper bound
            "p_bins":[[135, 145]], # set of momentum bins; for now really it is just a lower and upper bound
            "p_tot_ds_low":80., # downstream momentum cut lower bound
            "p_tot_ds_high":200., # downstream momentum cut upper bound
            "reco_files":file_names([8643,], "MAUS-v2.6.5"), # list of strings to be handed to glob
            "name":"3-140+M3-Test2", # appears on plots
            "color":4, # not used
            "pid":-13, # assume pid of tracks following TOF cut
            "do_amplitude":False,
            "do_extrapolation":False,
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

    # z position of detectors (used for track extrapolation)
    z_tof0 = 5287.2
    z_tof1 = 12929.6
    z_tku_1 = 15062.0
    z_tkd_1 = 18849.5
    z_tkd_2 = 19049.1
    z_tkd_3 = 19299.0
    z_tkd_4 = 19599.0
    z_tkd_5 = 19948.9
    z_tof2 = 21139.4
    z_tof1 = 12929.7
    z_tof2 = 21137.5


