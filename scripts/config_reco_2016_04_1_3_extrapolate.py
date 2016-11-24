def file_names(run_number_list, maus="MAUS-v2.6.5"):
    file_list = []
    for run in run_number_list:
        run = str(run).rjust(5, '0')
        a_file = "/home/cr67/MAUS/work/reco/"+maus+"/"+run+"/"+run+"_recon.root"
        file_list.append(a_file)
    return file_list

class Config(object):
    geometry = "Step_IV_Coils_2016_04_1.3.dat" # "scripts/fit_transfer_matrix/Step_IV_Coils.dat" #
    info_file = "geometry_08445/Maus_Information.gdml"
    data_dir = "output/2016-04_1.3/"
    # prerequisite for data to load
    will_require_tof1 = True
    will_require_tof2 = False
    tk_station = 1
    # prerequisite for space point cut
    will_require_triplets = False #True # require triplet space points
    cuts_active = {
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
          "tof_2_sp":True,
    }
    extrapolation_cuts = {
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
          "p_tot_ds":False,
          "tof_0_sp":True,
          "tof_1_sp":True,
          "tof_2_sp":False,
    }

    analyses = [
        {
            "plot_dir":data_dir+"plots_140MeVc-6_test_extrapolate/",
            "tof12_cut_low":32.,
            "tof12_cut_high":39.,
            "tof01_cut_low":28.,
            "tof01_cut_high":33.,
            "p_bins":[[137.5, 142.5]],
            "p_tot_ds_low":80.,
            "p_tot_ds_high":200.,
            "reco_files":file_names([8593,], "MAUS-v2.6.5"), # 8595, 8596, 8597, 8598, 8599, 8600, 8602
            "name":"6-140+M3-Test2_test_extrapolate",
            "color":4,
            "pid":-13, # tof cut -> muon peak
            "do_amplitude":False,
            "do_extrapolation":True
        },
    ]
    required_trackers = [0, 1] # for space points
    required_number_of_track_points = 12 # doesnt do anything
    #eloss_forwards_model = "bethe_bloch_forwards"
    #eloss_backwards_model = "bethe_bloch_backwards"
    global_min_step_size = 1.
    global_max_step_size = 100.
    will_load_tk_space_points = False
    will_load_tk_track_points = True
    number_of_spills = 100000
    momentum_from_tracker = True

    residuals_plots_nbins = 100

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


