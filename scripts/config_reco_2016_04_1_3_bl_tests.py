def file_names(run_number_list, maus="MAUS-v2.6.5"):
    file_list = []
    for run in run_number_list:
        run = str(run).rjust(5, '0')
        a_file = "/home/cr67/MAUS/work/reco/"+maus+"/"+run+"/"+run+"_recon.root"
        file_list.append(a_file)
    return file_list

class Config(object):
    geometry = "Test.dat" # "geometry_08445/ParentGeometryFile.dat" # "scripts/fit_transfer_matrix/Step_IV_Coils.dat" #
    data_dir = "output/2016-04_1.3/"
    # prerequisite for data to load
    will_require_tof1 = True
    will_require_tof2 = False
    tk_station = 5
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
          "tof_2_sp":False,
    }

    analyses = [
        {
            "plot_dir":data_dir+"plots_140MeVc-6/",
            "tof12_cut_low":32.,
            "tof12_cut_high":39.,
            "tof01_cut_low":28.,
            "tof01_cut_high":33.,
            "p_bins":[[137.5, 142.5]],
            "p_tot_ds_low":80.,
            "p_tot_ds_high":200.,
            "reco_files":file_names([8593, 8595, 8596, 8597, 8598, 8599, 8600, 8602], "MAUS-v2.6.5"),
            "name":"6-140+M3-Test2",
            "color":4,
            "pid":-13, # tof cut -> muon peak
            "do_amplitude":False,
            "do_extrapolate":False,
        },
    ]
    void = [{
            "plot_dir":data_dir+"plots_140MeVc-10-140+M3-Test1-diff4-8585/",
            "tof12_cut_low":32.,
            "tof12_cut_high":39.,
            "tof01_cut_low":28.,
            "tof01_cut_high":33.,
            "p_bins":[[135, 145]],
            "p_tot_ds_low":80.,
            "p_tot_ds_high":200.,
            "reco_files":file_names([8585,], "MAUS-v2.6.5"),
            "name":"10-140+M3-Test1-diff4-8585",
            "color":4,
            "pid":-13, # tof cut -> muon peak
        },
        {
            "plot_dir":data_dir+"plots_140MeVc-10-140+M3-Test1-diff4-8584/",
            "tof12_cut_low":32.,
            "tof12_cut_high":39.,
            "tof01_cut_low":28.,
            "tof01_cut_high":33.,
            "p_bins":[[135, 145]],
            "p_tot_ds_low":80.,
            "p_tot_ds_high":200.,
            "reco_files":file_names([8584,], "MAUS-v2.6.5"),
            "name":"10-140+M3-Test1-diff4-8584",
            "color":4,
            "pid":-13, # tof cut -> muon peak
        },
        {
            "plot_dir":data_dir+"plots_140MeVc-10-140+M3-Test1-diff4-8583/",
            "tof12_cut_low":32.,
            "tof12_cut_high":39.,
            "tof01_cut_low":28.,
            "tof01_cut_high":33.,
            "p_bins":[[135, 145]],
            "p_tot_ds_low":80.,
            "p_tot_ds_high":200.,
            "reco_files":file_names([8583,], "MAUS-v2.6.5"),
            "name":"10-140+M3-Test1-diff4-8583",
            "color":4,
            "pid":-13, # tof cut -> muon peak
        },
        {
            "plot_dir":data_dir+"plots_140MeVc-10-140+M3-Test1-diff4-8582/",
            "tof12_cut_low":32.,
            "tof12_cut_high":39.,
            "tof01_cut_low":28.,
            "tof01_cut_high":33.,
            "p_bins":[[135, 145]],
            "p_tot_ds_low":80.,
            "p_tot_ds_high":200.,
            "reco_files":file_names([8582,], "MAUS-v2.6.5"),
            "name":"10-140+M3-Test1-diff4-8582",
            "color":4,
            "pid":-13, # tof cut -> muon peak
        },
        {
            "plot_dir":data_dir+"plots_140MeVc-10-140+M3-Test1-diff4-8581/",
            "tof12_cut_low":32.,
            "tof12_cut_high":39.,
            "tof01_cut_low":28.,
            "tof01_cut_high":33.,
            "p_bins":[[135, 145]],
            "p_tot_ds_low":80.,
            "p_tot_ds_high":200.,
            "reco_files":file_names([8581,], "MAUS-v2.6.5"),
            "name":"10-140+M3-Test1-diff4-8581",
            "color":4,
            "pid":-13, # tof cut -> muon peak
        },
    ]
    void = [
        {
            "plot_dir":data_dir+"plots_140MeVc-10-q789+150pc/",
            "tof12_cut_low":32.,
            "tof12_cut_high":39.,
            "tof01_cut_low":28.,
            "tof01_cut_high":33.,
            "p_bins":[[135, 145]],
            "p_tot_ds_low":80.,
            "p_tot_ds_high":200.,
            "reco_files":file_names([8555,]),
            "name":"10-140+M3-Test1-q789+150pc",
            "color":4,
            "pid":-13, # tof cut -> muon peak
        },
        {
            "plot_dir":data_dir+"plots_140MeVc-10-q789+100pc/",
            "tof12_cut_low":32.,
            "tof12_cut_high":39.,
            "tof01_cut_low":28.,
            "tof01_cut_high":33.,
            "p_bins":[[135, 145]],
            "p_tot_ds_low":80.,
            "p_tot_ds_high":200.,
            "reco_files":file_names([8554,]),
            "name":"10-140+M3-Test1-q789+100pc",
            "color":4,
            "pid":-13, # tof cut -> muon peak
        },
        {
            "plot_dir":data_dir+"plots_140MeVc-10-q789+50pc/",
            "tof12_cut_low":32.,
            "tof12_cut_high":39.,
            "tof01_cut_low":28.,
            "tof01_cut_high":33.,
            "p_bins":[[135, 145]],
            "p_tot_ds_low":80.,
            "p_tot_ds_high":200.,
            "reco_files":file_names([8553,]),
            "name":"10-140+M3-Test1-q789+50pc",
            "color":4,
            "pid":-13, # tof cut -> muon peak
        },
        {
            "plot_dir":data_dir+"plots_140MeVc-10-q789-25pc/",
            "tof12_cut_low":32.,
            "tof12_cut_high":39.,
            "tof01_cut_low":28.,
            "tof01_cut_high":33.,
            "p_bins":[[135, 145]],
            "p_tot_ds_low":80.,
            "p_tot_ds_high":200.,
            "reco_files":file_names([8552,]),
            "name":"10-140+M3-Test1-q789-25pc",
            "color":4,
            "pid":-13, # tof cut -> muon peak
        },
        {
            "plot_dir":data_dir+"plots_140MeVc-10-q789-50pc/",
            "tof12_cut_low":32.,
            "tof12_cut_high":39.,
            "tof01_cut_low":28.,
            "tof01_cut_high":33.,
            "p_bins":[[135, 145]],
            "p_tot_ds_low":80.,
            "p_tot_ds_high":200.,
            "reco_files":file_names([8551,]),
            "name":"10-140+M3-Test1-q789-50pc",
            "color":4,
            "pid":-13, # tof cut -> muon peak
        },
        {
            "plot_dir":data_dir+"plots_140MeVc-10-q789-10pc/",
            "tof12_cut_low":32.,
            "tof12_cut_high":39.,
            "tof01_cut_low":28.,
            "tof01_cut_high":33.,
            "p_bins":[[135, 145]],
            "p_tot_ds_low":80.,
            "p_tot_ds_high":200.,
            "reco_files":file_names([8547,]),
            "name":"10-140+M3-Test1-q789-10pc",
            "color":4,
            "pid":-13, # tof cut -> muon peak
        },
        {
            "plot_dir":data_dir+"plots_140MeVc-10-q789+10pc/",
            "tof12_cut_low":32.,
            "tof12_cut_high":39.,
            "tof01_cut_low":28.,
            "tof01_cut_high":33.,
            "p_bins":[[135, 145]],
            "p_tot_ds_low":80.,
            "p_tot_ds_high":200.,
            "reco_files":file_names([8546,]),
            "name":"10-140+M3-Test1-q789+10pc",
            "color":4,
            "pid":-13, # tof cut -> muon peak
        },
        {
            "plot_dir":data_dir+"plots_140MeVc-10-var1/",
            "tof12_cut_low":32.,
            "tof12_cut_high":39.,
            "tof01_cut_low":28.,
            "tof01_cut_high":33.,
            "p_bins":[[135, 145]],
            "p_tot_ds_low":80.,
            "p_tot_ds_high":200.,
            "reco_files":file_names([8545,]),
            "name":"10-140+M3-Test1-D1-diff",
            "color":4,
            "pid":-13, # tof cut -> muon peak
        },
        {
            "plot_dir":data_dir+"plots_140MeVc-10/",
            "tof12_cut_low":32.,
            "tof12_cut_high":39.,
            "tof01_cut_low":28.,
            "tof01_cut_high":33.,
            "p_bins":[[135, 145]],
            "p_tot_ds_low":80.,
            "p_tot_ds_high":200.,
            "reco_files":file_names([8526, 8549,]),
            "name":"10-140+M3-Test1",
            "color":4,
            "pid":-13, # tof cut -> muon peak
        },]
    void = [{
            "plot_dir":data_dir+"plots_140MeVc-6/",
            "tof12_cut_low":32.,
            "tof12_cut_high":39.,
            "tof01_cut_low":28.,
            "tof01_cut_high":33.,
            "p_bins":[[137.5, 142.5]],
            "p_tot_ds_low":80.,
            "p_tot_ds_high":200.,
            "reco_files":file_names([8519, 8520, 8521, 8522]),
            "name":"6-140+M3-Test1",
            "color":4,
            "pid":-13, # tof cut -> muon peak
        },
        {
            "plot_dir":data_dir+"plots_140MeVc-3/",
            "tof12_cut_low":32.,
            "tof12_cut_high":39.,
            "tof01_cut_low":28.,
            "tof01_cut_high":33.,
            "p_bins":[[135, 145]],
            "p_tot_ds_low":80.,
            "p_tot_ds_high":200.,
            "reco_files":file_names([8518, 8529, 8530, 8531, 8532, 8537, 8538, 8541, 8542, 8543, 8544]),
            "name":"3-140+M3-Test2",
            "color":4,
            "pid":-13, # tof cut -> muon peak
        },
    ]
    required_trackers = [0, 1] # for space points
    required_number_of_track_points = 12 # doesnt do anything
    eloss_forwards_model = "bethe_bloch_forwards"
    eloss_backwards_model = "bethe_bloch_backwards"
    will_load_tk_space_points = False
    will_load_tk_track_points = True
    number_of_spills = 100000
    momentum_from_tracker = True

    z_tof1 = 12929.7
    z_tof2 = 21137.5
    z_tku = 15062.0
    z_tkd = 18848.5


