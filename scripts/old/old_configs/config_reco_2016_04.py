import copy

def mc_file_names(job_name, datasets):
    file_list = ["/home/cr67/work/reco/mc/"+job_name+"/"+datasets+"/*_sim.root"]
    return file_list

def reco_file_names(run_number_list, maus, do_globals):
    file_list = []
    for run in run_number_list:
        run = str(run).rjust(5, '0')
        a_file = "/home/cr67/work/reco/"+maus+"/"+run+"/"+run+"_recon"
        if not do_globals:
            a_file += ".root"
        else:
            a_file += "_global_test.root"
        file_list.append(a_file)
    print file_list
    return file_list

def get_analysis(run_list, name, tof01_min_max, maus_version, data_dir, amplitude_source, p_bins, tkd_cut, do_globals):
    plot_dir = data_dir+"plots_"+name+"/"
    plot_dir = plot_dir.replace(" ", "_")
    min_p = min([min(a_bin) for a_bin in p_bins])
    max_p = max([max(a_bin) for a_bin in p_bins])
    return {
            "plot_dir":plot_dir, # makedirs and then put plots in this directory. Removes any old plots there!!!
            "tof12_cut_low":32., # TOF12 cut lower bound
            "tof12_cut_high":39., # TOF12 cut upper bound
            "delta_tof01_lower":-2., # Delta TOF01 cut lower bound 
            "delta_tof01_upper":1.5, # Delta TOF01 cut upper bound 
            "delta_tof12_lower":-5., # Delta TOF01 cut lower bound 
            "delta_tof12_upper":5., # Delta TOF01 cut upper bound 
            "tof01_cut_low":tof01_min_max[0], # TOF01 cut lower bound
            "tof01_cut_high":tof01_min_max[1], # TOF01 cut upper bound
            "p_bins":p_bins, # set of momentum bins; for now really it is just a lower and upper bound
            "p_tot_ds_low":tkd_cut[0], # downstream momentum cut lower bound
            "p_tot_ds_high":tkd_cut[1], # downstream momentum cut upper bound
            "reco_files":reco_file_names(run_list, maus_version, do_globals), # list of strings to be handed to glob
            "name":name, # appears on plots
            "color":4, # not used
            "pid":-13, # assume pid of tracks following TOF cut
            "pvalue_threshold":0.02, # minimum allowed pvalue for pvalue cut
            "chi2_threshold":10.0, # maximum allowed chi2/dof for chi2 cut
            "amplitude_source":"output/2016-04_1.2_mc/plots_"+amplitude_source+"/amplitude.json",
            "field_uncertainty":0.02,
            "amplitude_chi2":False,
            "do_magnet_alignment":False,
            "do_amplitude":False,
            "do_extrapolation":False,
            "do_globals":do_globals,
            "do_mc":False,
            "do_plots":True,
            "csv_output_detectors":["tof1", "diffuser_us", "diffuser_mid", "diffuser_ds"], # write data at listed detector locations
            "csv_output_filename":"test", #"8590_mc_extrapolated_tracks.csv", # write a summary output of data in flat text format to listed filename; set to None to do nothing
            "extrapolation_source":"tku_tp"
        }


class Config(object):
    geometry = "Test.dat" #"geometry_08681/ParentGeometryFile.dat" #
    # location to which data and plots will be dumped following analysis
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
          "pvalue_us":False,
          "chi2_us":True,
          "aperture_ds":False,
          "scifi_tracks_ds":False,
          "scifi_nan_ds":False,
          "scifi_track_points_ds":False,
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
    downstream_cuts["tof2_sp"] = True
    downstream_cuts["pvalue_ds"] = False
    downstream_cuts["chi2_ds"] = True
    downstream_cuts["tof12"] = False # if TOF12 is out of range chuck it (but ignore "no TOF2" events)
    downstream_cuts["scifi_nan_ds"] = True
    downstream_cuts["scifi_tracks_ds"] = True
    extrapolation_cuts = upstream_cuts

    data_dir = "output/2016-04_1.2_globals/"

    analyses = []
    for run in [8681]:
        analyses.append(get_analysis([run], "3-140-"+str(run), [27, 32], "MAUS-v2.9.1", data_dir, "3-140_MC", [[135, 145]], [100, 200], True))
    for run in [8644]:
        analyses.append(get_analysis([run], "6-140-"+str(run), [27, 31], "MAUS-v2.9.1", data_dir, "10-140_MC", [[135, 145]], [100, 200], True))
    for run in [8645]:#8645]:
        analyses.append(get_analysis([run], "10-140-"+str(run), [27, 30], "MAUS-v2.9.1", data_dir, "10-140_MC", [[135, 145]], [100, 200], True))
    amplitude_bin_width = 5
    amplitude_max = 25

    required_trackers = [0, 1] # for space points
    required_number_of_track_points = 12 # doesnt do anything
    global_min_step_size = 1. # for extrapolation, set the extrapolation step size
    global_max_step_size = 100. # for extrapolation, set the extrapolation step size
    will_load_tk_space_points = True # determines whether data loader will attempt to load tracker space points
    will_load_tk_track_points = True # determines whether data loader will attempt to load tracker track points
    number_of_spills = None # if set to an integer, limits the number of spills loaded for each sub-analysis
    preanalysis_number_of_spills = 1000 # number of spills to analyse during "pre-analysis"
    analysis_number_of_spills = 1000 # number of spills to analyse during each "analysis" step
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
    # 13692: diffuser us edge
    # 13782: diffuser ds edge
    # 16639: lH2 window mount us (us edge)
    # 16710: lH2 window mount us (ds edge)
    # 16953.5: LiH centre
    # 17166: lH2 window mount ds (us edge)
    # 17242: lH2 window mount ds (ds edge)
    # 17585: SSU aperture
    # 18733: SSU He window
    tkd_offset = 0.
    detectors = [
        (5287.2, None, "tof0"),
        (12929.6-25., None, "tof1_us"),
        (12929.6, None, "tof1"),
        (12929.6+25., None, "tof1_ds"),
        (13968.0, None, "tku_5"),
        (14318.0, None, "tku_4"),
        (14618.0, None, "tku_3"),
        (14867.0, None, "tku_2"),
        (15068.0, None, "tku_tp"),
        (18836.8+tkd_offset, None, "tkd_tp"),
        (18855.+tkd_offset, None, "tkd_2"),
        (19205.+tkd_offset, None, "tkd_3"),
        (19505.+tkd_offset, None, "tkd_4"),
        (19855.+tkd_offset, None, "tkd_5"),
        (21114.4, None, "tof2"),
        (21139.4, None, "tof2"),
        (21159.4, None, "tof2"),
        (21208., None, "tof2"),
        (21214.4, None, "tof2"),
        (21220.4, None, "tof2"),

    ]

    plot_virtual_stations = [47, 62, 63, 72, 73, 75, 81]
    virtual_detectors = [500.*i for i in range(51)]+[det[0] for det in detectors]
    virtual_detectors += [12975, 13240., 13278, 13306.5, 13326, 13346.4, 13771.2, 13380, 13582, 13592, 13682, 13728, 13704.4, 13709.5, 13755, 13782, 13805, 13850, 13867, 16639, 16710, 16953, 16974, 17166, 17242, 17585,
                          18733, 19937.9, 16803.7, 16919.6, 16941, 16985.4, 17101.3, 18460, 18540., 18600., 19308, 19328,
                          19036.8, 19078, 19286.7, 19586.7, 20720]
    virtual_detectors = [(z, None, "virtual_"+str(i)) for i, z in enumerate(sorted(virtual_detectors))]
    for z, dummy, plane in virtual_detectors:
        print z, plane
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


