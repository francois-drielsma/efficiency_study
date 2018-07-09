import copy

def mc_file_names(job_name, datasets):
    file_list = ["/home/cr67/work/reco/mc/"+job_name+"/"+datasets+"/*_sim.root"]
    return file_list

def reco_file_names(run_number_list, maus, do_globals):
    file_list = []
    for run in run_number_list:
        run = str(run).rjust(5, '0')
        a_file = "/work/ast/cr67/reco/"+maus+"/"+run+"/"+run+"_recon"
        if not do_globals:
            a_file += ".root"
        else:
            a_file += ".root"
        file_list.append(a_file)
    print file_list
    return file_list

def get_amplitude_source(emittance, tk, name):
    if emittance == False or emittance == None:
        return None
    amplitude_source = "output/2017-02-Systematics/plots_Simulated_2017-2.7_"
    amplitude_source += str(emittance)+"-140_lH2_empty_Systematics_"+tk+"_"
    amplitude_source += name+"/amplitude/amplitude.json"
    return amplitude_source


def get_analysis(run_list, name, tof01_min_max, maus_version, data_dir, amplitude_source, p_bins, tkd_cut, do_globals):
    plot_dir = data_dir+"/plots_"+name+"/"
    plot_dir = plot_dir.replace(" ", "_")
    min_p = min([min(a_bin) for a_bin in p_bins])
    max_p = max([max(a_bin) for a_bin in p_bins])
    amp_systematic_sources = {}
    for tracker in ["tku", "tkd"]:
        for variable in ["pos_plus", "rot_plus", "scale_E1_plus", "scale_C_plus", "scale_E2_plus"]:
            this_source = get_amplitude_source(amplitude_source, tracker, variable)
            if this_source != None:
                amp_systematic_sources[this_source] = 1.

    return {
            "plot_dir":plot_dir, # makedirs and then put plots in this directory. Removes any old plots there!!!
            "tof12_cut_low":32., # TOF12 cut lower bound
            "tof12_cut_high":39., # TOF12 cut upper bound
            "delta_tof01_lower":-1., # Delta TOF01 cut lower bound 
            "delta_tof01_upper":+1.5, # Delta TOF01 cut upper bound 
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
            "chi2_threshold":4.0, # maximum allowed chi2/dof for chi2 cut
            "amplitude_source":get_amplitude_source(amplitude_source, "tku", "base"),
            "amplitude_systematic_reference":get_amplitude_source(amplitude_source, "tku", "base"),
            "amplitude_systematic_sources":amp_systematic_sources,
            "field_uncertainty":0.02,
            "csv_output_detectors":["tof1", "diffuser_us", "diffuser_mid", "diffuser_ds"], # write data at listed detector locations
            "csv_output_filename":"test", #"8590_mc_extrapolated_tracks.csv", # write a summary output of data in flat text format to listed filename; set to None to do nothing
            "extrapolation_source":"tku_tp",
            "amplitude_chi2":False,
            "amplitude_mc":False,
            "weight_tof01_source":None,
            "weight_tof01_target":plot_dir+"tof01_weights",
            "weight_tof01_mode":"build_distribution",
            "tracker_fiducial_radius":150.,
            "cov_fixed_us":None, #cov_us,
            "cov_fixed_ds":None, #cov_ds,
            "amplitude_algorithm":"binned",

            "do_magnet_alignment":False,
            "do_extrapolation":False,
            "do_globals":do_globals,
            "do_mc":False,
            "do_amplitude":True,
            "do_plots":True,
            "do_cuts_plots":True,
            "do_tof01_weighting":False,
            "do_optics":True,
            "do_data_recorder":True,
    }


class Config(object):
    # location to which data and plots will be dumped following analysis
    info_file = "geometry_08681/Maus_Information.gdml"
    will_require_tof1 = True # require at least one TOF1 Space point to even load the data
    will_require_tof2 = False # require at least one TOF2 Space point to even load the data
    tk_station = 1 # load data from a particular tracker station
    tk_plane = 0
    # prerequisite for space point cut
    will_require_triplets = False #True # require triplet space points
    global_through_cuts = True
    upstream_cuts = { # Set to true to make data_plotter and amplitude_analysis use these cuts; False to ignore the cut
          "any_cut":None,
          "scifi_space_clusters":False,
          "scifi_space_points":False,
          "scifi_tracks_us":True,
          "scifi_nan_us":False,
          "scifi_track_points_us":False,
          "scifi_fiducial_us":True,
          "aperture_us":False,
          "pvalue_us":False,
          "chi2_us":True,
          "aperture_ds":False,
          "scifi_tracks_ds":False,
          "scifi_nan_ds":False,
          "scifi_track_points_ds":False,
          "scifi_fiducial_ds":False,
          "pvalue_ds":False,
          "chi2_ds":False,
          "tof01":True,
          "tof12":False,
          "p_tot_us":True,
          "p_tot_ds":False,
          "tof_0_sp":True,
          "tof_1_sp":True,
          "tof_2_sp":False,
          "upstream_aperture_cut":True,
          "downstream_aperture_cut":False,
          "delta_tof01":False, #True, #extrapolatedtof01 compared to recon tof01
          "delta_tof12":False, #extrapolatedtof12 compared to recon tof12
          "global_through_tof0":global_through_cuts,
          "global_through_tof1":False,
          "global_through_tku_tp":False,
          "global_through_tkd_tp":False,
          "global_through_tof2":False,
          "tof01_selection":False,
          "mc_muon_us":False,
          "mc_stations_us":False,
          "mc_scifi_fiducial_us":False,
          "mc_muon_ds":False,
          "mc_stations_ds":False,
          "mc_scifi_fiducial_ds":False,
    }
    downstream_cuts = copy.deepcopy(upstream_cuts)
    downstream_cuts["p_tot_ds"] = True
    downstream_cuts["tof2_sp"] = False
    downstream_cuts["pvalue_ds"] = False
    downstream_cuts["chi2_ds"] = True
    downstream_cuts["scifi_fiducial_ds"] = True
    downstream_cuts["tof12"] = False # if TOF12 is out of range chuck it (but ignore "no TOF2" events)
    downstream_cuts["scifi_tracks_ds"] = True
    extrapolation_cuts = copy.deepcopy(downstream_cuts)
    extrapolation_cuts["downstream_aperture_cut"] = True
    extrapolation_cuts["tof2_sp"] = True
    extrapolation_cuts["global_through_tkd_tp"] = True
    extrapolation_cuts["global_through_tof2"] = True
    mc_true_us_cuts = copy.deepcopy(upstream_cuts)
    mc_true_ds_cuts = copy.deepcopy(upstream_cuts)
    cut_report  = [[], [], []]
    cut_report[0] = ["hline", "all events", "hline",]
    cut_report[0] += ["tof_1_sp", "tof_0_sp", "scifi_tracks_us", "chi2_us", "scifi_fiducial_us", "hline",]
    cut_report[0] += ["tof01", "p_tot_us", "hline",]
    if global_through_cuts:
        cut_report[0] += ["global_through_tof0"] #, "delta_tof01",]
    cut_report[0] += ["upstream_aperture_cut", "hline",]
    cut_report[0] += ["upstream_cut", "hline",]
    cut_report[1] += ["hline", "upstream_cut", "hline",]
    cut_report[1] += ["scifi_tracks_ds", "chi2_ds", "scifi_fiducial_ds", "p_tot_ds", "hline",]
    cut_report[1] += ["downstream_cut", "hline",]
    cut_report[2] =  ["hline", "downstream_cut", "hline",]
    cut_report[2] += ["downstream_aperture_cut", "tof_2_sp", "global_through_tkd_tp", "global_through_tof2", "hline",]
    cut_report[2] += ["extrapolation_cut", "hline"]


    data_dir = "output/2017-02-Test/" # to which data is written
    src_dir = "MAUS-Test-4-npe"
    correct_amplitude = True
    analyses = []
    # BUG need to fix the systematics "emittance"
    analyses.append(get_analysis([10069], "2017-2.7 3-140 lH2 empty", [1.5, 6.5], src_dir, data_dir, 6, [[135, 145]], [90, 170], True))
    analyses.append(get_analysis([9971],  "2017-2.7 3-140 lH2 full",  [1.5, 6.5], src_dir, data_dir, 6, [[135, 145]], [90, 170], True))
    analyses.append(get_analysis([10483], "2017-2.7 3-140 LiH",       [1.5, 6.5], src_dir, data_dir, 6, [[135, 145]], [90, 170], True))
    analyses.append(get_analysis([10444], "2017-2.7 3-140 None",      [1.5, 6.5], src_dir, data_dir, 6, [[135, 145]], [90, 170], True))

    analyses.append(get_analysis([10064], "2017-2.7 4-140 lH2 empty", [1.5, 6.0], src_dir, data_dir, 6, [[135, 145]], [90, 170], True))
    analyses.append(get_analysis([9962],  "2017-2.7 4-140 lH2 full",  [1.5, 6.0], src_dir, data_dir, 6, [[135, 145]], [90, 170], True))
    analyses.append(get_analysis([10484], "2017-2.7 4-140 LiH",       [1.5, 6.0], src_dir, data_dir, 6, [[135, 145]], [90, 170], True))
    analyses.append(get_analysis([10445], "2017-2.7 4-140 None",      [1.5, 6.0], src_dir, data_dir, 6, [[135, 145]], [90, 170], True))

    analyses.append(get_analysis([10051], "2017-2.7 6-140 lH2 empty", [1.5, 5.5], src_dir, data_dir, 6, [[135, 145]], [90, 170], True))
    analyses.append(get_analysis([9966],  "2017-2.7 6-140 lH2 full",  [1.5, 5.5], src_dir, data_dir, 6, [[135, 145]], [90, 170], True))
    analyses.append(get_analysis([10485], "2017-2.7 6-140 LiH",       [1.5, 5.5], src_dir, data_dir, 6, [[135, 145]], [90, 170], True))
    analyses.append(get_analysis([10446], "2017-2.7 6-140 None",      [1.5, 5.5], src_dir, data_dir, 6, [[135, 145]], [90, 170], True))

    analyses.append(get_analysis([10052], "2017-2.7 10-140 lH2 empty", [1.5, 4.5], src_dir, data_dir, 10, [[135, 145]], [90, 170], True))
    analyses.append(get_analysis([9970],  "2017-2.7 10-140 lH2 full",  [1.5, 4.5], src_dir, data_dir, 10, [[135, 145]], [90, 170], True))
    analyses.append(get_analysis([10486], "2017-2.7 10-140 LiH",       [1.5, 4.5], src_dir, data_dir, 10, [[135, 145]], [90, 170], True))
    analyses.append(get_analysis([10447], "2017-2.7 10-140 None",      [1.5, 4.5], src_dir, data_dir, 10, [[135, 145]], [90, 170], True))
    amplitude_bin_width = 5
    amplitude_max = 25

    required_trackers = [0, 1] # for space points
    required_number_of_track_points = 12 # doesnt do anything
    global_min_step_size = 1. # for extrapolation, set the extrapolation step size
    global_max_step_size = 100. # for extrapolation, set the extrapolation step size
    will_load_tk_space_points = True # determines whether data loader will attempt to load tracker space points
    will_load_tk_track_points = True # determines whether data loader will attempt to load tracker track points
    number_of_spills = None # if set to an integer, limits the number of spills loaded for each sub-analysis
    preanalysis_number_of_spills = 100 # number of spills to analyse during "pre-analysis"
    analysis_number_of_spills = 100 # number of spills to analyse during each "analysis" step
    momentum_from_tracker = True # i.e. not from TOFs
    time_from = "tof1"
    tof0_offset = 25.4
    tof1_offset = 0.
    tof2_offset = -27.9
    #z_tof2 - ztof1 = 21139.4-12929.6 = 8209.800000000001
    #dt_tof2 = 8209.8/299.8 = 27.384256170780517
    #compare with position of the electron peak dt_tof2 = 27.793

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
    tku_offset = -3.
    tkd_offset = 8.
    detectors = [
        (5287.2-25., None, "tof0_us"),
        (5287.2, None, "tof0"),
        (5287.2+25., None, "tof0_ds"),
        (12929.6-25., None, "tof1_us"),
        (12929.6, None, "tof1"),
        (12929.6+25., None, "tof1_ds"),
        (15068.0-1100.+tku_offset, None, "tku_5"),
        (15068.0-750.+tku_offset, None, "tku_4"),
        (15068.0-450.+tku_offset, None, "tku_3"),
        (15068.0-200.+tku_offset, None, "tku_2"),
        (15068.0+tku_offset, None, "tku_tp"),
        (18836.8+tkd_offset, None, "tkd_tp"),
        (18836.8+200.+tkd_offset, None, "tkd_2"),
        (18836.8+450.+tkd_offset, None, "tkd_3"),
        (18836.8+750.+tkd_offset, None, "tkd_4"),
        (18836.8+1100.+tkd_offset, None, "tkd_5"),
        (21114.4, None, "tof2_us"),
        (21139.4, None, "tof2"),
        (21164.4, None, "tof2_ds"),
        (21208., None, "cal"),
        (21214.4, None, "cal"),
        (21220.4, None, "cal"),
    ]

    virtual_detectors = [
        (12904., None, "virtual_tof1_us"),
        (12929.6, None, "virtual_tof1"),
        (12954., None, "virtual_tof1_ds"),
        (13205., None, "virtual_us_pry_us"),
        (13305., None, "virtual_us_pry_ds"),
        (13403., None, "virtual_tku_butterfly_us"),
        (13486., None, "virtual_tku_butterfly_ds"),
        (13578., None, "virtual_diffuser_shield"),
        (13698., None, "virtual_diffuser_us"),
        (13782., None, "virtual_diffuser_ds"),
        (16639., None, "virtual_lh2_us_window_flange_1"),
        (16710., None, "virtual_lh2_us_window_flange_2"),
        (16953.5, None, "virtual_absorber_centre"),
        (17166., None, "virtual_lh2_ds_window_flange_1"),
        (17242., None, "virtual_lh2_ds_window_flange_2"),
        (17585., None, "virtual_ssd_aperture"),
        (18733., None, "virtual_ssd_he"),
        (20447., None, "virtual_ds_pry"),
    ]
    plot_virtual_stations = [station for dummy, dummy, station in virtual_detectors]
    virtual_detectors += [(500.*i, None, "virtual_"+str(i)) for i in range(51)]
    virtual_detectors += [(z, None, "virtual_"+det) for z, dummy, det in detectors]
    virtual_detectors = sorted(virtual_detectors)
    for z, dummy, plane in virtual_detectors:
        print z, plane
    upstream_aperture_cut = {
        "global_through_virtual_diffuser_us":100.,
        "global_through_virtual_diffuser_ds":100.,
    }
    downstream_aperture_cut = {
        "global_through_virtual_lh2_us_window_flange_1":100.,
        "global_through_virtual_lh2_us_window_flange_2":100.,
        "global_through_virtual_absorber_centre":100.,
        "global_through_virtual_lh2_ds_window_flange_1":100.,
        "global_through_virtual_lh2_ds_window_flange_2":100.,
        "global_through_virtual_ssd_aperture":100.,
        "global_through_virtual_ssd_he":100.,
        "global_ds_virtual_ds_pry":100.,
    }

    mc_plots = {
        "mc_stations" : {
            "tku":"virtual_tku_tp",
            "tkd":"virtual_tkd_tp",
        }
    }
    bz_tku = 3e-3
    bz_tkd = -2e-3
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

