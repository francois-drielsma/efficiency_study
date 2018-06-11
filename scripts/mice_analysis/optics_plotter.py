import copy
import sys
import numpy
import ROOT
import xboa.common

from analysis_base import AnalysisBase

class OpticsPlotter(AnalysisBase):
    def __init__(self, config, config_anal, data_loader):
        super(OpticsPlotter, self).__init__(config, config_anal, data_loader)
        self.data_loader = data_loader
        self.ellipse_dict = {}
        self.var_out = {
            "tku_tp":{"ds":[],},
            "tkd_tp":{"ds":[],},
        }
        self.min_z_us = 12900. # BUG should come from config
        self.min_z_ds = 18800. # BUG should come from config
        for z_pos, dummy, detector in self.config.virtual_detectors:
            if z_pos < self.min_z_us:
                continue
            self.var_out["global_through_"+detector] = {"ds":[]}
            if z_pos < self.min_z_ds:
                continue
            self.var_out["global_ds_"+detector] = {"ds":[]}
        self.refaffed_ellipse_dict = {}

    def birth(self):
        self.set_plot_dir("optics_plots")
        self.get_var_list()
        for detector in self.var_out:
            for cut in self.var_out[detector]:
                self.birth_ellipse(detector, cut)
                self.process_ellipse(detector, cut)

    def process(self):
        self.get_var_list()
        for detector in self.var_out:
            for cut in self.var_out[detector]:
                self.process_ellipse(detector, cut)

    def death(self):
        us_lambda = lambda detector_name: detector_name == "tku_tp" or "global_through" in detector_name
        ds_lambda = lambda detector_name: detector_name == "tkd_tp" or "global_ds" in detector_name
        for prefix, detector_lambda, color in ("source_tku", us_lambda, ROOT.kBlue), ("source_tkd", ds_lambda, ROOT.kRed):
            self.refaff_ellipse_dict(detector_lambda) # rotate the data structure view
            for var, sub_var in [
                            ("mean", 5),
                            ("mean", 4),
                            ("beta_4d", None),
                            ("beta_x", None),
                            ("beta_y", None),
                            ("l_kin", None),
                            ("nevents", None),
                            ("sigma", 0),
                            ("sigma", 2),
                        ]:
                try:
                    self.plot("ds", var, sub_var, color, prefix) # mean z
                except Exception:
                    sys.excepthook(*sys.exc_info())
        self.base_death()
        self.print_plots()

    def will_cut_us(self, event):
        return event["upstream_cut"]

    def will_cut_ds(self, event):
        return event["downstream_cut"]

    def will_cut_ex(self, event):
        return event["extrapolation_cut"]

    def get_ellipse(self, detector, cut):
        return self.ellipse_dict[detector][cut]

    def clear_var_out(self):
        for det in self.var_out:
            for cut in self.var_out[det]:
                self.var_out[det][cut] = []

    def get_var_list(self):
        self.clear_var_out()
        top_detector_list = set(self.var_out.keys())
        for event in self.data_loader.events:
            detector_list = copy.deepcopy(top_detector_list)
            for detector_hit in event["data"]:
                detector = detector_hit["detector"]
                if detector not in detector_list:
                    continue
                detector_list.remove(detector) # only add once for each detector
                cuts = self.var_out[detector].keys()
                hit = detector_hit["hit"]
                this_var_list = [hit[var] for var in self.ellipse_variables]
                if "all" in cuts and not self.will_cut_us(event):
                    self.var_out[detector]["all"].append(this_var_list)
                if "us" in cuts and not self.will_cut_ds(event):
                    self.var_out[detector]["us"].append(this_var_list)
                if "ds" in cuts and not self.will_cut_ds(event):
                    self.var_out[detector]["ds"].append(this_var_list)

    def birth_ellipse(self, detector, cut):
        n_var = len(self.ellipse_variables)
        ellipse = {
        }
        ellipse["covariance"] = [[0. for j in range(n_var)] for i in range(n_var)]  
        ellipse["mean"] = [0. for i in range(n_var)]
        ellipse["sigma"] = [0. for i in range(n_var)]
        ellipse["nevents"] = 0.
        ellipse["emit_4d"] = 0.
        ellipse["beta_4d"] = 0.
        ellipse["alpha_4d"] = 0.
        ellipse["l_kin"] = 0.
        ellipse["z"] = 0.
        for axis in ["x", "y"]:
            ellipse["emit_"+axis] = 0.
            ellipse["beta_"+axis] = 0.
            ellipse["alpha_"+axis] = 0.
        if detector not in self.ellipse_dict:
            self.ellipse_dict[detector] = {}
        self.ellipse_dict[detector][cut] = ellipse

    @classmethod
    def print_ellipse(cls, ellipse):
        for row in ellipse:
          for cell in row:
            print str(round(cell)).rjust(8),
          print

    def process_ellipse(self, detector, cut):
        ellipse = self.get_ellipse(detector, cut)
        mean = ellipse["mean"]
        rms_ellipse = ellipse["covariance"]
        n_events = ellipse["nevents"]
        n_var = len(self.ellipse_variables)
        m_events = 0
        this_matrix = [[0. for j in range(n_var)] for i in range(n_var)]       
        this_mean = [0. for i in range(n_var)]
        all_data = self.var_out[detector][cut]
        for data in all_data:
            m_events += 1
            for i in range(n_var):
                this_mean[i] += data[i]
                for j in range(i, n_var):
                    this_matrix[i][j] += data[i]*data[j]
        if m_events+n_events == 0:
            return

        # update the main ellipse
        for i in range(n_var):
            for j in range(i, n_var):
                rms_ellipse[i][j] += mean[i]*mean[j]
 
        for i in range(n_var):
            mean[i] = mean[i]*n_events/(n_events+m_events) + \
                      this_mean[i]/(n_events+m_events)
            for j in range(i, n_var):
                rms_ellipse[i][j] = rms_ellipse[i][j]*n_events/(n_events+m_events) + \
                                this_matrix[i][j]/(n_events+m_events)
                rms_ellipse[i][j] -= mean[i]*mean[j]
                rms_ellipse[j][i] = rms_ellipse[i][j]

        ellipse["mean"] = mean
        ellipse["sigma"] = [abs(rms_ellipse[i][i])**0.5 for i in range(n_var)] # abs(sigma) as sigma_z == 0 causes nan error
        ellipse["covariance"] = rms_ellipse
        ellipse["nevents"] += m_events
        matrix = numpy.array(rms_ellipse)[0:4, 0:4]
        det = numpy.linalg.det(matrix)
        if det < 1e-9:
            # right at the back, near EMR, beam becomes very upright and det < 0
            return
        mu_mass = xboa.common.pdg_pid_to_mass[13]
        emit = det**0.25/mu_mass
        if abs(emit) > 1e-9:
            beta = (matrix[0][0]+matrix[2][2])*mean[4]/2./mu_mass/emit
            alpha = (matrix[0][1]+matrix[2][3])/2./mu_mass/emit
            ellipse["beta_4d"] = beta
            ellipse["alpha_4d"] = alpha
        l_kin = matrix[0][3]-matrix[1][2] # XBOA
        ellipse["emit_4d"] = emit
        for axis, index in [("x", 0), ("y", 2)]:
            matrix = numpy.array(rms_ellipse)[index:index+2, index:index+2]
            det = numpy.linalg.det(matrix)
            if det < 1e-9:
                # right at the back, near EMR, beam becomes very upright and det < 0
                continue
            emit = det**0.5/mu_mass
            beta = (matrix[0][0])*mean[4]/mu_mass/emit
            alpha = (matrix[0][1])/mu_mass/emit
            ellipse["beta_"+axis] = beta
            ellipse["alpha_"+axis] = alpha
            ellipse["emit_"+axis] = emit

    def death_ellipse(self, tracker, cut, print_to_screen):
        return

    def refaff_ellipse_dict(self, detector_lambda):
        """
        detector_lambda: lambda function that returns true if detector name is okay
        """
        self.refaffed_ellipse_dict = {
            "all":[],
            "ex":[],
            "us":[],
            "ds":[]
        }
        for detector in self.ellipse_dict:
            if not detector_lambda(detector):
                continue
            for cut in self.ellipse_dict[detector]:
                ellipse = self.ellipse_dict[detector][cut]
                self.refaffed_ellipse_dict[cut].append(ellipse)
        for cut in self.refaffed_ellipse_dict:
            z_var = self.ellipse_variables.index("z")
            refaff = sorted(self.refaffed_ellipse_dict[cut],
                            key = lambda ellipse: ellipse["mean"][z_var])
            self.refaffed_ellipse_dict[cut] = refaff

    @classmethod
    def name_lookup(cls, var, sub_var, cut):
        if sub_var == None:
            name = var+"_"+cut
            axis = var+" ("+cut+")"
        else:
            name = str(var)+"_"+str(sub_var)+"_"+cut
            axis = str(var)+" "+str(sub_var)+" ("+cut+")"
        return name, axis

    def plot(self, cut, var, sub_var, color, prefix):
        ellipse_list = self.refaffed_ellipse_dict[cut]
        z_var = self.ellipse_variables.index("z")
        pred =  lambda ellipse: ellipse["mean"][z_var] >= self.min_z_us
        z_list = [ellipse["mean"][z_var]*1e-3 for ellipse in ellipse_list if pred(ellipse)]
        if sub_var == None:
            var_list = [ellipse[var] for ellipse in ellipse_list if pred(ellipse)]
        else:
            var_list = [ellipse[var][sub_var] for ellipse in ellipse_list if pred(ellipse)]
        if "beta" in var:
            var_list = [ellipse[var][sub_var]*1e-3 for ellipse in ellipse_list if pred(ellipse)]
        name, axis = self.name_lookup(var, sub_var, cut)
        hist, graph = self.make_root_graph(name, name+"_"+prefix,
                      z_list, "z [m]", var_list, axis, True,
                      None, None, None, None)
        if len(self.get_plot(name)["histograms"]) == 1:
            hist.Draw()
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(color)
        graph.Draw("p l same")
        det_list = self.config.detectors
        det_list += [virt for virt in self.config.virtual_detectors if virt[2] == "virtual_absorber_centre"]
        for z_det, dummy, detector in self.config.detectors:
            var_min, var_max = min(var_list), max(var_list)
            delta = var_max-var_min+1
            var_max += delta
            var_min -= delta
            hist, graph = self.make_root_graph(name, name+"_"+detector,
                      [z_det*1e-3, z_det*1e-3], "z [m]", [var_min, var_max], axis, True,
                      None, None, None, None)
            if "tku" in detector or "tkd" in detector:
                line_color = ROOT.kBlue
            elif "tof" in detector or "cal" in detector:
                line_color = ROOT.kRed
            else:
                line_color = ROOT.kGreen
            graph.SetLineColor(line_color)
            graph.SetLineStyle(2)
            graph.Draw("same l")

    ellipse_variables = ["x", "px", "y", "py", "p", "z"]

