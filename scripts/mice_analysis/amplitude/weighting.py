import xboa.common
import ROOT

import utilities.root_style

class Weighting(object):
    def __init__(self, reco_data, truth_data, plot_dir):
        self.reco_data = reco_data
        self.truth_data = truth_data
        self.plot_dir = plot_dir

    def calculate_weighting(self):
        pass

    def apply_weighting(self, psv_list):
        pass

    def plot_slice(self):
        """For each bin, plot the ratio of kde reco vs truth with given third var"""
        pass

    def plot_sum(self, x_axis, nx_bins, x_min_max, y_axis, ny_bins, y_min_max):
        """For each bin, plot the number of events reco vs truth"""
        x_index = self.axis_lookup.index(x_axis)
        y_index = self.axis_lookup.index(y_axis)
        x_label = self.label_lookup[x_axis]
        y_label = self.label_lookup[y_axis]
        title = "efficiency_"+x_axis+"_vs_"+y_axis
        canvas = xboa.common.make_root_canvas(title)
        reco_hist = ROOT.TH2D("reco", "",
                                nx_bins, x_min_max[0], x_min_max[1],
                                ny_bins, y_min_max[0], y_min_max[1])
        truth_hist = ROOT.TH2D("truth", "",
                                nx_bins, x_min_max[0], x_min_max[1],
                                ny_bins, y_min_max[0], y_min_max[1])
        ratio_hist = ROOT.TH2D("ratio", ";"+x_label+";"+y_label,
                                nx_bins, x_min_max[0], x_min_max[1],
                                ny_bins, y_min_max[0], y_min_max[1])
        ratio_hist.SetStats(False)
        for sample in range(self.truth_data.n_samples):
            for a_bin in range(self.truth_data.n_bins):
                for run, spill, event, psv, amp in self.truth_data.retrieve(a_bin, sample):
                    truth_hist.Fill(psv[x_index], psv[y_index])
                for run, spill, event, psv, amp in self.reco_data.retrieve(a_bin, sample):
                    reco_hist.Fill(psv[x_index], psv[y_index])
        for i in range(truth_hist.GetNbinsX()):
            for j in range(truth_hist.GetNbinsY()):
                reco_count = reco_hist.GetBinContent(i, j)
                truth_count = truth_hist.GetBinContent(i, j)
                if truth_count > 0:
                    ratio_hist.SetBinContent(i, j, reco_count/truth_count)
        ratio_hist.Draw("COLZ")
        canvas.SetFrameFillColor(utilities.root_style.get_frame_fill())
        canvas.Update()
        for fmt in ["root", "png", "pdf"]:
            canvas.Print(self.plot_dir+"/weighting/"+title+"."+fmt)

    axis_lookup = ["x", "px", "y", "py"]
    label_lookup = {"x":"x [mm]", "y":"y [mm]", "px":"p_{x} [MeV/c]", "py":"p_{y} [MeV/c]"}
    root_objects = []
