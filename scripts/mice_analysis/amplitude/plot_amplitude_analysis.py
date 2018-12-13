import sys

import ROOT
import xboa.common as common

import utilities
from mice_analysis.amplitude.weighting import Weighting


class PlotAmplitudeAnalysis(object):
    def __init__(self, amplitude_analysis):
        self.amp = amplitude_analysis
        self.config = self.amp.config
        self.config_anal = self.amp.config_anal
        self.plot_dir = self.amp.plot_dir
        self.us_color = ROOT.kOrange+4
        self.ds_color = ROOT.kGreen+3


    def amplitude_scatter_plot(self, suffix):
        """
        Make a scatter plot of upstream versus downstream data
        """
        amp_dict_ds = self.amp.amplitudes[suffix]["amplitude_dict_downstream"]
        amp_dict_us = self.amp.amplitudes[suffix]["amplitude_dict_upstream"]
        amp_list_us = []
        amp_list_ds = []
        amp_list_delta = []

        suffix_label = self.get_suffix_label(suffix)
        for key, amp_ds in amp_dict_ds.iteritems():
            try:
                amp_us = amp_dict_us[key]
            except KeyError:
                sys.excepthook(*sys.exc_info())
                continue
            amp_delta = amp_us-amp_ds
            amp_list_us.append(amp_us)
            amp_list_ds.append(amp_ds)
            amp_list_delta.append(amp_delta)
        canvas = common.make_root_canvas("amplitude_residuals")
        canvas.Draw()
        n_points = min(len(amp_list_us), 10000) # no more than 10k points in the scatter
        hist, graph = common.make_root_graph("delta amplitude scatter",
                                             amp_list_us, "US "+suffix_label+" Amplitude [mm]",
                                             amp_list_delta, suffix_label+" US Amplitude - DS Amplitude [mm]", xmin=0., xmax=100., ymin=-50., ymax=50.)
        hist.SetTitle(self.config_anal['name'])
        hist.Draw()
        graph.Draw("P")
        canvas.Update()
        for format in ["eps", "png", "root"]:
            canvas.Print(self.plot_dir+"amplitude_delta_"+suffix+"_scatter."+format)

        canvas = common.make_root_canvas("delta_amplitude_hist")
        canvas.Draw()
        canvas.SetFrameFillColor(utilities.utilities.get_frame_fill())
        hist = common.make_root_histogram("delta amplitude hist",
                                          amp_list_us, suffix_label+" US Amplitude [mm]", 100,
                                          amp_list_delta, suffix_label+" US Amplitude - DS Amplitude [mm]", 100,
                                          xmin=0., xmax=100., ymin=-50., ymax=50.)
        hist.SetTitle(self.config_anal['name'])
        hist.Draw("COLZ")
        canvas.Update()
        for format in ["eps", "png", "root"]:
            canvas.Print(self.plot_dir+"amplitude_delta_"+suffix+"_hist."+format)


    def amplitude_residuals_plot(self, us_ds):
        """
        Make a plot of MC true - MC reco vs MC true
        """
        dict_name = "amplitude_dict_"+us_ds
        amp_true_dict = self.amp.amplitudes["reco_mc"][dict_name]
        amp_reco_dict = self.amp.amplitudes["reco"][dict_name]
        nx_bins = len(self.amp.bin_edge_list)-1
        xmax = max(self.amp.bin_edge_list)
        name = "amplitude residuals "+us_ds
        tk = {"upstream":"TKU", "downstream":"TKD"}[us_ds]
        title = self.config_anal["name"]
        ny_bins = 100
        res_hist = ROOT.TH2D(name, 
                          title+";"+tk+" A(true);"+tk+" A(reco) - A(true) [mm]",
                          nx_bins, 0., xmax,
                          ny_bins, -20., 20.)
        self.root_objects.append(res_hist)

        for key, amp_true in amp_true_dict.iteritems():
            amp_reco = amp_reco_dict[key]
            res_hist.Fill(amp_true, amp_reco-amp_true)
        # normalise the columns
        for i in range(1, res_hist.GetXaxis().GetNbins()+1):
            bin_sum = 0
            for j in range(res_hist.GetYaxis().GetNbins()+2):
                bin_sum += res_hist.GetBinContent(i, j)
            if abs(bin_sum) < 1e-6:
                continue
            for j in range(res_hist.GetYaxis().GetNbins()+2):
                bin_content = res_hist.GetBinContent(i, j)
                res_hist.SetBinContent(i, j, bin_content/bin_sum)

        canvas = common.make_root_canvas(name)
        #canvas.SetFrameFillColor(utilities.utilities.get_frame_fill())
        res_hist.SetStats(False)
        res_hist.Draw("COLZ")
        name = name.replace(" ", "_")
        for format in ["eps", "png", "root"]:
            canvas.Print(self.plot_dir+name+"."+format)

    def cdf_plot(self, suffix):
        us_data = self.amp.amplitudes[suffix]["all_upstream"]
        ds_data = self.amp.amplitudes[suffix]["all_downstream"]
        name = "amplitude_cdf_"+suffix
        canvas = common.make_root_canvas(name)
        canvas.SetName(name)

        upstream_graph = self.get_asymm_error_graph(us_data["corrected_cdf"],
                                                    us_data["cdf_stats_errors"],
                                                    style=24, color=self.us_color,
                                                    name = "Upstream CDF stats")
        downstream_graph = self.get_asymm_error_graph(ds_data["corrected_cdf"],
                                                      ds_data["cdf_stats_errors"],
                                                    style=26, color=self.ds_color,
                                                    name = "Downstream CDF stats")
        upstream_graph_sys = self.get_asymm_error_graph(us_data["corrected_cdf"],
                                                    us_data["cdf_sys_errors"],
                                                    fill=self.us_color,
                                                    name = "Upstream CDF sys")
        downstream_graph_sys = self.get_asymm_error_graph(ds_data["corrected_cdf"],
                                                      ds_data["cdf_sys_errors"],
                                                    fill=self.ds_color, name = "Downstream CDF sys")
        graph_list = [upstream_graph, downstream_graph, upstream_graph_sys, downstream_graph_sys]
        draw_list = ["p", "p", "2", "2"]

        hist = self.get_hist(graph_list, self.get_suffix_label(suffix)+" Amplitude [mm]", "Cumulative Number")
        hist.Draw()
        for i, graph in enumerate(graph_list):
            graph.Draw("SAME "+draw_list[i])
        self.text_box(graph_list)

        canvas.Update()
        for a_format in ["eps", "pdf", "root", "png"]:
            canvas.Print(self.plot_dir+"amplitude_cdf_"+suffix+"."+a_format)
        return canvas

    def cdf_ratio_plot(self, suffix):
        data = self.amp.amplitudes[suffix]["ratio"]
        name = "cdf_ratio_"+suffix
        canvas = common.make_root_canvas(name)
        canvas.SetName(name)

        cdf_graph_stats = self.get_asymm_error_graph(data["corrected_cdf"],
                                               data["cdf_stats_errors"],
                                               style=20, name = "CDF Ratio stats")
        cdf_graph_sys = self.get_asymm_error_graph(data["corrected_cdf"],
                                               data["cdf_sys_errors"], fill=ROOT.kGray, 
                                               name = "CDF Ratio sys")
        print "CDF"
        print data["corrected_cdf"]
        print data["cdf_stats_errors"]
        print data["cdf_sys_errors"]
        graph_list = [cdf_graph_stats, cdf_graph_sys]
        draw_list = ["p", "2"]
        
        hist = self.get_hist(graph_list, self.get_suffix_label(suffix)+" Amplitude [mm]", "Cumulative Number Ratio")
        hist.Draw()
        for i, graph in enumerate(graph_list):
            graph.Draw("SAME "+draw_list[i])
        self.text_box(graph_list)

        canvas.Update()
        for a_format in ["eps", "pdf", "root", "png"]:
            canvas.Print(self.plot_dir+canvas.GetName()+"."+a_format)
        return canvas

    def pdf_ratio_plot(self, suffix):
        data = self.amp.amplitudes[suffix]["ratio"]
        name = "pdf_ratio_"+suffix
        canvas = common.make_root_canvas(name)
        canvas.SetName(name)

        pdf_graph_stats = self.get_asymm_error_graph(data["corrected_pdf"],
                                               data["pdf_stats_errors"],
                                               style=20, name = "PDF Ratio stats")
        pdf_graph_sys = self.get_asymm_error_graph(data["corrected_pdf"],
                                               data["pdf_sys_errors"], fill=ROOT.kGray,
                                               name = "PDF Ratio sys")
        graph_list = [pdf_graph_stats, pdf_graph_sys]
        draw_list = ["p", "2"]

        hist = self.get_hist(graph_list, self.get_suffix_label(suffix)+" Amplitude [mm]", "Number Ratio")
        hist.Draw()
        for i, graph in enumerate(graph_list):
            graph.Draw("SAME "+draw_list[i])
        self.text_box(graph_list)
        print "PDF"
        print data["corrected_pdf"]
        print data["pdf_stats_errors"]
        print data["pdf_sys_errors"]

        canvas.Update()
        for a_format in ["eps", "pdf", "root", "png"]:
            canvas.Print(self.plot_dir+canvas.GetName()+"."+a_format)
        return canvas


    def pdf_plot(self, suffix):
        data = self.amp.amplitudes[suffix]
        us_data = self.amp.amplitudes[suffix]["all_upstream"]
        ds_data = self.amp.amplitudes[suffix]["all_downstream"]
        name = "amplitude_pdf_"+suffix
        canvas = common.make_root_canvas(name)
        canvas.SetName(name)
        raw_upstream_graph = self.get_asymm_error_graph(us_data["pdf"],
                                                    style=24, color=self.us_color, name="Raw upstream")
        raw_downstream_graph = self.get_asymm_error_graph(ds_data["pdf"],
                                                    style=26, color=self.ds_color, name="Raw downstream")
        scraped_graph = self.get_asymm_error_graph(data["upstream_scraped"]["pdf"],
                                                    style=25, color=ROOT.kViolet+2, name="Raw scraped") # 25 for not raw
        
        upstream_graph_stats = self.get_asymm_error_graph(us_data["corrected_pdf"],
                                                    us_data["pdf_stats_errors"],
                                                    style=20, color=self.us_color,
                                                    name = "Upstream stats")
        upstream_graph_sys = self.get_asymm_error_graph(us_data["corrected_pdf"],
                                                    us_data["pdf_sys_errors"],
                                                    fill=self.us_color,
                                                    name = "Upstream sys")
        downstream_graph_stats = self.get_asymm_error_graph(ds_data["corrected_pdf"],
                                                      ds_data["pdf_stats_errors"],
                                                    style=22, color=self.ds_color,
                                                    name = "Downstream stats")
        downstream_graph_sys = self.get_asymm_error_graph(ds_data["corrected_pdf"],
                                                      ds_data["pdf_sys_errors"],
                                                    fill=self.ds_color,
                                                    name = "Downstream sys")
        print "Plotting us sys", us_data["pdf_sys_errors"]
        print "Plotting ds sys", ds_data["pdf_sys_errors"]
        graph_list = [upstream_graph_sys, downstream_graph_sys, upstream_graph_stats, downstream_graph_stats, scraped_graph, raw_upstream_graph, raw_downstream_graph]
        draw_list = ["2", "2", "p", "p", "p", "p", "p"]
        hist = self.get_hist(graph_list, self.get_suffix_label(suffix)+" Amplitude [mm]", "Number")
        hist.Draw()
        same = "SAME "
        for i, graph in enumerate(graph_list):
            graph.Draw(same+draw_list[i])
        self.text_box(graph_list)

        if self.config_anal["amplitude_chi2"]:
            upstream_chi2 = self.chi2_graph(suffix, "all_upstream")
            upstream_chi2.SetLineColor(self.us_color) 
            upstream_chi2.Draw("SAMEL")
            downstream_chi2 = self.chi2_graph(suffix, "all_downstream")
            downstream_chi2.SetLineColor(self.ds_color)
            downstream_chi2.Draw("SAMEL")

        canvas.Update()
        for a_format in ["eps", "pdf", "root", "png"]:
            canvas.Print(self.plot_dir+"amplitude_pdf_"+suffix+"."+a_format)
        return canvas

    def get_hist(self, graph_list, x_axis, y_axis, min_x = None, max_x = None, min_y = None, max_y = None):
        if min_x == None:
            min_x = min([graph.GetXaxis().GetXmin() for graph in graph_list])
        if max_x == None:
            max_x = max([graph.GetXaxis().GetXmax() for graph in graph_list])
        if min_y == None:
            min_y = min([graph.GetYaxis().GetXmin() for graph in graph_list])
        if max_y == None:
            max_y = max([graph.GetYaxis().GetXmax() for graph in graph_list])
        print "HIST", min_x, max_x, min_y, max_y
        hist = ROOT.TH2D(graph_list[0].GetName()+"_hist", ";"+x_axis+";"+y_axis,
                          1000, min_x, max_x, 1000, min_y, max_y)
        hist.SetStats(False)
        self.root_objects.append(hist)
        return hist

    def matrix_plot(self, matrix, title, axis_1, axis_2):
        """
        Plot the matrix entries as a TH2D
        
        * matrix: matrix data to plot; list n_bins by n_bins
        * title: the title of the hist (used in the file name also)
        * axis_1: x axis title
        * axis_2: y axis title
        """
        bin_centre_list = self.amp.amplitudes["reco"]["all_upstream"]["bin_centre_list"]
        xmax = max(self.amp.bin_edge_list)
        nbins = len(bin_centre_list)
        matrix_hist = ROOT.TH2D(title, title+";"+axis_1+";"+axis_2, nbins, 0., xmax, nbins, 0., xmax)
        for i in range(nbins):
            for j in range(nbins):
                matrix_hist.Fill(bin_centre_list[i], bin_centre_list[j], matrix[i][j])
        self.root_objects.append(matrix_hist)
        canvas = common.make_root_canvas(title)
        canvas.SetFrameFillColor(utilities.utilities.get_frame_fill())
        matrix_hist.SetTitle(self.config_anal['name'])
        matrix_hist.SetStats(False)
        matrix_hist.Draw("COLZ")
        title = title.replace(" ", "_")
        for format in ["eps", "png", "root"]:
            canvas.Print(self.plot_dir+title+"."+format)

    def efficiency_plot(self):
        weighting = Weighting(self.amp.reco_mc_data_ds, self.amp.all_mc_data_ds, self.plot_dir)
        weighting.plot_sum("x", 20, [-200., 200.], "y", 20, [-200., 200.])
        weighting.plot_sum("x", 20, [-200., 200.], "px", 20, [-100., 100.])
        weighting.plot_sum("x", 20, [-200., 200.], "py", 20, [-100., 100.])
        weighting.plot_sum("y", 20, [-200., 200.], "px", 20, [-100., 100.])
        weighting.plot_sum("y", 20, [-200., 200.], "py", 20, [-100., 100.])
        weighting.plot_sum("px", 20, [-100., 100.], "py", 20, [-100., 100.])

    def systematics_plot(self, target):
        """
        Plot systematic corrections:
        * Momentum correction matrix
        * Crossing probability due to reconstruction
        * Inefficiency/impurity
        """
        if "upstream" in target:
            label = "US"
        elif "downstream" in target:
            label = "DS"
        matrix = self.amp.amplitudes["crossing_probability"][target]["migration_matrix"]
        title = ("crossing_probability_"+target).replace("_all", "").replace("_", " ")
        self.matrix_plot(matrix, title,
                         label+" "+self.get_suffix_label("reco")+" amplitude [mm]",
                         label+" "+self.get_suffix_label("reco_mc")+" amplitude [mm]")

        inefficiency = self.amp.amplitudes["inefficiency"][target]["pdf_ratio"]
        canvas = common.make_root_canvas("inefficiency_"+target)

        hist, graph = common.make_root_graph("inefficiency",
                                             self.amp.get_bin_centre_list(), label+" MC Truth Amplitude [mm]",
                                             inefficiency[:-1], label+" "+self.get_suffix_label("all_mc")+" number/"+self.get_suffix_label("reco_mc")+" number")
        hist.GetYaxis().SetRangeUser(0., 2.0)
        hist.SetTitle(self.config_anal['name'])
        hist.Draw()
        graph.SetMarkerStyle(24)
        graph.Draw("p")
        for format in "eps", "root", "png":
            canvas.Print(self.plot_dir+"amplitude_inefficiency_"+target+"."+format)
        hist.GetYaxis().SetRangeUser(0.9, 1.1)
        hist.Draw()
        graph.SetMarkerStyle(24)
        graph.Draw("p")
        for format in "eps", "root", "png":
            canvas.Print(self.plot_dir+"amplitude_inefficiency_zoom_"+target+"."+format)

    def chi2_graph(self, suffix, distribution):
        data = self.amp.amplitudes[suffix]
        pdf_list_tku = data[distribution]["pdf"]
        bin_centre_list = data[distribution]["bin_centre_list"]
        n_bins = len(bin_centre_list)
        weighted_integral = [bin_centre_list[i]*pdf_list_tku[i] for i in range(n_bins)]
        integral = sum(pdf_list_tku)#*(bin_centre_list[1]-bin_centre_list[0])
        max_bin = max(pdf_list_tku)
        emittance = sum(weighted_integral)/integral/4.
        dummy, graph = utilities.chi2_distribution.chi2_graph(emittance, max_bin, 100, 0., 100.)
        return graph

    def get_asymm_error_graph(self, points, errors=None, norm=1., style=None, color=None, fill=None, name="Graph"):
        graph = ROOT.TGraphAsymmErrors(len(points)-1)
        for i, low_edge in enumerate(self.amp.bin_edge_list[:-1]):
            high_edge = self.amp.bin_edge_list[i+1]
            centre = (low_edge+high_edge)/2.
            graph.SetPoint(i, centre, points[i]/norm)
            if errors != None:
                graph.SetPointError(i, centre-low_edge, high_edge-centre, errors[i]/norm, errors[i]/norm)
        if style != None:
            graph.SetMarkerStyle(style)
        if color != None:
            graph.SetMarkerColor(color)
        if fill != None:
            graph.SetFillColor(fill)
            graph.SetFillStyle(3001);
        graph.SetName(name)
        self.root_objects.append(graph)
        return graph

    def text_box(self, graph_list):
        legend = ROOT.TLegend(0.6, 0.65, 0.89, 0.89)
        for graph in graph_list:
            legend.AddEntry(graph, graph.GetName(), "lep")
        legend.SetBorderSize(0)
        legend.Draw()
        self.root_objects.append(legend)
        text_box = ROOT.TPaveText(0.6, 0.55, 0.89, 0.65, "NDC")
        text_box.SetFillColor(0)
        text_box.SetBorderSize(0)
        text_box.SetTextSize(0.04)
        text_box.SetTextAlign(12)
        text_box.AddText("MICE INTERNAL")
        text_box.Draw()
        self.root_objects.append(text_box)
        text_box = ROOT.TPaveText(0.6, 0.45, 0.89, 0.55, "NDC")
        text_box.SetFillColor(0)
        text_box.SetBorderSize(0)
        text_box.SetTextSize(0.03)
        text_box.SetTextAlign(12)
        text_box.AddText(self.config_anal['name'])
        text_box.Draw()
        self.root_objects.append(text_box)

    def get_suffix_label(self, suffix):
        return {"all_mc":"MC truth (all)",
                "reco_mc":"MC truth (recon)",
                "reco":"Reconstructed"}[suffix]

    root_objects = []