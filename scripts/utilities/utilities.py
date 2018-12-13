import math

import numpy
import ROOT

import xboa.common

ROOT_OBJECTS = []

def default_units(var):
    units_dict = {"x":"mm", "y":"mm", "z":"mm", "r":"mm",
                  "SP Res(x)":"mm", "SP Res(y)":"mm", 
                  "px":"MeV/c", "py":"MeV/c", "pz":"MeV/c", "pt":"MeV/c", "p":"MeV/c",
                  "energy":"MeV", "kinetic_energy":"MeV", "e_dep":"MeV",
                  "time":"ns", "delta_tof":"ns",
                  "tof01":"ns", "tof12":"ns", "tof02":"ns"}
    if var not in units_dict.keys():
        return ""
    return units_dict[var]

def detector_position(detector, config):
    for test_detector in config.detectors:
        if test_detector[2] == detector:
            return test_detector[0]
    raise KeyError("Could not find detector "+str(detector)+" in config")

def electron_tof(detector_1, detector_2, config):
    z_1 = detector_position(detector_1, config)
    z_2 = detector_position(detector_2, config)
    c_light = xboa.common.constants['c_light']
    electron_tof = (z_2-z_1)/c_light
    return electron_tof

def fractional_axis_range(data, fraction):
    data = sorted(data)
    data_len = len(data)
    if data_len == 0:
        return 0, 0
    half_fraction = (1.-fraction)/2.
    min_index = math.floor(half_fraction*data_len)
    max_index = math.ceil((1.-half_fraction)*data_len)
    try:
        xmin = data[int(min_index)]
        xmax = data[int(max_index)]
    except IndexError:
        xmin = data[0]
        xmax = data[-1]
    return xmin, xmax

def fit_peak_data(data, nsigma=3, frac_guess=0.99):
    xmin, xmax = fractional_axis_range(data, frac_guess)
    canvas = xboa.common.make_root_canvas("dummy")
    hist = xboa.common.make_root_histogram("dummy", data, "", 100, xmin=xmin, xmax=xmax)
    fit = fit_peak(hist, nsigma)
    return fit

def fit_peak(hist, nsigma=3, fit_option="Q", draw_option="", fit_range=None):
    mean = hist.GetBinCenter(hist.GetMaximumBin())
    sigma = hist.GetRMS()
    name = hist.GetName()
    draw_color = (8, 4, 2)
    # try to guess a good fit range
    for i in range(3):
        if fit_range != None:
            continue
        fit = ROOT.TF1(name+" fit", "gaus")
        fit.SetLineColor(draw_color[i])
        hist.Fit(fit, "QN", "", mean-nsigma*sigma, mean+nsigma*sigma)
        mean = fit.GetParameter(1)
        sigma = fit.GetParameter(2)
        
    fit = ROOT.TF1(name+"_fit", "gaus")
    fit.SetLineWidth(2)
    fit.SetLineColor(hist.GetLineColor())
    ROOT_OBJECTS.append(fit)
    if fit_range == None:
        hist.Fit(fit, fit_option, draw_option, mean-nsigma*sigma, mean+nsigma*sigma)
    else:
        hist.Fit(fit, fit_option, draw_option, fit_range[0], fit_range[1])
    return fit

text_boxes = []
def get_text_box(config, config_anal, data_list = None, fit = None, hist = None):
    text_box = ROOT.TPaveText(0.6, 0.4, 0.9, 0.9, "NDC")
    text_box.SetFillColor(0)
    text_box.SetBorderSize(0)
    text_box.SetTextSize(0.04)
    text_box.SetTextAlign(12)
    text_box.SetTextSize(0.03)
    mean = None
    sigma = None
    n_events = None
    if hist != None:
        mean = round(hist.GetMean())
        sigma = round(hist.GetRMS())
        n_events = hist.GetEntries()
    if data_list != None:
        mean = round(numpy.mean(data_list), 2)
        sigma = round(numpy.std(data_list), 2)
        n_events = len(data_list)
    if fit != None:
        mean = round(fit.GetParameter(1), 2)
        sigma = round(fit.GetParameter(2), 2)
    text_box.AddText("MICE Internal")
    text_box.AddText(config_anal["name"])
    text_box.AddText("Recon: "+config.maus_version)
    text_box.AddText("All events (black)")
    text_box.AddText("Number: "+str(n_events))
    text_box.AddText("Mean:   "+str(mean))
    text_box.AddText("Std:    "+str(sigma))
    text_box.SetBorderSize(1)
    text_box.Draw()
    text_boxes.append(text_box)
    return text_box

def set_root_verbosity(verbose_level):
    #verb_map = [ROOT.kPrint, ROOT.kInfo, ROOT.kWarning, ROOT.kError, ROOT.kBreak, ROOT.kSysError, ROOT.kFatal]
    ROOT.gErrorIgnoreLevel = 6000

def setup_gstyle():
    stops = [0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000]
    red   = [0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764]
    green = [0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832]
    blue  = [0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539]
    s = numpy.array(stops)
    r = numpy.array(red)
    g = numpy.array(green)
    b = numpy.array(blue)

    ncontours = 255
    npoints = len(s)
    ROOT.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    ROOT.gStyle.SetNumberContours(ncontours)

    # axes and labels
    ROOT.gStyle.SetPadBottomMargin(0.15)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    ROOT.gStyle.SetPadRightMargin(0.15)
    for axis in "X", "Y":
        ROOT.gStyle.SetNdivisions(505, axis)
        ROOT.gStyle.SetLabelSize(0.05, axis)
        ROOT.gStyle.SetTitleSize(0.06, axis)
        ROOT.gStyle.SetTitleOffset(1.10, axis)

def get_frame_fill():
    level = 0.9
    return ROOT.TColor.GetColor(0.2082*level, 0.1664*level, 0.5293*level)

import random
root_objects = []

def test_hist_2d():
    canvas = ROOT.TCanvas("test 2d", "test 2d")
    hist = ROOT.TH2D("test 2d", ";x axis [units];y axis [units]", 50, -1, 1, 50, -1, 1)
    hist.SetStats(False)
    for i in range(100000):
        x, y = random.gauss(0, 0.3), random.gauss(0, 0.3)
        hist.Fill(x, y)
    canvas.Draw()
    # if you want to emphasise the tails, leave the frame white
    hist.Draw("COLZ")
    canvas.Update()
    canvas.Print("test_hist_2d_no_fill.eps")
    canvas.Print("test_hist_2d_no_fill.png")

    # if you want to de-emphasise the tails, use a frame fill
    canvas.SetFrameFillColor(get_frame_fill())

    do_axes(hist)
    hist.Draw("COLZ")
    canvas.Update()
    canvas.Print("test_hist_2d_fill.eps")
    canvas.Print("test_hist_2d_fill.png")
    # these two lines force python to keep the ROOT stuff in memory
    root_objects.append(canvas)
    root_objects.append(hist)

def test_hist_1d():
    canvas = ROOT.TCanvas("test 1d", "test 1d")
    hist_data = ROOT.TH1D("test 1d data", ";x axis [units];y axis [units]", 100, -1, 1)
    hist_mc = ROOT.TH1D("test 1d mc", ";x axis [units];y axis [units]", 100, -1, 1)
    hist_data.SetStats(False)
    hist_mc.SetStats(False)
    for i in range(100000):
        x_data = random.gauss(0, 0.3)
        hist_data.Fill(x_data)
        x_mc = random.gauss(0.01, 0.27)
        hist_mc.Fill(x_mc)
    canvas.Draw()
    hist_mc.SetFillColor(ROOT.kOrange-2)
    hist_mc.Draw()

    hist_data.SetMarkerStyle(20)
    hist_data.Draw("e1 p same")
    canvas.Update()
    canvas.Print("test_hist_1d.eps")
    canvas.Print("test_hist_1d.png")
    # these two lines force python to keep the ROOT stuff in memory
    root_objects.append(canvas)
    root_objects.append(hist_data)
    root_objects.append(hist_mc)

    # if you are comparing MC graph with data graph, then following style is recommended:
    hist_mc.SetMarkerStyle(26)
    hist_mc.SetMarkerColor(ROOT.kRed)
    hist_mc.SetLineColor(ROOT.kRed)
    hist_mc.Draw("e1 p")

    hist_data.SetMarkerStyle(20)
    hist_data.Draw("e1 p same")
    canvas.Update()
    canvas.Print("test_graph_1d.eps")
    canvas.Print("test_graph_1d.png")

if __name__ == "__main__":
    setup_gstyle()
    test_hist_2d()
    test_hist_1d()
    raw_input("Finished")

