import math

import numpy
import ROOT

ROOT_OBJECTS = []

def fractional_axis_range(data, fraction):
    data = sorted(data)
    data_len = len(data)
    half_fraction = (1.-fraction)/2.
    min_index = math.floor(half_fraction*data_len)
    xmin = data[int(min_index)]
    max_index = math.ceil((1.-half_fraction)*data_len)
    try:
        xmax = data[int(max_index)]
    except IndexError:
        xmax = data[-1]
    return xmin, xmax

def fit_peak(hist, nsigma=3):
    mean = hist.GetBinCenter(hist.GetMaximumBin())
    sigma = hist.GetRMS()
    name = hist.GetName()
    draw_color = (8, 4, 2)
    for i in range(3):
        fit = ROOT.TF1(name+" fit", "gaus")
        fit.SetLineColor(draw_color[i])
        hist.Fit(fit, "QN", "", mean-nsigma*sigma, mean+nsigma*sigma)
        mean = fit.GetParameter(1)
        sigma = fit.GetParameter(2)
        
    fit = ROOT.TF1(name+"_fit", "gaus")
    fit.SetLineColor(1)
    #fit.SetLineStyle(3)
    fit.SetLineWidth(2)
    fit.SetLineColor(hist.GetLineColor())
    ROOT_OBJECTS.append(fit)
    hist.Fit(fit, "Q", "", mean-nsigma*sigma, mean+nsigma*sigma)
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

def set_palette(name = None, ncontours=999):
    """Set a color palette from a given RGB list
    stops, red, green and blue should all be lists of the same length
    see set_decent_colors for an example"""

    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.34, 0.00]
    # elif name == "whatever":
        # (define more palettes)
    else:
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]

    s = numpy.array(stops)
    r = numpy.array(red)
    g = numpy.array(green)
    b = numpy.array(blue)

    npoints = len(s)
    ROOT.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    ROOT.gStyle.SetNumberContours(ncontours)



