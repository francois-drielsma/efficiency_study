import random
import numpy
import ROOT
import os

ROOT_OBJECTS = []

def set_root_verbosity(verbose_level = 6000):
    ROOT.gErrorIgnoreLevel = verbose_level

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
    ROOT.gStyle.SetHistLineWidth(1)
    for axis in "X", "Y":
        ROOT.gStyle.SetNdivisions(505, axis)
        ROOT.gStyle.SetLabelSize(0.04, axis)
        ROOT.gStyle.SetTitleSize(0.06, axis)
        ROOT.gStyle.SetTitleOffset(1.10, axis)

def get_frame_fill():
    level = 0.9
    return ROOT.TColor.GetColor(0.2082*level, 0.1664*level, 0.5293*level)

def text_box():
    text_box = ROOT.TPaveText(0.6, 0.8, 0.85, 0.899, "NDC")
    text_box.SetBorderSize(0)
    text_box.SetFillColor(0)
    text_box.SetTextSize(0.04)
    text_box.SetTextAlign(12)
    text_box.AddText("MICE preliminary")
    text_box.Draw()
    ROOT_OBJECTS.append(text_box)
    text_box = ROOT.TPaveText(0.6, 0.7, 0.85, 0.8, "NDC")
    text_box.SetBorderSize(0)
    text_box.SetFillColor(0)
    text_box.SetTextSize(0.03)
    text_box.SetTextAlign(12)
    text_box.AddText("ISIS Cycle 2017/03")
    text_box.AddText("Run 7469 MAUS v1.2.3")
    text_box.Draw()
    ROOT_OBJECTS.append(text_box)

def test_hist_2d():
    canvas = ROOT.TCanvas("test 2d", "test 2d")
    hist = ROOT.TH2D("test 2d", ";x axis [units];y axis [units]", 50, -1, 1, 70, -1, 1.8)
    hist.SetStats(False)
    for i in range(100000):
        x, y = random.gauss(0, 0.3), random.gauss(0, 0.3)
        hist.Fill(x, y)
    canvas.Draw()
    # if you want to emphasise the tails, leave the frame white
    hist.Draw("COLZ")
    text_box()
    canvas.Update()
    canvas.Print("test_hist_2d_no_fill.eps")
    canvas.Print("test_hist_2d_no_fill.png")

    # if you want to de-emphasise the tails, use a frame fill
    canvas.SetFrameFillColor(get_frame_fill())

    hist.Draw("COLZ")
    text_box()
    canvas.Update()
    canvas.Print("test_hist_2d_fill.eps")
    canvas.Print("test_hist_2d_fill.png")
    # these two lines force python to keep the ROOT stuff in memory
    ROOT_OBJECTS.append(canvas)
    ROOT_OBJECTS.append(hist)

def test_hist_1d():
    canvas = ROOT.TCanvas("test 1d", "test 1d")
    hist_data = ROOT.TH1D("data", ";x axis [units];y axis [units]", 100, -1, 1)
    hist_mc = ROOT.TH1D("simulation", ";x axis [units];y axis [units]", 100, -1, 1)
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
    text_box()
    canvas.Update()
    canvas.Print("test_hist_1d.eps")
    canvas.Print("test_hist_1d.png")
    # these two lines force python to keep the ROOT stuff in memory
    ROOT_OBJECTS.append(canvas)
    ROOT_OBJECTS.append(hist_data)
    ROOT_OBJECTS.append(hist_mc)

    # if you are comparing MC graph with data graph, then following style is recommended:
    hist_mc.SetMarkerStyle(26)
    hist_mc.SetMarkerColor(ROOT.kRed)
    hist_mc.SetLineColor(ROOT.kRed)
    hist_mc.Draw("e1 p")

    hist_data.SetMarkerStyle(20)
    hist_data.Draw("e1 p same")
    text_box()
    canvas.Update()
    canvas.Print("test_graph_1d.eps")
    canvas.Print("test_graph_1d.png")

if __name__ == "__main__":
    test_dir = "style_test"
    try:
        os.mkdir(test_dir)
    except OSError:
        pass
    os.chdir(test_dir)
    setup_gstyle()
    test_hist_2d()
    test_hist_1d()
    raw_input("Finished - output in "+test_dir)

