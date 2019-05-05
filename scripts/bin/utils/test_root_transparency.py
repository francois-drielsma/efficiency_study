import random
import ROOT

canvas = ROOT.TCanvas("test transparency")
canvas.Draw()

hist1 = ROOT.TH1D("", ";x", 100, 0., 100.)
for i in range(1000):
    hist1.Fill(random.gauss(50., 20.))
hist1.SetFillColorAlpha(ROOT.kBlue, 0.5)
#hist1.SetFillColor(ROOT.kBlue)
#hist1.SetFillStyle(4050)
hist1.Draw()

hist2 = ROOT.TH1D("", ";x", 100, 0., 100.)
for i in range(1000):
    hist2.Fill(random.gauss(60., 10.))
hist2.SetFillColorAlpha(ROOT.kRed, 0.5)
#hist2.SetFillColor(ROOT.kRed)
#hist2.SetFillStyle(4050)
hist2.Draw("SAME")

graph1 = ROOT.TGraphAsymmErrors(10)
for i in range(10):
    x = i*10.+10
    y = i+0.5
    graph1.SetPoint(i, x, y)
    graph1.SetPointEXhigh(i, 3.0)
    graph1.SetPointEXlow(i, 3.0)
    graph1.SetPointEYhigh(i, 5.0)
    graph1.SetPointEYlow(i, 3.0)
    graph1.SetMarkerStyle(20)
    graph1.SetFillColorAlpha(ROOT.kGray, 0.5)
    #graph1.SetFillStyle(4100)
    #graph1.SetMarkerColorAlpha(ROOT.kGray, 0.5)
graph1.Draw("same p 3")

canvas.Print("test_transparency.png")
