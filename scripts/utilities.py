import ROOT

ROOT_OBJECTS = []

def fit_peak(hist):
    mean = hist.GetBinCenter(hist.GetMaximumBin())
    sigma = hist.GetRMS()
    name = hist.GetName()
    draw_color = (8, 4, 2)
    for i in range(3):
        print "FIT", i, "mean in", mean, "sigma in", sigma, "range", mean-3*sigma, mean+3*sigma
        fit = ROOT.TF1(name+" fit", "gaus")
        fit.SetLineColor(draw_color[i])
        hist.Fit(fit, "Q", "", mean-3*sigma, mean+3*sigma)
        mean = fit.GetParameter(1)
        sigma = fit.GetParameter(2)
        print "          mean out", mean, "sigma out", sigma
        
    fit = ROOT.TF1(name+"_fit", "gaus")
    fit.SetLineColor(1)
    ROOT_OBJECTS.append(fit)
    hist.Fit(fit, "Q", "", mean-3*sigma, mean+3*sigma)
    return fit

