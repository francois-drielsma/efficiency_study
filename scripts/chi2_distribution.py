import numpy
import scipy.stats
import xboa.common

def chi2_graph(emittance, max_bin, n_points, a_min, a_max):
    n_dof = 4
    if emittance == 0.:
        emittance = 1.
    bin_width = (a_max-a_min)/n_points
    x_list = [float(i)*bin_width+a_min for i in range(n_points)]
    chi2_list = [scipy.stats.chi2.pdf(x/(emittance), n_dof) for x in x_list]
    normalisation = max_bin/max(chi2_list)
    chi2_list = [chi2*normalisation for chi2 in chi2_list]
    hist, graph = xboa.common.make_root_graph("chi2", x_list, "\Chi^{2}", chi2_list, "")
    # mean = sum(x f(x))/sum(f(x))
    bin_list = [x_list[i]*chi2_list[i] for i in range(n_points)]
    print "Distribution mean", sum(bin_list)/sum(chi2_list)
    return hist, graph

def chi2_plot(emittance, sum_of_bins, n_points, a_min, a_max):
    hist, graph = chi2_graph(emittance, sum_of_bins, n_points, a_min, a_max)
    canvas = xboa.common.make_root_canvas("chi2")
    canvas.Draw()
    hist.Draw()
    graph.SetMarkerStyle(24)
    graph.Draw("SAMEP")
    canvas.Update()
    return canvas, hist, graph

if __name__ == "__main__":
    chi2_plot(6, 1000, 10, 0., 100.)
    raw_input()

