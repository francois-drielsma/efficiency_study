import ROOT
import numpy
import scipy.special

def bins_equal_number(n_bins, max_bin, dim, determinant, min_bin_width, max_bin_width):
    max_bin /= determinant
    min_bin_width /= determinant
    max_bin_width /= determinant
    bins = [0.]
    while True:
        gamma = scipy.special.gammainc(dim/2., bins[-1]/2.)
        bins.append(2.*scipy.special.gammaincinv(dim/2., gamma+1./n_bins))
        if bins[-1] != bins[-1]:
            bins[-1] = max_bin
        if bins[-1] - bins[-2] < min_bin_width:
            bins[-1] = bins[-2] + min_bin_width
        elif bins[-1] - bins[-2] > max_bin_width:
            bins[-1] = bins[-2] + max_bin_width
        if bins[-1] >= max_bin:
            break
    bins = [a_bin*determinant for a_bin in bins]
    bins = numpy.array(bins, dtype='float64')
    return bins

def bins_equal_volume(n_bins, max_bin, dim, min_width, bin_volume):
    # bin volume = r^T V^{-1} r which goes proportional to |r|^2
    bin_vol = max_bin**dim/n_bins
    bins = [max_bin]
    bins = [0.]+[(bin_volume*i)**(0.5) for i in range(1, n_bins)]
    for i in range(n_bins-1):
        if bins[i+1] - bins[i] < min_width:
            break
    n_constant_bins = int((max_bin - bins[i])/min_width)
    bins = bins[:i]+[bins[i]+min_width*j for j in range(n_constant_bins)]
    bins = numpy.array(bins, dtype='float64')
    return bins

def mv_normal_data(dim, determinant, n_entries):
    mean = numpy.array([0. for i in range(dim)])
    cov = [[0. for i in range(dim)] for j in range(dim)]
    for i in range(dim):
        cov[i][i] = determinant**1./dim
    cov = numpy.array(cov)
    cov_inv = numpy.linalg.inv(cov)

    entries = numpy.random.multivariate_normal(mean, cov, n_entries)
    amplitude = [determinant*numpy.dot(numpy.dot(numpy.transpose(entry), cov_inv), entry) for entry in entries]
    return amplitude

def uniform_data(dim, determinant, n_entries):
    low = numpy.array([-determinant for i in range(dim)])
    high = numpy.array([determinant for i in range(dim)])

    cov = [[0. for i in range(dim)] for j in range(dim)]
    for i in range(dim):
        cov[i][i] = determinant**1./dim
    cov = numpy.array(cov)
    cov_inv = numpy.linalg.inv(cov)

    entries = numpy.random.uniform(-determinant, determinant, dim*n_entries)
    entries = [entries[i:i+dim] for i in range(0, n_entries, dim)]
    print entries[0]
    amplitude = [determinant*numpy.dot(numpy.dot(numpy.transpose(entry), cov_inv), entry) for entry in entries]
    return amplitude


def test_error_function():
    determinant = 6.
    n_entries = 100000
    dim = 4
    n_bins = 5
    max_bin = determinant*4.*dim
    min_width = max_bin/20.
    canvas = ROOT.TCanvas("test_canvas", "test_canvas")
    hist = ROOT.TH1D("test_h1", "test_h1", 20, 0., max_bin)
    #bins = bins_equal_volume(n_bins, max_bin, dim, 5., 300.)
    bins = bins_equal_number(n_bins, max_bin, dim, determinant, 5., 15.)
    print bins
    hist_var = ROOT.TH1D("test_h_var", "test_h_var", len(bins)-1, bins)
    hist_var.SetLineColor(ROOT.kGreen+2)
    hist_var.SetMarkerStyle(26)
    hist_var.Draw("P")
    #data = uniform_data(dim, determinant, n_entries)
    data = mv_normal_data(dim, determinant, n_entries)

    for amplitude in data:
        hist.Fill(amplitude)
        hist_var.Fill(amplitude)
    hist_stack = ROOT.THStack()
    hist_stack.Add(hist)
    hist_stack.Add(hist_var)
    hist_stack.Draw("nostack")

    #chi2 = scipy.special.gammainc(a, x)
    canvas.Update()
    raw_input()
    
if __name__ == "__main__":
    test_error_function()