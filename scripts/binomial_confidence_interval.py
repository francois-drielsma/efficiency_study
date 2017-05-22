import xboa.common
import numpy
import ROOT

def test_binomial_confidence_interval(number_of_events, probability):
    print "Running test for n:", number_of_events, "p:", probability
    cx_list = range(number_of_events+1)
    c_list = [ROOT.TMath.BinomialI(probability, number_of_events, x) for x in cx_list]
    x_list = [i+0.5 for i in range(number_of_events)]
    f_list = [(c_list[i] - c_list[i+1])/2. for i in range(number_of_events)]

    canvas = xboa.common.make_root_canvas("binomial")
    hist, graph = xboa.common.make_root_graph("binomial", x_list, "x", f_list, "f(x)")
    hist.SetTitle("Binomial with n = "+str(number_of_events)+" and p = "+str(probability))
    canvas.Draw()
    hist.Draw()
    graph.Draw("SAME")
    canvas.Update()

    for interval, color in [(0.682, ROOT.kGray+1), (0.95, ROOT.kGray+2), (0.99, ROOT.kGray+3)]:
        lower, upper = BinomialConfidenceInterval.binomial_confidence_interval(int(probability*number_of_events), number_of_events, interval)
        print "    ", interval, "CI - lower", lower, "p", probability*number_of_events, "upper", upper
        assert(lower < probability*number_of_events)
        assert(upper > probability*number_of_events)
        hist, graph = xboa.common.make_root_graph("binomial", [lower, lower], "x", [-1., 10.], "f(x)")
        graph.Draw("SAME")
        graph.SetLineColor(color)
        hist, graph = xboa.common.make_root_graph("binomial", [upper, upper], "x", [-1., 10.], "f(x)")
        graph.Draw("SAME")
        graph.SetLineColor(color)
    canvas.Update()

class BinomialConfidenceInterval(object):
    def __init__(self, n_in_bin, n_total, interval):
        self.probability = n_in_bin/float(n_total)
        self.number_of_events = n_total
        self.y0 = ROOT.TMath.BinomialI(self.probability, n_total, n_in_bin)
        self.lower_bound = self.binary_search((1.-interval)/2., self.lower_dir)
        self.upper_bound = self.binary_search(1.-(1.-interval)/2., self.upper_dir)

    def binary_search(self, target_value, direction):
        # positive
        x0 = 0
        x1 = self.number_of_events/2
        x2 = self.number_of_events
        y0 = 0
        y2 = 1.
        while x0 != x1 and x2 != x1:
            y1 = 1-ROOT.TMath.BinomialI(self.probability, self.number_of_events, x1)
            if y1 > target_value:
                x2 = x1
                y2 = y1
            else:
                x0 = x1
                y0 = y1
            x1 = (x2 + x0)/2
            #print x0, y0, "**", x1, y1, "**", x2, y2, "**", target_value
        if direction == self.upper_dir:
            return x2
        else:
            return x0

    @classmethod
    def binomial_confidence_interval(cls, n_in_bin, n_total, interval):
        ci = BinomialConfidenceInterval(n_in_bin, n_total, interval)
        return ci.lower_bound, ci.upper_bound

    upper_dir = 1
    lower_dir = -1


if __name__ == "__main__":
    test_binomial_confidence_interval(819, 129./819.)
    test_binomial_confidence_interval(50, 0.1)
    test_binomial_confidence_interval(50, 0.9)
    test_binomial_confidence_interval(5000, 0.1)
    raw_input()
