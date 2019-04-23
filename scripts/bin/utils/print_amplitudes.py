import glob
import numpy
import xboa.common
from xboa.hit import Hit

from mice_analysis.amplitude.amplitude_data_binned import AmplitudeDataBinned
from mice_analysis.amplitude.plot_amplitude_data import PlotAmplitudeData

class AmplitudeMangler(object):
    def __init__(self):
        self.data = []
        self.canvas = xboa.common.make_root_canvas("amplitude")
        self.min_bin = 19
        self.amp_calc = None

    def load_target(self, target):
        self.data = None
        ps_matrix = None
        amp_vector = None
        for a_bin in range(21):
            for sample in range(2):
                ps_name = target+"ps_"+str(a_bin)+"_"+str(sample)
                ps = numpy.memmap(ps_name, dtype="float32", mode="r")
                if type(ps_matrix) == type(None):
                    ps_matrix = ps
                    print ps_matrix[0:4]
                else:
                    ps_matrix = numpy.concatenate((ps_matrix, ps))
                amp_name = target+"amp_"+str(a_bin)+"_"+str(sample)
                amp = numpy.memmap(amp_name, dtype="float32", mode="r")
                if type(amp_vector) == type(None):
                    amp_vector = amp
                else:
                    amp_vector = numpy.concatenate((amp_vector, amp) )
        n_rows = ps_matrix.shape[0]/4
        ps_matrix = numpy.reshape(ps_matrix, (n_rows, 4))
        amp_vector = numpy.reshape(amp_vector, (n_rows, 1))
        self.data = numpy.concatenate((ps_matrix, amp_vector), 1)
        print self.data[:3, :]

    def print_data(self, output_file_name):
        fout = open(output_file_name, "w")
        for item in self.data_keys:
            print >> fout, item.rjust(12),
        print >> fout
        for datum in self.data:
            for key in self.data_keys:
                print >>fout, format(datum[key], "12.8g"),
            print >> fout

    def get_amplitude(self):
        pass

    def plot_amplitude(self):
        self.canvas.cd()
        amp_data = self.data[:, 4:5]
        print amp_data.shape, amp_data[:5]
        hist = xboa.common.make_root_histogram("amplitude", amp_data, "A_{\perp}", 20,
                                               xmin=0., xmax=100.)
        hist.SetLineColor(4)
        hist.SetFillColorAlpha(4, 0.2)
        hist.Draw()

        psv_data = self.data[:, :4]
        cov = numpy.cov(psv_data, rowvar=False)
        cov_inv = numpy.linalg.inv(cov)
        const = numpy.linalg.det(cov)**0.25/xboa.common.pdg_pid_to_mass[13]
        amp_calc = lambda vec: const*numpy.dot(numpy.dot(numpy.transpose(vec), cov_inv), vec)
        amp_hack = [amp_calc(vec) for vec in psv_data]
        hist = xboa.common.make_root_histogram("amplitude", amp_hack, "A_{\perp}", 20,
                                               xmin=0., xmax=100.)
        hist.SetLineColor(2)
        hist.SetFillColorAlpha(2, 0.2)
        hist.Draw("SAME")

        amp_regen = self.regenerate_amplitudes()
        hist = xboa.common.make_root_histogram("amplitude", amp_regen, "A_{\perp}", 20,
                                               xmin=0., xmax=100.)
        hist.SetLineColor(8)
        hist.SetFillColorAlpha(8, 0.2)
        hist.Draw("SAME")


        self.canvas.Print("amplitude_mangler_"+str(self.min_bin)+".png")

    def regenerate_amplitudes(self):
        hit = lambda i, x: Hit.new_from_dict(
                  {'event_number':i, 'x':x[0], 'px':x[1], 'y':x[2], 'py':x[3]}
              )
        hit_list = [hit(i, item) for i, item in enumerate(self.data)]
        print len(hit_list)
        self.amp_calc = AmplitudeDataBinned('data/test_file',
                                       [i*5. for i in range(21)],
                                       xboa.common.pdg_pid_to_mass[13],
                                       False,
                                       self.min_bin,
                                       2000,
                                       )
        self.amp_calc.clear()
        self.amp_calc.append_hits(hit_list)
        amp_dict = self.amp_calc.fractional_amplitude()
        return amp_dict.values()

    def plot_data(self):
        plotter = PlotAmplitudeData(self.amp_calc, "./", "test")
        plotter.plot()

    data_keys = ["x", "px", "y", "py", "amp_4d"]
    root_objects = []


def main():
    a_dir = "../output/2017-02-7-v7/plots_2017-2.7_10-140_LiH/amplitude/data/"
    target = a_dir+"amp_data_recon_us_"
    mangle = AmplitudeMangler()
    mangle.load_target(target)
    #mangle.print_data("10-140_LiH_ps_us")
    for min_bin in [0]:#, 10, 6, 4, 2, 0]:
        mangle.min_bin = min_bin
        mangle.plot_amplitude()
    mangle.plot_data()

if __name__ == "__main__":
    main()