import bisect

import numpy.random
import xboa
import xboa.common
import xboa.bunch

class ToyAmplitude(object):
    def __init__(self):
        self.mass = xboa.common.pdg_pid_to_mass[13]
        self.bin_edges = [i*5 for i in range(21)]
        self.beam = None
        self.amp_pdf = [0 for i in range(21)]
        self.cdf = [0 for i in range(21)]

    def make_beam(self, n_particles, pz_ref, emit_ref, pz_actual):
        cov = xboa.bunch.Bunch.build_penn_ellipse(emit_ref,
                                                  self.mass,
                                                  600,
                                                  0.,
                                                  pz_ref,
                                                  0.,
                                                  0.,
                                                  1.)
        cov[1, 1] = cov[1,1]*(pz_actual/pz_ref)**2
        cov[3, 3] = cov[3,3]*(pz_actual/pz_ref)**2
        print cov
        mean = numpy.array([0., 0., 0., 0.])
        self.beam = numpy.random.multivariate_normal(mean, cov, n_particles)
        print "Generated beam with shape", self.beam.shape, "first event", self.beam[0]

    def get_amplitude_pdf(self):
        cov = numpy.cov(self.beam, rowvar=False)
        emittance = numpy.linalg.det(cov)**0.25/self.mass
        print "Generated beam had emittance", emittance, "cov matrix:"
        print cov
        cov_inv = numpy.linalg.inv(cov)
        for row in self.beam:
            amplitude = emittance*numpy.dot(row.transpose(), numpy.dot(cov_inv, row))
            amp_bin = bisect.bisect_left(self.bin_edges, amplitude)
            self.amp_pdf[amp_bin-1] += 1
        self.cdf[0] = self.amp_pdf[0]
        for i, n_amp in enumerate(self.amp_pdf[1:]):
            self.cdf[i+1] = self.cdf[i]+n_amp
        print self.cdf

    @classmethod
    def make_ratio_plots(cls, toy_1, toy_2):
        for i, cdf_1 in enumerate(toy_1.cdf):
            cdf_2 = toy_2.cdf[i]
            print i, cdf_2/float(cdf_1)

def main():
    toy_140 = ToyAmplitude()
    toy_140.make_beam(1000000, 140., 10., 140.)
    toy_140.get_amplitude_pdf()
    toy_125 = ToyAmplitude()
    toy_125.make_beam(1000000, 140., 10., 125.)
    toy_125.get_amplitude_pdf()
    ToyAmplitude.make_ratio_plots(toy_140, toy_125)
    
if __name__ == "__main__":
    main()