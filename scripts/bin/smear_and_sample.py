import json

import numpy
import numpy.random
import scipy.stats

from xboa.hit import Hit
from xboa.bunch import Bunch
import xboa.common

class SmearAndSample(object):
    def __init__(self, input_file, output_file, fmt, n_events):
        self.fsrc = open(input_file)
        self.output_file = output_file
        self.format = fmt
        self.n_events = n_events
        self.tgt_dist = None
        self.smear()
        self.write()

    def smear(self):
        src_dist = numpy.array([json.loads(line) for line in self.fsrc]).transpose()
        kernel = scipy.stats.gaussian_kde(src_dist)
        self.tgt_dist = kernel.resample(self.n_events).transpose()

    def write(self):
        hit_list = []
        for item in self.tgt_dist:
            my_dict = {"mass":self.mu_mass, "pid":-13}
            for i, key in enumerate(self.keys):
                my_dict[key] = item[i]
            hit_list.append(Hit.new_from_dict(my_dict, "energy"))
        bunch = Bunch.new_from_hits(hit_list)
        bunch.hit_write_builtin(self.format, self.output_file)
        print "Writing bunch of length", len(bunch)

    keys = ("x", "px", "y", "py", "pz")
    mu_mass = xboa.common.pdg_pid_to_mass[13]

def smear_test():
    src_fname = '/tmp/smear_test.in'
    tgt_fname = '/tmp/smear_test.out'
    tgt_nevents = 10000

    cov = numpy.array([
              [1000.,   0.,     0.,     0.,     0.,],
              [0.,    600.,     0.,     0.,     0.,],
              [0.,      0.,  1000.,     0.,     0.,],
              [0.,      0.,     0.,   600.,     0.,],
              [0.,      0.,     0.,     0.,   100.,],
          ])
    mean = numpy.array([0., 0., 0., 0., 0.])
    src_dist = numpy.random.multivariate_normal(mean, cov, 1000)
    fsrc = open(src_fname, 'w')
    for item in src_dist:
        fsrc.write(json.dumps(item.tolist())+'\n')
    fsrc.close()

    SmearAndSample(src_fname, tgt_fname, "icool_for003", tgt_nevents)

    fin = open(tgt_fname)
    tgt_bunch = Bunch.new_from_read_builtin("icool_for003", tgt_fname)
    tgt_dist = [[hit[key] for key in SmearAndSample.keys] for hit in tgt_bunch]

    if len(tgt_dist) != tgt_nevents:
        raise RuntimeError("Fail")

    tgt_dist = numpy.array(tgt_dist).transpose()
    numpy.set_printoptions(precision=8, linewidth=1000)
    print numpy.cov(src_dist.transpose())
    print
    print numpy.cov(tgt_dist)


if __name__ == "__main__":
    smear_test()