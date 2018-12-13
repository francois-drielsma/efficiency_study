import unittest
import importlib
import sys
import numpy

from mice_analysis.amplitude_analysis import AmplitudeAnalysis

def load_config():
    global config_mod
    if config_mod == None:
       print "Usage:"
       print "echo 'scripts/config/config_mc.py' | "+\
             "python scripts/mice_analysis/amplitude/test_amplitude_analysis.py"
    config_mod = config_mod.replace(".py", "")
    config_mod = config_mod.replace("scripts/", "")
    config_mod = config_mod.replace("/", ".")
    print "Using configuration module", config_mod
    config_file = importlib.import_module(config_mod)
    return config_file.Config()

class DataLoaderMockup(object):
    def __init__(self):
        self.events = []

class TestAmplitudeAnalysis(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = load_config()

    def setUp(self):
        self.data_loader_mockup = DataLoaderMockup()
        self.amps = AmplitudeAnalysis(self.config,
                                      self.config.analyses[0],
                                      self.data_loader_mockup)

    def matrix_calc(self, n_matrix, a_in, a_out):
        n_rows = 21
        row_sum = [0. for i in range(n_rows)]
        col_sum = [0. for i in range(n_rows)]
        for i in range(n_rows):
            for j in range(n_rows):
                row_sum[i] += n_matrix[i][j]
                col_sum[j] += n_matrix[i][j]
        for i in range(n_rows):
            for j in range(n_rows):
                n_matrix[i][j] /= row_sum[i]
        matrix_transpose = numpy.transpose(numpy.array(n_matrix))
        a_out_2 = numpy.dot(matrix_transpose, numpy.array(a_in)).tolist()


        print self.amps.matrix_str(n_matrix)
        print "row sum", row_sum
        print "a in   ", a_in
        print "col sum", col_sum
        print "a out  ", a_out
        print "a out 2", a_out_2

    def test_get_pdfs(self):
        print "\n"
        events_in = {}
        for i in range(100000):
            events_in[i] = numpy.random.uniform(0, 104.9)
        pdfs_in = self.amps.get_pdfs(events_in)["pdf"]
        for bin in pdfs_in:
            self.assertLess(bin, 6000)
            self.assertGreater(bin, 4000)
        for i in range(100000):
            events_in[i] = numpy.random.uniform(0, 4.9)
        pdfs_in = self.amps.get_pdfs(events_in)["pdf"]
        self.assertAlmostEqual(pdfs_in[0], 100000)
        for i in range(100000):
            events_in[i] = numpy.random.uniform(85.0001, 89.9999)
        pdfs_in = self.amps.get_pdfs(events_in)["pdf"]
        self.assertAlmostEqual(pdfs_in[17], 100000)

    def test_migration_matrix(self):
        print "\n"
        events_in = {}
        for i in range(10000):
            events_in[i] = numpy.random.uniform(0, 50.9)

        events_out = {}
        for i in range(10000):
            events_out[i] = events_in[i] + numpy.random.uniform(-1e1, 1e1)
            while events_out[i] < 0.:
                events_out[i] += 10.

        pdfs_in = numpy.array(self.amps.get_pdfs(events_in)["pdf"])
        pdfs_out = numpy.array(self.amps.get_pdfs(events_out)["pdf"])
        migration_matrix = numpy.array(self.amps.migration_matrix(events_in, events_out, True))
        pdfs_recovered = numpy.dot(numpy.transpose(migration_matrix), pdfs_in)

        print pdfs_in.tolist()
        print pdfs_out.tolist()
        print self.amps.matrix_str(migration_matrix.tolist())
        print pdfs_recovered.tolist()

    def test_performance_systematics(self):
        self.amps.amplitudes["reco"] = {
              "all_upstream":{"pdf":[1000/float(i) for i in range(8)]+[0., 0.]},
              "all_downstream":{"pdf":[1000/float(i) for i in range(8)]+[0., 0.]},
              "migration_matrix":[[float(i == j) for i in range(10)] for j in range(10)],
        }
        self.amps.amplitudes["systematics_performance_reference"] = None
        sys = self.amps.calculate_performance_systematics("all_upstream")
        [self.assertAlmostEqual(item, 0.) for item in sys]
        test_matrix = self.get_matrix()


if __name__ == "__main__":
    config_mod = raw_input()
    unittest.main()



