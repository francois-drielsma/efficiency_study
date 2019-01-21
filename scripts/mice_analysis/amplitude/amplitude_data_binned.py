import sys
import numpy
import bisect
import copy
import itertools
import math
import time
import resource

import ROOT

from xboa.hit import Hit
import xboa.common

class AmplitudeDataBinned(object):
    def __init__(self, file_name, bin_edges, mass, cov_fixed):
        """
        Initialise the AmplitudeData class for the basic amplitude calculation
        * File name is the name of the (temporary) file to which the amplitude data is written
        * bin_edges is a list of bin edges for grouping the amplitude data; must have a regular bin step size
        * mass is the particle mass used for calculating emittance
        * cov_fixed is not used

        Internally, we store using numpy.memmap (which is a buffered file-based 
        array thingy). We have several memmaps for each bin:
        * run_array stores the run number of each event in the bin
        * spill_array stores the spill number of each event in the bin
        * event_array stores the event number of each event in the bin
        * ps_matrix stores the 4D phase space vector of each event in the bin
        Each one is stored in a memmap file, based on "file_name"+suffix
        During the amplitude calculation, for events whose amplitude we have
        not yet finished counting:
        * mean, cov (covariance matrix) and cov_inv (inverse covariance matrix)
        * number of events in the covariance matrix
        And:
        * number of events in each bin
        We build, during the calculation, an "amplitude_list"
        * tuple of id, amplitude like ((run, spill, event), amplitude)
        """
        self.file_name = file_name
        self.bins = bin_edges
        self.min_events = 1000 # minimum number of events in the cov calculation
        self.min_bin = 4 # minimum bin to apply amplitude iterative procedure
        self.n_bins = len(self.bins) # we have an overflow bin
        self.n_samples = 2
        for i in range(self.n_bins-1):
            if self.bins[i+1] - self.bins[i] < 1e-6:
                print "Found bins:", self.bins
                raise ValueError("Bin spacing is too small for bin "+str(i))
        self.mass = mass # for emittance calculation
        self.state_list = [[], []]
        file_limit = resource.getrlimit(resource.RLIMIT_NOFILE)
        if file_limit[0] < 2048:
            resource.setrlimit(resource.RLIMIT_NOFILE, (2048, file_limit[1]))
            print "Adjusted number of allowed file handles from (soft, hard)", \
                  file_limit, "to", resource.getrlimit(resource.RLIMIT_NOFILE)

    def __del__(self):
        pass

    def clear(self):
        n_var = 4
        self.run_array = [[numpy.memmap(self.file_name+"_run_"+str(i)+"_"+str(j), dtype='int32', mode='w+', shape=(1)) for i in range(self.n_bins)] for j in range(2)]
        self.spill_array = [[numpy.memmap(self.file_name+"_spill_"+str(i)+"_"+str(j), dtype='int32', mode='w+', shape=(1)) for i in range(self.n_bins)] for j in range(2)]
        self.event_array = [[numpy.memmap(self.file_name+"_event_"+str(i)+"_"+str(j), dtype='int32', mode='w+', shape=(1)) for i in range(self.n_bins)] for j in range(2)]
        self.ps_matrix = [[numpy.memmap(self.file_name+"_ps_"+str(i)+"_"+str(j), dtype='float32', mode='w+', shape=(1, n_var)) for i in range(self.n_bins)] for j in range(2)]
        self.amp_array = [[numpy.memmap(self.file_name+"_amp_"+str(i)+"_"+str(j), dtype='float32', mode='w+', shape=(1)) for i in range(self.n_bins)] for j in range(2)]
        self.n_events = [[0 for i in range(self.n_bins)] for j in range(2)]
        self.n_cov_events = 0.
        self.amplitude_list = []
        self.state_list = [[], []]
        self.mean = numpy.array([0. for i in range(n_var)])
        self.cov = numpy.array([[0. for i in range(n_var)] for j in range(n_var)])
        self.cov_inv = numpy.array([[0. for i in range(n_var)] for j in range(n_var)])

    def append_hits(self, hit_list):
        var_list = ['x', 'px', 'y', 'py']
        for index in [0, 1]:
            tmp = [hit for i, hit in enumerate(hit_list) if i%2 == index]
            run_array = numpy.array([hit['particle_number'] for hit in tmp])
            spill_array = numpy.array([hit['spill'] for hit in tmp])
            event_array = numpy.array([hit['event_number'] for hit in tmp])
            ps_array = numpy.array([[hit[var] for var in var_list] for hit in tmp])
            amp_array = numpy.array([0. for hit in tmp])
            a_bin = self.n_bins-1
            self.append(run_array, spill_array, event_array, ps_array, amp_array, a_bin, index)

    def append(self, run_array, spill_array, event_array, ps_matrix, amp_array, a_bin, sample):
        if spill_array.shape[0] == 0:
              return # nothing to do
        self.n_events[sample][a_bin] += spill_array.shape[0]
        self.get_data(a_bin, sample, 'r+')
        self.run_array[sample][a_bin][-run_array.shape[0]:] = run_array[:]
        self.spill_array[sample][a_bin][-spill_array.shape[0]:] = spill_array[:]
        self.event_array[sample][a_bin][-event_array.shape[0]:] = event_array[:]
        self.ps_matrix[sample][a_bin][-ps_matrix.shape[0]:] = ps_matrix[:]
        self.amp_array[sample][a_bin][-amp_array.shape[0]:] = amp_array[:]
        self.run_array[sample][a_bin].flush()
        self.spill_array[sample][a_bin].flush()
        self.event_array[sample][a_bin].flush()
        self.ps_matrix[sample][a_bin].flush()
        self.amp_array[sample][a_bin].flush()

    def delete(self, index_list, a_bin, sample):
        if len(index_list) == 0:
            return
        run_array, spill_array, event_array, ps_matrix, amp_array = self.get_data(a_bin, sample, 'r+')
        self.run_array[sample][a_bin] = delete_elements(run_array, index_list)
        self.spill_array[sample][a_bin] = delete_elements(spill_array, index_list)
        self.event_array[sample][a_bin] = delete_elements(event_array, index_list)
        self.ps_matrix[sample][a_bin] = delete_elements(ps_matrix, index_list)
        self.amp_array[sample][a_bin] = delete_elements(amp_array, index_list)
        self.n_events[sample][a_bin] -= len(index_list)

    def retrieve(self, a_bin, sample):
        run_array, spill_array, event_array, ps_matrix, amp_array = self.get_data(a_bin, sample, 'r')
        for i, element in enumerate(ps_matrix):
            yield run_array[i], spill_array[i], event_array[i], element, amp_array[i]

    def get_data(self, a_bin, sample, mode):
        n_events = self.n_events[sample][a_bin]
        suffix = str(a_bin)+"_"+str(sample)
        self.run_array[sample][a_bin] = numpy.memmap(self.file_name+"_run_"+suffix, dtype='int32', mode=mode, shape=(n_events))
        self.spill_array[sample][a_bin] = numpy.memmap(self.file_name+"_spill_"+suffix, dtype='int32', mode=mode, shape=(n_events))
        self.event_array[sample][a_bin] = numpy.memmap(self.file_name+"_event_"+suffix, dtype='int32', mode=mode, shape=(n_events))
        self.ps_matrix[sample][a_bin] = numpy.memmap(self.file_name+"_ps_"+suffix, dtype='float32', mode=mode, shape=(n_events, 4))
        self.amp_array[sample][a_bin] = numpy.memmap(self.file_name+"_amp_"+suffix, dtype='float32', mode=mode, shape=(n_events))
        return self.run_array[sample][a_bin], self.spill_array[sample][a_bin], self.event_array[sample][a_bin], self.ps_matrix[sample][a_bin], self.amp_array[sample][a_bin]

    def set_amplitudes(self, a_bin, sample, amp_array):
        self.get_data(a_bin, sample, 'r+')
        #print "Set ampplitudes for bin, sample", a_bin, sample, "Before"
        #print "    ", self.amp_array[sample][a_bin]
        #print "New data"
        #print "    ", amp_array
        self.amp_array[sample][a_bin][:] = amp_array[:]
        self.amp_array[sample][a_bin].flush()
        #print "After"
        #print "    ", self.amp_array[sample][a_bin]
        if self.amp_array[sample][a_bin].shape[0] != self.n_events[sample][a_bin]:
            raise ValueError("Failed to set amplitudes")

    def get_cov_matrix(self, sample):
        print "Calculating cov_matrix"
        self.n_cov_events = 0.
        n_var = 4
        var_list = range(4)
        self.mean = numpy.array([0. for i in range(n_var)])
        self.cov = numpy.array([[0. for i in range(n_var)] for j in range(n_var)])
        for a_bin in range(self.n_bins):
            for run, spill, event, psv, amp in self.retrieve(a_bin, sample):
                self.n_cov_events += 1.
                if self.n_cov_events % 1000 == 0:
                    print self.n_cov_events,
                    sys.stdout.flush()
                for i in var_list:
                    delta = (self.n_cov_events-1.)/self.n_cov_events
                    self.mean[i] = self.mean[i]*delta + psv[i]/self.n_cov_events
                    for j in range(i, n_var):
                        self.cov[i][j] = self.cov[i][j]*delta + \
                                         psv[i]*psv[j]/self.n_cov_events
        for i in var_list:
            for j in range(i, n_var):
                self.cov[i][j] -= self.mean[i]*self.mean[j]
                self.cov[j][i] = self.cov[i][j]
        print "Done"

    def remove_from_cov_matrix(self, update_psv_array):
        if update_psv_array.shape[0] == 0:
            return
        max_to_remove = int(self.n_cov_events-self.min_events)
        if update_psv_array.shape[0] >= max_to_remove:
            if self.n_cov_events <= self.min_events:
                return
            else:
                update_psv_array = update_psv_array[:max_to_remove]
        #print "Removing", update_psv_array.shape[0], "values from cov matrix"
        n_var = 4
        var_list = range(4)
        # convert to not centred moments
        if update_psv_array.shape[0] >= int(round(self.n_cov_events)):
            for i in var_list:
                self.mean[i] = 0.
                for j in var_list:
                    self.cov[i][j] = 0.
            return

        for i in var_list:
            for j in var_list:
                self.cov[i][j] += self.mean[i]*self.mean[j]

        for psv in update_psv_array:
            delta = self.n_cov_events/(self.n_cov_events-1)            
            for i in var_list:
                self.mean[i] = self.mean[i]*delta - psv[i]/self.n_cov_events
                for j in range(i, n_var):
                    self.cov[i][j] = self.cov[i][j]*delta - \
                                     psv[i]*psv[j]/self.n_cov_events
            self.n_cov_events -= 1.

        for i in var_list:
            for j in range(i, n_var):
                self.cov[i][j] -= self.mean[i]*self.mean[j]
                self.cov[j][i] = self.cov[i][j]

    def get_rebins(self, a_bin, sample, will_store_amplitude):
        self.emittance = self.get_emittance()
        self.cov_inv = numpy.linalg.inv(self.cov)
        rebins = [[] for i in range(self.n_bins)]
        amplitudes = [0 for i in range(self.n_events[sample][a_bin])]
        emittance = self.get_emittance()
        for ev_number, event in enumerate((self.retrieve(a_bin, sample))):
            psv = [ui - self.mean[i] for i, ui in enumerate(event[3])]
            psv_t = numpy.transpose(psv)
            amplitudes[ev_number] = emittance*numpy.dot(numpy.dot(psv_t, self.cov_inv), psv)
            try:
                new_bin = bisect.bisect_right(self.bins, amplitudes[ev_number])-1
            except ValueError:
                print "Error calculating amplitude for bin", a_bin, "event", i, \
                       "data", event[0], event[1], event[2], event[3], "amplitude", amplitudes[ev_number], \
                       "det(cov)", numpy.linalg.det(self.cov), "cov:"
                print self.cov
                raise
            if new_bin != a_bin:
                rebins[new_bin].append(ev_number)
        if will_store_amplitude:
            self.set_amplitudes(a_bin, sample, numpy.array(amplitudes))
        return rebins

    def get_emittance(self):
        emittance = numpy.linalg.det(self.cov)**0.25/self.mass
        return emittance

    def queue_for_amplitude(self, runs, spills, events, psvs):
        new_list = [None]*spills.shape[0]
        emittance = self.get_emittance()
        for i in range(spills.shape[0]):
            psv = [ui - self.mean[i] for i, ui in enumerate(psvs[i])]
            psv_t = numpy.transpose(psv)
            amplitude = emittance*numpy.dot(numpy.dot(psv_t, self.cov_inv), psv)
            new_list[i] = ((runs[i], spills[i], events[i]), amplitude)
        self.amplitude_list += new_list

    def apply_rebins(self, source_bin, rebins, cut_bin, sample):
        """
        Rebin data. If data is moved from inside cut_bin to outside cut_bin, 
        remove it from the covariance matrix calculation.
        """
        all_rebins = []
        for rebin_list in rebins:
            all_rebins += rebin_list
        if len(all_rebins) == 0:
            return
        source_run_array, source_spill_array, source_event_array, source_ps_matrix, source_amp_array = self.get_data(source_bin, sample, 'r+')
        for target_bin, rebin_list in enumerate(rebins):
            runs = numpy.array([source_run_array[i] for i in rebin_list])
            spills = numpy.array([source_spill_array[i] for i in rebin_list])
            events = numpy.array([source_event_array[i] for i in rebin_list])
            psvs = numpy.array([source_ps_matrix[i] for i in rebin_list])
            amps = numpy.array([source_amp_array[i] for i in rebin_list])
            self.append(runs, spills, events, psvs, amps, target_bin, sample)
            # remove from covariance matrix if we are moving things to outside the cut
            if target_bin >= cut_bin and source_bin < cut_bin:
                self.remove_from_cov_matrix(psvs)
        self.delete(all_rebins, source_bin, sample)

    def fractional_amplitude_sample(self, ref_sample, test_sample):
        """
        calculate the covariance matrix
        bin the data
        loop over "cut_bin" from the outer bin inwards
            get data in the cut_bin
            calculate amplitudes in the cut bin
            remove from cov matrix anything in the cut_bin
            while events are still being cut
                loop over bin from the "cut bin"-1 inwards
                    calculate amplitudes
                    rebin events
                    remove from cov matrix anything in the cut_bin
        """
        self.get_cov_matrix(ref_sample)
        for cut_bin in range(self.n_bins):
            rebins = self.get_rebins(cut_bin, ref_sample, False)
            self.apply_rebins(cut_bin, rebins, self.n_bins+1, ref_sample)
        for cut_bin in range(self.n_bins):
            rebins = self.get_rebins(cut_bin, test_sample, True)
            self.apply_rebins(cut_bin, rebins, self.n_bins+1, test_sample)
        self.save_state(ref_sample)
        print "Fractional_amplitudes for ref/test sample", ref_sample, "/", test_sample
        print "Contained emittance", self.get_emittance(), "and", self.n_cov_events, "events, cov matrix"
        print self.cov
        print "Initially"
        print "    Test bins:", self.n_events[test_sample]
        print "    Ref bins:", self.n_events[ref_sample]
        print "===  Bin Emit N_cov Time  ==="
        for cut_bin in reversed(range(self.min_bin, self.n_bins)):
            now = time.time()
            run_array, spill_array, event_array, ps_matrix, amp_array = self.get_data(cut_bin, ref_sample, 'r')
            # remove events from cov matrix that are now in the ref bin
            self.remove_from_cov_matrix(ps_matrix)
            will_rebin = True
            # recalculate amplitudes and rebin the reference sample
            get_rebin_time, apply_rebin_time = 0., 0.
            while will_rebin:
                will_rebin = False
                for a_bin in range(cut_bin):
                    rebins = self.get_rebins(a_bin, ref_sample, False)
                    self.apply_rebins(a_bin, rebins, cut_bin, ref_sample)
                    this_will_rebin = sum([len(a_rebin_list) for a_rebin_list in rebins[cut_bin:]])
                    will_rebin = will_rebin or this_will_rebin
            # recalculate amplitudes and rebin the test sample
            for a_bin in range(cut_bin):
                rebins = self.get_rebins(a_bin, test_sample, True)
                self.apply_rebins(a_bin, rebins, a_bin, test_sample)
            self.save_state(ref_sample)
            print "   ", cut_bin, round(self.get_emittance(), 3), self.n_cov_events, round(time.time()-now)
            now = time.time()
        print "Finally"
        print "    Test bins:", self.n_events[test_sample]
        print "    Ref bins:", self.n_events[ref_sample]

    def fractional_amplitude(self):
        self.fractional_amplitude_sample(0, 1)
        self.fractional_amplitude_sample(1, 0)
        data = []
        for sample in range(2):
            for a_bin in range(self.n_bins):
                my_data = [((run, spill, event), amp) for run, spill, event, psv, amp in self.retrieve(a_bin, sample)]
                data += my_data
        data_dict = dict(data)
        return data_dict

    def save_state(self, sample):
        """
        Record internal state for plotting/diagnostics etc

        Append to self.state_list a dictionary containing:
          - emittance
          - mean
          - cov
          - number of events surviving in the reference sample
        """
        my_state = {
            "emittance":copy.deepcopy(self.emittance),
            "mean":copy.deepcopy(self.mean),
            "cov":copy.deepcopy(self.cov),
            "n_events":self.n_cov_events,
        }
        self.state_list[sample].append(my_state)

    def get_n_events(self):
        n_events = 0
        for bins in self.n_events:
            n_events += sum([a_bin for a_bin in bins])
        return n_events
            

def delete_elements(source_map, elements):
    if len(elements) == 0:
        return source_map
    filename = source_map.filename
    n_events = source_map.shape[0]-len(elements)
    if n_events == 0:
        new_shape = (1,)+source_map.shape[1:]
        new_map = numpy.memmap(filename, dtype=source_map.dtype, mode='w+', shape=new_shape)
        return new_map

    elements = sorted(elements)
    new_shape = (n_events,)+source_map.shape[1:]
    new_map = numpy.memmap(filename+'.tmp', dtype=source_map.dtype, mode='w+', shape=new_shape)
    new_map[:elements[0]] = source_map[:elements[0]]
    for i, next_element in enumerate(elements[1:]):
        prev_element = elements[i]
        new_map[prev_element-i:next_element-i-1] = source_map[prev_element+1:next_element]
    new_map[elements[-1]-len(elements)+1:] = source_map[elements[-1]+1:]
    very_new_map = numpy.memmap(filename, dtype=source_map.dtype, mode='w+', shape=new_shape)
    very_new_map[:] = new_map[:]
    #del new_map
    return very_new_map

def test_fractional_amplitude(n_spills):
    import xboa.common
    import xboa.bunch
    import time
    numpy.random.seed(0)
    mean_0 = numpy.array([0., 0., 0., 0.])
    mean_1 = numpy.array([75., 10., 0., 0.])
    mean_2 = numpy.array([-75., 10., -10., 0.])
    n_events_per_spill = 25
    mass = xboa.common.pdg_pid_to_mass[13]
    cov_1 = xboa.bunch.Bunch.build_penn_ellipse(6., mass, 100, 0., 200, 0., 0.0, 1)
    cov_2 = xboa.bunch.Bunch.build_penn_ellipse(6., mass, 333, -0.5, 200, 0., 0.0, 1)
    amp = AmplitudeDataBinned("/tmp/amplitude_analysis.tmp", [i * 5. for i in range(21)], mass, False)
    if n_spills%2 == 1:
        print "n_spills must be even; changing from", n_spills,
        n_spills += 1
        print "to", n_spills
    for spill in range(n_spills):
        for mean, cov in (mean_0, cov_2), (mean_0, cov_2), (mean_0, cov_2), (mean_0, cov_2):
            ps_data = []
            for event in range(n_events_per_spill):
                ps_vector = numpy.random.multivariate_normal(mean, cov)
                event_data = ps_vector.tolist()
                ps_data.append(event_data)
            run_data = numpy.array([1 for i in range(n_events_per_spill)])
            spill_data = numpy.array([spill for i in range(n_events_per_spill)])
            event_data = numpy.array(range(n_events_per_spill))
            ps_data = numpy.array(ps_data)
            ps_data = numpy.array(ps_data)
            amp_data = numpy.array([0. for i in range(n_events_per_spill)])
            amp.append(run_data, spill_data, event_data, ps_data, amp_data, amp.n_bins-1, spill%2)
    n_events = 4*n_spills*n_events_per_spill
    start_time = time.time()
    amp.min_events = max(n_events/10, 50)
    amp.save_frequency = n_events/10
    amp.fractional_amplitude()
    plot_fractional_amplitude(amp)
    return (time.time() - start_time)/60

def plot_fractional_amplitude(amp_data):
    import plot_amplitude_data
    import utilities.root_style
    plotter = plot_amplitude_data.PlotAmplitudeData(amp_data, "./", "test")
    utilities.root_style.setup_gstyle()
    plotter.plot()
    return
    canvas = ROOT.TCanvas("amplitude test 2d", "amplitude test 2d", 1400, 1000)
    canvas.Divide(2, 1)
    xpx_hist = ROOT.TH2D("amplitude test 2d", ";x [mm];p_{x} [MeV/c]", 50, -100, 100, 70, -50, 50)
    xpx_hist.SetStats(False)
    amp_hist = ROOT.TH1D("amplitude test 1d", ";amplitude [mm];", 12, 0., 60.)
    amp_hist.SetStats(False)
    amp_hist.SetFillColor(ROOT.kOrange-2)
    amp_hist.SetLineColor(ROOT.kBlack)
    for sample in range(2):
        for a_bin in range(21):
            for run, spill, evt, psv, amp in amp_data.retrieve(a_bin, sample):
                xpx_hist.Fill(psv[0], psv[1])
                amp_hist.Fill(amp)
    canvas.cd(1).SetFrameFillColor(utilities.root_style.get_frame_fill())
    xpx_hist.Draw("COLZ")
    print "get ellipse"
    delta_list = [-10, -7, -4, 1, 3]
    for color, sample in [(ROOT.kGreen, 0), (ROOT.kRed, 1)]: 
        step = len(amp_data.state_list[sample])/len(delta_list)+1
        print amp_data.state_list[sample][::step]
        for i, ellipse in enumerate(amp_data.state_list[sample][::step]):
            graph = plot_ellipse(ellipse)
            graph.SetLineColor(color+delta_list[i])
    canvas.cd(2)
    amp_hist.Draw("")
    canvas.Update()
    for fmt in ["root", "png", "pdf"]:
        canvas.Print("TestAmplitudeBinned."+fmt)
    root_objects.append(amp_hist)
    root_objects.append(xpx_hist)
    root_objects.append(canvas)


def plot_ellipse(ellipse):
    mean = ellipse["mean"][0:2]
    cov = ellipse["cov"][0:2, 0:2]
    try:
        points = xboa.common.make_shell(41, cov)
    except Exception:
        graph = ROOT.TGraph()
        return graph
    graph = ROOT.TGraph(len(points)+1)
    points = [(a_point[0, 0], a_point[0, 1]) for a_point in points]
    points = sorted(points, key = lambda points: math.atan2(points[1], points[0]))
    points.append(points[0])
    print ellipse
    for i, a_point in enumerate(points):
        graph.SetPoint(i, a_point[0]+mean[0], a_point[1]+mean[1])
    graph.SetLineWidth(2)
    graph.Draw("L")
    root_objects.append(graph)
    return graph


def test_memmap():
    test_map = numpy.memmap('test', dtype='int32', mode='w+', shape=(20))
    for i in range(20):
        test_map[i] = i
    print test_map
    del test_map
    test_map = numpy.memmap('test', dtype='int32', mode='r+', shape=(20))
    print test_map
    target_elements = [0, 4, 6, 12, 16, 19]
    test_map = delete_elements(test_map, target_elements)
    print test_map
    del test_map
    test_map = numpy.memmap('test', dtype='int32', mode='r+', shape=(20-len(target_elements)))
    print test_map

import cProfile, pstats, StringIO
if __name__ == "__main__":
    #test_memmap()
    results_list = []
    #plot_ellipse({"mean":numpy.array([-1., -1.]), "cov":numpy.array([[1., 0.9,], [0.9, 1.]])})
    #raw_input()
    for n_spills in [10]: #[2, 10, 100]:
        pr = cProfile.Profile()
        pr.enable()
        delta_t = test_fractional_amplitude(n_spills)
        pr.disable()
        s = StringIO.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print s.getvalue()
        results_list.append((n_spills, delta_t))
    print "N spills".ljust(10), "Time [mins]".rjust(10)
    for result in results_list:
        print str(result[0]).ljust(10), str(result[1]).rjust(10)
    raw_input()


