import sys
import numpy
import bisect
import time
import copy

from xboa.hit import Hit

class AmplitudeDataBinless(object):
    def __init__(self, file_name, bin_edges, mass, cov_fixed):
        """
        Initialise the AmplitudeData class for the basic amplitude calculation
        * File name is the name of the (temporary) file to which the amplitude data is written
        * bin_edges is a list of bin edges for grouping the amplitude data; must have a regular bin step size
        * mass is the particle mass used for calculating emittance

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
        self.mass = mass # for emittance calculation
        self.min_events = 1000 # minimum number of events in the cov calculation
        self.save_frequency = 100 # frequency with which we save the state
        self.cov_fixed = cov_fixed
        self.state_list = []
        self.clear()

    def clear(self):
        """
        Set up the data, with empty arrays. See __init__ doc for definitions
        """
        n_var = 4
        n_samples = 2
        self.run_array = [numpy.memmap(self.file_name+"_run_"+str(i), dtype='int32', mode='w+', shape=(1)) for i in range(n_samples)]
        self.spill_array = [numpy.memmap(self.file_name+"_spill_"+str(i), dtype='int32', mode='w+', shape=(1)) for i in range(n_samples)]
        self.event_array = [numpy.memmap(self.file_name+"_event_"+str(i), dtype='int32', mode='w+', shape=(1)) for i in range(n_samples)]
        self.ps_matrix = [numpy.memmap(self.file_name+"_ps_"+str(i), dtype='float32', mode='w+', shape=(1, n_var)) for i in range(n_samples)]
        self.amplitude = [numpy.memmap(self.file_name+"_amp_"+str(i), dtype='float32', mode='w+', shape=(1)) for i in range(n_samples)]
        self.n_events = [0 for i in range(n_samples)]
        self.n_cov_events = 0.
        self.amplitude_list = []
        self.mean = numpy.array([0. for i in range(n_var)])
        self.cov = numpy.array([[0. for i in range(n_var)] for j in range(n_var)])
        self.cov_inv = numpy.array([[0. for i in range(n_var)] for j in range(n_var)])
        self.state_list = []
        self.emittance = 0.

    def append_hits(self, hit_list):
        """
        Add a set of data to the amplitude_data structure
        * hit_list: list of hits, each formatted as an xboa.Hit
        
        Splits the list into two samples and fills data.
        """
        var_list = ['x', 'px', 'y', 'py']
        if self.n_events[0] > self.n_events[1]:
            bin_1, bin_2 = 1, 0
        else:
            bin_1, bin_2 = 0, 1

        # half the events
        temp_hit_list = hit_list[0::2]
        run_array = numpy.array([hit['particle_number'] for hit in temp_hit_list])
        spill_array = numpy.array([hit['spill'] for hit in temp_hit_list])
        event_array = numpy.array([hit['event_number'] for hit in temp_hit_list])
        ps_array = numpy.array([[hit[var] for var in var_list] for hit in temp_hit_list])
        amp_array = numpy.array([0. for hit in temp_hit_list])
        self.append(run_array, spill_array, event_array, ps_array, amp_array, bin_1)

        # other half of the events
        temp_hit_list = hit_list[1::2]
        run_array = numpy.array([hit['particle_number'] for hit in temp_hit_list])
        spill_array = numpy.array([hit['spill'] for hit in temp_hit_list])
        event_array = numpy.array([hit['event_number'] for hit in temp_hit_list])
        ps_array = numpy.array([[hit[var] for var in var_list] for hit in temp_hit_list])
        amp_array = numpy.array([0. for hit in temp_hit_list])
        self.append(run_array, spill_array, event_array, ps_array, amp_array, bin_2)

    def append(self, run_array, spill_array, event_array, ps_matrix, amp_array, sample):
        """
        Append data to a particular sample
        * run_array: array of run numbers (ints)
        * spill_array: array of spill numbers (ints)
        * event_array: array of particle event numbers (ints)
        * ps_matrix: array of 4-vectors like (x, px, y, py)
        * amp_array: array of floats with particle amplitudes (this gets 
          overwritten during fractional_amplitude method call)
        * sample: the sample to which the data should be stored
        """
        if spill_array.shape[0] == 0:
              return # nothing to do
        self.n_events[sample] += spill_array.shape[0]
        self.get_data(sample, 'r+')
        self.run_array[sample][-run_array.shape[0]:] = run_array[:]
        self.spill_array[sample][-spill_array.shape[0]:] = spill_array[:]
        self.event_array[sample][-event_array.shape[0]:] = event_array[:]
        self.ps_matrix[sample][-ps_matrix.shape[0]:] = ps_matrix[:]
        self.amplitude[sample][-amp_array.shape[0]:] = amp_array[:]
        self.run_array[sample].flush()
        self.spill_array[sample].flush()
        self.event_array[sample].flush()
        self.ps_matrix[sample].flush()
        self.amplitude[sample].flush()

    def delete(self, index_list, sample):
        """
        Clear some data from the sample
        * index_list: indices of events to clear; integer indices corresponding
          to the position in the run_array of the elements that will be removed
        * sample: the sample from which the events will be cleared
        Clears events from run_array, spill_array, etc and updates n_events
        """
        if len(index_list) == 0:
            return
        run_array, spill_array, event_array, ps_matrix = self.get_data(sample, 'r+')
        self.run_array[sample] = delete_elements(run_array, index_list)
        self.spill_array[sample] = delete_elements(spill_array, index_list)
        self.event_array[sample] = delete_elements(event_array, index_list)
        self.ps_matrix[sample] = delete_elements(ps_matrix, index_list)
        self.amplitude[sample] = delete_elements(amplitude, index_list)
        self.n_events[sample] -= len(index_list)

    def retrieve(self, sample):
        """
        Returns a generator for events from the given sample
        * sample: integer corresponding to teh event from which data should be
        retrieved.
        Does not change the sample (i.e. not a "pop" type operation).
        """
        run_array, spill_array, event_array, ps_matrix, amplitude = self.get_data(sample, 'r')
        for i, element in enumerate(ps_matrix):
            yield run_array[i], spill_array[i], event_array[i], element, amplitude[i]

    def retrieve_one(self, sample, i):
        """
        Returns a single event from the given sample
        * sample: integer corresponding to the sample from which data should be
        retrieved.
        * i: index of the event to be retrieved.
        Does not change the sample (i.e. not a "pop" type operation).
        """
        run_array, spill_array, event_array, ps_matrix, amplitude = self.get_data(sample, 'r')
        return run_array[i], spill_array[i], event_array[i], element, amplitude[i]

    def get_data(self, sample, mode):
        """
        Get data from the given sample
        * sample: integer corresponding to the sample
        * mode: mode in which the data should be returned; "r" read; "r+" 
          read and write; "w" overwrite existing data; "w+" write exisiting data
        Returns a tuple of the run_array, spill_array, event_array, ps_matrix, amplitude array
        """
        n_events = self.n_events[sample]
        self.run_array[sample] = numpy.memmap(self.file_name+"_run_"+str(sample), dtype='int32', mode=mode, shape=(n_events))
        self.spill_array[sample] = numpy.memmap(self.file_name+"_spill_"+str(sample), dtype='int32', mode=mode, shape=(n_events))
        self.event_array[sample] = numpy.memmap(self.file_name+"_event_"+str(sample), dtype='int32', mode=mode, shape=(n_events))
        self.ps_matrix[sample] = numpy.memmap(self.file_name+"_ps_"+str(sample), dtype='float32', mode=mode, shape=(n_events, 4))
        self.amplitude[sample] = numpy.memmap(self.file_name+"_amp_"+str(sample), dtype='float32', mode=mode, shape=(n_events))

        return self.run_array[sample], self.spill_array[sample], self.event_array[sample], self.ps_matrix[sample], self.amplitude[sample]

    def get_cov_matrix(self, sample):
        """
        Calculate the covariance matrix and emittance for the given sample; fill
        self.emittance, self.cov, self.cov_inv, self.n_cov_events with the
        relevant data.
        * sample: integer corresponding to the sample
        Returns None
        """
        print "Calculating cov_matrix",
        self.n_cov_events = 0.
        n_var = 4
        var_list = range(4)
        self.mean = numpy.array([0. for i in range(n_var)])
        self.cov = numpy.array([[0. for i in range(n_var)] for j in range(n_var)])
        for run, spill, event, psv, amplitude in self.retrieve(sample):
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
        print "... found", self.n_cov_events, "events"
        sys.stdout.flush()
        self.cov_inv = numpy.linalg.inv(self.cov)
        self.get_emittance()

    def remove_from_cov_matrix(self, update_psv_array):
        """
        Remove elements in the psv_array from the covariance matrix
        * update_psv_array: array of 4-vectors like (x, px, y, py) with data
          that should be removed from the covariance matrix and associated
          variables
        Updates self.emittance, self.cov, self.cov_inv, self.n_cov_events
        Returns None
        """
        if update_psv_array.shape[0] == 0:
            return
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
        self.cov_inv = numpy.linalg.inv(self.cov)
        self.get_emittance()

    def get_emittance(self):
        """
        Update the emittance calculation
        """
        self.emittance = numpy.linalg.det(self.cov)**0.25/self.mass

    def get_amplitude(self, sample, i):
        """
        Calculate the amplitude of a given event
        * sample: the sample from which the event should be taken
        * i: index of th event in the sample whose amplitude should be 
             calculated
        Uses the self.cov_inv to calculate amplitude (as calculated in
        *_cov_matrix functions)
        """
        psv = [uj-self.mean[j] for j, uj in enumerate(self.ps_matrix[sample][i])]
        psv_t = numpy.transpose(psv)
        amplitude = self.emittance*numpy.dot(numpy.dot(psv_t, self.cov_inv), psv)
        return amplitude

    def fractional_amplitude_sample(self, ref_sample, test_sample):
        """
        calculate the covariance matrix
        * ref_sample: reference sample, used for beam ellipse calculation
        * test_sample: test sample, whose amplitude will be calculated
        The idea of having two samples is that we can't get a bias where the
        reference sample is pulled by the statistical fluctuations from e.g.
        scattering in the test sample
        
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
        calculate covariance matrix of "ref bin"
        while number of events in "ref bin" > 10
            calculate amplitudes in "ref bin"
            designate highest amplitude in the "ref bin" as "amp cut"
            remove highest amplitude event from the "ref bin"
            update covariance matrix
            loop over "test bin"
                calculate amplitudes
                if amplitude > "amp cut", remove event from "test bin" and store the amplitude
        swap the "ref bin" and "test bin" designation and repeat
        """
        self.get_cov_matrix(ref_sample)
        test_used = [False]*self.n_events[test_sample]
        ref_used = [False]*self.n_events[ref_sample]
        amplitude_list = [None]*self.n_events[test_sample]
        then = time.time()
        sample_count = 0
        if self.n_cov_events <= self.min_events:
            self.min_events = self.n_cov_events - 1
        while self.n_cov_events > self.min_events:
            amplitude_cut = -1
            max_event = -1
            for i in range(self.n_events[ref_sample]):
                if ref_used[i]:
                    continue
                amplitude = self.get_amplitude(ref_sample, i)
                if amplitude > amplitude_cut:
                    amplitude_cut = amplitude
                    max_event = i
            try:
                ref_used[max_event] = True
                sample_count += 1
                self.remove_from_cov_matrix(numpy.array([self.ps_matrix[ref_sample][max_event]]))
            except Exception:
                sys.excepthook(*sys.exc_info())
                break
            if sample_count % self.save_frequency == 0:
                self.save_state()

            for i in range(self.n_events[test_sample]):
                if test_used[i]:
                    continue
                amplitude_list[i] = self.get_amplitude(test_sample, i)
                if amplitude_list[i] > amplitude_cut:
                    test_used[i] = True # freeze the amplitude
            
            if time.time() - then > 60:
                then = time.time()
                print "    ...", self.n_cov_events, "events remaining"
                sys.stdout.flush()
        print "Finished sample execution with ref", ref_sample, "test", \
              test_sample, "and", self.n_cov_events, "remaining"
        self.get_data(test_sample, 'r+')  # go into write mode
        self.amplitude[test_sample][:] = numpy.array(amplitude_list, dtype='float32')[:]

    def set_fixed_cov_matrix(self):
        n_var = 4
        self.n_cov_events = 0.
        self.mean = numpy.array([0. for i in range(n_var)])
        self.cov = numpy.array(self.cov_fixed)
        self.cov_inv = numpy.linalg.inv(self.cov)
        self.get_emittance()

    def fractional_amplitude_fixed(self):
        self.set_fixed_cov_matrix()
        for sample in [0, 1]:
            amplitude_list = [None]*self.n_events[sample]
            for i in range(self.n_events[sample]):
                amplitude_list[i] = self.get_amplitude(sample, i)
            self.get_data(sample, 'r+')  # go into write mode
            self.amplitude[sample][:] = numpy.array(amplitude_list, dtype='float32')[:]

    def fractional_amplitude(self):
        """
        Calculate the fractional amplitude for both samples, using the algorithm
        described in fractional_amplitude_sample, alternating test and reference
        samples.

        Returns a dictionary mapping (run, spill, event) tuples to amplitude
        """
        print "Calculating amplitudes with number/sample:", self.n_events
        if self.cov_fixed:
            self.fractional_amplitude_fixed()
        else:
            self.fractional_amplitude_sample(0, 1)
            self.fractional_amplitude_sample(1, 0)
        amplitude_list = []
        for sample in range(len(self.n_events)):
            runs, spills, events, psvs, amplitudes = self.get_data(sample, 'r')
            event_id_list = zip(runs, spills, events)
            amplitude_list += zip(event_id_list, amplitudes)
        amplitude_dict = dict(amplitude_list)
        print "Registered", len(amplitude_dict), "amplitudes"
        return amplitude_dict

    def save_state(self):
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
        self.state_list.append(my_state)

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
    mean = numpy.array([10., 10., 0., 0.])
    n_events_per_spill = 100
    cov_1 = xboa.bunch.Bunch.build_penn_ellipse(6., xboa.common.pdg_pid_to_mass[13], 333, 0., 200, 0., 0.004, 1)
    cov_2 = xboa.bunch.Bunch.build_penn_ellipse(20., xboa.common.pdg_pid_to_mass[13], 666, -1., 200, 0., 0.004, 1)
    amp = AmplitudeData("/tmp/amplitude_analysis.tmp", [i * 5. for i in range(21)], 105.658)
    for cov in cov_1,:# cov_2:
        for spill in range(n_spills):
            ps_data = []
            for event in range(n_events_per_spill):
                ps_vector = numpy.random.multivariate_normal(mean, cov)
                event_data = ps_vector.tolist()
                ps_data.append(event_data)
            spill_data = numpy.array([spill for i in range(n_events_per_spill)])
            event_data = numpy.array(range(n_events_per_spill))
            ps_data = numpy.array(ps_data)
            amp.append(spill_data, event_data, ps_data, 20)
    start_time = time.time()
    amp.fractional_amplitude()
    return (time.time() - start_time)/60


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

if __name__ == "__main__":
    #test_memmap()
    results_list = []
    for n_spills in range(10, 11, 10):
        delta_t = test_fractional_amplitude(n_spills)
        results_list.append((n_spills, delta_t))
    print "N spills".ljust(10), "Time [mins]".rjust(10)
    for result in results_list:
        print str(result[0]).ljust(10), str(result[1]).rjust(10)



