import sys
import numpy

class AmplitudeData(object):
    def __init__(self, file_name, bin_edges, mass):
        """
        Initialise the AmplitudeData class
        * File name is the name of the (temporary) file to which the amplitude data is written
        * bin_edges is a list of bin edges for grouping the amplitude data
        * mass is the particle mass used for calculating emittance
        """
        self.file_name = file_name
        self.bins = bin_edges
        self.mass = mass # for emittance calculation
        self.clear()

    def clear(self):
        self.n_bins = len(self.bins) # we have an overflow bin
        self.spill_array = [numpy.memmap(self.file_name+"_spill_"+str(i), dtype='int32', mode='w+', shape=(1)) for i in range(self.n_bins)]
        self.event_array = [numpy.memmap(self.file_name+"_event_"+str(i), dtype='int32', mode='w+', shape=(1)) for i in range(self.n_bins)]
        self.ps_matrix = [numpy.memmap(self.file_name+"_ps_"+str(i), dtype='float32', mode='w+', shape=(1, 4)) for i in range(self.n_bins)]
        self.n_events = [0 for i in range(self.n_bins)]
        self.n_cov_events = 0.

    def append(self, spill_array, event_array, ps_matrix, bin):
        if spill_array.shape[0] == 0:
              return # nothing to do
        self.n_events[bin] += spill_array.shape[0]
        self.get_data(bin, 'r+')
        self.spill_array[bin][-spill_array.shape[0]:] = spill_array[:]
        self.event_array[bin][-event_array.shape[0]:] = event_array[:]
        self.ps_matrix[bin][-ps_matrix.shape[0]:] = ps_matrix[:]
        self.spill_array[bin].flush()
        self.event_array[bin].flush()
        self.ps_matrix[bin].flush()

    def delete(self, index_list, bin):
        if len(index_list) == 0:
            return
        spill_array, event_array, ps_matrix = self.get_data(bin, 'r+')
        self.spill_array[bin] = delete_elements(spill_array, index_list)
        self.event_array[bin] = delete_elements(event_array, index_list)
        self.ps_matrix[bin] = delete_elements(ps_matrix, index_list)
        self.n_events[bin] -= len(index_list)

    def retrieve(self, bin):
        spill_array, event_array, ps_matrix = self.get_data(bin, 'r')
        for i, element in enumerate(ps_matrix):
            yield spill_array[i], event_array[i], element

    def get_data(self, bin, mode):
        n_events = self.n_events[bin]
        self.spill_array[bin] = numpy.memmap(self.file_name+"_spill_"+str(bin), dtype='int32', mode=mode, shape=(n_events))
        self.event_array[bin] = numpy.memmap(self.file_name+"_event_"+str(bin), dtype='int32', mode=mode, shape=(n_events))
        self.ps_matrix[bin] = numpy.memmap(self.file_name+"_ps_"+str(bin), dtype='float32', mode=mode, shape=(n_events, 4))
        return self.spill_array[bin], self.event_array[bin], self.ps_matrix[bin]

    def get_cov_matrix(self):
        print "Calculating cov_matrix"
        self.n_cov_events = 0.
        n_var = 4
        var_list = range(4)
        self.mean = numpy.array([0. for i in range(n_var)])
        self.cov = numpy.array([[0. for i in range(n_var)] for j in range(n_var)])
        for bin in range(self.n_bins):
            for spill, event, psv in self.retrieve(bin):
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
        #print "Removing", update_psv_array.shape[0], "values from cov matrix"
        n_var = 4
        var_list = range(4)
        # convert to not centred moments
        if update_psv_array.shape[0] >= int(round(self.n_cov_events)):
            for i in var_list:
                self.mean[i] = 0
                for j in var_list:
                    self.cov[i][j] = 0
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

    def get_rebins(self, bin):
        emittance = self.get_emittance()
        bin_step = (self.bins[1]-self.bins[0])/emittance
        cov_inv = numpy.linalg.inv(self.cov)
        rebins = [[] for i in range(self.n_bins)]
        for i, event in enumerate((self.retrieve(bin))):
            psv = event[2]
            psv_t = numpy.transpose(psv)
            amplitude = numpy.dot(numpy.dot(psv_t, cov_inv), psv)
            try:
                new_bin = int(amplitude/bin_step)
            except ValueError:
                print "Error calculating amplitude for bin", bin, "event", i, \
                       "data", event[0], event[1], event[2], "amplitude", amplitude, \
                       "step", bin_step, "det(cov)", numpy.linalg.det(self.cov), "cov:"
                print self.cov
                raise
            new_bin = min(self.n_bins-1, new_bin)
            if new_bin != bin:
                rebins[new_bin].append(i)
        return rebins

    def get_emittance(self):
        emittance = numpy.linalg.det(self.cov)**0.25/self.mass
        return emittance

    def apply_rebins(self, source_bin, rebins, cut_bin):
        """
        Events with rebin >= cut_bin will be removed from self.cov
        """
        all_rebins = []
        for rebin_list in rebins:
            all_rebins += rebin_list
        if len(all_rebins) == 0:
            return

        #print "Apply rebin from", source_bin, "to", [len(a_bin) for a_bin in rebins]
        source_spill_array, source_event_array, source_ps_matrix = self.get_data(source_bin, 'r+')
        for target_bin, rebin_list in enumerate(rebins):
            spills = numpy.array([source_spill_array[i] for i in rebin_list])
            events = numpy.array([source_event_array[i] for i in rebin_list])
            psvs = numpy.array([source_ps_matrix[i] for i in rebin_list])
            self.append(spills, events, psvs, target_bin)
            # remove from covariance matrix if we are moving things to outside the cut
            if target_bin >= cut_bin and source_bin < cut_bin:
                self.remove_from_cov_matrix(psvs)
            # add to covariance matrix if we are moving things to inside the cut
            #elif target_bin < cut_bin and source_bin >= cut_bin:
            #    self.add_to_cov_matrix(psvs)
        #print "Removing", len(all_rebins), "events from bin", source_bin, "which has currently", self.n_events[source_bin], "events"
        #print sorted(all_rebins)
        self.delete(all_rebins, source_bin)

    def fractional_amplitude(self):
        """
        bunch: all particles, on which amplitude is calculated
        amplitude_list: the set of bin edges over which amplitude is calculated
        bunch_delta: the set of particles which make it downstream
        """
        self.get_cov_matrix()
        for bin in reversed(range(self.n_bins)):
            rebins = self.get_rebins(bin)
            self.apply_rebins(bin, rebins, self.n_bins+1)
        for cut_bin in reversed(range(self.n_bins)):
            print "Fractional_amplitudes for cut bin", cut_bin, "and events", self.n_events, "total", sum(self.n_events)
            print "Contained emittance", self.get_emittance(), "and", self.n_cov_events, "events"
            print "cov matrix"
            print self.cov
            spill_array, event_array, ps_matrix = self.get_data(cut_bin, 'r')
            # remove events from cov matrix that are now in the cut bin
            self.remove_from_cov_matrix(ps_matrix)
            will_rebin = True
            # recalculate amplitudes and rebin
            while will_rebin:
                will_rebin = False
                for bin in reversed(range(cut_bin)): # we can't rebin from outside cut to inside cut
                    rebins = self.get_rebins(bin)
                    self.apply_rebins(bin, rebins, cut_bin)
                    this_will_rebin = sum([len(a_rebin_list) for a_rebin_list in rebins])
                    will_rebin = will_rebin or this_will_rebin

    def test(self):
        import xboa.common
        import xboa.bunch
        mean = numpy.array([10., 10., 0., 0.])
        n_events_per_spill = 100
        cov_1 = xboa.bunch.Bunch.build_penn_ellipse(1., xboa.common.pdg_pid_to_mass[13], 333, 0., 200, 0., 0.004, 1)
        cov_2 = xboa.bunch.Bunch.build_penn_ellipse(20., xboa.common.pdg_pid_to_mass[13], 666, -1., 200, 0., 0.004, 1)
        for cov in cov_1, cov_2:
            for spill in range(10):
                ps_data = []
                for event in range(n_events_per_spill):
                    ps_vector = numpy.random.multivariate_normal(mean, cov)
                    event_data = ps_vector.tolist()
                    ps_data.append(event_data)
                spill_data = numpy.array([spill for i in range(n_events_per_spill)])
                event_data = numpy.array(range(n_events_per_spill))
                ps_data = numpy.array(ps_data)
                self.append(spill_data, event_data, ps_data, 20)
        self.fractional_amplitude()

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
    amp = AmplitudeData("/tmp/amplitude_analysis.tmp", [i * 5. for i in range(21)], 105.658)
    amp.test()



