import numpy

class DensityData(object):
    def __init__(self, file_name):
        """
        Initialise the DensityData class for the density estimation
        * File name is the name of the (temporary) file to which the density data is written

        Internally, we store using numpy.memmap (which is a buffered file-based 
        array thingy). We have several memmaps:
        * run_array stores the run number of each event
        * spill_array stores the spill number of each event
        * event_array stores the event number of each event
        * ps_matrix stores the 4D phase space vector of each event
        Each one is stored in a memmap file, based on "file_name"+suffix
        During the density calculation, for events we have
        not yet finished counting:
        * number of events in the covariance matrix
        """
        self.file_name = file_name
        self.save_frequency = 100 # frequency with which we save the state
        self.clear()

    def clear(self):
        """
        Set up the data, with empty arrays. See __init__ doc for definitions
        """
        n_var = 4
        self.run_array = numpy.memmap(self.file_name+"_run", dtype='int32', mode='w+', shape=(1))
        self.spill_array = numpy.memmap(self.file_name+"_spill", dtype='int32', mode='w+', shape=(1))
        self.event_array = numpy.memmap(self.file_name+"_event", dtype='int32', mode='w+', shape=(1))
        self.ps_matrix = numpy.memmap(self.file_name+"_ps", dtype='float32', mode='w+', shape=(1, n_var))
        self.n_events = 0
        self.state_list = []

    def append_hits(self, hit_list):
        """
        Add a set of data to the density_data structure
        * hit_list: list of hits, each formatted as an xboa.Hit
        
        Fills data.
        """

        # All the events
        run_array = numpy.array([hit['particle_number'] for hit in hit_list])
        spill_array = numpy.array([hit['spill'] for hit in hit_list])
        event_array = numpy.array([hit['event_number'] for hit in hit_list])

        var_list = ['x', 'px', 'y', 'py']
        ps_array = numpy.array([[hit[var] for var in var_list] for hit in hit_list])

        self.append(run_array, spill_array, event_array, ps_array)

    def append(self, run_array, spill_array, event_array, ps_matrix):
        """
        Append data
        * run_array: array of run numbers (ints)
        * spill_array: array of spill numbers (ints)
        * event_array: array of particle event numbers (ints)
        * ps_matrix: array of 4-vectors like (x, px, y, py)
        """
        if spill_array.shape[0] == 0:
              return # nothing to do
        self.n_events += spill_array.shape[0]
        self.get_data('r+')
        self.run_array[-run_array.shape[0]:] = run_array[:]
        self.spill_array[-spill_array.shape[0]:] = spill_array[:]
        self.event_array[-event_array.shape[0]:] = event_array[:]
        self.ps_matrix[-ps_matrix.shape[0]:] = ps_matrix[:]
        self.run_array.flush()
        self.spill_array.flush()
        self.event_array.flush()
        self.ps_matrix.flush()

    def delete(self, index_list):
        """
        Clear some data
        * index_list: indices of events to clear; integer indices corresponding
          to the position in the run_array of the elements that will be removed
        Clears events from run_array, spill_array, etc and updates n_events
        """
        if len(index_list) == 0:
            return
        run_array, spill_array, event_array, ps_matrix = self.get_data('r+')
        self.run_array = delete_elements(run_array, index_list)
        self.spill_array = delete_elements(spill_array, index_list)
        self.event_array = delete_elements(event_array, index_list)
        self.ps_matrix = delete_elements(ps_matrix, index_list)
        self.n_events -= len(index_list)

    def delete_elements(self, source_map, elements):
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

    def retrieve(self):
        """
        Returns a generator for events
        Does not change the data (i.e. not a "pop" type operation).
        """
        run_array, spill_array, event_array, ps_matrix = self.get_data('r')
        for i, element in enumerate(ps_matrix):
            yield run_array[i], spill_array[i], event_array[i], element

    def retrieve_one(self, i):
        """
        Returns a single event
        * i: index of the event to be retrieved.
        Does not change the data (i.e. not a "pop" type operation).
        """
        run_array, spill_array, event_array, ps_matrix = self.get_data('r')
        return run_array[i], spill_array[i], event_array[i], element

    def get_data(self, mode):
        """
        Get data from the given sample
        * mode: mode in which the data should be returned; "r" read; "r+" 
          read and write; "w" overwrite existing data; "w+" write exisiting data
        Returns a tuple of the run_array, spill_array, event_array, ps_matrix
        """
        n_events = self.n_events
        self.run_array = numpy.memmap(self.file_name+"_run", dtype='int32', mode=mode, shape=(n_events))
        self.spill_array = numpy.memmap(self.file_name+"_spill", dtype='int32', mode=mode, shape=(n_events))
        self.event_array = numpy.memmap(self.file_name+"_event", dtype='int32', mode=mode, shape=(n_events))
        self.ps_matrix = numpy.memmap(self.file_name+"_ps", dtype='float32', mode=mode, shape=(n_events, 4))

        return self.run_array, self.spill_array, self.event_array, self.ps_matrix

    def save_state(self):
        """
        Record internal state for plotting/diagnostics etc

        Append to self.state_list an array containing:
          - "n_events" number of events surviving
        """
        my_state = {
            "n_events":self.n_cov_events,
        }
        self.state_list.append(my_state)

    def get_n_events(self):
        """
        Returns the total number of event stored
        """
        return self.n_events
