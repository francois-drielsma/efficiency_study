import glob
import json
import numpy
import xboa.common
numpy.set_printoptions( linewidth=1000, precision=4)


class Hacking(object):
    def __init__(self):
        self.bias_json = None
        self.list_of_lists = None
        self.max_bin = 12
        self.canvas_split = (3, 4)

    def setup_lists(self, bias_json, test_keys, test_targets):
        list_of_lists = {}
        for target in test_targets:
            for key_1, key_2 in test_keys:
                item = numpy.array(bias_json[key_1][target][key_2])
                list_of_lists[key_1] = {}
                list_of_lists[key_1][target] = {}
                new_item = numpy.zeros(item.shape).tolist()
                if len(item.shape) == 1:
                    for i in range(item.shape[0]):
                        new_item[i] = []
                elif len(item.shape) == 2:
                    for i in range(item.shape[0]):
                        for j in range(item.shape[1]):
                            new_item[i][j] = []
                list_of_lists[key_1][target][key_2] = new_item
        self.list_of_lists = list_of_lists
      
    def append_item(self, item, lists):
        item = numpy.array(item)
        if len(item.shape) == 2:
            for i in range(item.shape[0]):
                for j in range(item.shape[1]):
                    lists[i][j].append(item[i][j])
        else:
            for i in range(item.shape[0]):
                lists[i].append(item[i])

    def get_n_events(self, bias_json, sample):
        for name, target in [('us', 'all_upstream'), ('ds', 'all_downstream')]:
            print name, bias_json[sample][target]['weight'],

    def print_corrections(self, test_keys, test_targets, a_file_list):
        print "\nPrinting corrections for files:"
        for a_file in a_file_list:
            print "   ", a_file
        print

        for a_file in a_file_list:
            print "\n"+a_file,
            bias_json = json.loads(open(a_file).read())
            for key in 'all_mc', 'reco':
                print '...', key,
                self.get_n_events(bias_json, key)
            print
            for target in test_targets:
                for key_1, key_2 in test_keys:
                    print bias_json.keys()
                    item = bias_json[key_1][target][key_2]
                    item = numpy.array(item)
                    print key_1, target, key_2
                    if len(item.shape) == 2:
                        print numpy.array(item)[:self.max_bin, :self.max_bin]
                    elif len(item.shape) == 1:
                        print numpy.array(item)[:self.max_bin]
        
    def accumulate_corrections(self, test_keys, test_targets, a_file_list):
        print "\nAccumulating corrections for files:"
        for a_file in a_file_list:
            print "   ", a_file
        print
        for a_file in a_file_list:
            bias_json = json.loads(open(a_file).read())
            if self.list_of_lists == None:
                self.setup_lists(bias_json, test_keys, test_targets)
            for target in test_targets:
                for key_1, key_2 in test_keys:
                    item = bias_json[key_1][target][key_2]
                    lists = self.list_of_lists[key_1][target][key_2]
                    self.append_item(item, lists)
        print self.list_of_lists

    def plot_corrections_1d(self, lists, name, axis):
        canvas = xboa.common.make_root_canvas(name)
        canvas.Divide(*self.canvas_split)
        for i, a_list in enumerate(lists):
            if i >= self.max_bin:
                continue
            a_name = name+" "+str(i)
            print name, a_list
            canvas.cd(i+1)
            hist = xboa.common.make_root_histogram(a_name, a_list, axis, 10)
            hist.SetTitle(name)
            hist.SetStats(True)
            hist.Draw()
        canvas.Update()

    def plot_corrections_2d(self, lists, name, axis):
        canvas = xboa.common.make_root_canvas(name)
        canvas.Divide(*self.canvas_split)
        for i, a_list in enumerate(lists):
            if i >= self.max_bin:
                continue
            canvas.cd(i+1)
            hist = xboa.common.make_root_histogram(name, a_list[i], axis, 10)
            hist.SetTitle(name)
            hist.SetStats(True)
            hist.Draw()

    def plot_corrections(self, test_keys, test_targets):
        for target in test_targets:
            for key_1, key_2 in test_keys:
                lists = self.list_of_lists[key_1][target][key_2]
                name = key_1+" "+target+" "+key_2
                axis = key_2
                if len(numpy.array(lists).shape) == 3:
                    self.plot_corrections_2d(lists, name, axis)
                elif len(numpy.array(lists).shape) == 2:
                    self.plot_corrections_1d(lists, name, axis)


def main():
    a_file_list = glob.glob("output/2017-02-Systematics-6/plots_Simulated_2017-2.7_6-140_lH2_empty_Systematics_*/amplitude/amplitude.json")
    a_file_list = sorted(a_file_list)

    my_hacking = Hacking()
    my_hacking.accumulate_corrections([('crossing_probability', 'migration_matrix'), ('inefficiency','pdf_ratio')],
                      ['all_downstream'],
                      a_file_list)
    my_hacking.plot_corrections([('crossing_probability', 'migration_matrix'), ('inefficiency','pdf_ratio')],
                      ['all_downstream'],)
    raw_input()


def old():
    a_file_list = glob.glob("output/2017-02-Systematics-3/plots_Simulated_2017-2.7_6-140_lH2_empty_Systematics_*/amplitude/amplitude.json")
    a_file_list = sorted(a_file_list)
    Hacking().print_corrections([('crossing_probability', 'migration_matrix'), ('inefficiency','pdf_ratio')],
                      ['all_downstream'],
                      a_file_list)


    a_file_list = glob.glob("output/2017-02-Systematics-5/plots_Simulated_2017-2.7_6-140_lH2_empty_Systematics_tku*/amplitude/amplitude.json")
    a_file_list = sorted(a_file_list)

    Hacking().print_corrections([('crossing_probability', 'migration_matrix'), ('inefficiency','pdf_ratio')],
                      ['all_downstream'],
                      a_file_list)



if __name__ == "__main__":
    main()
