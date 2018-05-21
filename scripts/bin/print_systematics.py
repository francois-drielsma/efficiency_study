import os
import shutil
import glob
import json
import numpy
import xboa.common
import ROOT
numpy.set_printoptions( linewidth=1000, precision=4)


class Hacking(object):
    def __init__(self):
        self.bias_json = None
        self.list_of_lists = None
        self.max_bin = 12
        self.canvas_split = (3, 4)
        self.output_dir = "output/2017-02-Systematics/systematics_summary/"
    root_objects = []

    def clean_output_dir(self):
        try:
            shutil.rmtree(self.output_dir)
        except OSError:
            pass
        os.makedirs(self.output_dir)

    def get_bin_centre(self, i):
        return (i+0.5)*5.

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

    def names(self, a_file_list):
        self.names = []
        for a_file in a_file_list:
            name = a_file.split("/amplitude")[0]
            name = name.split("_Systematics_")[1]
            self.names.append(name)
            print name.ljust(12), a_file

    def accumulate_corrections(self, test_keys, test_targets, a_file_list):
        self.names(a_file_list)
        for a_file in a_file_list:
            bias_json = json.loads(open(a_file).read())
            if self.list_of_lists == None:
                self.setup_lists(bias_json, test_keys, test_targets)
            for target in test_targets:
                for key_1, key_2 in test_keys:
                    item = bias_json[key_1][target][key_2]
                    lists = self.list_of_lists[key_1][target][key_2]
                    self.append_item(item, lists)

    def plot_corrections_hist(self, lists, name, axis):
        canvas = xboa.common.make_root_canvas(name)
        canvas.Divide(*self.canvas_split)
        for i, a_list in enumerate(lists):
            if i >= self.max_bin:
                continue
            a_name = name+" "+str(i)
            canvas.cd(i+1)
            hist = xboa.common.make_root_histogram(a_name, a_list, axis, 10)
            hist.SetTitle(name)
            hist.SetStats(True)
            hist.Draw()
        canvas.Update()
        plot_name = name.replace(" ", "_")
        for fmt in ['.png', '.root', '.pdf']:
            canvas.Print(self.output_dir+plot_name+"_hist"+fmt)

    def plot_corrections_multigraph(self, lists, name, axis):
        canvas = xboa.common.make_root_canvas(name)
        graph_list = [ROOT.TGraph(self.max_bin) for i in range(len(lists[0]))]
        print name, axis, "multigraph"
        y_min, y_max = None, None
        for i, a_list in enumerate(lists):
            if i > self.max_bin:
                continue
            for j, value in enumerate(a_list):
                value = value - a_list[0]
                graph_list[j].SetPoint(i, self.get_bin_centre(i), value)
                y_min = max(-value, y_min)
                y_max = max(value, y_max)
        y_min *= -1.
        y_min -= (y_max-y_min)*0.1
        y_max += (y_max-y_min)*0.1
        draw_option = "A P L"
        for i, graph in enumerate(graph_list):
            if self.names[i] == "tku_base":
                continue
            graph.GetYaxis().SetRangeUser(y_min, y_max)
            graph.SetName(self.names[i])
            graph.SetTitle(self.names[i])
            graph.SetMarkerStyle(26)
            color = i/(len(graph_list)-1.)*ROOT.gStyle.GetNumberOfColors()
            color = ROOT.gStyle.GetColorPalette(int(color))
            graph.SetMarkerColor(color)
            graph.SetLineColor(color)
            graph.SetFillColor(10)
            graph.SetLineStyle(i)
            graph.Draw(draw_option)
            draw_option = "SAME P L"
        legend = canvas.BuildLegend(0.2, 0.5, 0.5, 0.9)
        self.root_objects.append(graph_list)
        plot_name = name.replace(" ", "_")
        canvas.Update()
        for fmt in ['.png', '.root', '.pdf']:
            canvas.Print(self.output_dir+plot_name+"_multigraph"+fmt)

    def plot_corrections_graph(self, lists, name, axis):
        canvas = xboa.common.make_root_canvas(name)
        graph = ROOT.TGraphErrors(self.max_bin)
        print name, axis, "graph"
        for i, a_list in enumerate(lists):
            if i > self.max_bin:
                continue
            mean, std = numpy.mean(a_list), numpy.std(a_list)
            print "   ", i, mean, std
            graph.SetPoint(i, self.get_bin_centre(i), mean)
            graph.SetPointError(i, 0, std/len(a_list)**0.5)
        graph.SetMarkerStyle(21)
        graph.SetTitle(name)
        graph.Draw("AP")
        self.root_objects.append(graph)
        plot_name = name.replace(" ", "_")
        canvas.Update()
        for fmt in ['.png', '.root', '.pdf']:
            canvas.Print(self.output_dir+plot_name+"_graph"+fmt)


    def plot_corrections(self, test_keys, test_targets, plot_routine_str):
        for target in test_targets:
            for key_1, key_2 in test_keys:
                lists = self.list_of_lists[key_1][target][key_2]
                name = key_1+" "+target+" "+key_2
                axis = key_2
                if len(numpy.array(lists).shape) == 3:
                    lists = [a_list[i] for i, a_list in enumerate(lists)]
                plot_routine = {
                    "hist":self.plot_corrections_hist,
                    "graph":self.plot_corrections_graph,
                    "multigraph":self.plot_corrections_multigraph,
                }[plot_routine_str]
                plot_routine(lists, name, axis)

def main():
    a_file_list = glob.glob("output/2017-02-Systematics/plots_Simulated_2017-2.7_10-140_lH2_empty_Systematics_?/amplitude/amplitude.json")
    a_file_list = sorted(a_file_list)

    my_hacking = Hacking()
    my_hacking.clean_output_dir()
    my_hacking.accumulate_corrections([('crossing_probability', 'migration_matrix')],
                      ['all_upstream'],
                      a_file_list)
    my_hacking.plot_corrections([('crossing_probability', 'migration_matrix')],
                      ['all_upstream'], "graph")

    my_hacking = Hacking()
    my_hacking.accumulate_corrections([('crossing_probability', 'migration_matrix'), ('inefficiency','pdf_ratio')],
                      ['all_downstream'],
                      a_file_list)
    my_hacking.plot_corrections([('crossing_probability', 'migration_matrix'), ('inefficiency','pdf_ratio')],
                      ['all_downstream'], "graph")


    a_file_list = ["output/2017-02-Systematics/plots_Simulated_2017-2.7_10-140_lH2_empty_Systematics_tku_base/amplitude/amplitude.json"]
    a_file_list += glob.glob("output/2017-02-Systematics/plots_Simulated_2017-2.7_10-140_lH2_empty_Systematics_tku*/amplitude/amplitude.json")
    my_hacking = Hacking()
    my_hacking.accumulate_corrections([('crossing_probability', 'migration_matrix'), ('inefficiency','pdf_ratio')],
                      ['all_upstream'],
                      a_file_list)
    my_hacking.plot_corrections([('crossing_probability', 'migration_matrix'), ('inefficiency','pdf_ratio')],
                      ['all_upstream'], "multigraph")

    a_file_list = ["output/2017-02-Systematics/plots_Simulated_2017-2.7_10-140_lH2_empty_Systematics_tku_base/amplitude/amplitude.json"]
    a_file_list += glob.glob("output/2017-02-Systematics/plots_Simulated_2017-2.7_10-140_lH2_empty_Systematics_tkd*/amplitude/amplitude.json")
    my_hacking = Hacking()
    my_hacking.accumulate_corrections([('crossing_probability', 'migration_matrix'), ('inefficiency','pdf_ratio')],
                      ['all_downstream'],
                      a_file_list)
    my_hacking.plot_corrections([('crossing_probability', 'migration_matrix'), ('inefficiency','pdf_ratio')],
                      ['all_downstream'], "multigraph")




if __name__ == "__main__":
    main()
