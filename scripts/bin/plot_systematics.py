import os
import shutil
import glob
import json
import numpy
import xboa.common
import ROOT
numpy.set_printoptions( linewidth=1000, precision=4)
import utilities.root_style

class Hacking(object):
    def __init__(self, output_dir):
        self.list_of_lists = None
        self.names = []
        self.max_bin = 12
        self.canvas_split = (3, 4)
        self.output_dir = output_dir+"/"
        self.run_configuration = ""
        self.stats_errors = {} # fractional stats error
    root_objects = []

    def clean_output_dir(self):
        try:
            shutil.rmtree(self.output_dir)
        except OSError:
            pass
        os.makedirs(self.output_dir)

    def get_run_configuration(self, file_name):
        self.run_configuration = file_name.split("_")[3]

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

    def get_names(self, a_file_list):
        self.names = []
        for a_file in a_file_list:
            name = a_file.split("/amplitude")[0]
            name = name.split("_Systematics_")[1]
            self.names.append(name)
            print name.ljust(12), a_file

    def accumulate_corrections(self, test_keys, test_targets, a_file_list):
        if len(a_file_list) == 0:
            raise KeyError("No files found")
        self.list_of_lists = None
        self.get_names(a_file_list)
        self.get_run_configuration(a_file_list[0])
        for i, a_file in enumerate(a_file_list):
            bias_json = json.loads(open(a_file).read())
            if i == 0:
                self.setup_lists(bias_json, test_keys, test_targets)
            for target in test_targets:
                for key_1, key_2 in test_keys:
                    item = bias_json[key_1][target][key_2]
                    lists = self.list_of_lists[key_1][target][key_2]
                    self.append_item(item, lists)
    
    def get_lists(self, key_1, target, key_2):
        name = key_1+" "+target+" "+key_2
        lists = self.list_of_lists[key_1]
        lists = lists[target]
        lists = lists[key_2]
        return name, lists

    def calculate_uncertainty(self, test_keys, test_targets):
        for target in test_targets:
            for key_1, key_2 in test_keys:
                name, lists = self.get_lists(key_1, target, key_2)
                # fill the map tree with stats error
                self.stats_errors[name] = [None]*min(self.max_bin+1, len(lists))
                if len(numpy.array(lists).shape) == 3:
                    lists = [lists[i][i] for i, a_list in enumerate(lists)]
                for i, a_list in enumerate(lists):
                    if i <= self.max_bin:
                        self.stats_errors[name][i] = numpy.std(a_list)/(len(a_list)-1.)**0.5
                print "Uncertainties:", name, self.stats_errors[name]

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
        graph_list = [ROOT.TGraphErrors(self.max_bin) for i in range(len(lists[0]))]
        x_min, x_max, y_min, y_max = 0., None, None, None
        for i, a_list in enumerate(lists):
            if i > self.max_bin:
                continue
            err = self.stats_errors[name][i]*(2.**0.5) #factor of sqrt(2) as we are looking at difference
            print i, str(self.get_bin_centre(i)).ljust(6), str(round(err, 1)).ljust(6),
            for j, value in enumerate(a_list):
                value = (value - a_list[0])
                graph_list[j].SetPoint(i, self.get_bin_centre(i), value)
                graph_list[j].SetPointError(i, 0., err)
                print str(round(value, 1)).ljust(6),
                y_min = max(-value+err, y_min) # use max(-y) to handle initial None
                y_max = max(value+err, y_max)
                x_max = max(self.get_bin_centre(i), x_max)
            print str(round(y_min, 2)).ljust(6)
        y_min *= -1.
        y_min -= (y_max-y_min)*0.1
        y_max += (y_max-y_min)*0.1
        x_max *= 1.5
        if name in self.multigraph_axis_range:
            [y_min, y_max] = self.multigraph_axis_range[name]
 
        draw_option = "SAME P L"
        a_type = ""
        if "migration_matrix" in name:
            a_type = "Matrix "
        elif "inefficiency" in name:
            a_type = "Efficiency "
        hist = ROOT.TH2D("", ";Amplitude [mm];Change in "+a_type+"Correction", 1000, x_min, x_max, 1000, y_min, y_max)
        hist.SetStats(False)
        hist.SetTitle(self.run_configuration)
        hist.Draw()
        self.root_objects.append(hist)
        legend = ROOT.TLegend(0.65, 0.5, 0.85, 0.9)
        self.root_objects.append(legend)
        for i, graph in enumerate(graph_list):
            if self.names[i] == "tku_base":
                continue
            graph.SetName(self.names[i])
            graph.SetTitle()
            graph.SetMarkerStyle(26+i)
            color = i/(len(graph_list)-1.)*ROOT.gStyle.GetNumberOfColors()
            color = ROOT.gStyle.GetColorPalette(int(color))
            graph.SetMarkerColor(color)
            graph.SetLineColor(color)
            graph.SetFillColor(10)
            graph.SetLineStyle(i)
            graph.Draw(draw_option)
            draw_option = "SAME P L"
            legend.AddEntry(graph, self.name_map[self.names[i]], "P L")
        legend.Draw()

        #legend = canvas.BuildLegend(0.65, 0.5, 0.85, 0.9)
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
            graph.SetPointError(i, 0, self.stats_errors[name][i])
        if name in self.graph_axis_range:
            [ymin, ymax] = self.graph_axis_range[name]
            graph.GetYaxis().SetRangeUser(ymin, ymax)
        if "migration_matrix" in name:
            axis = "Probability (true bin) = (recon bin)"
        elif "inefficiency" in name:
            axis = "Efficiency Correction"
        graph.GetYaxis().SetTitle(axis)
        graph.GetXaxis().SetTitle("Amplitude [mm]")
        graph.SetMarkerStyle(21)
        graph.SetTitle(self.run_configuration)
        graph.Draw("AP")
        self.root_objects.append(graph)
        plot_name = name.replace(" ", "_")
        canvas.Update()
        for fmt in ['.png', '.root', '.pdf']:
            canvas.Print(self.output_dir+plot_name+"_graph"+fmt)


    def plot_corrections(self, test_keys, test_targets, plot_routine_str):
        for target in test_targets:
            for key_1, key_2 in test_keys:
                name, lists = self.get_lists(key_1, target, key_2)
                axis = key_2
                if len(numpy.array(lists).shape) == 3:
                    lists = [lists[i][i] for i, a_list in enumerate(lists)]
                plot_routine = {
                    "hist":self.plot_corrections_hist,
                    "graph":self.plot_corrections_graph,
                    "multigraph":self.plot_corrections_multigraph,
                }[plot_routine_str]
                plot_routine(lists, name, axis)

    name_map = {
        "tku_density_plus":"TKU Density",
        "tku_pos_plus":"TKU Position",
        "tku_rot_plus":"TKU Rotation",
        "tku_scale_C_plus":"TKU Centre Coil",
        "tku_scale_E1_plus":"TKU End1 Coil",
        "tku_scale_E2_plus":"TKU End2 Coil",
        "tkd_density_plus":"TKD Density",
        "tkd_pos_plus":"TKD Position",
        "tkd_rot_plus":"TKD Rotation",
        "tkd_scale_C_plus":"TKD Centre Coil",
        "tkd_scale_E1_plus":"TKD End1 Coil",
        "tkd_scale_E2_plus":"TKD End2 Coil",
    }

    graph_axis_range = {
        "crossing_probability all_downstream migration_matrix":[0.6, 1.0],
        "crossing_probability all_upstream migration_matrix":[0.6, 1.0],
        "inefficiency all_downstream pdf_ratio":[1.0, 1.4],
    }


    multigraph_axis_range = {
        "crossing_probability all_downstream migration_matrix":[-0.12, 0.12],
        "crossing_probability all_upstream migration_matrix":[-0.12, 0.12],
        "inefficiency all_downstream pdf_ratio":[-0.12, 0.12],
    }


def file_list(src_dir, emittance, absorber, suffix):
    name = "from_scarf/"+src_dir+"/plots_Simulated_2017-2.7_"
    name += str(emittance)+"-140_"+absorber+"_Systematics_"+suffix
    name += "/amplitude/amplitude.json"
    name_list = sorted(glob.glob(name))
    if len(name_list) == 0:
        print name
        raise KeyError("Failed to find names for "+str(name))
    return name_list

def do_upstream(input_dir, emittance_list, output_dir):
    for emittance in emittance_list:
        a_file_list = file_list(input_dir, emittance, "lH2_empty", "?")+\
                      file_list(input_dir, emittance, "lH2_empty", "??")
        a_file_list = sorted(a_file_list)
        my_hacking = Hacking(output_dir+emittance+"-upstream")
        my_hacking.clean_output_dir()
        my_hacking.accumulate_corrections([('crossing_probability', 'migration_matrix')],
                          ['all_upstream'],
                          a_file_list)
        my_hacking.calculate_uncertainty([('crossing_probability', 'migration_matrix')],
                          ['all_upstream'])
        my_hacking.plot_corrections([('crossing_probability', 'migration_matrix')],
                          ['all_upstream'], "graph")

        a_file_list = file_list(input_dir, emittance, "lH2_empty", "tku_base")+\
                      file_list(input_dir, emittance, "lH2_empty", "tku_*")
        my_hacking.accumulate_corrections([('crossing_probability', 'migration_matrix')],
                          ['all_upstream'],
                          a_file_list)
        my_hacking.plot_corrections([('crossing_probability', 'migration_matrix')],
                          ['all_upstream'], "multigraph")



def do_downstream(input_dir, emittance_list, output_dir):
    for emittance in emittance_list:
        a_file_list = file_list(input_dir, emittance, "lH2_empty", "?")+\
                      file_list(input_dir, emittance, "lH2_empty", "??")
        a_file_list = sorted(a_file_list)
        my_hacking = Hacking(output_dir+emittance+"-downstream")
        my_hacking.clean_output_dir()
        my_hacking.accumulate_corrections([('crossing_probability', 'migration_matrix'), ('inefficiency','pdf_ratio')],
                          ['all_downstream'],
                          a_file_list)
        my_hacking.calculate_uncertainty([('crossing_probability', 'migration_matrix'), ('inefficiency','pdf_ratio')],
                          ['all_downstream'])
        my_hacking.plot_corrections([('crossing_probability', 'migration_matrix'), ('inefficiency','pdf_ratio')],
                          ['all_downstream'], "graph")

        a_file_list = file_list(input_dir, emittance, "lH2_empty", "tku_base")+\
                      file_list(input_dir, emittance, "lH2_empty", "tkd_*")
        my_hacking.accumulate_corrections([('crossing_probability', 'migration_matrix'), ('inefficiency','pdf_ratio')],
                          ['all_downstream'],
                          a_file_list)
        my_hacking.plot_corrections([('crossing_probability', 'migration_matrix'), ('inefficiency','pdf_ratio')],
                          ['all_downstream'], "multigraph")

def do_copy(input_dir, emittance_list, output_dir):
    for emittance in emittance_list:
        a_file_list = file_list(input_dir, emittance, "lH2_empty", "tku_base")
        plot_root = os.path.split(a_file_list[0])[0]
        for stream in "upstream", "downstream":
            plot_files = glob.glob(plot_root+"/crossing_probability_"+stream+".*")
            for src_file in plot_files:
                target_file = output_dir+"/"+emittance+"-"+stream+"/"+os.path.split(src_file)[1]
                print "Copying", src_file, "to", target_file
                shutil.copy(src_file, target_file)


def copy_more(input_dir, output_dir):
    a_dir_list = glob.glob("from_scarf/"+input_dir+"/*tku_base/amplitude/weighting")
    for a_dir in a_dir_list:
        print a_dir
        my_bl = a_dir.split("2017-2.7_")[1]
        my_bl = my_bl.split("_lH2_empty")[0]
        print my_bl
        target_dir = output_dir+"/"+my_bl+"_efficiency"
        try:
            shutil.rmtree(target_dir)
        except OSError:
            pass
        shutil.copytree(a_dir, target_dir)

def main():
    utilities.root_style.setup_gstyle()
    output_dir = "from_scarf/2017-02-Systematics/systematics_summary/"
    
    #do_copy("2017-02-Systematics", ["3", "4", "6", "10"], output_dir)
    copy_more("2017-02-Systematics", output_dir)
    #do_upstream("2017-02-Systematics", ["3", "4", "6", "10"], output_dir)
    #do_downstream("2017-02-Systematics", ["3", "4", "6", "10"], output_dir)

    return



if __name__ == "__main__":
    main()
