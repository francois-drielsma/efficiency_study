import xboa

class AnalysisBase(object):
    def __init__(self, config, config_anal, data_loader):
        self.data_loader = data_loader
        self.config = config
        self.config_anal = config_anal
        self.plot_dir = config_anal["plot_dir"]
        self.plots = {}
        self.process_list = {}
        self.death_list = {}

    def get_plot(self, name):
        if name not in self.plots:
            new_plot = {}
            new_plot["canvas"] = xboa.common.make_root_canvas(name)
            new_plot["histograms"] = {}
            new_plot["graphs"] = {}
            new_plot["misc"] = {}
            self.plots[name] = new_plot
            print "Adding canvas", name
        self.plots[name]["canvas"].cd()
        return self.plots[name]

    def make_root_histogram(self, name, *args):
        my_plot = self.get_plot(name)
        hist = xboa.common.make_root_histogram(*args)
        my_plot["histograms"][args[0]] = hist
        return hist

    def make_root_graph(self, name, *args):
        my_plot = self.get_plot(name)
        hist, graph = xboa.common.make_root_graph(*args)
        my_plot["graphs"][args[0]] = graph
        my_plot["histograms"][args[0]] = hist
        return hist, graph

    def print_plots(self):
        for name, my_plot in self.plots.iteritems():
            print "Printing canvas", name
            my_plot["canvas"].cd()
            my_plot["canvas"].Draw()
            my_plot["canvas"].Update()
            plot_title = name.replace(" ", "_")
            for format in ["png", "root", "eps"]:
                my_plot["canvas"].Print(self.plot_dir+"/"+plot_title+"."+format)

