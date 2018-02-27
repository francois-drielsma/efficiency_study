import ROOT
import utilities

class ConglomerateMerge(object):
    """
    Merge individual plots of the same data into one uber-canvas

    Each plot is intended to have the same axes, so we can suppress axes to
    save white space
    """
    def __init__(self, conglomerate_list):
        """
        conglomerate_list: list of ConglomerateBase (each of which is a merge
        of MC and data that pertains to several plots for a particular beam 
        setting)
        """
        self.conglomerate_list = conglomerate_list
        self.merge_list = []
        self.rows = 0
        self.cols = 0

    def merge_all(self, rows, cols):
        self.rows = rows
        self.cols = cols
        canvas_names = set()
        for cong_base in self.conglomerate_list:
            for cong in cong_base.conglomerations:
                canvas_names.add(cong.options["canvas_name"])
        for name in canvas_names:
            self.merge_one(name)

    def do_legend(self, i, j, legend, legend_size, pad):
        if legend == None:
            return
        # weird axis bounds!
        if j == 0:
            # x goes from 0.0 to 1.0/0.9?
            legend.SetX1(legend_size[0]*1.0/0.9)
            legend.SetX2(legend_size[2]*1.0/0.9)
        else:
            # x goes from -0.1 to 1.0/0.9?
            legend.SetX1((legend_size[0]-0.1)*1.0/0.9)
            legend.SetX2((legend_size[2])*1.0/0.9)
        if i == self.rows-1:
            # y goes from 0.0 to 1.0/0.9
            legend.SetY1(legend_size[1]*1.0/0.9)
            legend.SetY2(legend_size[3]*1.0/0.9)
        else:
            # y goes from -0.1 to 1.0/0.9
            legend.SetY1(legend_size[1]*1.0/0.9)
            legend.SetY2(legend_size[3]*1.0/0.9)
        legend.Draw()

    def draw_cong(self, pad, cong, i, j):
        pad.cd()
        legend_size = None
        if cong.legend != None:
            legend_size = cong.options["legend"]["pos"]
        cong.redraw(pad, cong.hist_list, cong.graph_list)

        #for graph in cong.graph_list:
        #    if graph.GetN() > 2:
        #        graph.Draw("LP")
        #    else:
        #        graph.Draw("SAME L")
        self.do_legend(i, j, cong.legend, legend_size, pad)

    def merge_one(self, canvas_name):
        """
        Merge one plot
        
        Draw data from all elements in self.conglomerate_list; find the "ith"
        plot and merge into one uber-canvas
        """
        print "Merging", canvas_name
        merge_dict = {}
        source_hist_list = []
        source_graph_list = []
        source_legend_list = []

        merge_canvas = ROOT.TCanvas(canvas_name+"_c1", canvas_name+"_c1", 1400, 1000)
        merge_dict["canvas"] = merge_canvas
        merge_canvas.Draw()
        merge_canvas.Divide(1, 1, 0., 0.)
        merge_pad = merge_canvas.GetPad(1)
        y_min, x_min = 0.1, 0.1
        merge_pad.SetPad(x_min, y_min, 0.99, 0.99)
        hist_index = 0
        merge_pad.Divide(self.cols, self.rows, 0., 0.)
        for i in range(self.rows):
            for j in range(self.cols):
                cong_base = self.conglomerate_list[hist_index]
                for cong in cong_base.conglomerations:
                    if cong.options["canvas_name"] == canvas_name:
                        break
                source_x_axis = cong.x_axis
                source_y_axis = cong.y_axis
                target_dir = cong.config.target_plot_dir
                self.draw_cong(merge_pad.GetPad(hist_index+1), cong, i, j)
                hist_index += 1

        merge_canvas.cd()
        if source_x_axis:
            source_x_axis.Draw()
        if source_y_axis:
            source_y_axis.Draw()
        merge_canvas.Update()
        for fmt in ["png", "root", "eps", "pdf"]:
            merge_canvas.Print(target_dir+"/"+canvas_name+"."+fmt)
        self.merge_list.append(merge_dict)