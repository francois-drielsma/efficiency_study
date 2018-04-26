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

    def do_mice_logo(self, do_isis):
        print "DOING MICE LOGO"
        y_min=0.895
        if do_isis:
            y_min = 0.8
        text_box = ROOT.TPaveText(0.2, y_min, 0.99, 0.99, "NDC")
        text_box.SetTextSize(0.05)
        text_box.SetFillColor(0)
        text_box.SetBorderSize(0)
        text_box.SetTextAlign(33)
        text_box.AddText("MICE Preliminary")
        if do_isis:
            text_box.AddText("ISIS User Runs 2017/02 and 2017/03")
        text_box.Draw()
        self.root_objects.append(text_box)
        
    def draw_cong(self, pad, cong, i, j):
        pad.cd()
        legend_size = None
        if cong.legend != None:
            legend_size = cong.options["legend"]["pos"]
        if j != 0:
            for hist in cong.hist_list:
                hist.GetYaxis().SetLabelOffset(100)
                hist.GetYaxis().SetRangeUser(0, 100000)

        cong.label_size = 0.1
        cong.murgle_many(pad, cong.hist_list, cong.graph_list)

        #for graph in cong.graph_list:
        #    if graph.GetN() > 2:
        #        graph.Draw("LP")
        #    else:
        #        graph.Draw("SAME L")
        if i == 0 and j == 2:
            self.do_legend(i, j, cong.legend, legend_size, pad)

    def extra_labels(self, cong_list):
        extra_labels = cong_list[0].conglomerations[0].options["extra_labels"]
        left_m = 0.128
        right_m = 0.12
        x_step = (1.0-left_m-right_m)/self.cols
        for i, item in enumerate(extra_labels["top"]):
            text_box = ROOT.TPaveText(left_m+x_step*i, 0.87, left_m+x_step*(i+1), 0.95, "NDC")
            text_box.SetTextSize(0.04)
            text_box.SetFillColor(0)
            text_box.SetBorderSize(0)
            print "Setting text", item
            text_box.AddText(item)
            text_box.Draw()
            self.root_objects.append(text_box)
        top_m = 0.14
        bottom_m = 0.14
        y_step = (1.0-top_m-bottom_m)/self.rows
        for i, item in enumerate(extra_labels["right"]):
            text_box = ROOT.TPaveText(1.0-right_m, 1.0-top_m-y_step*(i+1),
                                      0.99, 1.0-top_m-y_step*(i), "NDC")
            text_box.SetTextSize(0.04)
            text_box.SetTextAlign(12)
            text_box.SetFillColor(0)
            text_box.SetBorderSize(0)
            print "Setting text", item
            text_box.AddText(item)
            text_box.Draw()
            self.root_objects.append(text_box)
 
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
        hist_index = 0
        merge_canvas.Divide(1, 1, 0.1, 0.1)
        pad = merge_canvas.cd(1)
        pad.Divide(self.cols, self.rows, 0., 0.)
        for i in range(self.rows):
            for j in range(self.cols):
                cong_base = self.conglomerate_list[hist_index]
                for cong in cong_base.conglomerations:
                    if cong.options["canvas_name"] == canvas_name:
                        break
                source_x_axis = cong.x_axis
                source_y_axis = cong.y_axis
                target_dir = cong.config.target_plot_dir
                self.draw_cong(pad.GetPad(hist_index+1), cong, i, j)
                hist_index += 1
                self.do_mice_logo(i == 0 and j == 2)

        merge_canvas.cd()
        if source_x_axis:
            source_x_axis.Draw()
        if source_y_axis:
            source_y_axis.Draw()
        self.extra_labels(self.conglomerate_list)
        merge_canvas.Update()
        for fmt in ["png", "root", "eps", "pdf"]:
            merge_canvas.Print(target_dir+"/"+canvas_name+"."+fmt)
        self.merge_list.append(merge_dict)
        
    root_objects = []