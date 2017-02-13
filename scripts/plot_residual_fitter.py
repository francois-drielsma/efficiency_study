import json

import xboa.common as common
import config_reco_2016_04_magnet_alignment as config_file

def load_data(fname):
    fin = open(fname)
    data = [json.loads(line) for line in fin.readlines()]
    return data

def plot_alignment_convergence(data):
    x_residual_list = []
    y_residual_list = []
    iteration_list = []
    legend_list = []
    for i, item in enumerate(data):
        iteration_list.append(i)
        x_residual_list.append(item["score"]["x"])
        y_residual_list.append(item["score"]["y"])
    canvas = common.make_root_canvas("alignment convergence")
    hist, zero_graph = common.make_root_graph("x residual",
                                             [-iteration_list[-1], 2*iteration_list[-1]], "Iteration",
                                             [0., 0.], "Residual [mm]")
    hist, x_graph = common.make_root_graph("x residual", iteration_list, "Iteration", x_residual_list, "Residual [mm]")
    x_graph.SetMarkerStyle(26)
    x_graph.SetLineColor(10)
    x_graph.SetMarkerColor(4)
    canvas.Draw()
    hist.Draw()
    zero_graph.Draw("L")
    x_graph.Draw("P")
    
    hist, y_graph = common.make_root_graph("y residual", iteration_list, "Iteration", y_residual_list, "Residual [mm]")
    y_graph.SetLineColor(10)
    y_graph.SetMarkerStyle(26)
    y_graph.SetMarkerColor(8)
    y_graph.Draw("SAMEP")
    leg = common.make_root_legend(canvas, [x_graph, y_graph])
    leg.Draw()
    for format in ["pdf", "png", "root"]:
        canvas.Print("alignment_convergence."+format)

    
def do_plots(config):
    for anal in config.analyses:
        data = load_data(anal["plot_dir"]+"/alignment.json")
        plot_alignment_convergence(data)
    
def main():
    do_plots(config_file.Config())

if __name__ == "__main__":
    main()    
    raw_input()

