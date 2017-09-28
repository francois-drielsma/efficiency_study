import json

import ROOT
import xboa.common as common
#import config_reco_2016_04_magnet_alignment as config_file

ROOT_OBJECTS = []

def load_data(fname):
    fin = open(fname)
    data = [json.loads(line) for line in fin.readlines()]
    return data

def opt_key(optimisation_element):
    return "_".join(optimisation_element["axes"].keys())

def plot_alignment_convergence(data, plot_dir):
    del data[0]
    optimisation_set = set([opt_key(item) for item in data])
    for opt in optimisation_set:
        name = opt.replace("Coils", "")
        results = [item for item in data if opt_key(item) == opt]
        canvas = common.make_root_canvas("alignment convergence - "+opt)
        x_err_graph = ROOT.TGraphErrors(len(results))
        x_err_graph.SetMarkerColor(4)
        x_err_graph.SetMarkerStyle(24)
        y_err_graph = ROOT.TGraphErrors(len(results))
        y_err_graph.SetMarkerColor(8)
        y_err_graph.SetMarkerStyle(24)
        ROOT_OBJECTS.append(x_err_graph)
        ROOT_OBJECTS.append(y_err_graph)
        ymin, ymax = 0., 0.
        for j in range(4):
            for i, item in enumerate(results):
                result = item["score"]["results"][j]
                p_bin = result['p_bin']
                n = result["n"]
                x_err_graph.SetPoint(i, i, result["x"])
                x_err_graph.SetPointError(i, 0., result["x_err"]/n**0.5)
                y_err_graph.SetPoint(i, i+0.3, result["y"])
                y_err_graph.SetPointError(i, 0., result["y_err"]/n**0.5)
                ymin = min(ymin, result["x"]-2*result["x_err"]/n**0.5)
                ymax = max(ymax, result["x"]+2*result["x_err"]/n**0.5)
                ymin = min(ymin, result["y"]-2*result["y_err"]/n**0.5)
                ymax = max(ymax, result["y"]+2*result["y_err"]/n**0.5)
            hist, zero_graph = common.make_root_graph("x residual",
                                                     [0, len(results)], "Iteration",
                                                     [0., 0.], "Residual [mm]",
                                                      ymin=ymin, ymax=ymax)
            hist.SetTitle(name+" p: "+str(p_bin[0])+" to "+str(p_bin[1])+" MeV/c")
            hist.Draw()
            zero_graph.Draw("SAMEL")
            x_err_graph.Draw("SAMEP")
            y_err_graph.Draw("SAMEP")
            for format in ["pdf", "png", "root"]:
                canvas.Print(plot_dir+"/alignment_convergence_"+name+"_"+str(p_bin[0])+"_"+str(p_bin[1])+"."+format)


def plot_alignment(data, plot_dir):
    del data[0]
    optimisation_set = set([opt_key(item) for item in data])
    for opt in optimisation_set:
        name = opt.replace("Coils", "")
        results = [item for item in data if opt_key(item) == opt]
        canvas = common.make_root_canvas("alignment - "+name)
        n_list = range(len(results))
        x_list = [item["axes"][opt]["x"]["value"] for item in results]
        y_list = [item["axes"][opt]["y"]["value"] for item in results]
        xp_list = [item["axes"][opt]["xp"]["value"]*1000. for item in results]
        yp_list = [item["axes"][opt]["yp"]["value"]*1000. for item in results]

        hist, zero_graph = common.make_root_graph("x residual",
                                                 [0, len(results)], "Iteration",
                                                 [0., 0.], "Residual [mm]")
        hist, x_graph = common.make_root_graph("x alignment", n_list, "iteration", x_list, "misalignment [mm]")
        hist, y_graph = common.make_root_graph("y alignment", n_list, "iteration", y_list, "misalignment [mm]")
        hist, all_graph = common.make_root_graph("all alignment", n_list*2, "iteration", x_list+y_list, "misalignment [mm]")
        hist.SetTitle(name+" alignment")
        hist.Draw()
        zero_graph.Draw("SAMEL")
        x_graph.Draw("SAMEP")
        x_graph.SetMarkerStyle(24)
        x_graph.SetMarkerColor(4)
        y_graph.Draw("SAMEP")
        y_graph.SetMarkerStyle(24)
        y_graph.SetMarkerColor(8)
        canvas.Update()
        for format in ["pdf", "png", "root"]:
            canvas.Print(plot_dir+"/alignment_"+name+"."+format)


        canvas = common.make_root_canvas("tilt - "+name)
        hist, zero_graph = common.make_root_graph("x residual",
                                                 [0, len(results)], "Iteration",
                                                 [0., 0.], "Residual [mm]")
        hist, xp_graph = common.make_root_graph("x alignment", n_list, "iteration", xp_list, "tilt [mrad]")
        hist, yp_graph = common.make_root_graph("y alignment", n_list, "iteration", yp_list, "tilt [mrad]")
        hist, all_graph = common.make_root_graph("all alignment", n_list*2, "iteration", xp_list+yp_list, "titl [mrad]")
        hist.SetTitle(name+" tilt")
        hist.Draw()
        zero_graph.Draw("SAMEL")
        xp_graph.Draw("SAMEP")
        xp_graph.SetMarkerStyle(24)
        xp_graph.SetMarkerColor(4)
        yp_graph.Draw("SAMEP")
        yp_graph.SetMarkerStyle(24)
        yp_graph.SetMarkerColor(8)

        canvas.Update()
        for format in ["pdf", "png", "root"]:
            canvas.Print(plot_dir+"/tilt_"+name+"."+format)


    
def do_plots(config):
    for anal in config.analyses:
        data = load_data(anal["plot_dir"]+"/alignment.json")
        plot_alignment(data, anal["plot_dir"])
        plot_alignment_convergence(data, anal["plot_dir"])
    
def main():
    config = config_file.Config()
    config.analyses[0]["plot_dir"] = "iteration_4_magnet-alignment_first-mc/2016-04_1.2_magnet_alignment_pry_v3/plots_Mixed_momenta,_MAUS-v2.6.5/"
    do_plots(config_file.Config())
    

if __name__ == "__main__":
    main()    
    raw_input()

