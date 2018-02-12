import xboa.common
#E     mean dEdx    most prob dEdx   ref dEdx

fin = open("scripts/dedx/dedx_data")
fin.readline()
ref_de, mean_de, mpb_de, energy = [], [], [], []
for line in fin.readlines():
    words = line.split()
    energy.append(float(words[0])+xboa.common.pdg_pid_to_mass[13])
    mean_de.append(float(words[1]))
    mpb_de.append(float(words[2]))
    ref_de.append(float(words[3]))

canvas = xboa.common.make_root_canvas("e loss")

hist, ref_graph = xboa.common.make_root_graph("PDG Energy Loss", energy, "Energy [MeV]", ref_de, "dE/dx [MeV/mm]", ymax=-0.025, ymin=-0.05)
hist.Draw()
ref_graph.Draw("SAMEC")

hist, mean_graph = xboa.common.make_root_graph("Mean Energy Loss", energy, "Energy [MeV]", mean_de, "dE/dx [MeV/mm]")
mean_graph.SetLineColor(4)
mean_graph.Draw("C")

hist, mpb_graph = xboa.common.make_root_graph("Most Prob Energy Loss", energy, "Energy [MeV]", mpb_de, "dE/dx [MeV/mm]")
mpb_graph.SetLineColor(8)
mpb_graph.Draw("C")

xboa.common.make_root_legend(canvas, [ref_graph, mean_graph, mpb_graph])


for format in "eps", "root", "png":
    canvas.Print("plots/dEdx."+format)
    
raw_input()