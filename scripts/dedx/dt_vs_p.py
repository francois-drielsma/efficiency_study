import xboa.common


def dt(z, p):
    mass = xboa.common.pdg_pid_to_mass[13]
    energy = (p**2+mass**2)**0.5
    beta = p/energy
    c_light = xboa.common.constants['c_light']
    dt = z/(beta*c_light)
    return dt

def get_dtdp(dt_list, dp):
    deriv_list = []
    for i in range(len(dt_list)-1):
        dt1 = dt_list[i+1]
        dt0 = dt_list[i]
        deriv_list.append(dp/(dt1-dt0))
    return deriv_list

z_tof0 = 5287.2
z_tof1 = 12929.6
z_tof2 = 21114.4

dp = 1
momentum = [float(p) for p in range(100, 201, dp)]
dt_tof01 = [dt(z_tof1-z_tof0, p) for p in momentum]
dt_tof12 = [dt(z_tof2-z_tof1, p) for p in momentum]
d01_list = get_dtdp(dt_tof01, dp)
d12_list = get_dtdp(dt_tof12, dp)
momentum = [p+dp/2. for p in momentum[:-1]]

canvas = xboa.common.make_root_canvas("p vs dt")
hist, graph = xboa.common.make_root_graph("tof01", momentum*2, "#mu p [MeV/c]", d01_list+d12_list, "dp/dt [ns/MeV/c]")
hist01, graph01 = xboa.common.make_root_graph("tof01", momentum, "#mu p [MeV/c]", d01_list, "dp/dt [ns/MeV/c]")
hist12, graph12 = xboa.common.make_root_graph("tof12", momentum, "#mu p [MeV/c]", d12_list, "dp/dt [ns/MeV/c]")

hist.Draw()
graph01.SetLineColor(4)
graph01.SetLineStyle(2)
graph01.Draw("SAMEL")

graph12.SetLineColor(8)
graph12.SetLineStyle(3)
graph12.Draw("SAMEL")

xboa.common.make_root_legend(canvas, [graph01, graph12])

canvas.Update()
for fmt in "png", "eps", "root":
    canvas.Print("plots/dpdt_vs_p."+fmt)

raw_input()


