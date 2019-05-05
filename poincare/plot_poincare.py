from array import array
import ROOT
import InfoBox

# Phase space variables names and labels
psv_names = ["x", "px", "y", "py"]
psv_labels = [
    "x [mm]",
    "p_{x}  [MeV/c]",
    "y [mm]",
    "p_{y}  [MeV/c]",
]
nvars = len(psv_names)

# Phase space and amplitude dictionaries
phase_space = {}
amps = {}
locs = ['us', 'ds']
for loc in locs:
    phase_space[loc] = {}
    for var in psv_names:
        phase_space[loc][var] = []
    amps[loc] = []

# Imports data from the requested file
def import_data(loc, file_name):
    # Read the file line by line, append the dictionary
    data_file = open(file_name)
    for line in data_file.readlines():
         psv = {}
         psv['x'], psv['px'], psv['y'], psv['py'], amp = line.split()
         for var in psv_names:
             phase_space[loc][var].append(float(psv[var]))
         amps[loc].append(float(amp))
         
# Initalizes the canvas that will contain all of the phase space sections
def initialize_canvas(min_amp, max_amp):
    # Initialize main canvas
    canv = ROOT.TCanvas("c", "c", 1650, 1200)
    ROOT.gPad.SetRightMargin(.8)
    ROOT.gPad.SetLogz()
    
    # Initialize underlying TH2F (z scale)
    hempty = ROOT.TH2F("empty", "", 2, 0, 1, 2, 0, 1)
    hempty.Fill(0., min_amp, max_amp)
    hempty.Draw("COLZ")
    hempty.GetXaxis().SetLabelOffset(999)
    hempty.GetYaxis().SetLabelOffset(999)
    hempty.GetZaxis().SetTitle("A_{#perp}  [mm]  ")
    hempty.GetZaxis().SetTitleSize(.02)
    hempty.GetZaxis().SetTitleOffset(.9)
    hempty.GetZaxis().SetLabelSize(.02)

    # Set the position of the color palette
    ROOT.gPad.Update()
    palette = hempty.GetListOfFunctions().FindObject("palette")
    palette.SetX1NDC(.94)
    palette.SetX2NDC(.96)
    palette.SetY1NDC(.29)
    palette.SetY2NDC(.98)

    # Initialize the subpads
    pads = []
    for i in range(nvars):
        for j in range(nvars):
            cid = i*nvars+j
            lpadd = 0.
            if j == 0 or cid == 1:
                lpadd = 0.049
            rpadd = 0.
            bpadd = 0.
            if i == 3 or cid == 11:
                bpadd = 0.049
            tpadd = 0.
            pads.append(ROOT.TPad("pad%d" % cid, "",
                .05+.225*j-lpadd, .75-.23*i-bpadd, .275+.225*j+rpadd, .98-.23*i+tpadd))
            pads[cid].SetMargin(lpadd/.27, rpadd/.27, bpadd/.275, tpadd/.275)
            pads[cid].Draw()
  
    return canv, pads, hempty

# Draws a phase space 2D section
def plot_section(loc, x_var, y_var, max_n=100000):
    # Initialize the canvas
    x_name = psv_names[x_var]
    y_name = psv_names[y_var]
    x_label = psv_labels[x_var]+" "
    y_label = psv_labels[y_var]+" "
    title = "amplitude_phase_space_scatter_"+\
                loc+"_"+psv_names[x_var]+"_"+psv_names[y_var]
                
    # Fill the data
    data = []
    for i in range(min(max_n, len(phase_space[loc][x_name]))):
        data.append([phase_space[loc][x_name][i], phase_space[loc][y_name][i], amps[loc][i]])

    # Sort by descending order of amplitude (low amplitude points drawn last)
    data.sort(key=lambda el: el[2], reverse=True)

    # Initalize and draw the graph
    x_data = [el[0] for el in data]
    y_data = [el[1] for el in data]
    amp_list = [el[2] for el in data]
    graph = ROOT.TGraph2D(len(x_data), array('d', x_data), array('d', y_data), array('d', amp_list))
    graph.SetName(title)
    graph.SetTitle(";;;A_{#perp}  [mm]")
    ROOT.gStyle.SetPalette(1)
    graph.Draw("PCOLZ")

    # Set the view right (from above), set the style
    ROOT.gPad.SetTheta(90)
    ROOT.gPad.SetPhi(0)

    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(.5)
    
    # Set axis ranges
    ROOT.gPad.SetLogz()
    ROOT.gStyle.SetOptStat(0)
    hdummy = ROOT.TH2F(title+"_dummy", "", 1, -125+1e-3, 125-1e-3, 1, -125+1e-3, 125-1e-3)
    graph.SetHistogram(hdummy)
    graph.SetMinimum(0.1)
    graph.SetMaximum(100)

    # Remove default x and y axes, replace them by correctly placed TGaxis
    graph.GetXaxis().SetTitleOffset(999)
    graph.GetXaxis().SetLabelOffset(999)
    graph.GetXaxis().SetTickSize(0)

    graph.GetYaxis().SetTitleOffset(999)
    graph.GetYaxis().SetLabelOffset(999)
    graph.GetYaxis().SetTickSize(0)

    wmin = graph.GetXaxis().GetXmin()
    wmax = graph.GetXaxis().GetXmax()
    xaxis = ROOT.TGaxis(-.5775, -.5775, .5775, -.5775, wmin, wmax, 505)
    xaxis.SetTitle(x_label)
    xaxis.SetLabelFont(42)
    xaxis.SetTextFont(42)
    xaxis.Draw("SAME")

    wmin = graph.GetYaxis().GetXmin()
    wmax = graph.GetYaxis().GetXmax()
    yaxis = ROOT.TGaxis(-.5775, -.5775, -.5775, .5775, wmin, wmax, 505)
    yaxis.SetTitle(y_label)
    yaxis.SetLabelFont(42)
    yaxis.SetTextFont(42)
    yaxis.Draw("SAME")
    
    # Remove the z axis
    ROOT.gPad.Update()
    palette = graph.GetHistogram().GetListOfFunctions().FindObject("palette")
    palette.SetX1NDC(100.)
    palette.SetX2NDC(100.)
    palette.SetY1NDC(0.)
    palette.SetY2NDC(1.)
    
    # Set title and label sizes
    cid = int(ROOT.gPad.GetName()[3:])
    size = 0.08
    if cid < 11:
        size = 0.09
    xaxis.SetTitleSize(size)
    xaxis.SetLabelSize(size)
    xaxis.SetTitleOffset(1.)
    yaxis.SetTitleSize(size)
    yaxis.SetLabelSize(size)
    yaxis.SetTitleOffset(.9)
    graph.GetZaxis().SetTickLength(0)
    
    return graph, xaxis, yaxis
    
def draw_info_box(pad):
    pad.cd()
    pad.SetFillColorAlpha(0, 0.)
    info = InfoBox.InfoBox("Internal", "3.2.0", "6-140 LH_{2} full", "2017/03")
    info.SetCoordinates(.175, .25, .875, .75)
    info.Draw()
    return info
    
def draw_diagonal(pad):
    pad.cd()
    pad.SetFillColorAlpha(0, 0.)
    diag = ROOT.TLine(0., 1., 1., 0.)
    diag.Draw()
    return diag
    
def draw_tags(pad):
    pad.cd()
    uptext = ROOT.TLatex(.75, .55, "Upstream")
    uptext.SetTextAngle(90)
    uptext.SetTextFont(62)
    uptext.SetTextSize(.125)
    uptext.SetTextAlign(23)
    uptext.Draw("SAME")

    downtext = ROOT.TLatex(.5, .15, "Downstream")
    downtext.SetTextFont(62)
    downtext.SetTextSize(.125)
    downtext.SetTextAlign(23)
    downtext.Draw("SAME")
    
    return uptext, downtext

if __name__ == "__main__":
    # Import the data
    ROOT.gROOT.SetBatch()
    data_set = "6-140_lH2_full"
    for loc in locs:
        import_data(loc, data_set+"_ps_"+loc)
        
    # Initialize canvas, subpads
    min_amp = 0.1
    max_amp = 100
    canv, pads, palette = initialize_canvas(min_amp, max_amp)
    
    # Draw the upstream distributions in the top right corner
    # Draw the downstream distributions in the bottom left corner
#    shades = [kGray, kGray+1, kWhite, kGray+2, kGray+3, kBlack]
    graphs, xaxes, yaxes = {}, {}, {}
    for loc in locs:
#         shadeid = 0
        for i in range(nvars-1):
            for j in range(i+1, nvars):
                # Draw the graph
                padid = {"us":i*nvars+j, "ds":j*nvars+i}[loc]
                pads[padid].cd()
                ROOT.gPad.SetLogz()
                tag = loc+"_"+psv_names[i]+psv_names[j]
                graphs[tag], xaxes[tag], yaxes[tag] = plot_section(loc, i, j)

	        # If it is the top right corner, invert the axes
                if loc == "us":
                    print "flipping"
                    xtag = xaxes[tag].GetTitle()
                    xaxes[tag].SetTitle(yaxes[tag].GetTitle())
                    yaxes[tag].SetTitle(xtag)
	            for k in range(graphs[tag].GetN()):
	                x = graphs[tag].GetX()[k]
	                y = graphs[tag].GetY()[k]
	                graphs[tag].GetX()[k] = y
	                graphs[tag].GetY()[k] = x
	                
    # Draw tags
    tag_us, tag_ds = draw_tags(pads[0])
	                
    # Draw info box
    info = draw_info_box(pads[15])
    
    # Draw diagonal
    diag1 = draw_diagonal(pads[5])
    diag2 = draw_diagonal(pads[10])

    # Save
    canv.SaveAs("poincare_sections.pdf")
    del canv
