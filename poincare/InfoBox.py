import ROOT

class InfoBox:

    text = ROOT.TPaveText()	# Info box
    typ = ""			# Type of data
    maus = ""			# Version of MAUS
    run = ""			# Run tag
    cycle = ""			# Type of plot    

    def __init__(self, typ, maus, run="", cycle=""):
	# Record the info
        self.typ = typ
	self.maus = maus
	self.run = run
	self.cycle = cycle

	# Initialize the text box in the top right corner by default
	self.text = ROOT.TPaveText(.7, .7, .875, .875, "NDC")
  	self.text.SetLineColorAlpha(0, 0)
  	self.text.SetFillStyle(0)
  	self.text.SetTextAlign(12)
	self.text.SetTextFont(42)

    def SetPosition(self, pos, alpha=.02, beta=.16):
        # Get the limits of the text box from the pos variable
	if "r" in pos:
	    xmin = .9-alpha-beta
	    xmax = .9-alpha
	else: 
	    xmin = .1+alpha
	    xmax = .1+alpha+beta

	if "t" in pos:
	    ymin = .9-alpha-beta
	    ymax = .9-alpha
	else: 
	    ymin = .1+alpha
	    ymax = .1+alpha+beta

	# Set the new coordinates
	self.SetCoordinates(xmin, ymin, xmax, ymax)

    def SetCoordinates(self, xmin, ymin, xmax, ymax):
  	if self.text:
            del self.text
	self.text = ROOT.TPaveText(xmin, ymin, xmax, ymax, "NDC")
  	self.text.SetLineColorAlpha(0, 0);
  	self.text.SetFillStyle(0);
  	self.text.SetTextAlign(12);
  	self.text.SetTextFont(42);

    def Draw(self):
	# Clear the TPaveText
  	self.text.Clear()

	# Add the relevent info
        t = self.text.AddText("#scale[1.5]{MICE} %s" % self.typ)
        t.SetTextFont(62)
        if len(self.cycle):
            self.text.AddText("ISIS Cycle %s" % self.cycle)
  	if len(self.run):
            self.text.AddText("Run %s" % self.run)
   	if len(self.maus):
      	    self.text.AddText("MAUS v%s" % self.maus)

  	self.text.Draw("SAME");
