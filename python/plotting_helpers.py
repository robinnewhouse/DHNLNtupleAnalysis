import ROOT
from ROOT import * 
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True
import atlas_style


#get note
def getNote(size=12):
	n = ROOT.TLatex()
	n.SetNDC()
	n.SetTextFont(63)
	n.SetTextColor(1)
	n.SetTextSize(size)
	return n

def drawNote(note, size=14, ax=0.25, ay=0.82):
	a = getNote(size)
	a.DrawLatex(ax,ay,note)
	
def drawNotes(VtxConfig,lumi,truth_plot, channel=None):
	a = getNote(size=20)
	b = getNote(size=20)
	c = getNote(size=20)
	
	ax = 0.22
	ay = 0.82
	# if lumi != "":
	# 	a.DrawLatex(ax,ay,"\sqrt{s}  = 13 TeV, \int Ldt = %s fb^{-1}"%int(lumi))
	# else: 
	# 	a.DrawLatex(ax,ay,"\sqrt{s}  = 13 TeV")
	if not truth_plot:
		if VtxConfig == "VSI_Leptons": 
			b.DrawLatex(ax,ay-0.05,'VSI Leptons')
		else: 
			b.DrawLatex(ax,ay-0.05,'%s'%(VtxConfig))
	if channel != None:
		if  "uuu" in channel:
			c.DrawLatex(ax,ay-0.1,'channel: \mu\mu\mu')
		elif "ueu" in channel:
			c.DrawLatex(ax,ay-0.1,'channel: \mue\mu')
		elif "uee" in channel:
			c.DrawLatex(ax,ay-0.1,'channel: \muee')
		elif "eee" in channel:
			c.DrawLatex(ax,ay-0.1,'channel: eee')
		elif "eeu" in channel:
			c.DrawLatex(ax,ay-0.1,'channel: ee\mu')
		elif "euu" in channel:
			c.DrawLatex(ax,ay-0.1,'channel: e\mu\mu')

	atlas_style.ATLASLabel(0.22,0.87,"Internal")


def drawNotesMC(MC_campaign,Vertextype, channel,mass,lifetime):
	a = getNote()
	b = getNote()
	c = getNote()
	d = getNote()
	e = getNote()
	ax = 0.82
	ay = 0.82
	
	# a.DrawLatex(ax,ay,'%s'%MC_campaign) 
	a.DrawLatex(ax,ay,"\sqrt{s}  = 13 TeV, \int Ldt = %s fb^{-1}"%140) 
	b.DrawLatex(ax,ay-0.05,'(m_{HNL}, c\\tau) = (%s GeV,  %s mm) '%(mass.split("G")[0],lifetime.split("mm")[0]) )
	# b.DrawLatex(ax,ay-0.05,'HNL mass: %s GeV'%mass.split("G")[0])
	# c.DrawLatex(ax,ay-0.10,'HNL lifetime: %s mm'%lifetime.split("mm")[0])
	if  "uuu" in channel:
		d.DrawLatex(ax,ay-0.10,'channel: \mu\mu\mu')
	elif "ueu" in channel:
		d.DrawLatex(ax,ay-0.10,'channel: \mue\mu')
	elif "uee" in channel:
		d.DrawLatex(ax,ay-0.10,'channel: \muee')
	elif "eee" in channel:
		d.DrawLatex(ax,ay-0.10,'channel: eee')
	elif "eeu" in channel:
		d.DrawLatex(ax,ay-0.10,'channel: ee\mu')
	elif "euu" in channel:
		d.DrawLatex(ax,ay-0.10,'channel: e\mu\mu')
	e.DrawLatex(ax,ay-0.15,Vertextype)
	atlas_style.ATLASLabel(0.82,0.87,"Internal")

def drawNotesData(datarun,Vertextype,lumi,channel):
	a = getNote()
	b = getNote()
	c = getNote()
	
	ax = 0.82
	ay = 0.82

	
	a.DrawLatex(ax,ay,datarun)
	b.DrawLatex(ax,ay-0.05,"\sqrt{13} TeV, %s fb^{-1}"%lumi)
	if  "uuu" in channel:
		c.DrawLatex(ax,ay-0.1,'channel: \mu\mu\mu')
	elif "ueu" in channel:
		c.DrawLatex(ax,ay-0.1,'channel: \mue\mu')
	elif "uee" in channel:
		c.DrawLatex(ax,ay-0.1,'channel: \muee')
	elif "eee" in channel:
		c.DrawLatex(ax,ay-0.1,'channel: eee')
	elif "eeu" in channel:
		c.DrawLatex(ax,ay-0.1,'channel: ee\mu')
	elif "euu" in channel:
		c.DrawLatex(ax,ay-0.1,'channel: e\mu\mu')
	c.DrawLatex(ax,ay-.15,Vertextype)
	atlas_style.ATLASLabel(0.82,0.87,"Internal")


def drawNotesVertextype(Vertextype, lumi=1):
	a = getNote()
	b = getNote()
	
	ax = 0.25
	ay = 0.82
	a.DrawLatex(ax,ay,"\sqrt{s}= 13 TeV, %s fb^{-1}"%lumi)
	b.DrawLatex(ax,ay-0.05,Vertextype)
	atlas_style.ATLASLabel(0.25,0.87,"Internal")


x_label_names = {
	"DV_mass": "DV mass [GeV]",
	"trk_pt": "track p_{T} [GeV]",
	"trk_eta": "track \eta",
	"trk_phi": "track \phi",
	"trk_d0": "track d_{0}",
	"mvis": "Visible mass (m_{lll}) [GeV]",
	"dpt": "\Deltap_{T} between tracks in DV [GeV]",
	"deta": "\Delta\eta between tracks in DV",
	"dphi": "\Delta\phi between tracks in DV",
	"dR": "\DeltaR between tracks in DV",
	"mtrans": "m_{T} [GeV]",
	"HNLm": "HNL mass [GeV]",
	"HNLpt": "HNL p_{T} [GeV]",
	"HNLphi": "HNL \phi",
	"HNLeta": "HNL \eta",
	"DV_r": "DV r [mm]",
	"redmass": "reduced DV mass [GeV]",
	"redmassvis": "reduced visible mass [GeV]",
	"redmassHNL": "reduced HNL mass [GeV]",
	"DV_trk_d0_wrtSV": "d0 wrt SV",
	"DV_trk_z0_wrtSV": "z0 wrt SV"
}

 
def get_x_label(hist_name):
	if hist_name in x_label_names:
		return x_label_names[hist_name]
	else:
		return hist_name


def xlabelhistograms(hist):
	if "DV_r" in hist:
		if  ("redmassvis" in hist):
			return "reduced visible mass [GeV]"
		elif  ("redmass" in hist):
			if "redmassHNL" in hist:
				return "reduced HNL mass [GeV]"
			else:
				return "reduced DV mass [GeV]"
		else:
			return "DV r [mm]"
	if "DV_mass" in hist:
		return "DV mass [GeV]"
	if "trk_pt" in hist:
		return "track p_{T} [GeV]"
	if "trk_eta" in hist:
		return "track \eta"
	if "trk_phi" in hist:
		return "track \phi"
	if "trk_d0" in hist:
		return "track d_{0}"
	if "mvis" in hist:
		return "Visible mass (m_{lll}) [GeV]"
	if "dpt" in hist:
		return "\Deltap_{T} between tracks in DV [GeV]"
	if "deta" in hist:
		return "\Delta\eta between tracks in DV"
	if "dphi" in hist:
		return "\Delta\phi between tracks in DV"
	if "dR" in hist:
		return "\DeltaR between tracks in DV"
	if "mtrans" in hist:
		return "m_{T} [GeV]"
	if "HNLm" in hist:
		return "HNL mass [GeV]"
	if "HNLpt" in hist:
		return "HNL p_{T} [GeV]"
	if "HNLphi" in hist:
		return "HNL \phi"
	if "HNLeta" in hist:
		return "HNL \eta"
	else:
		return ""

def histColours(nhist): 
	color = ROOT.TColor()
	if nhist== 5:
		ncolor = color.GetColor("#FF1B28") # red 
	if nhist== 4:
		ncolor = color.GetColor("#FF7A43") # orange
	if nhist== 3:	
		ncolor = color.GetColor("#FFE425") # yellow
	if nhist== 2:	
		ncolor = color.GetColor("#39A230")# green
	if nhist== 0:
		ncolor = color.GetColor("#3ACEFF") # light blue
		# ncolor = color.GetColor("#2479B5") # dark blue
	if nhist== 1:
		ncolor = color.GetColor("#B576FF") # purple
	
	
	return ncolor

def bkgColours(nhist): 
	color = ROOT.TColor()
	if nhist== 0:
		ncolor = color.GetColor("#3F4FC0") #deep blue 
		# ncolor = color.GetColor("#c2a5cf") #light purple
	if nhist== 1:
		ncolor = color.GetColor("#a6dba0")# light green
	if nhist== 2:
		ncolor = color.GetColor("#7b3294") # dark purple
	if nhist== 3:
		ncolor = color.GetColor("#008837") # dark green
	return ncolor
	



