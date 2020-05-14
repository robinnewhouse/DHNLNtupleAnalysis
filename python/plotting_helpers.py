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

	
def drawNotes(channel,VtxConfig,lumi):
	a = getNote()
	b = getNote()
	c = getNote()
	d = getNote()
	e = getNote()
	f = getNote()
	ax = 0.22
	ay = 0.82

	b.DrawLatex(ax,ay,"\sqrt{s} = 13 TeV, \int Ldt = %s fb^{-1}"%lumi)
	if  "uuu" in channel:
		a.DrawLatex(ax,ay-0.05,'channel: \mu\mu\mu')
	elif "ueu" in channel:
		a.DrawLatex(ax,ay-0.05,'channel: \mue\mu')
	elif "uee" in channel:
		a.DrawLatex(ax,ay-0.05,'channel: \muee')
	elif "eee" in channel:
		a.DrawLatex(ax,ay-0.05,'channel: eee')
	elif "eeu" in channel:
		a.DrawLatex(ax,ay-0.05,'channel: ee\mu')
	elif "euu" in channel:
		a.DrawLatex(ax,ay-0.05,'channel: e\mu\mu')

	# c.DrawLatex(ax,ay-0.15,'(m_{HNL}, c\\tau) = (10, 100)')
	
	c.DrawLatex(ax,ay-0.1,'%s'%(VtxConfig))

	# else: 
	# 	a.DrawLatex(ax,ay,'%s'%MC_campaign) 
	# 	b.DrawLatex(ax,ay-0.05,'mass: %s GeV'%mass)
	# 	c.DrawLatex(ax,ay-0.10,'lifetime: %s mm'%lifetime)
	# 	if plepton == "muon":
	# 		d.DrawLatex(ax,ay-0.15,'Prompt muon')
	# 	if plepton == "electron":
	# 		d.DrawLatex(ax,ay-0.15,'Prompt electron')
	# 	if DV_type == "mumu":
	# 		e.DrawLatex(ax,ay-0.20,'DV type: \mu\mu\\nu')
	# 	if DV_type == "emu":
	# 		e.DrawLatex(ax,ay-0.20,'DV type: e\mu\\nu')
	# 	if DV_type == "ee":
	# 		e.DrawLatex(ax,ay-0.20,'DV type: ee\\nu')
	# 	f.DrawLatex(ax,ay-0.25,'%s'%(VtxConfig))
	# 	# else:
		# 	e.DrawLatex(ax,ay-0.20,'VSI Leptons')
	atlas_style.ATLASLabel(0.22,0.87,"Internal")


def drawNotesMC(MC_campaign,Vertextype, channel,mass,lifetime):
	a = getNote()
	b = getNote()
	c = getNote()
	d = getNote()
	e = getNote()
	ax = 0.82
	ay = 0.82
	
	a.DrawLatex(ax,ay,'%s'%MC_campaign) 
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
	if nhist== 0:
		return kAzure+6
	if nhist== 1:
		return kViolet+8
	if nhist== 2:
		return kRed
	if nhist== 3:
		return kGreen+1
	if nhist== 4:
		return kOrange -3
	else: 
		return kBlack

