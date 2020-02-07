# Plotting Script
import argparse, os, math, ROOT, glob, uproot, time


import numpy as np
from ROOT import*
from pylab import *
gROOT.LoadMacro("AtlasStyle.C")
gROOT.LoadMacro("AtlasUtils.C")
gROOT.LoadMacro("AtlasLabels.C")

#############################################################################################################################################
# globals
histos_savepath = '/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/output/VSI_DerivationStudies/' # change path here to save your histograms

#############################################################################################################################################

ROOT.gROOT.SetBatch(True)
SetAtlasStyle()

#get note
def getNote(size=14):
	n = ROOT.TLatex()
	n.SetNDC()
	n.SetTextFont(43)
	n.SetTextColor(1)
	n.SetTextSize(size)
	return n

def drawNotesMC(MC_campaign,Vertextype, DV_type,mass,lifetime):
	a = getNote()
	b = getNote()
	c = getNote()
	d = getNote()
	e = getNote()
	ax = 0.50
	ay = 0.87
	if MC_campaign == "merged": 
		a.DrawLatex(ax,ay,'all MC campaigns')
	else:
		a.DrawLatex(ax,ay,'%s'%MC_campaign) 
	b.DrawLatex(ax,ay-0.05,'mass: %s GeV'%mass)
	c.DrawLatex(ax,ay-0.10,'lifetime: %s mm'%lifetime)
	if DV_type == "emu":
		d.DrawLatex(ax,ay-0.15,'DV type: e\mu')
	elif DV_type == "mumu":
		d.DrawLatex(ax,ay-0.15,'DV type: \mu\mu')
	# if DV_Default == True:
	# 	e.DrawLatex(ax,ay-0.20,'VSI')
	# else:
	e.DrawLatex(ax,ay-0.20,Vertextype)
	ATLASLabel(0.25,0.87,"Internal")

def drawNotesData(datarun,Vertextype):
	a = getNote()
	b = getNote()
	
	ax = 0.5
	ay = 0.87

	a.DrawLatex(ax,ay,datarun)
	b.DrawLatex(ax,ay-0.05,Vertextype)
	ATLASLabel(0.25,0.87,"Internal")




def plot_cutflow(histos, ch_name, vertextype,savefilename):

	file = ROOT.TFile(histos)
  	hcutflow = file.Get('CutFlow_'+ch_name)

	MyC01= ROOT.TCanvas("MyC01","cutflow",600,400)
  	MyC01.Divide(1,1)
 	MyC01.cd(1)

 	ymax_cutflow = hcutflow.GetMaximum()
 	hcutflow.GetYaxis().SetRangeUser(0,ymax_cutflow*1.05)
 	hcutflow.SetFillColor(kAzure-4)
 	hcutflow.SetLineWidth(0)
  	hcutflow.Draw("HIST TEXT0 SAME")

  	if options.data == True: 
  		drawNotesData("data17 test",vertextype) 
  	else: 
  		drawNotesMC("",'VSI Leptons',"emu",'20','10') 
  	MyC01.SaveAs(histos_savepath +'Cutflow_'+savefilename+'.pdf')


def compare2(histos, h1name, h1label, h2name, h2label, xlabel,savefilename):
	############################################################################
	nRebin = 5 # set to 1 if you dont want to rebin.
	scaleymax = 1.6 # use this to scale height of y axis for asthetics
	############################################################################


	# get 2 histograms from the file 
	file = ROOT.TFile(histos)
  	h1 = file.Get(h1name)
  	h2 = file.Get(h2name)

  	print h1.GetEntries()
  	print h2.GetEntries()

  	#define your canvas
	MyC01= ROOT.TCanvas("MyC01","",600,400)
	MyC01.Divide(1,1)
	MyC01.cd(1)

	#format legend
	leg01 = ROOT.TLegend(0.60,0.7,0.91,0.92)
  	leg01.SetTextSize(0.035)
  	leg01.SetBorderSize(0)
  	leg01.SetFillColor(kWhite)
  	leg01.SetShadowColor(kWhite)

  	leg01.AddEntry(h1,h1label,"l")
	leg01.AddEntry(h2,h2label,"l")


	#rebin histograms
	h1.Rebin(nRebin)
	h2.Rebin(nRebin)

	# find the max of the 2 histograms
	y1_max = h1.GetMaximum()
	y2_max = h1.GetMaximum()
	y_max = max(y1_max,y2_max) # scale the max for asthetics 

	h1_binxmax = h1.FindLastBinAbove(0,1)
	h2_binxmax = h2.FindLastBinAbove(0,1)
	bin_xmax = max(h1_binxmax,h2_binxmax)

	h1_binxmin = h1.FindFirstBinAbove(0,1)
	h2_binxmin = h2.FindFirstBinAbove(0,1)
	bin_xmin = min(h1_binxmin,h2_binxmin)


	h1.SetLineColor(kAzure+6)
	h1.GetXaxis().SetTitle(xlabel)
  	h1.GetYaxis().SetTitle("entries")
  	h1.GetYaxis().SetRangeUser(0,y_max*scaleymax)
	h1.GetXaxis().SetRangeUser(0,3000)
	h1.GetXaxis().SetRange(bin_xmin-1,bin_xmax+1)
	h1.Draw("HIST")

	
 	h2.SetLineColor(kViolet+8)
 	h2.SetLineStyle(3)
 	h2.Draw("HIST SAME")

	leg01.Draw()
	ATLASLabel(0.25,0.87,"Internal")
	
	MyC01.SaveAs(histos_savepath +savefilename+'.pdf')




if __name__ == '__main__':
	import argparse
	class AppendActionCleanDefault(argparse._AppendAction):
		def __init__(self, *args, **kwargs):
			super(argparse._AppendAction, self).__init__(*args,**kwargs)
			self.index = 0
		def __call__(self, parser, namespace, values, option_string = None):
			items = argparse._copy.copy(argparse._ensure_value(namespace, self.dest, [])) if self.index else []
			if values:
				self.index += 1
				items.append(values)
				setattr(namespace, self.dest, items)


	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

	parser.add_argument("-d", "--isData",
						action="store_true",
						dest="data",
						help="Is this data?")

	parser.add_argument("-f", "--file",
						dest="file", default=[""],
						action = AppendActionCleanDefault,
						type = str,
						help="Input file",
						metavar="FILE")


	parent_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter, parents = [parser]) 

	options = parent_parser.parse_args()


	####################################################################################################################################
	# Here's where you configure what histograms to plot

	plot_cutflow(options.file[0], ch_name = "data_SS_main",
								  vertextype = "Standard VSI",
								  savefilename = "stdVSI")

	plot_cutflow(options.file[0], ch_name ="data_SS_deri",
								  vertextype = "Standard VSI in Derivation",
								  savefilename = "deriVSI")

	compare2(options.file[0], h1name="DV_r_data_SS_main", 
							  h1label="Standard VSI", 	
							  h2name="DV_r_data_SS_deri",
							  h2label="Standard VSI in Derivation", 
							  xlabel='r DV [mm]',
							  savefilename='hrDV_VSIcompare')

	####################################################################################################################################








