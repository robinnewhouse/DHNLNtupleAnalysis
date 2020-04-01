# Plotting Script
import argparse, os, math, ROOT, glob, uproot, time

import atlas_style
import numpy as np
from ROOT import *
from ROOT import gPad
from pylab import *
gROOT.LoadMacro("AtlasStyle.C")
gROOT.LoadMacro("AtlasUtils.C")
gROOT.LoadMacro("AtlasLabels.C")


#############################################################################################################################################
# globals
histos_savepath = '/home/dtrischuk/HNLAnalysis/SS_bkgStudies/plots/periodB/' # change path here to save your histograms

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
	
	ax = 0.25
	ay = 0.82

	a.DrawLatex(ax,ay,datarun)
	b.DrawLatex(ax,ay-0.05,Vertextype)
	ATLASLabel(0.25,0.87,"Internal")




def plot_cutflow(histos, ch_name, vertextype,savefilename):

	file = ROOT.TFile(histos)
  	hcutflow = file.Get('CutFlow_'+ch_name)

	MyC01= ROOT.TCanvas("MyC01","cutflow",600,400)
	# ROOT.gPad.SetLogy()
  # 	MyC01.Divide(1,1)
 	# MyC01.cd(1)

 	ymax_cutflow = hcutflow.GetMaximum()
 	hcutflow.GetYaxis().SetRangeUser(0,ymax_cutflow*1.05)
 	hcutflow.SetFillColor(kAzure-4)
 	hcutflow.SetLineWidth(0)
  	hcutflow.Draw("HIST TEXT0 SAME")

  	if options.data == True: 
  		drawNotesData("data18 period B",vertextype) 
  	else: 
  		drawNotesMC("",'VSI Leptons',"emu",'20','10') 
  	MyC01.SaveAs(histos_savepath +'Cutflow_'+savefilename+'.pdf')


def compare2(histos, h1name, h1label, h2name, h2label, xlabel,savefilename):
	############################################################################
	nRebin = 5 # set to 1 if you dont want to rebin.
	scaleymax = 1.6 # use this to scale height of y axis for asthetics
	############################################################################


	# get 2 histograms from input file

	file1 = ROOT.TFile(histos)
	file2 = ROOT.TFile(histos)
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


def compare3(histos, nRebin, h1name, h1label, h3name, h2name, h2label, h3label, vertextype, xlabel,savefilename):
	############################################################################
	# nRebin = 5 # set to 1 if you dont want to rebin.
	scaleymax = 1.2 # use this to scale height of y axis for asthetics
	############################################################################


	# get 2 histograms from input file

	file = ROOT.TFile(histos)
  	h1 = file.Get(h1name)
  	h2 = file.Get(h2name)
  	h3 = file.Get(h3name)


  	# print h1.GetEntries()
  	# print h2.GetEntries()

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
	leg01.AddEntry(h3,h3label,"l")


	#rebin histograms
	h1.Rebin(nRebin)
	h2.Rebin(nRebin)
	h3.Rebin(nRebin)

	# find the max of the 2 histograms
	y1_max = h1.GetMaximum()
	y2_max = h2.GetMaximum()
	y3_max = h3.GetMaximum()
	y_max = max(y1_max,y2_max,y3_max) # scale the max for asthetics 

	h1_binxmax = h1.FindLastBinAbove(0,1)
	h2_binxmax = h2.FindLastBinAbove(0,1)
	h3_binxmax = h3.FindLastBinAbove(0,1)
	bin_xmax = max(h1_binxmax,h2_binxmax,h3_binxmax)

	h1_binxmin = h1.FindFirstBinAbove(0,1)
	h2_binxmin = h2.FindFirstBinAbove(0,1)
	h3_binxmin = h3.FindFirstBinAbove(0,1)
	bin_xmin = min(h1_binxmin,h2_binxmin,h3_binxmin)


	h1.SetLineColor(kAzure+6)
	h1.GetXaxis().SetTitle(xlabel)
  	h1.GetYaxis().SetTitle("entries")
  	h1.GetYaxis().SetRangeUser(0,y_max*scaleymax)
	h1.GetXaxis().SetRangeUser(0,40)
	# h1.GetXaxis().SetRange(bin_xmin-1,bin_xmax+1)
	h1.Draw("HIST")

	
 	h2.SetLineColor(kRed)
 	# h2.SetLineStyle(3)
 	h2.Draw("HIST SAME")

 	h3.SetLineColor(kViolet+8)
 	# h2.SetLineStyle(3)
 	h3.Draw("HIST SAME")

 	mDVcut=TLine(4,0,4,y_max)
 	mDVcut.SetLineStyle(3)
	mDVcut.Draw("SAME")

	leg01.Draw()
	ATLASLabel(0.25,0.87,"Internal")

	if options.data == True: 
  		drawNotesData("data18",vertextype) 
	
	MyC01.SaveAs(histos_savepath +savefilename+'.pdf')


def compare4(histos, nRebin, h1name, h1label, h2name, h2label,h3name, h3label,h4name, h4label, vertextype, xlabel,savefilename):
	############################################################################
	# nRebin = 5 # set to 1 if you dont want to rebin.
	scaleymax = 10 # use this to scale height of y axis for asthetics
	############################################################################


	# get 2 histograms from input file

	file = ROOT.TFile(histos)
  	h1 = file.Get(h1name)
  	h2 = file.Get(h2name)
  	h3 = file.Get(h3name)
  	h4 = file.Get(h4name)


  	# print h1.GetEntries()
  	# print h2.GetEntries()

  	#define your canvas
  	
	MyC01= ROOT.TCanvas("MyC01","",600,400)
	ROOT.gPad.SetLogy()
	# MyC01.Divide(1,1)
	# ROOT.gPad.SetLogy()
	# MyC01.cd(1)


	#format legend
	leg01 = ROOT.TLegend(0.60,0.7,0.91,0.92)
  	leg01.SetTextSize(0.035)
  	leg01.SetBorderSize(0)
  	leg01.SetFillColor(kWhite)
  	leg01.SetShadowColor(kWhite)

  	leg01.AddEntry(h1,h1label,"l")
	leg01.AddEntry(h2,h2label,"l")
	leg01.AddEntry(h3,h3label,"l")
	leg01.AddEntry(h4,h4label,"l")


	#rebin histograms
	h1.Rebin(nRebin)
	h2.Rebin(nRebin)
	h3.Rebin(nRebin)
	h4.Rebin(nRebin)

	# find the max of the 2 histograms
	y1_max = h1.GetMaximum()
	y2_max = h2.GetMaximum()
	y3_max = h3.GetMaximum()
	y4_max = h4.GetMaximum()
	y_max = max(y1_max,y2_max,y3_max,y4_max ) # scale the max for asthetics 

	h1_binxmax = h1.FindLastBinAbove(0,1)
	h2_binxmax = h2.FindLastBinAbove(0,1)
	h3_binxmax = h3.FindLastBinAbove(0,1)
	h4_binxmax = h4.FindLastBinAbove(0,1)
	bin_xmax = max(h1_binxmax,h2_binxmax,h3_binxmax,h4_binxmax)

	h1_binxmin = h1.FindFirstBinAbove(0,1)
	h2_binxmin = h2.FindFirstBinAbove(0,1)
	h3_binxmin = h3.FindFirstBinAbove(0,1)
	h4_binxmin = h4.FindFirstBinAbove(0,1)
	bin_xmin = min(h1_binxmin,h2_binxmin,h3_binxmin,h4_binxmin)


	h1.SetLineColor(kAzure+6)
	# h1.SetFillColor(kAzure+6)
	h1.GetXaxis().SetTitle(xlabel)
  	h1.GetYaxis().SetTitle("entries")
  	h1.GetYaxis().SetRangeUser(1,y_max*scaleymax)
	h1.GetXaxis().SetRangeUser(0,50)
	# h1.GetXaxis().SetRange(bin_xmin-1,bin_xmax+1)
	h1.Draw("HIST")

	
 	h2.SetLineColor(kRed)
 	# h2.SetFillColor(kRed)
 	# h2.SetLineStyle(3)
 	h2.Draw("HIST SAME")

 	h3.SetLineColor(kViolet+8)
 	# h3.SetFillColor(kViolet+8)

 	# h2.SetLineStyle(3)
 	h3.Draw("HIST SAME")

 	h4.SetLineColor(kGreen+1)
 	# h4.SetFillColor(kGreen+1)
 	h4.SetLineStyle(3)
 	h4.Draw("HIST SAME")

 	mDVcut=TLine(4,0,4,y_max*10)
 	mDVcut.SetLineStyle(3)
	mDVcut.Draw("SAME")

	leg01.Draw()
	ATLASLabel(0.25,0.87,"Internal")

	if options.data == True: 
  		drawNotesData("data18",vertextype) 
	
	MyC01.SaveAs(histos_savepath +savefilename+'.pdf')



def compare2_wRatio(histos1, histos2, h1name, h1label, h2name, h2label, xlabel,savefilename):
	############################################################################
	nRebin = 10 # set to 1 if you dont want to rebin.
	scaleymax = 1.1 # use this to scale height of y axis for asthetics
	############################################################################


	# get 2 histograms from the file 
	file1 = ROOT.TFile(histos1)
	file2 = ROOT.TFile(histos2)
	h1 = file1.Get(h1name)
	h2 = file2.Get(h2name)


	print h1.GetEntries()
	print h2.GetEntries()

	#define your canvas
	# MyC01 = ROOT.TCanvas("MyC01","",600,400)
	MyC01 = ROOT.TCanvas("MyC01", "canvas", 800, 800)
	MyC01.Divide(1,1)
	MyC01.cd(1)
	pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
	pad1.SetBottomMargin(0)
	pad1.SetGridx()
	pad1.Draw()
	pad1.cd() 

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
	# h1.GetXaxis().SetRangeUser(0,5)
	h1.GetXaxis().SetRange(bin_xmin-1,bin_xmax+1)
	h1.Draw("HIST")

	
	h2.SetLineColor(kViolet+8)
	h2.SetLineStyle(3)
	h2.Draw("HIST SAME")

	# h1.GetYaxis().SetLabelSize(0.)
	# axis = ROOT.TGaxis( -5, 20, -5, 220, 20,220,510,"")
	# axis.SetLabelFont(43) # Absolute font size in pixel (precision 3)
	# axis.SetLabelSize(15)
	# axis.Draw()

	leg01.Draw()
	ATLASLabel(0.25,0.87,"Internal")

	MyC01.cd(1)
	pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
	pad2.SetTopMargin(0)
	pad2.SetBottomMargin(0.2)
	pad2.SetGridx() 
	pad2.Draw()
	pad2.cd() 

	hratio = h1.Clone("hratio")
	hratio.SetLineColor(kBlack)
	hratio.SetMinimum(0.8)  
	hratio.SetMaximum(1.35) 
	hratio.Sumw2()
	hratio.SetStats(0)   
	hratio.Divide(h2)
	hratio.SetMarkerStyle(21)
	hratio.GetYaxis().SetTitle("Orig. VSI / Rerun VSI")
	hratio.GetYaxis().SetNdivisions(505)
   	hratio.GetYaxis().SetTitleSize(20)
   	hratio.GetYaxis().SetTitleFont(43)
   	hratio.GetYaxis().SetTitleOffset(1.55)
   	hratio.GetYaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
   	hratio.GetYaxis().SetLabelSize(15)
	hratio.GetXaxis().SetTitleSize(20)
   	hratio.GetXaxis().SetTitleFont(43)
   	hratio.GetXaxis().SetTitleOffset(4.)
   	hratio.GetXaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
   	hratio.GetXaxis().SetLabelSize(15)


	hratio.Draw("ep")

	
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


	parser.add_argument("-f2", "--file2",
						dest="file2", default=[""],
						action = AppendActionCleanDefault,
						type = str,
						help="Input file 2",
						metavar="FILE")


	parent_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter, parents = [parser]) 

	options = parent_parser.parse_args()


	####################################################################################################################################
	# Here's where you configure what histograms to plot

	# plot_cutflow(options.file[0], ch_name = "emu",
	# 							  vertextype = "Run1",
	# 							  savefilename = "master")


	channel_VSI = ["pel_ee_VSI","pel_emu_VSI","pel_mumu_VSI","pmu_ee_VSI","pmu_emu_VSI","pmu_mumu_VSI"]
	channel_VSILep = ["pel_ee_VSILep","pel_emu_VSILep","pel_mumu_VSILep","pmu_ee_VSILep","pmu_emu_VSILep","pmu_mumu_VSILep"]

	for i in range(len(channel_VSI)):
		plot_cutflow(options.file[0], ch_name =channel_VSI[i],
								  vertextype = "VSI",
								  savefilename = "selLep_periodB_%s"%channel_VSI[i])

		plot_cutflow(options.file[0], ch_name =channel_VSILep[i],
								  vertextype = "VSI Leptons",
								  savefilename = "selLep_periodB_%s"%channel_VSILep[i])


	# plot_cutflow(options.file[0], ch_name ="SS",
	# 							  vertextype = "VSI",
	# 							  savefilename = "VSILep_02_fullSel_SS")

	# plot_cutflow(options.file[0], ch_name ="SS_1mu",
	# 							  vertextype = "VSI",
	# 							  savefilename = "VSILep_02_fullSel_SS_1mu")

	# plot_cutflow(options.file[0], ch_name ="SS_1el",
	# 							  vertextype = "VSI",
	# 							  savefilename = "VSILep_02_fullSel_SS_1el")

	# plot_cutflow(options.file[0], ch_name ="SS_2lep",
	# 							  vertextype = "VSI",
	# 							  savefilename = "VSILep_02_fullSel_SS_2lep")

	# compare4(options.file[0], nRebin=5,
	# 						  h1name="selDV_r_SS",
	# 						  h1label="SS", 
	# 						  h2name="selDV_r_SS_1mu",
	# 						  h2label="SS + 1 muon", 
	# 						  h3name="selDV_r_SS_1el",
	# 						  h3label="SS + 1 electron", 
	# 						  h4name="selDV_r_SS_2lep",
	# 						  h4label="SS + 2 leptons", 
	# 						  vertextype="VSI Leptons",
	# 						  xlabel='r DV [mm]',
	# 						  savefilename='hrDV_selLeptons_VSILep_02')

	compare4(options.file[0], nRebin = 2,
							  h1name="selDV_mass_pel_emu_VSI",
							  h1label="pel emu", 
							  h2name="selDV_mass_pel_mumu_VSI",
							  h2label="pel mumu", 
							  h3name="selDV_mass_pmu_emu_VSI",
							  h3label="pmu emu", 
							  h4name="selDV_mass_pmu_mumu_VSI",
							  h4label="pmu mumu", 
							 
							  vertextype="VSI",
							  xlabel='r DV [mm]',
							  savefilename='hmassDV_up2toLeptons_VSI')



	compare4(options.file[0], nRebin = 2,
							  h1name="selDV_r_pel_emu_VSILep",
							  h1label="pel emu", 
							  h2name="selDV_r_pel_mumu_VSILep",
							  h2label="pel mumu", 
							  h3name="selDV_r_pmu_emu_VSILep",
							  h3label="pmu emu", 
							  h4name="selDV_r_pmu_mumu_VSILep",
							  h4label="pmu mumu", 
							 
							  vertextype="VSI",
							  xlabel='r DV [mm]',
							  savefilename='hrDV_up2toLeptons_VSILeptons')



	# compare2_wRatio(options.file[0],options.file2[0], h1name="DV_r_SS", 
	# 												  h1label="VSI", 	
	# 												  h2name="DV_r_SS",
	# 												  h2label="VSI Leptons", 
	# 												  xlabel='r DV [mm]',
	# 												  savefilename='hrDV_compare2deri')

	####################################################################################################################################








