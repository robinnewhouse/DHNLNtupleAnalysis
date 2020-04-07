# Plotting Script
import argparse, os, math, ROOT, glob, uproot, time, json

import atlas_style
import numpy as np
import helpers
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
		helpers.drawNotesData("data18 period B",vertextype) 
	else: 
		helpers.drawNotesMC("",'VSI Leptons',"emu",'20','10') 
	MyC01.SaveAs(histos_savepath +'Cutflow_'+savefilename+'.pdf')

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





def compareN(file, hname, hlabel,savefilename,vertextype,setxrange="",scaleymax=1,nRebin=1,setlogy=False):
	############################################################################
	# nRebin = 5 # set to 1 if you dont want to rebin.
	# scaleymax = 1.6 # use this to scale height of y axis for asthetics
	############################################################################


	# get 2 histograms from input file
	nhist = len(hname)
	f = ROOT.TFile(file)
	# file2 = ROOT.TFile(histos)
	print file
	
	h = {}
	for i in xrange(nhist): 
		# print hname[i]
		h[hname[i]] = f.Get(hname[i])
	# exit()


	#define your canvas
	MyC01= ROOT.TCanvas("MyC01","",600,400)
	# MyC01.Divide(1,1)
	# MyC01.cd(1)

	#format legend
	leg01 = ROOT.TLegend(0.60,0.7,0.91,0.92)
	leg01.SetTextSize(0.035)
	leg01.SetBorderSize(0)
	leg01.SetFillColor(kWhite)
	leg01.SetShadowColor(kWhite)


	ymax_list = []
	h_binxmax_list = []
	h_binxmin_list = []
	for i in xrange(nhist):
		leg01.AddEntry(h[hname[i]],hlabel[i],"lp")
		h[hname[i]].Rebin(nRebin)

		ymax_list.append(h[hname[i]].GetMaximum())
		h_binxmax_list.append(h[hname[i]].FindLastBinAbove(0,1))
		h_binxmin_list.append(h[hname[i]].FindFirstBinAbove(0,1))

	y_max = max(ymax_list)
	bin_xmax = max(h_binxmax_list)
	bin_xmin = min(h_binxmin_list)

	shapelist = [20,22,21,33,29]
	for i in xrange(nhist): 
		h[hname[i]].SetLineColor(histColours(i))
		h[hname[i]].SetMarkerColor(histColours(i))
		h[hname[i]].SetMarkerSize(0.7)
		h[hname[i]].SetMarkerStyle(shapelist[i])
		if i == 0:
			h[hname[i]].GetXaxis().SetTitle(xlabelhistograms(hname[i]))
			h[hname[i]].GetYaxis().SetTitle("entries")
			h[hname[i]].GetYaxis().SetRangeUser(0,y_max*scaleymax)
			if setxrange != "":
				X_minmax = [item for item in setxrange.split(' ')]
				h[hname[i]].GetXaxis().SetRangeUser(int(X_minmax[0]),int(X_minmax[1]))
			else:
				h[hname[i]].GetXaxis().SetRange(bin_xmin-1,bin_xmax+1)
			h[hname[i]].Draw("E0 HIST")
		else: 
			h[hname[i]].Draw("E0 HIST SAME")
	
	if setlogy: 
		gPad.SetLogy()



	leg01.Draw()
	ATLASLabel(0.25,0.87,"Internal")

	if "pmu" in hname[0]: 
		if "mumu" in hname[0]:
			helpers.drawNotesMC("",vertextype, "mumumu","10","10")
		# if "mumu" in hname[0]:
		# 	helpers.drawNotesMC("",vertextype, "mumu","10","10")
		if "emu" in hname[0]:
			helpers.drawNotesMC("",vertextype, "emumu","10","10")

	if "DV_mass" in hname[0]:
		mDVcut=TLine(4,0,4,y_max)
		mDVcut.SetLineStyle(3)
		mDVcut.Draw("SAME")

	if "mvis" in hname[0]:
		min_mlll=TLine(50,0,50,y_max)
		min_mlll.SetLineStyle(3)
		min_mlll.Draw("SAME")

		max_mlll=TLine(84,0,84,y_max)
		max_mlll.SetLineStyle(3)
		max_mlll.Draw("SAME")

	if "DV_r" in hname[0]:
		if  "redmass" in hname[0]:
			pass
		else: 
			matlayers = [33.25,50.5,88.5,122.5,299]
			nmatlayers = len(matlayers)
			matlay={}
			for i in range(nmatlayers):
				print "plotting material layers!"
				matlay[i]=TLine(matlayers[i],0,matlayers[i],y_max)
				matlay[i].SetLineStyle(3)
				matlay[i].Draw("SAME")
			leg01.AddEntry(matlay[0],"material layers","l")

	

	
	
	MyC01.SaveAs(histos_savepath +savefilename+'.pdf')



def compare_dataMC(datafile,mcfile, nRebin, hdataname, hdatalabel, hmcname, hmclabel, setrange,scaleymax, vertextype,savefilename,setlogy=False):
	############################################################################
	# nRebin = 5 # set to 1 if you dont want to rebin.
	# scaleymax = 1.2 # use this to scale height of y axis for asthetics
	############################################################################


	# get 2 histograms from input file

	data = ROOT.TFile(datafile)
	mc = ROOT.TFile(mcfile)
	hdata = data.Get(hdataname)
	hmc = mc.Get(hmcname)


	#define your canvas
	MyC01= ROOT.TCanvas("MyC01","",600,400)

	#format legend
	leg01 = ROOT.TLegend(0.60,0.7,0.91,0.92)
	leg01.SetTextSize(0.035)
	leg01.SetBorderSize(0)
	leg01.SetFillColor(kWhite)
	leg01.SetShadowColor(kWhite)

	leg01.AddEntry(hdata,hdatalabel,"lp")
	leg01.AddEntry(hmc,hmclabel,"lp")



	#rebin histograms
	hdata.Rebin(nRebin)
	hmc.Rebin(nRebin)

	norm = 1
	scale_data = norm/(hdata.Integral())

	scale_mc = norm/(hmc.Integral())

	hdata.Scale(scale_data)

	hmc.Scale(scale_mc)	

	# find the max of the 2 histograms
	y1_max = hdata.GetMaximum()
	y2_max = hmc.GetMaximum()
	y_max = max(y1_max,y2_max) # scale the max for asthetics 

	hdata_binxmax = hdata.FindLastBinAbove(0,1)
	hmc_binxmax = hmc.FindLastBinAbove(0,1)
	bin_xmax = max(hdata_binxmax,hmc_binxmax)
	hdata_binxmin = hdata.FindFirstBinAbove(0,1)
	hmc_binxmin = hmc.FindFirstBinAbove(0,1)
	bin_xmin = min(hdata_binxmin,hmc_binxmin)


	hdata.SetLineColor(kBlack)
	hdata.GetXaxis().SetTitle(xlabelhistograms(hdataname))
	hdata.GetYaxis().SetTitle("entries")
	hdata.GetYaxis().SetRangeUser(0,y_max*scaleymax)


	if setrange == "":
		hdata.GetXaxis().SetRange(bin_xmin-1,bin_xmax+1)
	else: 

		minmax = [item for item in setrange.split(' ')]
		hdata.GetXaxis().SetRangeUser(int(minmax[0]),int(minmax[1]))
	
	hdata.SetMarkerSize(0.7)
	hmc.SetMarkerSize(0.7)
	hdata.Draw("E0 HIST")

	
	hmc.SetLineColor(kAzure+6)
	hmc.SetMarkerColor(kAzure+6)
	hmc.SetMarkerStyle(kFullTriangleUp)
	# hmc.SetLineStyle(3)
	hmc.Draw("E0 HIST SAME")

	if setlogy: 
		gPad.SetLogy()

	if "DV_mass" in hdataname:
		mDVcut=TLine(4,0,4,y_max)
		mDVcut.SetLineStyle(3)
		mDVcut.Draw("SAME")

	if "mvis" in hdataname:
		min_mlll=TLine(50,0,50,y_max)
		min_mlll.SetLineStyle(3)
		min_mlll.Draw("SAME")

		max_mlll=TLine(84,0,84,y_max)
		max_mlll.SetLineStyle(3)
		max_mlll.Draw("SAME")

	if "DV_r" in hdataname:
		if  "redmass" in hdataname:
			pass
		else: 
			matlayers = [33.25,50.5,88.5,122.5,299]
			nmatlayers = len(matlayers)
			matlay={}
			for i in range(nmatlayers):
				print "plotting material layers!"
				matlay[i]=TLine(matlayers[i],0,matlayers[i],y_max)
				matlay[i].SetLineStyle(3)
				matlay[i].Draw("SAME")
			leg01.AddEntry(matlay[0],"material layers","l")

	


	leg01.Draw()
	ATLASLabel(0.25,0.87,"Internal")
	helpers.drawNotesData("data 2018 period B",vertextype)
	# if "pmu" in hdataname: 
	# 	if "mumu" in hdataname:

	# 		helpers.drawNotes("mumu","muon",vertextype) 
	# 	if "emu" in hdataname:
	# 		helpers.drawNotes("emu","muon",vertextype) 
	
	MyC01.SaveAs(histos_savepath +savefilename+'.pdf')


def compare_dataMC_multi(datafile,mcfiles, nRebin, hdataname, hdatalabel, hmcnames, hmclabels, setrange,scaleymax, vertextype,savefilename,setlogy=False):
	############################################################################
	# nRebin = 5 # set to 1 if you dont want to rebin.
	# scaleymax = 1.2 # use this to scale height of y axis for asthetics
	############################################################################


	# get 2 histograms from input file

	data = ROOT.TFile(datafile)
	hdata = data.Get(hdataname)
	
	nmc_files = len(mcfiles)
	hmc ={}
	mc = {}


	for i in range(nmc_files): 
		print mcfiles[i]
		print hmcnames[i]
		mc[hmcnames[i]] = ROOT.TFile(mcfiles[i])
		hmc[hmcnames[i]] = mc[hmcnames[i]].Get(hmcnames[i])
		print hmc[hmcnames[i]]

	# exit()
	# h = {}
	# for i in xrange(nhist): 
	# 	# print hname[i]
	# 	h[hname[i]] = f.Get(hname[i])


	#define your canvas
	MyC01= ROOT.TCanvas("MyC01","",600,400)

	#format legend
	leg01 = ROOT.TLegend(0.50,0.7,0.91,0.92)
	leg01.SetTextSize(0.035)
	leg01.SetBorderSize(0)
	leg01.SetFillColor(kWhite)
	leg01.SetShadowColor(kWhite)

	leg01.AddEntry(hdata,hdatalabel,"lp")

	#rebin histograms
	hdata.Rebin(nRebin)

	for i in range(nmc_files): 
		leg01.AddEntry(hmc[hmcnames[i]],hmclabels[i],"lp")
		hmc[hmcnames[i]].Rebin(nRebin)

	norm = 1
	scale_data = norm/(hdata.Integral())
	hdata.Scale(scale_data)

	for i in range(nmc_files): 
		scale_mc = norm/(hmc[hmcnames[i]].Integral())
		hmc[hmcnames[i]].Scale(scale_mc)	

	# find the max of the n histograms
	ydata_max = hdata.GetMaximum()
	hdata_binxmax = hdata.FindLastBinAbove(0,1)
	hdata_binxmin = hdata.FindFirstBinAbove(0,1)

	ymax_list = []
	h_binxmax_list = []
	h_binxmin_list = []

	ymax_list.append(ydata_max)
	h_binxmax_list.append(hdata_binxmax)
	h_binxmin_list.append(hdata_binxmin)

	for i in xrange(nmc_files):
		ymax_list.append(hmc[hmcnames[i]].GetMaximum())
		h_binxmax_list.append(hmc[hmcnames[i]].FindLastBinAbove(0,1))
		h_binxmin_list.append(hmc[hmcnames[i]].FindFirstBinAbove(0,1))

	y_max = max(ymax_list)
	bin_xmax = max(h_binxmax_list)
	bin_xmin = min(h_binxmin_list)

	hdata.SetLineColor(kBlack)
	hdata.GetXaxis().SetTitle(xlabelhistograms(hdataname))
	hdata.GetYaxis().SetTitle("entries")
	hdata.GetYaxis().SetRangeUser(0,y_max*scaleymax)
	hdata.SetMarkerSize(0.7)
	
	if setrange == "":
		hdata.GetXaxis().SetRange(bin_xmin-1,bin_xmax+1)
	else: 

		minmax = [item for item in setrange.split(' ')]
		hdata.GetXaxis().SetRangeUser(int(minmax[0]),int(minmax[1]))

	hdata.Draw("E0 HIST")

	shapelist = [22,21,33,29]
	for i in xrange(nmc_files): 
		hmc[hmcnames[i]].SetMarkerSize(0.7)
		hmc[hmcnames[i]].SetLineColor(histColours(i))
		hmc[hmcnames[i]].SetMarkerColor(histColours(i))
		hmc[hmcnames[i]].SetMarkerStyle(shapelist[i])
		hmc[hmcnames[i]].Draw("E0 HIST SAME")


	if setlogy: 
		gPad.SetLogy()

	if "DV_mass" in hdataname:
		mDVcut=TLine(4,0,4,y_max)
		mDVcut.SetLineStyle(3)
		mDVcut.Draw("SAME")

	if "mvis" in hdataname:
		min_mlll=TLine(50,0,50,y_max)
		min_mlll.SetLineStyle(3)
		min_mlll.Draw("SAME")

		max_mlll=TLine(84,0,84,y_max)
		max_mlll.SetLineStyle(3)
		max_mlll.Draw("SAME")

	if "DV_r" in hdataname:
		if  "redmass" in hdataname:
			pass
		else: 
			matlayers = [33.25,50.5,88.5,122.5,299]
			nmatlayers = len(matlayers)
			matlay={}
			for i in range(nmatlayers):
				print "plotting material layers!"
				matlay[i]=TLine(matlayers[i],0,matlayers[i],y_max)
				matlay[i].SetLineStyle(3)
				matlay[i].Draw("SAME")
			leg01.AddEntry(matlay[0],"material layers","l")

	


	leg01.Draw()
	ATLASLabel(0.25,0.87,"Internal")
	helpers.drawNotesData("",vertextype)
	# if "pmu" in hdataname: 
	# 	if "mumu" in hdataname:

	# 		helpers.drawNotes("mumu","muon",vertextype) 
	# 	if "emu" in hdataname:
	# 		helpers.drawNotes("emu","muon",vertextype) 
	
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


def massplot2D(file, hname, hlabel,savefilename,vertextype, setxrange="",setyrange="",rebinx=1,rebiny=1):
	############################################################################
	# nRebin = 5 # set to 1 if you dont want to rebin.
	# scaleymax = 1.6 # use this to scale height of y axis for asthetics
	############################################################################


	# get 2 histograms from input file
	nhist = len(hname)
	f = ROOT.TFile(file)
	# file2 = ROOT.TFile(histos)
	
	h = f.Get(hname)
	# exit()


	#define your canvas
	MyC01= ROOT.TCanvas("MyC01","",800,500)
	MyC01.SetRightMargin(20)
	# MyC01.Divide(1,1)
	# MyC01.cd(1)


	#format legend
	leg01 = ROOT.TLegend(0.550,0.8,0.89,0.92)
	leg01.SetTextSize(0.035)
	leg01.SetBorderSize(0)
	leg01.SetFillColor(kWhite)
	leg01.SetShadowColor(kWhite)
	leg01.AddEntry(h,hlabel,"")
	# rebin
	hold = h.Clone()
	hold.SetDirectory(0)
	nbinsx = hold.GetXaxis().GetNbins()
	nbinsy = hold.GetYaxis().GetNbins()
	xmin  = hold.GetXaxis().GetXmin()
	xmax  = hold.GetXaxis().GetXmax()
	ymin  = hold.GetYaxis().GetXmin()
	ymax  = hold.GetYaxis().GetXmax()
	nx = nbinsx/rebinx
	ny = nbinsy/rebiny
	h.SetBins(nx,xmin,xmax,ny,ymin,ymax)
	
	#loop on all bins to reset contents and errors
	for biny in range(1, nbinsy): 
		if (biny <= nbinsy): 
			for binx in range(1, nbinsx):
				if binx <= nbinsx: 
					ibin = h.GetBin(binx,biny)
					h.SetBinContent(ibin,0)
			  
   
   #//loop on all bins and refill
	for biny in range(1, nbinsy): 
		if (biny <= nbinsy): 
			by  = hold.GetYaxis().GetBinCenter(biny)
			iy  = h.GetYaxis().FindBin(by)
			for binx in range(1, nbinsx):
				if binx <= nbinsx: 
					bx = hold.GetXaxis().GetBinCenter(binx)
					ix  = h.GetXaxis().FindBin(bx)
					Bin = hold.GetBin(binx,biny)
					ibin= h.GetBin(ix,iy)
					cu  = hold.GetBinContent(Bin)
					h.AddBinContent(ibin,cu)
  

	if setxrange != "": 
		X_minmax = [item for item in setxrange.split(' ')]
		h.GetXaxis().SetRangeUser(int(X_minmax[0]),int(X_minmax[1]))

	if setyrange != "": 
		Y_minmax = [item for item in setyrange.split(' ')]
		h.GetYaxis().SetRangeUser(int(Y_minmax[0]),int(Y_minmax[1]))

	if "DVmass_mvis" in hname: 
		h.GetXaxis().SetTitle("DV mass [GeV]")
		h.GetYaxis().SetTitle("Visible Mass [GeV]")
	if "DVmass_mhnl" in hname: 
		h.GetXaxis().SetTitle("DV mass [GeV]")
		h.GetYaxis().SetTitle("HNL Mass [GeV]")
	if "DVmass_mtrans" in hname: 
		h.GetXaxis().SetTitle("DV mass [GeV]")
		h.GetYaxis().SetTitle("Visible m_{T} [GeV]")
	if "mvis_mhnl" in hname: 
		h.GetXaxis().SetTitle("Visible mass [GeV]")
		h.GetYaxis().SetTitle("HNL mass [GeV]")
	if "mvis_mtrans" in hname: 
		h.GetXaxis().SetTitle("Visible Mass [GeV]")
		h.GetYaxis().SetTitle("Visible m_{T} [GeV]")
	if "mhnl_mtrans" in hname: 
		h.GetXaxis().SetTitle("HNL mass [GeV]")
		h.GetYaxis().SetTitle("Visible m_{T} [GeV]")

	h.GetZaxis().SetTitleOffset(-1.5)
	h.Draw("colz")


	leg01.Draw()
	ATLASLabel(0.25,0.87,"Internal")
	if "pmu" in hname: 
		if "mumu" in hname:
			helpers.drawNotes("mumu","muon",vertextype) 
		if "emu" in hname:
			helpers.drawNotes("emu","muon",vertextype) 
		
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
						help="data inputfile ",
						metavar="FILE")


	parser.add_argument("-f2", "--file2",
						dest="file2", default=[""],
						action = AppendActionCleanDefault,
						type = str,
						help="Input file 2",
						metavar="FILE")

	# parser.add_argument("--config",
	# 					dest="config",
	# 					type = str,
	# 					required = True,
	# 					help="Input config file for plotHisotgrams.py.")


	parent_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter, parents = [parser]) 

	options = parent_parser.parse_args()



	# with open(options.config, 'r') as json_config:
	# 	config_file = json.load(json_config) # load JSON config file that contains a channel name mapped to a list of selections

	# print options.file[0]

	files = [item for item in options.file[0].split(',')]
	# print files 

	# compare_dataMC(files[3],files[1],
	# 			nRebin = 1, 
	# 			setrange= "0 20",
	# 			scaleymax = 1.2,				
	# 			hdataname = "charge_DV_mass_pmu_emu_VSILep", 
	# 			hdatalabel= "data 2018 period B",
	# 			hmcname = "charge_DV_mass_pmu_emu_VSILep",
	# 			hmclabel = "(m_{HNL}, c\\tau, DVtype) = (10, 10, e\mu)",
	# 			vertextype = "VSI Leptons",
	# 			savefilename = "TEST_EMU"	)

	# compare_dataMC_multi(datafile=files[3],
	# 					mcfiles=[files[0],files[1]], 
	# 					nRebin=1, 
	# 					hdataname="charge_DV_mass_pmu_mumu_VSILep", 
	# 					hdatalabel= "data 2018 period B",
	# 					hmcnames=["charge_DV_mass_pmu_mumu_VSILep","charge_DV_mass_pmu_emu_VSILep"], 
	# 					hmclabels=["(m_{HNL}, c\\tau, DVtype) = (10, 10, \mu\mu)", "(m_{HNL}, c\\tau, DVtype) = (10, 10, e\mu)"],
	# 					setrange= "0 20",
	# 					scaleymax= 1.2, 
	# 					vertextype="VSI Leptons",
	# 					savefilename= "Test_multidatacompare" )


	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "0 20",
					scaleymax = 1.2,				
					hdataname = "charge_DV_mass_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_mass_pmu_mumu_VSILep","charge_DV_mass_pmu_emu_VSILep","charge_DV_mass_pel_ee_VSILep"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)" ],
					vertextype = "VSI Leptons",
					savefilename = "DV_mass_MCdatacompare_VSILep"	)


	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "0 40",
					scaleymax = 1.2,				
					hdataname = "charge_DV_mass_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_mass_pmu_mumu_VSI","charge_DV_mass_pmu_emu_VSI","charge_DV_mass_pel_ee_VSI"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI",
					savefilename = "DV_mass_MCdatacompare_VSI"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 5, 
					setrange= "",
					scaleymax = 1.2,
					hdataname = "charge_DV_r_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_r_pmu_mumu_VSILep","charge_DV_r_pmu_emu_VSILep","charge_DV_r_pel_ee_VSILep"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI Leptons",
					savefilename = "DV_r_MCdatacompare_VSILep"	)


	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 5, 
					setrange= "",
					scaleymax = 1.2,
					hdataname = "charge_DV_r_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_r_pmu_mumu_VSI","charge_DV_r_pmu_emu_VSI","charge_DV_r_pel_ee_VSI"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI",
					savefilename = "DV_r_MCdatacompare_VSI"	)


	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 2, 
					setrange= "0 100",
					scaleymax = 1.2,					
					hdataname = "charge_DV_trk_pt_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_trk_pt_pmu_mumu_VSILep","charge_DV_trk_pt_pmu_emu_VSILep","charge_DV_trk_pt_pel_ee_VSILep"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI Leptons",
					savefilename = "trk_pt_MCdatacompare_VSILep"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 2, 
					setrange= "0 100",
					setlogy= False,
					scaleymax = 1.2,					
					hdataname = "charge_DV_trk_pt_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_trk_pt_pmu_mumu_VSI","charge_DV_trk_pt_pmu_emu_VSI","charge_DV_trk_pt_pel_ee_VSI"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI",
					savefilename = "trk_pt_MCdatacompare_VSI"	)



	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "",
					scaleymax = 2.2,
					hdataname = "charge_DV_trk_phi_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_trk_phi_pmu_mumu_VSILep","charge_DV_trk_phi_pmu_emu_VSILep","charge_DV_trk_phi_pel_ee_VSILep"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI Leptons",
					savefilename = "trk_phi_MCdatacompare_VSILep"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "",
					scaleymax = 2.2,
					hdataname = "charge_DV_trk_phi_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_trk_phi_pmu_mumu_VSI","charge_DV_trk_phi_pmu_emu_VSI","charge_DV_trk_phi_pel_ee_VSI"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI",
					savefilename = "trk_phi_MCdatacompare_VSI"	)


	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "",
					scaleymax = 2.2,
					hdataname = "charge_DV_trk_eta_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_trk_eta_pmu_mumu_VSILep","charge_DV_trk_eta_pmu_emu_VSILep","charge_DV_trk_eta_pel_ee_VSILep"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI Leptons",
					savefilename = "trk_eta_MCdatacompare_VSILep"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "",
					scaleymax = 2.2,
					hdataname = "charge_DV_trk_eta_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_trk_eta_pmu_mumu_VSI","charge_DV_trk_eta_pmu_emu_VSI","charge_DV_trk_eta_pel_ee_VSI"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI",
					savefilename = "trk_eta_MCdatacompare_VSI"	)


	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "-10 10",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_d0_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_trk_d0_pmu_mumu_VSILep","charge_DV_trk_d0_pmu_emu_VSILep","charge_DV_trk_d0_pel_ee_VSILep"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI Leptons",
					savefilename = "trk_d0_MCdatacompare_VSILep"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "-10 10",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_d0_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_trk_d0_pmu_mumu_VSI","charge_DV_trk_d0_pmu_emu_VSI","charge_DV_trk_d0_pel_ee_VSI"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI",
					savefilename = "trk_d0_MCdatacompare_VSI"	)


	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "-10 10",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_d0_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_trk_d0_pmu_mumu_VSILep","charge_DV_trk_d0_pmu_emu_VSILep","charge_DV_trk_d0_pel_ee_VSILep"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI Leptons",
					savefilename = "trk_d0_MCdatacompare_VSILep"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "-10 10",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_d0_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_trk_d0_pmu_mumu_VSI","charge_DV_trk_d0_pmu_emu_VSI","charge_DV_trk_d0_pel_ee_VSI"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI",
					savefilename = "trk_d0_MCdatacompare_VSI"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 2, 
					setrange= "0 20",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_dpt_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_trk_dpt_pmu_mumu_VSILep","charge_DV_trk_dpt_pmu_emu_VSILep","charge_DV_trk_dpt_pel_ee_VSILep"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI Leptons",
					savefilename = "trk_dpt_MCdatacompare_VSILep"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 2, 
					setrange= "0 20",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_dpt_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_trk_dpt_pmu_mumu_VSI","charge_DV_trk_dpt_pmu_emu_VSI","charge_DV_trk_dpt_pel_ee_VSI"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI",
					savefilename = "trk_dpt_MCdatacompare_VSI"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 2, 
					setrange= "0 4",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_deta_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_trk_deta_pmu_mumu_VSILep","charge_DV_trk_deta_pmu_emu_VSILep","charge_DV_trk_deta_pel_ee_VSILep"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI Leptons",
					savefilename = "trk_deta_MCdatacompare_VSILep"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 2, 
					setrange= "0 4",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_deta_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_trk_deta_pmu_mumu_VSI","charge_DV_trk_deta_pmu_emu_VSI","charge_DV_trk_deta_pel_ee_VSI"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI",
					savefilename = "trk_deta_MCdatacompare_VSI"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "0 8",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_dphi_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_trk_dphi_pmu_mumu_VSILep","charge_DV_trk_dphi_pmu_emu_VSILep","charge_DV_trk_dphi_pel_ee_VSILep"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI Leptons",
					savefilename = "trk_dphi_MCdatacompare_VSILep"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "0 8",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_dphi_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_trk_dphi_pmu_mumu_VSI","charge_DV_trk_dphi_pmu_emu_VSI","charge_DV_trk_dphi_pel_ee_VSI"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI",
					savefilename = "trk_dphi_MCdatacompare_VSI"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 10, 
					setrange= "0 2",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_dR_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_trk_dR_pmu_mumu_VSILep","charge_DV_trk_dR_pmu_emu_VSILep","charge_DV_trk_dR_pel_ee_VSILep"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI Leptons",
					savefilename = "trk_dR_MCdatacompare_VSILep"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 10, 
					setrange= "0 10",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_dR_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_trk_dR_pmu_mumu_VSI","charge_DV_trk_dR_pmu_emu_VSI","charge_DV_trk_dR_pel_ee_VSI"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI",
					savefilename = "trk_dR_MCdatacompare_VSI"	)


	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 5, 
					setrange= "0 200",
					scaleymax = 1.5,
					hdataname = "charge_mvis_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_mvis_pmu_mumu_VSILep","charge_mvis_pmu_emu_VSILep","charge_mvis_pel_ee_VSILep"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI Leptons",
					savefilename = "mvis_MCdatacompare_VSILep"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 5, 
					setrange= "0 200",
					scaleymax = 1.5,
					hdataname = "charge_mvis_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_mvis_pmu_mumu_VSI","charge_mvis_pmu_emu_VSI","charge_mvis_pel_ee_VSI"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI",
					savefilename = "mvis_MCdatacompare_VSI"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 5, 
					setrange= "0 200",
					scaleymax = 1.5,
					hdataname = "charge_mtrans_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_mtrans_pmu_mumu_VSILep","charge_mtrans_pmu_emu_VSILep","charge_mtrans_pel_ee_VSILep"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI Leptons",
					savefilename = "mtrans_MCdatacompare_VSILep"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 5, 
					setrange= "0 200",
					scaleymax = 1.5,
					hdataname = "charge_mtrans_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_mtrans_pmu_mumu_VSI","charge_mtrans_pmu_emu_VSI","charge_mtrans_pel_ee_VSI"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI",
					savefilename = "mtrans_MCdatacompare_VSI"	)



	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "0 20",
					scaleymax = 1.5,
					hdataname = "charge_DV_redmass_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_redmass_pmu_mumu_VSILep","charge_DV_redmass_pmu_emu_VSILep","charge_DV_redmass_pel_ee_VSILep"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI Leptons",
					savefilename = "DV_redmass_MCdatacompare_VSILep"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "0 20",
					scaleymax = 1.5,
					hdataname = "charge_DV_redmass_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_redmass_pmu_mumu_VSI","charge_DV_redmass_pmu_emu_VSI","charge_DV_redmass_pel_ee_VSI"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI",
					savefilename = "DV_redmass_MCdatacompare_VSI"	)


	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 5, 
					setrange= "0 200",
					scaleymax = 1.5,
					hdataname = "charge_DV_redmassvis_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_redmassvis_pmu_mumu_VSILep","charge_DV_redmassvis_pmu_emu_VSILep","charge_DV_redmassvis_pel_ee_VSILep"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI Leptons",
					savefilename = "DV_redmassvis_MCdatacompare_VSILep"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 5, 
					setrange= "0 200",
					scaleymax = 1.5,
					hdataname = "charge_DV_redmassvis_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_redmassvis_pmu_mumu_VSI","charge_DV_redmassvis_pmu_emu_VSI","charge_DV_redmassvis_pel_ee_VSI"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI",
					savefilename = "DV_redmassvis_MCdatacompare_VSI"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "0 20",
					scaleymax = 1.5,
					hdataname = "charge_HNLm_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_HNLm_pmu_mumu_VSILep","charge_HNLm_pmu_emu_VSILep","charge_HNLm_pel_ee_VSILep"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI Leptons",
					savefilename = "HNLm_MCdatacompare_VSILep"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "0 20",
					scaleymax = 1.5,
					hdataname = "charge_HNLm_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_HNLm_pmu_mumu_VSI","charge_HNLm_pmu_emu_VSI","charge_HNLm_pel_ee_VSI"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI",
					savefilename = "HNLm_MCdatacompare_VSI"	)


	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "0 50",
					scaleymax = 1.5,
					hdataname = "charge_DV_redmassHNL_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_redmassHNL_pmu_mumu_VSILep","charge_DV_redmassHNL_pmu_emu_VSILep","charge_DV_redmassHNL_pel_ee_VSILep"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI Leptons",
					savefilename = "redmassHNL_MCdatacompare_VSILep"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "0 50",
					scaleymax = 1.5,
					hdataname = "charge_DV_redmassHNL_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_DV_redmassHNL_pmu_mumu_VSI","charge_DV_redmassHNL_pmu_emu_VSI","charge_DV_redmassHNL_pel_ee_VSI"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI",
					savefilename = "redmassHNL_MCdatacompare_VSI"	)



	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 5, 
					setrange= "0 200",
					scaleymax = 1.5,
					hdataname = "charge_HNLpt_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_HNLpt_pmu_mumu_VSILep","charge_HNLpt_pmu_emu_VSILep","charge_HNLpt_pel_ee_VSILep"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI Leptons",
					savefilename = "HNLpt_MCdatacompare_VSILep"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 5, 
					setrange= "0 200",
					scaleymax = 1.5,
					hdataname = "charge_HNLpt_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_HNLpt_pmu_mumu_VSI","charge_HNLpt_pmu_emu_VSI","charge_HNLpt_pel_ee_VSI"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI",
					savefilename = "HNLpt_MCdatacompare_VSI"	)


	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "-4 4",
					scaleymax = 1.5,
					hdataname = "charge_HNLphi_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_HNLphi_pmu_mumu_VSILep","charge_HNLphi_pmu_emu_VSILep","charge_HNLphi_pel_ee_VSILep"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI Leptons",
					savefilename = "HNLphi_MCdatacompare_VSILep"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "-4 4",
					scaleymax = 1.5,
					hdataname = "charge_HNLphi_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_HNLphi_pmu_mumu_VSI","charge_HNLphi_pmu_emu_VSI","charge_HNLphi_pel_ee_VSI"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI",
					savefilename = "HNLphi_MCdatacompare_VSI"	)


	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "-3 3",
					scaleymax = 1.5,
					hdataname = "charge_HNLeta_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_HNLeta_pmu_mumu_VSILep","charge_HNLeta_pmu_emu_VSILep","charge_HNLeta_pel_ee_VSILep"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI Leptons",
					savefilename = "HNLeta_MCdatacompare_VSILep"	)

	compare_dataMC_multi(datafile=files[3],
						mcfiles=[files[0],files[1],files[2]], 
					nRebin = 1, 
					setrange= "-3 3",
					scaleymax = 1.5,
					hdataname = "charge_HNLeta_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcnames = ["charge_HNLeta_pmu_mumu_VSI","charge_HNLeta_pmu_emu_VSI","charge_HNLeta_pel_ee_VSI"],
					hmclabels = ["(m_{HNL}, c\\tau, chan.) = (10, 10, \mu\mu\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, \mue\mu)", "(m_{HNL}, c\\tau, chan.) = (10, 10, eee)"],
					vertextype = "VSI",
					savefilename = "HNLeta_MCdatacompare_VSI"	)


	# compareN(files[0],
	# 		setxrange= "0 100",
	# 		scaleymax=1.5,
	# 		nRebin=5,
	# 		# hname = ["sel_DV_mass_pmu_mumu_VSILep","sel_mvis_pmu_mumu_VSILep","sel_HNLm_pmu_mumu_VSILep","sel_mtrans_pmu_mumu_VSILep","sel_DV_redmass_pmu_mumu_VSILep"],
	# 		hname = ["sel_mvis_pmu_mumu_VSILep"],
	# 		# hlabel = ["DV mass","Visible mass","m_{HNL}", "m_{T}","Reduced mass"], 
	# 		hlabel = ["Visible mass"], 
	# 		vertextype= "VSI Leptons",
	# 		savefilename="selMC_massvis_VSILep" )


	compareN(files[3],
			setxrange= "0 100",
			scaleymax=1.6,
			nRebin=2,
			hname = ["charge_DV_mass_pmu_mumu_VSILep","charge_mvis_pmu_mumu_VSILep","charge_HNLm_pmu_mumu_VSILep","charge_mtrans_pmu_mumu_VSILep","charge_DV_redmass_pmu_mumu_VSILep"],
			hlabel = ["DV mass","Visible mass","m_{HNL}", "m_{T}","Reduced mass"], 
			vertextype= "VSI Leptons",
			savefilename="chargeData_compareMassDef_VSILep" )


	compareN(files[3],
			setxrange= "0 100",
			scaleymax=1.2,
			nRebin=2,
			hname = ["charge_DV_mass_pmu_mumu_VSI","charge_mvis_pmu_mumu_VSI","charge_HNLm_pmu_mumu_VSI","charge_mtrans_pmu_mumu_VSI","charge_DV_redmass_pmu_mumu_VSI"],
			hlabel = ["DV mass","Visible mass","m_{HNL}", "m_{T}","Reduced mass"], 
			vertextype= "VSI",
			savefilename="chargeData_compareMassDef_VSI" )



	exit()
	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "0 20",
					scaleymax = 1.2,				
					hdataname = "charge_DV_mass_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_mass_pmu_mumu_VSILep",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI Leptons",
					savefilename = "DV_mass_MCdatacompare_VSILep"	)


	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "0 40",
					scaleymax = 1.2,				
					hdataname = "charge_DV_mass_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_mass_pmu_mumu_VSI",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI",
					savefilename = "DV_mass_MCdatacompare_VSI"	)

	compare_dataMC(files[1],files[0],
					nRebin = 5, 
					setrange= "",
					scaleymax = 1.2,
					hdataname = "charge_DV_r_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_r_pmu_mumu_VSILep",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI Leptons",
					savefilename = "DV_r_MCdatacompare_VSILep"	)


	compare_dataMC(files[1],files[0],
					nRebin = 5, 
					setrange= "",
					scaleymax = 1.2,
					hdataname = "charge_DV_r_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_r_pmu_mumu_VSI",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI",
					savefilename = "DV_r_MCdatacompare_VSI"	)


	compare_dataMC(files[1],files[0],
					nRebin = 2, 
					setrange= "0 100",
					scaleymax = 1.2,					
					hdataname = "charge_DV_trk_pt_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_trk_pt_pmu_mumu_VSILep",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI Leptons",
					savefilename = "trk_pt_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 2, 
					setrange= "0 100",
					setlogy= False,
					scaleymax = 1.2,					
					hdataname = "charge_DV_trk_pt_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_trk_pt_pmu_mumu_VSI",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI",
					savefilename = "trk_pt_MCdatacompare_VSI"	)



	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "",
					scaleymax = 2.2,
					hdataname = "charge_DV_trk_phi_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_trk_phi_pmu_mumu_VSILep",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI Leptons",
					savefilename = "trk_phi_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "",
					scaleymax = 2.2,
					hdataname = "charge_DV_trk_phi_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_trk_phi_pmu_mumu_VSI",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI",
					savefilename = "trk_phi_MCdatacompare_VSI"	)


	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "",
					scaleymax = 2.2,
					hdataname = "charge_DV_trk_eta_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_trk_eta_pmu_mumu_VSILep",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI Leptons",
					savefilename = "trk_eta_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "",
					scaleymax = 2.2,
					hdataname = "charge_DV_trk_eta_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_trk_eta_pmu_mumu_VSI",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI",
					savefilename = "trk_eta_MCdatacompare_VSI"	)


	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "-10 10",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_d0_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_trk_d0_pmu_mumu_VSILep",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI Leptons",
					savefilename = "trk_d0_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "-10 10",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_d0_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_trk_d0_pmu_mumu_VSI",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI",
					savefilename = "trk_d0_MCdatacompare_VSI"	)


	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "-10 10",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_d0_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_trk_d0_pmu_mumu_VSILep",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI Leptons",
					savefilename = "trk_d0_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "-10 10",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_d0_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_trk_d0_pmu_mumu_VSI",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI",
					savefilename = "trk_d0_MCdatacompare_VSI"	)

	compare_dataMC(files[1],files[0],
					nRebin = 2, 
					setrange= "0 20",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_dpt_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_trk_dpt_pmu_mumu_VSILep",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI Leptons",
					savefilename = "trk_dpt_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 2, 
					setrange= "0 20",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_dpt_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_trk_dpt_pmu_mumu_VSI",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI",
					savefilename = "trk_dpt_MCdatacompare_VSI"	)

	compare_dataMC(files[1],files[0],
					nRebin = 2, 
					setrange= "0 4",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_deta_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_trk_deta_pmu_mumu_VSILep",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI Leptons",
					savefilename = "trk_deta_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 2, 
					setrange= "0 4",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_deta_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_trk_deta_pmu_mumu_VSI",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI",
					savefilename = "trk_deta_MCdatacompare_VSI"	)

	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "0 8",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_dphi_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_trk_dphi_pmu_mumu_VSILep",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI Leptons",
					savefilename = "trk_dphi_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "0 8",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_dphi_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_trk_dphi_pmu_mumu_VSI",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI",
					savefilename = "trk_dphi_MCdatacompare_VSI"	)

	compare_dataMC(files[1],files[0],
					nRebin = 10, 
					setrange= "0 2",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_dR_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_trk_dR_pmu_mumu_VSILep",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI Leptons",
					savefilename = "trk_dR_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 10, 
					setrange= "0 10",
					scaleymax = 1.2,
					hdataname = "charge_DV_trk_dR_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_trk_dR_pmu_mumu_VSI",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI",
					savefilename = "trk_dR_MCdatacompare_VSI"	)


	compare_dataMC(files[1],files[0],
					nRebin = 5, 
					setrange= "0 200",
					scaleymax = 1.5,
					hdataname = "charge_mvis_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_mvis_pmu_mumu_VSILep",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI Leptons",
					savefilename = "mvis_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 5, 
					setrange= "0 200",
					scaleymax = 1.5,
					hdataname = "charge_mvis_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_mvis_pmu_mumu_VSI",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI",
					savefilename = "mvis_MCdatacompare_VSI"	)

	compare_dataMC(files[1],files[0],
					nRebin = 5, 
					setrange= "0 200",
					scaleymax = 1.5,
					hdataname = "charge_mtrans_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_mtrans_pmu_mumu_VSILep",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI Leptons",
					savefilename = "mtrans_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 5, 
					setrange= "0 200",
					scaleymax = 1.5,
					hdataname = "charge_mtrans_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_mtrans_pmu_mumu_VSI",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI",
					savefilename = "mtrans_MCdatacompare_VSI"	)



	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "0 20",
					scaleymax = 1.5,
					hdataname = "charge_DV_redmass_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_redmass_pmu_mumu_VSILep",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI Leptons",
					savefilename = "DV_redmass_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "0 20",
					scaleymax = 1.5,
					hdataname = "charge_DV_redmass_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_redmass_pmu_mumu_VSI",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI",
					savefilename = "DV_redmass_MCdatacompare_VSI"	)


	compare_dataMC(files[1],files[0],
					nRebin = 5, 
					setrange= "0 200",
					scaleymax = 1.5,
					hdataname = "charge_DV_redmassvis_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_redmassvis_pmu_mumu_VSILep",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI Leptons",
					savefilename = "DV_redmassvis_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 5, 
					setrange= "0 200",
					scaleymax = 1.5,
					hdataname = "charge_DV_redmassvis_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_redmassvis_pmu_mumu_VSI",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI",
					savefilename = "DV_redmassvis_MCdatacompare_VSI"	)

	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "0 20",
					scaleymax = 1.5,
					hdataname = "charge_HNLm_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_HNLm_pmu_mumu_VSILep",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI Leptons",
					savefilename = "HNLm_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "0 20",
					scaleymax = 1.5,
					hdataname = "charge_HNLm_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_HNLm_pmu_mumu_VSI",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI",
					savefilename = "HNLm_MCdatacompare_VSI"	)


	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "0 50",
					scaleymax = 1.5,
					hdataname = "charge_DV_redmassHNL_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_redmassHNL_pmu_mumu_VSILep",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI Leptons",
					savefilename = "redmassHNL_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "0 0",
					scaleymax = 1.5,
					hdataname = "charge_DV_redmassHNL_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_DV_redmassHNL_pmu_mumu_VSI",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI",
					savefilename = "redmassHNL_MCdatacompare_VSI"	)



	compare_dataMC(files[1],files[0],
					nRebin = 5, 
					setrange= "0 200",
					scaleymax = 1.5,
					hdataname = "charge_HNLpt_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_HNLpt_pmu_mumu_VSILep",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI Leptons",
					savefilename = "HNLpt_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 5, 
					setrange= "0 200",
					scaleymax = 1.5,
					hdataname = "charge_HNLpt_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_HNLpt_pmu_mumu_VSI",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI",
					savefilename = "HNLpt_MCdatacompare_VSI"	)


	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "-4 4",
					scaleymax = 1.5,
					hdataname = "charge_HNLphi_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_HNLphi_pmu_mumu_VSILep",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI Leptons",
					savefilename = "HNLphi_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "-4 4",
					scaleymax = 1.5,
					hdataname = "charge_HNLphi_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_HNLphi_pmu_mumu_VSI",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI",
					savefilename = "HNLphi_MCdatacompare_VSI"	)


	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "-3 3",
					scaleymax = 1.5,
					hdataname = "charge_HNLeta_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_HNLeta_pmu_mumu_VSILep",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI Leptons",
					savefilename = "HNLeta_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "-3 3",
					scaleymax = 1.5,
					hdataname = "charge_HNLeta_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "charge_HNLeta_pmu_mumu_VSI",
					hmclabel = "(m_{HNL}, c\\tau) = (10, 10)",
					vertextype = "VSI",
					savefilename = "HNLeta_MCdatacompare_VSI"	)




	# massplot2D(files[0], 
	# 		hname="charge_DVmass_mvis_pmu_mumu_VSI", 
	# 		hlabel="(m_{HNL}, c\\tau) = (10, 10)",
	# 		vertextype="VSI",
	# 		savefilename="MC_charge_DVmass_mvis_VSI",
	# 		setxrange="0 30",
	# 		setyrange="0 300")

	# massplot2D(files[0], 
	# 		hname="charge_DVmass_mvis_pmu_mumu_VSILep", 
	# 		hlabel="(m_{HNL}, c\\tau) = (10, 10)",
	# 		vertextype="VSI Leptons",
	# 		savefilename="MC_charge_DVmass_mvis_VSILep",
	# 		setxrange="0 30",
	# 		setyrange="0 300")

	# massplot2D(files[0], 
	# 		rebinx=2,
	# 		rebiny=5,
	# 		hname="charge_DVmass_mtrans_pmu_mumu_VSI", 
	# 		hlabel="(m_{HNL}, c\\tau) = (10, 10)",
	# 		vertextype="VSI",
	# 		savefilename="MC_charge_DVmass_mtrans_VSI",
	# 		setxrange="0 30",
	# 		setyrange="0 300")

	# massplot2D(files[0], 
	# 		rebinx=2,
	# 		rebiny=5,
	# 		hname="charge_DVmass_mtrans_pmu_mumu_VSILep", 
	# 		hlabel="(m_{HNL}, c\\tau) = (10, 10)",
	# 		vertextype="VSI Leptons",
	# 		savefilename="MC_charge_DVmass_mtrans_VSILep",
	# 		setxrange="0 30",
	# 		setyrange="0 300")

	# massplot2D(files[0], 
	# 		rebinx=5,
	# 		rebiny=5,
	# 		hname="charge_mvis_mtrans_pmu_mumu_VSI", 
	# 		hlabel="(m_{HNL}, c\\tau) = (10, 10)",
	# 		vertextype="VSI",
	# 		savefilename="MC_charge_mvis_mtrans_VSI",
	# 		setxrange="0 100",
	# 		setyrange="0 300")

	# massplot2D(files[0], 
	# 		rebinx=5,
	# 		rebiny=5,
	# 		hname="charge_mvis_mtrans_pmu_mumu_VSILep", 
	# 		hlabel="(m_{HNL}, c\\tau) = (10, 10)",
	# 		vertextype="VSI Leptons",
	# 		savefilename="MC_charge_mvis_mtrans_VSILep",
	# 		setxrange="0 100",
	# 		setyrange="0 300")

	# massplot2D(files[0], 
	# 		hname="charge_DVmass_mhnl_pmu_mumu_VSI", 
	# 		hlabel="(m_{HNL}, c\\tau) = (10, 10)",
	# 		vertextype="VSI",
	# 		savefilename="MC_charge_DVmass_mhnl_VSI",
	# 		setxrange="0 50",
	# 		setyrange="0 50")

	# massplot2D(files[0], 
	# 		hname="charge_DVmass_mhnl_pmu_mumu_VSILep", 
	# 		hlabel="(m_{HNL}, c\\tau) = (10, 10)",
	# 		vertextype="VSI Leptons",
	# 		savefilename="MC_charge_DVmass_mhnl_VSILep",
	# 		setxrange="0 50",
	# 		setyrange="0 50")


	# massplot2D(files[0], 
	# 		rebinx=5,
	# 		rebiny=1,
	# 		hname="charge_mvis_mhnl_pmu_mumu_VSI", 
	# 		hlabel="(m_{HNL}, c\\tau) = (10, 10)",
	# 		vertextype="VSI",
	# 		savefilename="MC_charge_mvis_mhnl_VSI",
	# 		setxrange="0 300",
	# 		setyrange="0 50")

	# massplot2D(files[0], 
	# 		rebinx=5,
	# 		rebiny=1,
	# 		hname="charge_mvis_mhnl_pmu_mumu_VSILep", 
	# 		hlabel="(m_{HNL}, c\\tau) = (10, 10)",
	# 		vertextype="VSI Leptons",
	# 		savefilename="MC_charge_mvis_mhnl_VSILep",
	# 		setxrange="0 300",
	# 		setyrange="0 50")


	# massplot2D(files[0], 
	# 		rebinx=1,
	# 		rebiny=5,
	# 		hname="charge_mhnl_mtrans_pmu_mumu_VSI", 
	# 		hlabel="(m_{HNL}, c\\tau) = (10, 10)",
	# 		vertextype="VSI",
	# 		savefilename="MC_charge_mhnl_mtrans_VSI",
	# 		setxrange="0 50",
	# 		setyrange="0 300")

	# massplot2D(files[0], 
	# 		rebinx=1,
	# 		rebiny=5,
	# 		hname="charge_mhnl_mtrans_pmu_mumu_VSILep", 
	# 		hlabel="(m_{HNL}, c\\tau) = (10, 10)",
	# 		vertextype="VSI Leptons",
	# 		savefilename="MC_charge_mhnl_mtrans_VSILep",
	# 		setxrange="0 50",
	# 		setyrange="0 300")




 #   # data 2D mass 
	# massplot2D(files[1], 
	# 		hname="charge_DVmass_mvis_pmu_mumu_VSI", 
	# 		hlabel="data 2018 period B",
	# 		vertextype="VSI",
	# 		savefilename="Data_charge_DVmass_mvis_VSI",
	# 		setxrange="0 30",
	# 		setyrange="0 300")

	# massplot2D(files[1], 
	# 		hname="charge_DVmass_mvis_pmu_mumu_VSILep", 
	# 		hlabel="data 2018 period B",
	# 		vertextype="VSI Leptons",
	# 		savefilename="Data_charge_DVmass_mvis_VSILep",
	# 		setxrange="0 30",
	# 		setyrange="0 300")

	# massplot2D(files[1], 
	# 		rebinx=2,
	# 		rebiny=5,
	# 		hname="charge_DVmass_mtrans_pmu_mumu_VSI", 
	# 		hlabel="data 2018 period B",
	# 		vertextype="VSI",
	# 		savefilename="Data_charge_DVmass_mtrans_VSI",
	# 		setxrange="0 30",
	# 		setyrange="0 300")

	# massplot2D(files[1], 
	# 		rebinx=2,
	# 		rebiny=5,
	# 		hname="charge_DVmass_mtrans_pmu_mumu_VSILep", 
	# 		hlabel="data 2018 period B",
	# 		vertextype="VSI Leptons",
	# 		savefilename="Data_charge_DVmass_mtrans_VSILep",
	# 		setxrange="0 30",
	# 		setyrange="0 300")

	# massplot2D(files[1], 
	# 		rebinx=5,
	# 		rebiny=5,
	# 		hname="charge_mvis_mtrans_pmu_mumu_VSI", 
	# 		hlabel="data 2018 period B",
	# 		vertextype="VSI",
	# 		savefilename="Data_charge_mvis_mtrans_VSI",
	# 		setxrange="0 100",
	# 		setyrange="0 300")

	# massplot2D(files[1], 
	# 		rebinx=5,
	# 		rebiny=5,
	# 		hname="charge_mvis_mtrans_pmu_mumu_VSILep", 
	# 		hlabel="data 2018 period B",
	# 		vertextype="VSI Leptons",
	# 		savefilename="Data_charge_mvis_mtrans_VSILep",
	# 		setxrange="0 100",
	# 		setyrange="0 300")

	# massplot2D(files[1], 
	# 		hname="charge_DVmass_mhnl_pmu_mumu_VSI", 
	# 		hlabel="data 2018 period B",
	# 		vertextype="VSI",
	# 		savefilename="data_charge_DVmass_mhnl_VSI",
	# 		setxrange="0 50",
	# 		setyrange="0 50")

	# massplot2D(files[1], 
	# 		hname="charge_DVmass_mhnl_pmu_mumu_VSILep", 
	# 		hlabel="data 2018 period B",
	# 		vertextype="VSI Leptons",
	# 		savefilename="data_charge_DVmass_mhnl_VSILep",
	# 		setxrange="0 50",
	# 		setyrange="0 50")


	# massplot2D(files[1], 
	# 		rebinx=5,
	# 		rebiny=1,
	# 		hname="charge_mvis_mhnl_pmu_mumu_VSI", 
	# 		hlabel="data 2018 period B",
	# 		vertextype="VSI",
	# 		savefilename="data_charge_mvis_mhnl_VSI",
	# 		setxrange="0 300",
	# 		setyrange="0 50")

	# massplot2D(files[1], 
	# 		rebinx=5,
	# 		rebiny=1,
	# 		hname="charge_mvis_mhnl_pmu_mumu_VSILep", 
	# 		hlabel="data 2018 period B",
	# 		vertextype="VSI Leptons",
	# 		savefilename="data_charge_mvis_mhnl_VSILep",
	# 		setxrange="0 300",
	# 		setyrange="0 50")


	# massplot2D(files[1], 
	# 		rebinx=1,
	# 		rebiny=5,
	# 		hname="charge_mhnl_mtrans_pmu_mumu_VSI", 
	# 		hlabel="data 2018 period B",
	# 		vertextype="VSI",
	# 		savefilename="data_charge_mhnl_mtrans_VSI",
	# 		setxrange="0 50",
	# 		setyrange="0 300")

	# massplot2D(files[1], 
	# 		rebinx=1,
	# 		rebiny=5,
	# 		hname="charge_mhnl_mtrans_pmu_mumu_VSILep", 
	# 		hlabel="data 2018 period B",
	# 		vertextype="VSI Leptons",
	# 		savefilename="data_charge_mhnl_mtrans_VSILep",
	# 		setxrange="0 50",
	# 		setyrange="0 300")





	# compareN(files[0],
	# 		setxrange= "0 100",
	# 		scaleymax=1.2,
	# 		nRebin=2,
	# 		hname = ["sel_DV_mass_pmu_mumu_VSILep","sel_mvis_pmu_mumu_VSILep","sel_HNLm_pmu_mumu_VSILep","sel_mtrans_pmu_mumu_VSILep","sel_DV_redmass_pmu_mumu_VSILep"],
	# 		hlabel = ["DV mass","Visible mass","m_{HNL}", "m_{T}","Reduced mass"], 
	# 		vertextype= "VSI Leptons",
	# 		savefilename="selMC_compareMassDef_VSILep" )


	# compareN(files[0],
	# 		setxrange= "0 100",
	# 		scaleymax=1.2,
	# 		nRebin=2,
	# 		hname = ["sel_DV_mass_pmu_mumu_VSI","sel_mvis_pmu_mumu_VSI","sel_HNLm_pmu_mumu_VSI","sel_mtrans_pmu_mumu_VSI","sel_DV_redmass_pmu_mumu_VSI"],
	# 		hlabel = ["DV mass","Visible mass","m_{HNL}", "m_{T}","Reduced mass"], 
	# 		vertextype= "VSI",
	# 		savefilename="selMC_compareMassDef_VSI" )

	# compareN(files[0],
	# 		setxrange= "0 100",
	# 		scaleymax=1.6,
	# 		nRebin=2,
	# 		hname = ["charge_DV_mass_pmu_mumu_VSILep","charge_mvis_pmu_mumu_VSILep","charge_HNLm_pmu_mumu_VSILep","charge_mtrans_pmu_mumu_VSILep","charge_DV_redmass_pmu_mumu_VSILep"],
	# 		hlabel = ["DV mass","Visible mass","m_{HNL}", "m_{T}","Reduced mass"], 
	# 		vertextype= "VSI Leptons",
	# 		savefilename="chargeMC_compareMassDef_VSILep" )


	# compareN(files[0],
	# 		setxrange= "0 100",
	# 		scaleymax=1.2,
	# 		nRebin=2,
	# 		hname = ["charge_DV_mass_pmu_mumu_VSI","charge_mvis_pmu_mumu_VSI","charge_HNLm_pmu_mumu_VSI","charge_mtrans_pmu_mumu_VSI","charge_DV_redmass_pmu_mumu_VSI"],
	# 		hlabel = ["DV mass","Visible mass","m_{HNL}", "m_{T}","Reduced mass"], 
	# 		vertextype= "VSI",
	# 		savefilename="chargeMC_compareMassDef_VSI" )


	# compareN(files[1],
	# 		setxrange= "0 100",
	# 		scaleymax=1.6,
	# 		nRebin=2,
	# 		hname = ["charge_DV_mass_pmu_mumu_VSILep","charge_mvis_pmu_mumu_VSILep","charge_HNLm_pmu_mumu_VSILep","charge_mtrans_pmu_mumu_VSILep","charge_DV_redmass_pmu_mumu_VSILep"],
	# 		hlabel = ["DV mass","Visible mass","m_{HNL}", "m_{T}","Reduced mass"], 
	# 		vertextype= "VSI Leptons",
	# 		savefilename="chargeData_compareMassDef_VSILep" )


	# compareN(files[1],
	# 		setxrange= "0 100",
	# 		scaleymax=1.2,
	# 		nRebin=2,
	# 		hname = ["charge_DV_mass_pmu_mumu_VSI","charge_mvis_pmu_mumu_VSI","charge_HNLm_pmu_mumu_VSI","charge_mtrans_pmu_mumu_VSI","charge_DV_redmass_pmu_mumu_VSI"],
	# 		hlabel = ["DV mass","Visible mass","m_{HNL}", "m_{T}","Reduced mass"], 
	# 		vertextype= "VSI",
	# 		savefilename="chargeData_compareMassDef_VSI" )





	exit()
# def compareN(file, hname, hlabel,savefilename):

	# for k, configs in config_file.items():
	# 	for key in config_file[k]["channels"]: 
	# 		print key
	# 	# print config_file[k]["channels"].keys()


	# compareN(options.file[0],["DVtypeDV_r_pel_ee_VSILep","DVtypeDV_r_pmu_ee_VSILep","DVtypeDV_r_pel_mumu_VSILep"],["eee","muee","mumu"], savefilename="testfile_MC" )
	####################################################################################################################################
	# Here's where you configure what histograms to plot

	# plot_cutflow(options.file[0], ch_name = "emu",
	# 							  vertextype = "Run1",
	# 							  savefilename = "master")


	# channel_VSI = ["pel_ee_VSI","pel_emu_VSI","pel_mumu_VSI","pmu_ee_VSI","pmu_emu_VSI","pmu_mumu_VSI"]
	# channel_VSILep = ["pel_ee_VSILep","pel_emu_VSILep","pel_mumu_VSILep","pmu_ee_VSILep","pmu_emu_VSILep","pmu_mumu_VSILep"]

	# for i in range(len(channel_VSI)):
	# 	plot_cutflow(options.file[0], ch_name =channel_VSI[i],
	# 							  vertextype = "VSI",
	# 							  savefilename = "selLep_periodB_%s"%channel_VSI[i])

	# 	plot_cutflow(options.file[0], ch_name =channel_VSILep[i],
	# 							  vertextype = "VSI Leptons",
	# 							  savefilename = "selLep_periodB_%s"%channel_VSILep[i])


	# plot_cutflow(options.file[0], ch_name ="charge",
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

	# compare4(options.file[0], nRebin = 2,
	# 						  h1name="selDV_mass_pel_emu_VSI",
	# 						  h1label="pel emu", 
	# 						  h2name="selDV_mass_pel_mumu_VSI",
	# 						  h2label="pel mumu", 
	# 						  h3name="selDV_mass_pmu_emu_VSI",
	# 						  h3label="pmu emu", 
	# 						  h4name="selDV_mass_pmu_mumu_VSI",
	# 						  h4label="pmu mumu", 
							 
	# 						  vertextype="VSI",
	# 						  xlabel='r DV [mm]',
	# 						  savefilename='hmassDV_up2toLeptons_VSI')



	# compare4(options.file[0], nRebin = 2,
	# 						  h1name="selDV_r_pel_emu_VSILep",
	# 						  h1label="pel emu", 
	# 						  h2name="selDV_r_pel_mumu_VSILep",
	# 						  h2label="pel mumu", 
	# 						  h3name="selDV_r_pmu_emu_VSILep",
	# 						  h3label="pmu emu", 
	# 						  h4name="selDV_r_pmu_mumu_VSILep",
	# 						  h4label="pmu mumu", 
							 
	# 						  vertextype="VSI",
	# 						  xlabel='r DV [mm]',
	# 						  savefilename='hrDV_up2toLeptons_VSILeptons')



	# compare2_wRatio(options.file[0],options.file2[0], h1name="DV_r_SS", 
	# 												  h1label="VSI", 	
	# 												  h2name="DV_r_SS",
	# 												  h2label="VSI Leptons", 
	# 												  xlabel='r DV [mm]',
	# 												  savefilename='hrDV_compare2deri')

	####################################################################################################################################








