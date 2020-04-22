# plotting functions
import argparse, os, math, ROOT, glob, uproot, time, json

import atlas_style
import numpy as np
import helpers
import plotting_helpers
from ROOT import *
from ROOT import gPad
from pylab import *


logger = helpers.getLogger('dHNLAnalysis.plotting')


def plot_cutflow(file, vertextype, outputDir="../output/"):

	Tfile = ROOT.TFile(file)
	if vertextype == "VSI":
		hcutflow = Tfile.Get('CutFlow_VSI')
	elif vertextype == "VSI Leptons":
		hcutflow = Tfile.Get('CutFlow_VSI_Leptons')

	MyC01= ROOT.TCanvas("MyC01","cutflow",1200,400)

	gPad.SetLogy()
	ymax_cutflow = hcutflow.GetMaximum()
	hcutflow.GetYaxis().SetRangeUser(0.1,ymax_cutflow*1000000000000000)
	hcutflow.SetFillColor(kAzure-4)
	hcutflow.SetLineWidth(0)
	hcutflow.GetYaxis().SetTickLength(0.)
	hcutflow.SetMarkerSize(2.2)
	hcutflow.Draw("HIST TEXT45 SAME")


	channel = file.split("histograms_")[1].split(".")[0]
	print channel
	print file
	fileInfo = helpers.File_info(file,channel)
	if "data" in file: 
		plotting_helpers.drawNotesData("data 2018",vertextype,lumi=60,channel=channel) 
	else: 
		plotting_helpers.drawNotesMC(fileInfo.MC_campaign,vertextype,channel,fileInfo.mass_str,fileInfo.ctau_str) 
		# MC_campaign,Vertextype, channel,mass,lifetime
	# elif "uuu" in file: 
	# 	# plotting_helpers.drawNotes("mumu","muon",vertextype) 
	# 	plotting_helpers.drawNotesMC("",vertextype,"uuu", "10", "1" )
	# elif "ueu" in file: 
	# 	plotting_helpers.drawNotes("emu","muon",vertextype) 
	# elif "uee" in file: 
	# 	plotting_helpers.drawNotes("ee","muon",vertextype) 
	# elif "eee" in file: 
	# 	plotting_helpers.drawNotes("ee","electron",vertextype) 
	# elif "eeu" in file: 
	# 	plotting_helpers.drawNotes("emu","electron",vertextype) 
	# elif "euu" in file: 
	# 	plotting_helpers.drawNotes("mumu","electron",vertextype) 


	if vertextype == "VSI":
		savefilename= "CutFlow_VSI_" + channel
	elif vertextype == "VSI Leptons":
		savefilename= "CutFlow_VSI_Leptons_" + channel

	MyC01.SaveAs(outputDir +savefilename+'.pdf')





def compareN(file, hname, hlabel,savefilename,vertextype,setxrange="",scaleymax=1,nRebin=1,setlogy=False,outputDir="../output/",normalize=False,lumi=1):
	############################################################################
	# nRebin = 5 # set to 1 if you dont want to rebin.
	# scaleymax = 1.6 # use this to scale height of y axis for asthetics
	############################################################################


	# get multiple histograms from input file
	nhist = len(hname)
	f = ROOT.TFile(file)
	
	h = {}
	if vertextype == "VSI":
		vertexSuffix = "_VSI"
	elif vertextype == "VSI Leptons":
		vertexSuffix = "_VSI_Leptons"

	for i in xrange(nhist): 
		h[i] = f.Get(hname[i]+vertexSuffix)
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


	h_binxmax_list = []
	h_binxmin_list = []
	for i in xrange(nhist):
		leg01.AddEntry(h[i],hlabel[i],"lp")
		h[i].Rebin(nRebin)
		h_binxmax_list.append(h[i].FindLastBinAbove(0,1))
		h_binxmin_list.append(h[i].FindFirstBinAbove(0,1))

	bin_xmax = max(h_binxmax_list)
	bin_xmin = min(h_binxmin_list)

	shapelist = [20,22,21,33,29]

	# normalize the histograms
	if normalize: 
		norm = 1
		for i in range(nhist):
			if (h[i].Integral() != 0):
				scale_mc = norm/(h[i].Integral(bin_xmin, bin_xmax))
			else:
				scale_mc = norm
			h[i].Scale(scale_mc)

	else: # use MC scaling to compare with data
		for i in range(nhist):
			h[i].Scale(lumi) # scale mc to a given luminosity in inverse fb

	# find the common min and max for y axis
	y_max = h[0].GetMaximum()
	for i in xrange(nhist):
		if h[i].GetMaximum() > y_max: y_max = h[i].GetMaximum()

	for i in xrange(nhist): 
		h[i].SetLineColor(plotting_helpers.histColours(i))
		h[i].SetMarkerColor(plotting_helpers.histColours(i))
		h[i].SetMarkerSize(0.7)
		h[i].SetMarkerStyle(shapelist[i])
		if i == 0:
			h[i].GetXaxis().SetTitle(plotting_helpers.xlabelhistograms(hname[i]))
			h[i].GetYaxis().SetTitle("entries")
			h[i].GetYaxis().SetRangeUser(0.0001,y_max*scaleymax)
			if setxrange != "":
				X_minmax = [item for item in setxrange.split(' ')]
				h[i].GetXaxis().SetRangeUser(int(X_minmax[0]),int(X_minmax[1]))
			else:
				h[i].GetXaxis().SetRange(bin_xmin-1,bin_xmax+1)
			h[i].Draw("E0 HIST")
		else: 
			h[i].Draw("E0 HIST SAME")
	
	if setlogy: 
		gPad.SetLogy()



	leg01.Draw()
	atlas_style.ATLASLabel(0.25,0.87,"Internal")

	# if "pmu" in hname[0]: 
	# 	if "mumu" in hname[0]:
	# 		plotting_helpers.drawNotesMC("",vertextype, "mumumu","10","10")
	# 	# if "mumu" in hname[0]:
	# 	# 	plotting_helpers.drawNotesMC("",vertextype, "mumu","10","10")
	# 	if "emu" in hname[0]:
	# 		plotting_helpers.drawNotesMC("",vertextype, "emumu","10","10")

	if "data" in file: 
		plotting_helpers.drawNotesData("data18 period B",vertextype) 
	if "uuu" in file: 
		plotting_helpers.drawNotes("mumu","muon",vertextype) 
	elif "ueu" in file: 
		plotting_helpers.drawNotes("emu","muon",vertextype) 
	elif "uee" in file: 
		plotting_helpers.drawNotes("ee","muon",vertextype) 
	elif "eee" in file: 
		plotting_helpers.drawNotes("ee","electron",vertextype) 
	elif "eeu" in file: 
		plotting_helpers.drawNotes("emu","electron",vertextype) 
	elif "euu" in file: 
		plotting_helpers.drawNotes("mumu","electron",vertextype) 

	if "HNLpt" in hname[0]:
		min_HNLptcut=TLine(20,0,20,y_max)
		min_HNLptcut.SetLineStyle(3)
		min_HNLptcut.Draw("SAME")

		max_HNLptcut=TLine(60,0,60,y_max)
		max_HNLptcut.SetLineStyle(3)
		max_HNLptcut.Draw("SAME")

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
				matlay[i]=TLine(matlayers[i],0,matlayers[i],y_max)
				matlay[i].SetLineStyle(3)
				matlay[i].Draw("SAME")
			leg01.AddEntry(matlay[0],"material layers","l")


	channel = file.split("_")[1].split(".")[0]

	if vertextype == "VSI":
		name = savefilename +"_"+ channel
	elif vertextype == "VSI Leptons":
		name= savefilename +"_"+ channel

	
	MyC01.SaveAs(outputDir +savefilename+'.pdf')





def compare_dataMC(datafile, mcfiles, hname, hdatalabel, hmclabels, vertextype, setrange = "", 
                   scaleymax=1.2, nRebin=1, setlogy=False, outputDir="../output/", save_name = "",
                   normalize=True, lumi=1):
	
	if setlogy: 
		scaleymax = scaleymax*150 #change default scale if log scale to give room for text & legend
	
	# get data histograms from input file
	
	data = ROOT.TFile(datafile)

	if vertextype == "VSI" and data.GetListOfKeys().Contains(hname + "_VSI"):
		hdata = data.Get(hname + "_VSI") 
	elif vertextype == "VSI Leptons" and data.GetListOfKeys().Contains(hname + "_VSI_Leptons"): 
		hdata = data.Get(hname + "_VSI_Leptons") 
	else: 
		logger.error("Couldn't find histogram: %s with vertextype: %s that you requested!"%(hname, vertextype ) ) 
		logger.error("Check that file %s has the histogram you are looking for!"%datafile)
		exit() 
	
	nmc_files = len(mcfiles)
	hmc ={}
	mc = {}

	# get MC histograms from input file

	for i in range(nmc_files): 
		mc[i] = ROOT.TFile(mcfiles[i])
		if vertextype == "VSI"and mc[i].GetListOfKeys().Contains(hname + "_VSI"):
			hmc[i] = mc[i].Get(hname + "_VSI")
		elif vertextype == "VSI Leptons" and mc[i].GetListOfKeys().Contains(hname + "_VSI_Leptons"):
			hmc[i] = mc[i].Get(hname + "_VSI_Leptons")
		else: 
			logger.error("Couldn't find histogram: %s with vertextype: %s that you requested!"%(hname, vertextype ))
			logger.error("Check that file %s has the histogram you are looking for!"%mcfiles[i])
			exit() 


	#define your canvas
	MyC01= ROOT.TCanvas("MyC01","",600,400)

	#format legend
	leg01 = ROOT.TLegend(0.50,0.7,0.91,0.92)
	leg01.SetTextSize(0.035)
	leg01.SetBorderSize(0)
	leg01.SetFillColor(kWhite)
	leg01.SetShadowColor(kWhite)

	# add data entry to legend 
	leg01.AddEntry(hdata,hdatalabel,"lp")

	#rebin histograms (if rebin =1 (default) it will not rebin the histogram)
	hdata.Rebin(nRebin)

	for i in range(nmc_files): 
		leg01.AddEntry(hmc[i],hmclabels[i],"lp")
		hmc[i].Rebin(nRebin)





	# find the common min and max for x axis
	bin_xmax = hdata.FindLastBinAbove(0,1)
	bin_xmin = hdata.FindFirstBinAbove(0,1)
	for i in xrange(nmc_files):
		if hmc[i].FindLastBinAbove(0,1) > bin_xmax: bin_xmax = hmc[i].FindLastBinAbove(0,1)
		if hmc[i].FindFirstBinAbove(0,1) > bin_xmin: bin_xmin = hmc[i].FindFirstBinAbove(0,1)
	if setrange == "":
		x_min = bin_xmin-1
		x_max = bin_xmax+1
	else: 
		x_min, x_max = [float(item) for item in setrange.split(' ')]
		bin_xmax = hdata.GetXaxis().FindBin(x_max)
		bin_xmin = hdata.GetXaxis().FindBin(x_min)
	hdata.GetXaxis().SetRangeUser(x_min, x_max)
	for i in xrange(nmc_files):
		hmc[i].GetXaxis().SetRangeUser(x_min, x_max)

	# normalize the histograms
	if normalize: 
		norm = 1
		if (hdata.Integral() != 0):
			# Normalize in specified range
			scale_data = norm/(hdata.Integral(bin_xmin, bin_xmax)) 
		else:
			scale_data = norm
		hdata.Scale(scale_data)

		for i in range(nmc_files):
			if (hmc[i].Integral() != 0):
				scale_mc = norm/(hmc[i].Integral(bin_xmin, bin_xmax))
			else:
				scale_mc = norm
			hmc[i].Scale(scale_mc)

	else: # use MC scaling to compare with data
		for i in range(nmc_files):
			hmc[i].Scale(lumi) # scale mc to a given luminosity in inverse fb


	# calculate some things for cuts
	if "DV_mass" in hname:
		b = hdata.Integral(hdata.FindFixBin(4), hdata.FindFixBin(20), "")
		print "background ", b
		for i in xrange(nmc_files):
			print hmc[i].FindFixBin(4)
			print hmc[i].FindFixBin(20)
			s = hmc[i].Integral(hmc[i].FindFixBin(4), hmc[i].FindFixBin(20), "")
			print "signal: ",  round(s,5)

			print "s/sqrt(b) ", round(s/np.sqrt(b), 5)


	if "HNLpt" in hname:
		b = hdata.Integral(hdata.FindFixBin(25), hdata.FindFixBin(50), "")
		print "background ", b
		for i in xrange(nmc_files):
			s = hmc[i].Integral(hmc[i].FindFixBin(25), hmc[i].FindFixBin(50), "")
			print "signal: ",  round(s,5)

			print "s/sqrt(b) ", round(s/np.sqrt(b), 5)


	# find the common min and max for y axis
	y_max = hdata.GetMaximum()
	for i in xrange(nmc_files):
		if hmc[i].GetMaximum() > y_max: y_max = hmc[i].GetMaximum()


	hdata.SetLineColor(kBlack)
	hdata.GetXaxis().SetTitle(plotting_helpers.xlabelhistograms(hname))
	hdata.GetYaxis().SetTitle("entries")
	hdata.GetYaxis().SetRangeUser(0.0001,y_max*scaleymax)
	hdata.SetMarkerSize(0.7)

	hdata.Draw("E0 HIST")

	shapelist = [22,21,33,29]
	for i in xrange(nmc_files): 
		hmc[i].SetMarkerSize(0.7)
		hmc[i].SetLineColor(plotting_helpers.histColours(i))
		hmc[i].SetMarkerColor(plotting_helpers.histColours(i))
		hmc[i].SetMarkerStyle(shapelist[i])
		hmc[i].GetYaxis().SetRangeUser(0.0001,y_max*scaleymax)
		hmc[i].Draw("E0 HIST SAME")


	if setlogy: 
		gPad.SetLogy()
		gPad.Update()

	# draw some lines 
	if "DV_mass" in hname:
		mDVcut=TLine(4,0,4,y_max)
		mDVcut.SetLineStyle(3)
		mDVcut.Draw("SAME")


	if "HNLpt" in hname:
		min_HNLptcut=TLine(20,0,20,y_max)
		min_HNLptcut.SetLineStyle(3)
		min_HNLptcut.Draw("SAME")

		max_HNLptcut=TLine(60,0,60,y_max)
		max_HNLptcut.SetLineStyle(3)
		max_HNLptcut.Draw("SAME")

	if "mvis" in hname:
		min_mlll=TLine(50,0,50,y_max)
		min_mlll.SetLineStyle(3)
		min_mlll.Draw("SAME")

		max_mlll=TLine(84,0,84,y_max)
		max_mlll.SetLineStyle(3)
		max_mlll.Draw("SAME")

	if "DV_r" in hname:
		if  "redmass" in hname:
			pass
		else: 
			matlayers = [33.25,50.5,88.5,122.5,299]
			nmatlayers = len(matlayers)
			matlay={}
			for i in range(nmatlayers):
				matlay[i]=TLine(matlayers[i],0,matlayers[i],y_max)
				matlay[i].SetLineStyle(3)
				matlay[i].Draw("SAME")
			leg01.AddEntry(matlay[0],"material layers","l")

	leg01.Draw()
	atlas_style.ATLASLabel(0.25,0.87,"Internal")
	

	if normalize: 
		plotting_helpers.drawNotesVertextype(vertextype, "")
	else: 
		plotting_helpers.drawNotesVertextype(vertextype, lumi)

	
	save_name = hname if save_name == "" else save_name
	if vertextype == "VSI":
		savefilename= save_name + "_compare_dataMC_VSI"
	elif vertextype == "VSI Leptons":
		savefilename= save_name + "_compare_dataMC_VSILep"
	else: 
		savefilename= save_name + "_compare_dataMC_VSILep"

	MyC01.SaveAs(outputDir +savefilename+'.pdf')




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
	atlas_style.ATLASLabel(0.25,0.87,"Internal")

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


def CorrPlot2D(file, hname, hlabel,vertextype, setxrange="",setyrange="",rebinx=1,rebiny=1, outputDir="../output/"):

	f = ROOT.TFile(file)
	if "data" in file:
		if vertextype == "VSI":
			h = f.Get(hname+"_VSI") # get histogram
			# h = f.Get(hname+"_pmu_mumu_VSI") # get histogram
		elif vertextype == "VSI Leptons": 
			h = f.Get(hname+"_VSI_Leptons") # get histogram
			# h = f.Get(hname+"_pmu_mumu_VSILep") # get histogram
	else:
		if vertextype == "VSI":
			h = f.Get(hname+"_VSI") # get histogram
		elif vertextype == "VSI Leptons": 
			h = f.Get(hname+"_VSI_Leptons") # get histogram


	#define your canvas
	MyC01= ROOT.TCanvas("MyC01","",800,500)
	MyC01.SetRightMargin(20)


	#format legend
	leg01 = ROOT.TLegend(0.50,0.8,0.89,0.89)
	leg01.SetTextSize(0.035)
	leg01.SetBorderSize(0)
	leg01.SetFillColor(kWhite)
	leg01.SetShadowColor(kWhite)
	leg01.AddEntry(h,hlabel,"")


	# Rebin 2D histogram!
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
	if "mhnl_hnlpt" in hname:
		h.GetXaxis().SetTitle("HNL mass [GeV]")
		h.GetYaxis().SetTitle("HNL p_{T} [GeV]")
	if "mhnl2D" in hname:
		h.GetXaxis().SetTitle("HNL mass (negative root) [GeV]")
		h.GetYaxis().SetTitle("HNL mass (postive root) [GeV]")



	h.GetZaxis().SetTitleOffset(-1.5)
	h.Draw("colz")


	leg01.Draw()
	atlas_style.ATLASLabel(0.25,0.87,"Internal")
	
	if "uuu" in file: 
		plotting_helpers.drawNotes("mumu","muon",vertextype) 
	elif "ueu" in file: 
		plotting_helpers.drawNotes("emu","muon",vertextype) 
	elif "uee" in file: 
		plotting_helpers.drawNotes("ee","muon",vertextype) 
	elif "eee" in file: 
		plotting_helpers.drawNotes("ee","electron",vertextype) 
	elif "eeu" in file: 
		plotting_helpers.drawNotes("emu","electron",vertextype) 
	elif "euu" in file: 
		plotting_helpers.drawNotes("mumu","electron",vertextype) 
	else: 
		plotting_helpers.drawNotesVertextype(vertextype)

	if vertextype == "VSI":
		savefilename = hname + "_2Dmass_VSI"
	elif vertextype == "VSI Leptons":
		savefilename = hname + "_2Dmass_VSILep"
	else: 
		savefilename = hname + "_2Dmass_VSILep"

	MyC01.SaveAs(outputDir +savefilename+'.pdf')
	