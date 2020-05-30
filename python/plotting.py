# plotting functions
from __future__ import division
import argparse
import os
import math
import ROOT
import glob
import uproot
import time
import json

import atlas_style
import numpy as np
import plotting_helpers
import helpers
from ROOT import kAzure, kWhite, kViolet, kBlack, TLine, TCanvas, TFile, TLegend, TPad, TGaxis
from ROOT import gPad
from pylab import *



logger = helpers.getLogger('dHNLAnalysis.plotting')


def plot_cutflow(file, vertextype, output_dir="../output/"):
	print file
	Tfile = ROOT.TFile(file)
	hcutflow = Tfile.Get('{}/CutFlow/CutFlow_{}'.format(vertextype,vertextype))

	MyC01= ROOT.TCanvas("MyC01","cutflow",1200,400)

	gPad.SetLogy()
	ymax_cutflow = hcutflow.GetMaximum()
	hcutflow.GetYaxis().SetRangeUser(0.1,ymax_cutflow*1000)
	# hcutflow.GetYaxis().SetRangeUser(0.1,ymax_cutflow)
	hcutflow.SetFillColor(kAzure-4)
	hcutflow.SetLineWidth(0)
	hcutflow.GetYaxis().SetTickLength(0.)
	hcutflow.SetMarkerSize(2.2)
	hcutflow.Draw("HIST TEXT35 SAME")


	channel = file.split("histograms_")[1].split(".")[0]

	fileInfo = helpers.FileInfo(file,channel)
	if "data" in file:
		plotting_helpers.drawNotesData("data 2018",vertextype,lumi=60,channel=channel)
	else:
		plotting_helpers.drawNotesMC(fileInfo.MC_campaign,vertextype,channel,fileInfo.mass_str,fileInfo.ctau_str)


	savefilename= "CutFlow_" + channel

	output_dir = os.path.join(os.path.abspath(output_dir), '{}/'.format(vertextype))
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	MyC01.SaveAs(output_dir +savefilename+'.pdf')


def compare(hist_channels, variable="", setrange=None, scaleymax=1.2, nRebin=1, setlogy=False, output_dir="../output",
			save_name="", vertical_lines=[], labels=[], normalize=True, lumi=1, **kwargs):
	histograms = []
	filenames = []
	labels = []
	channels = []  # TODO: Why do we need channels? It seems to be never used.
	tfiles = []  # root is stupid and will close the file if you're not careful
	# for key, val in hist_channels.items():
	for nhist in range(len(hist_channels)):
		filename = hist_channels[nhist][0]
		label = hist_channels[nhist][1]
		vtx_alg = hist_channels[nhist][2]
		selection = hist_channels[nhist][3]

		if 'data' in filename:
			channel = filename.split("data_")[1].split(".")[0]
		else:
			channel = filename.split("mm_")[1].split(".")[0]
		tfiles.append(ROOT.TFile(filename))  # get file

		if 'use_ntuple' in kwargs and kwargs['use_ntuple']:
			print('Using ntuple')
			if setrange is None or 'ntup_nbins' not in kwargs:
				raise ValueError('to use the ntuple, you must supply the range (e.g. setrange=(0,1000) '
								 'and the number of bins (e.g. ntup_nbins=40) in the arguments)')
			tmp_hist_name = "{}_{}_{}".format(vtx_alg, selection, variable)
			ntup_hist = ROOT.TH1D(tmp_hist_name, tmp_hist_name, kwargs['ntup_nbins'], setrange[0], setrange[1])  # create empty histogram
			ttree = tfiles[nhist].Get("{}/ntuples/{}".format(vtx_alg, selection))  # get TTree
			if not ttree:
				raise KeyError('Cannot find {}/ntuples/{}'.format(vtx_alg, selection))
			ttree.Draw(tmp_hist_name)  # fill histogram with data from ttree
			histograms.append(ntup_hist)
		else:
			hist_path = "{}/{}/{}".format(vtx_alg, selection, variable)
			histogram = tfiles[nhist].Get(hist_path)
			if not histogram:  # no histogram object. don't even try
				print("cannot find {}. Exiting".format(variable))
				return
			histograms.append(histogram)  # get variable with suffix
		channels.append(channel)
		filenames.append(filename)
		labels.append(label)

	n_h = len(histograms)
	h_idx = range(len(histograms))

	# define your canvas
	c = ROOT.TCanvas("canvas", "", 1200, 800)

	# format legend
	leg01 = ROOT.TLegend(0.57, 0.71, 0.92, 0.92)
	leg01.SetTextSize(0.035)
	leg01.SetBorderSize(0)
	leg01.SetFillColor(kWhite)
	leg01.SetShadowColor(kWhite)

	# rebin histograms
	if nRebin != 1:
		for i in h_idx:
			histograms[i].Rebin(nRebin)

	# find the common min and max for x axis	
	if setrange is None: # if the range is not defined by the user
		bin_xmin, bin_xmax = (np.inf, -np.inf)
		for i in h_idx:
			if histograms[i].FindFirstBinAbove(0, 1) > bin_xmin:
				bin_xmin = histograms[i].FindFirstBinAbove(0, 1)
			if histograms[i].FindLastBinAbove(0, 1) < bin_xmax:
				bin_xmax = histograms[i].FindLastBinAbove(0, 1)
		x_min = bin_xmin-1
		x_max = bin_xmax+1
	else: # if the range is designed by the user
		x_min, x_max = setrange
		bin_xmax = histograms[0].GetXaxis().FindBin(x_max)
		bin_xmin = histograms[0].GetXaxis().FindBin(x_min)

	# normalize the histograms
	if normalize:
		norm = 1
		if 'norm' in kwargs:
			norm = kwargs['norm']

		for i in h_idx:
			if histograms[i].Integral() != 0:
				scale_mc = norm / (histograms[i].Integral(bin_xmin, bin_xmax))
			else:
				scale_mc = norm
			histograms[i].Scale(scale_mc)

	# default scale histograms to a given luminosity
	else:
		for i in h_idx:
			if 'data' not in labels[i]:  # don't scale data histograms!
				histograms[i].Scale(lumi)

	# calculate yields
	mc_yield = {}
	for i in h_idx:
		mc_yield[i] = round(histograms[i].Integral(histograms[i].FindFirstBinAbove(0,1), histograms[i].FindLastBinAbove(0,1)),2)
		# if data in list, add it to the legend first
		if 'data' in labels[i]:
			if normalize:
				leg01.AddEntry(histograms[i],"\\bf{%s} )"%(labels[i]),"lp")
			else:
				leg01.AddEntry(histograms[i],"\\bf{%s}, \\bf{%s)}"%(labels[i],mc_yield[i]),"lp")
	#add non-data histograms to legend
	leg01.AddEntry(histograms[0],"","")
	for i in h_idx:
		if 'data' not in labels[i]: 
			if normalize: 
				leg01.AddEntry(histograms[i],"\\bf{%s} )"%(labels[i]),"lp")
			else:
				leg01.AddEntry(histograms[i],"\\bf{%s}, \\bf{%s)}"%(labels[i],mc_yield[i]),"lp")
		

	# set the common x limits for all histograms
	for i in h_idx:
		histograms[i].GetXaxis().SetRangeUser(x_min, x_max)

	# find the common min and max for y axis
	if 'y_max' in kwargs:
		y_max = kwargs['y_max']
	else:
		y_max = -np.inf
		for i in h_idx:
			if histograms[i].GetMaximum() > y_max:
				y_max = histograms[i].GetMaximum()

	shapelist = [22, 21, 33, 29, 30, 31, 32, 34, 35]
	for i in h_idx:
		if 'bin_labels' in kwargs:
			bin_labels = kwargs['bin_labels']
			for j, label in enumerate(bin_labels):
				histograms[i].GetXaxis().SetBinLabel(j+1, label)
		histograms[i].SetMarkerSize(1.5)
		if 'data' in labels[i]: 
			histograms[i].SetLineColor(kBlack)
			histograms[i].SetMarkerColor(kBlack)
			histograms[i].SetMarkerStyle(20)
			histograms[i].SetLineWidth(2)
		else: 
			histograms[i].SetLineWidth(2)
			histograms[i].SetLineColor(plotting_helpers.histColours(i))
			histograms[i].SetMarkerColor(plotting_helpers.histColours(i))
			histograms[i].SetMarkerStyle(shapelist[i])
		histograms[i].GetXaxis().SetTitle(plotting_helpers.get_x_label(variable))
		if not variable: histograms[i].GetXaxis().SetTitle(save_name)
		histograms[i].GetYaxis().SetTitle("entries")
		histograms[i].GetYaxis().SetRangeUser(0.00001 if setlogy else 0, y_max*10**scaleymax if setlogy else y_max*scaleymax)
		histograms[i].Draw("E0 HIST SAME")
	
	if setlogy:
		gPad.SetLogy()

	# draw vertical lines
	lines = []
	for i, x in enumerate(vertical_lines):
		lines.append(TLine(x, 0, x, y_max))
		lines[i].SetLineStyle(3)
		lines[i].SetLineWidth(3)
		lines[i].Draw("SAME")

	if 'vertical_legend' in kwargs:
		leg01.AddEntry(lines[0], "\\bf{%s}"%kwargs['vertical_legend'], "l")

	leg01.Draw()

	l = ROOT.TLatex()
	l.SetNDC()
	l.SetTextFont(43)
	l.SetTextColor(1)
	l.SetTextSize(20)
 	l.DrawLatex(0.65,0.855,"(     m_{HNL},      c\\tau,      chan.,      yield)")

	notes_x = 0.25
	notes_y = 0.87

	if normalize: 
		plotting_helpers.drawNotes(vtx_alg, "")
	else: 
		plotting_helpers.drawNotes(vtx_alg, lumi)

	if 'notes' in kwargs:
		for note in kwargs['notes']:
			notes_y -= 0.07
			plotting_helpers.drawNote(note, size=40, ax=notes_x, ay=notes_y)

	# plotting_helpers.getNote(35).DrawLatex(notes_x, notes_y-.05, vertextype)

	save_file_name = variable if save_name == "" else save_name
	output_dir = os.path.join(os.path.abspath(output_dir), 'plots/{}/{}/'.format(vtx_alg, selection))

	if not os.path.exists(output_dir): os.mkdir(output_dir)
	c.SaveAs(output_dir + save_file_name + '.pdf')
	# c.SaveAs(output_dir + save_file_name + '.eps')


def compare2_wRatio(histos1, histos2, h1name, h1label, h2name, h2label, xlabel,savefilename,outputDir):
	############################################################################
	nRebin = 10 # set to 1 if you dont want to rebin.
	scaleymax = 1.6 # use this to scale height of y axis for asthetics
	############################################################################

	# get 2 histograms from the file
	file1 = ROOT.TFile(histos1)
	file2 = ROOT.TFile(histos2)
	h1 = file1.Get(h1name)
	h2 = file2.Get(h2name)

	#define your canvas
	MyC01 = ROOT.TCanvas("MyC01", "canvas", 800, 800)
	MyC01.Divide(1,1)
	MyC01.cd(1)
	pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
	pad1.SetBottomMargin(0)
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

	shapelist = [22,21,33,29]
	h1.SetLineColor(plotting_helpers.histColours(0))
	h1.SetMarkerSize(0.9)
	h1.SetLineColor(plotting_helpers.histColours(0))
	h1.SetMarkerColor(plotting_helpers.histColours(0))
	h1.SetMarkerStyle(shapelist[0])
	h1.GetXaxis().SetTitle(xlabel)
	h1.GetYaxis().SetTitle("entries")
	h1.GetYaxis().SetRangeUser(0,y_max*scaleymax)
	h1.GetXaxis().SetRange(bin_xmin-1,bin_xmax+1)
	h1.Draw("HIST")



	h2.SetLineColor(plotting_helpers.histColours(1))
	h2.SetMarkerSize(0.9)
	h2.SetLineColor(plotting_helpers.histColours(1))
	h2.SetMarkerColor(plotting_helpers.histColours(1))
	h2.SetMarkerStyle(shapelist[1])
	h2.SetLineStyle(2)
	h2.Draw("HIST SAME")

	if "DV_r" in h1name or "DV_r" in h2name:
		if  "redmass" in h1name or "redmass" in h2name:
			pass
		else:
			matlayers = [33.25,50.5,88.5,122.5,299]
			nmatlayers = len(matlayers)
			matlay={}
			for i in range(nmatlayers):
				matlay[i]=TLine(matlayers[i],0,matlayers[i],y_max*1.05)
				matlay[i].SetLineStyle(3)
				matlay[i].Draw("SAME")
			leg01.AddEntry(matlay[0],"material layers","l")


	leg01.Draw()
	plotting_helpers.drawNotesVertextype("VSI")

	MyC01.cd(1)
	pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
	pad2.SetTopMargin(0)
	pad2.SetBottomMargin(0.2)
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
	hratio.GetYaxis().SetRangeUser(0.8,1.23)
	# line1=TLine(0,1,500,1)
	# line1.SetLineStyle(1)

	hratio.SetMarkerColor(plotting_helpers.histColours("ratio"))

	hratio.GetXaxis().SetRange(bin_xmin-1,bin_xmax+1)


	hratio.Draw("ep")
	# line1.Draw("SAME")

	if "DV_r" in h1name or "DV_r" in h2name:
		if  "redmass" in h1name or "redmass" in h2name:
			pass
		else:
			matlayers_ratio = [33.25,50.5,88.5,122.5,299]
			nmatlayers_ratio = len(matlayers_ratio)
			matlay_ratio={}
			for i in range(nmatlayers):
				matlay_ratio[i]=TLine(matlayers_ratio[i],0.8,matlayers_ratio[i],y_max)
				matlay_ratio[i].SetLineStyle(3)
				matlay_ratio[i].Draw("SAME")



	MyC01.SaveAs(outputDir +savefilename+'.pdf')


def CorrPlot2D(file, hname, hlabel,vertextype, setxrange="",setyrange="",rebinx=1,rebiny=1, lumi=1, outputDir="../output/"):

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
	leg01 = ROOT.TLegend(0.65,0.82,0.82,0.86)
	leg01.SetTextSize(0.035)
	leg01.SetBorderSize(0)
	leg01.SetFillColor(kWhite)
	leg01.SetShadowColor(kWhite)
	leg01.AddEntry(h,"\\bf{%s  ) }"%hlabel,"")

	if "data" in file:
		pass
	else:
		h.Scale(lumi)

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
	if "DVmass_hnlpt" in hname:
		h.GetXaxis().SetTitle("DV mass [GeV]")
		h.GetYaxis().SetTitle("HNL p_{T} [GeV]")
	if "DVmass_mtrans" in hname:
		h.GetXaxis().SetTitle("DV mass [GeV]")
		h.GetYaxis().SetTitle("Visible m_{T} [GeV]")
	if "mvis_mhnl" in hname:
		h.GetXaxis().SetTitle("Visible mass (m_{lll}) [GeV]")
		h.GetYaxis().SetTitle("HNL mass [GeV]")
	if "mvis_mtrans" in hname:
		h.GetXaxis().SetTitle("Visible Mass (m_{lll}) [GeV]")
		h.GetYaxis().SetTitle("Visible m_{T} [GeV]")
	if "mvis_hnlpt" in hname:
		h.GetXaxis().SetTitle("Visible Mass (m_{lll}) [GeV]")
		h.GetYaxis().SetTitle("HNL p_{T} [GeV]")
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

	channel = file.split("histograms_")[1].split(".")[0]

	if "data" in file:
		plotting_helpers.drawNotesVertextype(vertextype,lumi=60)
	else:
		plotting_helpers.drawNotesVertextype(vertextype,lumi=lumi)

	l = ROOT.TLatex()
	l.SetNDC()
	l.SetTextFont(43)
	l.SetTextColor(1)
	l.SetTextSize(14)
	l.DrawLatex(0.69,0.88,"(  m_{HNL},    c\\tau,    chan.  )")




	if vertextype == "VSI":
		savefilename = hname + "_2Dmass_VSI"
	elif vertextype == "VSI Leptons":
		savefilename = hname + "_2Dmass_VSILep"
	else:
		savefilename = hname + "_2Dmass_VSILep"

	MyC01.SaveAs(outputDir +savefilename+'.pdf')


