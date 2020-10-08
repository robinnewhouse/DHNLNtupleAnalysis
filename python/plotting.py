# plotting functions
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
shapelist = [22, 21, 33, 29, 30, 31, 32, 34, 35]

def plot_cutflow(file, selection,vertextype, output_dir="../output/"):
	Tfile = ROOT.TFile(file)
	
	if selection == "mixed": 
		hcutflow1 = Tfile.Get('{}/CutFlow/CutFlow_{}'.format(vertextype,"LNC"))
		hcutflow2 = Tfile.Get('{}/CutFlow/CutFlow_{}'.format(vertextype,"LNV"))
		if not hcutflow1 or not hcutflow2:  # no histogram object. don't even try
			print("Cannot find LNC &LNV cutflows for {} vertex configuration. Exiting".format(vertextype))
			return
		hcutflow1.Scale(0.5)
		hcutflow2.Scale(0.5)
		hcutflow = hcutflow1 + hcutflow2
	else:
		if selection == "all":
			hcutflow = Tfile.Get('{}/CutFlow/CutFlow'.format(vertextype))
		else: 
			hcutflow = Tfile.Get('{}/CutFlow/CutFlow_{}'.format(vertextype,selection))
		if not hcutflow:  # no histogram object. don't even try
			print("Cannot find cutflow for {} vertex configuration. Exiting".format(vertextype))
			return

	MyC01= ROOT.TCanvas("MyC01","cutflow",1200,400)
	# MyC01= ROOT.TCanvas("MyC01","cutflow",900,400)
	gPad.SetLogy()
	ymax_cutflow = hcutflow.GetMaximum()
	hcutflow.GetYaxis().SetRangeUser(0.1,ymax_cutflow*10000)
	hcutflow.SetFillColor(kAzure-4)
	# color = ROOT.TColor()
	# ncolor = color.GetColor("#FF7E00") # orange for VSI vs VSI Lep
	# hcutflow.SetFillColor(ncolor)
	hcutflow.SetLineWidth(0)
	hcutflow.GetYaxis().SetTickLength(0.)
	hcutflow.SetMarkerSize(2.2)
	# hcutflow.Scale(140)
	hcutflow.GetXaxis().SetLabelSize(0.05)
	# hcutflow.GetXaxis().SetRange(1,10)
	# hcutflow.LabelsOption("v")

	hcutflow.Draw("HIST SAME TEXT35")
	
	# Print rounded numbers on histogram
	# text = ROOT.TLatex()
	# text.SetTextAlign(21)
	# text.SetTextColor(hcutflow.GetMarkerColor());
	# text.SetTextSize(0.03*hcutflow.GetMarkerSize());
	# # text.SetTextAngle(35)
	# for j in xrange(1, hcutflow.GetNbinsX()+1):
	# 	x  = hcutflow.GetXaxis().GetBinCenter(j)
	# 	y  = hcutflow.GetBinContent(j)
	# 	if y != 0.0: 
	# 		text.DrawLatex(x, y+0.5, '%.1f' % y)

	channel = file.split(".root")[0].split("histograms_")[1]

	fileInfo = helpers.FileInfo(file,channel)
	if "data" in file:
		plotting_helpers.drawNotesData("data 2018",vertextype,lumi="",channel=channel)
	else:
		plotting_helpers.drawNotesMC(fileInfo.MC_campaign,vertextype,channel,fileInfo.mass_str,fileInfo.ctau_str)

	savefilename= selection+"_CutFlow_" + channel

	output_dir = os.path.join(os.path.abspath(output_dir),"Cutflows/", '{}/'.format(vertextype) )
	if not os.path.exists(output_dir): os.mkdir(output_dir)

	MyC01.SaveAs(output_dir +savefilename+'.pdf')


def compare(hist_channels, 
			variable, 
			setrange=None, 
			scaleymax=1.2, 
			nRebin=1, 
			setlogy=False, 
			output_dir="../output/",
			save_name="", 
			vertical_lines=[], 
			labels=[], 
			normalize=False, 
			drawRatio=False, 
			customVariable = False, 
			scalelumi=1.0,
			datalumi = 140.0, 
			do_cut_significance=False,
			debug = False, # enable debug printouts 
			**kwargs):
	
	draw_channel_info = not normalize
	histograms = []
	filenames = []
	labels = []
	tfiles = []  # root is stupid and will close the file if you're not careful

	ROOT.gStyle.SetOptTitle(0)
	for nhist in range(len(hist_channels)):
		filename = hist_channels[nhist][0]
		label = hist_channels[nhist][1]
		vtx_alg = hist_channels[nhist][2]
		selection = hist_channels[nhist][3]
		if len(hist_channels[nhist]) > 4: 
			MCtype = hist_channels[nhist][4]
		else: 
			MCtype = None
		tfiles.append(ROOT.TFile(filename))  # get file
	
		if  MCtype != None: 
			if "truth" in selection:  #truth has this extra folder to get to all...small hack -DT
				hist_path = "{}/{}/all/{}/{}".format(vtx_alg, selection,MCtype, variable)
			else: 
				hist_path = "{}/{}/{}/{}".format(vtx_alg, selection,MCtype, variable)
		else: 
			hist_path = "{}/{}/{}".format(vtx_alg, selection, variable)
	
		###################################################################################################
		# Plot with micro-ntuples 
		###################################################################################################
		if 'use_ntuple' in kwargs and kwargs['use_ntuple']:
			print('Using ntuple')
			if setrange is None or 'ntup_nbins' not in kwargs:
				raise ValueError('to use the ntuple, you must supply the range (e.g. setrange=(0,1000) '
								 'and the number of bins (e.g. ntup_nbins=40) in the arguments)')
			if MCtype != None: tmp_hist_name = "{}_{}_{}_{}".format(vtx_alg, MCtype,selection, variable)
			else: tmp_hist_name = "{}_{}_{}".format(vtx_alg, selection, variable)
			ntup_hist = ROOT.TH1D(tmp_hist_name, tmp_hist_name, kwargs['ntup_nbins'], setrange[0], setrange[1])  # create empty histogram
			ROOT.gStyle.SetOptTitle(0)
			if customVariable: 
				if not 'variable2' in kwargs or not 'customVariableAction'in kwargs:
					raise ValueError('To use a custom variable you must give variable, variable2 and a customVariableAction. Supported actions include "divide" and "multiply".')
				if MCtype != None: chain = ROOT.TChain('{}_ntuples_{}_{}'.format(vtx_alg, MCtype,selection)) 
				else: chain = ROOT.TChain('{}_ntuples_{}'.format(vtx_alg, selection)) 
				if not chain:
					raise KeyError('Cannot find {}_ntuples_{} micro-ntuple in file {}'.format(vtx_alg, selection, tfiles[nhist]))
				chain.Add(filename) 

				for event in chain: 
					var1 = event.GetLeaf(variable).GetValue() 
					var2 = event.GetLeaf(kwargs['variable2']).GetValue() 
					if kwargs['customVariableAction'] == "divide":
						customVar = var1/var2
						ntup_hist.GetXaxis().SetTitle("{}/{}".format(variable, kwargs['variable2'] ))
					if kwargs['customVariableAction'] == "multiply":
						customVar = var1*var2
						ntup_hist.GetXaxis().SetTitle("{}*{}}".format(variable, kwargs['variable2'] ))
					ntup_hist.Fill(customVar)
			else: 
				if MCtype != None: ttree = tfiles[nhist].Get('{}_ntuples_{}_{}'.format(vtx_alg, MCtype, selection))  # get TTree
				else: ttree = tfiles[nhist].Get('{}_ntuples_{}'.format(vtx_alg, selection))  # get TTree
				if not ttree:
					raise KeyError('Cannot find {}_ntuples_{} micro-ntuple in file {}'.format(vtx_alg, selection, tfiles[nhist]))
			
				ttree.Draw(variable+'>>'+tmp_hist_name, 'DV_weight')  # fill histogram with data from ttree. weighted with DV_weight
				# ttree.Draw(variable+'>>'+tmp_hist_name, 'DV_weight*(DV_mass > 2)')  # fill histogram with data from ttree. weighted with DV_weight
			ntup_hist.SetTitle("")
			histograms.append(ntup_hist)

		###################################################################################################
		# Plot with pre-binned histograms 
		###################################################################################################	
		else:
			histogram = tfiles[nhist].Get(hist_path)
			if not histogram:  # no histogram object. don't even try
				print('cannot find {}. Exiting'.format(variable))
				return
			histogram.SetTitle("")
			histograms.append(histogram)  # get variable with suffix

		filenames.append(filename)
		labels.append(label)

	if do_cut_significance: 
		makeAsimov(histograms[0],histograms[1],variable,selection, vtx_alg, scalelumi,datalumi, output_dir)
	
	n_h = len(histograms)
	h_idx = range(len(histograms))
	# define your canvas
	c = ROOT.TCanvas("canvas", "", 1200, 800)
	if drawRatio:
		c.Divide(1,1)
		c.cd(1)
		pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
		pad1.SetBottomMargin(0)
		pad1.Draw()
		pad1.cd()

	# format legend
	leg01 = ROOT.TLegend(0.57, 0.71, 0.92, 0.92)
	leg01.SetTextSize(0.025)
	leg01.SetBorderSize(0)
	leg01.SetFillColor(kWhite)
	leg01.SetShadowColor(kWhite)

	#rebin histograms
	if nRebin != 1:
		for i in h_idx:
			histograms[i].Rebin(nRebin)

			
	# find the common min and max for x axis	
	if setrange is None:
		bin_xmin, bin_xmax = (-np.inf, np.inf)
		for i in h_idx:
			if histograms[i].FindFirstBinAbove(0, 1) > bin_xmin:
				bin_xmin = histograms[i].FindFirstBinAbove(0, 1)
			if histograms[i].FindLastBinAbove(0, 1) < bin_xmax:
				bin_xmax = histograms[i].FindLastBinAbove(0, 1)
		x_min = bin_xmin-1
		x_max = bin_xmax+1
	else:
		x_min, x_max = setrange
		bin_xmax = histograms[0].GetXaxis().FindBin(x_max)
		bin_xmin = histograms[0].GetXaxis().FindBin(x_min)


	# normalize the histograms
	if normalize:
		if 'norm' in kwargs:
			norm = kwargs['norm']
		else: 
			norm = 1
		for i in h_idx:
			if (histograms[i].Integral() != 0):
				scale_mc = norm/(histograms[i].Integral(bin_xmin, bin_xmax))
			else:
				scale_mc = norm
			histograms[i].Scale(scale_mc)

	#default scale historams to a given luminosity 
	else: 
		for i in h_idx:
			if 'SS bkg' in labels[i]:
				#scale data histograms to scalelumi
				scale_data = scalelumi / datalumi
				histograms[i].Scale(scale_data)
				
			else:
				#scale MC histograms to scalelumi
				histograms[i].Scale(scalelumi)

	# calculate yields
	Yield = {}
	for i in h_idx:
		Yield[i] = round(histograms[i].Integral(histograms[i].FindFirstBinAbove(0,1), -1),2)
		# if data in list, add it to the legend first
		if 'SS bkg' in labels[i]: 
			if normalize:
				leg01.AddEntry(histograms[i],"\\bf{%s)}"%(labels[i]),"f")
			else: 
				leg01.AddEntry(histograms[i],"\\bf{%s}, \\bf{%s)}"%(labels[i],Yield[i]),"f")
	#add non-data histograms to legend 
	if draw_channel_info: leg01.AddEntry(histograms[0],"","")
	for i in h_idx:
		if 'SS bkg' not in labels[i]: 
			if normalize: 
				leg01.AddEntry(histograms[i],"\\bf{%s)}"%(labels[i]),"l")
			else:
				leg01.AddEntry(histograms[i],"\\bf{%s}, \\bf{%s)}"%(labels[i],Yield[i]),"l")
		

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

	for i in h_idx:
		if 'bin_labels' in kwargs:
			bin_labels = kwargs['bin_labels']
			for j, label in enumerate(bin_labels):
				histograms[i].GetXaxis().SetBinLabel(j+1, label)
		histograms[i].SetMarkerSize(1.5)
		if 'SS bkg' in labels[i]: 
			histograms[i].SetLineColor(plotting_helpers.bkgColours(i))
			histograms[i].SetFillColor(plotting_helpers.bkgColours(i))
			histograms[i].SetMarkerColor(plotting_helpers.bkgColours(i))
			histograms[i].SetLineColor(ROOT.kBlack)
			histograms[i].SetLineWidth(1)
			# histograms[i].SetMarkerStyle(20)
		else: 
			histograms[i].SetLineWidth(2)
			histograms[i].SetLineColor(plotting_helpers.histColours(i))
			# histograms[i].SetLineColor(ROOT.kBlack)
			histograms[i].SetMarkerColor(plotting_helpers.histColours(i))
			histograms[i].SetMarkerStyle(shapelist[i])
		if customVariable == False:
			histograms[i].GetXaxis().SetTitle(plotting_helpers.get_x_label(variable))
		if not variable: histograms[i].GetXaxis().SetTitle(save_name)
		histograms[i].GetYaxis().SetTitle("entries")
		histograms[i].GetYaxis().SetRangeUser(0.00001 if setlogy else 0, y_max*10**scaleymax if setlogy else y_max*scaleymax)
		# histograms[i].GetYaxis().SetRangeUser(0,400)
		histograms[i].Draw("HIST SAME")
		# histograms[i].Draw("E0 HIST SAME") # plot with errors

	if debug: 
		print("Making histogram for {}: ".format(variable))
		for i in h_idx:
			print("	" + labels[i])
			print( "		num. of events: " + str( histograms[i].GetEntries() ) ) 
			print( "		mean: " +  str( histograms[i].GetMean()) ) 
			print( "		std dev: " + str( histograms[i].GetStdDev()) )



		# get rounded numbers to draw on histogram
		# text = ROOT.TLatex()
		# text.SetTextAlign(21)
		# text.SetTextColor(histograms[i].GetMarkerColor());
		# text.SetTextSize(0.02*histograms[i].GetMarkerSize());

		# for j in xrange(1, histograms[i].GetNbinsX()+1):
		# 	x  = histograms[i].GetXaxis().GetBinCenter(j)
		# 	y  = histograms[i].GetBinContent(j)
		# 	if y != 0.0: 
		# 		text.DrawLatex(x, y+0.3, '%.1f' % y)
	
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
	l.SetTextSize(15)

	if draw_channel_info:
		l.DrawLatex(0.7,0.82,"(m_{HNL} ,  c\\tau ,  chan.  ,   yield)")

	notes_x = 0.25
	notes_y = 0.87

	if "truth" in selection: 
		truth_plot = True
	else: 
		truth_plot = False

	if normalize: 
		plotting_helpers.drawNotes(vtx_alg, "", truth_plot)
	else: 
		plotting_helpers.drawNotes(vtx_alg, scalelumi, truth_plot)

	if 'notes' in kwargs:
		for note in kwargs['notes']:
			notes_y -= 0.07
			plotting_helpers.drawNote(note, size=40, ax=notes_x, ay=notes_y)


	if drawRatio:
		c.cd(1)
		pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
		pad2.SetTopMargin(0)
		pad2.SetBottomMargin(0.2)
		pad2.Draw()
		pad2.cd()

		hist_numerator =  histograms[1]
		hist_denomenator =  histograms[0]
		hratio = hist_numerator.Clone("hratio")
		hratio.SetLineColor(kBlack)
		hratio.SetMinimum(0.8)
		hratio.SetMaximum(1.35)
		# hratio.Sumw2()
		hratio.SetStats(0)
		hratio.Divide(hist_denomenator)
		hratio.SetMarkerStyle(21)
		hratio.SetLineStyle(1)
		hratio.SetMarkerStyle(20)
		hratio.SetMarkerSize(1)
		if 'ratioLabel' in kwargs:
			hratio.GetYaxis().SetTitle(kwargs['ratioLabel'][0])
		else: 
			hratio.GetYaxis().SetTitle("ratio")
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
		hratio.GetYaxis().SetRangeUser(-1 ,3)
		line1=TLine(x_min,1,x_max,1)
		line1.SetLineStyle(1)


		hratio.SetMarkerColor(plotting_helpers.histColours("ratio"))
		hratio.Draw("ep")
		line1.Draw("SAME")

	if customVariable == True: 
		if save_name == "": 
			save_file_name = "{}_{}_{}_{}".format(selection, variable ,kwargs["variable2"],kwargs["customVariableAction"])
		else: 
			save_file_name = "{}_{}".format(selection, save_name)
	else:
		save_file_name = "{}_{}".format(selection, variable if save_name == "" else save_name)

	# Clean output directory
	# if vtx_alg == "VSI": 
	if "VSI_LRT" in vtx_alg: 
		output_dir = os.path.join(os.path.abspath(output_dir), 'plots/VSI/')
	# elif vtx_alg == "VSI_Leptons": 
	elif "VSI_Leptons_LRT" in vtx_alg: 
		output_dir = os.path.join(os.path.abspath(output_dir), 'plots/VSI_Leptons/')
	else:
		output_dir = os.path.join(os.path.abspath(output_dir), 'plots/')
	
	if not os.path.exists(output_dir): os.mkdir(output_dir)
	if not os.path.exists(output_dir+"eps_files/"): os.mkdir(output_dir+"eps_files/")
	
	c.SaveAs(output_dir + save_file_name + '.pdf')
	c.SaveAs(output_dir +"eps_files/" +save_file_name + '.eps')

	if do_cut_significance: 
		if 'ncuts' in kwargs:
			ncuts = kwargs['ncuts']
			cut_significance(variable,vtx_alg,scalelumi,histograms,filenames,labels,x_min,x_max, ncuts = ncuts)
		else:
			cut_significance(variable,vtx_alg,scalelumi,histograms,filenames,labels,x_min,x_max,output_dir,save_file_name)


def cut_significance(variable,vtx_alg,lumi,histograms,filenames,labels, x_min,x_max, output_dir, savefilename, ncuts=0):
	h_idx = range(len(histograms))
	
	#get the number of cuts you are going to make 
	if ncuts == 0: 
		nbins = histograms[0].GetNbinsX()
		step_size = (x_max- x_min) / nbins 
	else: 
		nbins = ncuts
		step_size = (x_max- x_min) / nbins 

	b = []
	s = {}
	s_label = {}
	cut_list = []
	# loop over all the histograms 
	for i in h_idx:
		#get the channel from the filename
		channel = filenames[i].split(".root")[0].split("histograms_")[1]
		s_label[channel] = labels[i]
		#loop over the cut you made
		for j in range(0,nbins+1):
			cut = x_min + j*step_size
			cut_list.append(cut)
			cut_yield = histograms[i].Integral(histograms[i].FindFixBin(cut),-1)
			full_yield = histograms[i].Integral(histograms[i].FindFirstBinAbove(0,1), -1)

			# check if you have a background estimate (need a better way to do this) -DT
			if 'data' in channel:
				b.append(cut_yield)
			else: 
				if j == 0:
					s[channel] = [cut_yield]
				else:
					s[channel].append(cut_yield)
	if len(b) == 0:
		print "Please provide a background estimate, if you want to calculate the significance. Exiting."
		return

	g_sig = {}
	c2 = {}
	leg = {}
	for key in s:
		g_sig[key] = ROOT.TGraph()
		n_cuts = len(s[key])
		yvals = []

		if len(b) != n_cuts: 
			print "Different number of signal and bkg cuts were applied. Cannot calculate significance. Exiting."
			return

		for npoint in range(n_cuts):
			if b[npoint] != 0: 
				significance_bkgONLY = s[key][npoint]/np.sqrt(b[npoint])
				significance_sPLUSb = s[key][npoint]/np.sqrt(s[key][npoint] + b[npoint])
				g_sig[key].SetPoint(npoint, cut_list[npoint], significance_bkgONLY)
				yvals.append(significance_bkgONLY)


		index = list(s).index(key)
		g_sig[key].SetLineColor(plotting_helpers.histColours(index))
		g_sig[key].SetMarkerColor(plotting_helpers.histColours(index))
		g_sig[key].GetYaxis().SetTitle("s/\surdb")
		g_sig[key].GetXaxis().SetTitle(plotting_helpers.get_x_label(variable))
		g_sig[key].SetLineWidth(2)
		y_max = max(yvals)
		g_sig[key].GetYaxis().SetRangeUser(0.00001, y_max*1.5)
		c2[key] = ROOT.TCanvas("sig_canvas_{}".format(index), "", 1200, 800)
		# format legend
		leg[key] = ROOT.TLegend(0.65, 0.71, 0.92, 0.92)
		leg[key].SetTextSize(0.035)
		leg[key].SetBorderSize(0)
		leg[key].SetFillColor(kWhite)
		leg[key].SetShadowColor(kWhite)

		leg[key].AddEntry(g_sig[key],"\\bf{%s}  )"%(s_label[key]),"l")
		g_sig[key].Draw("AL")
		leg[key].Draw()
		
		plotting_helpers.drawNotes(vtx_alg, lumi)

		l = ROOT.TLatex()
		l.SetNDC()
		l.SetTextFont(43)
		l.SetTextColor(1)
		l.SetTextSize(20)
		l.DrawLatex(0.715,0.855,"(     m_{HNL},      c\\tau,      chan.)")

		output_path = output_dir +"Cut_significance/"
		if not os.path.exists(output_path): os.mkdir(output_path)	
		c2[key].SaveAs(output_path + savefilename + "_{}".format(key) + '.pdf')


def makeAsimov(hbkg,hsignal,variable,selection, vtx_alg, scalelumi, datalumi, output_dir): 
	#scale MC histogram to luminosity
	c_asimov = ROOT.TCanvas("canvas_asimov_{}_{}".format(variable,vtx_alg), "", 1200, 800)
	c_asimov.cd()

	
	outFile = ROOT.TFile(output_dir + "Asimov_dataset.root","update")
	
	asimov_name = "asimov_{}".format(variable)
	signal_name = "signal_{}".format(variable)
	bkg_name =  "bkg_{}".format(variable)

	leg01 = ROOT.TLegend(0.57, 0.71, 0.92, 0.92)
	leg01.SetTextSize(0.035)
	leg01.SetBorderSize(0)
	leg01.SetFillColor(kWhite)
	leg01.SetShadowColor(kWhite)
	hsignal_clone = hsignal.Clone()
	hsignal_clone.Scale(scalelumi)


	hbkg_clone = hbkg.Clone()
	hbkg_clone.Scale(scalelumi/datalumi)

	hasimov = hsignal_clone.Clone(asimov_name)

	hasimov.Add(hsignal_clone,hbkg_clone,1,1)


	nbinsx = hasimov.GetXaxis().GetNbins()



	for binx in range(1, nbinsx):
		bin_content = hasimov.GetBinContent(binx)
		bin_error = np.sqrt(bin_content)
		hasimov.SetBinError(binx,bin_error)


	# Move ROOT toe base directory
	outFile.cd()
	# Make a subdirectory for vertex type. May be other channels in the future.
	if not outFile.FindObject(vtx_alg):
		outFile.mkdir(vtx_alg, "Analysis Channel " + vtx_alg)
	# Move ROOT to the channel subdirectory
	outFile.cd(vtx_alg)


	# Store saved histograms to file
	sel_dir = vtx_alg + '/' + selection
	if not outFile.GetDirectory(sel_dir):  # make TDirectory if necessary
		outFile.mkdir(sel_dir, "Analysis Selection " + selection)
	outFile.cd(sel_dir)  # change to TDirectory

	hsignal.Write(signal_name)
	hbkg.Write(bkg_name)
	hasimov.Write(asimov_name) 
	

	hsignal_clone.SetFillColor(632)
	hsignal_clone.SetLineColor(632)
	color = ROOT.TColor()
	hbkg_clone.SetFillColor(color.GetColor("#c2a5cf"))
	hbkg_clone.SetLineColor(color.GetColor("#c2a5cf"))
	hasimov.SetFillColor(0)
	hasimov.SetLineColor(ROOT.kBlack)
	hasimov.SetMarkerColor(ROOT.kBlack)
	hasimov.SetMarkerStyle(20)
	
	leg01.AddEntry(hsignal_clone,"HNL signal","f")
	leg01.AddEntry(hbkg_clone,"SS bkg","f")
	leg01.AddEntry(hasimov,"Asimov data set")
	a = ROOT.THStack("a_{}".format(variable),"Stacked 2D histograms")
	a.Add(hbkg_clone)
	a.Add(hsignal_clone)
	
	

	# hsignal_clone.Draw("HIST")
	# hbkg_clone.Draw("HIST SAME")
	a.Draw("HIST")
	a.GetXaxis().SetRangeUser(0,30)
	a.GetXaxis().SetTitle("m_{HNL} [GeV]")
	a.GetYaxis().SetRangeUser(0,a.GetMaximum()*2.5)

	hasimov.Draw("E0 SAME")
	
	plotting_helpers.drawNotes(vtx_alg, scalelumi)
	leg01.Draw("SAME")
	c_asimov.SaveAs(output_dir + "Asimov_plot_{}_{}".format(vtx_alg,variable) + '.pdf')

	outFile.Close()

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


