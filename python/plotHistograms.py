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
	if "mlll" in hist:
		return "tri-lepton mass [GeV]"
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
		return kBlue





def compareN(file, hname, hlabel,savefilename):
	############################################################################
	nRebin = 5 # set to 1 if you dont want to rebin.
	scaleymax = 1.6 # use this to scale height of y axis for asthetics
	############################################################################


	# get 2 histograms from input file
	nhist = len(hname)
	f = ROOT.TFile(file)
	# file2 = ROOT.TFile(histos)
  	
  	h = {}
  	for i in xrange(nhist): 
	  	h[hname[i]] = f.Get(hname[i])


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


  	ymax_list = []
  	h_binxmax_list = []
  	h_binxmin_list = []
  	for i in xrange(nhist):
	  	leg01.AddEntry(h[hname[i]],hlabel[i],"l")
	  	h[hname[i]].Rebin(nRebin)

		ymax_list.append(h[hname[i]].GetMaximum())
		h_binxmax_list.append(h[hname[i]].FindLastBinAbove(0,1))
		h_binxmin_list.append(h[hname[i]].FindFirstBinAbove(0,1))

	y_max = max(ymax_list)
	bin_xmax = max(h_binxmax_list)
	bin_xmin = min(h_binxmin_list)


	for i in xrange(nhist): 
		if i == 0:
			h[hname[i]].SetLineColor(histColours(i))
			h[hname[i]].GetXaxis().SetTitle(xlabelhistograms(hname[i]))
		  	h[hname[i]].GetYaxis().SetTitle("entries")
		  	h[hname[i]].GetYaxis().SetRangeUser(0,y_max*scaleymax)
			h[hname[i]].GetXaxis().SetRangeUser(0,3000)
			h[hname[i]].GetXaxis().SetRange(bin_xmin-1,bin_xmax+1)
			h[hname[i]].Draw("HIST")
		else: 
			h[hname[i]].SetLineColor(histColours(i))
			h[hname[i]].Draw("HIST SAME")


	leg01.Draw()
	ATLASLabel(0.25,0.87,"Internal")
	
	MyC01.SaveAs(histos_savepath +savefilename+'.pdf')



def compare_dataMC(datafile,mcfile, nRebin, hdataname, hdatalabel, hmcname, hmclabel, setrange,scaleymax, vertextype,savefilename):
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
	MyC01.Divide(1,1)
	MyC01.cd(1)

	#format legend
	leg01 = ROOT.TLegend(0.60,0.7,0.91,0.92)
  	leg01.SetTextSize(0.035)
  	leg01.SetBorderSize(0)
  	leg01.SetFillColor(kWhite)
  	leg01.SetShadowColor(kWhite)

  	leg01.AddEntry(hdata,hdatalabel,"l")
	leg01.AddEntry(hmc,hmclabel,"l")



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




	hdata.SetLineColor(kAzure+6)
	hdata.GetXaxis().SetTitle(xlabelhistograms(hdataname))
  	hdata.GetYaxis().SetTitle("entries")
  	hdata.GetYaxis().SetRangeUser(0,y_max*scaleymax)


  	if setrange == "":
		hdata.GetXaxis().SetRange(bin_xmin-1,bin_xmax+1)
	else: 

		minmax = [item for item in setrange.split(' ')]
		hdata.GetXaxis().SetRangeUser(int(minmax[0]),int(minmax[1]))
	
	hdata.Draw("HIST")

	
 	hmc.SetLineColor(kRed)
 	# hmc.SetLineStyle(3)
 	hmc.Draw("HIST SAME")

 	if "mass" in hdataname:
	 	mDVcut=TLine(4,0,4,y_max)
	 	mDVcut.SetLineStyle(3)
		mDVcut.Draw("SAME")
	if "mlll" in hdataname:
		min_mlll=TLine(50,0,50,y_max)
	 	min_mlll.SetLineStyle(3)
		min_mlll.Draw("SAME")

		max_mlll=TLine(84,0,84,y_max)
	 	max_mlll.SetLineStyle(3)
		max_mlll.Draw("SAME")

	leg01.Draw()
	ATLASLabel(0.25,0.87,"Internal")
	if "pmu" in hdataname: 
		if "mumu" in hdataname:
  			helpers.drawNotes("mumu","muon",vertextype) 
  		if "emu" in hdataname:
  			helpers.drawNotes("emu","muon",vertextype) 
	
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
						help="data inputfile ",
						metavar="FILE")


	parser.add_argument("-f2", "--file2",
						dest="file2", default=[""],
						action = AppendActionCleanDefault,
						type = str,
						help="Input file 2",
						metavar="FILE")

	parser.add_argument("--config",
						dest="config",
						type = str,
						required = True,
						help="Input config file for plotHisotgrams.py.")


	parent_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter, parents = [parser]) 

	options = parent_parser.parse_args()



	with open(options.config, 'r') as json_config:
		config_file = json.load(json_config) # load JSON config file that contains a channel name mapped to a list of selections

	# print options.file[0]

	files = [item for item in options.file[0].split(',')]
	# print files 


	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "0 20",
					scaleymax = 1.2,				
					hdataname = "SSDV_mass_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "SSDV_mass_pmu_mumu_VSILep",
					hmclabel = "MC 10G 10mm",
					vertextype = "VSI Leptons",
					savefilename = "DV_mass_MCdatacompare_VSILep"	)


	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "0 40",
					scaleymax = 1.2,				
					hdataname = "SSDV_mass_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "SSDV_mass_pmu_mumu_VSI",
					hmclabel = "MC 10G 10mm",
					vertextype = "VSI",
					savefilename = "DV_mass_MCdatacompare_VSI"	)

	compare_dataMC(files[1],files[0],
					nRebin = 5, 
					setrange= "",
					scaleymax = 1.2,
					hdataname = "SSDV_r_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "SSDV_r_pmu_mumu_VSILep",
					hmclabel = "MC 10G 10mm",
					vertextype = "VSI Leptons",
					savefilename = "DV_r_MCdatacompare_VSILep"	)


	compare_dataMC(files[1],files[0],
					nRebin = 5, 
					setrange= "",
					scaleymax = 1.2,
					hdataname = "SSDV_r_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "SSDV_r_pmu_mumu_VSI",
					hmclabel = "MC 10G 10mm",
					vertextype = "VSI",
					savefilename = "DV_r_MCdatacompare_VSI"	)


	compare_dataMC(files[1],files[0],
					nRebin = 2, 
					setrange= "0 100",
					scaleymax = 1.2,					
					hdataname = "SSDV_trk_pt_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "SSDV_trk_pt_pmu_mumu_VSILep",
					hmclabel = "MC 10G 10mm",
					vertextype = "VSI Leptons",
					savefilename = "trk_pt_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 2, 
					setrange= "0 100",
					scaleymax = 1.2,					
					hdataname = "SSDV_trk_pt_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "SSDV_trk_pt_pmu_mumu_VSI",
					hmclabel = "MC 10G 10mm",
					vertextype = "VSI",
					savefilename = "trk_pt_MCdatacompare_VSI"	)



	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "",
					scaleymax = 2.2,
					hdataname = "SSDV_trk_phi_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "SSDV_trk_phi_pmu_mumu_VSILep",
					hmclabel = "MC 10G 10mm",
					vertextype = "VSI Leptons",
					savefilename = "trk_phi_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "",
					scaleymax = 2.2,
					hdataname = "SSDV_trk_phi_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "SSDV_trk_phi_pmu_mumu_VSI",
					hmclabel = "MC 10G 10mm",
					vertextype = "VSI",
					savefilename = "trk_phi_MCdatacompare_VSI"	)


	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "",
					scaleymax = 2.2,
					hdataname = "SSDV_trk_eta_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "SSDV_trk_eta_pmu_mumu_VSILep",
					hmclabel = "MC 10G 10mm",
					vertextype = "VSI Leptons",
					savefilename = "trk_eta_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "",
					scaleymax = 2.2,
					hdataname = "SSDV_trk_eta_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "SSDV_trk_eta_pmu_mumu_VSI",
					hmclabel = "MC 10G 10mm",
					vertextype = "VSI",
					savefilename = "trk_eta_MCdatacompare_VSI"	)


	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "-10 10",
					scaleymax = 1.2,
					hdataname = "SSDV_trk_d0_pmu_mumu_VSILep", 
					hdatalabel= "data 2018 period B",
					hmcname = "SSDV_trk_d0_pmu_mumu_VSILep",
					hmclabel = "MC 10G 10mm",
					vertextype = "VSI Leptons",
					savefilename = "trk_d0_MCdatacompare_VSILep"	)

	compare_dataMC(files[1],files[0],
					nRebin = 1, 
					setrange= "-10 10",
					scaleymax = 1.2,
					hdataname = "SSDV_trk_d0_pmu_mumu_VSI", 
					hdatalabel= "data 2018 period B",
					hmcname = "SSDV_trk_d0_pmu_mumu_VSI",
					hmclabel = "MC 10G 10mm",
					vertextype = "VSI",
					savefilename = "trk_d0_MCdatacompare_VSI"	)


	# compare_dataMC(files[1],files[0],
	# 				nRebin = 5, 
	# 				setrange= "0 200",
	# 				scaleymax = 1.5,
	# 				hdataname = "precutDV_mlll_pmu_mumu_VSILep", 
	# 				hdatalabel= "data 2018 period B",
	# 				hmcname = "precutDV_mlll_pmu_mumu_VSILep",
	# 				hmclabel = "MC 10G 10mm",
	# 				vertextype = "VSI Leptons",
	# 				savefilename = "precutmlll_MCdatacompare_VSILep"	)

	# compare_dataMC(files[1],files[0],
	# 				nRebin = 5, 
	# 				setrange= "0 200",
	# 				scaleymax = 1.5,
	# 				hdataname = "precutDV_mlll_pmu_mumu_VSI", 
	# 				hdatalabel= "data 2018 period B",
	# 				hmcname = "precutDV_mlll_pmu_mumu_VSI",
	# 				hmclabel = "MC 10G 10mm",
	# 				vertextype = "VSI",
	# 				savefilename = "precutmlll_MCdatacompare_VSI"	)








	exit()
# def compareN(file, hname, hlabel,savefilename):

	for k, configs in config_file.items():
		for key in config_file[k]["channels"]: 
			print key
		# print config_file[k]["channels"].keys()


	compareN(options.file[0],["DVtypeDV_r_pel_ee_VSILep","DVtypeDV_r_pmu_ee_VSILep","DVtypeDV_r_pel_mumu_VSILep"],["eee","muee","mumu"], savefilename="testfile_MC" )
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








