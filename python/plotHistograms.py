# Signal Efficiency Plotting Script 
import argparse, os, math, ROOT, glob, uproot, time


import numpy as np
from ROOT import*
from pylab import *
gROOT.LoadMacro("AtlasStyle.C")
gROOT.LoadMacro("AtlasUtils.C")
gROOT.LoadMacro("AtlasLabels.C")


histos  = '/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/python/histograms'
histos_savepath = '/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/python/'
file = ROOT.TFile(histos)
ROOT.gROOT.SetBatch(True)
SetAtlasStyle()

def plot():
	hpt = ROOT.TH1D("hpt","hpt",300,0,150)
  	hpt.GetXaxis().SetTitle("p_{T} [GeV]")
  	hpt.GetYaxis().SetTitle("entries")
  	hpt.Sumw2()

  	hpt_wrtSV= ROOT.TH1D("hpt_wrtSV","hpt_wrtSV",300,0,150)
  	hpt_wrtSV.Sumw2()

  	hpt = file.Get("track_pt")
  	hpt_wrtSV = file.Get("track_pt_wrtSV")


  	MyC01= ROOT.TCanvas("MyC01","",600,400)
	MyC01.Divide(1,1)
	MyC01.cd(1)

	leg01 = ROOT.TLegend(0.70,0.8,0.91,0.92)
  	leg01.SetTextSize(0.035)
  	leg01.SetBorderSize(0)
  	leg01.SetFillColor(kWhite)
  	leg01.SetShadowColor(kWhite)

  	leg01.AddEntry(hpt,"track_pt","l")
	leg01.AddEntry(hpt_wrtSV,"track_pt_wrtSV","l")
	
	
	hpt.SetTitle("")
	hpt.SetStats(0)

	
	# hrdv.Rebin(2)
	
	hpt.SetLineColor(kRed)
	hpt.Rebin(0.5)
	hpt.GetXaxis().SetTitle("p_{T} [GeV]")
  	hpt.GetYaxis().SetTitle("entries")
	hpt.GetXaxis().SetRange(0,40)
	hpt.Draw("HIST")
 	# hpt.SetLineWidth(0)

 	hpt_wrtSV.SetLineColor(kBlue)
 	hpt_wrtSV.Rebin(0.5)
 	# hpt.SetLineWidth(0)

	leg01.Draw()
	hpt_wrtSV.Draw("HIST SAME")
	MyC01.SaveAs(histos_savepath +'trk_pt'+'.pdf')
	# hrdv.GetXaxis().SetRange(0,hrdv_bin_xmax+1) # add one bin for aesthetics
	# hrdv.GetYaxis().SetRangeUser(0,hrdv_ymax*1.7)
	# drawNotes(MC_campaign,DV_type,mass,lifetime)




if __name__ == '__main__':



	plot()