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

def plot_cutflow(file, vertextype, outputDir="/eos/home-r/rnewhous/public/plots/DHNLNtupleAnalysis/"):

    Tfile = ROOT.TFile(file)
    if vertextype == "VSI":
        hcutflow = Tfile.Get('CutFlow_VSI')
    elif vertextype == "VSI Leptons":
        hcutflow = Tfile.Get('CutFlow_VSI_Leptons')
    elif vertextype == "VSI Kaons":
        hcutflow = Tfile.Get('CutFlow_VSI_2')

    if not hcutflow: return

    # hcutflow.GetXaxis().SetBinLabel(1, "all")
    # hcutflow.GetXaxis().SetBinLabel(2, "PV")
    # hcutflow.GetXaxis().SetBinLabel(3, "Kshort mass")
    # hcutflow.GetXaxis().SetBinLabel(4, "alpha")

    MyC01 = ROOT.TCanvas("MyC01", "cutflow", 1200, 800)

    ymax_cutflow = hcutflow.GetMaximum()
    hcutflow.GetYaxis().SetRangeUser(0, ymax_cutflow*1.15)
    hcutflow.GetXaxis().SetRangeUser(0, hcutflow.FindLastBinAbove(0, 1))
    hcutflow.GetXaxis().SetLabelSize(0.05)

    hcutflow.SetFillColor(kAzure-4)
    hcutflow.SetLineWidth(0)
    hcutflow.Draw("HIST TEXT0 SAME")

    if "data" in file:
        plotting_helpers.drawNotesData("data18 period B", vertextype, 1, 'KShort')
    if "uuu" in file:
        plotting_helpers.drawNotes("mumu", "muon", vertextype)
    elif "ueu" in file:
        plotting_helpers.drawNotes("emu", "muon", vertextype)
    elif "uee" in file:
        plotting_helpers.drawNotes("ee", "muon", vertextype)
    elif "eee" in file:
        plotting_helpers.drawNotes("ee", "electron", vertextype)
    elif "eeu" in file:
        plotting_helpers.drawNotes("emu", "electron", vertextype)
    elif "euu" in file:
        plotting_helpers.drawNotes("mumu", "electron", vertextype)
    else:
        plotting_helpers.drawNote(vertextype, size=50, ax=0.7)

    channel = file.split("_")[1].split(".")[0]

    if vertextype == "VSI":
        savefilename = "CutFlow_VSI_" + channel
    elif vertextype == "VSI Leptons":
        savefilename = "CutFlow_VSI_Leptons_" + channel
    elif vertextype == "VSI Kaons":
        savefilename = "CutFlow_VSI_Kaons_" + channel

    # MyC01.SaveAs(outputDir +savefilename+'.pdf')
    MyC01.SaveAs(outputDir + savefilename+'.png')


def compareN(file, hname, hlabel, savefilename, vertextype, setxrange="", scaleymax=1, nRebin=1, setlogy=False, outputDir="../output/"):
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
    elif vertextype == "VSI Kaons":
        vertexSuffix = "_VSI_2"

    for i in xrange(nhist):
        h[i] = f.Get(hname[i]+vertexSuffix)
    # exit()

    # define your canvas
    MyC01 = ROOT.TCanvas("MyC01", "", 600, 400)
    # MyC01.Divide(1,1)
    # MyC01.cd(1)

    # format legend
    leg01 = ROOT.TLegend(0.60, 0.7, 0.91, 0.92)
    leg01.SetTextSize(0.035)
    leg01.SetBorderSize(0)
    leg01.SetFillColor(kWhite)
    leg01.SetShadowColor(kWhite)

    ymax_list = []
    h_binxmax_list = []
    h_binxmin_list = []
    for i in xrange(nhist):
        leg01.AddEntry(h[i], hlabel[i], "lp")
        h[i].Rebin(nRebin)

        ymax_list.append(h[i].GetMaximum())
        h_binxmax_list.append(h[i].FindLastBinAbove(0, 1))
        h_binxmin_list.append(h[i].FindFirstBinAbove(0, 1))

    y_max = max(ymax_list)
    bin_xmax = max(h_binxmax_list)
    bin_xmin = min(h_binxmin_list)

    shapelist = [20, 22, 21, 33, 29]
    for i in xrange(nhist):
        h[i].SetLineColor(plotting_helpers.histColours(i))
        h[i].SetMarkerColor(plotting_helpers.histColours(i))
        h[i].SetMarkerSize(0.7)
        h[i].SetMarkerStyle(shapelist[i])
        if i == 0:
            h[i].GetXaxis().SetTitle(plotting_helpers.xlabelhistograms(hname[i]))
            h[i].GetYaxis().SetTitle("entries")
            h[i].GetYaxis().SetRangeUser(0, y_max*scaleymax)
            if setxrange != "":
                X_minmax = [item for item in setxrange.split(' ')]
                h[i].GetXaxis().SetRangeUser(
                    int(X_minmax[0]), int(X_minmax[1]))
            else:
                h[i].GetXaxis().SetRange(bin_xmin-1, bin_xmax+1)
            h[i].Draw("E0 HIST")
        else:
            h[i].Draw("E0 HIST SAME")

    if setlogy:
        gPad.SetLogy()

    leg01.Draw()
    atlas_style.ATLASLabel(0.25, 0.87, "Internal")

    # if "pmu" in hname[0]:
    #     if "mumu" in hname[0]:
    #         plotting_helpers.drawNotesMC("",vertextype, "mumumu","10","10")
    #     # if "mumu" in hname[0]:
    #     #     plotting_helpers.drawNotesMC("",vertextype, "mumu","10","10")
    #     if "emu" in hname[0]:
    #         plotting_helpers.drawNotesMC("",vertextype, "emumu","10","10")

    if "data" in file:
        plotting_helpers.drawNotesData("data18 period B", vertextype)
    if "uuu" in file:
        plotting_helpers.drawNotes("mumu", "muon", vertextype)
    elif "ueu" in file:
        plotting_helpers.drawNotes("emu", "muon", vertextype)
    elif "uee" in file:
        plotting_helpers.drawNotes("ee", "muon", vertextype)
    elif "eee" in file:
        plotting_helpers.drawNotes("ee", "electron", vertextype)
    elif "eeu" in file:
        plotting_helpers.drawNotes("emu", "electron", vertextype)
    elif "euu" in file:
        plotting_helpers.drawNotes("mumu", "electron", vertextype)
    else:
        plotting_helpers.drawNotesVertextype(vertextype)

    if "DV_mass" in hname[0]:
        mDVcut = TLine(4, 0, 4, y_max)
        mDVcut.SetLineStyle(3)
        mDVcut.Draw("SAME")

    if "mvis" in hname[0]:
        min_mlll = TLine(50, 0, 50, y_max)
        min_mlll.SetLineStyle(3)
        min_mlll.Draw("SAME")

        max_mlll = TLine(84, 0, 84, y_max)
        max_mlll.SetLineStyle(3)
        max_mlll.Draw("SAME")

    if "DV_r" in hname[0]:
        if "redmass" in hname[0]:
            pass
        else:
            matlayers = [33.25, 50.5, 88.5, 122.5, 299]
            nmatlayers = len(matlayers)
            matlay = {}
            for i in range(nmatlayers):
                matlay[i] = TLine(matlayers[i], 0, matlayers[i], y_max)
                matlay[i].SetLineStyle(3)
                matlay[i].Draw("SAME")
            leg01.AddEntry(matlay[0], "material layers", "l")

    channel = file.split("_")[1].split(".")[0]

    if vertextype == "VSI":
        name = savefilename + "_" + channel
    elif vertextype == "VSI Leptons":
        name = savefilename + "_" + channel
    elif vertextype == "VSI Kaons":
        name = savefilename + "_" + channel

    # MyC01.SaveAs(outputDir +savefilename+'.pdf')
    MyC01.SaveAs(outputDir + savefilename+'.png')


def compare(hist_channels, variable="", setrange=None, scaleymax=1.9, nRebin=1, setlogy=False, outputDir="/eos/home-r/rnewhous/public/plots/DHNLNtupleAnalysis/",
            save_name="", vertical_lines=[], labels=[], norm=0, **kwargs):

    histograms = []
    tfiles = {} # root is stupid and will close the file if you're not careful
    for key, val in hist_channels.items():
        tfiles[key] = ROOT.TFile(val[0]) # get file
        if not variable: 
            histname = "{}_{}".format(val[3], val[1])
        else: 
            histname = "{}_{}".format(variable, val[1])
            if val[2]:  histname = "{}_{}".format(val[2], histname)
        histograms.append(tfiles[key].Get(histname)) # get variable with suffix

    n_h = len(histograms)
    h_idx = range(len(histograms))
    hmc = {}
    mc = {}

    # define your canvas
    MyC01 = ROOT.TCanvas("MyC01", "", 1200, 800)

    # format legend
    leg01 = ROOT.TLegend(0.50, 0.7, 0.91, 0.92)
    leg01.SetTextSize(0.035)
    leg01.SetBorderSize(0)
    leg01.SetFillColor(kWhite)
    leg01.SetShadowColor(kWhite)

    for i in h_idx:
        if not histograms[i]: # no histogram object. don't even try
            print( "cannot find {}. Exiting".format(histname))
            return
        leg01.AddEntry(histograms[i], hist_channels.keys()[i], "lp")
        histograms[i].Rebin(nRebin)

    # find the common min and max for x axis
    if setrange == None:
        bin_xmin, bin_xmax = (np.inf, -np.inf)
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

    # set the common x limits for all histograms
    for i in h_idx:
        histograms[i].GetXaxis().SetRangeUser(x_min, x_max)

    # normalize the histograms

    if norm > 0:
        for i in h_idx:
            if (histograms[i].Integral() != 0):
                scale_mc = norm/(histograms[i].Integral(bin_xmin, bin_xmax))
            else:
                scale_mc = norm
            histograms[i].Scale(scale_mc)

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
        histograms[i].SetLineColor(plotting_helpers.histColours(i))
        histograms[i].SetLineWidth(3)
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
        leg01.AddEntry(lines[0], kwargs['vertical_legend'], "l")

    leg01.Draw()
    notes_x = 0.25
    notes_y = 0.87
    atlas_style.ATLASLabel(notes_x, notes_y, "Internal")
    if 'notes' in kwargs:
        for note in kwargs['notes']:
            notes_y -= 0.07
            plotting_helpers.drawNote(note, size=40, ax=notes_x, ay=notes_y)

    # plotting_helpers.getNote(35).DrawLatex(notes_x, notes_y-.05, vertextype)

    save_file_name = "{}_{}".format(val[2], variable if save_name == "" else save_name)
    if not os.path.exists(outputDir): os.mkdir(outputDir)
    # MyC01.SaveAs(outputDir +savefilename+'.pdf')
    MyC01.SaveAs(outputDir + save_file_name + '.png')


def compare_dataMC(datafile, mcfiles, hname, hdatalabel, hmclabels, vertextype, setrange=None,
                   scaleymax=1.2, nRebin=1, setlogy=False, outputDir="../output/", save_name="", vertical_lines=[], **kwargs):

    # get 2 histograms from input file

    if datafile:
        data = ROOT.TFile(datafile)
        if vertextype == "VSI":
            # need to remake histograms with new format ;)
            hdata = data.Get(hname + "_VSI")
            # hdata = data.Get(hname + "_pmu_mumu_VSI") # need to remake histograms with new format ;)
        elif vertextype == "VSI Leptons":
            hdata = data.Get(hname + "_VSI_Leptons")  # need to fix this ;)
            # hdata = data.Get(hname + "_pmu_mumu_VSILep") # need to fix this ;)
        elif vertextype == "VSI Kaons":
            hdata = data.Get(hname + "_VSI_2")  # need to fix this ;)
            # hdata = data.Get(hname + "_pmu_mumu_VSILep") # need to fix this ;)
        else:
            logger.error("Couldn't find the data  histogram you requested!")
            logger.error(
                "Check file %s has the histogram you are looking for!" % datafile)
            exit()

    nmc_files = len(mcfiles)
    hmc = {}
    mc = {}

    for i in range(nmc_files):
        mc[i] = ROOT.TFile(mcfiles[i])
        if vertextype == "VSI":
            hmc[i] = mc[i].Get(hname + "_VSI")
        elif vertextype == "VSI Leptons":
            hmc[i] = mc[i].Get(hname + "_VSI_Leptons")
        elif vertextype == "VSI Kaons":
            hmc[i] = mc[i].Get(hname + "_VSI_2")
        else:
            logger.error("Couldn't find the mc histogram you requested!")
            logger.error(
                "Check file %s has the histogram you are looking for!" % mcfiles[i])
            exit()

    # define your canvas
    MyC01 = ROOT.TCanvas("MyC01", "", 1200, 800)

    # format legend
    leg01 = ROOT.TLegend(0.50, 0.7, 0.91, 0.92)
    leg01.SetTextSize(0.035)
    leg01.SetBorderSize(0)
    leg01.SetFillColor(kWhite)
    leg01.SetShadowColor(kWhite)

    if (datafile):
        leg01.AddEntry(hdata, hdatalabel, "lp")

        # rebin histograms
        hdata.Rebin(nRebin)

    for i in range(nmc_files):
        leg01.AddEntry(hmc[i], hmclabels[i], "lp")
        hmc[i].Rebin(nRebin)

    # find the common min and max for x axis
    if setrange == None:
        bin_xmin, bin_xmax = (-np.inf, np.inf)
        if (datafile):
            bin_xmax = hdata.FindLastBinAbove(0, 1)
            bin_xmin = hdata.FindFirstBinAbove(0, 1)
        for i in xrange(nmc_files):
            if hmc[i].FindLastBinAbove(0, 1) < bin_xmax:
                bin_xmax = hmc[i].FindLastBinAbove(0, 1)
            if hmc[i].FindFirstBinAbove(0, 1) > bin_xmin:
                bin_xmin = hmc[i].FindFirstBinAbove(0, 1)
            x_min = bin_xmin-1
            x_max = bin_xmax+1
    else:
        x_min, x_max = setrange
        if (datafile):
            bin_xmax = hdata.GetXaxis().FindBin(x_max)
            bin_xmin = hdata.GetXaxis().FindBin(x_min)
        else:
            bin_xmax = hmc[0].GetXaxis().FindBin(x_max)
            bin_xmin = hmc[0].GetXaxis().FindBin(x_min)

    # set the common x limits for all histograms
    if (datafile):
        hdata.GetXaxis().SetRangeUser(x_min, x_max)
    for i in xrange(nmc_files):
        hmc[i].GetXaxis().SetRangeUser(x_min, x_max)

    # normalize the histograms
    norm = 1
    if (datafile):
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

    # find the common min and max for y axis
    if 'y_max' in kwargs:
        y_max = kwargs['y_max']
    else:
        y_max = -666666
        if (datafile):
            y_max = hdata.GetMaximum()
        for i in xrange(nmc_files):
            if hmc[i].GetMaximum() > y_max:
                y_max = hmc[i].GetMaximum()

    if (datafile):
        hdata.SetLineColor(kBlack)
        hdata.SetLineWidth(3)
        hdata.GetXaxis().SetTitle(plotting_helpers.xlabelhistograms(hname))
        hdata.GetYaxis().SetTitle("entries")
        hdata.GetYaxis().SetRangeUser(0, y_max*scaleymax)
        hdata.SetMarkerSize(1.5)

        hdata.Draw("E0 HIST")

    shapelist = [22, 21, 33, 29]
    for i in xrange(nmc_files):
        hmc[i].SetMarkerSize(1.5)
        hmc[i].SetLineColor(plotting_helpers.histColours(i))
        hmc[i].SetLineWidth(3)
        hmc[i].SetMarkerColor(plotting_helpers.histColours(i))
        hmc[i].SetMarkerStyle(shapelist[i])
        hmc[i].GetXaxis().SetTitle(plotting_helpers.xlabelhistograms(hname))
        hmc[i].GetYaxis().SetTitle("entries")
        hmc[i].GetYaxis().SetRangeUser(0, y_max*scaleymax)
        hmc[i].Draw("E0 HIST SAME")

    if setlogy:
        gPad.SetLogy()

    # draw vertical lines
    lines = []
    if 'material_layers' in kwargs and kwargs['material_layers']:
        vertical_lines.extend([33.25, 50.5, 88.5, 122.5, 299])

    for i, x in enumerate(vertical_lines):
        lines.append(TLine(x, 0, x, y_max))
        lines[i].SetLineStyle(3)
        lines[i].SetLineWidth(3)
        lines[i].Draw("SAME")

    if 'vertical_legend' in kwargs:
        leg01.AddEntry(lines[0], kwargs['vertical_legend'], "l")

    if 'material_layers' in kwargs and kwargs['material_layers']:
        leg01.AddEntry(lines[0], "material layers", "l")

    leg01.Draw()
    notes_x = 0.25
    notes_y = 0.87
    atlas_style.ATLASLabel(notes_x, notes_y, "Internal")
    plotting_helpers.getNote(35).DrawLatex(notes_x, notes_y-.05, vertextype)

    save_name = hname if save_name == "" else save_name
    if vertextype == "VSI":
        savefilename = save_name + "_compare_dataMC_VSI"
    elif vertextype == "VSI Leptons":
        savefilename = save_name + "_compare_dataMC_VSILep"
    elif vertextype == "VSI Kaons":
        savefilename = save_name + "_compare_dataMC_VSIKaons"
    else:
        savefilename = save_name + "_compare_dataMC"

    # MyC01.SaveAs(outputDir +savefilename+'.pdf')
    MyC01.SaveAs(outputDir + savefilename+'.png')


def compare2_wRatio(histos1, histos2, h1name, h1label, h2name, h2label, xlabel, savefilename):
    ############################################################################
    nRebin = 10  # set to 1 if you dont want to rebin.
    scaleymax = 1.1  # use this to scale height of y axis for asthetics
    ############################################################################

    # get 2 histograms from the file
    file1 = ROOT.TFile(histos1)
    file2 = ROOT.TFile(histos2)
    h1 = file1.Get(h1name)
    h2 = file2.Get(h2name)

    print h1.GetEntries()
    print h2.GetEntries()

    # define your canvas
    # MyC01 = ROOT.TCanvas("MyC01","",600,400)
    MyC01 = ROOT.TCanvas("MyC01", "canvas", 800, 800)
    MyC01.Divide(1, 1)
    MyC01.cd(1)
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)
    pad1.SetGridx()
    pad1.Draw()
    pad1.cd()

    # format legend
    leg01 = ROOT.TLegend(0.60, 0.7, 0.91, 0.92)
    leg01.SetTextSize(0.035)
    leg01.SetBorderSize(0)
    leg01.SetFillColor(kWhite)
    leg01.SetShadowColor(kWhite)

    leg01.AddEntry(h1, h1label, "l")
    leg01.AddEntry(h2, h2label, "l")

    # rebin histograms
    h1.Rebin(nRebin)
    h2.Rebin(nRebin)

    # find the max of the 2 histograms
    y1_max = h1.GetMaximum()
    y2_max = h1.GetMaximum()
    y_max = max(y1_max, y2_max)  # scale the max for asthetics

    h1_binxmax = h1.FindLastBinAbove(0, 1)
    h2_binxmax = h2.FindLastBinAbove(0, 1)
    bin_xmax = max(h1_binxmax, h2_binxmax)

    h1_binxmin = h1.FindFirstBinAbove(0, 1)
    h2_binxmin = h2.FindFirstBinAbove(0, 1)
    bin_xmin = min(h1_binxmin, h2_binxmin)

    h1.SetLineColor(kAzure+6)
    h1.GetXaxis().SetTitle(xlabel)
    h1.GetYaxis().SetTitle("entries")
    h1.GetYaxis().SetRangeUser(0, y_max*scaleymax)
    # h1.GetXaxis().SetRangeUser(0,5)
    h1.GetXaxis().SetRange(bin_xmin-1, bin_xmax+1)
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
    atlas_style.ATLASLabel(0.25, 0.87, "Internal")

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
    hratio.GetYaxis().SetLabelFont(43)  # Absolute font size in pixel (precision 3)
    hratio.GetYaxis().SetLabelSize(15)
    hratio.GetXaxis().SetTitleSize(20)
    hratio.GetXaxis().SetTitleFont(43)
    hratio.GetXaxis().SetTitleOffset(4.)
    hratio.GetXaxis().SetLabelFont(43)  # Absolute font size in pixel (precision 3)
    hratio.GetXaxis().SetLabelSize(15)

    hratio.Draw("ep")

    # MyC01.SaveAs(histos_savepath +savefilename+'.pdf')
    MyC01.SaveAs(histos_savepath + savefilename+'.png')


def CorrPlot2D(file, hname, hlabel, vertextype, setxrange="", setyrange="", rebinx=1, rebiny=1, outputDir="../output/"):

    f = ROOT.TFile(file)
    if "data" in file:
        if vertextype == "VSI":
            h = f.Get(hname+"_VSI")  # get histogram
            # h = f.Get(hname+"_pmu_mumu_VSI") # get histogram
        elif vertextype == "VSI Leptons":
            h = f.Get(hname+"_VSI_Leptons")  # get histogram
        elif vertextype == "VSI Kaons":
            h = f.Get(hname+"_VSI_2")  # get histogram
            # h = f.Get(hname+"_pmu_mumu_VSILep") # get histogram
    else:
        if vertextype == "VSI":
            h = f.Get(hname+"_VSI")  # get histogram
        elif vertextype == "VSI Leptons":
            h = f.Get(hname+"_VSI_Leptons")  # get histogram
        elif vertextype == "VSI Kaons":
            h = f.Get(hname+"_VSI_2")  # get histogram

    # define your canvas
    MyC01 = ROOT.TCanvas("MyC01", "", 800, 500)
    MyC01.SetRightMargin(20)

    # format legend
    leg01 = ROOT.TLegend(0.550, 0.8, 0.89, 0.92)
    leg01.SetTextSize(0.035)
    leg01.SetBorderSize(0)
    leg01.SetFillColor(kWhite)
    leg01.SetShadowColor(kWhite)
    leg01.AddEntry(h, hlabel, "")

    # Rebin 2D histogram!
    hold = h.Clone()
    hold.SetDirectory(0)
    nbinsx = hold.GetXaxis().GetNbins()
    nbinsy = hold.GetYaxis().GetNbins()
    xmin = hold.GetXaxis().GetXmin()
    xmax = hold.GetXaxis().GetXmax()
    ymin = hold.GetYaxis().GetXmin()
    ymax = hold.GetYaxis().GetXmax()
    nx = int(nbinsx/rebinx)
    ny = int(nbinsy/rebiny)
    h.SetBins(nx, xmin, xmax, ny, ymin, ymax)

    # loop on all bins to reset contents and errors
    for biny in range(1, nbinsy):
        if (biny <= nbinsy):
            for binx in range(1, nbinsx):
                if binx <= nbinsx:
                    ibin = h.GetBin(binx, biny)
                    h.SetBinContent(ibin, 0)

   # //loop on all bins and refill
    for biny in range(1, nbinsy):
        if (biny <= nbinsy):
            by = hold.GetYaxis().GetBinCenter(biny)
            iy = h.GetYaxis().FindBin(by)
            for binx in range(1, nbinsx):
                if binx <= nbinsx:
                    bx = hold.GetXaxis().GetBinCenter(binx)
                    ix = h.GetXaxis().FindBin(bx)
                    Bin = hold.GetBin(binx, biny)
                    ibin = h.GetBin(ix, iy)
                    cu = hold.GetBinContent(Bin)
                    h.AddBinContent(ibin, cu)

    if setxrange != "":
        X_minmax = [item for item in setxrange.split(' ')]
        h.GetXaxis().SetRangeUser(int(X_minmax[0]), int(X_minmax[1]))

    if setyrange != "":
        Y_minmax = [item for item in setyrange.split(' ')]
        h.GetYaxis().SetRangeUser(int(Y_minmax[0]), int(Y_minmax[1]))

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
    atlas_style.ATLASLabel(0.25, 0.87, "Internal")

    if "uuu" in file:
        plotting_helpers.drawNotes("mumu", "muon", vertextype)
    elif "ueu" in file:
        plotting_helpers.drawNotes("emu", "muon", vertextype)
    elif "uee" in file:
        plotting_helpers.drawNotes("ee", "muon", vertextype)
    elif "eee" in file:
        plotting_helpers.drawNotes("ee", "electron", vertextype)
    elif "eeu" in file:
        plotting_helpers.drawNotes("emu", "electron", vertextype)
    elif "euu" in file:
        plotting_helpers.drawNotes("mumu", "electron", vertextype)
    else:
        plotting_helpers.drawNotesVertextype(vertextype)

    if vertextype == "VSI":
        savefilename = hname + "_2Dmass_VSI"
    elif vertextype == "VSI Leptons":
        savefilename = hname + "_2Dmass_VSILep"
    elif vertextype == "VSI Kaons":
        savefilename = hname + "_2Dmass_VSIKaons"
    else:
        savefilename = hname + "_2Dmass"

    # MyC01.SaveAs(outputDir +savefilename+'.pdf')
    MyC01.SaveAs(outputDir + savefilename+'.png')
