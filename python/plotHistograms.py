# Plotting Script
import argparse, os, math, ROOT, glob, uproot, time, json
import numpy as np
import helpers
from ROOT import gROOT
from ROOT import gPad
from pylab import *
import plotting
import atlas_style

ROOT.gROOT.SetBatch(True)

# trying to  set ATLAS style a million different ways
gROOT.SetStyle("ATLAS")  # might have to change how you set atlas style like this, depends how you have setup python
# atlas_style.AtlasStyle()	
# gROOT.LoadMacro("AtlasStyle.C")
# gROOT.LoadMacro("AtlasUtils.C")
# gROOT.LoadMacro("AtlasLabels.C")
# SetAtlasStyle()

logger = helpers.getLogger('dHNLAnalysis.plotHisotgrams')

#############################################################################################################################################
# globals
outputDir = '../output/plots/'  # change path here to save your histograms somewhere else!
MATERIAL_LAYERS = [33.25, 50.5, 88.5, 122.5, 299]


#############################################################################################################################################


def compareMCdata(config_file, hist_type="all"):
    # add plots here if you want to compare different distributions with with data & MC
    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_mass",
                            vertical_lines=[4],
                            nRebin=100,
                            scaleymax=1,
                            setrange=(0, 20),
                            vertextype="VSI",
                            outputDir=outputDir)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_mass",
                            vertical_lines=[4],
                            nRebin=100,
                            setrange=(0, 20),
                            vertextype="VSI Kaons",
                            outputDir=outputDir)

    # zoomed x axis
    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_mass",
                            save_name=hist_type + "_DV_mass_zoomed",
                            scaleymax=2.2,
                            vertical_lines=[.4977 - .01, .4977 + .01],
                            setrange=(0.4, 0.6),
                            vertextype="VSI",
                            outputDir=outputDir,
                            vertical_legend="K_{S}^{0} mass cut",
                            y_max=0.3)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_mass",
                            save_name=hist_type + "_DV_mass_zoomed",
                            scaleymax=2.2,
                            vertical_lines=[.4977 - .01, .4977 + .01],
                            setrange=(0.4, 0.6),
                            vertextype="VSI Kaons",
                            outputDir=outputDir,
                            vertical_legend="K_{S}^{0} mass cut",
                            y_max=0.3)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_alpha",
                            setrange=(-0.5, 0.5),
                            scaleymax=2.2,
                            vertextype="VSI",
                            outputDir=outputDir)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_alpha",
                            setrange=(-0.5, 0.5),
                            scaleymax=2.2,
                            vertextype="VSI Kaons",
                            outputDir=outputDir)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_r",
                            nRebin=5,
                            vertextype="VSI Kaons",
                            outputDir=outputDir,
                            material_layers=True)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_r",
                            nRebin=5,
                            vertextype="VSI",
                            outputDir=outputDir,
                            material_layers=True)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_trk_pt",
                            nRebin=2,
                            setrange=(0, 100),
                            vertextype="VSI Kaons",
                            outputDir=outputDir)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_trk_pt",
                            nRebin=2,
                            setrange=(0, 100),
                            vertextype="VSI",
                            outputDir=outputDir)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_trk_phi",
                            scaleymax=2.2,
                            vertextype="VSI Kaons",
                            outputDir=outputDir)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_trk_phi",
                            scaleymax=2.2,
                            vertextype="VSI",
                            outputDir=outputDir)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_trk_eta",
                            scaleymax=2.2,
                            vertextype="VSI Kaons",
                            outputDir=outputDir)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_trk_eta",
                            scaleymax=2.2,
                            vertextype="VSI",
                            outputDir=outputDir)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_trk_d0",
                            setrange=(-10, 10),
                            vertextype="VSI Kaons",
                            outputDir=outputDir)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_trk_d0",
                            setrange=(-10, 10),
                            vertextype="VSI",
                            outputDir=outputDir)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_trk_dpt",
                            nRebin=2,
                            setrange=(0, 20),
                            vertextype="VSI Kaons",
                            outputDir=outputDir)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_trk_dpt",
                            nRebin=2,
                            setrange=(0, 20),
                            vertextype="VSI",
                            outputDir=outputDir)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_trk_deta",
                            nRebin=2,
                            setrange=(0, 4),
                            vertextype="VSI Kaons",
                            outputDir=outputDir)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_trk_deta",
                            nRebin=2,
                            setrange=(0, 4),
                            vertextype="VSI",
                            outputDir=outputDir)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_trk_dphi",
                            setrange=(0, 8),
                            vertextype="VSI Kaons",
                            outputDir=outputDir)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_trk_dphi",
                            setrange=(0, 8),
                            vertextype="VSI",
                            outputDir=outputDir)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_trk_dR",
                            nRebin=10,
                            setrange=(0, 2),
                            vertextype="VSI Kaons",
                            outputDir=outputDir)

    plotting.compare_dataMC(datafile=config_file["dataFile"],
                            mcfiles=config_file["mcFile"],
                            hdatalabel=config_file["dataLabel"],
                            hmclabels=config_file["mcLabel"],
                            hname=hist_type + "_DV_trk_dR",
                            nRebin=10,
                            setrange=(0, 10),
                            vertextype="VSI",
                            outputDir=outputDir)


# plotting.compare_dataMC(datafile=config_file["dataFile"],
# 						mcfiles=config_file["mcFile"],
# 						hdatalabel= config_file["dataLabel"],
# 						hmclabels = config_file["mcLabel"],
# 						hname = hist_type+"_mvis",
# 						vertical_lines=[50, 84],
# 						nRebin = 5,
# 						setrange= (0, 200),
# 						scaleymax = 1.5,
# 						vertextype = "VSI Kaons",
# 						outputDir= outputDir)

# plotting.compare_dataMC(datafile=config_file["dataFile"],
# 						mcfiles=config_file["mcFile"],
# 						hdatalabel= config_file["dataLabel"],
# 						hmclabels = config_file["mcLabel"],
# 						hname = hist_type+"_mvis",
# 						vertical_lines=[50, 84],
# 						nRebin = 5,
# 						setrange= (0, 200),
# 						scaleymax = 1.5,
# 						vertextype = "VSI",
# 						outputDir= outputDir)

# plotting.compare_dataMC(datafile=config_file["dataFile"],
# 						mcfiles=config_file["mcFile"],
# 						hdatalabel= config_file["dataLabel"],
# 						hmclabels = config_file["mcLabel"],
# 						hname = hist_type+"_mtrans",
# 						nRebin = 5,
# 						setrange= (0, 200),
# 						scaleymax = 1.5,
# 						vertextype = "VSI Kaons",
# 						outputDir= outputDir)

# plotting.compare_dataMC(datafile=config_file["dataFile"],
# 						mcfiles=config_file["mcFile"],
# 						hdatalabel= config_file["dataLabel"],
# 						hmclabels = config_file["mcLabel"],
# 						hname = hist_type+"_mtrans",
# 						nRebin = 5,
# 						setrange= (0, 200),
# 						scaleymax = 1.5,
# 						vertextype = "VSI",
# 						outputDir= outputDir)


# plotting.compare_dataMC(datafile=config_file["dataFile"],
# 						mcfiles=config_file["mcFile"],
# 						hdatalabel= config_file["dataLabel"],
# 						hmclabels = config_file["mcLabel"],
# 						hname = "all_DV_redmass",
# 						setrange= (0, 20),
# 						scaleymax = 1.5,
# 						vertextype = "VSI Kaons",
# 						outputDir= outputDir)

# plotting.compare_dataMC(datafile=config_file["dataFile"],
# 						mcfiles=config_file["mcFile"],
# 						hdatalabel= config_file["dataLabel"],
# 						hmclabels = config_file["mcLabel"],
# 						hname = "all_DV_redmass",
# 						setrange= (0, 20),
# 						scaleymax = 1.5,
# 						vertextype = "VSI",
# 						outputDir= outputDir)


# plotting.compare_dataMC(datafile=config_file["dataFile"],
# 						mcfiles=config_file["mcFile"],
# 						hdatalabel= config_file["dataLabel"],
# 						hmclabels = config_file["mcLabel"],
# 						hname = "all_DV_redmassvis",
# 						nRebin = 5,
# 						setrange= (0, 200),
# 						scaleymax = 1.5,
# 						vertextype = "VSI Kaons",
# 						outputDir= outputDir)

# plotting.compare_dataMC(datafile=config_file["dataFile"],
# 						mcfiles=config_file["mcFile"],
# 						hdatalabel= config_file["dataLabel"],
# 						hmclabels = config_file["mcLabel"],
# 						hname = "all_DV_redmassvis",
# 						nRebin = 5,
# 						setrange= (0, 200),
# 						scaleymax = 1.5,
# 						vertextype = "VSI",
# 						outputDir= outputDir)

# plotting.compare_dataMC(datafile=config_file["dataFile"],
# 						mcfiles=config_file["mcFile"],
# 						hdatalabel= config_file["dataLabel"],
# 						hmclabels = config_file["mcLabel"],
# 						hname = "all_HNLm",
# 						setrange= (0, 20),
# 						scaleymax = 1.5,
# 						vertextype = "VSI Kaons",
# 						outputDir= outputDir)

# plotting.compare_dataMC(datafile=config_file["dataFile"],
# 						mcfiles=config_file["mcFile"],
# 						hdatalabel= config_file["dataLabel"],
# 						hmclabels = config_file["mcLabel"],
# 						hname = "all_HNLm",
# 						setrange= (0, 20),
# 						scaleymax = 1.5,
# 						vertextype = "VSI",
# 						outputDir= outputDir)


# plotting.compare_dataMC(datafile=config_file["dataFile"],
# 						mcfiles=config_file["mcFile"],
# 						hdatalabel= config_file["dataLabel"],
# 						hmclabels = config_file["mcLabel"],
# 						hname = "all_DV_redmassHNL",
# 						setrange= (0, 50),
# 						scaleymax = 1.5,
# 						vertextype = "VSI Kaons",
# 						outputDir= outputDir)

# plotting.compare_dataMC(datafile=config_file["dataFile"],
# 						mcfiles=config_file["mcFile"],
# 						hdatalabel= config_file["dataLabel"],
# 						hmclabels = config_file["mcLabel"],
# 						hname = "all_DV_redmassHNL",
# 						nRebin = 1,
# 						setrange= (0, 50),
# 						scaleymax = 1.5,
# 						vertextype = "VSI",
# 						outputDir= outputDir)


# plotting.compare_dataMC(datafile=config_file["dataFile"],
# 						mcfiles=config_file["mcFile"],
# 						hdatalabel= config_file["dataLabel"],
# 						hmclabels = config_file["mcLabel"],
# 						hname = hist_type+"_HNLpt",
# 						nRebin = 5,
# 						setrange= (0, 200),
# 						scaleymax = 1.5,
# 						vertextype = "VSI Kaons",
# 						outputDir= outputDir)

# plotting.compare_dataMC(datafile=config_file["dataFile"],
# 						mcfiles=config_file["mcFile"],
# 						hdatalabel= config_file["dataLabel"],
# 						hmclabels = config_file["mcLabel"],
# 						hname = hist_type+"_HNLpt",
# 						nRebin = 5,
# 						setrange= (0, 200),
# 						scaleymax = 1.5,
# 						vertextype = "VSI",
# 						outputDir= outputDir)


# plotting.compare_dataMC(datafile=config_file["dataFile"],
# 						mcfiles=config_file["mcFile"],
# 						hdatalabel= config_file["dataLabel"],
# 						hmclabels = config_file["mcLabel"],
# 						hname = hist_type+"_HNLphi",
# 						setrange= (-4, 4),
# 						scaleymax = 1.5,
# 						vertextype = "VSI Kaons",
# 						outputDir= outputDir)

# plotting.compare_dataMC(datafile=config_file["dataFile"],
# 						mcfiles=config_file["mcFile"],
# 						hdatalabel= config_file["dataLabel"],
# 						hmclabels = config_file["mcLabel"],
# 						hname = hist_type+"_HNLphi",
# 						setrange= (-4, 4),
# 						scaleymax = 1.5,
# 						vertextype = "VSI",
# 						outputDir= outputDir)


# plotting.compare_dataMC(datafile=config_file["dataFile"],
# 						mcfiles=config_file["mcFile"],
# 						hdatalabel= config_file["dataLabel"],
# 						hmclabels = config_file["mcLabel"],
# 						hname = hist_type+"_HNLeta",
# 						setrange= (-3, 3),
# 						scaleymax = 1.5,
# 						vertextype = "VSI Kaons",
# 						outputDir= outputDir)

# plotting.compare_dataMC(datafile=config_file["dataFile"],
# 						mcfiles=config_file["mcFile"],
# 						hdatalabel= config_file["dataLabel"],
# 						hmclabels = config_file["mcLabel"],
# 						hname = hist_type+"_HNLeta",
# 						setrange= (-3, 3),
# 						scaleymax = 1.5,
# 						vertextype = "VSI",
# 						outputDir= outputDir)


def make2Dmassplots(config_file):
    # MC FILES
    plotting.CorrPlot2D(file=config_file["mcFile"][0],
                        hname="charge_DVmass_mvis",
                        hlabel=config_file["mcLabel"][0],
                        rebinx=2,
                        rebiny=5,
                        vertextype="VSI",
                        setxrange="0 30",
                        setyrange="0 300",
                        outputDir=outputDir + "2Dmassplots/mc/")

    plotting.CorrPlot2D(file=config_file["mcFile"][0],
                        hlabel=config_file["mcLabel"][0],
                        hname="charge_DVmass_mvis",
                        rebinx=2,
                        rebiny=5,
                        vertextype="VSI Kaons",
                        setxrange="0 30",
                        setyrange="0 300",
                        outputDir=outputDir + "2Dmassplots/mc/")

    plotting.CorrPlot2D(file=config_file["mcFile"][0],
                        rebinx=2,
                        rebiny=5,
                        hname="charge_DVmass_mtrans",
                        hlabel=config_file["mcLabel"][0],
                        vertextype="VSI",
                        setxrange="0 30",
                        setyrange="0 300",
                        outputDir=outputDir + "2Dmassplots/mc/")

    plotting.CorrPlot2D(file=config_file["mcFile"][0],
                        rebinx=2,
                        rebiny=5,
                        hname="charge_DVmass_mtrans",
                        hlabel=config_file["mcLabel"][0],
                        vertextype="VSI Kaons",
                        setxrange="0 30",
                        setyrange="0 300",
                        outputDir=outputDir + "2Dmassplots/mc/")

    plotting.CorrPlot2D(file=config_file["mcFile"][0],
                        rebinx=5,
                        rebiny=5,
                        hname="charge_mvis_mtrans",
                        hlabel=config_file["mcLabel"][0],
                        vertextype="VSI",
                        setxrange="0 100",
                        setyrange="0 300",
                        outputDir=outputDir + "2Dmassplots/mc/")

    plotting.CorrPlot2D(file=config_file["mcFile"][0],
                        rebinx=5,
                        rebiny=5,
                        hname="charge_mvis_mtrans",
                        hlabel=config_file["mcLabel"][0],
                        vertextype="VSI Kaons",
                        setxrange="0 100",
                        setyrange="0 300",
                        outputDir=outputDir + "2Dmassplots/mc/")

    plotting.CorrPlot2D(file=config_file["mcFile"][0],
                        hname="charge_DVmass_mhnl",
                        hlabel=config_file["mcLabel"][0],
                        vertextype="VSI",
                        setxrange="0 50",
                        setyrange="0 50",
                        outputDir=outputDir + "2Dmassplots/mc/")

    plotting.CorrPlot2D(file=config_file["mcFile"][0],
                        hname="charge_DVmass_mhnl",
                        hlabel=config_file["mcLabel"][0],
                        vertextype="VSI Kaons",
                        setxrange="0 50",
                        setyrange="0 50",
                        outputDir=outputDir + "2Dmassplots/mc/")

    plotting.CorrPlot2D(file=config_file["mcFile"][0],
                        rebinx=5,
                        rebiny=1,
                        hname="charge_mvis_mhnl",
                        hlabel=config_file["mcLabel"][0],
                        vertextype="VSI",
                        setxrange="0 300",
                        setyrange="0 50",
                        outputDir=outputDir + "2Dmassplots/mc/")

    plotting.CorrPlot2D(file=config_file["mcFile"][0],
                        rebinx=5,
                        rebiny=1,
                        hname="charge_mvis_mhnl",
                        hlabel=config_file["mcLabel"][0],
                        vertextype="VSI Kaons",
                        setxrange="0 300",
                        setyrange="0 50",
                        outputDir=outputDir + "2Dmassplots/mc/")

    plotting.CorrPlot2D(file=config_file["mcFile"][0],
                        rebinx=1,
                        rebiny=5,
                        hname="charge_mhnl_mtrans",
                        hlabel=config_file["mcLabel"][0],
                        vertextype="VSI",
                        setxrange="0 50",
                        setyrange="0 300",
                        outputDir=outputDir + "2Dmassplots/mc/")

    plotting.CorrPlot2D(file=config_file["mcFile"][0],
                        rebinx=1,
                        rebiny=5,
                        hname="charge_mhnl_mtrans",
                        hlabel=config_file["mcLabel"][0],
                        vertextype="VSI Kaons",
                        setxrange="0 50",
                        setyrange="0 300",
                        outputDir=outputDir + "2Dmassplots/mc/")

    # DATA FILES
    plotting.CorrPlot2D(config_file["dataFile"],
                        hname="charge_DVmass_mvis",
                        hlabel=config_file["dataLabel"],
                        rebinx=2,
                        rebiny=5,
                        vertextype="VSI",
                        setxrange="0 30",
                        setyrange="0 300",
                        outputDir=outputDir + "2Dmassplots/data/")

    plotting.CorrPlot2D(config_file["dataFile"],
                        hname="charge_DVmass_mvis",
                        hlabel=config_file["dataLabel"],
                        rebinx=2,
                        rebiny=5,
                        vertextype="VSI Kaons",
                        setxrange="0 30",
                        setyrange="0 300",
                        outputDir=outputDir + "2Dmassplots/data/")

    plotting.CorrPlot2D(config_file["dataFile"],
                        rebinx=2,
                        rebiny=5,
                        hname="charge_DVmass_mtrans",
                        hlabel=config_file["dataLabel"],
                        vertextype="VSI",
                        setxrange="0 30",
                        setyrange="0 300",
                        outputDir=outputDir + "2Dmassplots/data/")

    plotting.CorrPlot2D(config_file["dataFile"],
                        rebinx=2,
                        rebiny=5,
                        hname="charge_DVmass_mtrans",
                        hlabel=config_file["dataLabel"],
                        vertextype="VSI Kaons",
                        setxrange="0 30",
                        setyrange="0 300",
                        outputDir=outputDir + "2Dmassplots/data/")

    plotting.CorrPlot2D(config_file["dataFile"],
                        rebinx=5,
                        rebiny=5,
                        hname="charge_mvis_mtrans",
                        hlabel=config_file["dataLabel"],
                        vertextype="VSI",
                        setxrange="0 100",
                        setyrange="0 300",
                        outputDir=outputDir + "2Dmassplots/data/")

    plotting.CorrPlot2D(config_file["dataFile"],
                        rebinx=5,
                        rebiny=5,
                        hname="charge_mvis_mtrans",
                        hlabel=config_file["dataLabel"],
                        vertextype="VSI Kaons",
                        setxrange="0 100",
                        setyrange="0 300",
                        outputDir=outputDir + "2Dmassplots/data/")

    plotting.CorrPlot2D(config_file["dataFile"],
                        hname="charge_DVmass_mhnl",
                        hlabel=config_file["dataLabel"],
                        vertextype="VSI",
                        setxrange="0 50",
                        setyrange="0 50",
                        outputDir=outputDir + "2Dmassplots/data/")

    plotting.CorrPlot2D(config_file["dataFile"],
                        hname="charge_DVmass_mhnl",
                        hlabel=config_file["dataLabel"],
                        vertextype="VSI Kaons",
                        setxrange="0 50",
                        setyrange="0 50",
                        outputDir=outputDir + "2Dmassplots/data/")

    plotting.CorrPlot2D(config_file["dataFile"],
                        rebinx=5,
                        rebiny=1,
                        hname="charge_mvis_mhnl",
                        hlabel=config_file["dataLabel"],
                        vertextype="VSI",
                        setxrange="0 300",
                        setyrange="0 50",
                        outputDir=outputDir + "2Dmassplots/data/")

    plotting.CorrPlot2D(config_file["dataFile"],
                        rebinx=5,
                        rebiny=1,
                        hname="charge_mvis_mhnl",
                        hlabel=config_file["dataLabel"],
                        vertextype="VSI Kaons",
                        setxrange="0 300",
                        setyrange="0 50",
                        outputDir=outputDir + "2Dmassplots/data/")

    plotting.CorrPlot2D(config_file["dataFile"],
                        rebinx=1,
                        rebiny=5,
                        hname="charge_mhnl_mtrans",
                        hlabel=config_file["dataLabel"],
                        vertextype="VSI",
                        setxrange="0 50",
                        setyrange="0 300",
                        outputDir=outputDir + "2Dmassplots/data/")

    plotting.CorrPlot2D(config_file["dataFile"],
                        rebinx=1,
                        rebiny=5,
                        hname="charge_mhnl_mtrans",
                        hlabel=config_file["dataLabel"],
                        vertextype="VSI Kaons",
                        setxrange="0 50",
                        setyrange="0 300",
                        outputDir=outputDir + "2Dmassplots/data/")


def makeCutflows(config_file):
    # make Cutflow plots here
    plotting.plot_cutflow(file=config_file["mcFile"][0],
                          vertextype="VSI Kaons",
                          outputDir=outputDir + "Cutflows/")


def compareHistograms(config_file):
    # compare different histograms from the same file here
    plotting.compareN(file=config_file["mcFile"][0],
                      hname=["sel_DV_mass", "sel_mvis", "sel_HNLm", "sel_mtrans", "sel_DV_redmass"],
                      hlabel=["DV mass", "Visible mass", "m_{HNL}", "m_{T}", "Reduced mass"],
                      setxrange="0 100",
                      scaleymax=1.5,
                      nRebin=5,
                      vertextype="VSI Kaons",
                      savefilename="selMC_mass",
                      outputDir=outputDir)


def compare_histograms(config_file, hist_type="all"):
    hist_channels = {}
    # hist_channels[<Legend name>] = (<filename>, <hist suffix>)
    hist_channels["VSI"] = (config_file["mcFile"][0], "VSI", hist_type)
    hist_channels["VSI Kaons"] = (config_file["mcFile"][0], "VSI_2", hist_type)
    # hist_channels["VSI Leptons"] = (config_file["mcFile"][1], "VSI_Leptons", hist_type)

    plotting.compare(hist_channels,
                     variable='DV_mass',
                     setrange=(0, 20),
                     nRebin=100,
                     )

    plotting.compare(hist_channels,
                     variable='DV_mass',
                     vertical_lines=[.4977 - .01, .4977 + .01],
                     save_name="DV_mass_zoomed_Kshort",
                     setrange=(0.4, 0.6),
                     vertical_legend="K_{S}^{0} mass cut",
                     )

    plotting.compare(hist_channels,
                     variable='DV_mass',
                     vertical_lines=[3.096 - .01, 3.096 + .01],
                     save_name="DV_mass_zoomed_JPsi",
                     setrange=(2.8, 3.2),
                     vertical_legend="J/\Psi mass cut",
                     )

    hist_channels = {}
    hist_channels["prompt electrons"] = (config_file["mcFile"][0], "VSI_2", hist_type, 'prompt_electron')
    hist_channels["prompt muons"] = (config_file["mcFile"][0], "VSI_2", hist_type, 'prompt_muon')
    hist_channels["prompt leptons"] = (config_file["mcFile"][0], "VSI_2", hist_type, 'prompt_lepton')

    plotting.compare(hist_channels,
                     save_name="prompt_leptons",
                     scaleymax=2.2,
                     setrange=(0, 30),
                     norm=0,
                     notes=["VSI Kaons", "Dijet JZ4W"],
                     )

    hist_channels = {}
    # hist_channels[<Legend name>] = (<filename>, <hist suffix>)
    hist_channels["VSI"] = (config_file["mcFile"][0], "VSI", hist_type)
    hist_channels["VSI Kaons"] = (config_file["mcFile"][0], "VSI_2", hist_type)

    plotting.compare(hist_channels,
                     variable='DV_r',
                     nRebin=15,
                     setrange=(0, 350),
                     vertical_lines=MATERIAL_LAYERS,
                     vertical_legend="Material Layers",
                     )

    plotting.compare(hist_channels,
                     variable='DV_eta',
                     setrange=(-3, 3),
                     )

    plotting.compare(hist_channels,
                     variable='DV_trk_d0',
                     setrange=(-10, 10),
                     )

    plotting.compare(hist_channels,
                     variable='DV_trk_z0',
                     setrange=(-10, 10),
                     )

    plotting.compare(hist_channels,
                     variable='DV_trk_pt',
                     setrange=(0, 20),
                     )

    plotting.compare(hist_channels,
                     variable='DV_sum_track_pt',
                     setrange=(0, 20),
                     )

    plotting.compare(hist_channels,
                     variable='DV_sum_track_pt_wrt_pv',
                     setrange=(0, 20),
                     )

    plotting.compare(hist_channels,
                     variable='DV_sum_track_pt_diff',
                     setlogy=True,
                     setrange=(0, 1),
                     )

    plotting.compare(hist_channels,
                     variable='DV_sum_track_charge',
                     # setrange=(-2, 2),
                     )

    plotting.compare(hist_channels,
                     variable='DV_pt',
                     setrange=(0, 20),
                     )

    hist_channels = {}
    # hist_channels[<Legend name>] = (<filename>, <hist suffix>)
    # hist_channels["VSI"] = (config_file["mcFile"][0], "VSI", '')
    hist_channels["VSI Kaons"] = (config_file["mcFile"][0], "VSI_2", '')

    plotting.compare(hist_channels,
                     variable='muon_pt',
                     nRebin=5,
                     setrange=(0, 100),
                     )

    plotting.compare(hist_channels,
                     variable='el_pt',
                     nRebin=5,
                     setrange=(0, 100),
                     )

    plotting.compare(hist_channels,
                     variable='muon_quality',
                     setlogy=True,
                     bin_labels=['None', 'Loose', 'Medium', 'Tight']
                     )

    plotting.compare(hist_channels,
                     variable='el_quality',
                     setlogy=True,
                     bin_labels=['None', 'LHLoose', 'LHMedium', 'LHTight']
                     )

if __name__ == '__main__':
    import argparse


    class AppendActionCleanDefault(argparse._AppendAction):
        def __init__(self, *args, **kwargs):
            super(argparse._AppendAction, self).__init__(*args, **kwargs)
            self.index = 0

        def __call__(self, parser, namespace, values, option_string=None):
            items = argparse._copy.copy(argparse._ensure_value(namespace, self.dest, [])) if self.index else []
            if values:
                self.index += 1
                items.append(values)
                setattr(namespace, self.dest, items)


    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

    parser.add_argument("--config",
                        dest="config",
                        type=str,
                        required=True,
                        help="Input config file for plotHisotgrams.py.")

    parent_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[parser])

    options = parent_parser.parse_args()

    if os.path.exists(outputDir):
        pass
    else:
        logger.info('Making output directories...')
        os.mkdir(outputDir)
        os.mkdir(outputDir + "Cutflows/")
        os.mkdir(outputDir + "2Dmassplots")
        os.mkdir(outputDir + "2Dmassplots/data/")
        os.mkdir(outputDir + "2Dmassplots/mc/")

    with open(options.config, 'r') as json_config:
        config_file = json.load(json_config)  # load JSON config file

    plotting.plot_cutflow(file=config_file["mcFile"][0],
                          vertextype="VSI Kaons",
                          )
    plotting.plot_cutflow(file=config_file["mcFile"][0],
                          vertextype="VSI"
                          )

    # execute plotting here, comment out functions in you dont want to plot them again.
    compare_histograms(config_file, "all")
    # compare_histograms(config_file, "alpha")
    compare_histograms(config_file, "sel")
# compareMCdata(config_file, "mass")
# compareMCdata(config_file, "sel")
# makeCutflows(config_file)
# compareHistograms(config_file)
# make2Dmassplots(config_file)
