# Plotting Script
import argparse, os, math, ROOT, glob, uproot, time, json
import numpy as np
import helpers
from ROOT import *
from ROOT import gPad
from pylab import *
import plotting
import atlas_style

ROOT.gROOT.SetBatch(True)

#trying to  set ATLAS style a million different ways
gROOT.SetStyle("ATLAS") #might have to change how you set atlas style like this, depends how you have setup python
# atlas_style.AtlasStyle()	
# gROOT.LoadMacro("AtlasStyle.C")
# gROOT.LoadMacro("AtlasUtils.C")
# gROOT.LoadMacro("AtlasLabels.C")
# SetAtlasStyle()

logger = helpers.getLogger('dHNLAnalysis.plotHisotgrams')

#############################################################################################################################################
# globals
outputDir = '../output/plots/' # change path here to save your histograms somewhere else!


#############################################################################################################################################



def compareMC(config_file):

	# add plots here if you want to compare different distributions with with data & MC
	plotting.compare_MC(
							mcfiles=config_file["mcFile"], 
							hmclabels = config_file["mcLabel"],
							hname = "charge_DV_mass", 
							scaleymax = 1.5,
							setrange= "0 20",					
							vertextype = "VSI",
							outputDir=outputDir)
	plotting.compare_MC(
							mcfiles=config_file["mcFile"], 
							hmclabels = config_file["mcLabel"],
							hname = "charge_DV_mass", 
							scaleymax = 1.5,
							setrange= "0 20",					
							vertextype = "VSI Leptons",
							outputDir=outputDir)
	plotting.compare_MC(
							mcfiles=config_file["mcFile"], 
							hmclabels = config_file["mcLabel"],
							hname = "charge_DV_r", 
							setrange= "0 300",					
							nRebin = 4, 
							vertextype = "VSI",
							outputDir=outputDir)
	plotting.compare_MC(
							mcfiles=config_file["mcFile"], 
							hmclabels = config_file["mcLabel"],
							hname = "charge_DV_r", 
							setrange= "0 300",					
							nRebin = 4, 
							vertextype = "VSI Leptons",
							outputDir=outputDir)
	plotting.compare_MC(
							mcfiles=config_file["mcFile"], 
							hmclabels = config_file["mcLabel"],
							hname = "charge_HNLm", 
							scaleymax = 1.5,
							setrange= "0 20",					
							vertextype = "VSI",
							outputDir=outputDir)
	plotting.compare_MC(
							mcfiles=config_file["mcFile"], 
							hmclabels = config_file["mcLabel"],
							hname = "charge_HNLm", 
							scaleymax = 1.5,
							setrange= "0 20",					
							vertextype = "VSI Leptons",
							outputDir=outputDir)



	plotting.compare_MC(
							mcfiles=config_file["mcFile"], 
							hmclabels = config_file["mcLabel"],
							hname = "sel_DV_mass", 
							scaleymax = 1.5,
							setrange= "0 20",					
							vertextype = "VSI",
							outputDir=outputDir)
	plotting.compare_MC(
							mcfiles=config_file["mcFile"], 
							hmclabels = config_file["mcLabel"],
							hname = "sel_DV_mass", 
							scaleymax = 1.5,
							setrange= "0 20",					
							vertextype = "VSI Leptons",
							outputDir=outputDir)
	plotting.compare_MC(
							mcfiles=config_file["mcFile"], 
							hmclabels = config_file["mcLabel"],
							hname = "sel_DV_r", 
							setrange= "0 300",					
							nRebin = 4, 
							vertextype = "VSI",
							outputDir=outputDir)
	plotting.compare_MC(
							mcfiles=config_file["mcFile"], 
							hmclabels = config_file["mcLabel"],
							hname = "sel_DV_r", 
							setrange= "0 300",					
							nRebin = 4, 
							vertextype = "VSI Leptons",
							outputDir=outputDir)
	plotting.compare_MC(
							mcfiles=config_file["mcFile"], 
							hmclabels = config_file["mcLabel"],
							hname = "sel_HNLm", 
							scaleymax = 1.5,
							setrange= "0 20",					
							vertextype = "VSI",
							outputDir=outputDir)
	plotting.compare_MC(
							mcfiles=config_file["mcFile"], 
							hmclabels = config_file["mcLabel"],
							hname = "sel_HNLm", 
							scaleymax = 1.5,
							setrange= "0 20",					
							vertextype = "VSI Leptons",
							outputDir=outputDir)




def compareMCdata(config_file):

	# add plots here if you want to compare different distributions with with data & MC
	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel=  config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = "charge_DV_mass", 
							scaleymax = 1.5,
							setrange= "0 20",					
							vertextype = "VSI",
							outputDir=outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel=  config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],				
							hname = "charge_DV_mass", 
						    setrange= "0 20",
							vertextype = "VSI Leptons",
							outputDir=outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = "charge_DV_r", 
							nRebin = 5, 
							vertextype = "VSI Leptons",
							outputDir= outputDir)


	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = "charge_DV_r", 
							nRebin = 5, 
							vertextype = "VSI",
							outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_DV_trk_pt",
	# 						nRebin = 2, 
	# 						setrange= "0 100",					
	# 						vertextype = "VSI Leptons",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_DV_trk_pt",
	# 						nRebin = 2, 
	# 						setrange= "0 100",
	# 						vertextype = "VSI",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_DV_trk_phi", 
	# 						scaleymax = 2.2,
	# 						vertextype = "VSI Leptons",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_DV_trk_phi", 
	# 						scaleymax = 2.2,
	# 						vertextype = "VSI",
	# 						outputDir= outputDir)


	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_DV_trk_eta", 
	# 						scaleymax = 2.2,
	# 						vertextype = "VSI Leptons",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"],
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"], 
	# 						hname = "charge_DV_trk_eta", 
	# 						scaleymax = 2.2,
	# 						vertextype = "VSI",
	# 						outputDir= outputDir)


	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_DV_trk_d0",
	# 						setrange= "-10 10",
	# 						vertextype = "VSI Leptons",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_DV_trk_d0", 
	# 						setrange= "-10 10",							
	# 						vertextype = "VSI",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_DV_trk_dpt",
	# 						nRebin = 2, 
	# 						setrange= "0 20",
	# 						vertextype = "VSI Leptons",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_DV_trk_dpt",
	# 						nRebin = 2, 
	# 						setrange= "0 20",						
	# 						vertextype = "VSI",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_DV_trk_deta", 
	# 						nRebin = 2, 
	# 						setrange= "0 4",
	# 						vertextype = "VSI Leptons",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 					  	mcfiles=config_file["mcFile"],
	# 					  	hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"], 
	# 						hname = "charge_DV_trk_deta", 
	# 						nRebin = 2, 
	# 						setrange= "0 4",
	# 						vertextype = "VSI",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_DV_trk_dphi",
	# 						setrange= "0 8",
	# 						vertextype = "VSI Leptons",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_DV_trk_dphi",  
	# 						setrange= "0 8",	
	# 						vertextype = "VSI",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_DV_trk_dR", 
	# 						nRebin = 10, 
	# 						setrange= "0 2",							
	# 						vertextype = "VSI Leptons",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_DV_trk_dR", 
	# 						nRebin = 10, 
	# 						setrange= "0 10",
	# 						vertextype = "VSI",
	# 						outputDir= outputDir)


	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_mvis", 
	# 						nRebin = 5, 
	# 						setrange= "0 200",
	# 						scaleymax = 1.5,
	# 						vertextype = "VSI Leptons",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_mvis", 
	# 						nRebin = 5, 
	# 						setrange= "0 200",
	# 						scaleymax = 1.5,							
	# 						vertextype = "VSI",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_mtrans",
	# 						nRebin = 5, 
	# 						setrange= "0 200",
	# 						scaleymax = 1.5,
	# 						vertextype = "VSI Leptons",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_mtrans", 
	# 						nRebin = 5, 
	# 						setrange= "0 200",
	# 						scaleymax = 1.5,
	# 						vertextype = "VSI",
	# 						outputDir= outputDir)



	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"],
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_DV_redmass",  
	# 						setrange= "0 20",
	# 						scaleymax = 1.5,
	# 						vertextype = "VSI Leptons",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"],
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_DV_redmass", 
	# 						setrange= "0 20",
	# 						scaleymax = 1.5,
	# 						vertextype = "VSI",
	# 						outputDir= outputDir)


	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_DV_redmassvis",
	# 						nRebin = 5, 
	# 						setrange= "0 200",
	# 						scaleymax = 1.5,
	# 						vertextype = "VSI Leptons",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_DV_redmassvis", 
	# 						nRebin = 5, 
	# 						setrange= "0 200",
	# 						scaleymax = 1.5,
	# 						vertextype = "VSI",
	# 						outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = "charge_HNLm",
							setrange= "0 20",
							scaleymax = 1.5,							 							
							vertextype = "VSI Leptons",
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = "charge_HNLm", 
							setrange= "0 20",
							scaleymax = 1.5,
							vertextype = "VSI",
							outputDir= outputDir)


	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"],
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"], 
	# 						hname = "charge_DV_redmassHNL",
	# 						setrange= "0 50",
	# 						scaleymax = 1.5,
	# 						vertextype = "VSI Leptons",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_DV_redmassHNL", 
	# 						nRebin = 1, 
	# 						setrange= "0 50",
	# 						scaleymax = 1.5,				
	# 						vertextype = "VSI",
	# 						outputDir= outputDir)



	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_HNLpt",
	# 						nRebin = 5, 
	# 						setrange= "0 200",
	# 						scaleymax = 1.5,
	# 						vertextype = "VSI Leptons",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_HNLpt",
	# 						nRebin = 5, 
	# 						setrange= "0 200",
	# 						scaleymax = 1.5,
	# 						vertextype = "VSI",
	# 						outputDir= outputDir)


	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_HNLphi", 
	# 						setrange= "-4 4",
	# 						scaleymax = 1.5,
	# 						vertextype = "VSI Leptons",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_HNLphi", 
	# 						setrange= "-4 4",
	# 						scaleymax = 1.5,
	# 						vertextype = "VSI",
	# 						outputDir= outputDir)


	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_HNLeta",
	# 						setrange= "-3 3",
	# 						scaleymax = 1.5,
	# 						vertextype = "VSI Leptons",
	# 						outputDir= outputDir)

	# plotting.compare_dataMC(datafile=config_file["dataFile"],
	# 						mcfiles=config_file["mcFile"], 
	# 						hdatalabel= config_file["dataLabel"],
	# 						hmclabels = config_file["mcLabel"],
	# 						hname = "charge_HNLeta",
	# 						setrange= "-3 3",
	# 						scaleymax = 1.5,			
	# 						vertextype = "VSI",
	# 						outputDir= outputDir)




def make2Dmassplots(config_file):

	# MC FILES
	plotting.CorrPlot2D( file= config_file["mcFile"][0], 
				hname="charge_DVmass_mvis", 
				hlabel=config_file["mcLabel"][0],
				rebinx=2,
				rebiny=5,
				vertextype="VSI",
				setxrange="0 30",
				setyrange="0 300",
				outputDir=outputDir + "2Dmassplots/mc/")

	plotting.CorrPlot2D( file= config_file["mcFile"][0], 
				hlabel=config_file["mcLabel"][0],
				hname="charge_DVmass_mvis",
				rebinx=2,
				rebiny=5,
				vertextype="VSI Leptons",
				setxrange="0 30",
				setyrange="0 300",
				outputDir=outputDir + "2Dmassplots/mc/")

	plotting.CorrPlot2D(file= config_file["mcFile"][0], 
			rebinx=2,
			rebiny=5,
			hname="charge_DVmass_mtrans", 
			hlabel=config_file["mcLabel"][0],
			vertextype="VSI",
			setxrange="0 30",
			setyrange="0 300",
			outputDir=outputDir + "2Dmassplots/mc/")

	plotting.CorrPlot2D(file= config_file["mcFile"][0], 
			rebinx=2,
			rebiny=5,
			hname="charge_DVmass_mtrans", 
			hlabel=config_file["mcLabel"][0],
			vertextype="VSI Leptons",
			setxrange="0 30",
			setyrange="0 300",
			outputDir=outputDir + "2Dmassplots/mc/")

	plotting.CorrPlot2D(file= config_file["mcFile"][0], 
			rebinx=5,
			rebiny=5,
			hname="charge_mvis_mtrans", 
			hlabel=config_file["mcLabel"][0],
			vertextype="VSI",
			setxrange="0 100",
			setyrange="0 300",
			outputDir=outputDir + "2Dmassplots/mc/")

	plotting.CorrPlot2D(file= config_file["mcFile"][0], 
			rebinx=5,
			rebiny=5,
			hname="charge_mvis_mtrans", 
			hlabel=config_file["mcLabel"][0],
			vertextype="VSI Leptons",
			setxrange="0 100",
			setyrange="0 300",
			outputDir=outputDir + "2Dmassplots/mc/")

	plotting.CorrPlot2D(file= config_file["mcFile"][0], 
			hname="charge_DVmass_mhnl", 
			hlabel=config_file["mcLabel"][0],
			vertextype="VSI",
			setxrange="0 50",
			setyrange="0 50",
			outputDir=outputDir + "2Dmassplots/mc/")

	plotting.CorrPlot2D(file= config_file["mcFile"][0], 
			hname="charge_DVmass_mhnl", 
			hlabel=config_file["mcLabel"][0],
			vertextype="VSI Leptons",
			setxrange="0 50",
			setyrange="0 50",
			outputDir=outputDir + "2Dmassplots/mc/")


	plotting.CorrPlot2D(file= config_file["mcFile"][0], 
			rebinx=5,
			rebiny=1,
			hname="charge_mvis_mhnl", 
			hlabel=config_file["mcLabel"][0],
			vertextype="VSI",
			setxrange="0 300",
			setyrange="0 50",
			outputDir=outputDir + "2Dmassplots/mc/")

	plotting.CorrPlot2D(file= config_file["mcFile"][0], 
			rebinx=5,
			rebiny=1,
			hname="charge_mvis_mhnl", 
			hlabel=config_file["mcLabel"][0],
			vertextype="VSI Leptons",
			setxrange="0 300",
			setyrange="0 50",
			outputDir=outputDir + "2Dmassplots/mc/")


	plotting.CorrPlot2D(file= config_file["mcFile"][0], 
			rebinx=1,
			rebiny=5,
			hname="charge_mhnl_mtrans", 
			hlabel=config_file["mcLabel"][0],
			vertextype="VSI",
			setxrange="0 50",
			setyrange="0 300",
			outputDir=outputDir + "2Dmassplots/mc/")

	plotting.CorrPlot2D(file= config_file["mcFile"][0], 
			rebinx=1,
			rebiny=5,
			hname="charge_mhnl_mtrans", 
			hlabel=config_file["mcLabel"][0],
			vertextype="VSI Leptons",
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
			outputDir=outputDir +"2Dmassplots/data/")

	plotting.CorrPlot2D(config_file["dataFile"], 
			hname="charge_DVmass_mvis", 
			hlabel=config_file["dataLabel"],
			rebinx=2,
			rebiny=5,
			vertextype="VSI Leptons",
			setxrange="0 30",
			setyrange="0 300",
			outputDir=outputDir +"2Dmassplots/data/")

	plotting.CorrPlot2D(config_file["dataFile"], 
			rebinx=2,
			rebiny=5,
			hname="charge_DVmass_mtrans", 
			hlabel=config_file["dataLabel"],
			vertextype="VSI",
			setxrange="0 30",
			setyrange="0 300",
			outputDir=outputDir +"2Dmassplots/data/")

	plotting.CorrPlot2D(config_file["dataFile"], 
			rebinx=2,
			rebiny=5,
			hname="charge_DVmass_mtrans", 
			hlabel=config_file["dataLabel"],
			vertextype="VSI Leptons",
			setxrange="0 30",
			setyrange="0 300",
			outputDir=outputDir +"2Dmassplots/data/")

	plotting.CorrPlot2D(config_file["dataFile"], 
			rebinx=5,
			rebiny=5,
			hname="charge_mvis_mtrans", 
			hlabel=config_file["dataLabel"],
			vertextype="VSI",
			setxrange="0 100",
			setyrange="0 300",
			outputDir=outputDir +"2Dmassplots/data/")

	plotting.CorrPlot2D(config_file["dataFile"], 
			rebinx=5,
			rebiny=5,
			hname="charge_mvis_mtrans", 
			hlabel=config_file["dataLabel"],
			vertextype="VSI Leptons",
			setxrange="0 100",
			setyrange="0 300",
			outputDir=outputDir +"2Dmassplots/data/")

	plotting.CorrPlot2D(config_file["dataFile"], 
			hname="charge_DVmass_mhnl", 
			hlabel=config_file["dataLabel"],
			vertextype="VSI",
			setxrange="0 50",
			setyrange="0 50",
			outputDir=outputDir +"2Dmassplots/data/")

	plotting.CorrPlot2D(config_file["dataFile"], 
			hname="charge_DVmass_mhnl", 
			hlabel=config_file["dataLabel"],
			vertextype="VSI Leptons",
			setxrange="0 50",
			setyrange="0 50",
			outputDir=outputDir +"2Dmassplots/data/")


	plotting.CorrPlot2D(config_file["dataFile"], 
			rebinx=5,
			rebiny=1,
			hname="charge_mvis_mhnl", 
			hlabel=config_file["dataLabel"],
			vertextype="VSI",
			setxrange="0 300",
			setyrange="0 50",
			outputDir=outputDir +"2Dmassplots/data/")

	plotting.CorrPlot2D(config_file["dataFile"], 
			rebinx=5,
			rebiny=1,
			hname="charge_mvis_mhnl", 
			hlabel=config_file["dataLabel"],
			vertextype="VSI Leptons",
			setxrange="0 300",
			setyrange="0 50",
			outputDir=outputDir +"2Dmassplots/data/")


	plotting.CorrPlot2D(config_file["dataFile"], 
			rebinx=1,
			rebiny=5,
			hname="charge_mhnl_mtrans", 
			hlabel=config_file["dataLabel"],
			vertextype="VSI",
			setxrange="0 50",
			setyrange="0 300",
			outputDir=outputDir +"2Dmassplots/data/")

	plotting.CorrPlot2D(config_file["dataFile"], 
			rebinx=1,
			rebiny=5,
			hname="charge_mhnl_mtrans", 
			hlabel=config_file["dataLabel"],
			vertextype="VSI Leptons",
			setxrange="0 50",
			setyrange="0 300",
			outputDir=outputDir +"2Dmassplots/data/")


def makeCutflows(config_file):
	# make Cutflow plots here
	plotting.plot_cutflow(file = config_file["dataFile"],
						  vertextype= "VSI",
						  outputDir=outputDir + "Cutflows/")
	plotting.plot_cutflow(file = config_file["dataFile"],
						  vertextype= "VSI Leptons",
						  outputDir=outputDir + "Cutflows/")
	plotting.plot_cutflow(file = config_file["mcFile"][0],
						  vertextype= "VSI",
						  outputDir=outputDir + "Cutflows/")
	plotting.plot_cutflow(file = config_file["mcFile"][0],
						  vertextype= "VSI Leptons",
						  outputDir=outputDir + "Cutflows/")
	plotting.plot_cutflow(file = config_file["mcFile"][1],
						  vertextype= "VSI",
						  outputDir=outputDir + "Cutflows/")
	plotting.plot_cutflow(file = config_file["mcFile"][1],
						  vertextype= "VSI Leptons",
						  outputDir=outputDir + "Cutflows/")
	plotting.plot_cutflow(file = config_file["mcFile"][2],
						  vertextype= "VSI",
						  outputDir=outputDir + "Cutflows/")
	plotting.plot_cutflow(file = config_file["mcFile"][2],
						  vertextype= "VSI Leptons",
						  outputDir=outputDir + "Cutflows/")

def compareHistograms(config_file):
	# compare different histograms from the same file here
	plotting.compareN(file=config_file["mcFile"][0],
			hname = ["sel_DV_mass","sel_HNLm"],
			hlabel = ["DV mass","m_{HNL}"], 
			setxrange= "0 20",
			scaleymax=1.5,
			nRebin=1,
			vertextype= "VSI",
			savefilename="selMC_mass_VSI",
			outputDir=outputDir)
	plotting.compareN(file=config_file["mcFile"][0],
			hname = ["sel_DV_mass","sel_HNLm"],
			hlabel = ["DV mass","m_{HNL}"], 
			setxrange= "0 20",
			scaleymax=1.5,
			nRebin=1,
			vertextype= "VSI Leptons",
			savefilename="selMC_mass_VSI_Leptons",
			outputDir=outputDir)
	plotting.compareN(file=config_file["mcFile"][1],
			hname = ["sel_DV_mass","sel_HNLm"],
			hlabel = ["DV mass","m_{HNL}"], 
			setxrange= "0 20",
			scaleymax=1.5,
			nRebin=1,
			vertextype= "VSI",
			savefilename="selMC_mass_VSI",
			outputDir=outputDir)
	plotting.compareN(file=config_file["mcFile"][1],
			hname = ["sel_DV_mass","sel_HNLm"],
			hlabel = ["DV mass","m_{HNL}"], 
			setxrange= "0 20",
			scaleymax=1.5,
			nRebin=1,
			vertextype= "VSI Leptons",
			savefilename="selMC_mass_VSI_Leptons",
			outputDir=outputDir)
	plotting.compareN(file=config_file["mcFile"][2],
			hname = ["sel_DV_mass","sel_HNLm"],
			hlabel = ["DV mass","m_{HNL}"], 
			setxrange= "0 20",
			scaleymax=1.5,
			nRebin=1,
			vertextype= "VSI",
			savefilename="selMC_mass_VSI",
			outputDir=outputDir)
	plotting.compareN(file=config_file["mcFile"][2],
			hname = ["sel_DV_mass","sel_HNLm"],
			hlabel = ["DV mass","m_{HNL}"], 
			setxrange= "0 20",
			scaleymax=1.5,
			nRebin=1,
			vertextype= "VSI Leptons",
			savefilename="selMC_mass_VSI_Leptons",
			outputDir=outputDir)
	# plotting.compareN(file=config_file["mcFile"][0],
	# 		hname = ["sel_DV_mass","sel_mvis","sel_HNLm","sel_mtrans","sel_DV_redmass"],
	# 		hlabel = ["DV mass","Visible mass","m_{HNL}", "m_{T}","Reduced mass"], 
	# 		setxrange= "0 100",
	# 		scaleymax=1.5,
	# 		nRebin=5,
	# 		vertextype= "VSI",
	# 		savefilename="selMC_mass_VSI",
	# 		outputDir=outputDir)
	# plotting.compareN(file=config_file["mcFile"][0],
	# 		hname = ["sel_DV_mass","sel_mvis","sel_HNLm","sel_mtrans","sel_DV_redmass"],
	# 		hlabel = ["DV mass","Visible mass","m_{HNL}", "m_{T}","Reduced mass"], 
	# 		setxrange= "0 100",
	# 		scaleymax=1.5,
	# 		nRebin=5,
	# 		vertextype= "VSI Leptons",
	# 		savefilename="selMC_mass_VSI_Leptons",
	# 		outputDir=outputDir)

	plotting.compareN(file=config_file["mcFile"][0],
			hname = ["charge_DV_mass","charge_HNLm"],
			hlabel = ["DV mass","m_{HNL}"], 
			setxrange= "0 20",
			scaleymax=1.5,
			nRebin=1,
			vertextype= "VSI",
			savefilename="chargeMC_mass_VSI",
			outputDir=outputDir)
	plotting.compareN(file=config_file["mcFile"][0],
			hname = ["charge_DV_mass","charge_HNLm"],
			hlabel = ["DV mass","m_{HNL}"], 
			setxrange= "0 20",
			scaleymax=1.5,
			nRebin=1,
			vertextype= "VSI Leptons",
			savefilename="chargeMC_mass_VSI_Leptons",
			outputDir=outputDir)
	plotting.compareN(file=config_file["mcFile"][1],
			hname = ["charge_DV_mass","charge_HNLm"],
			hlabel = ["DV mass","m_{HNL}"], 
			setxrange= "0 20",
			scaleymax=1.5,
			nRebin=1,
			vertextype= "VSI",
			savefilename="chargeMC_mass_VSI",
			outputDir=outputDir)
	plotting.compareN(file=config_file["mcFile"][1],
			hname = ["charge_DV_mass","charge_HNLm"],
			hlabel = ["DV mass","m_{HNL}"], 
			setxrange= "0 20",
			scaleymax=1.5,
			nRebin=1,
			vertextype= "VSI Leptons",
			savefilename="chargeMC_mass_VSI_Leptons",
			outputDir=outputDir)
	plotting.compareN(file=config_file["mcFile"][2],
			hname = ["charge_DV_mass","charge_HNLm"],
			hlabel = ["DV mass","m_{HNL}"], 
			setxrange= "0 20",
			scaleymax=1.5,
			nRebin=1,
			vertextype= "VSI",
			savefilename="chargeMC_mass_VSI",
			outputDir=outputDir)
	plotting.compareN(file=config_file["mcFile"][2],
			hname = ["charge_DV_mass","charge_HNLm"],
			hlabel = ["DV mass","m_{HNL}"], 
			setxrange= "0 20",
			scaleymax=1.5,
			nRebin=1,
			vertextype= "VSI Leptons",
			savefilename="chargeMC_mass_VSI_Leptons",
			outputDir=outputDir)
	# plotting.compareN(file=config_file["mcFile"][0],
	# 		hname = ["charge_DV_mass","charge_mvis","charge_HNLm","charge_mtrans","charge_DV_redmass"],
	# 		hlabel = ["DV mass","Visible mass","m_{HNL}", "m_{T}","Reduced mass"], 
	# 		setxrange= "0 100",
	# 		scaleymax=1.5,
	# 		nRebin=5,
	# 		vertextype= "VSI",
	# 		savefilename="chargeMC_mass_VSI",
	# 		outputDir=outputDir)
	# plotting.compareN(file=config_file["mcFile"][0],
	# 		hname = ["charge_DV_mass","charge_mvis","charge_HNLm","charge_mtrans","charge_DV_redmass"],
	# 		hlabel = ["DV mass","Visible mass","m_{HNL}", "m_{T}","Reduced mass"], 
	# 		setxrange= "0 100",
	# 		scaleymax=1.5,
	# 		nRebin=5,
	# 		vertextype= "VSI Leptons",
	# 		savefilename="chargeMC_mass_VSI_Leptons",
	# 		outputDir=outputDir)





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

	parser.add_argument("--config",
						dest="config",
						type = str,
						required = True,
						help="Input config file for plotHisotgrams.py.")

	parent_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter, parents = [parser]) 

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
		config_file = json.load(json_config) # load JSON config file


	#execute plotting here, comment out functions in you dont want to plot them again.	
	compareMCdata(config_file)
	makeCutflows(config_file)
	# compareHistograms(config_file)
	# make2Dmassplots(config_file)
	compareMC(config_file)

	
	
