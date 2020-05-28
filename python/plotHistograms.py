# Plotting Script
import ROOT
import json
import helpers
import os
from ROOT import gROOT
from pylab import *
import plotting

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
MATERIAL_LAYERS = [33.25, 50.5, 88.5, 122.5, 299]


#############################################################################################################################################




def compareMCdata(config_file):
	# configs for all plots!
	histcut = "DVtype"
	setlogy = True
	normalize = False
	lumi = 60

	#add plots here if you want to compare different distributions with with data & MC
	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel=  config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_mass", 
							setlogy = setlogy,
							scaleymax = 100,
							setrange= "0 20",					
							vertextype = "VSI",
							normalize = normalize,
							lumi = lumi,
							outputDir=outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel=  config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],				
							hname = histcut + "_DV_mass",
							scaleymax = 2.2, 
						    setlogy = setlogy,
						    setrange= "0 20",
							vertextype = "VSI Leptons",
							normalize = normalize,
							lumi = lumi,
							outputDir=outputDir)
	
	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_r",
							setlogy = setlogy,
							nRebin = 5, 
							setrange= "0 300", 						
							vertextype = "VSI Leptons",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)


	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_r", 
							setlogy = setlogy,
							nRebin = 5, 
							setrange= "0 300", 
							vertextype = "VSI",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_trk_pt",
							setlogy = setlogy,
							nRebin = 2, 
							setrange= "0 100",					
							vertextype = "VSI Leptons",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_trk_pt",
							setlogy = setlogy,
							nRebin = 2, 
							setrange= "0 100",
							vertextype = "VSI",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_trk_phi", 
							setlogy = setlogy,
							scaleymax = 2.2,
							vertextype = "VSI Leptons",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_trk_phi", 
							setlogy = setlogy,
							scaleymax = 2.2,
							vertextype = "VSI",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)


	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_trk_eta", 
							setlogy = setlogy,
							scaleymax = 2.2*100,
							vertextype = "VSI Leptons",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"],
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"], 
							hname = histcut + "_DV_trk_eta", 
							setlogy = setlogy,
							scaleymax = 2.2*100,
							vertextype = "VSI",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)


	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_trk_d0",
							setlogy = setlogy,
							scaleymax = 1.5,
							setrange= "-100 100",
							vertextype = "VSI Leptons",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_trk_d0", 
							setlogy = setlogy,
							scaleymax = 1.5,
							setrange= "-100 100",							
							vertextype = "VSI",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_trk_dpt",
							setlogy = setlogy,
							nRebin = 2, 
							setrange= "0 20",
							vertextype = "VSI Leptons",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_trk_dpt",
							setlogy = setlogy,
							nRebin = 2, 
							setrange= "0 20",						
							vertextype = "VSI",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_trk_deta", 
							setlogy = setlogy,
							nRebin = 2, 
							setrange= "0 4",
							vertextype = "VSI Leptons",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
						  	mcfiles=config_file["mcFile"],
						  	hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"], 
							hname = histcut + "_DV_trk_deta", 
							setlogy = setlogy,
							nRebin = 2, 
							setrange= "0 4",
							vertextype = "VSI",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_trk_dphi",
							setlogy = setlogy,
							setrange= "0 3.2",
							vertextype = "VSI Leptons",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_trk_dphi",  
							setlogy = setlogy,
							setrange= "0 3.2",	
							vertextype = "VSI",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_trk_dR", 
							setlogy = setlogy,
							nRebin = 10, 
							setrange= "0 5",							
							vertextype = "VSI Leptons",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_trk_dR", 
							setlogy = setlogy,
							nRebin = 10, 
							setrange= "0 10",
							vertextype = "VSI",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)


	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_mvis", 
							setlogy = setlogy,
							nRebin = 5, 
							setrange= "0 200",
							scaleymax = 1.5,
							vertextype = "VSI Leptons",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_mvis", 
							setlogy = setlogy,
							nRebin = 5, 
							setrange= "0 200",
							scaleymax = 1.5,							
							vertextype = "VSI",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_mtrans",
							setlogy = setlogy,
							nRebin = 5, 
							setrange= "0 200",
							scaleymax = 1.5,
							vertextype = "VSI Leptons",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_mtrans", 
							setlogy = setlogy,
							nRebin = 5, 
							setrange= "0 200",
							scaleymax = 1.5,
							vertextype = "VSI",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)



	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"],
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_redmass",  
							setlogy = setlogy,					
							setrange= "0 50",
							scaleymax = 1.5,
							vertextype = "VSI Leptons",
							normalize = normalize,
							lumi = lumi,
	 						outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"],
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_redmass", 
							setlogy = setlogy,						
							setrange= "0 50",
							scaleymax = 1.5,
							vertextype = "VSI",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)


	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_redmassvis",
							setlogy = setlogy,						
							nRebin = 5, 
							setrange= "0 200",
							scaleymax = 1.5,
							vertextype = "VSI Leptons",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_redmassvis", 
							setlogy = setlogy,
							nRebin = 5, 
							setrange= "0 200",
							scaleymax = 1.5,
							vertextype = "VSI",
							normalize = normalize,
							lumi = lumi,
 							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_HNLm",
							setlogy = setlogy,
							setrange= "0 30",
							scaleymax = 1.2,							 							
							vertextype = "VSI Leptons",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_HNLm", 
							setlogy = setlogy,
							setrange= "0 30",
							scaleymax = 1.2,
							vertextype = "VSI",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)


	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"],
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"], 
							hname = histcut + "_DV_redmassHNL",
							setlogy = setlogy,						
							setrange= "0 50",
							scaleymax = 1.5,
							vertextype = "VSI Leptons",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)
						

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_DV_redmassHNL", 
							setlogy = setlogy,						
							nRebin = 1, 
							setrange= "0 50",
							scaleymax = 1.5,				
							vertextype = "VSI",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)
	

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_HNLpt",
							setlogy = setlogy,
							nRebin = 5, 
							setrange= "0 200",
							scaleymax = 1.5,
							vertextype = "VSI Leptons",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_HNLpt",
							setlogy = setlogy,
							nRebin = 5, 
							setrange= "0 200",
							scaleymax = 1.5,
							vertextype = "VSI",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)


	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_HNLphi", 
							setlogy = setlogy,
							setrange= "-4 4",
							scaleymax = 1.5,
							vertextype = "VSI Leptons",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_HNLphi", 
							setlogy = setlogy,
							setrange= "-4 4",
							scaleymax = 1.5,
							vertextype = "VSI",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)


	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_HNLeta",
							setlogy = setlogy,
							setrange= "-3 3",
							scaleymax = 1.5,
							vertextype = "VSI Leptons",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)

	plotting.compare_dataMC(datafile=config_file["dataFile"],
							mcfiles=config_file["mcFile"], 
							hdatalabel= config_file["dataLabel"],
							hmclabels = config_file["mcLabel"],
							hname = histcut + "_HNLeta",
							setlogy = setlogy,
							setrange= "-3 3",
							scaleymax = 1.5,			
							vertextype = "VSI",
							normalize = normalize,
							lumi = lumi,
							outputDir= outputDir)





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
	plotting.plot_cutflow(file = config_file["mcFile"][0],
						  vertextype= "VSI Leptons",
						  outputDir=outputDir + "Cutflows/")

def compareHistograms(config_file):
	# compare different histograms from the same file here
	plotting.compareN(file=config_file["mcFile"][0],
			hname = ["sel_DV_mass","sel_mvis","sel_HNLm","sel_mtrans","sel_DV_redmass"],
			hlabel = ["DV mass","Visible mass","m_{HNL}", "m_{T}","Reduced mass"], 
			setxrange= "0 100",
			scaleymax=1.5,
			nRebin=5,
			vertextype= "VSI Leptons",
			savefilename="selMC_mass",
			outputDir=outputDir)

def compare_histograms(config_file, selection):
	# Example: comparing data and MC for VSI and VSI_Leptons
	hist_channels = {}
	# hist_channels[<Legend name>] = (<filename>, <vertex directory>, <selection directory>)
	hist_channels["MC VSI"] = (config_file["mcFile"][0], "VSI", selection)
	hist_channels["Data VSI"] = (config_file["dataFile"], "VSI", selection)
	hist_channels["MC VSI Leptons"] = (config_file["mcFile"][0], "VSI_Leptons", selection)
	hist_channels["Data VSI Leptons"] = (config_file["dataFile"], "VSI_Leptons", selection)

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
	# compareMCdata(config_file)
	# makeCutflows(config_file)
	# compareHistograms(config_file)
	# make2Dmassplots(config_file)
	compare_histograms(config_file, 'all')

	
	