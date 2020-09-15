# Plotting Script
import os, math, ROOT, json,sys
import numpy as np
from ROOT import *
from ROOT import gPad
from pylab import *
sys.path.append('../python/')
# Considering your module contains a function called my_func, you could import it:
import plotting
import helpers
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
outputDir = '../output/MCreweighting/' # change path here to save your histograms somewhere else!
MATERIAL_LAYERS = [33.25, 50.5, 88.5, 122.5, 299]
normalize = True
setlogy = False
drawRatio=False
draw_channel_info = True
do_cut_significane = True

#############################################################################################################################################


def makeCutflows(config_file):
	# vtx_channels = ["VSI", "VSI_Leptons"]
	vtx_channels = ["VSI"]
	for vtx_channel in vtx_channels:
		plotting.plot_cutflow(file = config_file["dataFile"],
							  vertextype= vtx_channel,
							  output_dir=outputDir)
		nMCfiles = len(config_file["mcFiles"])

		for i in range(nMCfiles): 
			plotting.plot_cutflow(file = config_file["mcFiles"][i],
								  vertextype= vtx_channel,
								  output_dir=outputDir)


def check_rerunningVSI(config_file, selection):
	hist_channels = []
	# hist_channels[i] = (<filename>, <legend label>,<vertex directory>, <selection directory>)
	hist_channels.append([config_file["rerunVSIfile"],"VSI", "VSI", selection])
	hist_channels.append([config_file["rerunVSIfile"],"VSI_2", "VSI_2", selection])

	plotting.compare(hist_channels,
						 variable='DV_r',
						 nRebin=10,
						 setrange=(0, 350),
						 setlogy = setlogy,
						 scaleymax=2.5,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 normalize = True,
						 vertical_lines=MATERIAL_LAYERS,
						 vertical_legend="Material Layers",
						 output_dir= outputDir
						 )



def compare_histograms(config_file, selection):
	# vtx_channels = ["VSI", "VSI_Leptons"]
	vtx_channels = ["VSI"]
	for vtx_channel in vtx_channels:
		hist_channels = []
		# hist_channels[i] = (<filename>, <legend label>,<vertex directory>, <selection directory>)
		# hist_channels.append([config_file["dataFile"],config_file["dataLabel"], vtx_channel, selection])
		hist_channels.append([config_file["mcFiles"][0],config_file["mcLabels"][0], vtx_channel, selection])
		hist_channels.append([config_file["mcFiles"][1], config_file["mcLabels"][1], vtx_channel, selection])
		
		#get integrated luminosity to scale MC files to (ideally this should come from a value in the nutple TD DO) - DT
		scalelumi = config_file["scaleLumi"] # luminosity you want to scale everything to 
		datalumi = config_file["dataLumi"] #  lumi of the data you are looking at
		
		if "ratioLabel" in config_file.keys():
			ratioLabel = config_file["ratioLabel"]
		else:
			ratioLabel = [""]


		# DV Variables
		plotting.compare(hist_channels,
						 variable='DV_r',
						 nRebin=10,
						 setrange=(0, 350),
						 # setrange=(0, 10),
						 setlogy = setlogy,
						 scaleymax=1.6,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 vertical_lines=MATERIAL_LAYERS,
						 vertical_legend="Material Layers",
						 output_dir= outputDir,
						 use_ntuple = False,
						 ntup_nbins=350,
						 )
		# DV Track Variables 
		plotting.compare(hist_channels,
						 variable='DV_trk_pt',
						 setrange=(0, 100),
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir = outputDir
						 )
		plotting.compare(hist_channels,
					 variable='DV_trk_0_pt',
					 setrange=(0, 100),
					 setlogy = setlogy,
					 scalelumi = scalelumi,
					 datalumi = datalumi,
					 drawRatio = drawRatio,
					 ratioLabel = ratioLabel,
					 normalize = normalize,
					 draw_channel_info= draw_channel_info,
					 nRebin=2,
					 output_dir = outputDir
					 )

		plotting.compare(hist_channels,
				 variable='DV_trk_1_pt',
				 setrange=(0, 100),
				 setlogy = setlogy,
				 scalelumi = scalelumi,
				 datalumi = datalumi,
				 drawRatio = drawRatio,
				 ratioLabel = ratioLabel,
				 normalize = normalize,
				 draw_channel_info= draw_channel_info,
				 nRebin=2,
				 output_dir = outputDir
				 )


		plotting.compare(hist_channels,
						 variable='DV_trk_eta',
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_phi',
						 nRebin = 2, 
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_eta',
						 setrange=(-3, 3),
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_d0',
						 setrange=(-10, 10),
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_z0',
						 setrange=(-500, 500),
						 setlogy = setlogy,
						 scaleymax=1.6,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 nRebin = 4,
						 draw_channel_info= draw_channel_info,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_absz0',
						 setrange=(0, 250),
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_dTheta',
						 setrange=(0, 3),
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_nSCTHoles',
						 setrange=(0, 3),
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_nSiHits',
						 setrange=(0, 23),
						 setlogy = setlogy,
						 scaleymax=1.6,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_nSCTHits',
						 setrange=(0, 14),
						 setlogy = setlogy,
						 scaleymax=1.9,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_isLRT',
						 setrange=(0, 2),
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_chi2',
						 setrange=(0, 10),
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 # nRebin =  4,
						 draw_channel_info= draw_channel_info,
						 use_ntuple = False,
						 ntup_nbins=200,
						 output_dir = outputDir
						 )


		plotting.compare(hist_channels,
						 variable='DV_trk_dpt',
						 setrange=(0, 20),
						 nRebin = 2, 
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_deta',
						 setrange=(0, 3.2),
						 nRebin = 2, 
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_dphi',
						 setrange=(0, 3.2),
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_dR',
						 setrange=(0, 10),
						 nRebin = 10, 
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir = outputDir
						 )

		# Mass Variables 
		plotting.compare(hist_channels,
						 variable='DV_mass',
						 setrange=(0, 30),
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 # vertical_lines=[50],
						 nRebin=2,
						 # vertical_legend="DV mass cut",
						 output_dir= outputDir
						 )

		plotting.compare(hist_channels,
						 variable='mvis',
						 nRebin=5,
						 setrange=(0, 200),
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir= outputDir
						 )

		plotting.compare(hist_channels,
						 variable='mtrans',
						 nRebin=5,
						 setrange=(0, 200),
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir= outputDir
						 )
		plotting.compare(hist_channels,
						 variable='HNLm',
						 setrange=(0, 30),
						 setlogy = setlogy,
						 scaleymax=1.4,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 nRebin =2,
						 draw_channel_info= draw_channel_info,
						 do_cut_significane= do_cut_significane,
						 output_dir= outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_redmass',
						 setrange=(0, 50),
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir= outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_redmassvis',
						 setrange=(0, 200),
						 nRebin=5,
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir= outputDir
						 )
		plotting.compare(hist_channels,
						 variable='DV_redmassHNL',
						 setrange=(0, 50),
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir= outputDir
						 )

		# HNL quantities
		plotting.compare(hist_channels,
						 variable='HNLpt',
						 setrange=(0, 200),
						 nRebin=5,
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir= outputDir
						 )

		plotting.compare(hist_channels,
						 variable='HNLeta',
						 setrange=(0, 50),
						 nRebin=5,
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir= outputDir
						 )

		plotting.compare(hist_channels,
						 variable='HNLphi',
						 setrange=(-3,3),
						 setlogy = setlogy,
						 scalelumi = scalelumi,
						 datalumi = datalumi,
						 drawRatio = drawRatio,
						 ratioLabel = ratioLabel,
						 normalize = normalize,
						 draw_channel_info= draw_channel_info,
						 output_dir= outputDir
						 )

		plotting.compare(hist_channels,
					 variable='HNLphi',
					 setrange=(-3,3),
					 setlogy = setlogy,
					 scalelumi = scalelumi,
					 datalumi = datalumi,
					 drawRatio = drawRatio,
					 ratioLabel = ratioLabel,
					 normalize = normalize,
					 draw_channel_info= draw_channel_info,
					 output_dir= outputDir
					 )

def compare_TRUTH_histograms(config_file, selection):
	# vtx_channels = ["VSI", "VSI_Leptons"]
	vtx_channels = ["VSI"]
	for vtx_channel in vtx_channels:
		hist_channels = []
		# hist_channels[i] = (<filename>, <legend label>,<vertex directory>, <selection directory>)
		# hist_channels.append([config_file["dataFile"],config_file["dataLabel"], vtx_channel, selection])
		hist_channels.append([config_file["mcFiles"][0],config_file["mcLabels"][0], vtx_channel, selection])
		hist_channels.append([config_file["mcFiles"][1], config_file["mcLabels"][1], vtx_channel, selection])
		
		#get integrated luminosity to scale MC files to (ideally this should come from a value in the nutple TD DO) - DT
		scalelumi = config_file["scaleLumi"] # luminosity you want to scale everything to 
		datalumi = config_file["dataLumi"] #  lumi of the data you are looking at
		
		if "ratioLabel" in config_file.keys():
			ratioLabel = config_file["ratioLabel"]
		else:
			ratioLabel = [""]



		plotting.compare(hist_channels,
					 variable='all_DV_mass',
					 setrange=(0,20),
					 setlogy = setlogy,
					 scalelumi = scalelumi,
					 datalumi = datalumi,
					 drawRatio = drawRatio,
					 ratioLabel = ratioLabel,
					 normalize = normalize,
					 draw_channel_info= draw_channel_info,
					 output_dir= outputDir
					 )
		plotting.compare(hist_channels,
				 variable='all_plep_pt',
				 #setrange=(-3,3),
				 setlogy = setlogy,
				 scalelumi = scalelumi,
				 datalumi = datalumi,
				 drawRatio = drawRatio,
				 ratioLabel = ratioLabel,
				 normalize = normalize,
				 # nRebin =2,
				 draw_channel_info= draw_channel_info,
				 output_dir= outputDir
				 )


		plotting.compare(hist_channels,
				 variable='all_lep1_trk_pt',
				 #setrange=(-3,3),
				 setlogy = setlogy,
				 scalelumi = scalelumi,
				 datalumi = datalumi,
				 drawRatio = drawRatio,
				 ratioLabel = ratioLabel,
				 normalize = normalize,
				 # nRebin =2,
				 draw_channel_info= draw_channel_info,
				 output_dir= outputDir
				 )

		plotting.compare(hist_channels,
				 variable='all_lep2_trk_pt',
				 #setrange=(-3,3),
				 setlogy = setlogy,
				 scalelumi = scalelumi,
				 datalumi = datalumi,
				 drawRatio = drawRatio,
				 ratioLabel = ratioLabel,
				 normalize = normalize,
				 # nRebin =2,
				 draw_channel_info= draw_channel_info,
				 output_dir= outputDir
				 )

		plotting.compare(hist_channels,
			 variable='all_m12',
			 setrange=(0,20),
			 setlogy = setlogy,
			 scalelumi = scalelumi,
			 datalumi = datalumi,
			 drawRatio = drawRatio,
			 ratioLabel = ratioLabel,
			 normalize = normalize,
			 draw_channel_info= draw_channel_info,
			 output_dir= outputDir
			 )
		plotting.compare(hist_channels,
			 variable='all_m12_sq',
			 setrange=(0,100),
			 setlogy = setlogy,
			 scalelumi = scalelumi,
			 datalumi = datalumi,
			 drawRatio = drawRatio,
			 ratioLabel = ratioLabel,
			 normalize = normalize,
			 draw_channel_info= draw_channel_info,
			 nRebin= 5,
			 output_dir= outputDir
			 )

		plotting.compare(hist_channels,
			 variable='all_m13',
			 setrange=(0,20),
			 setlogy = setlogy,
			 scalelumi = scalelumi,
			 datalumi = datalumi,
			 drawRatio = drawRatio,
			 ratioLabel = ratioLabel,
			 normalize = normalize,
			 draw_channel_info= draw_channel_info,
			 output_dir= outputDir
			 )
		plotting.compare(hist_channels,
			 variable='all_m13_sq',
			 setrange=(0,100),
			 setlogy = setlogy,
			 scalelumi = scalelumi,
			 datalumi = datalumi,
			 drawRatio = drawRatio,
			 ratioLabel = ratioLabel,
			 normalize = normalize,
			 draw_channel_info= draw_channel_info,
			 nRebin= 5,
			 output_dir= outputDir
			 )
		plotting.compare(hist_channels,
			 variable='all_m23',
			 setrange=(0,20),
			 setlogy = setlogy,
			 scalelumi = scalelumi,
			 datalumi = datalumi,
			 drawRatio = drawRatio,
			 ratioLabel = ratioLabel,
			 normalize = normalize,
			 draw_channel_info= draw_channel_info,
			 output_dir= outputDir
			 )
		plotting.compare(hist_channels,
			 variable='all_m23_sq',
			 setrange=(0,100),
			 setlogy = setlogy,
			 scalelumi = scalelumi,
			 datalumi = datalumi,
			 drawRatio = drawRatio,
			 ratioLabel = ratioLabel,
			 normalize = normalize,
			 draw_channel_info= draw_channel_info,
			 nRebin= 5,
			 output_dir= outputDir
			 )

		plotting.compare(hist_channels,
			 variable='all_nu_trk_pt',
			 # setrange=(0,20),
			 setlogy = setlogy,
			 scalelumi = scalelumi,
			 datalumi = datalumi,
			 drawRatio = drawRatio,
			 ratioLabel = ratioLabel,
			 normalize = normalize,
			 # nRebin =2,
			 draw_channel_info= draw_channel_info,
			 output_dir= outputDir
			 )


		plotting.compare(hist_channels,
			 variable='all_s12',
			 setrange=(0,5999),
			 setlogy = setlogy,
			 scalelumi = scalelumi,
			 datalumi = datalumi,
			 drawRatio = drawRatio,
			 ratioLabel = ratioLabel,
			 normalize = normalize,
			 nRebin =100,
			 draw_channel_info= draw_channel_info,
			 output_dir= outputDir
			 )

		plotting.compare(hist_channels,
			 variable='all_s13',
			 setrange=(0,5999),
			 setlogy = setlogy,
			 scalelumi = scalelumi,
			 datalumi = datalumi,
			 drawRatio = drawRatio,
			 ratioLabel = ratioLabel,
			 normalize = normalize,
			 nRebin =100,
			 draw_channel_info= draw_channel_info,
			 output_dir= outputDir
			 )

		plotting.compare(hist_channels,
			 variable='all_s14',
			 setrange=(0,5999),
			 setlogy = setlogy,
			 scalelumi = scalelumi,
			 datalumi = datalumi,
			 drawRatio = drawRatio,
			 ratioLabel = ratioLabel,
			 normalize = normalize,
			 nRebin =100,
			 draw_channel_info= draw_channel_info,
			 output_dir= outputDir
			 )

		plotting.compare(hist_channels,
			 variable='all_s23',
			 # setrange=(0,20),
			 setlogy = setlogy,
			 scalelumi = scalelumi,
			 datalumi = datalumi,
			 drawRatio = drawRatio,
			 ratioLabel = ratioLabel,
			 normalize = normalize,
			 # nRebin =2,
			 draw_channel_info= draw_channel_info,
			 output_dir= outputDir
			 )

		plotting.compare(hist_channels,
			 variable='all_s24',
			 # setrange=(0,20),
			 setlogy = setlogy,
			 scalelumi = scalelumi,
			 datalumi = datalumi,
			 drawRatio = drawRatio,
			 ratioLabel = ratioLabel,
			 normalize = normalize,
			 # nRebin =2,
			 draw_channel_info= draw_channel_info,
			 output_dir= outputDir
			 )

		plotting.compare(hist_channels,
			 variable='all_s34',
			 # setrange=(0,20),
			 setlogy = setlogy,
			 scalelumi = scalelumi,
			 datalumi = datalumi,
			 drawRatio = drawRatio,
			 ratioLabel = ratioLabel,
			 normalize = normalize,
			 # nRebin =2,
			 draw_channel_info= draw_channel_info,
			 output_dir= outputDir
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


	if not os.path.exists(outputDir): os.mkdir(outputDir)
	if not os.path.exists(outputDir + "plots/"): os.mkdir(outputDir + "plots/")
	if not os.path.exists(outputDir + "Cutflows/"): os.mkdir(outputDir + "Cutflows/")
	if not os.path.exists(outputDir + "2Dmassplots"): os.mkdir(outputDir + "2Dmassplots")
	if not os.path.exists(outputDir + "2Dmassplots/data/"): os.mkdir(outputDir + "2Dmassplots/data/")
	if not os.path.exists(outputDir + "2Dmassplots/data/"): os.mkdir(outputDir + "2Dmassplots/data/")




	with open(options.config, 'r') as json_config:
		config_file = json.load(json_config) # load JSON config file


	#execute plotting here, comment out functions in you dont want to plot them again.	
	makeCutflows(config_file)
	compare_histograms(config_file, 'DVtype')
	compare_TRUTH_histograms(config_file, 'truth')
	# check_rerunningVSI(config_file,"all")

	
	