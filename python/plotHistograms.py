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
outputDir = '../output/' # change path here to save your histograms somewhere else!
MATERIAL_LAYERS = [33.25, 50.5, 88.5, 122.5, 299]
lumi = 60
normalize = False
setlogy = False

#############################################################################################################################################


def makeCutflows(config_file):
	vtx_channels = ["VSI", "VSI_Leptons"]
	for vtx_channel in vtx_channels:
		plotting.plot_cutflow(file = config_file["dataFile"],
							  vertextype= vtx_channel,
							  output_dir=outputDir + "Cutflows/")

		plotting.plot_cutflow(file = config_file["mcFiles"][0],
							  vertextype= vtx_channel,
							  output_dir=outputDir + "Cutflows/")

		plotting.plot_cutflow(file = config_file["mcFiles"][1],
							  vertextype= vtx_channel,
							  output_dir=outputDir + "Cutflows/")

def compare_histograms(config_file, selection):
	vtx_channels = ["VSI_Leptons"]
	for vtx_channel in vtx_channels:
		hist_channels = []
		# hist_channels[i] = (<filename>, <legend label>,<vertex directory>, <selection directory>)
		hist_channels.append([config_file["dataFile"],config_file["dataLabel"], vtx_channel, selection])
		hist_channels.append([config_file["mcFiles"][0],config_file["mcLabels"][0], vtx_channel, selection])
		hist_channels.append([config_file["mcFiles"][1], config_file["mcLabels"][1], vtx_channel, selection])

		#get integrated luminosity corresponding to the data file
		lumi = config_file["dataLumi"]

		# DV Variables
		plotting.compare(hist_channels,
						 variable='DV_r',
						 nRebin=10,
						 setrange=(0, 350),
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
						 vertical_lines=MATERIAL_LAYERS,
						 vertical_legend="Material Layers",
						 output_dir= outputDir,
						 use_ntuple = True,
						 ntup_nbins=100
						 )

		# DV Track Variables 
		plotting.compare(hist_channels,
						 variable='DV_trk_pt',
						 setrange=(0, 100),
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_eta',
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_phi',
						 nRebin = 2, 
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_eta',
						 setrange=(-3, 3),
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_d0',
						 setrange=(-10, 10),
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_z0',
						 setrange=(-10, 10),
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_dpt',
						 setrange=(0, 20),
						 nRebin = 2, 
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_deta',
						 setrange=(0, 3.2),
						 nRebin = 2, 
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_dphi',
						 setrange=(0, 3.2),
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
						 output_dir = outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_trk_dR',
						 setrange=(0, 10),
						 nRebin = 10, 
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
						 output_dir = outputDir
						 )

		# Mass Variables 
		plotting.compare(hist_channels,
						 variable='DV_mass',
						 setrange=(0, 20),
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
						 vertical_lines=[4],
						 vertical_legend="DV mass cut",
						 output_dir= outputDir
						 )

		plotting.compare(hist_channels,
						 variable='mvis',
						 nRebin=5,
						 setrange=(0, 200),
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
						 output_dir= outputDir
						 )

		plotting.compare(hist_channels,
						 variable='mtrans',
						 nRebin=5,
						 setrange=(0, 200),
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
						 output_dir= outputDir
						 )
		plotting.compare(hist_channels,
						 variable='HNLm',
						 setrange=(0, 30),
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
						 output_dir= outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_redmass',
						 setrange=(0, 50),
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
						 output_dir= outputDir
						 )

		plotting.compare(hist_channels,
						 variable='DV_redmassvis',
						 setrange=(0, 200),
						 nRebin=5,
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
						 output_dir= outputDir
						 )
		plotting.compare(hist_channels,
						 variable='DV_redmassHNL',
						 setrange=(0, 50),
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
						 output_dir= outputDir
						 )

		# HNL quantities
		plotting.compare(hist_channels,
						 variable='HNLpt',
						 setrange=(0, 50),
						 nRebin=5,
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
						 output_dir= outputDir
						 )

		plotting.compare(hist_channels,
						 variable='HNLeta',
						 setrange=(0, 50),
						 nRebin=5,
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
						 output_dir= outputDir
						 )

		plotting.compare(hist_channels,
						 variable='HNLphi',
						 setrange=(-3,3),
						 setlogy = setlogy,
						 lumi = lumi,
						 normalize = normalize,
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

	if os.path.exists(outputDir):
		pass
	else: 
		logger.info('Making output directories...')
		os.mkdir(outputDir)
		os.mkdir(outputDir + "plots/")
		os.mkdir(outputDir + "Cutflows/")
		os.mkdir(outputDir + "2Dmassplots")
		os.mkdir(outputDir + "2Dmassplots/data/")
		os.mkdir(outputDir + "2Dmassplots/mc/")



	with open(options.config, 'r') as json_config:
		config_file = json.load(json_config) # load JSON config file


	#execute plotting here, comment out functions in you dont want to plot them again.	
	makeCutflows(config_file)
	compare_histograms(config_file, 'DVtype')

	
	