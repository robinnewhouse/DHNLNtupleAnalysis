#!/usr/bin/env python
import os
import helpers
import ROOT
import csv
import analysis
import treenames
import json
# import reweighting
# import systematics

logger = helpers.getLogger('dHNLAnalysis.makeHistograms')

blinded = True

def main():

	output_path ="../output/"
	if os.path.exists(output_path) == False:
		logger.info('Making output directory')
		os.mkdir(output_path)

	if options.update == False:
		if os.path.exists(output_path +"histograms.root"):
			logger.info('Removing histograms.root')
			os.remove(output_path + "histograms.root") # by default remove histrogram file that if you previously created it.

	with open(options.config, 'r') as json_config:
		config_file = json.load(json_config) # load JSON config file that contains a channel name mapped to a list of selections


	analysisCode = {}
	# Define that we're using a specific type of anaysis
	anaClass = getattr(analysis, "WmuHNL")

	file = options.file[0]
	treename = "outTree"

	#loop over all the configurations in the config file
	for k, configs in config_file.items():

		channels =  config_file[k]["channels"] # define channels for each section in config file
		vtx_container =  config_file[k]["vtx_container"] # define vertex container to be used 

		tree = treenames.Tree(file, treename, vtx_container) # define variables to be accessed from rootfile
		nentries = options.nevents or len(tree.dvmass)

		for channel, selections in channels.items():
			# Make instance of the analysis class
			if blinded: 
				if tree.isData and "OS" in selections:
					logger.fatal("You are running on data and you cannot look at OS verticies!!!")
					exit()
			ana = anaClass(channel, selections, output_path + "histograms_MC10G_10mm_el.root")
			
			# Loop over each event
			for ievt in xrange(nentries):
				if (ievt % 1000 == 0):
					print "Channel {}: processing event {}".format(channel, ievt)
				# Create an event instance to keep track of basic event properties
				evt = helpers.Event(tree=tree, ievt=ievt, idv=None)
				ndv = len(tree.dvx[ievt])

				# Run preselection cuts to avoid processing unnecessary events
				presel = ana.preSelection(evt)

				# Loop over each vertex in the event
				for idv in xrange(ndv):
					DVevt = helpers.Event(tree=tree, ievt=ievt, idv=idv)
					ana.DVSelection(DVevt)

				ana.unlock()
			# Call functions to finalize analysis
			ana.end()
			# Store analysis in dictionary for possible later use
			analysisCode[channel] = ana


if __name__ == "__main__":
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



	###########################################################################################################################
	# Parser arguments are not fully functional yet, but we should implement arguments like this to help with running the code.
	###########################################################################################################################

	# parser.add_argument("-d", "--data",
	# 					action="store_true",
	# 					dest="data",
	# 					help="Is this data?")

	parser.add_argument("-f", "--file",
						dest="file", default=["/eos/atlas/atlascerngroupdisk/phys-exotics/ueh/HNL/DHNLAlg_testNtuples/newframework_Ntuple_WmuHNL_20G_lt10dd_emu_wTrigMatch.root"],
						action = AppendActionCleanDefault,
						type = str,
						help="Input file",
						metavar="FILE")

	parser.add_argument("--config",
						dest="config",
						type = str,
						required = True,
						help="Input config file for makeHisotgrams.py.")

	parser.add_argument('--nevents',
						default = None,
						type = int,
						help='Number of events are going to be processed for test-only purpose.')
	parser.add_argument("-u", "--update",
						action="store_true",
						dest="update",
						help="Update histogram file? Default is to recreate?")





	parent_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter, parents = [parser])
	# subparsers = parent_parser.add_subparsers(title = 'SL/FH ttbar resonances anaylsis', dest = 'analysis')


	options = parent_parser.parse_args()
	logger.info("-> Calling main")
	# helpers.initialise_binds()
	main()
	logger.info("The end.")








