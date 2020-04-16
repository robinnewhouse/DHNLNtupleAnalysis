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

blinded = True # Dont dont change this flag! This ensures you do not accidentilly unblind when looking at data. 

def main():
	
	output_path ="../output/"
	if os.path.exists(output_path) == False:
		logger.info('Making output directory')
		os.mkdir(output_path)

	# if options.update == False:
	# 	if os.path.exists(output_path +"histograms.root"):
	# 		logger.info('Removing histograms.root')
	# 		os.remove(output_path + "histograms.root") # by default remove histrogram file that if you previously created it.

	with open(options.config, 'r') as json_config:
		config_file = json.load(json_config) # load JSON config file that contains a channel name mapped to a list of selections


	analysisCode = {}
	# Define that we're using a specific type of anaysis
	anaClass = getattr(analysis, "WmuHNL")

	file = options.input[0]
	treename = "outTree"

	#loop over all the channels in the config file
	for channel, configs in config_file.items():
		

		logger.info('Running on channel: %s'%channel)
		#create one output file per channel in your config file
		if "data" in options.config.split("config")[1]:
			outputfile = output_path + "histograms_data_%s.root"%channel
		else:
			outputfile = output_path + helpers.Output_filename(file, channel)

		if os.path.exists(outputfile):
			if options.force == False:
				logger.error("Output %s file already exists. Either re-run with -f/--force OR choose a different output path."%helpers.Output_filename(file, channel))
				exit()
			else:
				logger.info('Removing %s'%outputfile)
				os.remove(outputfile) # if force option is given then remove histrograms file that was previously created.

		# loop over the vertex containers in each channel (usually just VSI & VSI Leptons)
		for vtx_container in config_file[channel]["vtx_containers"]:
			
			selections =  config_file[channel]["selections"] # define selections for the channel from the config file
			tree = treenames.Tree(file, treename, vtx_container) # define variables in tree to be accessed from rootfile
			nentries = options.nevents or len(tree.dvmass)
			
			#blinding flag to prevent accidental unblinding in data
			if blinded:
				if tree.isData and "OS" in selections:
					logger.error("You are running on data and you cannot look at OS verticies!!!")
					exit()

			# Make instance of the analysis class
			ana = anaClass(vtx_container, selections, outputfile)
			
			# Loop over each event
			for ievt in xrange(nentries):
				if (ievt % 1000 == 0):
					logger.info("Channel {}: processing event {}".format("%s_%s"%(channel,vtx_container), ievt))
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
			analysisCode["%s_%s"%(channel,vtx_container)] = ana


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

	parser.add_argument("-i", "--input",
						dest="input",
						required = True,
						action = AppendActionCleanDefault,
						type = str,
						help="Input ntuple produced by DHNLAlgorithm.",
						metavar="INPUT")

	parser.add_argument("-f", "--force",
						action="store_true",
						dest="force",
						help="Overwrite previous histograms output file if it exists. (default: False)")

	parser.add_argument("--config",
						dest="config",
						type = str,
						required = True,
						help="Input config file for makeHisotgrams.py.")

	parser.add_argument('--nevents',
						default = None,
						type = int,
						help='Number of events are going to be processed for test-only purpose.')

	# parser.add_argument("-u", "--update",
	# 					action="store_true",
	# 					dest="update",
	# 					help="Update histogram file? Default is to recreate?")





	parent_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter, parents = [parser])
	# subparsers = parent_parser.add_subparsers(title = 'SL/FH ttbar resonances anaylsis', dest = 'analysis')


	options = parent_parser.parse_args()
	logger.info("-> Calling main")
	# helpers.initialise_binds()
	main()
	logger.info("The end.")








