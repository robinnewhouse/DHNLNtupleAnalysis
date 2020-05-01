#!/usr/bin/env python
import os,sys
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
		


	with open(options.config, 'r') as json_config:
		config_file = json.load(json_config) # load JSON config file that contains a channel name mapped to a list of selections


	analysisCode = {}
	# Define that we're using a specific type of anaysis
	# anaClass = getattr(analysis, "oldHNLanalysis")
	anaClass = getattr(analysis, options.analysis)

	file = options.input[0] # get file 
	treename = "outTree" # define tree name 

	#loop over all the channels in the config file
	for channel, configs in config_file.items():
   		
		logger.info('Running on channel: %s'%channel)
		file_info = helpers.File_info(file, channel) # If you are running on MC this will give info about signal mass and lifetime

		#create one output file per channel in your config file
		if "data" in options.config.split("config")[1]:
			outputfile = output_path + "histograms_data_%s.root"%channel
		else:
			outputfile = output_path + file_info.Output_filename
		if os.path.exists(outputfile):
			if options.force == False:
				if "data" in options.config.split("config")[1]:
					logger.error("Output histograms_data_%s.root file already exists. Either re-run with -f/--force OR choose a different output path."%channel)
				else:
					logger.error("Output %s file already exists. Either re-run with -f/--force OR choose a different output path."%file_info.Output_filename)
				exit()
			else:
				logger.info('Removing %s'%outputfile)
				os.remove(outputfile) # if force option is given then remove histrograms file that was previously created.

		# loop over the vertex containers in each channel (usually just VSI & VSI Leptons)
		for vtx_container in config_file[channel]["vtx_containers"]:
			
			selections =  config_file[channel]["selections"] # define selections for the channel from the config file
			# Try to load only the number of entries of you need
			nentries = options.nevents if options.nevents else None
			tree = treenames.Tree(file, treename, vtx_container, nentries) # define variables in tree to be accessed from rootfile
			if len(tree.dvmass) < nentries or nentries == None:
					nentries = len(tree.dvmass)

			#blinding flag to prevent accidental unblinding in data
			if blinded:
				if tree.isData: 
					if "CR" in selections:
						pass
					else: 
						if "OS" in selections:
							logger.error("You are running on data and you cannot look at OS verticies!!!")
							sys.exit(1)  # abort because of error
						if ("SS" in selections) == False:
							logger.error("You are running on data and you are not in the CR. You must only look at SS vertices!!!")
							sys.exit(1)  # abort because of error


			# Make instance of the analysis class
			ana = anaClass(vtx_container, selections, outputfile,isdata=tree.isData)

			
			# Loop over each event
			for ievt in range(nentries):
				if (ievt % 1000 == 0):
					logger.info("Channel {}: processing event {}".format("%s_%s"%(channel,vtx_container), ievt))
				# Create an event instance to keep track of basic event properties
				evt = helpers.Event(tree=tree, ievt=ievt,mass=file_info.mass,ctau=file_info.ctau)
				ndv = len(tree.dvx[ievt])

				# Run preselection cuts to avoid processing unnecessary events
				presel = ana.preSelection(evt)

				# Loop over each vertex in the event
				for idv in range(ndv):
					DVevt = helpers.Event(tree=tree, ievt=ievt, idv=idv,mass=file_info.mass,ctau=file_info.ctau)
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
	parser.add_argument('-a','--analysis',
						dest="analysis",
						default = "oldAnalysis",
						type = str,
						help='Name of the analysis you want to run. Default is the old 36fb dHNL analysis')

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








