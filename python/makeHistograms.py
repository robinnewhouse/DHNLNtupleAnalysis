#!/usr/bin/env python
import os, sys
import helpers
import analysis
import trees
import json
import time
# import reweighting
# import systematics



blinded = True  # Dont dont change this flag! This ensures you do not accidentally unblind when looking at data.


def main():
	start = time.clock()
	#set debug level
	helpers.logger_debug_level = helpers.get_debug_level(options.debug_level)
	# set up own logger
	logger = helpers.getLogger('dHNLAnalysis.makeHistograms', level=helpers.logger_debug_level)
	# set up logger for helper module
	helpers.logger.setLevel(helpers.logger_debug_level)
	logger.info("-> Calling main")

	with open(options.config, 'r') as json_config:
		# load JSON config file that contains a channel name mapped to a list of selections
		config_file = json.load(json_config)

	output_path = os.path.join(os.path.abspath("../output/"+ options.output), '')
	if options.output_file:
		output_path = os.path.join(os.path.dirname(options.output_file), '')  # override
		if len(config_file.items()) > 1:
			logger.error("Writing multiple channels to one file. This will overwrite previous channels. Please reconsider.")
			exit(1)

	if not os.path.exists(output_path):
		logger.info('Making output directory')
		os.mkdir(output_path)
	analysisCode = {}
	# Define that we're using a specific type of analysis
	anaClass = getattr(analysis, options.analysis)

	input_file = options.input[0]  # get file
	logger.info("Running event selection on: {}".format(input_file))
	treename = "outTree"  # define tree name
	

	# loop over all the channels in the config file
	for channel, configs in config_file.items():

		logger.info('Running on channel: {}'.format(channel))
		# If you are running on MC this will give info about signal mass and lifetime
		file_info = helpers.FileInfo(input_file, channel)
		# Try to load only the number of entries of you need
		entries = options.nevents if options.nevents else None
		# Create new Tree class using uproot
		tree = trees.Tree(input_file, treename, entries, mc_campaign=file_info.MC_campaign, mass=file_info.mass,
											channel=channel, ctau=file_info.ctau, not_hnl_mc=options.notHNLmc, skip_events=options.skipEvents, br=file_info.br)
		logger.info('Mass dependent BR: {}'.format(file_info.br))

		# create one output file per channel in your config file
		if "SSbkg" in options.config.split("config")[1]:
			output_file = output_path + "histograms_SSbkg_{}.root".format(channel)
		else:
			if "CR" in config_file[channel]["selections"]:
				if tree.is_data:
					output_file = output_path + "CR_histograms_data_{}.root".format(channel)
				else:
					output_file = output_path + "CR_" + file_info.output_filename
			elif  "CR_BE" in config_file[channel]["selections"]:
				if tree.is_data:
					output_file = output_path + "CR_BE_histograms_data_{}.root".format(channel)
				else:
					output_file = output_path + "CR_BE" + file_info.output_filename
			elif "BE" in config_file[channel]["selections"]:
				if tree.is_data:
					output_file = output_path + "BE_histograms_data_{}.root".format(channel)
				else:
					output_file = output_path + "CR_BE" + file_info.output_filename
			else:
				output_file = output_path + file_info.output_filename

		# override if available
		if options.output_file:output_file = os.path.abspath(options.output_file)

		if os.path.exists(output_file):
			if not options.force:
				logger.error("Output {} file already exists. Either re-run with -f/--force OR choose a different output path.".format(output_file))
				exit()
			else:
				logger.info('Removing {}'.format(output_file))
				os.remove(output_file)  # if force option is given then remove histograms file that was previously created.

		
		if entries is None or tree.numentries < entries:
			entries = tree.numentries
		# specify this to reduce number of entries loaded in each array
		if options.skipEvents is not None: entries = entries + options.skipEvents  # if skipping events then entries needs to be updated
		tree.max_entries = entries
		logger.info('Going to process {}  events'.format(entries))

		# loop over the vertex containers in each channel (usually just VSI & VSI Leptons)
		for vtx_container in config_file[channel]["vtx_containers"]:

			# define selections for the channel from the config file
			selections = config_file[channel]["selections"]

			# define variables in tree to be accessed from root file
			tree.vtx_container = vtx_container

			# blinding flag to prevent accidental unblinding in data
			do_CR = "CR" in selections or "CR_BE" in selections
			if blinded and tree.is_data and not do_CR :
				# if "CR_BE" in selections: 
				# 	pass 
				# else: 
				if "OS" in selections or "SS" not in selections:
					logger.error("You are running on data and you cannot look at OS vertices!!! "
								"Please include 'SS', not 'OS' in selections, "
								"or add 'CR' if you are trying to look in the control region.")
					sys.exit(1)  # abort because of error

			# Make instance of the analysis class
			ana = anaClass(options.analysis, tree, vtx_container, selections, output_file, options.saveNtuples, weight_override=options.weight)

			# Loop over each event
			while tree.ievt < entries:
				if tree.ievt % 1000 == 0:
					logger.info("Channel {}_{}: processing event {} / {}".format(channel, vtx_container, tree.ievt, entries))
				# Create an event instance to keep track of basic event properties
				# evt = helpers.Event(tree=tree, ievt=tree.ievt, mass=file_info.mass, ctau=file_info.ctau)

				# Run preselection cuts to avoid processing unnecessary events
				presel = ana.preSelection()

				# Loop over each vertex in the event
				while tree.idv < tree.ndv:
					# DVevt = helpers.Event(tree=tree, ievt=tree.ievt, idv=idv, mass=file_info.mass, ctau=file_info.ctau)
					ana.DVSelection()
					tree.increment_dv()

				tree.reset_dv()
				tree.increment_event()
				ana.unlock()

			# Call functions to finalize analysis
			tree.reset_event()
			ana.end()
			# Store analysis in dictionary for possible later use
			# This is a huge memory hog and will likely crash if too many histograms are declared
			# Recommended not to use unless necessary and unless a minimal number of histograms are written. # RN
			# analysisCode["%s_%s"%(channel,vtx_container)] = ana
	logger.info("The end.")
	logger.info("Time elapsed: {}".format(time.clock()-start))


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

	parser.add_argument("-o", "--output",
						dest="output",
						type=str,
						default = "",
						help="Output directory to store histograms.")

	parser.add_argument("--output_file",
						dest="output_file",
						type=str,
						default = "",
						help="Overrides the output filename and output directory.")

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

	parser.add_argument('--weight',
						default = None,
						type = float,
						help='Use this flag if you want to override the weight calculation for this sample.')

	parser.add_argument('-a','--analysis',
						dest="analysis",
						default = "run2Analysis",
						type = str,
						help='Name of the analysis you want to run. Default is the full run 2 139fb dHNL analysis')

	parser.add_argument('-s','--saveNtuples',
						dest="saveNtuples",
						default = "DVtype",
						type = str,
						help='Name of cut after which you want to save the micro-ntuples. Default is the save them after DV type cut')

	parser.add_argument('-d','--debug',
						dest="debug_level",
						default = 'INFO',
						type = str,
						help='debug level. Default is INFO. Options include are CRITICAL, ERROR, WARNING, INFO, DEBUG ')

	parser.add_argument('--notHNLmc',
						action="store_true",
						dest="notHNLmc",
						default = False,
						help='Not running on HNL mc. Default: False. Useful for running on mc that is not HNL mc. Turn HNL specific truth info storing off.')

	parser.add_argument('--skipEvents',
						dest="skipEvents",
						default = None,
						type = int,
						help='Skip this number of events when processing inputfile.')
	







	# parser.add_argument("-u", "--update",
	# 					action="store_true",
	# 					dest="update",
	# 					help="Update histogram file? Default is to recreate?")





	parent_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter, parents = [parser])
	# subparsers = parent_parser.add_subparsers(title = 'SL/FH ttbar resonances anaylsis', dest = 'analysis')


	options = parent_parser.parse_args()
	
	# helpers.initialise_binds()
	main()
	








