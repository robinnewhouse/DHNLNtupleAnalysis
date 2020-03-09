#!/usr/bin/env python
import os
import helpers
import ROOT
import csv
import analysis
import treenames
# import reweighting
# import systematics

logger = helpers.getLogger('dHNLAnalysis.makeHistograms')

def main():
	if options.output_file:
		output_path = os.path.abspath(os.path.join(options.output_file, os.pardir))
		output_file = options.output_file
	else:
		output_path = "../output/"
		output_file = output_path + "histograms.root"

	logger.info('Output file: ' + output_file)

	# Make output directory if it doesn't exist
	if not os.path.exists(output_path):
		logger.info('Making output directory')
		os.mkdirs(output_path)

	if not options.update and os.path.exists(output_file):
		logger.info('Removing output file')
		os.remove(output_file) # by default remove histogram file that if you previously created it.


	###########################################################################################################################
	# Put a map for a 1 one word key to a list of inputs for the selections.
	# histograms will be saved for each channel called histName_ch in histograms.root
	###########################################################################################################################
	channels = {
		# 'mumu': ['alltriggers', 'pmuon', '4-filter', 'nDV', 'fidvol', '2track', 'OS', 'mumu', '2-tight', 'cosmicveto', 'mlll', 'DVmass'],
		# 'mumu': ['alltriggers', 'pmuon', 'mumu-filter', 'nDV', 'fidvol', '2track', 'OS', 'mumu', '2-tight', 'cosmicveto', 'mlll', 'DVmass'],
		'mumu_TP': ['alltriggers', 'pmuon', 'mumu-filter', 'TP', 'nDV', 'fidvol', '2track', 'OS', 'mumu', '2-tight', 'cosmicveto', 'mlll', 'DVmass'],
		'mumu_FP': ['alltriggers', 'pmuon', 'mumu-filter', 'FP', 'nDV', 'fidvol', '2track', 'OS', 'mumu', '2-tight', 'cosmicveto', 'mlll', 'DVmass'],
		'mumu_FN': ['alltriggers', 'pmuon', 'mumu-filter', 'FN', 'nDV', 'fidvol', '2track', 'OS', 'mumu', '2-tight', 'cosmicveto', 'mlll', 'DVmass'],
		# 'mue': ['alltriggers', 'pmuon', '4-filter', 'nDV', 'fidvol', '2track', 'OS', 'mue', '2-tight', 'cosmicveto', 'mlll', 'DVmass'],
		# 'emu': ['alltriggers', 'pmuon', '4-filter', 'nDV', 'fidvol', '2track', 'OS', 'emu', '2-tight', 'cosmicveto', 'mlll', 'DVmass'],
		# 'OS' : ['alltriggers','pmuon', '4-filter', 'nDV', 'fidvol','2track','OS'],
		# 'SS' : ['alltriggers','pmuon', '4-filter', 'nDV', 'fidvol','2track','SS']
	}


	analysisCode = {}
	# Define that we're using a specific type of anaysis
	anaClass = getattr(analysis, "WmuHNL")

	logger.info('Loading ntuple tree from ' + options.file[0])
	tree = treenames.Tree(options.file[0],  "outTree")
	aod_tree = treenames.Tree(options.aod_file,  "outTree") if options.aod_file else None
	
	nentries = options.nevents or len(tree.dvmass)
	logger.info('Finished loading. Processing {} events'.format(nentries))

	for channel, selections in channels.items():
		# Make instance of the analysis class
		ana = anaClass(channel, selections, output_file)
		# Loop over each event
		for ievt in xrange(nentries):
			if (ievt % 1000 == 0):
				print "Channel {}: processing event {}".format(channel, ievt)
			# Create an event instance to keep track of basic event properties
			# aod_tree is temporarily necessary while doing filter mismatch tests 
			evt = helpers.Event(tree=tree, ievt=ievt, aod_tree=aod_tree) 
			ndv = len(tree.dvx[ievt])

			# Run preselection cuts to avoid processing unnecessary events
			presel = ana.preSelection(evt)

			# Loop over each vertex in the event
			for idv in xrange(ndv):
				DVevt = helpers.Event(tree=tree, ievt=ievt, idv=idv)
				ana.DVSelection(DVevt)

			ana.unlock()  # !!! Does this do anything now? RN
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
						action=AppendActionCleanDefault,
						type=str,
						help="Input file",
						metavar="FILE")
	parser.add_argument("--output_file",
						dest="output_file",
						default="",
						help="Path of output histograms file.")
	parser.add_argument('--nevents',
						default=None,
						type=int,
						help='Number of events are going to be processed for test-only purpose.')
	parser.add_argument("-u", "--update",
						action="store_true",
						dest="update",
						help="Update histogram file? Default is to recreate?")
	parser.add_argument("--aod_file",
						dest="aod_file",
						default="",
						help="AOD file to use in filter decision comparison studies.")

	parent_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, parents=[parser])
	# subparsers = parent_parser.add_subparsers(title = 'SL/FH ttbar resonances anaylsis', dest = 'analysis')

	options = parent_parser.parse_args()
	logger.info("-> Calling main")
	# helpers.initialise_binds()
	main()
	logger.info("The end.")
