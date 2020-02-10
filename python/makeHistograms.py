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
	output_path ="/eos/home-r/rnewhous/public/HNL/mc/21.0.77_full/histograms/"
	if os.path.exists(output_path) == False:
		logger.info('Making output directory')
		os.mkdir(output_path)

	if options.update == "False":
		if os.path.exists(output_path +"histograms.root"): 
			logger.info('Removing histograms.root')
			os.remove(output_path + "histograms.root") # by default remove histrogram file that if you previously created it.

	###########################################################################################################################
	# Put a map for a 1 one word key to a list of inputs for the selections.
	# histograms will be saved for each channel called histName_ch in histograms.root
	###########################################################################################################################
	channels = {
		# 'emu' : ['alltriggers','pmuon', '4-filter', 'nDV', 'fidvol','2track','OS', 'emu','2-tight','cosmicveto', 'mlll','DVmass'],
		#  'OS' : ['alltriggers','pmuon', '4-filter', 'nDV', 'fidvol','2track','OS'],
		#  'SS' : ['alltriggers','pmuon', '4-filter', 'nDV', 'fidvol','2track','SS']
		'TN': ['alltriggers', 'TN'],
		'FN': ['alltriggers', 'FN'],
		'FP': ['alltriggers', 'FP'],
		'TP': ['alltriggers', 'TP'],
	}


	analysisCode = {}
	anaClass = getattr(analysis, "WmuHNL")

	tree = treenames.Tree(options.file[0], "outTree")
	aod_tree = treenames.Tree(options.AODfile, "outTree")

	nentries = options.nevents or len(tree.dvmass)

	for k in channels:
		ch = k
		ana = anaClass(ch, channels[k],output_path+"histograms.root")
		for ievt in xrange(nentries):
			evt = helpers.Event(tree=tree, ievt = ievt , idv = None)
			aod_evt = helpers.Event(tree=aod_tree, ievt = ievt , idv = None) # for filter studies
			ndv = len(tree.dvx[ievt])

			presel = ana.preSelection(evt, aod_evt=aod_evt)

			for idv in xrange(ndv):
				DVevt = helpers.Event(tree=tree, ievt = ievt , idv = idv)
				ana.DVSelection(DVevt)

			ana.unlock()
		ana.end()

	# analysisCode["emu"] = anaClass("emu", channels["emu"],"histograms.root")

	# for ievt in xrange(nentries):
	# 	evt = helpers.Event(tree=tree, ievt = ievt , idv = None)
	# 	ndv = len(tree.dvx[ievt])

	# 	for ana in analysisCode.itervalues():
	# 		presel = ana.preSelection(evt)

	# 		for idv in xrange(ndv): 
	# 			DVevt = helpers.Event(tree=tree, ievt = ievt , idv = idv)
	# 			ana.DVSelection(DVevt)

	# 	ana.unlock()
	# analysisCode["emu"].end()



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
	parser.add_argument('--nevents',
                        default = None,
                        type = int,
                        help='Number of events are going to be processed for test-only purpose.')
	parser.add_argument("-u", "--update",
                        dest="update",
                        default="False",
                        help="Update histogram file? Default is to recreate.",
                        metavar="update")
	parser.add_argument("--AODfile",
						dest="AODfile",
						default="",
						help="AOD file to use in filter decision comparison studies.")


	parent_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter, parents = [parser]) 
    # subparsers = parent_parser.add_subparsers(title = 'SL/FH ttbar resonances anaylsis', dest = 'analysis')


	options = parent_parser.parse_args()
	logger.info("-> Calling main")
	# helpers.initialise_binds()
	main()
	logger.info("The end.")








