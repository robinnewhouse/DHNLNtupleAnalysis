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
	

	# put a map for a 1 one word key to a list of inputs for the selections
	channels = { 
			   'emu' : ['alltriggers','pmuon', '4-filter', 'nDV', 'fidvol','2track','OS', 'emu','2-tight','cosmicveto', 'mlll', 'DVmass'],   
			   'mumu'  : ['alltriggers','pmuon', '4-filter' 'mumu']}


	analysisCode = {}
	anaClass = getattr(analysis, "WmuHNL")

	# file = fileName
	file = "/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_20G_lt10dd_emu_wTrigMatch.root"

	treename = "outTree"
	tree = treenames.Tree(file, treename)

	nentries = len(tree.dvmass)

	# for k in channels:
	# 	ch = k
	# 	analysisCode[k] = anaClass(ch, channels[k])

	analysisCode["emu"] = anaClass("emu", channels["emu"],"histograms.root")

	for ievt in xrange(nentries):
		evt = helpers.Event(tree=tree, ievt = ievt , idv = None)
		ndv = len(tree.dvx[ievt])

		for ana in analysisCode.itervalues():
			presel = ana.preSelection(evt)
			# print presel

			for idv in xrange(ndv): 
				DVevt = helpers.Event(tree=tree, ievt = ievt , idv = idv)
				# if presel: # current analysis will only look at the DV if the presel is met (will help code run faster, but we can turn this off)
				ana.DVSelection(DVevt)

		ana.unlock()
	analysisCode["emu"].end()



if __name__ == "__main__":
	main()
	logger.info("The end.")









