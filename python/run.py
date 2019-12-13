# import treenames	
# import selections 
# import helpers
# import ROOT, time
import EventSel
# from ROOT import *
# gROOT.LoadMacro("AtlasStyle.C")
# gROOT.LoadMacro("AtlasUtils.C")
# gROOT.LoadMacro("AtlasLabels.C")



if __name__ == '__main__':
	
	# emu files: 
	# file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/fixElmatching_mc16_13TeV.311660.Pythia8EvtGen_A14NNPDF23LO_WmuHNL50_20G_lt10dd_el.merge.DAOD_RPVLL.e7422_e5984_a875_r10739_r10706.root"
	# mumu files:
	# file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_10G_lt10dd_mumu.root"


	file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_20G_lt10dd_emu_wTightmuel.root"
	EventSel.Event_Sel(file, "311660", "muon", "emu", "20", "10", "partial MC16d", "VSI")
	
	file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_20G_lt10dd_emu_VSILep.root"
	EventSel.Event_Sel(file, "311660", "muon", "emu", "20", "10", "partial MC16d", "VSILep")
	
	file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_10G_lt10dd_mumu_wTightmuel.root"
	EventSel.Event_Sel(file, "311633", "muon", "mumu", "20", "10", "partial MC16d", "VSI")

	file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_10G_lt10dd_mumu_VSILep.root"
	EventSel.Event_Sel(file, "311633", "muon", "mumu", "10", "10", "partial MC16d", "VSILep")

