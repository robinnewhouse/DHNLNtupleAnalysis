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


# ------------------------------------------------------------------------------------------------------------------------------------------------------
	# files with no trigger matching 

	# file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_20G_lt10dd_emu_wTightmuel.root"
	# EventSel.Event_Sel(file, "311660", "muon", "emu", "20", "10", "partial MC16d", "VSI")
	
	# file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_20G_lt10dd_emu_VSILep.root"
	# EventSel.Event_Sel(file, "311660", "muon", "emu", "20", "10", "partial MC16d", "VSILep")
	
	# file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_10G_lt10dd_mumu_wTightmuel.root"
	# EventSel.Event_Sel(file, "311633", "muon", "mumu", "20", "10", "partial MC16d", "VSI")

	# file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_10G_lt10dd_mumu_VSILep.root"
	# EventSel.Event_Sel(file, "311633", "muon", "mumu", "10", "10", "partial MC16d", "VSILep")
	

	# file = "/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WelHNL_10G_10mm_emu.root"
	# EventSel.Event_Sel(file, "312990", "electron", "emu", "10", "10", "partial MC16a", "VSI")

	# file = "/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WelHNL_10G_10mm_emu_VSILep.root"
	# EventSel.Event_Sel(file, "312990", "electron", "emu", "10", "10", "partial MC16a", "VSILep")

	# file = "/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WelHNL_10G_10mm_ee_VSILep.root"
	# EventSel.Event_Sel(file, "312987", "electron", "ee", "10", "10", "partial MC16a", "VSILep")
# ------------------------------------------------------------------------------------------------------------------------------------------------------
	# with trigger matching

	# file = "/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_20G_lt10dd_emu_wTrigMatch.root"
	# file = "/home/dtrischuk/HNLAnalysis/DHNL/run/testData/data-tree/data15_13TeV.00284427.physics_Main.merge.DAOD_RPVLL.r11761_r11764_p4054.root"
	# EventSel.Event_Sel(file, "311660", "muon", "emu", "20", "10", "partial MC16d", "VSI")
	
# ------------------------------------------------------------------------------------------------------------------------------------------------------
	# data runs
	file = "/data/hnl/dataprep/ntuples/VSI/data15_13TeV.00284427.physics_Main.merge.DAOD_RPVLL.r11761_r11764_p4054.root"
	EventSel.Event_Sel(file, "311660", "muon", "emu", "-1", "-1", "2015 test data", "VSI")
	
	file = "/data/hnl/dataprep/ntuples/VSI/data16_13TeV.00304178.physics_Main.merge.DAOD_RPVLL.r11761_r11764_p4054.root"
	EventSel.Event_Sel(file, "311660", "muon", "emu", "-1", "-1", "2016 test data", "VSI")
	
	file = "/data/hnl/dataprep/ntuples/VSI/data17_13TeV.00340368.physics_Main.merge.DAOD_RPVLL.r11761_r11764_p4054.root"
	EventSel.Event_Sel(file, "311660", "muon", "emu", "-1", "-1", "2017 test data", "VSI")

	file = "/data/hnl/dataprep/ntuples/VSI/data18_13TeV.00358031.physics_Main.merge.DAOD_RPVLL.r11760_r11764_p4054.root"
	EventSel.Event_Sel(file, "311660", "muon", "emu", "-1", "-1", "2018 test data", "VSI")



	file = "/data/hnl/dataprep/ntuples/VSILep/data15_13TeV.00284427.physics_Main.merge.DAOD_RPVLL.r11761_r11764_p4054.root"
	EventSel.Event_Sel(file, "311660", "muon", "emu", "-1", "-1", "2015 test data", "VSILep")
	
	file = "/data/hnl/dataprep/ntuples/VSILep/data16_13TeV.00304178.physics_Main.merge.DAOD_RPVLL.r11761_r11764_p4054.root"
	EventSel.Event_Sel(file, "311660", "muon", "emu", "-1", "-1", "2016 test data", "VSILep")
	
	file = "/data/hnl/dataprep/ntuples/VSILep/data17_13TeV.00340368.physics_Main.merge.DAOD_RPVLL.r11761_r11764_p4054.root"
	EventSel.Event_Sel(file, "311660", "muon", "emu", "-1", "-1", "2017 test data", "VSILep")

	file = "/data/hnl/dataprep/ntuples/VSILep/data18_13TeV.00358031.physics_Main.merge.DAOD_RPVLL.r11760_r11764_p4054.root"
	EventSel.Event_Sel(file, "311660", "muon", "emu", "-1", "-1", "2018 test data", "VSILep")



