# from ROOT import*
import ROOT
import numpy as np
import helpers
import logging
# logger = helpers.getLogger('dHNLAnalysis.selections',level = logging.WARNING)


class Trigger():
	def __init__(self, evt, trigger,invert=False):
		self.evt = evt
		self.invert = invert

		# trigger lists taken from https://acode-browser1.usatlas.bnl.gov/lxr/source/athena/PhysicsAnalysis/SUSYPhys/LongLivedParticleDPDMaker/share/PhysDESDM_HNL.py?v=21.0#0008
		apiSingleMuonTriggerlist = ["HLT_mu20_iloose_L1MU15", "HLT_mu24_iloose", "HLT_mu24_ivarloose", "HLT_mu24_ivarmedium", "HLT_mu26_imedium", "HLT_mu26_ivarmedium", "HLT_mu40", "HLT_mu50", "HLT_mu60_0eta105_msonly"]
		apiSingleElectronTriggerlist = ["HLT_e24_lhmedium_L1EM20VH", "HLT_e24_lhtight_nod0_ivarloose", "HLT_e26_lhtight_nod0", "HLT_e26_lhtight_nod0_ivarloose", "HLT_e60_lhmedium_nod0", "HLT_e60_lhmedium", "HLT_e60_medium","HLT_e120_lhloose", "HLT_e140_lhloose_nod0", "HLT_e300_etcut"]
		apiMultiMuonTriggerlist = ["HLT_2mu10",  "HLT_2mu10_nomucomb",  "HLT_2mu14",  "HLT_2mu14_nomucomb",  "HLT_mu20_nomucomb_mu6noL1_nscan03",  "HLT_mu20_mu8noL1",  "HLT_mu22_mu8noL1",  "HLT_mu22_mu8noL1_calotag_0eta010",  "HLT_3mu4",  "HLT_mu6_2mu4",  "HLT_3mu6",  "HLT_3mu6_msonly",  "HLT_mu20_2mu4noL1",  "HLT_4mu4",  "HLT_mu11_nomucomb_2mu4noL1_nscan03_L1MU11_2MU6",  "HLT_mu11_nomucomb_2mu4noL1_nscan03_L1MU11_2MU6_bTau",  "HLT_mu20_msonly_mu10noL1_msonly_nscan05_noComb",  "HLT_mu11_nomucomb_mu6noL1_nscan03_L1MU11_2MU6_bTau",  "HLT_mu6_nomucomb_2mu4_nomucomb_bTau_L1MU6_3MU4",  "HLT_2mu6_nomucomb_mu4_nomucomb_bTau_L12MU6_3MU4"]
		apiMultiElectronTriggerlist = ["HLT_2e12_lhloose_L12EM10VH", "HLT_2e15_lhvloose_nod0_L12EM13VH", "HLT_e17_lhloose_2e9_lhloose", "HLT_e17_lhloose_nod0_2e9_lhloose_nod0", "HLT_e17_lhloose_nod0_2e10_lhloose_nod0_L1EM15VH_3EM8VH", "HLT_2e17_lhvloose_nod0", "HLT_2e17_lhvloose_nod0_L12EM15VHI", "HLT_2e24_lhvloose_nod0", "HLT_e24_lhvloose_nod0_2e12_lhvloose_nod0_L1EM20VH_3EM10VH"]
		apiElectronMuonTriggerlist = ["HLT_e17_lhloose_mu14", "HLT_e17_lhloose_nod0_mu14", "HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1", "HLT_e26_lhmedium_nod0_mu8noL1", "HLT_e7_lhmedium_nod0_mu24", "HLT_e12_lhloose_nod0_2mu10", "HLT_2e12_lhloose_mu10", "HLT_2e12_lhloose_nod0_mu10", "HLT_e7_lhmedium_mu24", "HLT_e12_lhloose_2mu10"]	
		self.DAOD_RPVLLTriggerlist = ['HLT_mu26_ivarmedium', 'HLT_xe110_pufit_xe65_L1XE50', 'HLT_e26_lhtight_nod0_ivarloose', 'HLT_mu50', 'HLT_mu60_0eta105_msonly', 'HLT_mu80_msonly_3layersEC', 'HLT_mu22_mu8noL1', 'HLT_xe120_pufit_L1XE50', 'HLT_xe110_pufit_xe70_L1XE50', 'HLT_j0_perf_ds1_L1J100', 'HLT_g25_medium_L1EM22VHI_4j35_0eta490_invm1000', 'HLT_5j70_0eta240_L14J15', 'HLT_5j50_gsc70_boffperf_split_0eta240_L14J15', 'HLT_6j45_gsc55_boffperf_split_0eta240_L14J15', 'HLT_7j45_L14J15', 'HLT_7j35_gsc45_boffperf_split_L14J15', 'HLT_2mu14', 'HLT_mu6_dRl1_mu20_msonly_iloosems_mu6noL1_dRl1_msonly', 'HLT_3j200', 'HLT_4j120', 'HLT_4j85_gsc115_boffperf_split', 'HLT_ht1000_L1J100', 'HLT_2g20_tight_icalovloose_L12EM15VHI', 'HLT_2g22_tight_L12EM15VHI', 'HLT_2g25_tight_L12EM20VH', 'HLT_e26_lhtight_nod0_e15_etcut_L1EM7_Zee', 'HLT_2e17_lhvloose_nod0_L12EM15VHI', 'HLT_2e24_lhvloose_nod0', 'HLT_e24_lhmedium_nod0_L1EM20VH_g25_medium', 'HLT_g35_medium_g25_medium_L12EM20VH', 'HLT_j30_jes_cleanLLP_PS_llp_L1LLP-NOMATCH', 'HLT_j30_jes_cleanLLP_PS_llp_noiso_L1LLP-NOMATCH', 'HLT_6j50_gsc70_boffperf_split_L14J15', 'HLT_6j55_0eta240_L14J15', 'HLT_2j35_bmv2c1060_split_2j35_L14J15.0ETA25', 'HLT_2j35_bmv2c1070_split_2j35_bmv2c1085_split_L14J15.0ETA25', 'HLT_g140_loose', 'HLT_j420', 'HLT_j360_a10t_lcw_jes_60smcINF_j360_a10t_lcw_jes_L1SC111', 'HLT_j370_a10t_lcw_jes_35smcINF_j370_a10t_lcw_jes_L1SC111', 'HLT_e140_lhloose_nod0', 'HLT_tau160_medium1_tracktwoEF_L1TAU100', 'HLT_j30_jes_cleanLLP_PS_llp_L1LLP-RO', 'HLT_j30_jes_cleanLLP_PS_llp_noiso_L1LLP-RO', 'HLT_2g50_loose_L12EM20VH', 'HLT_j460_a10_lcw_subjes_L1J100', 'HLT_j460_a10_lcw_subjes_L1SC111', 'HLT_j460_a10r_L1J100', 'HLT_j460_a10r_L1SC111', 'HLT_j420_a10t_lcw_jes_35smcINF_L1J100', 'HLT_5j85_L14J15', 'HLT_5j60_gsc85_boffperf_split_L14J15', 'HLT_6j70_L14J15', 'HLT_mu4_j90_xe90_pufit_2dphi10_L1MU4_J50_XE50_DPHI-J20s2XE30', 'HLT_mu4_j90_xe90_pufit_2dphi10_L1MU4_XE60', 'HLT_e5_lhloose_nod0_j50_xe70_pufit_2dphi10_L1J40_XE50_DPHI-J20s2XE50', 'HLT_g80_loose_xe80noL1', 'HLT_mu20_msonly_iloosems_mu6noL1_msonly_nscan05_L1MU20_XE30', 'HLT_mu20_msonly_iloosems_mu6noL1_msonly_nscan05_L1MU20_J40', 'HLT_e60_lhmedium_nod0', 'HLT_tau200_medium1_tracktwoEF_L1TAU100', 'HLT_mu20_mu6noL1_bTau', 'HLT_e17_lhmedium_nod0_ivarloose_tau25_medium1_tracktwoEF', 'HLT_mu20_2mu4noL1', 'HLT_mu11_2mu4noL1_bNocut_L1MU11_2MU6', 'HLT_g85_tight_L1EM22VHI_3j50noL1', 'HLT_e5_lhvloose_nod0_mu4_j30_xe40_pufit_2dphi10_L1MU4_J30_XE40_DPHI-J20s2XE30', 'HLT_g300_etcut', 'HLT_e300_etcut', 'HLT_j225_gsc420_boffperf_split', 'HLT_j460_a10t_lcw_jes_L1J100', 'HLT_j460_a10t_lcw_jes_L1SC111', 'HLT_2j330_a10t_lcw_jes_35smcINF_L1J100', 'HLT_2j330_a10t_lcw_jes_35smcINF_L1SC111', 'HLT_g45_loose_6j45_0eta240', 'HLT_tau35_medium1_tracktwoEF_xe70_L1XE45', 'HLT_j55_gsc75_bmv2c1040_split_3j55_gsc75_boffperf_split', 'HLT_j60_gsc85_bmv2c1050_split_3j60_gsc85_boffperf_split', 'HLT_j45_gsc55_bmv2c1050_split_ht700_L1HT190-J15s5.ETA21', 'HLT_mu14_ivarloose_tau25_medium1_tracktwoEF', 'HLT_mu14_ivarloose_tau25_medium1_tracktwoEF_L1DR-MU10TAU12I_TAU12I-J25', 'HLT_mu14_tau25_medium1_tracktwoEF_xe50', 'HLT_mu14_ivarloose_tau25_medium1_tracktwoEF_xe50', 'HLT_j30_muvtx', 'HLT_j30_muvtx_noiso', 'HLT_g35_loose_L1EM24VHI_mu18', 'HLT_2j35_gsc45_bmv2c1050_split_2j35_gsc45_boffperf_split_L14J15.0ETA25', 'HLT_2j45_gsc55_bmv2c1060_split_2j45_gsc55_boffperf_split_L14J15.0ETA25', 'HLT_2j35_bmv2c1060_split_3j35', 'HLT_2j35_gsc45_bmv2c1060_split_3j35_gsc45_boffperf_split', 'HLT_3j50_gsc65_bmv2c1077_split_L13J35.0ETA23', 'HLT_3j35_bmv2c1070_split_j35_L14J15.0ETA25', 'HLT_4j35_bmv2c1077_split_L14J15.0ETA25', 'HLT_3j35_bmv2c1077_split_xe70_pufit_xe50_L13J15.0ETA25_XE40', 'HLT_e17_lhloose_nod0_mu14', 'HLT_e26_lhmedium_nod0_mu8noL1', 'HLT_e7_lhmedium_nod0_mu24', 'HLT_g25_medium_mu24', 'HLT_2e5_lhvloose_nod0_j40_xe70_pufit_2dphi10_L1J40_XE50_DPHI-J20s2XE50', 'HLT_g35_tight_icalotight_L1EM24VHI_mu18noL1', 'HLT_e17_lhmedium_nod0_ivarloose_tau25_medium1_tracktwoEF_xe50', 'HLT_e17_lhmedium_nod0_tau25_medium1_tracktwoEF_xe50', 'HLT_e17_lhmedium_nod0_ivarloose_tau25_medium1_tracktwoEF_L1DR-EM15TAU12I-J25', 'HLT_e24_lhmedium_nod0_ivarloose_tau35_medium1_tracktwoEF', 'HLT_e17_lhmedium_nod0_ivarloose_tau25_medium1_tracktwoEF_L1EM15VHI_2TAU12IM_4J12', 'HLT_g25_medium_L1EM22VHI_j35_0eta490_bmv2c1077_split_3j35_0eta490_invm700', 'HLT_g25_medium_L1EM22VHI_2j35_0eta490_bmv2c1077_split_2j35_0eta490', 'HLT_g20_tight_icaloloose_j35_bmv2c1077_split_3j35_0eta490_invm500', 'HLT_j85_gsc100_bmv2c1050_split_xe85_pufit_xe50_L1XE55', 'HLT_g0_hiptrt_L1EM22VHI', 'HLT_j30_jes_cleanLLP_PS_llp_L1TAU100', 'HLT_j30_jes_cleanLLP_PS_llp_noiso_L1TAU100', 'HLT_e5_lhloose_nod0_j40_xe70_pufit_2dphi10_L1XE60', 'HLT_j110_gsc150_boffperf_split_2j45_gsc55_bmv2c1070_split_L1J85_3J30', 'HLT_tau40_medium1_tracktwoEF_tau35_medium1_tracktwoEF', 'HLT_tau35_medium1_tracktwoEF_tau25_medium1_tracktwoEF_L1DR-TAU20ITAU12I-J25', 'HLT_tau35_medium1_tracktwoEF_tau25_medium1_tracktwoEF_03dR30_L1DR-TAU20ITAU12I-J25', 'HLT_tau80_medium1_tracktwoEF_L1TAU60_tau35_medium1_tracktwoEF_L1TAU12IM_L1TAU60_DR-TAU20ITAU12I', 'HLT_tau80_medium1_tracktwoEF_L1TAU60_tau60_medium1_tracktwoEF_L1TAU40', 'HLT_2j35_gsc45_bmv2c1070_split_xe85_pufit_xe50_L12J15_XE55', 'HLT_j175_gsc225_bmv2c1040_split', 'HLT_j225_gsc275_bmv2c1060_split', 'HLT_j225_gsc300_bmv2c1070_split', 'HLT_j225_gsc360_bmv2c1077_split', 'HLT_tau60_medium1_tracktwoEF_tau25_medium1_tracktwoEF_xe50', 'HLT_2mu4_invm1_j20_xe40_pufit_2dphi10_L12MU4_J20_XE30_DPHI-J20s2XE30', 'HLT_2mu4_invm1_j20_xe60_pufit_2dphi10_L12MU4_J40_XE50', 'HLT_mu14_ivarloose_tau35_medium1_tracktwoEF_L1MU10_TAU20IM_J25_2J20', 'HLT_mu14_ivarloose_tau35_medium1_tracktwoEF', 'HLT_mu11_mu6_bDimu2700_L1LFV-MU11', 'HLT_mu11_mu6_bTau_L1LFV-MU11', 'HLT_mu11_mu6_bBmumuxv2_L1LFV-MU11', 'HLT_mu11_mu6_bDimu2700', 'HLT_mu11_mu6_bBmumuxv2', 'HLT_mu11_mu6_bTau', 'HLT_e12_lhloose_nod0_2mu10', 'HLT_g15_loose_2mu10_msonly', 'HLT_e5_lhmedium_nod0_mu4_j30_xe65_pufit_2dphi10_L1MU4_XE60', 'HLT_g35_tight_icalotight_L1EM24VHI_mu15noL1_mu2noL1', 'HLT_g35_loose_L1EM24VHI_mu15_mu2noL1', 'HLT_3mu6_msonly', 'HLT_e70_lhloose_nod0_xe70noL1', 'HLT_2e5_lhvloose_nod0_j40_xe70_pufit_2dphi10_L1XE60', 'HLT_2j45_gsc55_bmv2c1050_split_ht300_L1HT190-J15s5.ETA21', 'HLT_mu11_mu6_bDimu_L1LFV-MU11', 'HLT_mu11_mu6_bDimu2700_Lxy0_L1LFV-MU11', 'HLT_mu11_mu6_bDimu_Lxy0_L1LFV-MU11', 'HLT_mu11_mu6_bDimu_Lxy0', 'HLT_mu11_mu6_bDimu2700_Lxy0', 'HLT_mu11_mu6_bDimu', 'HLT_ht300_2j40_0eta490_invm700_L1HT150-J20s5.ETA31_MJJ-400-CF_AND_2j35_gsc45_bmv2c1070_split', 'HLT_j55_gsc80_bmv2c1070_split_j45_gsc60_bmv2c1085_split_j45_320eta490', 'HLT_j80_0eta240_j60_j45_320eta490_AND_2j35_gsc45_bmv2c1070_split', 'HLT_j150_gsc175_bmv2c1060_split_j45_gsc60_bmv2c1060_split', 'HLT_mu6_2mu4_bDimu2700', 'HLT_mu10_mgonly_L1LATE-MU10_J50', 'HLT_mu10_mgonly_L1LATE-MU10_XE40', 'HLT_2g25_loose_g15_loose', 'HLT_mu11_mu6_bJpsimumu_L1LFV-MU11', 'HLT_mu11_mu6_bBmumux_BpmumuKp_L1LFV-MU11', 'HLT_mu11_mu6_bJpsimumu_Lxy0_L1LFV-MU11', 'HLT_mu11_mu6_bBmumux_BpmumuKp', 'HLT_mu11_mu6_bJpsimumu_Lxy0', 'HLT_mu11_mu6_bJpsimumu', 'HLT_mu11_mu6_bUpsimumu_L1LFV-MU11', 'HLT_mu11_mu6_bUpsimumu', 'HLT_3mu4_bTau', 'HLT_mu11_mu6_bBmumu_L1LFV-MU11', 'HLT_3mu6_bDimu', 'HLT_2mu6_bBmumu_Lxy0_L1BPH-2M9-2MU6_BPH-2DR15-2MU6', 'HLT_mu11_mu6_bBmumu', 'HLT_3mu6']
		if trigger == "muononly":
			self.allowed_trigger_list = apiSingleMuonTriggerlist
		elif trigger == "electrononly":
			self.allowed_trigger_list = apiSingleElectronTriggerlist
		elif trigger == "all":
			self.allowed_trigger_list = apiSingleMuonTriggerlist + apiSingleElectronTriggerlist
		elif trigger == "DAOD_RPVLL":
			self.allowed_trigger_list = self.DAOD_RPVLLTriggerlist
		else:
			self.allowed_trigger_list = list(trigger)
		# print self.allowed_trigger_list

	def passes(self):
		# evaluate whether the event passed the trigger
		# This method checks if there is any overlap between sets a and b
		# https://stackoverflow.com/questions/3170055/test-if-lists-share-any-items-in-python
		event_triggers = self.evt.tree.passedtriggers[self.evt.ievt]
		if self.invert: # invert trigger requirement
			# check if event_triggers includes a trigger in the DAOD_RPVLL trigger list
			if not set(event_triggers).isdisjoint(self.DAOD_RPVLLTriggerlist): 
				return set(event_triggers).isdisjoint(self.allowed_trigger_list)
			
			else: 
				return False
		else:	#default check if event_triggers includes a trigger on the allowed trigger list
			return not set(event_triggers).isdisjoint(self.allowed_trigger_list)



class Filter():
	def __init__(self, evt, filter_type):
		self.evt = evt
		self.filter_type = filter_type

	def passes(self):
		if self.filter_type == "mu-mu":
			return self.evt.tree.mumufilter[self.evt.ievt]

		if self.filter_type == "mu-el":
			return self.evt.tree.muelfilter[self.evt.ievt]

		if self.filter_type == "el-mu":
			return self.evt.tree.elmufilter[self.evt.ievt]

		if self.filter_type == "el-el":
			return self.evt.tree.elelfilter[self.evt.ievt]

		if self.filter_type == "4-filter":
			return (self.evt.tree.mumufilter[self.evt.ievt]
					or self.evt.tree.elmufilter[self.evt.ievt]
					or self.evt.tree.elelfilter[self.evt.ievt]
					or self.evt.tree.muelfilter[self.evt.ievt])

		if self.filter_type == "3-filter":
			return (self.evt.tree.mumufilter[self.evt.ievt]
					or self.evt.tree.elmufilter[self.evt.ievt]
					or self.evt.tree.elelfilter[self.evt.ievt])

		if self.filter_type == "2-filter":
			return (self.evt.tree.mumufilter[self.evt.ievt]
					or self.evt.tree.elmufilter[self.evt.ievt])

		if self.filter_type == "1-filter":
			return self.evt.tree.mumufilter[self.evt.ievt]


class PromptLepton():
	def __init__(self, evt, lepton="any", quality="tight", min_dR=0.05):
		self.evt = evt
		self.plepVec = ROOT.TLorentzVector(0,0,0,0)
		self.plepd0 = -2000
		self.plepz0 = -2000
		self.nPlep = 0
		quality = quality
		ndv = len(self.evt.tree.dvx[self.evt.ievt])	

		if lepton == "muon":
			if quality == "tight": #tight muon is requested
				lepquality = self.evt.tree.tightmu
			if quality == "medium":
				lepquality = self.evt.tree.mediummu
			if quality == "loose":
				lepquality = self.evt.tree.loosemu

			nleps = len(self.evt.tree.muonpt[self.evt.ievt])
			passPfilter = self.evt.tree.muonpassPfilter

		if lepton == "electron":
			if quality == "tight": #tight muon is requested
				lepquality = self.evt.tree.tightel
			if quality == "medium":
				lepquality = self.evt.tree.mediumel
			if quality == "loose":
				lepquality = self.evt.tree.looseel

			nleps = len(self.evt.tree.elpt[self.evt.ievt])
			passPfilter = self.evt.tree.elpassPfilter

		# initalize variable for prompt lepton				
		self.highestpt_lep = ROOT.TLorentzVector(0,0,0,0)
		self.highestpt_lep_d0 = -2000
		self.highestpt_lep_z0 = -2000


		for ilep in range(nleps): 
			overlap = False
			plepVec_i = ROOT.TLorentzVector()

			if lepton == "muon": 
				pt = self.evt.tree.muonpt[self.evt.ievt][ilep]
				eta = self.evt.tree.muoneta[self.evt.ievt][ilep]
				phi = self.evt.tree.muonphi[self.evt.ievt][ilep]
				mass = self.evt.tree.muonmass[self.evt.ievt][ilep]
				plepVec_i.SetPtEtaPhiM(pt,eta,phi,mass)

				lepd0 = self.evt.tree.muond0[self.evt.ievt][ilep]
				lepz0 = self.evt.tree.muonz0[self.evt.ievt][ilep]
				lepz0sintheta = self.evt.tree.muonz0sintheta[self.evt.ievt][ilep]

			if lepton == "electron":
				pt = self.evt.tree.elpt[self.evt.ievt][ilep]
				eta = self.evt.tree.eleta[self.evt.ievt][ilep]
				phi = self.evt.tree.elphi[self.evt.ievt][ilep]
				mass = self.evt.tree.elmass[self.evt.ievt][ilep]
				plepVec_i.SetPtEtaPhiM(pt,eta,phi,mass)

				lepd0 = self.evt.tree.eld0[self.evt.ievt][ilep]
				lepz0 = self.evt.tree.elz0[self.evt.ievt][ilep]
				lepz0sintheta = self.evt.tree.elz0sintheta[self.evt.ievt][ilep]

			# check if the plep passes the DRAW filter, quality and track significance cuts before looping over tracks in DV
			if passPfilter[self.evt.ievt][ilep] == True and lepquality[self.evt.ievt][ilep] == True and abs(lepd0) < 3 and abs(lepz0sintheta) < 0.5: # if lepton pass the track significance cuts : 
				for idv in range(ndv):
					leptracks = helpers.Tracks()
					trackevt = helpers.Event(self.evt.tree, self.evt.ievt, idv)
					leptracks.getTracks(trackevt)
					dlepVec = leptracks.lepVec
					ndtracks = len(dlepVec)
						
					for itr in range(ndtracks): # check overlap with DVs
						dR = dlepVec[itr].DeltaR(plepVec_i)
						if dR < min_dR:  # set overlap to true if muon overlaps with displaced track
							overlap = True
		
				if overlap == False: # if lepton doesnt overlap with and DV tracks
					self.nPlep = self.nPlep + 1 
					if (plepVec_i.Pt() > self.highestpt_lep.Pt()): # if pt is larger then the previous prompt lepton found 
						self.highestpt_lep = plepVec_i  #get highest pt prompt lepton!
						self.highestpt_lep_d0 = lepd0
						self.highestpt_lep_z0 = lepz0

					#for trigger matching
					# if self.evt.tree.muontrigmatched[self.evt.ievt][ilep] == 0:
					# 	print "is muon trig matched?", self.evt.tree.muontrigmatched[self.evt.ievt][ilep]
					# 	print "pt of the highest pt TIGHT muon: ", self.highestpt_plep.Pt() 
					# 	print "muon pt: ", self.evt.tree.muonpt[self.evt.ievt]
					# 	print "trigger matched: ", self.evt.tree.muontrigmatched[self.evt.ievt]
					# 	print "lepton quality: ", lepquality[self.evt.ievt]
					# 	print "muon type: ", self.evt.tree.muontype[self.evt.ievt]

	def passes(self):

		if self.nPlep > 0 and self.highestpt_lep.Pt() > 0: # check if you found a prompt lepton
			self.plepVec = self.highestpt_lep
			self.plepd0 = self.highestpt_lep_d0
			self.plepz0 = self.highestpt_lep_z0
			return True
		else: 
			return False

		

class InvertedPromptLepton():
	def __init__(self, evt, d0_cut=3.0, z0_sin_theta_cut=0.5):
		self.evt = evt
		self.n_prompt_leptons = 0 
		self.n_prompt_muons = 0 
		self.n_prompt_electrons = 0 

		n_muons = len(self.evt.tree.muonpt[self.evt.ievt])
		n_electrons = len(self.evt.tree.elpt[self.evt.ievt])
		self.n_prompt_muons = 0

		for imu in range(n_muons):
			mupt = self.evt.tree.muonpt[self.evt.ievt][imu]
			mud0 = self.evt.tree.muond0[self.evt.ievt][imu]
			mueta = self.evt.tree.muoneta[self.evt.ievt][imu]
			muphi = self.evt.tree.muonphi[self.evt.ievt][imu]
			mumass = self.evt.tree.muonmass[self.evt.ievt][imu]
			muz0sintheta = self.evt.tree.muonz0sintheta[self.evt.ievt][imu]
			mutype = self.evt.tree.muontype[self.evt.ievt][imu]
			muVec_i = ROOT.TLorentzVector()
			muVec_i.SetPtEtaPhiM(mupt,mueta,muphi,mumass)

			# check the muon has some quality (at least loose)
			if not ((self.evt.tree.tightmu[self.evt.ievt][imu] == 1) or
					(self.evt.tree.mediummu[self.evt.ievt][imu] == 1) or
					(self.evt.tree.loosemu[self.evt.ievt][imu]) == 1):
				# muon doesn't satisfy any quality, ignore it
				continue
			
			if abs(mud0) < d0_cut and abs(muz0sintheta) < z0_sin_theta_cut: # check muon passes prompt track cuts
				self.n_prompt_muons += 1 # count the number of prompt leptons 
				
		self.n_prompt_electrons = 0
		for iel in range(n_electrons): 
			elpt = self.evt.tree.elpt[self.evt.ievt][iel]
			eld0 = self.evt.tree.eld0[self.evt.ievt][iel]
			eleta = self.evt.tree.eleta[self.evt.ievt][iel]
			elphi = self.evt.tree.elphi[self.evt.ievt][iel]
			elmass = self.evt.tree.elmass[self.evt.ievt][iel]
			elz0sintheta = self.evt.tree.elz0sintheta[self.evt.ievt][iel]
			elVec_i = ROOT.TLorentzVector()
			elVec_i.SetPtEtaPhiM(elpt,eleta,elphi,elmass)

			# check the electron has some quality (at least loose)
			if not ((self.evt.tree.tightel[self.evt.ievt][iel] == 1) or
					(self.evt.tree.mediumel[self.evt.ievt][iel] == 1) or
					(self.evt.tree.looseel[self.evt.ievt][iel]) == 1):
				# electron doesn't satisfy any quality, ignore it
				continue
			
			if abs(eld0) < d0_cut and abs(elz0sintheta) < z0_sin_theta_cut: # check electron passes prompt track cuts
				# electron satisfies prompt lepton requirements
				self.n_prompt_electrons += 1 # count the number of prompt leptons 

		self.n_prompt_leptons = self.n_prompt_electrons + self.n_prompt_muons

	def passes(self):
		return self.n_prompt_leptons == 0





class nDV():
	def __init__(self, evt):
		self.evt = evt

	def passes(self):
		return len(self.evt.tree.dvx[self.evt.ievt]) > 0



class DVradius():
	def __init__(self, evt):
		self.evt = evt

		self.rdv = -1
		self.ntracks = len(self.evt.tree.dvx[self.evt.ievt])

		if self.ntracks  > 0:
			dx = self.evt.tree.dvx[self.evt.ievt][self.evt.idv]
			dy = self.evt.tree.dvy[self.evt.ievt][self.evt.idv]
			self.rdv = np.sqrt(dx**2 + dy**2)



	def passes(self, _min = 4,_max = 300):
		if (self.rdv > _min and self.rdv < _max):
			return True
		else: 
			return False




class DVntracks():
	def __init__(self, evt, ntrk=2,  decaymode="leptonic"):
		self.evt = evt
		self.ntrk = ntrk
		self.decaymode = decaymode

		self.ntracks = -1 

		if self.decaymode == "leptonic":
			self.ntracks = len(self.evt.tree.trackpt[self.evt.ievt][self.evt.idv])

	def passes(self):
		if self.ntracks == self.ntrk: 			
			return True
		else: 
			return False 




class ChargeDV(): 
	def __init__(self, evt, sel="OS",decaymode="leptonic"): 
		self.evt = evt
		self.decaymode = decaymode
		self.sel = sel
		self.ntracks = -1 
		self.charge_trk1 = -2 # dont make default -1 since that's a valid charge! :)
		self.charge_trk2 = -2 # dont make default -1 since that's a valid charge! :)

		if self.decaymode == "leptonic":
			self.ntracks = len(self.evt.tree.trackpt[self.evt.ievt][self.evt.idv])

			if self.ntracks == 2: 
				self.charge_trk1 = self.evt.tree.trackcharge[self.evt.ievt][self.evt.idv][0]
				self.charge_trk2 = self.evt.tree.trackcharge[self.evt.ievt][self.evt.idv][1]

	def passes(self): 
		if self.sel == 'OS':
			if self.charge_trk1 != self.charge_trk2: 
				return True
				
		elif self.sel == 'SS':
			if self.charge_trk1 == self.charge_trk2: 
				return True
		else:
			return False



class DVtype():
	def __init__(self, evt, dv_type, decaymode="leptonic"):
		self.evt = evt
		self.decaymode = decaymode
		self.dv_type = dv_type
		
		if self.decaymode == "leptonic": 
			self.ntracks = len(self.evt.tree.trackpt[self.evt.ievt][self.evt.idv])
			self.nel = -1
			self.nmu = -1


			self.muons = helpers.Tracks()
			self.muons.getMuons(self.evt)
			self.nmu = len(self.muons.lepVec)
			

			self.electrons = helpers.Tracks()
			self.electrons.getElectrons(self.evt)
			self.nel = len(self.electrons.lepVec)	


		


	def passes(self): 
		combined = 0 

		if self.dv_type == "emu": 
			if self.nel == 1 and self.nmu == 1: 
				mu1_type = self.evt.tree.muontype[self.evt.ievt][self.muons.lepIndex[0]]

				if mu1_type == combined:  # Only count combined muons 
					return True
				else:
					return False
			else:
				return False


		elif self.dv_type == "mumu":
			if self.nmu == 2: 
				mu1_type = self.evt.tree.muontype[self.evt.ievt][self.muons.lepIndex[0]]
				mu2_type = self.evt.tree.muontype[self.evt.ievt][self.muons.lepIndex[1]]

				if (mu1_type == combined and mu2_type == combined) :  # Only count combined muons 
					return True
				else:
					return False
			else:
				return False

		elif self.dv_type == "ee":
			if self.nel == 2: 
				return True
			else:
				return False

		elif self.dv_type == "mumu-notcomb":
			if self.nmu == 2: 

				return True
			else: 
				return False

		elif self.dv_type == "1-lep":
			if self.nmu > 0 or self.nel> 0: 					
				return True
			else: 
				return False
		elif self.dv_type == "2-lep":
			if self.nmu == 2 or (self.nmu == 1 and self.nel == 1) or self.nel ==2:
				return True
			else:
				return False






class Trackqual():
	def __init__(self,evt, decaymode="leptonic", quality="2-tight"):
		self.evt = evt
		self.decaymode = decaymode
		self.quality = quality 

		if self.decaymode == "leptonic": 
			self.ntracks = len(self.evt.tree.trackpt[self.evt.ievt][self.evt.idv])
			self.ntight = 0
			self.DVmuons = []
			self.DVelectrons = []

			muons = helpers.Tracks()
			muons.getMuons(self.evt)
		
			

			electrons = helpers.Tracks()
			electrons.getElectrons(self.evt)
			self.nel = len(electrons.lepVec)

			self.nmu_tight = 0
			self.nel_tight = 0

			self.ndvmu = len(muons.lepVec)
			self.ndvel = len(electrons.lepVec)
		
			for imu in range(self.ndvmu):
				muindex = muons.lepIndex[imu]
				muisTight = self.evt.tree.tightmu[self.evt.ievt][muindex]
				if muisTight: 
					self.nmu_tight = self.nmu_tight + 1

			for iel in range(self.ndvel):
				elindex = electrons.lepIndex[iel]
				elisTight = self.evt.tree.tightel[self.evt.ievt][elindex]
				if elisTight: 
					self.nel_tight = self.nel_tight + 1


	def passes(self):
			if self.quality == "2-tight": 
				# print self.nmu_tight
				if (self.nmu_tight == 2 or self.nel_tight == 2 or (self.nmu_tight == 1 and self.nel_tight == 1) ):
					return True
				else:
					return False

			if self.quality == "1-tight":	 
				if (self.nmu_tight > 0 or self.nel_tight > 0):
					return True
				else: 
					return False

class Cosmicveto():
	def __init__(self, evt, decaymode="leptonic",cosmicvetocut=0.05 ):
		self.evt = evt
		self.decaymode = decaymode
		self.cosmicvetocut = cosmicvetocut

		self.separation = -1 

		if self.decaymode == "leptonic": 
			ntracks = len(self.evt.tree.trackpt[self.evt.ievt][self.evt.idv])
			if ntracks == 2: 

				sumeta = self.evt.tree.tracketa[self.evt.ievt][self.evt.idv][0] + self.evt.tree.tracketa[self.evt.ievt][self.evt.idv][1]
				dphi = abs(self.evt.tree.trackphi[self.evt.ievt][self.evt.idv][0] - self.evt.tree.trackphi[self.evt.ievt][self.evt.idv][1])

				self.separation = np.sqrt( sumeta**2 + (np.pi -dphi)**2 )


	def passes(self):		
		if (self.separation > self.cosmicvetocut):
			return True
		else:
			return False
		


class Mlll():
	def __init__(self, dv_type, plep, dMu, dEl, decaymode="leptonic", _minmlll= 50 , _maxmlll = 84):
		self.decaymode = decaymode
		self.dv_type = dv_type
		self.plep = plep
		self.dMu = dMu
		self.dEl = dEl
		self._minmlll = _minmlll
		self._maxmlll = _maxmlll

		self.mlll = -1
		self.plll = ROOT.TLorentzVector(0,0,0,0)

		if self.decaymode == "leptonic":	
		
			if self.dv_type == "emu": 
				self.plll = self.plep + self.dEl[0] + self.dMu[0]	

			if self.dv_type == "mumu": 
				self.plll = self.plep + self.dMu[0] + self.dMu[1]
				

			if self.dv_type == "ee": 
				self.plll = self.plep + self.dEl[0] + self.dEl[1]
			
			self.mlll = self.plll.M()

	def passes(self):
		
		if (self.mlll> self._minmlll and self.mlll < self._maxmlll):
			return True
		else: 
			return False


class Mltt():
	
	def __init__(self, plep, trks, decaymode="leptonic", _minmltt= 50 , _maxmltt = 84):
		self.decaymode = decaymode
		self.plep = plep
		self.trks = trks
		self._minmltt = _minmltt
		self._maxmltt = _maxmltt

		self.mltt = -1
		self.plll = ROOT.TLorentzVector(0,0,0,0)

		if self.decaymode == "leptonic":	
	
			self.plll = self.plep + self.trks[0] + self.trks[1]
			self.mltt = self.plll.M()

	def passes(self):
		if (self.mltt> self._minmltt and self.mltt < self._maxmltt):
			return True
		else: 
			return False


class Mtrans():
	def __init__(self, plep, trks, decaymode="leptonic", _minmtrans= 50 , _maxmtrans = 84):
		self.decaymode = decaymode
		self.plep = plep
		self.trks = trks
		self._minmtrans = _minmtrans
		self._maxmtrans = _maxmtrans

		self.mtrans = -1
		self.mvis = -1 
		self.plll = ROOT.TLorentzVector(0,0,0,0)

		if self.decaymode == "leptonic":	
			self.plll = self.plep + self.trks[0] + self.trks[1]
			self.mvis = self.plll.M()
			self.mtrans = self.plll.Perp()

	def passes(self):
		
		if (self.mtrans> self._minmtrans and self.mtrans < self._maxmtrans):
			return True
		else: 
			return False



class DVmass():
	def __init__(self, evt, decaymode="leptonic",dvmasscut=4):
		self.evt = evt
		self.decaymode = decaymode
		self.dvmasscut = dvmasscut
		self.dvmass = self.evt.tree.dvmass[self.evt.ievt][self.evt.idv]

	def passes(self):
		if (self.dvmass > self.dvmasscut):
			return True
		else: 
			return False


class Mhnl():
	def __init__(self, evt, plep, trks,hnlmasscut=4):
		self.evt = evt
		self.plep = plep
		self.trks = trks
		self.hnlmasscut = hnlmasscut

		self.mvis = -1 
		self.mhnl = -1
		self.mhnl2 = -1
		self.hnlpt = -1
		self.hnleta = -99
		self.hnlphi = -99
		self.mtrans_rot = -1 
		self.neg_mhnl12 = -1
		self.neg_mhnl13 = -1
		self.neg_mhnl23 = -1
		self.pos_mhnl12 = -1
		self.pos_mhnl13 = -1
		self.pos_mhnl23 = -1

		self.nu_vec = ROOT.TLorentzVector()
		self.nu_vec2 = ROOT.TLorentzVector()

		def rotate_vector(r,v):
			r_new = ROOT.TVector3(v)
			rotation_axis = ROOT.TVector3(-1*r.Y(),r.X(),0.0)
			rotation_angle = -1*r.Theta() 
			r_new.Rotate(rotation_angle,rotation_axis)
			
			if (r== v and r_new.X() > 0.001): 
				#if r=v then you should end up with a vector all in the z component
				print logger.ERROR("Roatating vectors did not work!! Check HNL mass calculation.")
			return r_new

		def unrotate_vector(r,v):
			r_new = v

			rotation_axis = ROOT.TVector3(-1*r.Y(),r.X(),0.0)
			rotation_angle = r.Theta() 
			r_new.Rotate(rotation_angle,rotation_axis)
		
			return r_new

		#primary vertex vector
		pv_vec = ROOT.TVector3( self.evt.tree.pvx[self.evt.ievt],
								self.evt.tree.pvy[self.evt.ievt],
								self.evt.tree.pvz[self.evt.ievt]) 
		#displaced vertex vector
		dv_vec = ROOT.TVector3( self.evt.tree.dvx[self.evt.ievt][self.evt.idv],
							    self.evt.tree.dvy[self.evt.ievt][self.evt.idv],
								self.evt.tree.dvz[self.evt.ievt][self.evt.idv])

		# vector defining direction hnl trajectory
		hnl_vec =  dv_vec- pv_vec
		
		lepp_vec = ROOT.TVector3(plep.Px(),plep.Py(),plep.Pz()) 
		trkp_vec = []
		ntrk = len(self.trks)
		for i in xrange(ntrk):
			trkp_vec.append( ROOT.TVector3(self.trks[i].Px(),self.trks[i].Py(),self.trks[i].Pz()) )

		# print "Original DV vector: (", hnl_vec.X(),",", hnl_vec.Y(),",",hnl_vec.Z() , ")"

		#rotate coordinate system so hnl vector = z-axis
		lepp_vec_rot = rotate_vector(hnl_vec,lepp_vec)
		hnl_vec_rot = rotate_vector(hnl_vec,hnl_vec)
		trkp_vec_rot = []
		for i in xrange(ntrk):
			trkp_vec_rot.append(rotate_vector(hnl_vec,trkp_vec[i]))

		# check hnl_vec_rot is all in z component (x & y maybe be very small (~10^-15) due to precision when rotating vector)
		# print "New vector: (", hnl_vec_rot.X(),",", hnl_vec_rot.Y(),",",hnl_vec_rot.Z() , ")"

		# neutrino x-y momentum in rotated plane equals -1* sum of visible x-y momentum 
		pnu_rot = ROOT.TLorentzVector() 
		pnu_rot_x = -1*trkp_vec_rot[0].Px()-trkp_vec_rot[1].Px()
		pnu_rot_y = -1*trkp_vec_rot[0].Py()-trkp_vec_rot[1].Py()

		pnu_rot.SetPx(pnu_rot_x)
		pnu_rot.SetPy(pnu_rot_y)

		#make 4 vectors for visible particles in the rotated coordinate system (assume m = 0)
		plep_rot = ROOT.TLorentzVector(lepp_vec_rot,lepp_vec_rot.Mag())
		ptk1_rot = ROOT.TLorentzVector(trkp_vec_rot[0],trkp_vec_rot[0].Mag())
		ptk2_rot = ROOT.TLorentzVector(trkp_vec_rot[1],trkp_vec_rot[1].Mag())

		# make 4-vector for the visible particles in rotated coordinate system
		pvis_rot = plep_rot + ptk1_rot + ptk2_rot
		self.mvis = pvis_rot.M() # visible mass of the system
		self.mtrans_rot = pvis_rot.Perp() 

		# solve pw = pvis + pnu to find pnu_z
		m_w = 80.379 # mass of W boson in GeV
		K = (m_w**2 - self.mvis**2)/2 + pvis_rot.Px()*pnu_rot.Px() + pvis_rot.Py()*pnu_rot.Py()
		A = pvis_rot.Pz()**2 - pvis_rot.E()**2
		B = 2*K*pvis_rot.Pz()
		C = K**2 - pvis_rot.E()**2*(pnu_rot.Px()**2+ pnu_rot.Py()**2)

		# two solutions to the quadratic equation
		if (B**2 - 4*A*C) < 0: 
			#pw != pvis+ pnu so we can't solve quadratic
			return
		else: 
			pnu_rot_z1= (-B + np.sqrt(B**2 - 4*A*C) )/(2*A) 
			pnu_rot_z2= (-B - np.sqrt(B**2 - 4*A*C) )/(2*A)

			# two solutions for neutrino momentum vector
			pnu_rot1_vec = ROOT.TVector3(pnu_rot_x,pnu_rot_y,pnu_rot_z1) 
			pnu_rot2_vec = ROOT.TVector3(pnu_rot_x,pnu_rot_y,pnu_rot_z2) 

			# make a 4-vector for the neutrino
			pnu_rot1 = ROOT.TLorentzVector(pnu_rot1_vec,pnu_rot1_vec.Mag())
			pnu_rot2 = ROOT.TLorentzVector(pnu_rot2_vec,pnu_rot2_vec.Mag())

			# 2 solutions for HNL 4-vector
			pHNL_1 = ptk1_rot + ptk2_rot + pnu_rot1 # postive root
			pHNL_2 = ptk1_rot + ptk2_rot + pnu_rot2 # negative root
						#1        #2          #3
			
			neg_mhnl12_vec = ptk1_rot + ptk2_rot
			neg_mhnl13_vec = ptk1_rot + pnu_rot2
			neg_mhnl23_vec = ptk2_rot + pnu_rot2


			pos_mhnl12_vec = ptk1_rot + ptk2_rot
			pos_mhnl13_vec = ptk1_rot + pnu_rot1
			pos_mhnl23_vec = ptk2_rot + pnu_rot1



			# truth studies show pHNL_2 is the solution that gets us the HNL mass
			self.mhnl2 = pHNL_1.M()
			self.mhnl = pHNL_2.M()

			pHNL_2_lab = unrotate_vector(hnl_vec,pHNL_2)
			hnl_vec_unrot = unrotate_vector(hnl_vec,hnl_vec_rot)

			nu_vec_unrot2 = unrotate_vector(hnl_vec,pnu_rot1)

			nu_vec_unrot = unrotate_vector(hnl_vec,pnu_rot2)

			self.nu_vec = nu_vec_unrot
			self.nu_vec2 = nu_vec_unrot2

			# check that you rotated back to the original reference frame
			#print "Unrotated DV vector: (", hnl_vec_unrot.X(),",", hnl_vec_unrot.Y(),",",hnl_vec_unrot.Z() , ")"

			self.hnlpt = pHNL_2_lab.Pt()
			self.hnleta = pHNL_2_lab.Eta()
			self.hnlphi = pHNL_2_lab.Phi()


			self.neg_mhnl12 = neg_mhnl12_vec.M()
			self.neg_mhnl13 = neg_mhnl13_vec.M()
			self.neg_mhnl23 = neg_mhnl23_vec.M()

			self.pos_mhnl12 = pos_mhnl12_vec.M()
			self.pos_mhnl13 = pos_mhnl13_vec.M()
			self.pos_mhnl23 = pos_mhnl23_vec.M()

	def passes(self):
		
		if (self.mhnl > self.hnlmasscut ):
			return True
		else: 
			return False



class PV():
	def __init__(self, evt):
		self.evt = evt
		self.pv_x = self.evt.tree.pvx[self.evt.ievt]
		self.pv_y = self.evt.tree.pvy[self.evt.ievt]
		self.pv_z = self.evt.tree.pvz[self.evt.ievt]

	def passes(self):
		
		if (self.pv_x != -999.0 and self.pv_y != -999.0 and self.pv_z != -999.0 ): 
			return True
		else: 
			return False # no primary vertex in the event





