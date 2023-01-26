import ROOT
import numpy as np
import os
import sys
import helpers
import scale_factors
import selections
import observables
import ntuples
import systematics

UNLOCKED = 0
SELECTION_LOCKED = 1
FILL_LOCKED = 2


class Analysis(object):

	def __init__(self, name, tree, vtx_container, selection_list, output_file, save_ntuples, weight_override=None):
		# set up logger for self
		self.logger = helpers.getLogger('dHNLAnalysis.analysis', level=helpers.logger_debug_level)
		# set up logger for helper module
		selections.logger.setLevel(helpers.logger_debug_level)

		self.name = name
		self.weight_override = weight_override
		self.sel = selection_list
		self.output_file = output_file
		self.fi = ROOT.TFile.Open(output_file, 'update')
		self.ch = vtx_container
		self.save_ntuples = save_ntuples
		self.hist_suffixes = [self.ch]
		self.h = {}
		self.micro_ntuples = {}
		self.tree = tree
		self._locked = UNLOCKED
		# create an instance of Observables to store histograms
		self.observables = observables.Observables()

		self.jetVariables = {}
		self.jetVariables['pt'] = ROOT.std.vector('float')()
		self.jetVariables['eta'] = ROOT.std.vector('float')()
		self.jetVariables['phi'] = ROOT.std.vector('float')()
		self.jetVariables['E'] = ROOT.std.vector('float')()
		self.jetVariables['DL1dv00'] = ROOT.std.vector('float')()
		self.jetVariables['DL1dv01'] = ROOT.std.vector('float')()
		self.jetVariables['GN1'] = ROOT.std.vector('float')()
		self.triggerWeight = -999.

		self.events_with_trig_match_plep = 0
		self.events_with_trig_match_dlep = 0
		self.events_with_trig_match_both_pdlep =0

		# Calculate MC event weights as "L * xsec / num. MC events" (One number per file)
		# Single flavour mixing weights
		self.weight = helpers.MCEventWeight(self.tree, mixing_type ="single-flavour").get_mc_event_weight()
		self.mc_event_weight_one_dirac_hnl_single_flavour = helpers.MCEventWeight(self.tree, mixing_type ="single-flavour", dirac_limit = True).get_mc_event_weight()
		self.mc_event_weight_one_majorana_hnl_single_flavour = helpers.MCEventWeight(self.tree, mixing_type ="single-flavour").get_mc_event_weight()
		# Quasi-dirac pair "Majorana limit" with IH or NH mixing
		self.mc_event_weight_majorana_limit_ih = helpers.MCEventWeight(self.tree, mixing_type ="IH").get_mc_event_weight()
		self.mc_event_weight_majorana_limit_nh = helpers.MCEventWeight(self.tree, mixing_type ="NH").get_mc_event_weight()
		# Quasi-dirac pair "Majorana limit" with IH or NH mixing
		self.mc_event_weight_dirac_limit_ih = helpers.MCEventWeight(self.tree, mixing_type ="IH", dirac_limit = True).get_mc_event_weight()
		self.mc_event_weight_dirac_limit_nh = helpers.MCEventWeight(self.tree, mixing_type ="NH", dirac_limit = True).get_mc_event_weight()

		if self.tree.channel == "uue" or self.tree.channel == "eeu":
			self.mc_event_weight_majorana_limit_ih_flip_e_and_mu = helpers.MCEventWeight(self.tree, mixing_type ="IH", flip_e_and_mu = True).get_mc_event_weight()
			self.mc_event_weight_majorana_limit_nh_flip_e_and_mu = helpers.MCEventWeight(self.tree, mixing_type ="NH", flip_e_and_mu = True).get_mc_event_weight()
			self.mc_event_weight_dirac_limit_ih_flip_e_and_mu = helpers.MCEventWeight(self.tree, mixing_type ="IH", flip_e_and_mu = True, dirac_limit = True).get_mc_event_weight()
			self.mc_event_weight_dirac_limit_nh_flip_e_and_mu = helpers.MCEventWeight(self.tree, mixing_type ="NH", flip_e_and_mu = True, dirac_limit = True).get_mc_event_weight()

		# prepare systematics
		if self.tree.tree_name == 'nominal':
			self.lepton_reco_sf = {'nominal': 1,
								   'MUON_EFF_RECO_SYS_LOWPT__1down': 1,
								   'MUON_EFF_RECO_SYS_LOWPT__1up': 1,
								   'MUON_EFF_RECO_SYS__1down': 1,
								   'MUON_EFF_RECO_SYS__1up': 1,
								   'EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR__1down': 1,
								   'EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR__1up': 1,
								   'EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR__1down': 1,
								   'EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR__1up': 1, }
			self.lepton_trig_sf = {'nominal': 1,
								   'MUON_EFF_TrigSystUncertainty__1down': 1,
								   'MUON_EFF_TrigSystUncertainty__1up': 1,
								   'EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down': 1,
								   'EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up': 1,
								   }
		else:
			self.lepton_reco_sf = {'nominal': 1}
			self.lepton_trig_sf = {'nominal': 1}

		# setting all the relevant variables for the cuts based on the input selections
		# trigger cut
		self.do_trigger_cut = False
		if 'alltriggers' in self.sel:
			self.trigger = 'alltriggers'
			self.do_trigger_cut = True
		if 'electrononly' in self.sel:
			self.trigger = 'electrononly'
			self.do_trigger_cut = True
		if 'muononly' in self.sel:
			self.trigger = 'muononly'
			self.do_trigger_cut = True
		if (not self.do_trigger_cut) and ('CR' not in self.sel):
			self.logger.warn('You did not specify a trigger configuration for this channel. Skipping trigger selection.')
			
		# filter cut
		self.do_filter_cut = False
		if '4-filter' in self.sel:
			self.filter_type = '4-filter'
		elif '3-filter' in self.sel:
			self.filter_type = '3-filter'
		elif '2-filter' in self.sel:
			self.filter_type = '2-filter'
		elif '1-filter' in self.sel:
			self.filter_type = '1-filter'
		elif 'mumu-filter' in self.sel:
			self.filter_type = 'mu-mu'
		elif 'elmu-filter' in self.sel:
			self.filter_type = 'el-mu'
		elif 'elel-filter' in self.sel:
			self.filter_type = 'el-el'
		elif 'muel-filter' in self.sel:
			self.filter_type = 'mu-el'
		else:
			if 'CR' not in self.sel:
				self.logger.warn('You did not specify a filter configuration for this channel. Skipping filter selection.')
			self.do_filter_cut = False

		# prompt lepton cut
		if 'pmuon' in self.sel:
			self.plep = 'muon'
			self.do_prompt_lepton_cut = True
		elif 'pelectron' in self.sel:
			self.plep = 'electron'
			self.do_prompt_lepton_cut = True
		else:
			if 'CR' not in self.sel:
				self.logger.warn('You did not specify a prompt lepton for this channel. Skipping prompt lepton selection.')
			self.do_prompt_lepton_cut = False
		
		if "medium_plep" in self.sel: 
			self.plep_quality = "medium"
		else: 
			self.plep_quality =  "tight"

		self.do_same_event_cut =  "SE" in self.sel
		self.do_different_event_cut =  "DE" in self.sel
		
		# mass window cut
		self.do_mass_window_cut = 'mass_window' in self.sel

		# nDV cut
		self.do_ndv_cut = ('nDV' in self.sel)
		if not self.do_ndv_cut: self.logger.warn('You did not add nDV cut. Skipping nDV selection.')

		# fiducial volume
		self.do_fidvol_cut = 'fidvol' in self.sel
		if not self.do_fidvol_cut:
			self.logger.warn('You did not add DV in fiducial volume cut. Skipping DV in fiducial volume selection.')

		# 2 (or more) track cut
		self.do_ntrk_cut = True
		if '2track' in self.sel:
			self.ntrk = 2
		elif '3track' in self.sel:
			self.ntrk = 3
		elif '4track' in self.sel:
			self.ntrk = 4
		else:
			self.logger.warn('You did not add an ntrack cut. Skipping ntrack selection.')
			self.do_ntrk_cut = False

		# Opposite sign children vertex cut
		self.do_opposite_sign_cut = 'OS' in self.sel
		# Same sign children vertex cut
		self.do_same_sign_cut = 'SS' in self.sel
		if not (self.do_opposite_sign_cut or self.do_same_sign_cut):
			self.logger.warn('You did not add an SS or OS track cut. Skipping SS/OS track selection.')

		# DV type
		self.do_dv_type_cut = True
		if 'mumu' in self.sel:
			self.dv_type = "mumu"
		elif 'emu' in self.sel:
			self.dv_type = "emu"
		elif 'ee' in self.sel:
			self.dv_type = "ee"
		elif '1-lep' in self.sel:
			self.dv_type = "1-lep"
		elif '2-lep' in self.sel:
			self.dv_type = "2-lep"
		else:
			self.logger.warn('You did not specify a DV type for this channel. Skipping DV type selection.')
			self.do_dv_type_cut = False

		# Track quality
		self.do_track_quality_cut = True
		if '1-tight' in self.sel:
			self.track_quality = '1-tight'
		elif '1-medium' in self.sel:
			self.track_quality = '1-medium'
		elif '1-loose' in self.sel:
			self.track_quality = '1-loose'
		elif '2-tight' in self.sel:
			self.track_quality = '2-tight'
		elif '2-medium' in self.sel:
			self.track_quality = '2-medium'
		elif '2-loose' in self.sel:
			self.track_quality = '2-loose'
		elif 'tight-medium' in self.sel:
			self.track_quality = 'tight-medium'
		elif 'tight-loose' in self.sel:
			self.track_quality = 'tight-loose'
		elif 'medium-loose' in self.sel:
			self.track_quality = 'medium-loose'
		elif 'tight-veryloose' in self.sel:
			self.track_quality = 'tight-veryloose'
		elif 'medium-veryloose' in self.sel:
			self.track_quality = 'medium-veryloose'
		elif 'loose-veryloose' in self.sel:
			self.track_quality = 'loose-veryloose'
		elif 'tight-veryveryloose' in self.sel:
			self.track_quality = 'tight-veryveryloose'
		elif 'medium-veryveryloose' in self.sel:
			self.track_quality = 'medium-veryveryloose'
		elif 'loose-veryveryloose' in self.sel:
			self.track_quality = 'loose-veryveryloose'
		elif 'any-loose' in self.sel:
			self.track_quality = 'any-loose'
		elif 'any-veryveryloose' in self.sel:
			self.track_quality = 'any-veryveryloose'
		elif '2-any' in self.sel:
			self.track_quality = '2-any'
		elif '2-veryveryloose' in self.sel:
			self.track_quality = '2-veryveryloose'
		else:
			if "CR" not in self.sel:
				self.logger.warn('You did not specify a DV track quality for this channel. Skipping DV track quality selection.')
			self.do_track_quality_cut = False

		# cosmic veto cut
		self.do_cosmic_veto_cut = 'cosmicveto' in self.sel
		if not self.do_cosmic_veto_cut and 'CR' not in self.sel:
			self.logger.warn('You did not add a cosmic veto cut for this channel. Skipping cosmic veto selection.')

		# tri-lepton mass cut
		self.do_trilepton_mass_cut = 'mlll' in self.sel
		if not self.do_trilepton_mass_cut  and ("CR" not in self.sel or "inverted_mlll" not in self.sel):
			self.logger.warn('You did not add a mlll cut for this channel. Skipping tri-lepton mass selection.')
		# material veto cut
		self.do_mat_veto_cut = "matveto" in self.sel

		# DV mass cut
		self.do_dv_mass_cut = 'DVmass' in self.sel
		if not self.do_dv_mass_cut and "CR" not in self.sel:
			self.logger.warn('You did not add a DVmass cut for this channel. Skipping displaced vertex mass selection.')

		# HNL mass cut
		self.do_HNL_mass_cut = 'HNLmass' in self.sel
		
		# alpha cut
		self.do_alpha_cut = 'alpha' in self.sel
		
		# lepton pT cut
		self.do_lep_pt_cut = 'lep_pt' in self.sel

		# HNL pT cut
		self.do_HNL_pt_cut = 'HNLpt' in self.sel

		if 'CR' in self.sel:  # DO NOT CHANGE THESE CUTS OR YOU MIGHT UNBLIND DATA!!!
			self.do_CR = True
			self.tree.fake_aod = False
			self.do_trigger_cut = False  # do not apply trigger cut
			self.do_invert_trigger_cut = False  # do not apply inverted trigger cut
			self.do_filter_cut = False  # do not apply filter cut
			self.do_prompt_lepton_cut = False  # do not apply prompt lepton cut
			self.do_invert_prompt_lepton_cut = True  # invert prompt lepton cut
			self.do_inverted_mlll_cut = False
			self.do_inverted_mhnl_cut = False
			if 'ptrack' in self.sel: 
				self.do_prompt_track_cut = True # apply prompt track cut
			else: 
				self.do_prompt_track_cut = False # DO NOT apply prompt track cut
			self.logger.info('You are setup up to look in the inverted prompt lepton control region!')
		elif "CR_BE" in self.sel: # if running on fakeAOD that already has CR cuts applied (be careful with this setting!!!!!)
			self.tree.fake_aod = True
			self.do_CR = True
			self.do_trigger_cut = False  # do not apply trigger cut
			self.do_invert_trigger_cut = False  # do not apply inverted trigger cut
			self.do_filter_cut = False  # do not apply filter cut
			self.do_prompt_lepton_cut = False  # do not apply prompt lepton cut
			self.do_invert_prompt_lepton_cut = False  # do not apply inverted prompt lepton cut
			self.do_prompt_track_cut = False # do not apply prompt track cut
			self.do_inverted_mlll_cut = False
			self.do_inverted_mhnl_cut = False
			self.logger.info('You are running on a fake AOD created from events in the inverted prompt lepton control region!')
		elif "BE" in self.sel: # if running on fakeAOD without any CR cuts applied
			if "realDAOD" in self.sel:
				self.logger.info('You are running the background analysis on a real AOD created from events in the signal region!')
				self.tree.fake_aod = False
				self.do_trigger_cut = False  # apply a trigger cut
			else:
				self.logger.info('You are running the background analysis on a fake AOD created from events in the signal region!')
				self.tree.fake_aod = True
				self.do_trigger_cut = False  # apply a trigger cut
			self.do_CR = False
			self.do_invert_trigger_cut = False  # do not apply inverted trigger cut
			self.do_filter_cut = False  # do not apply filter cut, not sure if we need this for BE -DT
			self.do_prompt_lepton_cut = False  # do not apply prompt lepton cut, dont apply this cut no prompt leptons in the fake DAODs yet... -DT
			self.do_invert_prompt_lepton_cut = False  # no fake leptons in DAODs... -DT
			self.do_prompt_track_cut = False # do not apply prompt track cut
			self.do_inverted_mlll_cut = False
			self.do_inverted_mhnl_cut = False
			
			if "OS" in self.sel: raise ValueError("This analysis is blinded! You cannot look at OS DV from data events!!") # another blinded check -DT
		elif "inverted_mlll" in self.sel: 
			self.logger.warning('You are looking at a validation region with a SR pre-selection and an inverted m_lll cut!')
			self.do_CR = False
			self.tree.fake_aod = False
			self.do_invert_prompt_lepton_cut = False
			self.do_invert_trigger_cut = False
			self.do_prompt_track_cut = False 
			self.do_inverted_mlll_cut = True
			self.do_inverted_mhnl_cut = False
			self.do_trilepton_mass_cut = False
			self.save_ntuples = "mvis" # only save ntuples after mlll selection is applied!
		elif "inverted_mhnl" in self.sel:
			self.logger.warning('You are looking at a validation region with a SR pre-selection and an inverted m_hnl cut!')
			self.do_CR = False
			self.tree.fake_aod = False
			self.do_invert_prompt_lepton_cut = False
			self.do_invert_trigger_cut = False
			self.do_prompt_track_cut = False
			self.do_inverted_mlll_cut = False
			self.do_inverted_mhnl_cut = True
			self.do_HNL_mass_cut = False
			self.do_trilepton_mass_cut = False
			self.save_ntuples = "mhnl" # only save ntuples after mlll selection is applied!
		else:
			self.do_CR = False
			self.tree.fake_aod = False
			self.do_invert_prompt_lepton_cut = False
			self.do_invert_trigger_cut = False
			self.do_prompt_track_cut = False
			self.do_inverted_mlll_cut = False
			self.do_inverted_mhnl_cut = False

		if self.do_CR: 
			self.do_trigger_matching_cut = False
			self.do_zmass_veto = False
			if self.do_opposite_sign_cut == True and self.do_same_event_cut == True:
				self.be_region = "RegionAprime"
			elif self.do_same_sign_cut == True and self.do_same_event_cut == True: 
				self.be_region = "RegionBprime"
			elif self.do_opposite_sign_cut == True and self.do_different_event_cut == True: 
				self.be_region = "RegionCprime"
			elif self.do_same_sign_cut == True and self.do_different_event_cut == True: 
				self.be_region = "RegionDprime"
		else: 
			self.do_trigger_matching_cut = True
			self.do_zmass_veto = True
			if self.do_opposite_sign_cut == True and self.do_same_event_cut == True:
				self.be_region = "RegionA"
				raise ValueError("This analysis is blinded! You cannot look at OS DV from data events in the prompt lepton region")
			elif self.do_same_sign_cut == True and self.do_same_event_cut == True: 
				self.be_region = "RegionB"
			elif self.do_opposite_sign_cut == True and self.do_different_event_cut == True: 
				self.be_region = "RegionC"
				raise ValueError("This analysis is blinded! You cannot look at OS DV from data events in the prompt lepton region")
			elif self.do_same_sign_cut == True and self.do_different_event_cut == True: 
				self.be_region = "RegionD"

		self.check_input_consistency()


	# Getter helper functions
	def get_dv(self, key):
		return self.tree['secVtxLoose_{}_{}'.format(self.ch, key)]

	def get(self, key):
		return self.tree[key]


	def fill_ntuple(self, selection, ntuple_name, variable,full_name="",weight=None):
		"""
		A helper function for filling micro-ntuples.
		:param selection: the step of selection the analysis it at. May be "None" in which case there will be no prefix.
		:param ntuple_name: base name of the ntuple. When saved, a prefix and suffix will be appended.
		:param variable: variable you want to fill the histogram with.
		:param mc_type: defines whether the vertex is lepton-number conserving (LNC) or violating (LNV)
		:param full_name: override the automatic naming of the ntuple.
		"""
		#directory = '{ch}/{selection}/'.format(ch=self.ch, selection=selection)

		if self.MCEventType.isLNC:
			mc_types = ["LNC_","LNC_plus_LNV_"]
		elif self.MCEventType.isLNV:
			mc_types = ["LNV_","LNC_plus_LNV_"]
		else: mc_types = [""]

		if not selection:
			raise ValueError("You must indicate a selection in order to store the ntuple. Use 'all' if no selection.")
		prefixes = [mc_type + selection for mc_type in mc_types]

		# Retrieve the ntuple for this selection. If it doesn't exist, create it.
		for prefix in prefixes:
			if prefix not in self.micro_ntuples:
				# temp name. not written
				self.micro_ntuples[prefix] = ntuples.Ntuples('ntuples_{}_{}'.format(prefix, self.ch))
			# The name of the ntuple
			if not full_name:
				full_name = ntuple_name
			self.micro_ntuples[prefix][full_name] = variable

	def check_input_consistency(self):
		if self.do_trilepton_mass_cut or self.do_HNL_mass_cut or self.do_HNL_pt_cut:
			if self.do_invert_prompt_lepton_cut:
				self.logger.error("You are looking in the CR without prompt leptons so you cannot cut on mll, HNLpt or HNLm!!")
				sys.exit(1)  # abort because of error

			if not self.do_prompt_lepton_cut:
				if self.do_trilepton_mass_cut or self.do_HNL_mass_cut or self.do_HNL_pt_cut:
					self.logger.warning("You cannot cut on mlll, HNLpt or HNLm without first selecting a prompt lepton. Apply a prompt lepton cut!")
					sys.exit(1)  # abort because of error

		if self.do_opposite_sign_cut and self.do_same_sign_cut:
			self.logger.error("These cuts are mutually exclusive. You will get zero events!")
			sys.exit(1)  # abort because of error

		if self.do_inverted_mlll_cut and self.do_trilepton_mass_cut: 
			self.logger.error("Cannot invert mlll cut and also do mlll cut. These cuts are mutually exclusive!")
			sys.exit(1)  # abort because of error


	def unlock(self):
		self._locked = UNLOCKED

	def write(self):
		"""
		Write ntuples and histograms to root file
		"""
		# ____________________________________________________________
		# Write ntuples
		# Move ROOT to base directory
		self.fi.cd()
		# self.fi.mkdir(self.tree.tree_name)
		# self.fi.cd(self.tree.tree_name)
		for key, ntuple in self.micro_ntuples.items():
			ntuple.write(self.tree.tree_name + '_' + self.ch + '_ntuples_' + key)
		# ____________________________________________________________
		# Write histograms
		self.observables.write_histograms(root_file=self.fi, tree_name=self.tree.tree_name)
		self.logger.info("Histograms written to {}".format(self.output_file))

		self.fi.Close()

	def end(self):

		self.write()

		# Clean up memory
		del self.observables.histogram_dict
		del self.micro_ntuples



	# Protected function to create the selection object and return its success
	# The selection object may be used to fill additional histograms (see _prompt_lepton_cut)
	def _pv_cut(self):
		pv_sel = selections.PV(self.tree)
		return pv_sel.passes()

	def _trigger_cut(self):
		trigger_sel = selections.Trigger(self.tree, trigger=self.trigger)
		return trigger_sel.passes()

	def _invert_trigger_cut(self):
		trigger_sel = selections.Trigger(self.tree, trigger=self.trigger, invert=True)
		return trigger_sel.passes()

	def _nmuon_cut(self):
		n_muons = len(self.tree['muon_pt'])
		ncomb_muon = 0
		for imu in range(n_muons):
			mutype = self.tree['muon_type'][imu]
			if mutype == 0:
				ncomb_muon += 1
		return ncomb_muon != 3

	def _filter_cut(self):
		filter_sel = selections.Filter(self.tree, filter_type=self.filter_type)
		return filter_sel.passes()

	def _pass_prompt_lepton_overlap(self):
		prompt_lepton_overlap = selections.PromptLeptonOverlap(tree=self.tree, plep=self.plep, selected_plep = self.plep_sel)
		return prompt_lepton_overlap.passes()

	def _prompt_lepton_cut(self):
		self.found_plep = False # intitalize the plep each event 
		self.plep_sel = selections.PromptLepton(self.tree, lepton=self.plep,quality=self.plep_quality) # run plep selection 
		self.found_plep = self.plep_sel.found_plep # check if you found any prompt leptons 
		# Add to histogram all prompt leptons that pass selection.
		# If _prompt_lepton_cut() is run after trigger and filter cut then those cuts will also be applied.
		#if self.plep_sel.passes():
		#	# trig_match = selections.Lep_TriggerMatching(self.tree, self.plep, self.plep_sel.plep_Index)
		#	# if trig_match.lep_isTrigMatched:
		#	# 	self.events_with_trig_match_plep = self.events_with_trig_match_plep + 1
		#	
		#	# fill prompt lepton histograms
		#	self.fill_ntuple('all', 'prompt_lepton_pt', self.plep_sel.plepVec.Pt())
		#	self.fill_ntuple('all', 'prompt_lepton_eta', self.plep_sel.plepVec.Eta())
		#	self.fill_ntuple('all', 'prompt_lepton_phi', self.plep_sel.plepVec.Phi())
		#	self.fill_ntuple('all', 'prompt_lepton_d0', self.plep_sel.plep_d0)
		#	self.fill_ntuple('all', 'prompt_lepton_z0', self.plep_sel.plep_z0)
		return self.plep_sel.passes() # full plep selection find the highest pt plep that doesnt overlap with any DVs

	def _trigger_matched_medium_lepton_cut(self):
		"""
		Require that at least one lepton passes trigger matching.
		Also require that the lepton that passes trigger matching is of at least Medium quality.
		"""
		# get muons and electrons for this event
		muons = helpers.Muons(self.tree)
		electrons = helpers.Electrons(self.tree)
		self.trigger_matched_medium = selections.RequireMediumTriggerMatching(
			self.tree,
			prompt_lepton_index=self.plep_sel.plep_index,
			prompt_lepton_type=self.plep,
			muons=muons,
			electrons=electrons,
			dv_type=self.dv_type)

		return self.trigger_matched_medium.passes()

	def _prompt_track_cut(self):
		self.found_ptrk = False # initialize the plep each event
		self.ptrk_sel = selections.PromptTrack(self.tree) # run ptrack selection
		self.found_ptrk = self.ptrk_sel.found_trk # check if you found any prompt leptons
		# Add to histogram all prompt leptons that pass selection.
		# If _prompt_lepton_cut() is run after trigger and filter cut then those cuts will also be applied.
		if self.ptrk_sel.passes():
			self.fill_ntuple('all', 'ptrk_pt', self.ptrk_sel.trkVec.Pt())
			self.fill_ntuple('all', 'ptrk_eta', self.ptrk_sel.trkVec.Eta())
			self.fill_ntuple('all', 'ptrk_phi', self.ptrk_sel.trkVec.Phi())
			self.fill_ntuple('all', 'ptrk_d0', self.ptrk_sel.trkd0)
			self.fill_ntuple('all', 'ptrk_z0', self.ptrk_sel.trkz0)


		return self.ptrk_sel.passes() # full plep selection find the highest pt plep that doesnt overlap with any DVs

	def _invert_prompt_lepton_cut(self):
		self.invt_lep = selections.InvertedPromptLepton(self.tree)
		return self.invt_lep.passes()

	def _ndv_cut(self):
		return self.tree.ndv > 0

	def _fidvol_cut(self):
		fidvol_sel = selections.DVRadius(self.tree)
		return fidvol_sel.passes()

	def _ntrk_cut(self):
		ntracks_sel = selections.DVNTracks(self.tree, ntrk=self.ntrk)
		return ntracks_sel.passes()

	def _charge_cut(self):
		sign_pair = "SS" if self.do_same_sign_cut else "OS"
		charge_sel = selections.ChargeDV(self.tree, sel=sign_pair)
		return charge_sel.passes()
	
	def _be_event_type_cut(self): 
		if self.do_same_event_cut: return not self.tree.dv('shuffled')
		if self.do_different_event_cut: return self.tree.dv('shuffled')
		
	def _dv_type_cut(self):
		dv_sel = selections.DVType(self.tree, dv_type=self.dv_type)
		# if dv_sel.passes():
		# 	if not self.tree.fake_aod:
		# 		trig_match = selections.TriggerMatching_disp(self.tree, self.dv_type, dv_sel.dMu_Index, dv_sel.dEl_Index)
		# 		if trig_match.dlep_isTrigMatched:
		# 			self.events_with_trig_match_dlep = self.events_with_trig_match_dlep + 1
		# 		count_trig_match_disp_event = True

		return dv_sel.passes()

	def _track_quality_cut(self):
		track_quality_sel = selections.TrackQuality(self.tree, quality=self.track_quality)
		return track_quality_sel.passes()

	def _cosmic_veto_cut(self):
		cosmic_veto_sel = selections.CosmicVeto(self.tree)
		# self.h["DV_trk_sep"][self.ch].Fill(cosmic_veto_sel.separation)
		return cosmic_veto_sel.passes()

	def _trilepton_mass_cut(self):
		plep_vec = self.plep_sel.plepVec

		muons = helpers.Tracks(self.tree)
		muons.get_muons()
		muVec = muons.lepVec

		electrons = helpers.Tracks(self.tree)
		electrons.get_electrons()
		elVec = electrons.lepVec

		mlll_sel = selections.Mlll(dv_type=self.dv_type, plep=plep_vec, dMu=muVec, dEl=elVec,minmlll=40,maxmlll=90)
		return mlll_sel.passes()
	
	def _invert_trilepton_mass_cut(self):
		plep_vec = self.plep_sel.plepVec

		muons = helpers.Tracks(self.tree)
		muons.get_muons()
		muVec = muons.lepVec

		electrons = helpers.Tracks(self.tree)
		electrons.get_electrons()
		elVec = electrons.lepVec

		mlll_sel = selections.Mlll(dv_type=self.dv_type, plep=plep_vec, dMu=muVec, dEl=elVec,minmlll=25,maxmlll=125, invert=True)
		return mlll_sel.passes()

	def _mat_veto_cut(self):
		return selections.MaterialVeto(self.tree).passes()

	def _dv_mass_cut(self):
		dv_mass_sel = selections.DVMass(self.tree, dvmasscut=5.5)
		return dv_mass_sel.passes()
	
	def _Bhadron_veto(self):
		bhadron_sel = selections.BHadronVeto(self.tree, dv_type=self.dv_type, dv_mass_cut=5.5)
		return bhadron_sel.passes()

	def _Zmass_veto(self):
		zmass_veto = selections.ZMassVeto(self.tree, plep_vec = self.plep_sel.plepVec, plep=self.plep, plep_charge= self.plep_sel.plep_charge, dv_type= self.dv_type)
		return zmass_veto.passes()

	def _dv_lep_pt_cut(self):
		dv_lep_pt_sel = selections.DVLepPt(self.tree, self.dv_type)
		return dv_lep_pt_sel.passes()

	def _alpha_cut(self):
		alpha_sel = selections.Alpha(self.tree)
		return alpha_sel.passes()
	
	def _HNL_mass_cut(self):
		muons = helpers.Tracks(self.tree)
		muons.get_muons()
		muVec = muons.lepVec

		electrons = helpers.Tracks(self.tree)
		electrons.get_electrons()
		elVec = electrons.lepVec

		mHNL_sel = selections.Mhnl(self.tree, self.dv_type, plep=self.plep_sel.plepVec, dMu=muVec,dEl=elVec)
		return mHNL_sel.passes()

	def _invert_HNL_mass_cut(self):
		muons = helpers.Tracks(self.tree)
		muons.get_muons()
		muVec = muons.lepVec

		electrons = helpers.Tracks(self.tree)
		electrons.get_electrons()
		elVec = electrons.lepVec

		mHNL_sel = selections.Mhnl(self.tree, self.dv_type, plep=self.plep_sel.plepVec, dMu=muVec,dEl=elVec,hnlmasscut=20, invert= True)
		return mHNL_sel.passes()

	def _multitrk_2lep_cut(self):
		if self.tree.dv('ntrk') >= 2:  # 2+ trk vertex
			dv_type_sel = selections.DVType(self.tree, dv_type=self.dv_type)
			if dv_type_sel.passes():  # 2 leptons in the DV
				sign_pair = "SS" if self.do_same_sign_cut else "OS"
				charge_sel = selections.ChargeDV(self.tree, sel=sign_pair, trk_charge=dv_type_sel.lepton_charge)
				return charge_sel.passes()

	def _truth_match(self):
		"""
		Truth matching function for displaced vertices.
		linkTruth_score is calculated using DVAnalysisBase.
		pdgId 72 signifies a heavy neutral lepton parent particle.
		"""
		maxlinkTruth_score = self.tree.dv('maxlinkTruth_score')
		maxlinkTruth_parent_pdgId = abs(self.tree.dv('maxlinkTruth_parent_pdgId'))
		return self.tree.dv('maxlinkTruth_score') > 0.75 and abs(self.tree.dv('maxlinkTruth_parent_pdgId')) == 72

	def initialize_cut_bools(self):
		###########################################################################################################################
		# Cut bools that will be intialized in the pre selection for every event. These bools tell the code if the cutflow has already been filled for this event.
		# Default is to select the first event that passes the selection.
		###########################################################################################################################
		self.passed_preselection_cuts = False
		self.passed_fidvol_cut = False
		self.passed_ntrk_cut = False
		self.passed_charge_cut = False
		self.passed_dv_type_cut = False
		self.passed_dlep_pt_cut = False
		self.passed_be_event_type_cut = False
		self.passed_track_quality_cut = False
		self.passed_trig_matched_cut = False
		self.passed_cosmic_veto_cut = False
		self.passed_alpha_cut = False
		self.passed_trilepton_mass_cut = False
		self.passed_mat_veto_cut = False
		self.passed_dv_mass_cut = False
		self.passed_zmass_veto = False
		self.passed_HNL_mass_cut = False
		self.passed_HNL_pt_cut = False
		self.count_trig_match_disp_event = False

		self.selected_dv_index = -1


	def preSelection(self,use_truth=True):
		# if self.tree.max_entries == self.tree.ievt +1:
		# 	print "number of events with prompt lepton trigger matched: ", self.events_with_trig_match_plep
		# 	print "number of events with disp. lepton trigger matched: ", self.events_with_trig_match_dlep
		########################################
		# get the event weight for each event
		########################################
		self.calculate_event_weight()

		######################################################################################################
		# Preselection are all the cuts that are requied per event
		# Current cuts include: trigger, filter, plepton, DV cut
		######################################################################################################

		#initialize the cut bools for every event
		self.initialize_cut_bools()

		#if not self.tree.fake_aod:
		#	self._fill_leptons()

		if use_truth and not self.tree.is_data and not self.tree.is_bkg_mc:
			self._fill_truth_histos(sel='truth_all')

		self._fill_cutflow(0)

		######################################################################################################
		# Selection code is designed so that it will pass the selection only if the cut true or cut is unused
		# ex. passTrigger is true if the trigcut is true OR if trigcut is not used)
		######################################################################################################

		if self.do_trigger_cut:
			if self._trigger_cut():
				# Fill the cutflow plot at the specified bin
				self._fill_cutflow(1)
			else:
				return

		# Check to make sure event has a PV otherwise throw event away (this happens very rarely with data).
		if self._pv_cut():
			self._fill_cutflow(2)
		else:
			return

		# when running on data skip over any events without any DVs to speed up running
		if self.tree.is_data:
			if self.tree.ndv == 0: return

		if self.do_invert_trigger_cut:
			if self._invert_trigger_cut():
				self._fill_cutflow(2)
			else:
				return

		if self.do_filter_cut:
			if self._filter_cut():
				self._fill_cutflow(3)
			else:
				return

		if self.do_prompt_lepton_cut:
			plep_cut = self._prompt_lepton_cut()
			if self.found_plep:
				self._fill_cutflow(4)
			if plep_cut:
				if not self._pass_prompt_lepton_overlap():
					return
				self._fill_cutflow(5)
			else:
				return

		if self.do_invert_prompt_lepton_cut:
			if self._invert_prompt_lepton_cut():
				self._fill_cutflow(3)
			else:
				return

		if self.do_prompt_track_cut:
			ptrk_cut = self._prompt_track_cut()
			if self.found_ptrk:
				self._fill_cutflow(4)
			if ptrk_cut:
				self._fill_cutflow(5)
			else:
				return

		if self.do_ndv_cut:
			if self._ndv_cut():
				self._fill_cutflow(6)
			else:
				return

		# If you've made it here, preselection is passed
		self.passed_preselection_cuts = True
		#if use_truth and not self.tree.is_data and not self.tree.is_bkg_mc:
		#	self._fill_truth_histos(sel='truth/presel')

	def calculate_event_weight(self):
		"""
		Calculates and stores event weight only for event-level scale factors
		Do not use this for final signal weighting. These weights must be combined with
		the object-level scale factors.

		In general the weight is equal to (cross section weight * spin correlations weight)
		Due to the various models (single flavour mixing, inverted heiarchy (ih) and normal heiarchy (nh)),
		different weights are computed and saved to the mini-trees.

		These different weights can then be picked up by the limit setting framework to interpret the different models.
		"""
		# MC re-weighting to include spin correlations and fix lepton ordering bug (One number per event)
		self.MCEventType = selections.MCEventType(self.tree, mixing_type = "single-flavour")  # if data then MCEventType weight defaults to 1
		if self.weight_override is None:
			# #####################################################################################
			# Compute "model weight" as the product of mc event weight x spin corr weight
			# #####################################################################################
			weight_majorana_limit_ih = self.mc_event_weight_majorana_limit_ih * selections.MCEventType(self.tree, mixing_type = "IH").weight
			weight_majorana_limit_nh = self.mc_event_weight_majorana_limit_nh * selections.MCEventType(self.tree, mixing_type = "NH").weight
			weight_dirac_limit_ih = self.mc_event_weight_dirac_limit_ih * selections.MCEventType(self.tree, mixing_type = "IH").weight
			weight_dirac_limit_nh = self.mc_event_weight_dirac_limit_nh * selections.MCEventType(self.tree, mixing_type = "NH").weight

			# For uue and eeu channels, compute a second weight to reweight uue --> ueu  or eeu --> eue
			if self.tree.channel == "uue" or self.tree.channel == "eeu":	
				weight_majorana_limit_ih_2 = self.mc_event_weight_majorana_limit_ih_flip_e_and_mu * selections.MCEventType(self.tree, mixing_type = "IH",flip_e_and_mu = True).weight
				weight_majorana_limit_nh_2 = self.mc_event_weight_majorana_limit_nh_flip_e_and_mu * selections.MCEventType(self.tree, mixing_type = "NH",flip_e_and_mu = True).weight
				weight_dirac_limit_ih_2 = self.mc_event_weight_dirac_limit_ih_flip_e_and_mu * selections.MCEventType(self.tree, mixing_type = "IH",flip_e_and_mu = True).weight
				weight_dirac_limit_nh_2 = self.mc_event_weight_dirac_limit_nh_flip_e_and_mu * selections.MCEventType(self.tree, mixing_type = "NH",flip_e_and_mu = True).weight
				# Total weight is the sum of the two channels for models 
				self.model_weight_majorana_limit_ih = weight_majorana_limit_ih + weight_majorana_limit_ih_2
				self.model_weight_majorana_limit_nh = weight_majorana_limit_nh + weight_majorana_limit_nh_2
				self.model_weight_dirac_limit_ih = weight_dirac_limit_ih + weight_dirac_limit_ih_2
				self.model_weight_dirac_limit_nh = weight_dirac_limit_nh + weight_dirac_limit_nh_2
			# For all other channels the weight
			else:
				self.model_weight_majorana_limit_ih = weight_majorana_limit_ih
				self.model_weight_majorana_limit_nh = weight_majorana_limit_nh
				self.model_weight_dirac_limit_ih = weight_dirac_limit_ih
				self.model_weight_dirac_limit_nh = weight_dirac_limit_nh
			# Single flavour mixing models 
			self.model_weight_one_dirac_hnl_single_flavour = self.mc_event_weight_one_dirac_hnl_single_flavour * selections.MCEventType(self.tree, mixing_type = "single-flavour").weight
			self.model_weight_one_majorana_hnl_single_flavour = self.mc_event_weight_one_majorana_hnl_single_flavour * selections.MCEventType(self.tree, mixing_type = "single-flavour").weight
		else:
			self.model_weight_one_dirac_hnl_single_flavour = self.weight_override
			self.one_majorana_hnl_single_flavour = self.weight_override
			self.model_weight_diract_limit_ih = self.weight_override
			self.model_weight_diract_limit_nh = self.weight_override
			self.model_weight_majorana_limit_ih = self.weight_override
			self.model_weight_majorana_limit_nh = self.weight_override


	def DVSelection(self):
		raise NotImplementedError("Please implement this method in your own Analysis subclass")

	def _fill_cutflow(self, nbin):
		if not self.tree.is_data and not self.tree.is_bkg_mc:
			# store weighted cutflow with nominal scale factors
			scale_factor = self.lepton_reco_sf['nominal'] * self.lepton_trig_sf['nominal'] * self.tree['weight_pileup']
			if self.MCEventType.isLNC:
				self.CutFlow_LNC.Fill(nbin) # raw counts (only LNC events)
				self.CutFlow_LNC_weighted.Fill(nbin, self.model_weight_one_majorana_hnl_single_flavour * scale_factor)
				self.CutFlow_weighted_one_hnl_dirac.Fill(nbin, self.model_weight_one_dirac_hnl_single_flavour * scale_factor)
			if self.MCEventType.isLNV:
				self.CutFlow_LNV.Fill(nbin) # raw counts (only LNV events)
				self.CutFlow_LNV_weighted.Fill(nbin, self.model_weight_one_majorana_hnl_single_flavour * scale_factor)

			self.CutFlow.Fill(nbin) # raw counts (all events)
			# Models that count both LNC and LNV events!
			self.CutFlow_weighted_majorana_limit_ih.Fill(nbin, self.model_weight_majorana_limit_ih * scale_factor)
			self.CutFlow_weighted_majorana_limit_nh.Fill(nbin, self.model_weight_majorana_limit_nh * scale_factor)
			self.CutFlow_weighted_dirac_limit_ih.Fill(nbin, self.model_weight_dirac_limit_ih * scale_factor)
			self.CutFlow_weighted_dirac_limit_nh.Fill(nbin, self.model_weight_dirac_limit_nh * scale_factor)
			self.CutFlow_weighted_one_hnl_majorana.Fill(nbin, self.model_weight_one_majorana_hnl_single_flavour * scale_factor)
		else:
			self.CutFlow.Fill(nbin)
			self.CutFlow_weighted.Fill(nbin, self.weight)

	def _fill_leptons(self):
		sel = 'all'
		for imu in range(len(self.tree['muon_type'])):
			self.fill_ntuple(sel, 'muon_type', self.tree['muon_type'][imu])
			self.fill_ntuple(sel, 'muon_pt', self.tree['muon_pt'][imu])
			self.fill_ntuple(sel, 'muon_eta', self.tree['muon_eta'][imu])
			self.fill_ntuple(sel, 'muon_phi', self.tree['muon_phi'][imu])
			if self.tree['muon_isTight'][imu] == 1:  self.fill_ntuple(sel, 'muon_quality', 3)
			if self.tree['muon_isMedium'][imu] == 1: self.fill_ntuple(sel, 'muon_quality', 2)
			if self.tree['muon_isLoose'][imu] == 1:  self.fill_ntuple(sel, 'muon_quality', 1)
			else: self.fill_ntuple(sel, 'muon_quality', 0)
			self.fill_ntuple(sel, 'muon_isLRT', self.tree['muon_isLRT'][imu])

		for iel in range(len(self.tree['el_pt'])):
			self.fill_ntuple(sel, 'el_pt', self.tree['el_pt'][iel])
			self.fill_ntuple(sel, 'el_eta', self.tree['el_eta'][iel])
			self.fill_ntuple(sel, 'el_phi', self.tree['el_phi'][iel])
			if self.tree['el_LHTight'][iel] == 1:  self.fill_ntuple(sel, 'el_quality', 3)
			if self.tree['el_LHMedium'][iel] == 1: self.fill_ntuple(sel, 'el_quality', 2)
			if self.tree['el_LHLoose'][iel] == 1:  self.fill_ntuple(sel, 'el_quality', 1)
			else: self.fill_ntuple(sel, 'el_quality', 0)
			self.fill_ntuple(sel, 'el_isLRT', self.tree['el_isLRT'][imu])

	def _fill_all_dv_histos(self):
		sel = 'all'	
		self.fill_ntuple(sel, 'DV_num_trks', self.tree.dv('ntrk'))
		self.fill_ntuple(sel, 'DV_x', self.tree.dv('x'))
		self.fill_ntuple(sel, 'DV_y', self.tree.dv('y'))
		self.fill_ntuple(sel, 'DV_z', self.tree.dv('z'))
		self.fill_ntuple(sel, 'DV_r', self.tree.dv('r'))
		self.fill_ntuple(sel, 'DV_distFromPV', self.tree.dv('distFromPV'))
		self.fill_ntuple(sel, 'DV_mass', self.tree.dv('mass'))
		self.fill_ntuple(sel, 'DV_pt', self.tree.dv('pt'))
		self.fill_ntuple(sel, 'DV_eta', self.tree.dv('eta'))
		self.fill_ntuple(sel, 'DV_phi', self.tree.dv('phi'))
		self.fill_ntuple(sel, 'DV_minOpAng', self.tree.dv('minOpAng'))
		self.fill_ntuple(sel, 'DV_maxOpAng', self.tree.dv('maxOpAng'))
		self.fill_ntuple(sel, 'DV_charge', self.tree.dv('charge'))
		self.fill_ntuple(sel, 'DV_chi2', self.tree.dv('chi2'))
		# self.fill_ntuple(sel, 'DV_chi2_assoc', self.tree.dv('chi2_assoc'))
		
		# Write values to ntuple. TTree already created and array already filled by fill_ntuple
		if self.MCEventType.isLNC: 
			self.micro_ntuples["LNC_"+sel].fill()
			self.micro_ntuples["LNC_plus_LNV_"+sel].fill()
		elif self.MCEventType.isLNV: 
			self.micro_ntuples["LNV_"+sel].fill()
			self.micro_ntuples["LNC_plus_LNV_"+sel].fill()
		else: self.micro_ntuples[sel].fill()
			

	def _fill_truth_histos(self, sel):
		truth_info = helpers.Truth()
		truth_info.get_truth_particles(self.tree)
		
		# self.fill_ntuple(sel, 'event_type_MCweight', self.MCEventType.weight)  #if not weight_override else weight_override
		# self.fill_ntuple(sel, 'M2_spin_corr_MCweight', self.MCEventType.M2_spin_corr)  #if not weight_override else weight_override
		# self.fill_ntuple(sel, 'M2_nocorr_MCweight', self.MCEventType.M2_nocorr)  #if not weight_override else weight_override
		self.fill_ntuple(sel, 'W_pt', truth_info.W_pt)
		self.fill_ntuple(sel, 'W_eta', truth_info.W_eta)
		self.fill_ntuple(sel, 'W_phi', truth_info.W_phi)
		self.fill_ntuple(sel, 'W_mass', truth_info.W_mass)
		self.fill_ntuple(sel, 'pt_hnl', truth_info.N_pt)
		self.fill_ntuple(sel, 'eta_hnl', truth_info.N_eta)
		self.fill_ntuple(sel, 'phi_hnl', truth_info.N_phi)
		self.fill_ntuple(sel, 'mhnl', truth_info.mhnl)
		# self.fill_ntuple(sel, 'generated_mN', truth_info.mhnl)
		self.fill_ntuple(sel, 'DV_mass', truth_info.DV_mass)
		self.fill_ntuple(sel, 'mlll', truth_info.lll_mass)
		self.fill_ntuple(sel, 'DV_r', truth_info.DV_r)
		self.fill_ntuple(sel, 'DV_x', truth_info.DV_x)
		self.fill_ntuple(sel, 'DV_y', truth_info.DV_y)
		self.fill_ntuple(sel, 'DV_z', truth_info.DV_z)
		self.fill_ntuple(sel, 'prompt_lepton_pt', truth_info.plep_vec.Pt())
		self.fill_ntuple(sel, 'prompt_lepton_eta', truth_info.plep_vec.Eta())
		self.fill_ntuple(sel, 'prompt_lepton_phi', truth_info.plep_vec.Phi())
		self.fill_ntuple(sel, 'prompt_lepton_mass', truth_info.plep_vec.M())

		self.fill_ntuple(sel, 'displaced_lepton1_pt', truth_info.dv_track_1.Pt())
		self.fill_ntuple(sel, 'displaced_lepton1_eta', truth_info.dv_track_1.Eta())
		self.fill_ntuple(sel, 'displaced_lepton1_phi', truth_info.dv_track_1.Phi())
		self.fill_ntuple(sel, 'displaced_lepton1_mass', truth_info.dv_track_1.M())

		self.fill_ntuple(sel, 'displaced_lepton2_pt', truth_info.dv_track_2.Pt())
		self.fill_ntuple(sel, 'displaced_lepton2_eta', truth_info.dv_track_2.Eta())
		self.fill_ntuple(sel, 'displaced_lepton2_phi', truth_info.dv_track_2.Phi())
		self.fill_ntuple(sel, 'displaced_lepton2_mass', truth_info.dv_track_2.M())



		self.fill_ntuple(sel, 'properLifetime', truth_info.properLifetime)
		# if truth_info.W_charge == 1: 
		# 	self.fill_ntuple(sel, 'Wplus_HNLpt', truth_info.HNL_vec.Pt())
		# 	self.fill_ntuple(sel, 'Wplus_HNLeta', truth_info.HNL_vec.Eta())
		# 	self.fill_ntuple(sel, 'Wplus_HNLphi', truth_info.HNL_vec.Phi())
		# 	self.fill_ntuple(sel, 'Wplus_HNLE', truth_info.HNL_vec.E())
		# if truth_info.W_charge == -1: 
		# 	self.fill_ntuple(sel, 'Wminus_HNLpt', truth_info.HNL_vec.Pt())
		# 	self.fill_ntuple(sel, 'Wminus_HNLeta', truth_info.HNL_vec.Eta())
		# 	self.fill_ntuple(sel, 'Wminus_HNLphi', truth_info.HNL_vec.Phi())
		# 	self.fill_ntuple(sel, 'Wminus_HNLE', truth_info.HNL_vec.E())
		# if len(truth_info.trkVec) == 2: 
		# 	DV_4vec= truth_info.trkVec[1]+ truth_info.trkVec[0]
		# 	lep12 = truth_info.dLepVec[0] + truth_info.dLepVec[1] 
		# 	lep23 = truth_info.dLepVec[1] + truth_info.dLepVec[2] 
		# 	lep13 = truth_info.dLepVec[0] + truth_info.dLepVec[2] 
		# 	all_leptons = truth_info.plep_vec + truth_info.trkVec[0] + truth_info.trkVec[1]
		# 	self.fill_ntuple(sel, 'DV_mass', DV_4vec.M())
		# 	self.fill_ntuple(sel, 'm12', lep12.M())
		# 	self.fill_ntuple(sel, 'm23', lep23.M())
		# 	self.fill_ntuple(sel, 'm13', lep13.M())
		# 	self.fill_ntuple(sel, 'm12_sq', lep12.M()**2)
		# 	self.fill_ntuple(sel, 'm23_sq', lep23.M()**2)
		# 	self.fill_ntuple(sel, 'm13_sq', lep13.M()**2)
		# 	self.fill_ntuple(sel, 's12', self.MCEventType.s12) 
		# 	self.fill_ntuple(sel, 's13', self.MCEventType.s13) 
		# 	self.fill_ntuple(sel, 's14', self.MCEventType.s14) 
		# 	self.fill_ntuple(sel, 's23', self.MCEventType.s23) 
		# 	self.fill_ntuple(sel, 's24', self.MCEventType.s24) 
		# 	self.fill_ntuple(sel, 's34', self.MCEventType.s34) 
		# 	self.fill_ntuple(sel, 'lep1_trk_pt', self.MCEventType.p_2.Pt()) # topological ordered
		# 	self.fill_ntuple(sel, 'lep1_trk_eta', self.MCEventType.p_2.Eta())
		# 	self.fill_ntuple(sel, 'lep1_trk_phi', self.MCEventType.p_2.Phi())
		# 	self.fill_ntuple(sel, 'lep2_trk_pt', self.MCEventType.p_3.Pt())
		# 	self.fill_ntuple(sel, 'lep2_trk_eta', self.MCEventType.p_3.Eta())
		# 	self.fill_ntuple(sel, 'lep2_trk_phi', self.MCEventType.p_3.Phi())
		# 	self.fill_ntuple(sel, 'nu_trk_pt', self.MCEventType.p_4.Pt())
		# 	self.fill_ntuple(sel, 'nu_trk_eta', self.MCEventType.p_4.Eta())
		# 	self.fill_ntuple(sel, 'nu_trk_phi', self.MCEventType.p_4.Phi())
		# 	disp_lep = [self.MCEventType.p_2,self.MCEventType.p_3,self.MCEventType.p_4]
		# 	# pt order the displaced leptons
		# 	disp_lep.sort(key=lambda x: x.Pt(), reverse=True)
		# 	self.fill_ntuple(sel, 'dlep1_pt', disp_lep[0].Pt()) # pt ordered
		# 	self.fill_ntuple(sel, 'dlep1_eta', disp_lep[0].Eta())
		# 	self.fill_ntuple(sel, 'dlep1_phi', disp_lep[0].Phi())
		# 	self.fill_ntuple(sel, 'dlep2_pt', disp_lep[1].Pt())
		# 	self.fill_ntuple(sel, 'dlep2_eta', disp_lep[1].Eta())
		# 	self.fill_ntuple(sel, 'dlep2_phi', disp_lep[1].Phi())
		# 	self.fill_ntuple(sel, 'dlep3_pt', disp_lep[2].Pt())
		# 	self.fill_ntuple(sel, 'dlep3_eta', disp_lep[2].Eta())
		# 	self.fill_ntuple(sel, 'dlep3_phi', disp_lep[2].Phi())

		# 	if (abs(truth_info.dTrk_d0[0]) < 2 and abs(truth_info.dTrk_d0[1]) < 2):
		# 		self.fill_ntuple(sel, 'DV_d0_cut',1, fill_ntuple=False)
		# 	else:
		# 		self.fill_ntuple(sel, 'DV_d0_cut',0, fill_ntuple=False)

		# 	n_el = len(truth_info.dEl)
		# 	n_mu = len(truth_info.dMu)
		# 	for iel in range(n_el):
		# 		self.fill_ntuple(sel, 'DV_El_{}_pt'.format(iel), truth_info.dEl[iel].Pt(), fill_ntuple=False)
		# 		self.fill_ntuple(sel, 'DV_El_{}_eta'.format(iel), truth_info.dEl[iel].Eta(), fill_ntuple=False)
		# 		self.fill_ntuple(sel, 'DV_El_{}_phi'.format(iel), truth_info.dEl[iel].Phi(), fill_ntuple=False)
		# 		self.fill_ntuple(sel, 'DV_El_{}_d0'.format(iel), truth_info.dEl_d0[iel], fill_ntuple=False)
		# 		self.fill_ntuple(sel, 'DV_El_{}_charge'.format(iel), truth_info.dEl_charge[iel], fill_ntuple=False)
			
		# 	for imu in range(n_mu):
		# 		self.fill_ntuple(sel, 'DV_Mu_{}_pt'.format(imu), truth_info.dMu[imu].Pt(), fill_ntuple=False)
		# 		self.fill_ntuple(sel, 'DV_Mu_{}_eta'.format(imu), truth_info.dMu[imu].Eta(), fill_ntuple=False)
		# 		self.fill_ntuple(sel, 'DV_Mu_{}_phi'.format(imu), truth_info.dMu[imu].Phi(), fill_ntuple=False)
		# 		self.fill_ntuple(sel, 'DV_Mu_{}_d0'.format(imu), truth_info.dMu_d0[imu], fill_ntuple=False)
		# 		self.fill_ntuple(sel, 'DV_Ml_{}_charge'.format(imu), truth_info.dMu_charge[imu], fill_ntuple=False)

		# 	for itrk in range(2):
		# 		self.fill_ntuple(sel, 'DV_trk_pt', truth_info.trkVec[itrk].Pt(), fill_ntuple=False) # do the same here but also save charge for lepton truth matching later
		# 		self.fill_ntuple(sel, 'DV_trk_eta', truth_info.trkVec[itrk].Eta(), fill_ntuple=False)
		# 		self.fill_ntuple(sel, 'DV_trk_phi', truth_info.trkVec[itrk].Phi(), fill_ntuple=False)
		# 		self.fill_ntuple(sel, 'DV_trk_d0',truth_info.dTrk_d0[itrk], fill_ntuple=False)

		# Write values to ntuple. TTree already created and array already filled by fill_ntuple
		if self.MCEventType.isLNC:
			self.micro_ntuples["LNC_"+sel].fill()
			self.micro_ntuples["LNC_plus_LNV_"+sel].fill()
		elif self.MCEventType.isLNV:
			self.micro_ntuples["LNV_"+sel].fill()
			self.micro_ntuples["LNC_plus_LNV_"+sel].fill()
		else: self.micro_ntuples[sel].fill()

	def _fill_systematic_branches(self, sel):
		"""
		Fills the branches to be used for systematic variations.
		These systematic weights (including nominal) should be applied by multiplication with the existing event weight.
		This is done in the statistical interpretation framework.
		@param sel: selection step in cuts
		"""
		# fill nominal branch
		self.fill_ntuple(sel, 'SF_nominal', self.lepton_reco_sf['nominal'] * self.lepton_trig_sf['nominal'] * self.tree['weight_pileup'])
		# vary lepton reconstruction systematics, hold rest nominal
		for systematic in self.lepton_reco_sf.keys():
			if systematic == 'nominal': continue
			self.fill_ntuple(sel, 'SF_' + systematic, self.lepton_reco_sf[systematic] * self.lepton_trig_sf['nominal'] * self.tree['weight_pileup'])
		# vary lepton trigger systematics, hold rest nominal
		for systematic in self.lepton_trig_sf.keys():
			if systematic == 'nominal': continue
			self.fill_ntuple(sel, 'SF_' + systematic, self.lepton_reco_sf['nominal'] * self.lepton_trig_sf[systematic] * self.tree['weight_pileup'])
		# vary lepton trigger systematics, hold rest nominal
		for systematic in ['weight_pileup_up', 'weight_pileup_down']:
			if not self.tree.tree_name == 'nominal': continue # only do pileup systematics on nominal tree
			self.fill_ntuple(sel, 'SF_' + systematic, self.lepton_reco_sf['nominal'] * self.lepton_trig_sf['nominal'] * self.tree[systematic])


	def _fill_selected_dv_histos(self, sel, do_lock=True):
		if self._locked < FILL_LOCKED and do_lock:
			# These are the histograms you only want to fill ONCE per DV
			# ____________________________________________________________
			# Systematics
			#if self.save_ntuples == sel:
			self._fill_systematic_branches(sel)
			# ____________________________________________________________
			# Fill mc weights for the different HNL models
			# Generic weight to be used for MC bkg
			evt_weight = self.tree['mcEventWeight']
			self.fill_ntuple(sel, 'weight', self.weight * evt_weight)
			# Weights for LNC + LNV decays
			self.fill_ntuple(sel, 'model_weight_one_majorana_hnl_LNCplusLNV_single_flavour_mixing', self.model_weight_one_majorana_hnl_single_flavour * evt_weight)
			self.fill_ntuple(sel, 'model_weight_quasi_dirac_pair_LNCplusLNV_ih_mixing', self.model_weight_majorana_limit_ih * evt_weight)
			self.fill_ntuple(sel, 'model_weight_quasi_dirac_pair_LNCplusLNV_nh_mixing', self.model_weight_majorana_limit_nh * evt_weight)
			# LNC only weights
			self.fill_ntuple(sel, 'model_weight_one_dirac_hnl_LNC_single_flavour_mixing', self.model_weight_one_dirac_hnl_single_flavour * evt_weight)
			self.fill_ntuple(sel, 'model_weight_quasi_dirac_pair_LNC_ih_mixing', self.model_weight_dirac_limit_ih * evt_weight)
			self.fill_ntuple(sel, 'model_weight_quasi_dirac_pair_LNC_nh_mixing', self.model_weight_dirac_limit_nh * evt_weight)
			# ____________________________________________________________
			# Fill HNL cross sections for different models
			one_majorana_hnl_single_flavour_xsec = helpers.MCEventWeight(self.tree, mixing_type="single-flavour").hnl_xsec_generic_model(channel=self.tree.channel, mass=self.tree.mass, ctau=self.tree.ctau)
			ih_xsec = helpers.MCEventWeight(self.tree, mixing_type="IH").hnl_xsec_generic_model(channel=self.tree.channel, mass=self.tree.mass, ctau=self.tree.ctau)
			nh_xsec = helpers.MCEventWeight(self.tree, mixing_type="NH").hnl_xsec_generic_model(channel=self.tree.channel, mass=self.tree.mass, ctau=self.tree.ctau)
			self.fill_ntuple(sel, 'LNC_xsec_one_majorana_hnl_single_flavour', one_majorana_hnl_single_flavour_xsec)
			self.fill_ntuple(sel, 'LNC_xsec_one_dirac_hnl_single_flavour', one_majorana_hnl_single_flavour_xsec * 2)
			self.fill_ntuple(sel, 'NH_xsec', nh_xsec * 2)
			self.fill_ntuple(sel, 'IH_xsec', ih_xsec * 4)
			# ____________________________________________________________
			self.fill_ntuple(sel, 'event_is_LNC', self.MCEventType.isLNC)
			self.fill_ntuple(sel, 'event_is_LNV', self.MCEventType.isLNV)
			# ____________________________________________________________
			# add mc event weight
			self.fill_ntuple(sel, 'mcEventWeight', self.tree['mcEventWeight'])
			self.fill_ntuple(sel, 'runNumber',self.tree["runNumber"])
			if (self.tree["mcChannelNumber"] != self.tree.mcChannelNumber) and not self.tree.is_data and self.tree.is_bkg_mc:
				self.logger.error("DSID {} put in config does not match DSID {} read from the ntuple. Please check your configuration.".format(self.tree.mcChannelNumber,self.tree["mcChannelNumber"]))
				sys.exit(1)  # abort because of error
			self.fill_ntuple(sel, 'mcChannelNumber',self.tree["mcChannelNumber"])
			self.fill_ntuple(sel, 'eventNumber',self.tree['eventNumber'])
			self.AddExtraVariables(sel)

			# ____________________________________________________________
			tracks = helpers.Tracks(self.tree)
			tracks.get_tracks()

			muons = helpers.Muons(self.tree)
			mu_vec = muons.lepVec
			mu_index = muons.lepIndex

			electrons = helpers.Electrons(self.tree)
			el_vec = electrons.lepVec
			el_index = electrons.lepIndex

			# fill histograms that require a prompt lepton to be identified
			if self.do_prompt_lepton_cut:
				plep_vec = self.plep_sel.plepVec
				plep_d0 = self.plep_sel.plep_d0
				plep_z0 = self.plep_sel.plep_z0
				plep_charge = self.plep_sel.plep_charge
				plep_isTight = self.plep_sel.plep_isTight
				plep_index = self.plep_sel.plep_index
				plep_is_trigger_matched = selections.LeptonTriggerMatching(self.tree, self.plep, plep_index).is_trigger_matched

				self.fill_ntuple(sel, 'prompt_lepton_pt', plep_vec.Pt())
				self.fill_ntuple(sel, 'prompt_lepton_eta', plep_vec.Eta())
				self.fill_ntuple(sel, 'prompt_lepton_phi', plep_vec.Phi())
				self.fill_ntuple(sel, 'prompt_lepton_d0', plep_d0)
				self.fill_ntuple(sel, 'prompt_lepton_z0', plep_z0)
				self.fill_ntuple(sel, 'prompt_lepton_charge', plep_charge)
				self.fill_ntuple(sel, 'prompt_lepton_isTight', plep_isTight)
				self.fill_ntuple(sel, 'prompt_lepton_is_trigger_matched', plep_is_trigger_matched)

				if tracks.ntracks == 2:
					Mlll = selections.Mlll(dv_type=self.dv_type, plep=plep_vec, dMu=mu_vec, dEl=el_vec)
					Mhnl = selections.Mhnl(self.tree, self.dv_type, plep=plep_vec, dMu=mu_vec, dEl=el_vec)

					self.fill_ntuple(sel, 'mlll', Mlll.mlll)
					self.fill_ntuple(sel, 'mtransverse', Mlll.mtrans)
					self.fill_ntuple(sel, 'm_hnl', Mhnl.mhnl)
					self.fill_ntuple(sel, 'lifetime_hnl', Mhnl.lifetime_hnl)
					self.fill_ntuple(sel, 'alt_mhnl', Mhnl.alt_mhnl)
					self.fill_ntuple(sel, 'pt_hnl', Mhnl.hnlpt)
					self.fill_ntuple(sel, 'eta_hnl', Mhnl.hnleta)
					self.fill_ntuple(sel, 'phi_hnl', Mhnl.hnlphi)
					dR = tracks.lepVec[0].DeltaR(tracks.lepVec[1])
					if dR == 0.0:
						self.fill_ntuple(sel, 'DV_reduced_mass', -1)
						self.fill_ntuple(sel, 'DV_reduced_mlll', -1)
						self.fill_ntuple(sel, 'DV_reduced_mhnl', -1)
					else:
						self.fill_ntuple(sel, 'DV_reduced_mass', self.tree.dv('mass') / dR)
						self.fill_ntuple(sel, 'DV_reduced_mlll', Mlll.mlll / dR)
						self.fill_ntuple(sel, 'DV_reduced_mhnl', Mhnl.mhnl / dR)

				# get invariant mass between prompt lepton + same flavour displaced lepton for Z->ll veto
				if self.plep == 'muon' and len(mu_vec) > 0:
					zmass_veto_var = selections.ZMassVeto(self.tree, plep_vec=self.plep_sel.plepVec, plep=self.plep, plep_charge=self.plep_sel.plep_charge, dv_type=self.dv_type)
					self.fill_ntuple(sel, 'mll_dMu_plep_is_OS', zmass_veto_var.mll_dMu_plep_is_OS)
					self.fill_ntuple(sel, 'mll_dMu_plep_is_SS', zmass_veto_var.mll_dMu_plep_is_SS)
					self.fill_ntuple(sel, 'mll_dMu_plep', zmass_veto_var.mll_dMu_plep)

				if self.plep == 'muon' and len(el_vec) > 0:
					zmass_veto_var = selections.ZMassVeto(self.tree, plep_vec=self.plep_sel.plepVec, plep=self.plep, plep_charge=self.plep_sel.plep_charge, dv_type=self.dv_type)
					self.fill_ntuple(sel, 'mll_dEl_plep_is_OS', zmass_veto_var.mll_dEl_plep_is_OS)
					self.fill_ntuple(sel, 'mll_dEl_plep_is_SS', zmass_veto_var.mll_dEl_plep_is_SS)
					self.fill_ntuple(sel, 'mll_dEl_plep', zmass_veto_var.mll_dEl_plep)

				if self.plep == 'electron' and len(mu_vec) > 0:
					zmass_veto_var = selections.ZMassVeto(self.tree, plep_vec=self.plep_sel.plepVec, plep=self.plep, plep_charge=self.plep_sel.plep_charge, dv_type=self.dv_type)
					self.fill_ntuple(sel, 'mll_dMu_plep_is_OS', zmass_veto_var.mll_dMu_plep_is_OS)
					self.fill_ntuple(sel, 'mll_dMu_plep_is_SS', zmass_veto_var.mll_dMu_plep_is_SS)
					self.fill_ntuple(sel, 'mll_dMu_plep', zmass_veto_var.mll_dMu_plep)

				if self.plep == 'electron' and len(el_vec) > 0:
					zmass_veto_var = selections.ZMassVeto(self.tree, plep_vec=self.plep_sel.plepVec, plep=self.plep, plep_charge=self.plep_sel.plep_charge, dv_type=self.dv_type)
					self.fill_ntuple(sel, 'mll_dEl_plep_is_OS', zmass_veto_var.mll_dEl_plep_is_OS)
					self.fill_ntuple(sel, 'mll_dEl_plep_is_SS', zmass_veto_var.mll_dEl_plep_is_SS)
					self.fill_ntuple(sel, 'mll_dEl_plep', zmass_veto_var.mll_dEl_plep)

			if self.do_prompt_track_cut:
				ptrk_vec = self.ptrk_sel.trkVec
				ptrkd0 = self.ptrk_sel.trkd0
				ptrkz0 = self.ptrk_sel.trkz0
				self.fill_ntuple(sel, 'ptrk_pt', ptrk_vec.Pt())
				self.fill_ntuple(sel, 'ptrk_eta', ptrk_vec.Eta())
				self.fill_ntuple(sel, 'ptrk_phi', ptrk_vec.Phi())
				self.fill_ntuple(sel, 'ptrk_d0', ptrkd0)
				self.fill_ntuple(sel, 'ptrk_z0', ptrkz0)

				if tracks.ntracks == 2:
					Mhnl = selections.Mhnl(self.tree, self.dv_type, plep=ptrk_vec, dMu=mu_vec, dEl=el_vec, use_tracks=True, trks=tracks.lepVec)
					self.fill_ntuple(sel, 'HNLm', Mhnl.mhnl)
					self.fill_ntuple(sel, 'HNLpt', Mhnl.hnlpt)
					self.fill_ntuple(sel, 'HNLeta', Mhnl.hnleta)
					self.fill_ntuple(sel, 'HNLphi', Mhnl.hnlphi)

			if tracks.ntracks == 2:
				deta = abs(tracks.eta[0] - tracks.eta[1])
				dphi = abs(tracks.lepVec[0].DeltaPhi(tracks.lepVec[1]))
				dpt = abs(tracks.pt[0] - tracks.pt[1])
				dR = tracks.lepVec[0].DeltaR(tracks.lepVec[1])
				self.fill_ntuple(sel, 'DV_trk_deta', deta)
				self.fill_ntuple(sel, 'DV_trk_dphi', dphi)
				self.fill_ntuple(sel, 'DV_trk_dpt', dpt)
				self.fill_ntuple(sel, 'DV_trk_dR', dR)

				cosmic_veto = selections.CosmicVeto(self.tree)
				self.fill_ntuple(sel, 'DV_cosmic_sep', cosmic_veto.separation)

				self.fill_ntuple(sel, 'DV_trk_max_chi2_toSV', max(self.tree.dv('trk_chi2_toSV')[0], self.tree.dv('trk_chi2_toSV')[1]))
				self.fill_ntuple(sel, 'DV_trk_min_chi2_toSV', min(self.tree.dv('trk_chi2_toSV')[0], self.tree.dv('trk_chi2_toSV')[1]))
				self.fill_ntuple(sel, 'DV_trk_max_d0_wrtSV', max(self.tree.dv('trk_d0_wrtSV')[0], self.tree.dv('trk_d0_wrtSV')[1]))
				self.fill_ntuple(sel, 'DV_trk_min_d0_wrtSV', min(self.tree.dv('trk_d0_wrtSV')[0], self.tree.dv('trk_d0_wrtSV')[1]))
				self.fill_ntuple(sel, 'DV_trk_max_errd0_wrtSV', max(self.tree.dv('trk_errd0_wrtSV')[0], self.tree.dv('trk_errd0_wrtSV')[1]))
				self.fill_ntuple(sel, 'DV_trk_min_errd0_wrtSV', min(self.tree.dv('trk_errd0_wrtSV')[0], self.tree.dv('trk_errd0_wrtSV')[1]))
				self.fill_ntuple(sel, 'DV_trk_max_z0_wrtSV', max(self.tree.dv('trk_z0_wrtSV')[0], self.tree.dv('trk_z0_wrtSV')[1]))
				self.fill_ntuple(sel, 'DV_trk_min_z0_wrtSV', min(self.tree.dv('trk_z0_wrtSV')[0], self.tree.dv('trk_z0_wrtSV')[1]))
				self.fill_ntuple(sel, 'DV_trk_max_errz0_wrtSV', max(self.tree.dv('trk_errz0_wrtSV')[0], self.tree.dv('trk_errz0_wrtSV')[1]))
				self.fill_ntuple(sel, 'DV_trk_min_errz0_wrtSV', min(self.tree.dv('trk_errz0_wrtSV')[0], self.tree.dv('trk_errz0_wrtSV')[1]))

				DV_mumu = selections.DVType(self.tree, dv_type="mumu").passes()
				DV_ee = selections.DVType(self.tree, dv_type="ee").passes()
				DV_emu = selections.DVType(self.tree, dv_type="emu").passes()
				DV_1lep = (len(mu_vec) == 1 and len(el_vec) == 0) or (len(mu_vec) == 0 and len(el_vec) == 1)

				self.fill_ntuple(sel, 'DV_mumu', DV_mumu)
				self.fill_ntuple(sel, 'DV_ee', DV_ee)
				self.fill_ntuple(sel, 'DV_emu', DV_emu)
				self.fill_ntuple(sel, 'DV_1lepton', DV_1lep)
				# pass_el_mu_overlap = selections.ElectronMuonOverlapCheck(self.tree).passes()
				# self.fill_ntuple(sel, 'DV_pass_el_mu_overlap', pass_el_mu_overlap)

				pass_lep_pt_cut = selections.DVLepPt(self.tree, self.dv_type).pass_pt_cuts
				self.fill_ntuple(sel, 'DV_pass_lepton_pt', pass_lep_pt_cut)

				# calculate momentum parallel and perpendicular to the decay vector = DV-PV
				dv = ROOT.TVector3(self.tree.dv('x'), self.tree.dv('y'), self.tree.dv('z'))
				pv = ROOT.TVector3(self.tree['truth_PV_x'], self.tree['truth_PV_y'], self.tree['truth_PV_z'])
				decay_vector = dv - pv
				pvec_0 = ROOT.TVector3(tracks.lepVec[0].Px(), tracks.lepVec[0].Py(), tracks.lepVec[0].Pz())
				pvec_1 = ROOT.TVector3(tracks.lepVec[1].Px(), tracks.lepVec[1].Py(), tracks.lepVec[1].Pz())

				mom_perp_0 = helpers.mom_perp(pvec_0, decay_vector)
				mom_parall_0 = helpers.mom_parall(pvec_0, decay_vector)
				mom_frac_parall_0 = helpers.mom_frac_parall(pvec_0, decay_vector)
				pvec_0_mag = pvec_0.Mag()

				mom_perp_1 = helpers.mom_perp(pvec_1, decay_vector)
				mom_parall_1 = helpers.mom_parall(pvec_1, decay_vector)
				mom_frac_parall_1 = helpers.mom_frac_parall(pvec_1, decay_vector)
				pvec_1_mag = pvec_1.Mag()

				# pt order the visible leptons in the DV
				if self.tree.dv('trk_pt_wrtSV')[1] > self.tree.dv('trk_pt_wrtSV')[0]:
					self.fill_ntuple(sel, 'DV_trk_0_pt', self.tree.dv('trk_pt_wrtSV')[1])
					self.fill_ntuple(sel, 'DV_trk_0_eta', self.tree.dv('trk_eta_wrtSV')[1])
					self.fill_ntuple(sel, 'DV_trk_0_phi', self.tree.dv('trk_phi_wrtSV')[1])
					self.fill_ntuple(sel, 'DV_trk_0_d0', self.tree.dv('trk_d0')[1])
					self.fill_ntuple(sel, 'DV_trk_0_z0', self.tree.dv('trk_z0')[1])
					self.fill_ntuple(sel, 'DV_trk_0_charge', self.tree.dv('trk_charge')[1])
					self.fill_ntuple(sel, 'DV_trk_0_chi2', self.tree.dv('trk_chi2')[1])
					self.fill_ntuple(sel, 'DV_trk_0_isSelected', self.tree.dv('trk_isSelected')[1])
					self.fill_ntuple(sel, 'DV_trk_0_isAssociated', self.tree.dv('trk_isAssociated')[1])
					self.fill_ntuple(sel, 'DV_trk_0_mom_parall', mom_parall_1)
					self.fill_ntuple(sel, 'DV_trk_0_mom_perp', mom_perp_1)
					self.fill_ntuple(sel, 'DV_trk_0_mom_mag', pvec_1_mag)
					self.fill_ntuple(sel, 'DV_trk_0_mom_frac_parall', mom_frac_parall_1)

					self.fill_ntuple(sel, 'DV_trk_1_pt', self.tree.dv('trk_pt_wrtSV')[0])
					self.fill_ntuple(sel, 'DV_trk_1_eta', self.tree.dv('trk_eta_wrtSV')[0])
					self.fill_ntuple(sel, 'DV_trk_1_phi', self.tree.dv('trk_phi_wrtSV')[0])
					self.fill_ntuple(sel, 'DV_trk_1_d0', self.tree.dv('trk_d0')[0])
					self.fill_ntuple(sel, 'DV_trk_1_z0', self.tree.dv('trk_z0')[0])
					self.fill_ntuple(sel, 'DV_trk_1_charge', self.tree.dv('trk_charge')[0])
					self.fill_ntuple(sel, 'DV_trk_1_chi2', self.tree.dv('trk_chi2')[0])
					self.fill_ntuple(sel, 'DV_trk_1_isSelected', self.tree.dv('trk_isSelected')[0])
					self.fill_ntuple(sel, 'DV_trk_1_isAssociated', self.tree.dv('trk_isAssociated')[0])
					self.fill_ntuple(sel, 'DV_trk_1_mom_parall', mom_parall_0)
					self.fill_ntuple(sel, 'DV_trk_1_mom_perp', mom_perp_0)
					self.fill_ntuple(sel, 'DV_trk_1_mom_mag', pvec_0_mag)
					self.fill_ntuple(sel, 'DV_trk_1_mom_frac_parall', mom_frac_parall_0)
				else:
					self.fill_ntuple(sel, 'DV_trk_0_pt', self.tree.dv('trk_pt_wrtSV')[0])
					self.fill_ntuple(sel, 'DV_trk_0_eta', self.tree.dv('trk_eta_wrtSV')[0])
					self.fill_ntuple(sel, 'DV_trk_0_phi', self.tree.dv('trk_phi_wrtSV')[0])
					self.fill_ntuple(sel, 'DV_trk_0_d0', self.tree.dv('trk_d0')[0])
					self.fill_ntuple(sel, 'DV_trk_0_z0', self.tree.dv('trk_z0')[0])
					self.fill_ntuple(sel, 'DV_trk_0_charge', self.tree.dv('trk_charge')[0])
					self.fill_ntuple(sel, 'DV_trk_0_chi2', self.tree.dv('trk_chi2')[0])
					self.fill_ntuple(sel, 'DV_trk_0_isSelected', self.tree.dv('trk_isSelected')[0])
					self.fill_ntuple(sel, 'DV_trk_0_isAssociated', self.tree.dv('trk_isAssociated')[0])
					self.fill_ntuple(sel, 'DV_trk_0_mom_parall', mom_parall_0)
					self.fill_ntuple(sel, 'DV_trk_0_mom_perp', mom_perp_0)
					self.fill_ntuple(sel, 'DV_trk_0_mom_mag', pvec_0_mag)
					self.fill_ntuple(sel, 'DV_trk_0_mom_frac_parall', mom_frac_parall_0)

					self.fill_ntuple(sel, 'DV_trk_1_pt', self.tree.dv('trk_pt_wrtSV')[1])
					self.fill_ntuple(sel, 'DV_trk_1_eta', self.tree.dv('trk_eta_wrtSV')[1])
					self.fill_ntuple(sel, 'DV_trk_1_phi', self.tree.dv('trk_phi_wrtSV')[1])
					self.fill_ntuple(sel, 'DV_trk_1_d0', self.tree.dv('trk_d0')[1])
					self.fill_ntuple(sel, 'DV_trk_1_z0', self.tree.dv('trk_z0')[1])
					self.fill_ntuple(sel, 'DV_trk_1_charge', self.tree.dv('trk_charge')[1])
					self.fill_ntuple(sel, 'DV_trk_1_chi2', self.tree.dv('trk_chi2')[1])
					self.fill_ntuple(sel, 'DV_trk_1_isSelected', self.tree.dv('trk_isSelected')[1])
					self.fill_ntuple(sel, 'DV_trk_1_isAssociated', self.tree.dv('trk_isAssociated')[1])
					self.fill_ntuple(sel, 'DV_trk_1_mom_parall', mom_parall_1)
					self.fill_ntuple(sel, 'DV_trk_1_mom_perp', mom_perp_1)
					self.fill_ntuple(sel, 'DV_trk_1_mom_mag', pvec_1_mag)
					self.fill_ntuple(sel, 'DV_trk_1_mom_frac_parall', mom_frac_parall_1)

			# fill 3 different track calculations
			if self.dv_type == "mumu":
				assert(tracks.ntracks == 2)
				lep_matched_DV_vec = muons.lepmatched_lepVec[0] + muons.lepmatched_lepVec[1]
				trk_percent_diff_0 = helpers.pT_diff(muons.std_lepVec[0].Pt(), muons.lepmatched_lepVec[0].Pt() )
				trk_percent_diff_1 = helpers.pT_diff( muons.std_lepVec[1].Pt(), muons.lepmatched_lepVec[1].Pt() )
				trk_percent_diff_truth_lep_matched_0 = helpers.pT_diff(  muons.lepmatched_lepVec[0].Pt(), muons.std_lepVec[0].Pt() )
				trk_percent_diff_truth_lep_matched_1 = helpers.pT_diff(  muons.lepmatched_lepVec[1].Pt(), muons.std_lepVec[1].Pt() )
				dv_mass_diff = abs(lep_matched_DV_vec.M() - self.tree.dv('mass') ) / lep_matched_DV_vec.M()
				# trigger matching for displaced leptons
				dmu_0_is_trig_matched = selections.LeptonTriggerMatching(self.tree, "muon", mu_index[0]).is_trigger_matched
				dmu_1_is_trig_matched = selections.LeptonTriggerMatching(self.tree, "muon", mu_index[1]).is_trigger_matched

				self.fill_ntuple(sel, 'DV_mu_0_trk_pt_wrtSV', muons.lepVec[0].Pt())
				self.fill_ntuple(sel, 'DV_mu_1_trk_pt_wrtSV', muons.lepVec[1].Pt())
				self.fill_ntuple(sel, 'DV_mu_0_std_trk_pt', muons.std_lepVec[0].Pt())
				self.fill_ntuple(sel, 'DV_mu_1_std_trk_pt', muons.std_lepVec[1].Pt())
				self.fill_ntuple(sel, 'DV_mu_0_lepmatched_trk_pt', muons.lepmatched_lepVec[0].Pt())
				self.fill_ntuple(sel, 'DV_mu_1_lepmatched_trk_pt', muons.lepmatched_lepVec[1].Pt())
				self.fill_ntuple(sel, 'DV_mu_0_lepmatched_trk_eta', muons.lepmatched_lepVec[0].Eta())
				self.fill_ntuple(sel, 'DV_mu_1_lepmatched_trk_eta', muons.lepmatched_lepVec[1].Eta())
				self.fill_ntuple(sel, 'DV_mu_0_lepmatched_trk_phi', muons.lepmatched_lepVec[0].Phi())
				self.fill_ntuple(sel, 'DV_mu_1_lepmatched_trk_phi', muons.lepmatched_lepVec[1].Phi())
				self.fill_ntuple(sel, 'DV_mass_lepmatched', lep_matched_DV_vec.M())
				self.fill_ntuple(sel, 'DV_mass_diff', dv_mass_diff)
				self.fill_ntuple(sel, 'DV_mu_0_pt_diff', trk_percent_diff_0 )
				self.fill_ntuple(sel, 'DV_mu_1_pt_diff', trk_percent_diff_1 )
				self.fill_ntuple(sel, 'DV_mu_0_pt_diff_lep_matched', trk_percent_diff_truth_lep_matched_0 )
				self.fill_ntuple(sel, 'DV_mu_1_pt_diff_lep_matched', trk_percent_diff_truth_lep_matched_1 )
				self.fill_ntuple(sel, 'DV_mu_0_is_trigger_matched', dmu_0_is_trig_matched)
				self.fill_ntuple(sel, 'DV_mu_1_is_trigger_matched', dmu_1_is_trig_matched)
				self.fill_ntuple(sel, 'DV_mu_0_isMuon', 1)
				self.fill_ntuple(sel, 'DV_mu_1_isMuon', 1)
				self.fill_ntuple(sel, 'DV_mu_0_charge', muons.lepCharge[0])
				self.fill_ntuple(sel, 'DV_mu_1_charge', muons.lepCharge[1])
				self.fill_ntuple(sel, 'DV_mu_0_isElectron', 0)
				self.fill_ntuple(sel, 'DV_mu_1_isElectron', 0)
				self.fill_ntuple(sel, 'DV_mu_0_isLRT', muons.muon_isLRT[0])
				self.fill_ntuple(sel, 'DV_mu_1_isLRT', muons.muon_isLRT[1])
				self.fill_ntuple(sel, 'DV_mu_0_muon_isLoose', self.tree.get('muon_isLoose')[muons.lepIndex[0]])
				self.fill_ntuple(sel, 'DV_mu_1_muon_isLoose', self.tree.get('muon_isLoose')[muons.lepIndex[1]])
				self.fill_ntuple(sel, 'DV_mu_0_muon_isMedium', self.tree.get('muon_isMedium')[muons.lepIndex[0]])
				self.fill_ntuple(sel, 'DV_mu_1_muon_isMedium', self.tree.get('muon_isMedium')[muons.lepIndex[1]])
				self.fill_ntuple(sel, 'DV_mu_0_muon_isTight', self.tree.get('muon_isTight')[muons.lepIndex[0]])
				self.fill_ntuple(sel, 'DV_mu_1_muon_isTight', self.tree.get('muon_isTight')[muons.lepIndex[1]])

			if self.dv_type == "emu":
				assert(tracks.ntracks == 2)
				lep_matched_DV_vec = muons.lepmatched_lepVec[0] + electrons.lepmatched_lepVec[0]
				trk_percent_diff_0 = helpers.pT_diff( muons.std_lepVec[0].Pt(), muons.lepmatched_lepVec[0].Pt() )
				trk_percent_diff_1 = helpers.pT_diff( electrons.std_lepVec[0].Pt() , electrons.lepmatched_lepVec[0].Pt() )
				trk_percent_diff_truth_lep_matched_0 = helpers.pT_diff(  muons.lepmatched_lepVec[0].Pt(), muons.std_lepVec[0].Pt() )
				trk_percent_diff_truth_lep_matched_1 = helpers.pT_diff( electrons.lepmatched_lepVec[0].Pt(), electrons.std_lepVec[0].Pt() )
				dv_mass_diff = abs(lep_matched_DV_vec.M() - self.tree.dv('mass') ) / lep_matched_DV_vec.M()
				# trigger matching for displaced leptons
				dmu_is_trig_matched = selections.LeptonTriggerMatching(self.tree, "muon", mu_index[0]).is_trigger_matched
				del_is_trig_matched = selections.LeptonTriggerMatching(self.tree, "electron", el_index[0]).is_trigger_matched

				self.fill_ntuple(sel, 'DV_mu_0_trk_pt_wrtSV', muons.lepVec[0].Pt())
				self.fill_ntuple(sel, 'DV_el_1_trk_pt_wrtSV', electrons.lepVec[0].Pt())
				self.fill_ntuple(sel, 'DV_mu_0_std_trk_pt', muons.std_lepVec[0].Pt())
				self.fill_ntuple(sel, 'DV_el_1_std_trk_pt', electrons.std_lepVec[0].Pt())
				self.fill_ntuple(sel, 'DV_mu_0_lepmatched_trk_pt', muons.lepmatched_lepVec[0].Pt())
				self.fill_ntuple(sel, 'DV_el_1_lepmatched_trk_pt', electrons.lepmatched_lepVec[0].Pt())
				self.fill_ntuple(sel, 'DV_mu_0_lepmatched_trk_eta', muons.lepmatched_lepVec[0].Eta())
				self.fill_ntuple(sel, 'DV_el_1_lepmatched_trk_eta', electrons.lepmatched_lepVec[0].Eta())
				self.fill_ntuple(sel, 'DV_mu_0_lepmatched_trk_phi', muons.lepmatched_lepVec[0].Phi())
				self.fill_ntuple(sel, 'DV_el_1_lepmatched_trk_phi', electrons.lepmatched_lepVec[0].Phi())
				self.fill_ntuple(sel, 'DV_mass_lepmatched', lep_matched_DV_vec.M())
				self.fill_ntuple(sel, 'DV_mass_diff', dv_mass_diff)
				self.fill_ntuple(sel, 'DV_mu_0_pt_diff', trk_percent_diff_0 )
				self.fill_ntuple(sel, 'DV_el_1_pt_diff', trk_percent_diff_1 )
				self.fill_ntuple(sel, 'DV_mu_0_pt_diff_lep_matched', trk_percent_diff_truth_lep_matched_0 )
				self.fill_ntuple(sel, 'DV_el_1_pt_diff_lep_matched', trk_percent_diff_truth_lep_matched_1 )
				self.fill_ntuple(sel, 'DV_mu_0_is_trigger_matched', dmu_is_trig_matched)
				self.fill_ntuple(sel, 'DV_el_1_is_trigger_matched', del_is_trig_matched)
				self.fill_ntuple(sel, 'DV_mu_0_isMuon', 1)
				self.fill_ntuple(sel, 'DV_el_1_isMuon', 0)
				self.fill_ntuple(sel, 'DV_mu_0_isElectron', 0)
				self.fill_ntuple(sel, 'DV_el_1_isElectron', 1)
				self.fill_ntuple(sel, 'DV_mu_0_charge', muons.lepCharge[0])
				self.fill_ntuple(sel, 'DV_el_1_charge', electrons.lepCharge[0])

				self.fill_ntuple(sel, 'DV_mu_0_muon_isLoose', self.tree.get('muon_isLoose')[muons.lepIndex[0]])
				self.fill_ntuple(sel, 'DV_mu_0_muon_isMedium', self.tree.get('muon_isMedium')[muons.lepIndex[0]])
				self.fill_ntuple(sel, 'DV_mu_0_muon_isTight', self.tree.get('muon_isTight')[muons.lepIndex[0]])
				self.fill_ntuple(sel, 'DV_el_1_electron_LHTight', self.tree.get('el_LHTight')[electrons.lepIndex[0]])
				self.fill_ntuple(sel, 'DV_el_1_electron_LHMedium', self.tree.get('el_LHMedium')[electrons.lepIndex[0]])
				self.fill_ntuple(sel, 'DV_el_1_electron_LHLoose', self.tree.get('el_LHLoose')[electrons.lepIndex[0]])
				#self.fill_ntuple(sel, 'DV_el_1_electron_isLHVeryLoose', self.tree.get('el_isLHVeryLoose')[electrons.lepIndex[0]])
				#self.fill_ntuple(sel, 'DV_el_1_electron_VeryVeryLoose', self.tree.get('el_isLHVeryLoose_mod1')[electrons.lepIndex[0]])
				#self.fill_ntuple(sel, 'DV_el_1_electron_VeryVeryLooseSi', self.tree.get('el_isLHVeryLoose_modSi')[electrons.lepIndex[0]])

			if self.dv_type == "ee":
				assert(tracks.ntracks == 2)
				lep_matched_DV_vec = electrons.lepmatched_lepVec[0] + electrons.lepmatched_lepVec[1]
				trk_percent_diff_0 = helpers.pT_diff( electrons.std_lepVec[0].Pt() , electrons.lepmatched_lepVec[0].Pt() )
				trk_percent_diff_1 = helpers.pT_diff( electrons.std_lepVec[1].Pt(), electrons.lepmatched_lepVec[1].Pt() )
				trk_percent_diff_truth_lep_matched_0 = helpers.pT_diff(  electrons.lepmatched_lepVec[0].Pt(), electrons.std_lepVec[0].Pt() )
				trk_percent_diff_truth_lep_matched_1 = helpers.pT_diff(  electrons.lepmatched_lepVec[1].Pt(), electrons.std_lepVec[1].Pt() )
				dv_mass_diff = abs(lep_matched_DV_vec.M() - self.tree.dv('mass') ) / lep_matched_DV_vec.M()
				# trigger matching for displaced leptons
				del_0_is_trig_matched = selections.LeptonTriggerMatching(self.tree, "electron", el_index[0]).is_trigger_matched
				del_1_is_trig_matched = selections.LeptonTriggerMatching(self.tree, "electron", el_index[1]).is_trigger_matched

				self.fill_ntuple(sel, 'DV_el_0_trk_pt_wrtSV', electrons.lepVec[0].Pt())
				self.fill_ntuple(sel, 'DV_el_1_trk_pt_wrtSV', electrons.lepVec[1].Pt())
				self.fill_ntuple(sel, 'DV_el_0_std_trk_pt', electrons.std_lepVec[0].Pt())
				self.fill_ntuple(sel, 'DV_el_1_std_trk_pt', electrons.std_lepVec[1].Pt())
				self.fill_ntuple(sel, 'DV_el_0_lepmatched_trk_pt', electrons.lepmatched_lepVec[0].Pt())
				self.fill_ntuple(sel, 'DV_el_1_lepmatched_trk_pt', electrons.lepmatched_lepVec[1].Pt())
				self.fill_ntuple(sel, 'DV_el_0_lepmatched_trk_eta', electrons.lepmatched_lepVec[0].Eta())
				self.fill_ntuple(sel, 'DV_el_1_lepmatched_trk_eta', electrons.lepmatched_lepVec[1].Eta())
				self.fill_ntuple(sel, 'DV_el_0_lepmatched_trk_phi', electrons.lepmatched_lepVec[0].Phi())
				self.fill_ntuple(sel, 'DV_el_1_lepmatched_trk_phi', electrons.lepmatched_lepVec[1].Phi())
				self.fill_ntuple(sel, 'DV_mass_lepmatched', lep_matched_DV_vec.M())
				self.fill_ntuple(sel, 'DV_mass_diff', dv_mass_diff)
				self.fill_ntuple(sel, 'DV_el_0_pt_diff', trk_percent_diff_0 )
				self.fill_ntuple(sel, 'DV_el_1_pt_diff', trk_percent_diff_1 )
				self.fill_ntuple(sel, 'DV_el_0_pt_diff_lep_matched', trk_percent_diff_truth_lep_matched_0 )
				self.fill_ntuple(sel, 'DV_el_1_pt_diff_lep_matched', trk_percent_diff_truth_lep_matched_1 )
				self.fill_ntuple(sel, 'DV_el_0_is_trigger_matched', del_0_is_trig_matched)
				self.fill_ntuple(sel, 'DV_el_1_is_trigger_matched', del_1_is_trig_matched)
				self.fill_ntuple(sel, 'DV_el_0_isElectron', 1)
				self.fill_ntuple(sel, 'DV_el_1_isElectron', 1)
				self.fill_ntuple(sel, 'DV_el_0_isMuon', 0)
				self.fill_ntuple(sel, 'DV_el_1_isMuon', 0)
				self.fill_ntuple(sel, 'DV_el_0_charge', electrons.lepCharge[0])
				self.fill_ntuple(sel, 'DV_el_1_charge', electrons.lepCharge[1])
				self.fill_ntuple(sel, 'DV_el_0_electron_LHTight', self.tree.get('el_LHTight')[electrons.lepIndex[0]])
				self.fill_ntuple(sel, 'DV_el_0_electron_LHMedium', self.tree.get('el_LHMedium')[electrons.lepIndex[0]])
				self.fill_ntuple(sel, 'DV_el_0_electron_LHLoose', self.tree.get('el_LHLoose')[electrons.lepIndex[0]])
				#self.fill_ntuple(sel, 'DV_el_0_electron_isLHVeryLoose', self.tree.get('el_isLHVeryLoose')[electrons.lepIndex[0]])
				#self.fill_ntuple(sel, 'DV_el_0_electron_VeryVeryLoose', self.tree.get('el_isLHVeryLoose_mod1')[electrons.lepIndex[0]])
				#self.fill_ntuple(sel, 'DV_el_0_electron_VeryVeryLooseSi', self.tree.get('el_isLHVeryLoose_modSi')[electrons.lepIndex[0]])
				self.fill_ntuple(sel, 'DV_el_1_electron_LHTight', self.tree.get('el_LHTight')[electrons.lepIndex[1]])
				self.fill_ntuple(sel, 'DV_el_1_electron_LHMedium', self.tree.get('el_LHMedium')[electrons.lepIndex[1]])
				self.fill_ntuple(sel, 'DV_el_1_electron_LHLoose', self.tree.get('el_LHLoose')[electrons.lepIndex[1]])
				#self.fill_ntuple(sel, 'DV_el_1_electron_isLHVeryLoose', self.tree.get('el_isLHVeryLoose')[electrons.lepIndex[1]])
				#self.fill_ntuple(sel, 'DV_el_1_electron_VeryVeryLoose', self.tree.get('el_isLHVeryLoose_mod1')[electrons.lepIndex[1]])
				#self.fill_ntuple(sel, 'DV_el_1_electron_VeryVeryLooseSi', self.tree.get('el_isLHVeryLoose_modSi')[electrons.lepIndex[1]])

			# fill standard track variables for electrons
			for lep in range(len(el_vec)):
				self.fill_ntuple(sel, 'DV_el_trk_pt', el_vec[lep].Pt())
				self.fill_ntuple(sel, 'DV_el_trk_eta', el_vec[lep].Eta())
				self.fill_ntuple(sel, 'DV_el_trk_phi', el_vec[lep].Phi())

			# fill standard track variables for muons
			for lep in range(len(mu_vec)):
				self.fill_ntuple(sel, 'DV_mu_trk_pt', mu_vec[lep].Pt())
				self.fill_ntuple(sel, 'DV_mu_trk_eta', mu_vec[lep].Eta())
				self.fill_ntuple(sel, 'DV_mu_trk_phi', mu_vec[lep].Phi())

			# fill standard track variable histograms
			for i in range(tracks.ntracks):
				self.fill_ntuple(sel, 'DV_trk_pt', self.tree.dv('trk_pt_wrtSV')[i])
				self.fill_ntuple(sel, 'DV_trk_eta', self.tree.dv('trk_eta_wrtSV')[i])
				self.fill_ntuple(sel, 'DV_trk_phi', self.tree.dv('trk_phi_wrtSV')[i])
				self.fill_ntuple(sel, 'DV_trk_d0', self.tree.dv('trk_d0')[i])
				self.fill_ntuple(sel, 'DV_trk_z0', self.tree.dv('trk_z0')[i])
				self.fill_ntuple(sel, 'DV_trk_absz0', abs(self.tree.dv('trk_z0')[i]))
				self.fill_ntuple(sel, 'DV_trk_charge', self.tree.dv('trk_charge')[i])
				self.fill_ntuple(sel, 'DV_trk_chi2', self.tree.dv('trk_chi2')[i])
				# self.fill_ntuple(sel, 'DV_trk_isLRT', self.tree.dv('trk_isLRT')[i])
				self.fill_ntuple(sel, 'DV_trk_isSelected', self.tree.dv('trk_isSelected')[i])
				self.fill_ntuple(sel, 'DV_trk_isAssociated', self.tree.dv('trk_isAssociated')[i])
				self.fill_ntuple(sel, 'DV_trk_nPixelHits', self.tree.dv('trk_nPixelHits')[i])
				self.fill_ntuple(sel, 'DV_trk_nSCTHits', self.tree.dv('trk_nSCTHits')[i])
				# self.fill_ntuple(sel, 'DV_trk_nSCTHoles', self.tree.dv('trk_nSCTHoles')[i])
				self.fill_ntuple(sel, 'DV_trk_nSiHits', self.tree.dv('trk_nSCTHits')[i] + self.tree.dv('trk_nPixelHits')[i])
				# self.fill_ntuple(sel, 'DV_trk_dTheta', self.tree.dv('trk_dTheta')[i])
				self.fill_ntuple(sel, 'DV_trk_chi2_toSV'.format(i), self.tree.dv('trk_chi2_toSV')[i])
				self.fill_ntuple(sel, 'DV_trk_d0_wrtSV'.format(i), self.tree.dv('trk_d0_wrtSV')[i])
				self.fill_ntuple(sel, 'DV_trk_errd0_wrtSV'.format(i), self.tree.dv('trk_errd0_wrtSV')[i])
				self.fill_ntuple(sel, 'DV_trk_z0_wrtSV'.format(i), self.tree.dv('trk_z0_wrtSV')[i])
				self.fill_ntuple(sel, 'DV_trk_errz0_wrtSV'.format(i), self.tree.dv('trk_errz0_wrtSV')[i])

			# fill standard dv histograms
			self.fill_ntuple(sel, 'DV_num_trks', self.tree.dv('ntrk'))
			self.fill_ntuple(sel, 'DV_x', self.tree.dv('x'))
			self.fill_ntuple(sel, 'DV_y', self.tree.dv('y'))
			self.fill_ntuple(sel, 'DV_z', self.tree.dv('z'))
			self.fill_ntuple(sel, 'DV_r', self.tree.dv('r'))
			self.fill_ntuple(sel, 'PV_x', self.tree['truth_PV_x'])
			self.fill_ntuple(sel, 'PV_y', self.tree['truth_PV_y'])
			self.fill_ntuple(sel, 'PV_z', self.tree['truth_PV_z'])
			self.fill_ntuple(sel, 'DV_distFromPV', self.tree.dv('distFromPV'))
			self.fill_ntuple(sel, 'DV_mass', self.tree.dv('mass'))
			self.fill_ntuple(sel, 'DV_pt', self.tree.dv('pt'))
			self.fill_ntuple(sel, 'DV_eta', self.tree.dv('eta'))
			self.fill_ntuple(sel, 'DV_phi', self.tree.dv('phi'))
			self.fill_ntuple(sel, 'DV_minOpAng', self.tree.dv('minOpAng'))
			self.fill_ntuple(sel, 'DV_maxOpAng', self.tree.dv('maxOpAng'))
			self.fill_ntuple(sel, 'DV_charge', self.tree.dv('charge'))
			self.fill_ntuple(sel, 'DV_chi2', self.tree.dv('chi2'))
			# self.fill_ntuple(sel, 'DV_chi2_assoc', self.tree.dv('chi2_assoc'))
			self.fill_ntuple(sel, 'DV_max_dR', self.tree.dv('maxDR'))
			self.fill_ntuple(sel, 'DV_max_dR_wrtSV', self.tree.dv('maxDR_wrtSV'))
			self.fill_ntuple(sel, 'DV_maxd0', self.tree.dv('maxd0'))
			self.fill_ntuple(sel, 'DV_mind0', self.tree.dv('mind0'))
			self.fill_ntuple(sel, 'DV_ntrk', self.tree.dv('ntrk'))
			self.fill_ntuple(sel, 'DV_ntrk_lrt', self.tree.dv('ntrk_lrt'))
			self.fill_ntuple(sel, 'DV_ntrk_sel', self.tree.dv('ntrk_sel'))
			self.fill_ntuple(sel, 'DV_ntrk_assoc', self.tree.dv('ntrk_assoc'))
			self.fill_ntuple(sel, 'DV_pass_mat_veto', selections.MaterialVeto(self.tree).passes())

			# Christian: TODO: all of the things below have to be reimplemented:

			# if not self.tree.is_data:
			# 	# is truth matched:
			# 	# self.fill_ntuple(sel, 'DV_truth_matched', self._truth_match())

			# 	# add proper lifetime
			# 	truth_info = helpers.Truth()
			# 	truth_info.get_truth_particles(self.tree)
			# 	self.fill_ntuple(sel, 'properLifetime', truth_info.properLifetime)

			# 	self.fill_ntuple(sel, 'truth_DV_mass', truth_info.DV_mass)

			# 	if self.dv_type == "mumu":
			# 		# Get the truth index (truth matching by charge)
			# 		if not self.tree.is_bkg_mc and 'tt' not in self.tree.mc_ch_str:
			# 			if self.tree.dv('trk_charge')[0] == truth_info.dMu_charge[0]: trk_0_truth_index = 0
			# 			elif self.tree.dv('trk_charge')[0] == truth_info.dMu_charge[1]: trk_0_truth_index = 1
			# 			else: raise Exception("Can't truth match lepton by charge. Something is strange.")

			# 			if self.tree.dv('trk_charge')[1] == truth_info.dMu_charge[0]: trk_1_truth_index = 0
			# 			elif self.tree.dv('trk_charge')[1] == truth_info.dMu_charge[1]: trk_1_truth_index = 1
			# 			else: raise Exception("Can't truth match lepton by charge. Something is strange.")

			# 			self.fill_ntuple(sel, 'DV_trk_0_d0_truth', truth_info.dMu_d0[trk_0_truth_index])
			# 			self.fill_ntuple(sel, 'DV_trk_1_d0_truth', truth_info.dMu_d0[trk_1_truth_index])

			# 			for i in range(len(muons.lepVec)):
			# 				delta = muons.lepVec[i].Pt() - muons.lepmatched_lepVec[i].Pt()
			# 				self.fill_ntuple(sel, 'DV_trk_v_mu_pt', delta / muons.lepmatched_lepVec[i].Pt())

			# 	if self.dv_type == "emu":
			# 		# truth d0 # it looks like these are already matched so track0 is always the muon. Could that be?
			# 		if not self.tree.is_bkg_mc and 'tt' not in self.tree.mc_ch_str:
			# 			self.fill_ntuple(sel, 'DV_trk_0_d0_truth', truth_info.dMu_d0[0])
			# 			self.fill_ntuple(sel, 'DV_trk_1_d0_truth', truth_info.dEl_d0[0])

			# 	if self.dv_type == "ee":
			# 		if not self.tree.is_bkg_mc and 'tt' not in self.tree.mc_ch_str:
			# 			# Get the truth index (truth matching by charge)
			# 			if self.tree.dv('trk_charge')[0] == truth_info.dEl_charge[0]: trk_0_truth_index = 0
			# 			elif self.tree.dv('trk_charge')[0] == truth_info.dEl_charge[1]: trk_0_truth_index = 1
			# 			else:	raise Exception("Can't truth match lepton by charge. Something is strange.")

			# 			if self.tree.dv('trk_charge')[1] == truth_info.dEl_charge[0]: trk_1_truth_index = 0
			# 			elif self.tree.dv('trk_charge')[1] == truth_info.dEl_charge[1]: trk_1_truth_index = 1
			# 			else:	raise Exception("Can't truth match lepton by charge. Something is strange.")

			# 			self.fill_ntuple(sel, 'DV_trk_0_d0_truth', truth_info.dEl_d0[trk_0_truth_index])
			# 			self.fill_ntuple(sel, 'DV_trk_1_d0_truth', truth_info.dEl_d0[trk_1_truth_index])

			# 			for i in range(len(electrons.lepVec)):
			# 				delta = electrons.lepVec[i].Pt() - electrons.lepmatched_lepVec[i].Pt()
			# 				self.fill_ntuple(sel, 'DV_trk_v_el_pt', delta/electrons.lepmatched_lepVec[i].Pt() )



			self.fill_ntuple(sel, 'DV_alpha', selections.Alpha(self.tree).alpha)
			
			trk_quality = selections.TrackQuality(self.tree)

			self.fill_ntuple(sel, 'DV_2tight', trk_quality.DV_2tight)
			self.fill_ntuple(sel, 'DV_2medium', trk_quality.DV_2medium)
			self.fill_ntuple(sel, 'DV_2loose', trk_quality.DV_2loose)
			self.fill_ntuple(sel, 'DV_1tight', trk_quality.DV_1tight)
			self.fill_ntuple(sel, 'DV_1medium', trk_quality.DV_1medium)
			self.fill_ntuple(sel, 'DV_1loose', trk_quality.DV_1loose)
			self.fill_ntuple(sel, 'DV_tight_loose', trk_quality.DV_tight_loose)
			self.fill_ntuple(sel, 'DV_tight_medium', trk_quality.DV_tight_medium)
			self.fill_ntuple(sel, 'DV_medium_loose', trk_quality.DV_medium_loose)
            # Can't handle veryloose or veryveryloose WP yet
			#self.fill_ntuple(sel, 'DV_tight_veryloose', trk_quality.DV_tight_veryloose)
			#self.fill_ntuple(sel, 'DV_medium_veryloose', trk_quality.DV_medium_veryloose)
			#self.fill_ntuple(sel, 'DV_loose_veryloose', trk_quality.DV_loose_veryloose)
			#self.fill_ntuple(sel, 'DV_tight_veryveryloose', trk_quality.DV_tight_veryveryloose)
			#self.fill_ntuple(sel, 'DV_medium_veryveryloose', trk_quality.DV_medium_veryveryloose)
			#self.fill_ntuple(sel, 'DV_loose_veryveryloose', trk_quality.DV_loose_veryveryloose)
			#self.fill_ntuple(sel, 'DV_2veryveryloose', trk_quality.DV_2veryveryloose)
			#self.fill_ntuple(sel, 'DV_1veryveryloose', trk_quality.DV_1veryveryloose)

			# ____________________________________________________________
			# Trigger matching requirement
			if not self.do_CR:
				self.trigger_matched_medium = selections.RequireMediumTriggerMatching(
					self.tree,
					prompt_lepton_index=self.plep_sel.plep_index,
					prompt_lepton_type=self.plep,
					muons=muons,
					electrons=electrons,
					dv_type=self.dv_type)

				self.fill_ntuple(sel, 'n_trigger_matched_medium', self.trigger_matched_medium.n_trigger_matched_medium)
				# This will also be used for applying the lepton trigger matching scale factors

			if not self.do_CR and not self.tree.is_data:
				self.fill_systematics(sel)

			# ============================================================
			# Write values to ntuple. TTree already created and array already filled by fill_ntuple
			if self.MCEventType.isLNC:
				self.micro_ntuples["LNC_" + sel].fill()
				self.micro_ntuples["LNC_plus_LNV_" + sel].fill()
			elif self.MCEventType.isLNV:
				self.micro_ntuples["LNV_" + sel].fill()
				self.micro_ntuples["LNC_plus_LNV_" + sel].fill()
			else:
				if not self.tree.is_data: self.logger.warn("MC Ntuple being filled that is neither LNV or LNC. Please investigate this.")
				self.micro_ntuples[sel].fill()

			if sel == "sel":
				self._locked = FILL_LOCKED  # this only becomes unlocked after the event loop finishes in makeHistograms so you can only fill one DV from each event.

	def fill_systematics(self, sel):
		"""
		Calculate systematic variations store them to micro-ntuples
		"""

		# ============================================================
		# Systematics
		# ____________________________________________________________
		# Calculate weight for tracking/vertexing uncertainties
		vertexing_systematic = systematics.get_vertexing_systematic(self.tree.dv('r'), self.tree.dv('pt'))
		self.fill_ntuple(sel, 'vertexing_1DOWN', vertexing_systematic, weight=1)

		# ____________________________________________________________
		# Calculate weight for displaced lepton identification uncertainties
		lepton_0_type, lepton_1_type = '', ''
		if self.dv_type == "mumu": lepton_0_type, lepton_1_type = 'muon', 'muon'
		if self.dv_type == "emu": lepton_0_type, lepton_1_type = 'muon', 'electron'
		if self.dv_type == "ee": lepton_0_type, lepton_1_type = 'electron', 'electron'

		d0_extrapolation_systematic = systematics.get_combined_d0_extrapolation_systematic(
			self.tree.dv('trk_d0')[0], lepton_0_type,
			self.tree.dv('trk_d0')[1], lepton_1_type
		)
		# storing systematic with standard ATLAS "1-sigma down" notation
		self.fill_ntuple(sel, 'd0_extrapolation_1DOWN', d0_extrapolation_systematic, weight=1)

	def AddExtraVariables(self,sel):
		"""
		GUGLIELMO :: Function to add jet variables, we might put a flag to not call it by default.
		"""
		for jet_index in range(len(self.tree['jet_pt'])):
			self.jetVariables['pt'].push_back(self.tree['jet_pt'][jet_index])
			self.jetVariables['eta'].push_back(self.tree['jet_eta'][jet_index])
			self.jetVariables['phi'].push_back(self.tree['jet_phi'][jet_index])
			self.jetVariables['E'].push_back(self.tree['jet_E'][jet_index])
			self.jetVariables['DL1dv00'].push_back(self.tree['jet_DL1dv00'][jet_index])
			self.jetVariables['DL1dv01'].push_back(self.tree['jet_DL1dv01'][jet_index])
			if (self.tree['jet_GN1'][jet_index] < -999. ): self.jetVariables['GN1'].push_back(-999.)
			else: self.jetVariables['GN1'].push_back(self.tree['jet_GN1'][jet_index])
   
		self.fill_ntuple(sel, 'jet_pt', self.jetVariables['pt'])
		self.fill_ntuple(sel, 'jet_eta', self.jetVariables['eta'])
		self.fill_ntuple(sel, 'jet_phi', self.jetVariables['phi'])
		self.fill_ntuple(sel, 'jet_E', self.jetVariables['E'])
		self.fill_ntuple(sel, 'jet_DL1dv00', self.jetVariables['DL1dv00'])
		self.fill_ntuple(sel, 'jet_DL1dv01', self.jetVariables['DL1dv01'])
		self.fill_ntuple(sel, 'jet_GN1', self.jetVariables['GN1'])
	
#		if self.MCEventType.isLNC: 
#			self.micro_ntuples["LNC_"+sel].fill()
#			self.micro_ntuples["LNC_plus_LNV_"+sel].fill()
#		elif self.MCEventType.isLNV: 
#			self.micro_ntuples["LNV_"+sel].fill()
#			self.micro_ntuples["LNC_plus_LNV_"+sel].fill()
#		else: self.micro_ntuples[sel].fill()
		
		self.jetVariables['pt'].clear()
		self.jetVariables['eta'].clear()
		self.jetVariables['phi'].clear()
		self.jetVariables['E'].clear()
		self.jetVariables['DL1dv00'].clear()
		self.jetVariables['DL1dv01'].clear()
		self.jetVariables['GN1'].clear()


class run2Analysis(Analysis):
	def __init__(self, name, tree, vtx_container, selections, output_file, save_ntuples, weight_override=None,use_truth=True):

		Analysis.__init__(self, name, tree, vtx_container, selections, output_file, save_ntuples, weight_override)
		self.logger.info('Running  Full Run 2 Analysis cuts')
		self.logger.info('save_ntuples = {}'.format(self.save_ntuples))

		# Define cutflow histogram "by hand"
		self.cutflow_dir = 'CutFlow/'
		if not self.dv_type == "ee":
			self.observables.histogram_dict[self.cutflow_dir+ 'CutFlow'] = ROOT.TH1D('CutFlow_raw_counts', 'CutFlow_raw_counts', 21, -0.5, 20.5)
		else:
			self.observables.histogram_dict[self.cutflow_dir+ 'CutFlow'] = ROOT.TH1D('CutFlow_raw_counts', 'CutFlow_raw_counts', 22, -0.5, 21.5)
		self.CutFlow = self.observables.histogram_dict[self.cutflow_dir + 'CutFlow']
		# Bin labels are 1 greater than histogram bins
		self.CutFlow.GetXaxis().SetBinLabel(1, "all")
		if self.do_trigger_cut:
			if self.do_CR == False:
				self.CutFlow.GetXaxis().SetBinLabel(2, "trigger")
			else:
				self.CutFlow.GetXaxis().SetBinLabel(2, "DAOD_RPVLL triggers")
		if self.do_invert_trigger_cut:
			self.CutFlow.GetXaxis().SetBinLabel(2, "invert trigger")
		self.CutFlow.GetXaxis().SetBinLabel(3, "PV")
		if self.do_filter_cut:
			self.CutFlow.GetXaxis().SetBinLabel(4, "%s" % self.filter_type)
		if self.do_prompt_lepton_cut:
			self.CutFlow.GetXaxis().SetBinLabel(5, "{} prompt {}".format(self.plep_quality,self.plep))
			self.CutFlow.GetXaxis().SetBinLabel(6, "plep overlap")
		if self.do_invert_prompt_lepton_cut:
			self.CutFlow.GetXaxis().SetBinLabel(4, "invert prompt lepton")
			if self.do_prompt_track_cut:
				self.CutFlow.GetXaxis().SetBinLabel(5, "prompt track")
				self.CutFlow.GetXaxis().SetBinLabel(6, "no ptrk overlap with DV")	
		if self.do_ndv_cut:
			self.CutFlow.GetXaxis().SetBinLabel(7, "DV")
		if self.do_fidvol_cut:
			self.CutFlow.GetXaxis().SetBinLabel(8, "fiducial")
		if self.do_ntrk_cut:
			self.CutFlow.GetXaxis().SetBinLabel(9, "%s-track DV" % self.ntrk)
		if self.do_opposite_sign_cut:
			self.CutFlow.GetXaxis().SetBinLabel(10, "OS DV")
		if self.do_same_sign_cut:
			self.CutFlow.GetXaxis().SetBinLabel(10, "SS DV")
		if self.do_dv_type_cut:
			self.CutFlow.GetXaxis().SetBinLabel(11, "%s DV" % self.dv_type)
		if self.do_cosmic_veto_cut:
			self.CutFlow.GetXaxis().SetBinLabel(12, "cosmic veto")
		if self.do_lep_pt_cut:
			self.CutFlow.GetXaxis().SetBinLabel(13, "lepton pt")
		if self.do_mat_veto_cut:
			self.CutFlow.GetXaxis().SetBinLabel(14, "mat. veto")
		if self.do_track_quality_cut:
			if not self.dv_type == "ee": qual_bin = 14
			else: qual_bin = 15
			self.CutFlow.GetXaxis().SetBinLabel(qual_bin, "{}-lepton DV".format(self.track_quality))
		if self.do_trigger_matching_cut:
			if not self.dv_type == "ee": qual_bin = 15
			else: qual_bin = 16
			self.CutFlow.GetXaxis().SetBinLabel(qual_bin, "trig. match".format(self.track_quality))
		if self.do_trilepton_mass_cut:
			if not self.dv_type == "ee": mlll_bin = 16
			else: mlll_bin = 17
			self.CutFlow.GetXaxis().SetBinLabel(mlll_bin, "m_{lll}")
		if self.do_inverted_mlll_cut:
			if not self.dv_type == "ee": mlll_bin = 16
			else: mlll_bin = 17
			self.CutFlow.GetXaxis().SetBinLabel(mlll_bin, "inverted m_{lll}")
		if self.do_inverted_mhnl_cut:	
			if not self.dv_type == "ee": mhnl_bin = 16
			else: mhnl_bin = 17
			self.CutFlow.GetXaxis().SetBinLabel(mhnl_bin, "inverted m_{hnl}")
		if self.do_dv_mass_cut:
			if not self.dv_type == "ee": mdv_bin = 17
			else: mdv_bin = 18
			self.CutFlow.GetXaxis().SetBinLabel(mdv_bin, "B-hadron veto")
		if self.do_zmass_veto:
			if not self.dv_type == "ee": mdv_bin = 18
			else: mdv_bin = 19
			self.CutFlow.GetXaxis().SetBinLabel(mdv_bin, "Z mass veto")
		if self.do_HNL_mass_cut:
			if not self.dv_type == "ee": mhnl_bin = 19
			else: mhnl_bin = 20
			self.CutFlow.GetXaxis().SetBinLabel(mhnl_bin, "m_{HNL}")
		if self.do_alpha_cut:
			if not self.dv_type == "ee": alpha_bin = 20
			else: alpha_bin = 21
			self.CutFlow.GetXaxis().SetBinLabel(alpha_bin, "alpha")
		if use_truth and not self.tree.is_data:

			if not self.dv_type == "ee": truth_match_bin = 21
			else: truth_match_bin = 22
			self.CutFlow.GetXaxis().SetBinLabel(truth_match_bin, "truth matched")
		
		


		# Store LNC and LNV cutflows in the observables collection
		if not self.tree.is_data and not self.tree.is_bkg_mc:
			# #########################################################################################################
			# Raw cutflows
			# #########################################################################################################
			self.CutFlow_LNC = self.CutFlow.Clone()
			self.CutFlow_LNC.SetName("CutFlow_one_hnl_dirac"+"_"+self.ch)
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_LNC_raw_counts'] = self.CutFlow_LNC

			self.CutFlow_LNV = self.CutFlow.Clone()
			self.CutFlow_LNV.SetName("CutFlow_LNV_only"+"_"+self.ch)
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_LNV_raw_counts'] = self.CutFlow_LNV
			
			# #########################################################################################################
			# Weighted cutflows
			# #########################################################################################################

			# One HNL Majorana model (combined LNC and LNV decays)
			self.CutFlow_weighted_one_hnl_majorana = self.CutFlow.Clone()
			self.CutFlow_weighted_one_hnl_majorana.SetName("CutFlow_weighted_one_hnl_majorana"+"_"+self.ch)
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_weighted_one_hnl_majorana'] = self.CutFlow_weighted_one_hnl_majorana
			# LNC only decays for one Majorana HNL model (for comparing LNC and LNV acceptance)
			self.CutFlow_LNV_weighted = self.CutFlow.Clone()
			self.CutFlow_LNV_weighted.SetName("CutFlow_weighted_one_hnl_majorana_LNV_only"+"_"+self.ch)
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_weighted_one_hnl_majorana_LNV_only'] = self.CutFlow_LNV_weighted
			# LNC only decays for one Majorana HNL model (for comparing LNC and LNV acceptance)
			self.CutFlow_LNC_weighted = self.CutFlow.Clone()
			self.CutFlow_LNC_weighted.SetName("CutFlow_weighted_one_hnl_majorana_LNC_only"+"_"+self.ch)
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_weighted_one_hnl_majorana_LNC_only'] = self.CutFlow_LNC_weighted

			# One HNL Dirac model (LNC only decays)
			self.CutFlow_weighted_one_hnl_dirac = self.CutFlow.Clone()
			self.CutFlow_weighted_one_hnl_dirac.SetName("CutFlow_weighted_one_hnl_dirac"+"_"+self.ch)
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_weighted_one_hnl_dirac'] = self.CutFlow_weighted_one_hnl_dirac

			# Quasi-Dirac pair "Majorana limit" with inverted hierarchy mixing
			self.CutFlow_weighted_majorana_limit_ih = self.CutFlow.Clone()
			self.CutFlow_weighted_majorana_limit_ih.SetName("CutFlow_weighted_majorana_limit_ih"+"_"+self.ch)
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_weighted_majorana_limit_ih'] = self.CutFlow_weighted_majorana_limit_ih

			# Quasi-Dirac pair "Majorana limit" with normal hierarchy mixing
			self.CutFlow_weighted_majorana_limit_nh = self.CutFlow.Clone()
			self.CutFlow_weighted_majorana_limit_nh.SetName("CutFlow_weighted_majorana_limit_nh"+"_"+self.ch)
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_weighted_majorana_limit_nh'] = self.CutFlow_weighted_majorana_limit_nh

			# Quasi-Dirac pair "Dirac limit" with inverted hierarchy mixing
			self.CutFlow_weighted_dirac_limit_ih = self.CutFlow.Clone()
			self.CutFlow_weighted_dirac_limit_ih.SetName("CutFlow_weighted_dirac_limit_ih"+"_"+self.ch)
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_weighted_dirac_limit_ih'] = self.CutFlow_weighted_dirac_limit_ih

			# Quasi-Dirac pair "Dirac limit" with normal hierarchy mixing
			self.CutFlow_weighted_dirac_limit_nh = self.CutFlow.Clone()
			self.CutFlow_weighted_dirac_limit_nh.SetName("CutFlow_weighted_dirac_limit_nh"+"_"+self.ch)
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_weighted_dirac_limit_nh'] = self.CutFlow_weighted_dirac_limit_nh
		else:
			# Weighted cutflow
			self.CutFlow_weighted = self.CutFlow.Clone()
			self.CutFlow_weighted.SetName("CutFlow_weighted"+"_"+self.ch)
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_weighted'] = self.CutFlow_weighted


	def DVSelection(self,use_truth=False):
		######################################################################################################
		# DV Selection is any cuts that are done per DV
		# Current cuts include: fiducial vol, ntrack, OS, DVtype, track quality, cosmic veto, mlll, mDV
		######################################################################################################
		# reset lepton systematics for each DV
		for key in self.lepton_reco_sf.keys():
			self.lepton_reco_sf[key] = 1
		for key in self.lepton_trig_sf.keys():
			self.lepton_trig_sf[key] = 1

		# Fill all the histograms with ALL DVs (this could be more that 1 per event). Useful for vertexing efficiency studies.
		if self.save_ntuples == 'all':
			self._fill_all_dv_histos()

		# We only want to save one DV per event: Choose first DV that passes "DV type" selection as "the DV for the event"
		# If we have already selected our DV for this event, don't bother looking at the rest of the DVs!
		if self.selected_dv_index != -1: return

		# only do the DV selection if the preselction was passed for the event.
		if not self.passed_preselection_cuts:
			return

		# There is an extra bit of logic here since we look at several DVs
		# Do we want to use this cut?
		if self.do_fidvol_cut:
			# Does this cut pass?
			if self._fidvol_cut():
				# Has the cutflow already been filled for this event?
				if not self.passed_fidvol_cut:
					self._fill_cutflow(7)
					self.passed_fidvol_cut = True
			# If this cut doesn't pass, don't continue to check other cuts
			else:
				return

		if self.do_ntrk_cut:
			if self._ntrk_cut():
				if not self.passed_ntrk_cut:
					self._fill_cutflow(8)
					self.passed_ntrk_cut = True
			else:
				return

		if self.do_opposite_sign_cut or self.do_same_sign_cut:
			if self._charge_cut():
				if not self.passed_charge_cut:
					self._fill_cutflow(9)
					self.passed_charge_cut = True
			else:
				return

		if self.do_dv_type_cut:
			if self._dv_type_cut():
				if not self.passed_dv_type_cut:
					self._fill_cutflow(10)
					self.passed_dv_type_cut = True
					if not self.do_inverted_mlll_cut and not self.do_inverted_mhnl_cut:
						if self.save_ntuples == "DVtype":
							self._fill_selected_dv_histos(self.save_ntuples)
					# Select this DV as the DV for the event!
					self.selected_dv_index = self.tree.idv
			else:
				return

		if self.do_cosmic_veto_cut:
			if self._cosmic_veto_cut():
				if not self.passed_cosmic_veto_cut:
					self._fill_cutflow(11)
					self.passed_cosmic_veto_cut = True
					if not self.do_inverted_mlll_cut and not self.do_inverted_mhnl_cut:
						if self.save_ntuples == "cosmic":
							self._fill_selected_dv_histos(self.save_ntuples)
			else:
				return

		if self.do_lep_pt_cut:
			if self._dv_lep_pt_cut():
				if not self.passed_dlep_pt_cut:
					self._fill_cutflow(12)
					self.passed_lep_pt_cut = True
					if not self.do_inverted_mlll_cut and not self.do_inverted_mhnl_cut:
						if self.save_ntuples == "DVlep_pt":
							self._fill_selected_dv_histos(self.save_ntuples)
			else:
				return

		if self.do_mat_veto_cut:
			if self._mat_veto_cut():
				if not self.passed_mat_veto_cut:
					self._fill_cutflow(13)
					self.passed_mat_veto_cut = True
					if not self.do_inverted_mlll_cut and not self.do_inverted_mhnl_cut:
						if self.save_ntuples == "matveto":
							self._fill_selected_dv_histos(self.save_ntuples)
			else:
				return

		if self.do_track_quality_cut:
			if self._track_quality_cut():
				if not self.passed_track_quality_cut:

					# comment out to run on bw99 sample 
					
					# if not self.do_CR and not self.tree.is_data:
					# 	# update event weight
					# 	for systematic in self.lepton_reco_sf.keys():
					# 		self.lepton_reco_sf[systematic] *= scale_factors.get_reco_scale_factor(self, systematic)

					if not self.dv_type == "ee": self._fill_cutflow(13)
					else: self._fill_cutflow(13+1)
					self.passed_track_quality_cut = True
					if not self.do_inverted_mlll_cut and not self.do_inverted_mhnl_cut:
						if self.save_ntuples == "trkqual":
							self._fill_selected_dv_histos(self.save_ntuples)
			else:
				return

		if self.do_trigger_matching_cut:
			if self._trigger_matched_medium_lepton_cut():
				if not self.passed_trig_matched_cut:

					# if not self.do_CR and not self.tree.is_data:
						# update event weight TODO: fix this (Christian)
						# for systematic in self.lepton_trig_sf.keys():
						# 	self.lepton_trig_sf[systematic] *= scale_factors.get_trigger_scale_factor(self, systematic)

					if not self.dv_type == "ee": self._fill_cutflow(14)
					else: self._fill_cutflow(14+1)
					self.passed_trig_matched_cut = True
					if not self.do_inverted_mlll_cut and not self.do_inverted_mhnl_cut:
						if self.save_ntuples == "trig_match":
							self._fill_selected_dv_histos(self.save_ntuples)
			else:
				return

		if self.do_trilepton_mass_cut:
			if self._trilepton_mass_cut():
				if not self.passed_trilepton_mass_cut:
					if not self.dv_type == "ee": self._fill_cutflow(15)
					else: self._fill_cutflow(15+1)
					self.passed_trilepton_mass_cut = True
					if self.save_ntuples == "mvis":
						self._fill_selected_dv_histos(self.save_ntuples)
			else:
				return
		
		if self.do_inverted_mlll_cut:
			if self._invert_trilepton_mass_cut():
				if not self.passed_trilepton_mass_cut:
					if not self.dv_type == "ee": self._fill_cutflow(15)
					else: self._fill_cutflow(15+1)
					self.passed_trilepton_mass_cut = True
					if self.save_ntuples == "mvis":
						self._fill_selected_dv_histos(self.save_ntuples)
			else:
				return
		
		if self.do_inverted_mhnl_cut:
			if self._invert_HNL_mass_cut():
				if not self.passed_HNL_mass_cut:
					if not self.dv_type == "ee": self._fill_cutflow(15)
					else: self._fill_cutflow(15+1)
					self.passed_HNL_mass_cut = True
					if self.save_ntuples == "mhnl":
						self._fill_selected_dv_histos(self.save_ntuples)
			else:
				return

		if self.do_dv_mass_cut:
			# do b-hadron veto instead of just straight DV mass cut
			if self._Bhadron_veto():
				if not self.passed_dv_mass_cut:
					if not self.dv_type == "ee": self._fill_cutflow(16)
					else: self._fill_cutflow(16+1)
					self.passed_dv_mass_cut = True
					if self.save_ntuples == "mDV":
						self._fill_selected_dv_histos(self.save_ntuples)
			else:
				return

		if self.do_zmass_veto:
			# do z mass veto for all channels
			if self._Zmass_veto():
				if not self.passed_zmass_veto:
					if not self.dv_type == "ee": self._fill_cutflow(17)
					else: self._fill_cutflow(17+1)
					self.passed_zmass_veto = True
					if self.save_ntuples == "Zmass_veto":
						self._fill_selected_dv_histos(self.save_ntuples)
			else:
				return

		if self.do_HNL_mass_cut :
			if self._HNL_mass_cut():
				if not self.passed_HNL_mass_cut:
					if not self.dv_type == "ee": self._fill_cutflow(18)
					else: self._fill_cutflow(18+1)
					self.passed_HNL_mass_cut = True
					if self.save_ntuples == "mHNL":
						self._fill_selected_dv_histos(self.save_ntuples)
			else:
				return

		if self.do_alpha_cut:
			if self._alpha_cut():
				if not self.passed_alpha_cut:
					if not self.dv_type == "ee": self._fill_cutflow(19)
					else: self._fill_cutflow(19+1)
					self.passed_alpha_cut = True
					if self.save_ntuples == "alpha":
						self._fill_selected_dv_histos(self.save_ntuples)
			else:
				return

		# Fill histos of truth-matched DVs
		if use_truth and not self.tree.is_data and not self.tree.is_bkg_mc:
			if self._truth_match():
				if not self.dv_type == "ee": self._fill_cutflow(20)
				else: self._fill_cutflow(20+1)
				if self.save_ntuples == "match":
					self._fill_selected_dv_histos(self.save_ntuples)
