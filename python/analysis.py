import ROOT
import numpy as np
import os
import sys
import helpers
import selections
import observables
import ntuples
import systematics

UNLOCKED = 0
SELECTION_LOCKED = 1
FILL_LOCKED = 2


class Analysis(object):
	def __init__(self, name, tree, vtx_container, selection_list, outputFile, saveNtuples, weight_override=None):
		# set up logger for self
		self.logger = helpers.getLogger('dHNLAnalysis.analysis', level=helpers.logger_debug_level)
		# set up logger for helper module
		selections.logger.setLevel(helpers.logger_debug_level)

		self.name = name
		self.weight_override = weight_override
		self.sel = selection_list
		self.outputFile = outputFile
		self.fi = ROOT.TFile.Open(outputFile, 'update')
		self.ch = vtx_container
		self.saveNtuples = saveNtuples
		self.histSuffixes = [self.ch]
		self.h = {}
		self.micro_ntuples = {}
		self.tree = tree
		self._locked = UNLOCKED
		# create an instance of Observables to store histograms
		self.observables = observables.Observables()

		self.events_with_trig_match_plep = 0
		self.events_with_trig_match_dlep = 0
		self.events_with_trig_match_both_pdlep =0

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
		self.do_filter_cut = True
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
			self.saveNtuples = "mvis" # only save ntuples after mlll selection is applied!
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
			self.saveNtuples = "mhnl" # only save ntuples after mlll selection is applied!
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
		return self.tree['secVtx_{}_{}'.format(self.ch, key)]

	def get(self, key):
		return self.tree[key]

	# hist filling helper functions
	def fill_hist(self, selection, hist_name, variable_1, variable_2=None, fill_ntuple=True, weight=None):
		"""
		A helper function for filling registered histograms
		:param selection: the step of selection the analysis it at. May be "None" in which case there will be no prefix.
		:param hist_name: base name of the histogram. When saved, a prefix and suffix will be appended.
		:param variable_1: variable you want to fill the histogram with.
		:param variable_2: if histogram is 2d, variable you want to fill the second axis of the histogram with
		:param fill_ntuple: set to True if you want to simultaneously fill an ntuple with this variable
		:param weight: Used to override weight calculation. Useful for storing actual weights which should not be themselves weighted.
		"""
		# use calculated weight for this event unless a weight is specified
		weight_LNC_only = self.weight_LNC_only if weight is None else weight
		weight_LNC_plus_LNV = self.weight_LNC_plus_LNV if weight is None else weight

		# define here the directory structure where this histogram is stored.
		directory = '{ch}/{selection}/'.format(ch=self.ch, selection=selection)
		# fill LNC histograms. Not used for data.
		if self.MCEventType.isLNC: 
			self.observables.fill_hist(directory+'LNC/', hist_name, variable_1, variable_2, weight_LNC_only)
		# fill LNV histograms. Not used for data.
		if self.MCEventType.isLNV: 
			self.observables.fill_hist(directory+'LNV/', hist_name, variable_1, variable_2, weight_LNC_only)
		if self.MCEventType.isLNC or self.MCEventType.isLNV: 
			#fill LNC_plus_LNV histograms for every event
			# U, m, ctau relationship changes becuase twice as many decay channels avaliable if HNL can decay via LNC and LNV
			self.observables.fill_hist(directory+'LNC_plus_LNV/', hist_name, variable_1, variable_2, weight_LNC_plus_LNV)

		# Unless suppressed, fill the corresponding micro-ntuple with the variable
		# Will not fill variables from 2D histograms to prevent double-counting
		# TODO Can we clean this up in some way?
		save_sel = self.saveNtuples == selection or 'truth_'+self.saveNtuples == selection or self.saveNtuples == 'allcuts'
		if fill_ntuple and (variable_2 is None) and save_sel:
			# Need selection to define ntuple tree
			# TODO redo this method to use the directory correctly
			if self.MCEventType.isLNC: 
				self.fill_ntuple(selection, hist_name, variable_1, mc_type="LNC")
				self.fill_ntuple(selection, hist_name, variable_1, mc_type="LNC_plus_LNV")
			elif self.MCEventType.isLNV: 
				self.fill_ntuple(selection, hist_name, variable_1, mc_type="LNV")
				self.fill_ntuple(selection, hist_name, variable_1, mc_type="LNC_plus_LNV")
			else: self.fill_ntuple(selection, hist_name, variable_1)


	def fill_ntuple(self, selection, ntuple_name, variable, mc_type=None, full_name=""):
		"""
		A helper function for filling micro-ntuples. Often called from the fill_hist function.
		If you are using this in your analysis,
		please check that it is not also being called by fill_hist to prevent double-counting.
		:param selection: the step of selection the analysis it at. May be "None" in which case there will be no prefix.
		:param ntuple_name: base name of the ntuple. When saved, a prefix and suffix will be appended.
		:param variable: variable you want to fill the histogram with.
		:param mc_type: defines whether the vertex is lepton-number conserving (LNC) or violating (LNV)
		:param full_name: override the automatic naming of the ntuple.
		"""
		if not selection:
			raise ValueError("You must indicate a selection in order to store the ntuple. Use 'all' if no selection.")
		if mc_type is not None:
			selection = mc_type + "_" + selection

		# Retrieve the ntuple for this selection. If it doesn't exist, create it.
		if selection not in self.micro_ntuples:
			# temp name. not written
			self.micro_ntuples[selection] = ntuples.Ntuples('ntuples_{}_{}'.format(selection, self.ch))
		# The name of the ntuple
		if not full_name:
			full_name = ntuple_name
		self.micro_ntuples[selection][full_name] = variable

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
		# Move ROOT to base directory
		self.fi.cd()
		[ntuple.write(self.ch+'_ntuples_'+key) for key, ntuple in self.micro_ntuples.items()]
		self.observables.write_histograms(root_file=self.fi)
		self.logger.info("Histograms written to {}".format(self.outputFile))

		self.fi.Close()

	def end(self):
		# make acceptance Histograms
		# for hist in self.observables.histogram_dict.values():
		# 	hist.SetBinContent(hist.GetNbinsX(), hist.GetBinContent(hist.GetNbinsX()) + hist.GetBinContent(hist.GetNbinsX() + 1)) # merge overflow into last bin
		# 	hist.SetBinContent(1, hist.GetBinContent(1) + hist.GetBinContent(0)) # merge underflow into first bin

		# make acceptance Histograms
		# TOD: it doesnt looks like on data the acceptance histograms are working as expected. -DT
		if not self.tree.is_data and not self.tree.not_hnl_mc:
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_LNV_acceptance'] = self.CutFlow_LNV.Clone()
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_LNC_acceptance'] = self.CutFlow_LNC.Clone()
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_LNV_acceptance'].SetName("CutFlow_LNV_acceptance"+"_"+self.ch)
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_LNV_acceptance'].SetName("CutFlow_LNV_acceptance"+"_"+self.ch)
			if self.CutFlow_LNV.GetBinContent(1) != 0: # Protect against zero-division
				self.observables.histogram_dict[self.cutflow_dir+'CutFlow_LNV_acceptance'].Scale(1.0/self.CutFlow_LNV.GetBinContent(1))
			if self.CutFlow_LNC.GetBinContent(1) != 0: # Protect against zero-division
				self.observables.histogram_dict[self.cutflow_dir+'CutFlow_LNC_acceptance'].Scale(1.0/self.CutFlow_LNC.GetBinContent(1))

		self.observables.histogram_dict[self.cutflow_dir+'CutFlow_acceptance'] = self.CutFlow.Clone()
		self.observables.histogram_dict[self.cutflow_dir+'CutFlow_acceptance'].SetName("CutFlow_acceptance"+"_"+self.ch)
		if self.CutFlow.GetBinContent(1) != 0: # Protect against zero-division
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_acceptance'].Scale(1.0/self.CutFlow.GetBinContent(1))
		
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

	def _prompt_lepton_cut(self):
		self.found_plep = False # intitalize the plep each event 
		self.plep_sel = selections.PromptLepton(self.tree, lepton=self.plep,quality=self.plep_quality) # run plep selection 
		self.found_plep = self.plep_sel.found_plep # check if you found any prompt leptons 
		# Add to histogram all prompt leptons that pass selection.
		# If _prompt_lepton_cut() is run after trigger and filter cut then those cuts will also be applied.
		if self.plep_sel.passes():
			# trig_match = selections.Lep_TriggerMatching(self.tree, self.plep, self.plep_sel.plep_Index)
			# if trig_match.lep_isTrigMatched:
			# 	self.events_with_trig_match_plep = self.events_with_trig_match_plep + 1

			# fill prompt lepton histograms
			self.fill_hist('all', 'plep_pt', self.plep_sel.plepVec.Pt())
			self.fill_hist('all', 'plep_eta', self.plep_sel.plepVec.Eta())
			self.fill_hist('all', 'plep_phi', self.plep_sel.plepVec.Phi())
			self.fill_hist('all', 'plep_d0', self.plep_sel.plepd0)
			self.fill_hist('all', 'plep_z0', self.plep_sel.plepz0)
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
			prompt_lepton_index=self.plep_sel.plep_Index,
			prompt_lepton_type=self.plep,
			muons=muons,
			electrons=electrons,
			dv_type=self.dv_type)

		return self.trigger_matched_medium.passes()

	def _prompt_track_cut(self):
		self.found_ptrk = False # intitalize the plep each event
		self.ptrk_sel = selections.PromptTrack(self.tree) # run ptrack selection
		self.found_ptrk = self.ptrk_sel.found_trk # check if you found any prompt leptons
		# Add to histogram all prompt leptons that pass selection.
		# If _prompt_lepton_cut() is run after trigger and filter cut then those cuts will also be applied.
		if self.ptrk_sel.passes():
			self.fill_hist('all', 'ptrk_pt', self.ptrk_sel.trkVec.Pt())
			self.fill_hist('all', 'ptrk_eta', self.ptrk_sel.trkVec.Eta())
			self.fill_hist('all', 'ptrk_phi', self.ptrk_sel.trkVec.Phi())
			self.fill_hist('all', 'ptrk_d0', self.ptrk_sel.trkd0)
			self.fill_hist('all', 'ptrk_z0', self.ptrk_sel.trkz0)


		return self.ptrk_sel.passes() # full plep selection find the highest pt plep that doesnt overlap with any DVs

	def _invert_prompt_lepton_cut(self):
		self.invt_lep = selections.InvertedPromptLepton(self.tree)
		return self.invt_lep.passes()

	def _ndv_cut(self):
		return self.tree.ndv > 0

	def _fidvol_cut(self):
		fidvol_sel = selections.DVradius(self.tree)
		return fidvol_sel.passes()

	def _ntrk_cut(self):
		ntracks_sel = selections.DVntracks(self.tree, ntrk=self.ntrk)
		return ntracks_sel.passes()

	def _charge_cut(self):
		sign_pair = "SS" if self.do_same_sign_cut else "OS"
		charge_sel = selections.ChargeDV(self.tree, sel=sign_pair)
		return charge_sel.passes()
	
	def _be_event_type_cut(self): 
		if self.do_same_event_cut: return not self.tree.dv('shuffled')
		if self.do_different_event_cut: return self.tree.dv('shuffled')
		
	def _dv_type_cut(self):
		dv_sel = selections.DVtype(self.tree, dv_type=self.dv_type)
		# if dv_sel.passes():
		# 	if not self.tree.fake_aod:
		# 		trig_match = selections.TriggerMatching_disp(self.tree, self.dv_type, dv_sel.dMu_Index, dv_sel.dEl_Index)
		# 		if trig_match.dlep_isTrigMatched:
		# 			self.events_with_trig_match_dlep = self.events_with_trig_match_dlep + 1
		# 		count_trig_match_disp_event = True

		return dv_sel.passes()

	def _track_quality_cut(self):
		track_quality_sel = selections.Trackqual(self.tree, quality=self.track_quality)
		return track_quality_sel.passes()

	def _cosmic_veto_cut(self):
		cosmic_veto_sel = selections.Cosmicveto(self.tree)
		# self.h["DV_trk_sep"][self.ch].Fill(cosmic_veto_sel.separation)
		return cosmic_veto_sel.passes()

	def _trilepton_mass_cut(self):
		plep_vec = self.plep_sel.plepVec

		muons = helpers.Tracks(self.tree)
		muons.getMuons()
		muVec = muons.lepVec

		electrons = helpers.Tracks(self.tree)
		electrons.getElectrons()
		elVec = electrons.lepVec

		mlll_sel = selections.Mlll(dv_type=self.dv_type, plep=plep_vec, dMu=muVec, dEl=elVec,minmlll=40,maxmlll=90)
		return mlll_sel.passes()
	
	def _invert_trilepton_mass_cut(self):
		plep_vec = self.plep_sel.plepVec

		muons = helpers.Tracks(self.tree)
		muons.getMuons()
		muVec = muons.lepVec

		electrons = helpers.Tracks(self.tree)
		electrons.getElectrons()
		elVec = electrons.lepVec

		mlll_sel = selections.Mlll(dv_type=self.dv_type, plep=plep_vec, dMu=muVec, dEl=elVec,minmlll=25,maxmlll=125, invert=True)
		return mlll_sel.passes()

	def _mat_veto_cut(self):
		return selections.Mat_veto(self.tree).passes()

	def _dv_mass_cut(self):
		dv_mass_sel = selections.DVmass(self.tree, dvmasscut=5.5)
		return dv_mass_sel.passes()
	
	def _Bhadron_veto(self):
		bhadron_sel = selections.Bhadron_veto( self.tree, dv_type=self.dv_type, dvmasscut=5.5)
		return bhadron_sel.passes()

	def _Zmass_veto(self):
		zmass_veto = selections.Zmass_veto(self.tree, plep_vec = self.plep_sel.plepVec, plep=self.plep, plepcharge = self.plep_sel.plepcharge, dv_type= self.dv_type)
		return zmass_veto.passes()

	def _dv_lep_pt_cut(self):
		dv_lep_pt_sel = selections.DV_lep_pt(self.tree,self.dv_type)
		return dv_lep_pt_sel.passes()

	def _alpha_cut(self):
		alpha_sel = selections.Alpha(self.tree)
		return alpha_sel.passes()
	
	def _HNL_mass_cut(self):
		muons = helpers.Tracks(self.tree)
		muons.getMuons()
		muVec = muons.lepVec

		electrons = helpers.Tracks(self.tree)
		electrons.getElectrons()
		elVec = electrons.lepVec

		mHNL_sel = selections.Mhnl(self.tree, self.dv_type, plep=self.plep_sel.plepVec, dMu=muVec,dEl=elVec)
		return mHNL_sel.passes()

	def _invert_HNL_mass_cut(self):
		muons = helpers.Tracks(self.tree)
		muons.getMuons()
		muVec = muons.lepVec

		electrons = helpers.Tracks(self.tree)
		electrons.getElectrons()
		elVec = electrons.lepVec

		mHNL_sel = selections.Mhnl(self.tree, self.dv_type, plep=self.plep_sel.plepVec, dMu=muVec,dEl=elVec,hnlmasscut=20, invert= True)
		return mHNL_sel.passes()

	def _multitrk_2lep_cut(self):
		if self.tree.dv('ntrk') >= 2:  # 2+ trk vertex
			dv_type_sel = selections.DVtype(self.tree, dv_type=self.dv_type)
			if dv_type_sel.passes():  # 2 leptons in the DV
				sign_pair = "SS" if self.do_same_sign_cut else "OS"
				charge_sel = selections.ChargeDV(self.tree, sel=sign_pair, trk_charge=dv_type_sel.lepton_charge)
				return charge_sel.passes()

	def _truth_match(self):
		"""Truth matching function for displaced vertices.
		linkTruth_score is calculated using DVAnalysisBase.
		pdgId 50 signifies a heavy neutral lepton parent particle."""
		maxlinkTruth_score = self.tree.dv('maxlinkTruth_score')
		maxlinkTruth_parent_pdgId = abs(self.tree.dv('maxlinkTruth_parent_pdgId'))
		return self.tree.dv('maxlinkTruth_score') > 0.75 and abs(self.tree.dv('maxlinkTruth_parent_pdgId')) == 50

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


	def preSelection(self):
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

		if not self.tree.fake_aod:
			self._fill_leptons()

		if not self.tree.is_data and not self.tree.not_hnl_mc:
			self._fill_truth_histos(sel='truth/all')
		# 	if self.MCEventType.isLNC:
		# 		self.CutFlow_LNC.SetBinContent(1, self.tree.all_entries/2)  # all events
		# 	if self.MCEventType.isLNV:
		# 		self.CutFlow_LNV.SetBinContent(1, self.tree.all_entries/2)  # all events
		# self.CutFlow.SetBinContent(1, self.tree.all_entries)  # all events
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
		if not self.tree.is_data and not self.tree.not_hnl_mc:
			self._fill_truth_histos(sel='truth/presel')

	def calculate_event_weight(self):
		# MC re-weighting to include spin correlations and fix lepton ordering bug
		self.MCEventType = selections.MCEventType(self.tree) # if data then MCEventType weight defaults to 1

		# calculate mass lifetime weight 
		self.mass_lt_weight_LNC_only = helpers.get_mass_lt_weight(self.tree, lnc_plus_lnv=False)
		self.mass_lt_weight_LNC_plus_LNV = helpers.get_mass_lt_weight(self.tree, lnc_plus_lnv=True)

		# self.mass_lt_weight = helpers.get_mass_lt_weight(self.tree.mass, self.tree.ctau,lnv=self.MCEventType.isLNV)
		self.logger.debug('LNC only mass-lifetime weight for this signal sample is: {}'.format(self.mass_lt_weight_LNC_only))
		self.logger.debug('LNC+LNV mass-lifetime weight for this signal sample is: {}'.format(self.mass_lt_weight_LNC_plus_LNV))

		if self.weight_override is None:
			self.weight_LNC_only = self.mass_lt_weight_LNC_only * self.MCEventType.weight
			self.weight_LNC_plus_LNV = self.mass_lt_weight_LNC_plus_LNV * self.MCEventType.weight
		else:
			self.weight_LNC_only = self.weight_override
			self.weight_LNC_plus_LNV = self.weight_override

		# calculate the product of all scale factors
		self.scale_factor_product = 1
		try:
			self.scale_factor_product *= self.tree['weight_pileup']
		except KeyError:
			pass

		self.weight_LNC_only *= self.scale_factor_product
		self.weight_LNC_plus_LNV *= self.scale_factor_product

	def DVSelection(self):
		raise NotImplementedError("Please implement this method in your own Analysis subclass")

	def _fill_cutflow(self, nbin):
		if not self.tree.is_data and not self.tree.not_hnl_mc:
			if self.MCEventType.isLNC:
				self.CutFlow_LNC.Fill(nbin)
				self.CutFlow_LNC_weighted.Fill(nbin, self.weight_LNC_only) # weight LNC only
			if self.MCEventType.isLNV:
				self.CutFlow_LNV.Fill(nbin)
				self.CutFlow_LNV_weighted.Fill(nbin, self.weight_LNC_only) # weight LNC only since LNV events are scaled to 100% LNV (do not include extra factor of 2 for LNC+LNV model)
			self.CutFlow.Fill(nbin)
			self.CutFlow_LNC_plus_LNV.Fill(nbin, self.weight_LNC_plus_LNV)
		else:
			self.CutFlow.Fill(nbin)

	def _fill_multitrk_histos(self):
		self.fill_hist('2lepMultitrk', 'num_trks', self.tree.dv('ntrk'))
		muons = helpers.Tracks(self.tree)
		muons.getMuons()
		if muons.lepisAssoc[0] == 1 and muons.lepisAssoc[1] == 1:
			self.fill_hist('2lepMultitrk', 'bothmuon_isAssociated', 1)
		else:
			self.fill_hist('2lepMultitrk', 'bothmuon_isAssociated', 0)
		if muons.lepisAssoc[0] == 0 and muons.lepisAssoc[1] == 0:
			self.fill_hist('2lepMultitrk', 'nomuon_isAssociated', 1)
			tracks = helpers.Tracks(self.tree)
			tracks.getTracks()
			trk_assoc = tracks.lepisAssoc
			num_trk_assoc  = sum(trk_assoc)
			ntrk = [3,4,5,6,7,8,9,10]
			for i in range(len(ntrk)):
				if len(tracks.lepisAssoc) == ntrk[i]:
					self.fill_hist('2lepMultitrk', 'num_assoc_{}trk'.format(ntrk[i]), num_trk_assoc)
			for i in range(self.tree.ntrk):
				if self.tree.dv('trk_muonIndex')[i] >= 0:
					self.fill_hist('2lepMultitrk', 'lep_trk_pt', self.tree.dv('trk_pt_wrtSV')[i])
					self.fill_hist('2lepMultitrk', 'lep_trk_eta', self.tree.dv('trk_eta_wrtSV')[i])
					self.fill_hist('2lepMultitrk', 'lep_trk_phi', self.tree.dv('trk_phi_wrtSV')[i])
					self.fill_hist('2lepMultitrk', 'lep_trk_d0', self.tree.dv('trk_d0')[i])
					self.fill_hist('2lepMultitrk', 'lep_trk_z0', self.tree.dv('trk_z0')[i])
					self.fill_hist('2lepMultitrk', 'lep_trk_absz0', abs(self.tree.dv('trk_z0')[i]))
					self.fill_hist('2lepMultitrk', 'lep_trk_charge', self.tree.dv('trk_charge')[i])
					self.fill_hist('2lepMultitrk', 'lep_trk_chi2', self.tree.dv('trk_chi2')[i])
					self.fill_hist('2lepMultitrk', 'lep_trk_isLRT', self.tree.dv('trk_isLRT')[i])
					self.fill_hist('2lepMultitrk', 'lep_trk_isSelected', self.tree.dv('trk_isSelected')[i])
					self.fill_hist('2lepMultitrk', 'lep_trk_isAssociated', self.tree.dv('trk_isAssociated')[i])
					self.fill_hist('2lepMultitrk', 'lep_trk_nPixelHits', self.tree.dv('trk_nPixelHits')[i])
					self.fill_hist('2lepMultitrk', 'lep_trk_nSCTHits', self.tree.dv('trk_nSCTHits')[i])
					self.fill_hist('2lepMultitrk', 'lep_trk_nSCTHoles', self.tree.dv('trk_nSCTHoles')[i])
					self.fill_hist('2lepMultitrk', 'lep_trk_nSiHits', self.tree.dv('trk_nSCTHits')[i]+self.tree.dv('trk_nPixelHits')[i])
					# self.fill_hist('2lepMultitrk', 'lep_trk_dTheta', self.tree.dv('trk_dTheta')[i])
					self.fill_hist('2lepMultitrk', 'lep_trk_chi2_toSV'.format(i), self.tree.dv('trk_chi2_toSV')[i])
					self.fill_hist('2lepMultitrk', 'lep_trk_d0_wrtSV'.format(i), self.tree.dv('trk_d0_wrtSV')[i])
					self.fill_hist('2lepMultitrk', 'lep_trk_errd0_wrtSV'.format(i), self.tree.dv('trk_errd0_wrtSV')[i])
					self.fill_hist('2lepMultitrk', 'lep_trk_z0_wrtSV'.format(i), self.tree.dv('trk_z0_wrtSV')[i])
					self.fill_hist('2lepMultitrk', 'lep_trk_errz0_wrtSV'.format(i), self.tree.dv('trk_errz0_wrtSV')[i])
				else:
					self.fill_hist('2lepMultitrk', 'nonlep_trk_pt', self.tree.dv('trk_pt_wrtSV')[i])
					self.fill_hist('2lepMultitrk', 'nonlep_trk_eta', self.tree.dv('trk_eta_wrtSV')[i])
					self.fill_hist('2lepMultitrk', 'nonlep_trk_phi', self.tree.dv('trk_phi_wrtSV')[i])
					self.fill_hist('2lepMultitrk', 'nonlep_trk_d0', self.tree.dv('trk_d0')[i])
					self.fill_hist('2lepMultitrk', 'nonlep_trk_z0', self.tree.dv('trk_z0')[i])
					self.fill_hist('2lepMultitrk', 'nonlep_trk_absz0', abs(self.tree.dv('trk_z0')[i]))
					self.fill_hist('2lepMultitrk', 'nonlep_trk_charge', self.tree.dv('trk_charge')[i])
					self.fill_hist('2lepMultitrk', 'nonlep_trk_chi2', self.tree.dv('trk_chi2')[i])
					self.fill_hist('2lepMultitrk', 'nonlep_trk_isLRT', self.tree.dv('trk_isLRT')[i])
					self.fill_hist('2lepMultitrk', 'nonlep_trk_isSelected', self.tree.dv('trk_isSelected')[i])
					self.fill_hist('2lepMultitrk', 'nonlep_trk_isAssociated', self.tree.dv('trk_isAssociated')[i])
					self.fill_hist('2lepMultitrk', 'nonlep_trk_nPixelHits', self.tree.dv('trk_nPixelHits')[i])
					self.fill_hist('2lepMultitrk', 'nonlep_trk_nSCTHits', self.tree.dv('trk_nSCTHits')[i])
					self.fill_hist('2lepMultitrk', 'nonlep_trk_nSCTHoles', self.tree.dv('trk_nSCTHoles')[i])
					self.fill_hist('2lepMultitrk', 'nonlep_trk_nSiHits', self.tree.dv('trk_nSCTHits')[i]+self.tree.dv('trk_nPixelHits')[i])
					# self.fill_hist('2lepMultitrk', 'nonlep_trk_dTheta', self.tree.dv('trk_dTheta')[i])
					self.fill_hist('2lepMultitrk', 'nonlep_trk_chi2_toSV'.format(i), self.tree.dv('trk_chi2_toSV')[i])
					self.fill_hist('2lepMultitrk', 'nonlep_trk_d0_wrtSV'.format(i), self.tree.dv('trk_d0_wrtSV')[i])
					self.fill_hist('2lepMultitrk', 'nonlep_trk_errd0_wrtSV'.format(i), self.tree.dv('trk_errd0_wrtSV')[i])
					self.fill_hist('2lepMultitrk', 'nonlep_trk_z0_wrtSV'.format(i), self.tree.dv('trk_z0_wrtSV')[i])
					self.fill_hist('2lepMultitrk', 'nonlep_trk_errz0_wrtSV'.format(i), self.tree.dv('trk_errz0_wrtSV')[i])

		else:
			self.fill_hist('2lepMultitrk', 'nomuon_isAssociated', 0)
		if muons.lepisAssoc[0] != muons.lepisAssoc[1]:
			self.fill_hist('2lepMultitrk', 'onemuon_isAssociated', 1)
		else:
			self.fill_hist('2lepMultitrk', 'onemuon_isAssociated', 0)
		for i in range(len(muons.lepisAssoc)):
			self.fill_hist('2lepMultitrk', 'muon_isAssociated', muons.lepisAssoc[i])

		for i in range(self.tree.ntrk):
			self.fill_hist('2lepMultitrk', 'all_trk_pt', self.tree.dv('trk_pt_wrtSV')[i])
			self.fill_hist('2lepMultitrk', 'all_trk_eta', self.tree.dv('trk_eta_wrtSV')[i])
			self.fill_hist('2lepMultitrk', 'all_trk_phi', self.tree.dv('trk_phi_wrtSV')[i])
			self.fill_hist('2lepMultitrk', 'all_trk_d0', self.tree.dv('trk_d0')[i])
			self.fill_hist('2lepMultitrk', 'all_trk_z0', self.tree.dv('trk_z0')[i])
			self.fill_hist('2lepMultitrk', 'all_trk_absz0', abs(self.tree.dv('trk_z0')[i]))
			self.fill_hist('2lepMultitrk', 'all_trk_charge', self.tree.dv('trk_charge')[i])
			self.fill_hist('2lepMultitrk', 'all_trk_chi2', self.tree.dv('trk_chi2')[i])
			self.fill_hist('2lepMultitrk', 'all_trk_isLRT', self.tree.dv('trk_isLRT')[i])
			self.fill_hist('2lepMultitrk', 'all_trk_isSelected', self.tree.dv('trk_isSelected')[i])
			self.fill_hist('2lepMultitrk', 'all_trk_isAssociated', self.tree.dv('trk_isAssociated')[i])
			self.fill_hist('2lepMultitrk', 'all_trk_nPixelHits', self.tree.dv('trk_nPixelHits')[i])
			self.fill_hist('2lepMultitrk', 'all_trk_nSCTHits', self.tree.dv('trk_nSCTHits')[i])
			self.fill_hist('2lepMultitrk', 'all_trk_nSCTHoles', self.tree.dv('trk_nSCTHoles')[i])
			self.fill_hist('2lepMultitrk', 'all_trk_nSiHits', self.tree.dv('trk_nSCTHits')[i]+self.tree.dv('trk_nPixelHits')[i])
			# self.fill_hist('2lepMultitrk', 'all_trk_dTheta', self.tree.dv('trk_dTheta')[i])
			self.fill_hist('2lepMultitrk', 'all_trk_chi2_toSV', self.tree.dv('trk_chi2_toSV')[i])
			self.fill_hist('2lepMultitrk', 'all_trk_d0_wrtSV', self.tree.dv('trk_d0_wrtSV')[i])
			self.fill_hist('2lepMultitrk', 'all_trk_errd0_wrtSV', self.tree.dv('trk_errd0_wrtSV')[i])
			self.fill_hist('2lepMultitrk', 'all_trk_z0_wrtSV', self.tree.dv('trk_z0_wrtSV')[i])
			self.fill_hist('2lepMultitrk', 'all_trk_errz0_wrtSV', self.tree.dv('trk_errz0_wrtSV')[i])

	def _fill_leptons(self):
		sel = 'all'
		for imu in range(len(self.tree['muon_type'])):
			self.fill_hist(sel, 'muon_type', self.tree['muon_type'][imu])
			self.fill_hist(sel, 'muon_pt', self.tree['muon_pt'][imu])
			self.fill_hist(sel, 'muon_eta', self.tree['muon_eta'][imu])
			self.fill_hist(sel, 'muon_phi', self.tree['muon_phi'][imu])
			if self.tree['muon_isTight'][imu] == 1:  self.fill_hist(sel, 'muon_quality', 3)
			if self.tree['muon_isMedium'][imu] == 1: self.fill_hist(sel, 'muon_quality', 2)
			if self.tree['muon_isLoose'][imu] == 1:  self.fill_hist(sel, 'muon_quality', 1)
			else: self.fill_hist(sel, 'muon_quality', 0)

		for iel in range(len(self.tree['el_pt'])):
			self.fill_hist(sel, 'el_pt', self.tree['el_pt'][iel])
			self.fill_hist(sel, 'el_eta', self.tree['el_eta'][iel])
			self.fill_hist(sel, 'el_phi', self.tree['el_phi'][iel])
			if self.tree['el_LHTight'][iel] == 1:  self.fill_hist(sel, 'el_quality', 3)
			if self.tree['el_LHMedium'][iel] == 1: self.fill_hist(sel, 'el_quality', 2)
			if self.tree['el_LHLoose'][iel] == 1:  self.fill_hist(sel, 'el_quality', 1)
			else: self.fill_hist(sel, 'el_quality', 0)

	def _fill_all_dv_histos(self):
		sel = 'all'
		# self.fill_hist(sel, 'charge_ntrk', self.tree.dv('charge'), self.tree.dv('ntrk'))
		
		self.fill_hist(sel, 'DV_num_trks', self.tree.dv('ntrk'))
		self.fill_hist(sel, 'DV_x', self.tree.dv('x'))
		self.fill_hist(sel, 'DV_y', self.tree.dv('y'))
		self.fill_hist(sel, 'DV_z', self.tree.dv('z'))
		self.fill_hist(sel, 'DV_r', self.tree.dv('r'))
		self.fill_hist(sel, 'DV_distFromPV', self.tree.dv('distFromPV'))
		self.fill_hist(sel, 'DV_mass', self.tree.dv('mass'))
		self.fill_hist(sel, 'DV_pt', self.tree.dv('pt'))
		self.fill_hist(sel, 'DV_eta', self.tree.dv('eta'))
		self.fill_hist(sel, 'DV_phi', self.tree.dv('phi'))
		self.fill_hist(sel, 'DV_minOpAng', self.tree.dv('minOpAng'))
		self.fill_hist(sel, 'DV_maxOpAng', self.tree.dv('maxOpAng'))
		self.fill_hist(sel, 'DV_charge', self.tree.dv('charge'))
		self.fill_hist(sel, 'DV_chi2', self.tree.dv('chi2'))
		# self.fill_hist(sel, 'DV_chi2_assoc', self.tree.dv('chi2_assoc'))
		
		if sel == self.saveNtuples or self.saveNtuples == 'allcuts': 
			if self.MCEventType.isLNC: 
				self.micro_ntuples["LNC_"+sel].fill()
				self.micro_ntuples["LNC_plus_LNV_"+sel].fill()
			elif self.MCEventType.isLNV: 
				self.micro_ntuples["LNV_"+sel].fill()
				self.micro_ntuples["LNC_plus_LNV_"+sel].fill()
			else: self.micro_ntuples[sel].fill()
			

	def _fill_truth_histos(self, sel):
		truth_info = helpers.Truth()
		truth_info.getTruthParticles(self.tree)
		
		self.fill_hist(sel, 'event_type_MCweight', self.MCEventType.weight)  #if not weight_override else weight_override
		self.fill_hist(sel, 'M2_spin_corr_MCweight', self.MCEventType.M2_spin_corr)  #if not weight_override else weight_override
		self.fill_hist(sel, 'M2_nocorr_MCweight', self.MCEventType.M2_nocorr)  #if not weight_override else weight_override
		self.fill_hist(sel, 'W_pt', truth_info.W_vec.Pt())
		self.fill_hist(sel, 'W_eta', truth_info.W_vec.Eta())
		self.fill_hist(sel, 'W_phi', truth_info.W_vec.Phi())
		self.fill_hist(sel, 'W_mass', truth_info.W_vec.M())
		self.fill_hist(sel, 'HNLpt', truth_info.HNL_vec.Pt())
		self.fill_hist(sel, 'HNLeta', truth_info.HNL_vec.Eta())
		self.fill_hist(sel, 'HNLphi', truth_info.HNL_vec.Phi())
		self.fill_hist(sel, 'HNLm', truth_info.HNL_vec.M())

		self.fill_hist(sel, 'mHNLcalc', truth_info.mhnl)
		self.fill_hist(sel, 'DV_mass', truth_info.dvmass)

		self.fill_hist(sel, 'DV_r', truth_info.truth_dvr)
		
		self.fill_hist(sel, 'DV_x', truth_info.truth_dvx)
		self.fill_hist(sel, 'DV_y', truth_info.truth_dvy)
		self.fill_hist(sel, 'DV_z', truth_info.truth_dvz)

		self.fill_hist(sel, 'plep_pt', truth_info.plep_vec.Pt())
		self.fill_hist(sel, 'plep_eta', truth_info.plep_vec.Eta())
		self.fill_hist(sel, 'plep_phi', truth_info.plep_vec.Phi())
		self.fill_hist(sel, 'plep_mass', truth_info.plep_vec.M())
		self.fill_hist(sel, 'properLifetime', truth_info.properLifetime)
		
		if truth_info.W_charge == 1: 
			self.fill_hist(sel, 'Wplus_HNLpt', truth_info.HNL_vec.Pt())
			self.fill_hist(sel, 'Wplus_HNLeta', truth_info.HNL_vec.Eta())
			self.fill_hist(sel, 'Wplus_HNLphi', truth_info.HNL_vec.Phi())
			self.fill_hist(sel, 'Wplus_HNLE', truth_info.HNL_vec.E())


		if truth_info.W_charge == -1: 
			self.fill_hist(sel, 'Wminus_HNLpt', truth_info.HNL_vec.Pt())
			self.fill_hist(sel, 'Wminus_HNLeta', truth_info.HNL_vec.Eta())
			self.fill_hist(sel, 'Wminus_HNLphi', truth_info.HNL_vec.Phi())
			self.fill_hist(sel, 'Wminus_HNLE', truth_info.HNL_vec.E())
		if len(truth_info.trkVec) == 2: 
			DV_4vec= truth_info.trkVec[1]+ truth_info.trkVec[0]
			lep12 = truth_info.dLepVec[0] + truth_info.dLepVec[1] 
			lep23 = truth_info.dLepVec[1] + truth_info.dLepVec[2] 
			lep13 = truth_info.dLepVec[0] + truth_info.dLepVec[2] 
			all_leptons = truth_info.plep_vec + truth_info.trkVec[0] + truth_info.trkVec[1]
			self.fill_hist(sel, 'DV_mass', DV_4vec.M())
			self.fill_hist(sel, 'mvis', all_leptons.M())
			self.fill_hist(sel, 'm12', lep12.M())
			self.fill_hist(sel, 'm23', lep23.M())
			self.fill_hist(sel, 'm13', lep13.M())
			self.fill_hist(sel, 'm12_sq', lep12.M()**2)
			self.fill_hist(sel, 'm23_sq', lep23.M()**2)
			self.fill_hist(sel, 'm13_sq', lep13.M()**2)
			self.fill_hist(sel, 's12', self.MCEventType.s12) 
			self.fill_hist(sel, 's13', self.MCEventType.s13) 
			self.fill_hist(sel, 's14', self.MCEventType.s14) 
			self.fill_hist(sel, 's23', self.MCEventType.s23) 
			self.fill_hist(sel, 's24', self.MCEventType.s24) 
			self.fill_hist(sel, 's34', self.MCEventType.s34) 
			self.fill_hist(sel, 'lep1_trk_pt', self.MCEventType.p_2.Pt()) # topological ordered
			self.fill_hist(sel, 'lep1_trk_eta', self.MCEventType.p_2.Eta())
			self.fill_hist(sel, 'lep1_trk_phi', self.MCEventType.p_2.Phi())
			self.fill_hist(sel, 'lep2_trk_pt', self.MCEventType.p_3.Pt())
			self.fill_hist(sel, 'lep2_trk_eta', self.MCEventType.p_3.Eta())
			self.fill_hist(sel, 'lep2_trk_phi', self.MCEventType.p_3.Phi())
			self.fill_hist(sel, 'nu_trk_pt', self.MCEventType.p_4.Pt())
			self.fill_hist(sel, 'nu_trk_eta', self.MCEventType.p_4.Eta())
			self.fill_hist(sel, 'nu_trk_phi', self.MCEventType.p_4.Phi())
			if self.MCEventType.weight > 5: 
				self.fill_hist(sel, 'largew_plep_pt', truth_info.plep_vec.Pt())
				self.fill_hist(sel, 'largew_plep_eta', truth_info.plep_vec.Eta())
				self.fill_hist(sel, 'largew_plep_phi', truth_info.plep_vec.Phi())
				self.fill_hist(sel, 'largew_lep1_trk_pt', self.MCEventType.p_2.Pt())
				self.fill_hist(sel, 'largew_lep1_trk_eta', self.MCEventType.p_2.Eta())
				self.fill_hist(sel, 'largew_lep1_trk_phi', self.MCEventType.p_2.Phi())
				self.fill_hist(sel, 'largew_lep2_trk_pt', self.MCEventType.p_3.Pt())
				self.fill_hist(sel, 'largew_lep2_trk_eta', self.MCEventType.p_3.Eta())
				self.fill_hist(sel, 'largew_lep2_trk_phi', self.MCEventType.p_3.Phi())
				self.fill_hist(sel, 'largew_nu_trk_pt', self.MCEventType.p_4.Pt())
				self.fill_hist(sel, 'largew_nu_trk_eta', self.MCEventType.p_4.Eta())
				self.fill_hist(sel, 'largew_nu_trk_phi', self.MCEventType.p_4.Phi())

			disp_lep = [self.MCEventType.p_2,self.MCEventType.p_3,self.MCEventType.p_4]
			# pt order the displaced leptons
			disp_lep.sort(key=lambda x: x.Pt(), reverse=True)

			self.fill_hist(sel, 'dlep1_pt', disp_lep[0].Pt()) # pt ordered
			self.fill_hist(sel, 'dlep1_eta', disp_lep[0].Eta())
			self.fill_hist(sel, 'dlep1_phi', disp_lep[0].Phi())
			self.fill_hist(sel, 'dlep2_pt', disp_lep[1].Pt())
			self.fill_hist(sel, 'dlep2_eta', disp_lep[1].Eta())
			self.fill_hist(sel, 'dlep2_phi', disp_lep[1].Phi())
			self.fill_hist(sel, 'dlep3_pt', disp_lep[2].Pt())
			self.fill_hist(sel, 'dlep3_eta', disp_lep[2].Eta())
			self.fill_hist(sel, 'dlep3_phi', disp_lep[2].Phi())
	
			if (abs(truth_info.dTrk_d0[0]) < 2 and abs(truth_info.dTrk_d0[1]) < 2):
				self.fill_hist(sel, 'DV_d0_cut',1, fill_ntuple=False)
			else:
				self.fill_hist(sel, 'DV_d0_cut',0, fill_ntuple=False)

			n_el = len(truth_info.dEl)
			n_mu = len(truth_info.dMu)

		
			for iel in range(n_el):
				self.fill_hist(sel, 'DV_El_{}_pt'.format(iel), truth_info.dEl[iel].Pt(), fill_ntuple=False)
				self.fill_hist(sel, 'DV_El_{}_eta'.format(iel), truth_info.dEl[iel].Eta(), fill_ntuple=False)
				self.fill_hist(sel, 'DV_El_{}_phi'.format(iel), truth_info.dEl[iel].Phi(), fill_ntuple=False)
				self.fill_hist(sel, 'DV_El_{}_d0'.format(iel), truth_info.dEl_d0[iel], fill_ntuple=False)
				self.fill_hist(sel, 'DV_El_{}_charge'.format(iel), truth_info.dEl_charge[iel], fill_ntuple=False)
			
			for imu in range(n_mu):
				self.fill_hist(sel, 'DV_Mu_{}_pt'.format(imu), truth_info.dMu[imu].Pt(), fill_ntuple=False)
				self.fill_hist(sel, 'DV_Mu_{}_eta'.format(imu), truth_info.dMu[imu].Eta(), fill_ntuple=False)
				self.fill_hist(sel, 'DV_Mu_{}_phi'.format(imu), truth_info.dMu[imu].Phi(), fill_ntuple=False)
				self.fill_hist(sel, 'DV_Mu_{}_d0'.format(imu), truth_info.dMu_d0[imu], fill_ntuple=False)
				self.fill_hist(sel, 'DV_Ml_{}_charge'.format(imu), truth_info.dMu_charge[imu], fill_ntuple=False)

			for itrk in range(2):
				self.fill_hist(sel, 'DV_trk_pt', truth_info.trkVec[itrk].Pt(), fill_ntuple=False) # do the same here but also save charge for lepton truth matching later
				self.fill_hist(sel, 'DV_trk_eta', truth_info.trkVec[itrk].Eta(), fill_ntuple=False)
				self.fill_hist(sel, 'DV_trk_phi', truth_info.trkVec[itrk].Phi(), fill_ntuple=False)
				self.fill_hist(sel, 'DV_trk_d0',truth_info.dTrk_d0[itrk], fill_ntuple=False)


			# TODO: figure out a ntuple scheme that can store these variables as well
		if sel == self.saveNtuples or self.saveNtuples == 'allcuts': 
			if self.MCEventType.isLNC: 
				self.micro_ntuples["LNC_"+sel].fill()
				self.micro_ntuples["LNC_plus_LNV_"+sel].fill()
			elif self.MCEventType.isLNV: 
				self.micro_ntuples["LNV_"+sel].fill()
				self.micro_ntuples["LNC_plus_LNV_"+sel].fill()
			else: self.micro_ntuples[sel].fill()


	def _fill_selected_dv_histos(self, sel, do_lock=True):
		if self._locked < FILL_LOCKED and do_lock:
			# these are the histograms you only want to fill ONCE per DV

			# ____________________________________________________________
			# pileup
			try:
				self.fill_hist(sel, 'weight_pileup', self.tree['weight_pileup'])
				# pileup systematics
				self.fill_hist(sel, 'weight_pileup_1UP', self.tree['weight_pileup_up'])
				self.fill_hist(sel, 'weight_pileup_1DOWN', self.tree['weight_pileup_down'])
			except KeyError:
				pass

			# ____________________________________________________________
			# fill the DV weight for LNC or LNV only model assumption (Dirac neutrino)
			self.fill_hist(sel, 'DV_weight_LNC_only', self.weight_LNC_only)
			# fill the DV weight for LNC plus LNV signal model assumption (Majorana neutrino)
			# extra factor of 1/2 is from the U, m , ctau relationship changing due to there being twice
			# as many decay channels open if the HNL can decay both LNC and LNV
			self.fill_hist(sel, 'DV_weight_LNC_plus_LNV', self.weight_LNC_plus_LNV)
			self.fill_hist(sel, 'event_is_LNC', self.MCEventType.isLNC)
			self.fill_hist(sel, 'event_is_LNV', self.MCEventType.isLNV)

			# ____________________________________________________________
			tracks = helpers.Tracks(self.tree)
			tracks.getTracks()

			muons = helpers.Muons(self.tree)
			mu_vec = muons.lepVec
			mu_index = muons.lepIndex

			electrons = helpers.Electrons(self.tree)
			el_vec = electrons.lepVec
			el_index = electrons.lepIndex

			# fill histograms that require a prompt lepton to be identified
			if self.do_prompt_lepton_cut:
				plep_vec = self.plep_sel.plepVec
				plepd0 = self.plep_sel.plepd0
				plepz0 = self.plep_sel.plepz0
				plepcharge = self.plep_sel.plepcharge
				plep_isTight = self.plep_sel.plep_isTight
				plep_Index = self.plep_sel.plep_Index
				plep_is_trigger_matched = selections.LeptonTriggerMatching(self.tree, self.plep, plep_Index).is_trigger_matched

				self.fill_hist(sel, 'plep_pt', plep_vec.Pt())
				self.fill_hist(sel, 'plep_eta', plep_vec.Eta())
				self.fill_hist(sel, 'plep_phi', plep_vec.Phi())
				self.fill_hist(sel, 'plep_d0', plepd0)
				self.fill_hist(sel, 'plep_z0', plepz0)
				self.fill_hist(sel, 'plep_charge', plepcharge)
				self.fill_hist(sel, 'plep_isTight', plep_isTight)
				self.fill_hist(sel, 'plep_is_trigger_matched', plep_is_trigger_matched)

				if tracks.ntracks == 2:
					Mlll = selections.Mlll(dv_type=self.dv_type, plep=plep_vec, dMu=mu_vec, dEl=el_vec)
					Mhnl = selections.Mhnl(self.tree, self.dv_type, plep=plep_vec, dMu=mu_vec, dEl=el_vec)
					# Mhnl_old = selections.Mhnl_old(self.tree, plep=plep_vec, trks=trkVec)
					Mhnl_fixWmass = selections.Mhnl(self.tree, self.dv_type, plep=plep_vec, dMu=mu_vec, dEl=el_vec, fixWMass=True)

					self.fill_hist(sel, 'mvis', Mlll.mlll)
					self.fill_hist(sel, 'mtrans', Mlll.mtrans)
					self.fill_hist(sel, 'HNLm', Mhnl.mhnl)
					self.fill_hist(sel, 'HNLm_altbinning', Mhnl.mhnl)
					self.fill_hist(sel, 'alt_HNLm', Mhnl.alt_mhnl)
					self.fill_hist(sel, 'HNLm_fixWmass', Mhnl_fixWmass.mhnl)
					self.fill_hist(sel, 'HNLpt', Mhnl.hnlpt)
					self.fill_hist(sel, 'HNLeta', Mhnl.hnleta)
					self.fill_hist(sel, 'HNLphi', Mhnl.hnlphi)
					dR = tracks.lepVec[0].DeltaR(tracks.lepVec[1])
					if dR == 0.0:
						self.fill_hist(sel, 'DV_redmass', -1)
						self.fill_hist(sel, 'DV_redmassvis', -1)
						self.fill_hist(sel, 'DV_redmassHNL', -1)
					else:
						self.fill_hist(sel, 'DV_redmass', self.tree.dv('mass') / dR)
						self.fill_hist(sel, 'DV_redmassvis', Mlll.mlll / dR)
						self.fill_hist(sel, 'DV_redmassHNL', Mhnl.mhnl / dR)

				# get invariant mass between prompt lepton + same flavour displaced lepton for Z->ll veto
				if self.plep == 'muon' and len(mu_vec) > 0:
					zmass_veto_var = selections.Zmass_veto(self.tree, plep_vec=self.plep_sel.plepVec, plep=self.plep, plepcharge=self.plep_sel.plepcharge, dv_type=self.dv_type)
					self.fill_hist(sel, 'mll_dMu_plep_is_OS', zmass_veto_var.mll_dMu_plep_is_OS)
					self.fill_hist(sel, 'mll_dMu_plep_is_SS', zmass_veto_var.mll_dMu_plep_is_SS)
					self.fill_hist(sel, 'mll_dMu_plep', zmass_veto_var.mll_dMu_plep)

				if self.plep == 'electron' and len(el_vec) > 0:
					zmass_veto_var = selections.Zmass_veto(self.tree, plep_vec = self.plep_sel.plepVec, plep=self.plep, plepcharge = self.plep_sel.plepcharge, dv_type= self.dv_type)
					self.fill_hist(sel, 'mll_dEl_plep_is_OS', zmass_veto_var.mll_dEl_plep_is_OS)
					self.fill_hist(sel, 'mll_dEl_plep_is_SS', zmass_veto_var.mll_dEl_plep_is_SS)
					self.fill_hist(sel, 'mll_dEl_plep', zmass_veto_var.mll_dEl_plep )

			if self.do_prompt_track_cut:
				ptrk_vec = self.ptrk_sel.trkVec
				ptrkd0 = self.ptrk_sel.trkd0
				ptrkz0 = self.ptrk_sel.trkz0
				self.fill_hist(sel, 'ptrk_pt', ptrk_vec.Pt())
				self.fill_hist(sel, 'ptrk_eta', ptrk_vec.Eta())
				self.fill_hist(sel, 'ptrk_phi', ptrk_vec.Phi())
				self.fill_hist(sel, 'ptrk_d0', ptrkd0)
				self.fill_hist(sel, 'ptrk_z0', ptrkz0)

				if tracks.ntracks == 2:
					Mhnl = selections.Mhnl(self.tree, self.dv_type, plep=ptrk_vec, dMu=mu_vec, dEl=el_vec, use_tracks=True, trks=tracks.lepVec)
					self.fill_hist(sel, 'HNLm', Mhnl.mhnl)
					self.fill_hist(sel, 'HNLpt', Mhnl.hnlpt)
					self.fill_hist(sel, 'HNLeta', Mhnl.hnleta)
					self.fill_hist(sel, 'HNLphi', Mhnl.hnlphi)

			if tracks.ntracks == 2:
				deta = abs(tracks.eta[0] - tracks.eta[1])
				dphi = abs(tracks.lepVec[0].DeltaPhi(tracks.lepVec[1]))
				dpt = abs(tracks.pt[0] - tracks.pt[1])
				dR = tracks.lepVec[0].DeltaR(tracks.lepVec[1])
				self.fill_hist(sel, 'DV_trk_deta', deta)
				self.fill_hist(sel, 'DV_trk_dphi', dphi)
				self.fill_hist(sel, 'DV_trk_dpt', dpt)
				self.fill_hist(sel, 'DV_trk_dR', dR)

				cosmic_veto = selections.Cosmicveto(self.tree)
				self.fill_hist(sel, 'DV_cosmic_sep', cosmic_veto.separation)

				self.fill_hist(sel, 'DV_trk_max_chi2_toSV', max(self.tree.dv('trk_chi2_toSV')[0], self.tree.dv('trk_chi2_toSV')[1]))
				self.fill_hist(sel, 'DV_trk_min_chi2_toSV', min(self.tree.dv('trk_chi2_toSV')[0], self.tree.dv('trk_chi2_toSV')[1]))
				self.fill_hist(sel, 'DV_trk_max_d0_wrtSV', max(self.tree.dv('trk_d0_wrtSV')[0], self.tree.dv('trk_d0_wrtSV')[1]))
				self.fill_hist(sel, 'DV_trk_min_d0_wrtSV', min(self.tree.dv('trk_d0_wrtSV')[0], self.tree.dv('trk_d0_wrtSV')[1]))
				self.fill_hist(sel, 'DV_trk_max_errd0_wrtSV', max(self.tree.dv('trk_errd0_wrtSV')[0], self.tree.dv('trk_errd0_wrtSV')[1]))
				self.fill_hist(sel, 'DV_trk_min_errd0_wrtSV', min(self.tree.dv('trk_errd0_wrtSV')[0], self.tree.dv('trk_errd0_wrtSV')[1]))
				self.fill_hist(sel, 'DV_trk_max_z0_wrtSV', max(self.tree.dv('trk_z0_wrtSV')[0], self.tree.dv('trk_z0_wrtSV')[1]))
				self.fill_hist(sel, 'DV_trk_min_z0_wrtSV', min(self.tree.dv('trk_z0_wrtSV')[0], self.tree.dv('trk_z0_wrtSV')[1]))
				self.fill_hist(sel, 'DV_trk_max_errz0_wrtSV', max(self.tree.dv('trk_errz0_wrtSV')[0], self.tree.dv('trk_errz0_wrtSV')[1]))
				self.fill_hist(sel, 'DV_trk_min_errz0_wrtSV', min(self.tree.dv('trk_errz0_wrtSV')[0], self.tree.dv('trk_errz0_wrtSV')[1]))

				DV_mumu = selections.DVtype(self.tree, dv_type="mumu").passes()
				DV_ee = selections.DVtype(self.tree, dv_type="ee").passes()
				DV_emu = selections.DVtype(self.tree, dv_type="emu").passes()
				DV_1lep = (len(mu_vec) == 1 and len(el_vec) == 0) or (len(mu_vec) == 0 and len(el_vec) == 1)

				self.fill_hist(sel, 'DV_mumu', DV_mumu)
				self.fill_hist(sel, 'DV_ee', DV_ee)
				self.fill_hist(sel, 'DV_emu', DV_emu)
				self.fill_hist(sel, 'DV_1lep', DV_1lep)
				pass_el_mu_overlap = selections.electron_muon_overlap_check(self.tree).passes()
				self.fill_hist(sel, 'DV_pass_el_mu_overlap', pass_el_mu_overlap)

				pass_lep_pt_cut = selections.DV_lep_pt(self.tree, self.dv_type).pass_pt_cut
				self.fill_hist(sel, 'DV_pass_lep_pt', pass_lep_pt_cut)

				# calculate momentum parallel and perpendicular to the decay vector = DV-PV
				dv = ROOT.TVector3(self.tree.dv('x'), self.tree.dv('y'), self.tree.dv('z'))
				pv = ROOT.TVector3(self.tree['vertex_x'], self.tree['vertex_y'], self.tree['vertex_z'])
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
					self.fill_hist(sel, 'DV_trk_0_pt', self.tree.dv('trk_pt_wrtSV')[1])
					self.fill_hist(sel, 'DV_trk_0_eta', self.tree.dv('trk_eta_wrtSV')[1])
					self.fill_hist(sel, 'DV_trk_0_phi', self.tree.dv('trk_phi_wrtSV')[1])
					self.fill_hist(sel, 'DV_trk_0_d0', self.tree.dv('trk_d0')[1])
					self.fill_hist(sel, 'DV_trk_0_z0', self.tree.dv('trk_z0')[1])
					self.fill_hist(sel, 'DV_trk_0_charge', self.tree.dv('trk_charge')[1])
					self.fill_hist(sel, 'DV_trk_0_chi2', self.tree.dv('trk_chi2')[1])
					self.fill_hist(sel, 'DV_trk_0_isSelected', self.tree.dv('trk_isSelected')[1])
					self.fill_hist(sel, 'DV_trk_0_isAssociated', self.tree.dv('trk_isAssociated')[1])
					self.fill_hist(sel, 'DV_trk_0_mom_parall', mom_parall_1)
					self.fill_hist(sel, 'DV_trk_0_mom_perp', mom_perp_1)
					self.fill_hist(sel, 'DV_trk_0_mom_mag', pvec_1_mag)
					self.fill_hist(sel, 'DV_trk_0_mom_frac_parall', mom_frac_parall_1)

					self.fill_hist(sel, 'DV_trk_1_pt', self.tree.dv('trk_pt_wrtSV')[0])
					self.fill_hist(sel, 'DV_trk_1_eta', self.tree.dv('trk_eta_wrtSV')[0])
					self.fill_hist(sel, 'DV_trk_1_phi', self.tree.dv('trk_phi_wrtSV')[0])
					self.fill_hist(sel, 'DV_trk_1_d0', self.tree.dv('trk_d0')[0])
					self.fill_hist(sel, 'DV_trk_1_z0', self.tree.dv('trk_z0')[0])
					self.fill_hist(sel, 'DV_trk_1_charge', self.tree.dv('trk_charge')[0])
					self.fill_hist(sel, 'DV_trk_1_chi2', self.tree.dv('trk_chi2')[0])
					self.fill_hist(sel, 'DV_trk_1_isSelected', self.tree.dv('trk_isSelected')[0])
					self.fill_hist(sel, 'DV_trk_1_isAssociated', self.tree.dv('trk_isAssociated')[0])
					self.fill_hist(sel, 'DV_trk_1_mom_parall', mom_parall_0)
					self.fill_hist(sel, 'DV_trk_1_mom_perp', mom_perp_0)
					self.fill_hist(sel, 'DV_trk_1_mom_mag', pvec_0_mag)
					self.fill_hist(sel, 'DV_trk_1_mom_frac_parall', mom_frac_parall_0)
				else:
					self.fill_hist(sel, 'DV_trk_0_pt', self.tree.dv('trk_pt_wrtSV')[0])
					self.fill_hist(sel, 'DV_trk_0_eta', self.tree.dv('trk_eta_wrtSV')[0])
					self.fill_hist(sel, 'DV_trk_0_phi', self.tree.dv('trk_phi_wrtSV')[0])
					self.fill_hist(sel, 'DV_trk_0_d0', self.tree.dv('trk_d0')[0])
					self.fill_hist(sel, 'DV_trk_0_z0', self.tree.dv('trk_z0')[0])
					self.fill_hist(sel, 'DV_trk_0_charge', self.tree.dv('trk_charge')[0])
					self.fill_hist(sel, 'DV_trk_0_chi2', self.tree.dv('trk_chi2')[0])
					self.fill_hist(sel, 'DV_trk_0_isSelected', self.tree.dv('trk_isSelected')[0])
					self.fill_hist(sel, 'DV_trk_0_isAssociated', self.tree.dv('trk_isAssociated')[0])
					self.fill_hist(sel, 'DV_trk_0_mom_parall', mom_parall_0)
					self.fill_hist(sel, 'DV_trk_0_mom_perp', mom_perp_0)
					self.fill_hist(sel, 'DV_trk_0_mom_mag', pvec_0_mag)
					self.fill_hist(sel, 'DV_trk_0_mom_frac_parall', mom_frac_parall_0)

					self.fill_hist(sel, 'DV_trk_1_pt', self.tree.dv('trk_pt_wrtSV')[1])
					self.fill_hist(sel, 'DV_trk_1_eta', self.tree.dv('trk_eta_wrtSV')[1])
					self.fill_hist(sel, 'DV_trk_1_phi', self.tree.dv('trk_phi_wrtSV')[1])
					self.fill_hist(sel, 'DV_trk_1_d0', self.tree.dv('trk_d0')[1])
					self.fill_hist(sel, 'DV_trk_1_z0', self.tree.dv('trk_z0')[1])
					self.fill_hist(sel, 'DV_trk_1_charge', self.tree.dv('trk_charge')[1])
					self.fill_hist(sel, 'DV_trk_1_chi2', self.tree.dv('trk_chi2')[1])
					self.fill_hist(sel, 'DV_trk_1_isSelected', self.tree.dv('trk_isSelected')[1])
					self.fill_hist(sel, 'DV_trk_1_isAssociated', self.tree.dv('trk_isAssociated')[1])
					self.fill_hist(sel, 'DV_trk_1_mom_parall', mom_parall_1)
					self.fill_hist(sel, 'DV_trk_1_mom_perp', mom_perp_1)
					self.fill_hist(sel, 'DV_trk_1_mom_mag', pvec_1_mag)
					self.fill_hist(sel, 'DV_trk_1_mom_frac_parall', mom_frac_parall_1)

			# fill 3 different track calculations
			if self.dv_type == "mumu":
				lep_matched_DV_vec = muons.lepmatched_lepVec[0] + muons.lepmatched_lepVec[1]
				self.fill_hist(sel, 'DV_lep_0_trk_pt_wrtSV', muons.lepVec[0].Pt())
				self.fill_hist(sel, 'DV_lep_1_trk_pt_wrtSV', muons.lepVec[1].Pt())
				self.fill_hist(sel, 'DV_lep_0_std_trk_pt', muons.std_lepVec[0].Pt())
				self.fill_hist(sel, 'DV_lep_1_std_trk_pt', muons.std_lepVec[1].Pt())
				self.fill_hist(sel, 'DV_lep_0_lepmatched_trk_pt', muons.lepmatched_lepVec[0].Pt())
				self.fill_hist(sel, 'DV_lep_1_lepmatched_trk_pt', muons.lepmatched_lepVec[1].Pt())
				self.fill_hist(sel, 'DV_lep_0_lepmatched_trk_eta', muons.lepmatched_lepVec[0].Eta())
				self.fill_hist(sel, 'DV_lep_1_lepmatched_trk_eta', muons.lepmatched_lepVec[1].Eta())
				self.fill_hist(sel, 'DV_lep_0_lepmatched_trk_phi', muons.lepmatched_lepVec[0].Phi())
				self.fill_hist(sel, 'DV_lep_1_lepmatched_trk_phi', muons.lepmatched_lepVec[1].Phi())
				self.fill_hist(sel, 'DV_mass_lepmatched', lep_matched_DV_vec.M())

				# trigger matching for displaced leptons
				dmu_0_is_trig_matched = selections.LeptonTriggerMatching(self.tree, "muon", mu_index[0]).is_trigger_matched
				dmu_1_is_trig_matched = selections.LeptonTriggerMatching(self.tree, "muon", mu_index[1]).is_trigger_matched
				self.fill_hist(sel, 'DV_lep_0_is_trigger_matched', dmu_0_is_trig_matched)
				self.fill_hist(sel, 'DV_lep_1_is_trigger_matched', dmu_1_is_trig_matched)

				self.fill_hist(sel, 'DV_lep_0_isMuon', 1)
				self.fill_hist(sel, 'DV_lep_1_isMuon', 1)
				self.fill_hist(sel, 'DV_lep_0_isElectron', 0)
				self.fill_hist(sel, 'DV_lep_1_isElectron', 0)
				self.fill_hist(sel, 'DV_lep_0_muon_isLoose', self.tree.get('muon_isLoose')[muons.lepIndex[0]])
				self.fill_hist(sel, 'DV_lep_1_muon_isLoose', self.tree.get('muon_isLoose')[muons.lepIndex[1]])
				self.fill_hist(sel, 'DV_lep_0_muon_isMedium', self.tree.get('muon_isMedium')[muons.lepIndex[0]])
				self.fill_hist(sel, 'DV_lep_1_muon_isMedium', self.tree.get('muon_isMedium')[muons.lepIndex[1]])
				self.fill_hist(sel, 'DV_lep_0_muon_isTight', self.tree.get('muon_isTight')[muons.lepIndex[0]])
				self.fill_hist(sel, 'DV_lep_1_muon_isTight', self.tree.get('muon_isTight')[muons.lepIndex[1]])

			if self.dv_type == "emu":
				lep_matched_DV_vec = muons.lepmatched_lepVec[0] + electrons.lepmatched_lepVec[0]
				self.fill_hist(sel, 'DV_lep_0_trk_pt_wrtSV', muons.lepVec[0].Pt())
				self.fill_hist(sel, 'DV_lep_1_trk_pt_wrtSV', electrons.lepVec[0].Pt())
				self.fill_hist(sel, 'DV_lep_0_std_trk_pt', muons.std_lepVec[0].Pt())
				self.fill_hist(sel, 'DV_lep_1_std_trk_pt', electrons.std_lepVec[0].Pt())
				self.fill_hist(sel, 'DV_lep_0_lepmatched_trk_pt', muons.lepmatched_lepVec[0].Pt())
				self.fill_hist(sel, 'DV_lep_1_lepmatched_trk_pt', electrons.lepmatched_lepVec[0].Pt())
				self.fill_hist(sel, 'DV_lep_0_lepmatched_trk_eta', muons.lepmatched_lepVec[0].Eta())
				self.fill_hist(sel, 'DV_lep_1_lepmatched_trk_eta', electrons.lepmatched_lepVec[0].Eta())
				self.fill_hist(sel, 'DV_lep_0_lepmatched_trk_phi', muons.lepmatched_lepVec[0].Phi())
				self.fill_hist(sel, 'DV_lep_1_lepmatched_trk_phi', electrons.lepmatched_lepVec[0].Phi())
				self.fill_hist(sel, 'DV_mass_lepmatched', lep_matched_DV_vec.M())

				delta_el = electrons.lepVec[0].Pt() - electrons.lepmatched_lepVec[0].Pt()
				delta_mu = muons.lepVec[0].Pt() - muons.lepmatched_lepVec[0].Pt()
				self.fill_hist(sel, 'DV_trk_v_el_pt', delta_el / electrons.lepmatched_lepVec[0].Pt())
				self.fill_hist(sel, 'DV_trk_v_mu_pt', delta_mu / muons.lepmatched_lepVec[0].Pt())

				# trigger matching for displaced leptons
				dmu_is_trig_matched = selections.LeptonTriggerMatching(self.tree, "muon", mu_index[0]).is_trigger_matched
				del_is_trig_matched = selections.LeptonTriggerMatching(self.tree, "electron", el_index[0]).is_trigger_matched

				self.fill_hist(sel, 'DV_lep_0_is_trigger_matched', dmu_is_trig_matched)
				self.fill_hist(sel, 'DV_lep_1_is_trigger_matched', del_is_trig_matched)

				self.fill_hist(sel, 'DV_lep_0_isMuon', 1)
				self.fill_hist(sel, 'DV_lep_1_isMuon', 0)
				self.fill_hist(sel, 'DV_lep_0_isElectron', 0)
				self.fill_hist(sel, 'DV_lep_1_isElectron', 1)

				self.fill_hist(sel, 'DV_lep_0_muon_isLoose', self.tree.get('muon_isLoose')[muons.lepIndex[0]])
				self.fill_hist(sel, 'DV_lep_0_muon_isMedium', self.tree.get('muon_isMedium')[muons.lepIndex[0]])
				self.fill_hist(sel, 'DV_lep_0_muon_isTight', self.tree.get('muon_isTight')[muons.lepIndex[0]])

				self.fill_hist(sel, 'DV_lep_1_electron_LHTight', self.tree.get('el_LHTight')[electrons.lepIndex[0]])
				self.fill_hist(sel, 'DV_lep_1_electron_LHMedium', self.tree.get('el_LHMedium')[electrons.lepIndex[0]])
				self.fill_hist(sel, 'DV_lep_1_electron_LHLoose', self.tree.get('el_LHLoose')[electrons.lepIndex[0]])
				self.fill_hist(sel, 'DV_lep_1_electron_isLHVeryLoose', self.tree.get('el_isLHVeryLoose')[electrons.lepIndex[0]])
				self.fill_hist(sel, 'DV_lep_1_electron_VeryVeryLoose', self.tree.get('el_isLHVeryLoose_mod1')[electrons.lepIndex[0]])
				self.fill_hist(sel, 'DV_lep_1_electron_VeryVeryLooseSi', self.tree.get('el_isLHVeryLoose_modSi')[electrons.lepIndex[0]])

			if self.dv_type == "ee":
				lep_matched_DV_vec = electrons.lepmatched_lepVec[0] + electrons.lepmatched_lepVec[1]
				self.fill_hist(sel, 'DV_lep_0_trk_pt_wrtSV', electrons.lepVec[0].Pt())
				self.fill_hist(sel, 'DV_lep_1_trk_pt_wrtSV', electrons.lepVec[1].Pt())
				self.fill_hist(sel, 'DV_lep_0_std_trk_pt', electrons.std_lepVec[0].Pt())
				self.fill_hist(sel, 'DV_lep_1_std_trk_pt', electrons.std_lepVec[1].Pt())
				self.fill_hist(sel, 'DV_lep_0_lepmatched_trk_pt', electrons.lepmatched_lepVec[0].Pt())
				self.fill_hist(sel, 'DV_lep_1_lepmatched_trk_pt', electrons.lepmatched_lepVec[1].Pt())
				self.fill_hist(sel, 'DV_lep_0_lepmatched_trk_eta', electrons.lepmatched_lepVec[0].Eta())
				self.fill_hist(sel, 'DV_lep_1_lepmatched_trk_eta', electrons.lepmatched_lepVec[1].Eta())
				self.fill_hist(sel, 'DV_lep_0_lepmatched_trk_phi', electrons.lepmatched_lepVec[0].Phi())
				self.fill_hist(sel, 'DV_lep_1_lepmatched_trk_phi', electrons.lepmatched_lepVec[1].Phi())
				self.fill_hist(sel, 'DV_mass_lepmatched', lep_matched_DV_vec.M())

				# trigger matching for displaced leptons
				del_0_is_trig_matched = selections.LeptonTriggerMatching(self.tree, "electron", el_index[0]).is_trigger_matched
				del_1_is_trig_matched = selections.LeptonTriggerMatching(self.tree, "electron", el_index[1]).is_trigger_matched
				self.fill_hist(sel, 'DV_lep_0_is_trigger_matched', del_0_is_trig_matched)
				self.fill_hist(sel, 'DV_lep_1_is_trigger_matched', del_1_is_trig_matched)

				self.fill_hist(sel, 'DV_lep_0_isElectron', 1)
				self.fill_hist(sel, 'DV_lep_1_isElectron', 1)
				self.fill_hist(sel, 'DV_lep_0_isMuon', 0)
				self.fill_hist(sel, 'DV_lep_1_isMuon', 0)
				self.fill_hist(sel, 'DV_lep_0_electron_LHTight', self.tree.get('el_LHTight')[electrons.lepIndex[0]])
				self.fill_hist(sel, 'DV_lep_0_electron_LHMedium', self.tree.get('el_LHMedium')[electrons.lepIndex[0]])
				self.fill_hist(sel, 'DV_lep_0_electron_LHLoose', self.tree.get('el_LHLoose')[electrons.lepIndex[0]])
				self.fill_hist(sel, 'DV_lep_0_electron_isLHVeryLoose', self.tree.get('el_isLHVeryLoose')[electrons.lepIndex[0]])
				self.fill_hist(sel, 'DV_lep_0_electron_VeryVeryLoose', self.tree.get('el_isLHVeryLoose_mod1')[electrons.lepIndex[0]])
				self.fill_hist(sel, 'DV_lep_0_electron_VeryVeryLooseSi', self.tree.get('el_isLHVeryLoose_modSi')[electrons.lepIndex[0]])
				self.fill_hist(sel, 'DV_lep_1_electron_LHTight', self.tree.get('el_LHTight')[electrons.lepIndex[1]])
				self.fill_hist(sel, 'DV_lep_1_electron_LHMedium', self.tree.get('el_LHMedium')[electrons.lepIndex[1]])
				self.fill_hist(sel, 'DV_lep_1_electron_LHLoose', self.tree.get('el_LHLoose')[electrons.lepIndex[1]])
				self.fill_hist(sel, 'DV_lep_1_electron_isLHVeryLoose', self.tree.get('el_isLHVeryLoose')[electrons.lepIndex[1]])
				self.fill_hist(sel, 'DV_lep_1_electron_VeryVeryLoose', self.tree.get('el_isLHVeryLoose_mod1')[electrons.lepIndex[1]])
				self.fill_hist(sel, 'DV_lep_1_electron_VeryVeryLooseSi', self.tree.get('el_isLHVeryLoose_modSi')[electrons.lepIndex[1]])

			# fill standard track variables for electrons
			for lep in range(len(el_vec)):
				self.fill_hist(sel, 'DV_el_trk_pt', el_vec[lep].Pt())
				self.fill_hist(sel, 'DV_el_trk_eta', el_vec[lep].Eta())
				self.fill_hist(sel, 'DV_el_trk_phi', el_vec[lep].Phi())

			# fill standard track variables for muons
			for lep in range(len(mu_vec)):
				self.fill_hist(sel, 'DV_mu_trk_pt', mu_vec[lep].Pt())
				self.fill_hist(sel, 'DV_mu_trk_eta', mu_vec[lep].Eta())
				self.fill_hist(sel, 'DV_mu_trk_phi', mu_vec[lep].Phi())

			# fill standard track variable histograms
			for i in range(tracks.ntracks):
				self.fill_hist(sel, 'DV_trk_pt', self.tree.dv('trk_pt_wrtSV')[i])
				self.fill_hist(sel, 'DV_trk_eta', self.tree.dv('trk_eta_wrtSV')[i])
				self.fill_hist(sel, 'DV_trk_phi', self.tree.dv('trk_phi_wrtSV')[i])
				self.fill_hist(sel, 'DV_trk_d0', self.tree.dv('trk_d0')[i])
				self.fill_hist(sel, 'DV_trk_z0', self.tree.dv('trk_z0')[i])
				self.fill_hist(sel, 'DV_trk_absz0', abs(self.tree.dv('trk_z0')[i]))
				self.fill_hist(sel, 'DV_trk_charge', self.tree.dv('trk_charge')[i])
				self.fill_hist(sel, 'DV_trk_chi2', self.tree.dv('trk_chi2')[i])
				self.fill_hist(sel, 'DV_trk_isLRT', self.tree.dv('trk_isLRT')[i])
				self.fill_hist(sel, 'DV_trk_isSelected', self.tree.dv('trk_isSelected')[i])
				self.fill_hist(sel, 'DV_trk_isAssociated', self.tree.dv('trk_isAssociated')[i])
				self.fill_hist(sel, 'DV_trk_nPixelHits', self.tree.dv('trk_nPixelHits')[i])
				self.fill_hist(sel, 'DV_trk_nSCTHits', self.tree.dv('trk_nSCTHits')[i])
				# self.fill_hist(sel, 'DV_trk_nSCTHoles', self.tree.dv('trk_nSCTHoles')[i])
				self.fill_hist(sel, 'DV_trk_nSiHits', self.tree.dv('trk_nSCTHits')[i] + self.tree.dv('trk_nPixelHits')[i])
				# self.fill_hist(sel, 'DV_trk_dTheta', self.tree.dv('trk_dTheta')[i])
				self.fill_hist(sel, 'DV_trk_chi2_toSV'.format(i), self.tree.dv('trk_chi2_toSV')[i])
				self.fill_hist(sel, 'DV_trk_d0_wrtSV'.format(i), self.tree.dv('trk_d0_wrtSV')[i])
				self.fill_hist(sel, 'DV_trk_errd0_wrtSV'.format(i), self.tree.dv('trk_errd0_wrtSV')[i])
				self.fill_hist(sel, 'DV_trk_z0_wrtSV'.format(i), self.tree.dv('trk_z0_wrtSV')[i])
				self.fill_hist(sel, 'DV_trk_errz0_wrtSV'.format(i), self.tree.dv('trk_errz0_wrtSV')[i])

			# fill standard dv histograms
			self.fill_hist(sel, 'DV_num_trks', self.tree.dv('ntrk'))
			self.fill_hist(sel, 'DV_x', self.tree.dv('x'))
			self.fill_hist(sel, 'DV_y', self.tree.dv('y'))
			self.fill_hist(sel, 'DV_z', self.tree.dv('z'))
			self.fill_hist(sel, 'DV_r', self.tree.dv('r'))
			self.fill_hist(sel, 'PV_x', self.tree['vertex_x'])
			self.fill_hist(sel, 'PV_y', self.tree['vertex_y'])
			self.fill_hist(sel, 'PV_z', self.tree['vertex_z'])
			self.fill_hist(sel, 'DV_distFromPV', self.tree.dv('distFromPV'))
			self.fill_hist(sel, 'DV_mass', self.tree.dv('mass'))
			self.fill_hist(sel, 'DV_pt', self.tree.dv('pt'))
			self.fill_hist(sel, 'DV_eta', self.tree.dv('eta'))
			self.fill_hist(sel, 'DV_phi', self.tree.dv('phi'))
			self.fill_hist(sel, 'DV_minOpAng', self.tree.dv('minOpAng'))
			self.fill_hist(sel, 'DV_maxOpAng', self.tree.dv('maxOpAng'))
			self.fill_hist(sel, 'DV_charge', self.tree.dv('charge'))
			self.fill_hist(sel, 'DV_chi2', self.tree.dv('chi2'))
			# self.fill_hist(sel, 'DV_chi2_assoc', self.tree.dv('chi2_assoc'))
			self.fill_hist(sel, 'DV_max_dR', self.tree.dv('maxDR'))
			self.fill_hist(sel, 'DV_max_dR_wrtSV', self.tree.dv('maxDR_wrtSV'))
			self.fill_hist(sel, 'DV_maxd0', self.tree.dv('maxd0'))
			self.fill_hist(sel, 'DV_mind0', self.tree.dv('mind0'))
			self.fill_hist(sel, 'DV_ntrk', self.tree.dv('ntrk'))
			self.fill_hist(sel, 'DV_ntrk_lrt', self.tree.dv('ntrk_lrt'))
			self.fill_hist(sel, 'DV_ntrk_sel', self.tree.dv('ntrk_sel'))
			self.fill_hist(sel, 'DV_ntrk_assoc', self.tree.dv('ntrk_assoc'))
			self.fill_hist(sel, 'DV_pass_mat_veto', self.tree.dv('pass_mat'))

			
			if not self.tree.is_data:
				maxlinkTruth_score = self.tree.dv('maxlinkTruth_score')
				maxlinkTruth_parent_pdgId = abs(self.tree.dv('maxlinkTruth_parent_pdgId'))
				is_truth_matched =  self.tree.dv('maxlinkTruth_score') > 0.75 and abs(self.tree.dv('maxlinkTruth_parent_pdgId')) == 50
				#is truth matched: 
				self.fill_hist(sel, 'DV_truth_matched', is_truth_matched)

				# add proper lifetime
				truth_info = helpers.Truth()
				truth_info.getTruthParticles(self.tree)
				self.fill_hist(sel, 'properLifetime', truth_info.properLifetime)

				if self.dv_type == "mumu":
					# Get the truth index (truth matching by charge)
					if not self.tree.not_hnl_mc:
						if self.tree.dv('trk_charge')[0] == truth_info.dMu_charge[0]: trk_0_truth_index = 0
						elif self.tree.dv('trk_charge')[0] == truth_info.dMu_charge[1]: trk_0_truth_index = 1
						else:	raise Exception("Can't truth match lepton by charge. Something is strange.")

						if self.tree.dv('trk_charge')[1] == truth_info.dMu_charge[0]: trk_1_truth_index = 0
						elif self.tree.dv('trk_charge')[1] == truth_info.dMu_charge[1]: trk_1_truth_index = 1
						else:	raise Exception("Can't truth match lepton by charge. Something is strange.")

						self.fill_hist(sel, 'DV_trk_0_d0_truth', truth_info.dMu_d0[trk_0_truth_index])
						self.fill_hist(sel, 'DV_trk_1_d0_truth', truth_info.dMu_d0[trk_1_truth_index])

						for i in range(len(muons.lepVec)):
							delta = muons.lepVec[i].Pt() - muons.lepmatched_lepVec[i].Pt()
							self.fill_hist(sel, 'DV_trk_v_mu_pt', delta/muons.lepmatched_lepVec[i].Pt() )
				if self.dv_type == "emu":
					# truth d0 # it looks like these are already matched so track0 is always the muon. Could that be?
					if not self.tree.not_hnl_mc:
						self.fill_hist(sel, 'DV_trk_0_d0_truth', truth_info.dMu_d0[0])
						self.fill_hist(sel, 'DV_trk_1_d0_truth', truth_info.dEl_d0[0])

				if self.dv_type == "ee":
					if not self.tree.not_hnl_mc:
						# Get the truth index (truth matching by charge)
						if self.tree.dv('trk_charge')[0] == truth_info.dEl_charge[0]: trk_0_truth_index = 0
						elif self.tree.dv('trk_charge')[0] == truth_info.dEl_charge[1]: trk_0_truth_index = 1
						else:	raise Exception("Can't truth match lepton by charge. Something is strange.")

						if self.tree.dv('trk_charge')[1] == truth_info.dEl_charge[0]: trk_1_truth_index = 0
						elif self.tree.dv('trk_charge')[1] == truth_info.dEl_charge[1]: trk_1_truth_index = 1
						else:	raise Exception("Can't truth match lepton by charge. Something is strange.")

						self.fill_hist(sel, 'DV_trk_0_d0_truth', truth_info.dEl_d0[trk_0_truth_index])
						self.fill_hist(sel, 'DV_trk_1_d0_truth', truth_info.dEl_d0[trk_1_truth_index])

						for i in range(len(muons.lepVec)):
							delta = electrons.lepVec[i].Pt() - electrons.lepmatched_lepVec[i].Pt()
							self.fill_hist(sel, 'DV_trk_v_el_pt', delta/electrons.lepmatched_lepVec[i].Pt() )



			self.fill_hist(sel, 'DV_alpha', selections.Alpha(self.tree).alpha)
			
			trk_quality = selections.Trackqual(self.tree)

			self.fill_hist(sel, 'DV_2tight', trk_quality.DV_2tight)
			self.fill_hist(sel, 'DV_2medium', trk_quality.DV_2medium)
			self.fill_hist(sel, 'DV_2loose', trk_quality.DV_2loose)
			self.fill_hist(sel, 'DV_1tight', trk_quality.DV_1tight)
			self.fill_hist(sel, 'DV_1medium', trk_quality.DV_1medium)
			self.fill_hist(sel, 'DV_1loose', trk_quality.DV_1loose)
			self.fill_hist(sel, 'DV_tight_loose', trk_quality.DV_tight_loose)
			self.fill_hist(sel, 'DV_tight_medium', trk_quality.DV_tight_medium)
			self.fill_hist(sel, 'DV_medium_loose', trk_quality.DV_medium_loose)
			self.fill_hist(sel, 'DV_tight_veryloose', trk_quality.DV_tight_veryloose)
			self.fill_hist(sel, 'DV_medium_veryloose', trk_quality.DV_medium_veryloose)
			self.fill_hist(sel, 'DV_loose_veryloose', trk_quality.DV_loose_veryloose)
			self.fill_hist(sel, 'DV_tight_veryveryloose', trk_quality.DV_tight_veryveryloose)
			self.fill_hist(sel, 'DV_medium_veryveryloose', trk_quality.DV_medium_veryveryloose)
			self.fill_hist(sel, 'DV_loose_veryveryloose', trk_quality.DV_loose_veryveryloose)
			self.fill_hist(sel, 'DV_2veryveryloose', trk_quality.DV_2veryveryloose)
			self.fill_hist(sel, 'DV_1veryveryloose', trk_quality.DV_1veryveryloose)

			# ____________________________________________________________
			# Trigger matching requirement
			if not self.do_CR:
				self.trigger_matched_medium = selections.RequireMediumTriggerMatching(
					self.tree,
					prompt_lepton_index=self.plep_sel.plep_Index,
					prompt_lepton_type=self.plep,
					muons=muons,
					electrons=electrons,
					dv_type=self.dv_type)

				self.fill_hist(sel, 'n_trigger_matched_medium', self.trigger_matched_medium.n_trigger_matched_medium)
				# This will also be used for applying the lepton trigger matching uncertainties

			# ============================================================
			# Systematics

			# ____________________________________________________________
			# Calculate weight for tracking/vertexing uncertainties
			vertexing_systematic = systematics.get_vertexing_systematic(self.tree.dv('r'), self.tree.dv('pt'))
			self.fill_hist(sel, 'vertexing_1DOWN', vertexing_systematic, weight=1)

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
			self.fill_hist(sel, 'd0_extrapolation_1DOWN', d0_extrapolation_systematic, weight=1)

			# ____________________________________________________________
			# Calculate weight for lepton trigger scale factor uncertainty
			if not self.do_CR:
				if self.trigger_matched_medium.prompt_trigger_matched_medium:
					# trigger uncertainty from prompt lepton (TODO: needs to be implemented)
					pass
				elif self.trigger_matched_medium.displaced_0_trigger_matched_medium:
					# trigger uncertainty from displaced lepton 0 (TODO: needs to be implemented)
					pass
				elif self.trigger_matched_medium.displaced_1_trigger_matched_medium:
					# trigger uncertainty from displaced lepton 1 (TODO: needs to be implemented)
					pass

			# ============================================================
			# fill TTree with ntuple information. Already set by fill_hist
			if sel == self.saveNtuples or self.saveNtuples == 'allcuts':
				if self.MCEventType.isLNC: 
					self.micro_ntuples["LNC_"+sel].fill()
					self.micro_ntuples["LNC_plus_LNV_"+sel].fill()
				elif self.MCEventType.isLNV: 
					self.micro_ntuples["LNV_"+sel].fill()
					self.micro_ntuples["LNC_plus_LNV_"+sel].fill()
				else: self.micro_ntuples[sel].fill()

			if sel == "sel":
				self._locked = FILL_LOCKED  # this only becomes unlocked after the event loop finishes in makeHistograms so you can only fill one DV from each event.


class run2Analysis(Analysis):
	def __init__(self, name, tree, vtx_container, selections, outputFile, saveNtuples, weight_override=None):

		Analysis.__init__(self, name, tree, vtx_container, selections, outputFile, saveNtuples, weight_override)
		self.logger.info('Running  Full Run 2 Analysis cuts')

		# Define cutflow histogram "by hand"
		self.cutflow_dir = self.ch + '/CutFlow/'
		if not self.dv_type == "ee":
			self.observables.histogram_dict[self.cutflow_dir+ 'CutFlow'] = ROOT.TH1D('CutFlow', 'CutFlow', 21, -0.5, 20.5)
		else:
			self.observables.histogram_dict[self.cutflow_dir+ 'CutFlow'] = ROOT.TH1D('CutFlow', 'CutFlow', 22, -0.5, 21.5)
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
			self.CutFlow.GetXaxis().SetBinLabel(6, "no plep overlap with DV")
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
		if not self.tree.is_data:

			if not self.dv_type == "ee": truth_match_bin = 21
			else: truth_match_bin = 22
			self.CutFlow.GetXaxis().SetBinLabel(truth_match_bin, "truth matched")
		self.CutFlow_LNC_plus_LNV = self.CutFlow.Clone()
		self.CutFlow_LNC_plus_LNV.SetName("CutFlow_LNC_plus_LNV"+"_"+self.ch)
		self.observables.histogram_dict[self.cutflow_dir+'CutFlow_LNC_plus_LNV'] = self.CutFlow_LNC_plus_LNV
		# Store LNC and LNV cutflows in the observables collection
		if not self.tree.is_data and not self.tree.not_hnl_mc:
			self.CutFlow_LNV = self.CutFlow.Clone()
			self.CutFlow_LNC = self.CutFlow.Clone()
			self.CutFlow_LNV.SetName("CutFlow_LNV"+"_"+self.ch)
			self.CutFlow_LNC.SetName("CutFlow_LNC"+"_"+self.ch)
			self.CutFlow_LNV_weighted = self.CutFlow.Clone()
			self.CutFlow_LNC_weighted = self.CutFlow.Clone()
			self.CutFlow_LNV_weighted.SetName("CutFlow_weighted_LNV"+"_"+self.ch)
			self.CutFlow_LNC_weighted.SetName("CutFlow_weighted_LNC"+"_"+self.ch)

			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_LNV'] = self.CutFlow_LNV
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_LNC'] = self.CutFlow_LNC
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_weighted_LNV'] = self.CutFlow_LNV_weighted
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_weighted_LNC'] = self.CutFlow_LNC_weighted

	def DVSelection(self):
		######################################################################################################
		# DV Selection is any cuts that are done per DV
		# Current cuts include: fiducial vol, ntrack, OS, DVtype, track quality, cosmic veto, mlll, mDV
		######################################################################################################
				
		# Fill all the histograms with ALL DVs (this could be more that 1 per event). Useful for vertexing efficiency studies.
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
						self._fill_selected_dv_histos("DVtype")
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
						self._fill_selected_dv_histos("cosmic")
			else:
				return

		if self.do_lep_pt_cut:
			if self._dv_lep_pt_cut():
				if not self.passed_dlep_pt_cut:
					self._fill_cutflow(12)
					self.passed_lep_pt_cut = True
					if not self.do_inverted_mlll_cut and not self.do_inverted_mhnl_cut:
						self._fill_selected_dv_histos("DVlep_pt")
			else:
				return

		if self.do_mat_veto_cut:
			if self._mat_veto_cut():
				if not self.passed_mat_veto_cut:
					self._fill_cutflow(13)
					self.passed_mat_veto_cut = True
					if not self.do_inverted_mlll_cut and not self.do_inverted_mhnl_cut:
						self._fill_selected_dv_histos("matveto")
			else:
				return

		if self.do_track_quality_cut:
			if self._track_quality_cut():
				if not self.passed_track_quality_cut:
					if not self.dv_type == "ee": self._fill_cutflow(13)
					else: self._fill_cutflow(13+1)
					self.passed_track_quality_cut = True
					if not self.do_inverted_mlll_cut and not self.do_inverted_mhnl_cut:
						self._fill_selected_dv_histos("trkqual")
			else:
				return

		if self.do_trigger_matching_cut:
			if self._trigger_matched_medium_lepton_cut():
				if not self.passed_trig_matched_cut:
					if not self.dv_type == "ee": self._fill_cutflow(14)
					else: self._fill_cutflow(14+1)
					self.passed_trig_matched_cut = True
					if not self.do_inverted_mlll_cut and not self.do_inverted_mhnl_cut:
						self._fill_selected_dv_histos("trig_match")
			else:
				return

		if self.do_trilepton_mass_cut:
			if self._trilepton_mass_cut():
				if not self.passed_trilepton_mass_cut:
					if not self.dv_type == "ee": self._fill_cutflow(15)
					else: self._fill_cutflow(15+1)
					self.passed_trilepton_mass_cut = True
					self._fill_selected_dv_histos("mvis")
			else:
				return
		
		if self.do_inverted_mlll_cut:
			if self._invert_trilepton_mass_cut():
				if not self.passed_trilepton_mass_cut:
					if not self.dv_type == "ee": self._fill_cutflow(15)
					else: self._fill_cutflow(15+1)
					self.passed_trilepton_mass_cut = True
					self._fill_selected_dv_histos("mvis")
			else:
				return
		
		if self.do_inverted_mhnl_cut:
			if self._invert_HNL_mass_cut():
				if not self.passed_HNL_mass_cut:
					if not self.dv_type == "ee": self._fill_cutflow(15)
					else: self._fill_cutflow(15+1)
					self.passed_HNL_mass_cut = True
					self._fill_selected_dv_histos("mhnl")
			else:
				return

		if self.do_dv_mass_cut:
			# do b-hadron veto instead of just straight DV mass cut
			if self._Bhadron_veto():
				if not self.passed_dv_mass_cut:
					if not self.dv_type == "ee": self._fill_cutflow(16)
					else: self._fill_cutflow(16+1)
					self.passed_dv_mass_cut = True
					self._fill_selected_dv_histos("mDV")
			else:
				return

		if self.do_zmass_veto:
			# do z mass veto for all channels
			if self._Zmass_veto():
				if not self.passed_zmass_veto:
					if not self.dv_type == "ee": self._fill_cutflow(17)
					else: self._fill_cutflow(17+1)
					self.passed_zmass_veto = True
					self._fill_selected_dv_histos("Zmass_veto")
			else:
				return

		if self.do_HNL_mass_cut :
			if self._HNL_mass_cut():
				if not self.passed_HNL_mass_cut:
					if not self.dv_type == "ee": self._fill_cutflow(18)
					else: self._fill_cutflow(18+1)
					self.passed_HNL_mass_cut = True
					self._fill_selected_dv_histos("mHNL")
			else:
				return

		if self.do_alpha_cut:
			if self._alpha_cut():
				if not self.passed_alpha_cut:
					if not self.dv_type == "ee": self._fill_cutflow(19)
					else: self._fill_cutflow(19+1)
					self.passed_alpha_cut = True
					self._fill_selected_dv_histos("alpha")
			else:
				return

		# Fill histos of truth-matched DVs
		if not self.tree.is_data and not self.tree.not_hnl_mc:
			if self._truth_match():
				if not self.dv_type == "ee": self._fill_cutflow(20)
				else: self._fill_cutflow(20+1)
				self._fill_selected_dv_histos("match")



class BEAnalysis(Analysis):
	def __init__(self, name, tree, vtx_container, selections, outputFile, saveNtuples, weight_override=None):

		Analysis.__init__(self, name, tree, vtx_container, selections, outputFile, saveNtuples, weight_override)
		self.logger.info('Running Background Estimate Analysis Cuts')

		# Define cutflow histogram "by hand"
		self.cutflow_dir = self.ch + '/CutFlow/'
		self.observables.histogram_dict[self.cutflow_dir + 'CutFlow'] = ROOT.TH1D('CutFlow', 'CutFlow', 17, -0.5, 16.5)
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
			self.CutFlow.GetXaxis().SetBinLabel(5, "{} prompt {}".format(self.plep_quality, self.plep))
			self.CutFlow.GetXaxis().SetBinLabel(6, "no plep overlap with DV")
		if self.do_invert_prompt_lepton_cut:
			self.CutFlow.GetXaxis().SetBinLabel(6, "invert prompt lepton")
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
		if self.do_same_event_cut:
			self.CutFlow.GetXaxis().SetBinLabel(11, "Same-Event")
		if self.do_different_event_cut:
			self.CutFlow.GetXaxis().SetBinLabel(11, "Different-Events")
		if self.do_dv_type_cut:
			self.CutFlow.GetXaxis().SetBinLabel(12, "%s DV" % self.dv_type)
		if self.do_mat_veto_cut:
			self.CutFlow.GetXaxis().SetBinLabel(13, "mat. veto")
		if self.do_cosmic_veto_cut:
			self.CutFlow.GetXaxis().SetBinLabel(14, "cosmic veto")
		if self.do_dv_mass_cut:
			self.CutFlow.GetXaxis().SetBinLabel(15, "m_{DV}")
		if self.do_track_quality_cut:
			self.CutFlow.GetXaxis().SetBinLabel(16, "{}-lepton DV".format(self.track_quality))
		if self.do_trilepton_mass_cut:
			self.CutFlow.GetXaxis().SetBinLabel(17, "m_{lll}")

	def DVSelection(self):
		######################################################################################################
		# DV Selection is any cuts that are done per DV
		# Current cuts include: fiducial vol, ntrack, OS, DVtype, track quality, cosmic veto, mlll, mDV
		######################################################################################################
		# Fill all the histograms with ALL DVs (this could be more that 1 per event). Useful for vertexing efficiency studies.
		self._fill_all_dv_histos()

		# We only want to save one DV per event: Choose first DV that passes "DV type" selection as "the DV for the event"
		# If we have already selected our DV for this event, don't bother looking at the rest of the DVs!
		if self.selected_dv_index != -1: return

		# There is an extra bit of logic here since we look at several DVs
		# but only want to fill the cutflow once per event
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

		if self.do_same_event_cut or self.do_different_event_cut:
			if self._be_event_type_cut():
				if not self.passed_be_event_type_cut:
					self._fill_cutflow(10)
					self.passed_be_event_type_cut = True
			else:
				return

		if self.do_dv_type_cut:
			if self._dv_type_cut():
				if not self.passed_dv_type_cut:
					self._fill_cutflow(11)
					self.passed_dv_type_cut = True
					# save the first DV that passes the DVtype selection
					self._fill_selected_dv_histos("DVtype")
					# Select this DV as the DV for the event!
					self.selected_dv_index = self.tree.idv
			else:
				return

		if self.do_cosmic_veto_cut:
			if self._cosmic_veto_cut():
				if not self.passed_cosmic_veto_cut:
					self._fill_cutflow(13)
					self.passed_cosmic_veto_cut = True
					self._fill_selected_dv_histos("cosmic")
			else:
				return

		if self.do_dv_mass_cut:
			if self._dv_mass_cut():
				if not self.passed_dv_mass_cut:
					self._fill_cutflow(14)
					self.passed_dv_mass_cut = True
					self._fill_selected_dv_histos("mDV")
			else:
				return

		if self.do_track_quality_cut:
			if self._track_quality_cut():
				if not self.passed_track_quality_cut:
					self._fill_cutflow(15)
					self.passed_track_quality_cut = True
					self._fill_selected_dv_histos("trkqual")
			else:
				return

		if self.do_trilepton_mass_cut:
			if self._trilepton_mass_cut():
				if not self.passed_trilepton_mass_cut:
					self._fill_cutflow(16)
					self.passed_trilepton_mass_cut = True
			else:
				return
		# self._fill_selected_dv_histos("mlll")

		# Fill all the histograms with only selected DVs. (ie. the ones that pass the full selection)
		self._fill_selected_dv_histos("sel")


class KShort(Analysis):
	def __init__(self, name, tree, vtx_container, selections, outputFile, saveNtuples, weight_override=None):
		Analysis.__init__(self, name, tree, vtx_container, selections, outputFile, saveNtuples, weight_override)
		self.logger.info('Running KShort Analysis cuts')

		self.add('CutFlow', 17, -0.5, 16.5)
		# Bin labels are 1 greater than histogram bins
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(1, "all")
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(2, "PV")
		if self.do_invert_prompt_lepton_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(3, "invert prompt lepton")
		if self.do_alpha_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(4, "alpha")
		if self.do_mass_window_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(5, "K0 mass")
		# add other histograms

	def _fill_leptons(self):
		sel = 'all'
		prompt = selections.InvertedPromptLepton(self.tree)
		self.fill_hist(sel, 'prompt_muon', prompt.n_prompt_muons)
		self.fill_hist(sel, 'prompt_electron', prompt.n_prompt_electrons)
		self.fill_hist(sel, 'prompt_lepton', prompt.n_prompt_leptons)
		for imu in range(len(self.tree['muon_type'])):
			self.fill_hist(sel, 'muon_type', self.tree['muon_type'][imu])
			self.fill_hist(sel, 'muon_pt', self.tree['muon_pt'][imu])
			self.fill_hist(sel, 'muon_eta', self.tree['muon_eta'][imu])
			self.fill_hist(sel, 'muon_phi', self.tree['muon_phi'][imu])
			if self.tree['muon_isTight'][imu] == 1:  self.fill_hist(sel, 'muon_quality', 3)
			if self.tree['muon_isMedium'][imu] == 1: self.fill_hist(sel, 'muon_quality', 2)
			if self.tree['muon_isLoose'][imu] == 1:  self.fill_hist(sel, 'muon_quality', 1)
			else: self.fill_hist(sel, 'muon_quality', 0)

		for iel in range(len(self.tree['el_pt'])):
			self.fill_hist(sel, 'el_pt', self.tree['el_pt'][iel])
			self.fill_hist(sel, 'el_eta', self.tree['el_eta'][iel])
			self.fill_hist(sel, 'el_phi', self.tree['el_phi'][iel])
			if self.tree['el_LHTight'][iel] == 1:  self.fill_hist(sel, 'el_quality', 3)
			if self.tree['el_LHMedium'][iel] == 1: self.fill_hist(sel, 'el_quality', 2)
			if self.tree['el_LHLoose'][iel] == 1:  self.fill_hist(sel, 'el_quality', 1)
			else: self.fill_hist(sel, 'el_quality', 0)


	def _fill_selected_dv_histos(self, sel, do_lock=True):
		if not self.tree.is_data and not self.tree.not_hnl_mc:
			if self.MCEventType.isLNC: 
				sel =  sel + "_LNC" 
			if self.MCEventType.isLNV:
				sel =  sel + "_LNV" 


		if self._locked < FILL_LOCKED or not do_lock:
			# these are the histograms you only want to fill ONCE per DV

			# sel refers to the last selection that was applied
			for itrk in range(self.tree.ntrk):  # loop over tracks
				self.fill_hist(sel, 'DV_trk_pt', self.tree.dv('trk_pt_wrtSV')[itrk], fill_ntuple=False)
				self.fill_hist(sel, 'DV_trk_eta', self.tree.dv('trk_eta_wrtSV')[itrk], fill_ntuple=False)
				self.fill_hist(sel, 'DV_trk_phi', self.tree.dv('trk_phi_wrtSV')[itrk], fill_ntuple=False)
				self.fill_hist(sel, 'DV_trk_d0', self.tree.dv('trk_d0')[itrk], fill_ntuple=False)
				self.fill_hist(sel, 'DV_trk_z0', self.tree.dv('trk_z0')[itrk], fill_ntuple=False)
				self.fill_hist(sel, 'DV_trk_charge', self.tree.dv('trk_charge')[itrk], fill_ntuple=False)
				self.fill_hist(sel, 'DV_trk_chi2', self.tree.dv('trk_chi2_toSV')[itrk], fill_ntuple=False)


			self.fill_hist(sel, 'DV_num_trks', self.tree.dv('ntrk'))
			self.fill_hist(sel, 'DV_x', self.tree.dv('x'))
			self.fill_hist(sel, 'DV_y', self.tree.dv('y'))
			self.fill_hist(sel, 'DV_z', self.tree.dv('z'))
			self.fill_hist(sel, 'DV_r', self.tree.dv('r'))
			self.fill_hist(sel, 'DV_distFromPV', self.tree.dv('distFromPV'))
			self.fill_hist(sel, 'DV_mass', self.tree.dv('mass'))
			self.fill_hist(sel, 'DV_pt', self.tree.dv('pt'))
			self.fill_hist(sel, 'DV_eta', self.tree.dv('eta'))
			self.fill_hist(sel, 'DV_phi', self.tree.dv('phi'))
			self.fill_hist(sel, 'DV_minOpAng', self.tree.dv('minOpAng'))
			self.fill_hist(sel, 'DV_maxOpAng', self.tree.dv('maxOpAng'))
			self.fill_hist(sel, 'DV_charge', self.tree.dv('charge'))
			self.fill_hist(sel, 'DV_chi2', self.tree.dv('chi2'))
			# self.fill_hist(sel, 'DV_chi2_assoc', self.tree.dv('chi2_assoc'))
			self.fill_hist(sel, 'DV_alpha', selections.Alpha(self.tree).alpha)

			# kshort stuff
			track_sum = selections.SumTrack(self.tree)
			self.fill_hist(sel, 'DV_sum_track_pt', track_sum.sum_track_pt)
			self.fill_hist(sel, 'DV_sum_track_pt_wrt_pv', track_sum.sum_track_pt_wrt_pv)
			self.fill_hist(sel, 'DV_sum_track_pt_diff', track_sum.sum_track_pt_wrt_pv - track_sum.sum_track_pt)
			self.fill_hist(sel, 'DV_sum_track_charge', track_sum.sum_track_charge)

			if sel == self.saveNtuples or self.saveNtuples == 'allcuts': 
				if self.MCEventType.isLNC: 
					self.micro_ntuples["LNC_"+sel].fill()
					self.micro_ntuples["LNC_plus_LNV_"+sel].fill()
				elif self.MCEventType.isLNV: 
					self.micro_ntuples["LNV_"+sel].fill()
					self.micro_ntuples["LNC_plus_LNV_"+sel].fill()
				else: self.micro_ntuples[sel].fill()


	#########################################################################################################################
	# Define new cuts you want to apply here. This will overwrite whatever cuts are defined in the parent analysis class.
	#########################################################################################################################
	def _mass_window_cut(self):
		mass_min = .4977 - .01  # kshort mass +/- epsilon
		mass_max = .4977 + .01
		dv_mass = self.tree.dv('mass')
		return (dv_mass > mass_min) and (dv_mass < mass_max)

	def _alpha_cut(self):
		alpha_sel = selections.Alpha(self.tree)
		return alpha_sel.passes()

	def preSelection(self):
		self.passed_preselection_cuts = False
		self.passed_alpha_cut = False
		self.passed_mass_window_cut = False
		# fill cutflow bin
		self._fill_cutflow(0)

		# Check to make sure event has a PV otherwise throw event away (this happens very rarely with data).
		if self._pv_cut():
			self.h['CutFlow'][self.ch].Fill(1)
		else:
			return

		if self.do_invert_prompt_lepton_cut:
			if self._invert_prompt_lepton_cut():
				self.h['CutFlow'][self.ch].Fill(2)
			else:
				return

		self.passed_preselection_cuts = True
		self._fill_leptons()

	def DVSelection(self):
		######################################################################################################
		# DV Selection is any cuts that are done per DV
		# Current cuts include: fiducial vol, ntrack, OS, DVtype, track quality, cosmic veto, mlll, mDV
		######################################################################################################
		# Fill all the histograms with ALL DVs (this could be more that 1 per event). Useful for vertexing efficiency studies.
		self._fill_selected_dv_histos("all", do_lock=False)

		# only do the DV selection if the preselction was passed for the event.
		if not self.passed_preselection_cuts:
			return

		# self._fill_selected_dv_histos("presel", do_lock=False)

		if self.do_alpha_cut:
			if self._alpha_cut():
				self.h['CutFlow'][self.ch].Fill(3)
				self.passed_alpha_cut = True
			else:
				return
		# self._fill_selected_dv_histos("alpha", do_lock=False)

		if self.do_mass_window_cut:
			if self._mass_window_cut():
				self.h['CutFlow'][self.ch].Fill(4)
				self.passed_mass_window_cut = True
			else:
				return
		# self._fill_selected_dv_histos("mass", do_lock=False)

		self._fill_selected_dv_histos("sel", do_lock=False)
		# self._fill_selected_dv_ntuples("sel", do_lock=False)



class BEAnalysis(Analysis):
	def __init__(self, name, tree, vtx_container, selections, outputFile, saveNtuples, weight_override=None):

		Analysis.__init__(self, name, tree, vtx_container, selections, outputFile, saveNtuples, weight_override)
		self.logger.info('Running Background Estimate Analysis Cuts')

		# Define cutflow histogram "by hand"
		self.cutflow_dir = self.ch + '/CutFlow/'
		self.observables.histogram_dict[self.cutflow_dir+ 'CutFlow'] = ROOT.TH1D('CutFlow', 'CutFlow', 17, -0.5, 16.5)
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
			self.CutFlow.GetXaxis().SetBinLabel(6, "no plep overlap with DV")
		if self.do_invert_prompt_lepton_cut:
			self.CutFlow.GetXaxis().SetBinLabel(6, "invert prompt lepton")
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
		if self.do_same_event_cut:
			self.CutFlow.GetXaxis().SetBinLabel(11, "Same-Event")
		if self.do_different_event_cut:
			self.CutFlow.GetXaxis().SetBinLabel(11, "Different-Events")
		if self.do_dv_type_cut:
			self.CutFlow.GetXaxis().SetBinLabel(12, "%s DV" % self.dv_type)
		if self.do_mat_veto_cut:
			self.CutFlow.GetXaxis().SetBinLabel(13, "mat. veto")
		if self.do_cosmic_veto_cut:
			self.CutFlow.GetXaxis().SetBinLabel(14, "cosmic veto")
		if self.do_dv_mass_cut:
			self.CutFlow.GetXaxis().SetBinLabel(15, "m_{DV}")
		if self.do_track_quality_cut:
			self.CutFlow.GetXaxis().SetBinLabel(16, "{}-lepton DV".format(self.track_quality))
		if self.do_trilepton_mass_cut:
			self.CutFlow.GetXaxis().SetBinLabel(17, "m_{lll}")

	def DVSelection(self):
		######################################################################################################
		# DV Selection is any cuts that are done per DV
		# Current cuts include: fiducial vol, ntrack, OS, DVtype, track quality, cosmic veto, mlll, mDV
		######################################################################################################

		# Fill all the histograms with ALL DVs (this could be more that 1 per event). Useful for vertexing efficiency studies.
		self._fill_all_dv_histos()

		# There is an extra bit of logic here since we look at several DVs
		# but only want to fill the cutflow once per event
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

		if self.do_CR: # protect against saving OS DV when youre not looking in the CR
			self._fill_selected_dv_histos("2trk")
			OS_sel = selections.ChargeDV(self.tree, sel="OS").passes()
			SS_sel = selections.ChargeDV(self.tree, sel="SS").passes()
			if OS_sel:
				self._fill_selected_dv_histos("allOS") # save OS histograms
			elif SS_sel:
				self._fill_selected_dv_histos("allSS") # save SS histograms

		if self.do_opposite_sign_cut or self.do_same_sign_cut:

			if self._charge_cut():
				if not self.passed_charge_cut:
					self._fill_cutflow(9)
					self.passed_charge_cut = True
			else:
				return

		if self.do_same_event_cut or self.do_different_event_cut:

			if self._be_event_type_cut():
				if not self.passed_be_event_type_cut:
					self._fill_cutflow(10)
					self.passed_be_event_type_cut = True
			else:
				return

		self._fill_selected_dv_histos(self.be_region)

		if self.do_dv_type_cut:
			if self._dv_type_cut():
				if not self.passed_dv_type_cut:
					self._fill_cutflow(11)
					self.passed_dv_type_cut = True
			else:
				return

		self._fill_selected_dv_histos("DVtype")

		if self.do_mat_veto_cut:
			if self._mat_veto_cut():
				if not self.passed_mat_veto_cut:
					self._fill_cutflow(12)
					self.passed_mat_veto_cut = True
			else:
				return

		self._fill_selected_dv_histos("mat_veto")

		if self.do_cosmic_veto_cut:
			if self._cosmic_veto_cut():
				if not self.passed_cosmic_veto_cut:
					self._fill_cutflow(13)
					self.passed_cosmic_veto_cut = True
			else:
				return

		self._fill_selected_dv_histos("cosmic")

		if self.do_dv_mass_cut:
			if self._dv_mass_cut():
				if not self.passed_dv_mass_cut:
					self._fill_cutflow(14)
					self.passed_dv_mass_cut = True
			else:
				return
		self._fill_selected_dv_histos("mDV")


		if self.do_track_quality_cut:
			if self._track_quality_cut():
				if not self.passed_track_quality_cut:
					self._fill_cutflow(15)
					self.passed_track_quality_cut = True
			else:
				return

		self._fill_selected_dv_histos("trkqual")


		if self.do_trilepton_mass_cut:
			if self._trilepton_mass_cut():
				if not self.passed_trilepton_mass_cut:
					self._fill_cutflow(16)
					self.passed_trilepton_mass_cut = True
			else:
				return
		# self._fill_selected_dv_histos("mlll")


		# Fill all the histograms with only selected DVs. (ie. the ones that pass the full selection)
		self._fill_selected_dv_histos("sel")
