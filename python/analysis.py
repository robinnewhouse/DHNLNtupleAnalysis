# from ROOT import*
import ROOT
import numpy as np
import os
import sys
import helpers
import selections
import observables
import ntuples


UNLOCKED = 0
SELECTION_LOCKED = 1
FILL_LOCKED = 2


class Analysis(object):
	def __init__(self, name, tree, vtx_container, selection_list, outputFile, saveNtuples, debug_level,weight_override=None):
		self.logger = helpers.getLogger('dHNLAnalysis.analysis', level=debug_level)
		selections.set_debug_level(debug_level)
		self.name = name
		self.weight_override = weight_override
		self.sel = selection_list
		self.outputFile = outputFile
		self.fi = ROOT.TFile.Open(outputFile, 'update')
		self.ch = vtx_container
		self.histSuffixes = [self.ch]
		self.h = {}
		self.micro_ntuples = {}
		self.tree = tree
		# self.weight = 1
		self._locked = UNLOCKED
		# create an instance of Observables to store histograms
		self.observables = observables.Observables()

		self.events_with_trig_match_plep = 0
		self.events_with_trig_match_dlep = 0
		self.events_with_trig_match_both_pdlep =0

		# setting all the relevant variables for the cuts based on the input selections
		# trigger cut
		if 'alltriggers' in self.sel:
			self.trigger = 'all'
			self.do_trigger_cut = True
		else:
			if 'CR' not in self.sel:
				self.logger.warn('You did not specify a trigger configuration for this channel. Skipping trigger selection.')
			self.do_trigger_cut = False

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
	
			
		if 'CR' in self.sel:  # DO NOT CHANGE THESE CUTS OR YOU MIGHT UNBLIND DATA!!!
			self.do_CR = True
			self.fakeAOD = False
			self.do_trigger_cut = False  # do not apply trigger cut
			self.do_invert_trigger_cut = False  # do not apply inverted trigger cut
			self.do_filter_cut = False  # do not apply filter cut
			self.do_prompt_lepton_cut = False  # do not apply prompt lepton cut
			self.do_invert_prompt_lepton_cut = True  # invert prompt lepton cut
			self.logger.info('You are setup up to look in the inverted prompt lepton control region!')
		elif "CR_BE" in self.sel: # if running on fakeAOD that already has CR cuts applied (be careful with this setting!!!!!)
			self.fakeAOD = True
			self.do_CR = True
			self.do_trigger_cut = False  # do not apply trigger cut
			self.do_invert_trigger_cut = False  # do not apply inverted trigger cut
			self.do_filter_cut = False  # do not apply filter cut
			self.do_prompt_lepton_cut = False  # do not apply prompt lepton cut
			self.do_invert_prompt_lepton_cut = False  # do not apply inverted prompt lepton cut
			self.logger.info('You are running on a fakeAOD created from events in the invertex prompt lepton control region!')
		else:
			self.do_CR = False
			self.fakeAOD = False
			self.do_invert_prompt_lepton_cut = False
			self.do_invert_trigger_cut = False

		# alpha cut
		self.do_alpha_cut = 'alpha' in self.sel

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
		
		
		if self.do_CR: 
			if self.do_opposite_sign_cut == True and self.do_same_event_cut == True: 
				self.be_region = "RegionC"
			elif self.do_opposite_sign_cut == True and self.do_different_event_cut == True: 
				self.be_region = "RegionD"
			elif self.do_same_sign_cut == True and self.do_same_event_cut == True: 
				self.be_region = "RegionCprime"
			elif self.do_same_sign_cut == True and self.do_different_event_cut == True: 
				self.be_region = "RegionDprime"
		else: 
			if self.do_opposite_sign_cut == True and self.do_same_event_cut == True: 
				self.be_region = "RegionA"
			elif self.do_opposite_sign_cut == True and self.do_different_event_cut == True: 
				self.be_region = "RegionB"
			elif self.do_same_sign_cut == True and self.do_same_event_cut == True: 
				self.be_region = "RegionAprime"
			elif self.do_same_sign_cut == True and self.do_different_event_cut == True: 
				self.be_region = "RegionBprime"
		if "CR_BE" in self.sel: 
			self.saveNtuples = self.be_region # update saveNtuples selection
		else:
			self.saveNtuples = saveNtuples

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
		if not self.do_trilepton_mass_cut  and "CR" not in self.sel:
			self.logger.warn('You did not add a mlll cut for this channel. Skipping tri-lepton mass selection.')
		# material veto cut
		self.do_mat_veto_cut = "matveto" in self.sel

		# DV mass cut
		self.do_dv_mass_cut = 'DVmass' in self.sel
		if not self.do_dv_mass_cut and "CR" not in self.sel:
			self.logger.warn('You did not add a DVmass cut for this channel. Skipping displaced vertex mass selection.')

		# HNL mass cut
		self.do_HNL_mass_cut = 'HNLmass' in self.sel

		# HNL pT cut
		self.do_HNL_pt_cut = 'HNLpt' in self.sel

		self.check_input_consistency()

	# Getter helper functions
	def get_dv(self, key):
		return self.tree['secVtx_{}_{}'.format(self.ch, key)]

	def get(self, key):
		return self.tree[key]

	# hist filling helper functions
	def fill_hist(self, selection, hist_name, variable_1, variable_2=None, fill_ntuple=True,beAnalysis=False):
		"""
		A helper function for filling registered histograms
		:param selection: the step of selection the analysis it at. May be "None" in which case there will be no prefix.
		:param hist_name: base name of the histogram. When saved, a prefix and suffix will be appended.
		:param variable_1: variable you want to fill the histogram with.
		:param variable_2: if histogram is 2d, variable you want to fill the second axis of the histogram with
		:param fill_ntuple: set to True if you want to simultaneously fill an ntuple with this variable
		"""

		# define here the directory structure where this histogram is stored.
		directory = '{ch}/{selection}/'.format(ch=self.ch, selection=selection)
		# add LNC or LNV directory for simulation. Not used for data.
		if self.MCEventType.isLNC: directory += 'LNC/'
		if self.MCEventType.isLNV: directory += 'LNV/'

		self.observables.fill_hist(directory, hist_name, variable_1, variable_2, self.weight)

		# Unless suppressed, fill the corresponding micro-ntuple with the variable
		# Will not fill variables from 2D histograms to prevent double-counting
		# TODO Can we clean this up in some way?
		save_sel = self.saveNtuples == selection or 'truth_'+self.saveNtuples == selection or self.saveNtuples == 'allcuts'
		if fill_ntuple and variable_2 is None and save_sel:
			# Note: selection and hist_name will be overridden by full_name
			# Need selection to define ntuple tree
			# TODO redo this method to use the directory correctly
			if self.MCEventType.isLNC: self.fill_ntuple(selection, hist_name, variable_1,MCtype="LNC")
			elif self.MCEventType.isLNV: self.fill_ntuple(selection, hist_name, variable_1,MCtype="LNV")
			else: self.fill_ntuple(selection, hist_name, variable_1)
	def fill_ntuple(self, selection, ntuple_name, variable,MCtype=None, full_name=""):
		"""
		A helper function for filling micro-ntuples. Often called from the fill_hist function.
		If you are using this in you analysis,
		please check that it is not also being called by fill_hist to prevent double-counting.
		:param selection: the step of selection the analysis it at. May be "None" in which case there will be no prefix.
		:param ntuple_name: base name of the ntuple. When saved, a prefix and suffix will be appended.
		:param variable: variable you want to fill the histogram with.
		:param full_name: override the automatic naming of the ntuple.
		"""
		if not selection:
			raise ValueError("You must indicate a selection in order to store the ntuple. Use 'all' if no selection.")
		# Retrieve the ntuple for this selection. If it doesn't exist, create it.
		if MCtype !=None: 
			selection = MCtype + "_" + 	selection 

		if selection not in self.micro_ntuples:
			self.micro_ntuples[selection] = ntuples.Ntuples('ntuples_{}_{}'.format(selection, self.ch))  # temp name. not written
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
				if self.do_trilepton_mass_cut or self.do_HNL_mass_cut or do_HNL_pt_cut:
					self.logger.warning("You cannot cut on mlll, HNLpt or HNLm without first selecting a prompt lepton. Apply a prompt lepton cut!")
					sys.exit(1)  # abort because of error

		if self.do_opposite_sign_cut and self.do_same_sign_cut:
			self.logger.error("These cuts are mutually exclusive. You will get zero events!")
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
		for hist in self.observables.histogram_dict.values(): 
			hist.SetBinContent(hist.GetNbinsX(), hist.GetBinContent(hist.GetNbinsX()) + hist.GetBinContent(hist.GetNbinsX() + 1)) # merge overflow into last bin
			hist.SetBinContent(1, hist.GetBinContent(1) + hist.GetBinContent(0)) # merge underflow into first bin


		# make acceptance Histograms 
		# TOD: it doesnt looks like on data the acceptance histograms are working as expected. -DT
		if not self.tree.is_data and not self.tree.notHNLmc:
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
			trig_match = selections.TriggerMatching_prompt(self.tree, self.plep, self.plep_sel.plep_Index)
			if trig_match.plep_isTrigMatched:
				self.events_with_trig_match_plep = self.events_with_trig_match_plep + 1


			self.fill_hist('all', 'plep_pt', self.plep_sel.plepVec.Pt())
			self.fill_hist('all', 'plep_eta', self.plep_sel.plepVec.Eta())
			self.fill_hist('all', 'plep_phi', self.plep_sel.plepVec.Phi())
			self.fill_hist('all', 'plep_d0', self.plep_sel.plepd0)
			self.fill_hist('all', 'plep_z0', self.plep_sel.plepz0)
		return self.plep_sel.passes() # full plep selection find the highest pt plep that doesnt overlap with any DVs

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
		return True # fill this in when tree is ready!


	def _dv_type_cut(self):
		dv_sel = selections.DVtype(self.tree, dv_type=self.dv_type,fakeAOD = self.fakeAOD)
		if dv_sel.passes():
			if self.fakeAOD == False: 
				trig_match = selections.TriggerMatching_disp(self.tree, self.dv_type, dv_sel.dMu_Index, dv_sel.dEl_Index)
				if trig_match.dlep_isTrigMatched:
					self.events_with_trig_match_dlep = self.events_with_trig_match_dlep + 1
				count_trig_match_disp_event = True

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

		muons = helpers.Tracks(self.tree,self.fakeAOD)
		muons.getMuons()
		muVec = muons.lepVec

		electrons = helpers.Tracks(self.tree,self.fakeAOD)
		electrons.getElectrons()
		elVec = electrons.lepVec

		mlll_sel = selections.Mlll(dv_type=self.dv_type, plep=plep_vec, dMu=muVec, dEl=elVec)
		return mlll_sel.passes()

	def _mat_veto_cut(self):
		return selections.Mat_veto(self.tree).passes()


	def _dv_mass_cut(self):
		dv_mass_sel = selections.DVmass(self.tree, dvmasscut=2)
		return dv_mass_sel.passes()

	def _multitrk_2lep_cut(self):
		if self.tree.dv('ntrk') >= 2:  # 2+ trk vertex
			dv_type_sel = selections.DVtype(self.tree, dv_type=self.dv_type,fakeAOD = self.fakeAOD)
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
		self.passed_be_event_type_cut = False
		self.passed_track_quality_cut = False
		self.passed_cosmic_veto_cut = False
		self.passed_trilepton_mass_cut = False
		self.passed_mat_veto_cut = False
		self.passed_dv_mass_cut = False
		self.passed_HNL_mass_cut = False
		self.passed_HNL_pt_cut = False
		self.count_trig_match_disp_event = False


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

		self._fill_leptons()

		if not self.tree.is_data and not self.tree.notHNLmc:
			self._fill_truth_histos(sel='truth/all')
			if self.MCEventType.isLNC: 
				self.CutFlow_LNC.SetBinContent(1, self.tree.all_entries/2)  # all events
			if self.MCEventType.isLNV:
				self.CutFlow_LNV.SetBinContent(1, self.tree.all_entries/2)  # all events
		
		self.CutFlow.SetBinContent(1, self.tree.all_entries)  # all events
		

		######################################################################################################
		# Selection code is deisgned so that it will pass the selection only if the cut true or cut is unused
		# ex. passTrigger is true if the trigcut is true OR if trigcut is not used)
		######################################################################################################

		if self.do_trigger_cut:
			if self._trigger_cut():
				# Fill the plot at the specified bin
				self._fill_cutflow(1)
			else:
				return
		
		if self._pv_cut(): #Check to make sure event has a PV otherwise throw event away (this happens very rarely with data).
			self._fill_cutflow(2)
		else:
			return

		if self.tree.is_data: # when running on data skip over any events without any DVs to speed up running
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
		if not self.tree.is_data and not self.tree.notHNLmc:
			self._fill_truth_histos(sel='truth/presel')

	def calculate_event_weight(self):
		# MC re-weighting to include spin correlations and fix lepton ordering bug
		self.MCEventType = selections.MCEventType(self.tree) # if data then MCEventType weight defaults to 1
		# calculate mass lifetime weight 
		self.mass_lt_weight = helpers.get_mass_lt_weight(self.tree,logger=self.logger, both_lnc_lnv=False) 
		# self.mass_lt_weight = helpers.get_mass_lt_weight(self.tree.mass, self.tree.ctau,lnv=self.MCEventType.isLNV)  
		self.logger.debug('Event weight for this signal sample is: {}'.format(self.mass_lt_weight))
		if self.weight_override == None: 
			self.weight = self.mass_lt_weight*self.MCEventType.weight 
		else: 
			self.weight = self.weight_override
		

	def DVSelection(self):
		raise NotImplementedError("Please implement this method in your own Analysis subclass")

	def _fill_cutflow(self, nbin):
		if not self.tree.is_data and not self.tree.notHNLmc:
			if self.MCEventType.isLNC:
				self.CutFlow_LNC.Fill(nbin)
			if self.MCEventType.isLNV:
				self.CutFlow_LNV.Fill(nbin)
			self.CutFlow.Fill(nbin)
		else:
			self.CutFlow.Fill(nbin)

	def _fill_multitrk_histos(self):
		self.fill_hist('2lepMultitrk', 'num_trks', self.tree.dv('ntrk'))
		muons = helpers.Tracks(self.tree,self.fakeAOD)
		muons.getMuons()
		if muons.lepisAssoc[0] == 1 and muons.lepisAssoc[1] == 1:
			self.fill_hist('2lepMultitrk', 'bothmuon_isAssociated', 1)
		else:
			self.fill_hist('2lepMultitrk', 'bothmuon_isAssociated', 0)
		if muons.lepisAssoc[0] == 0 and muons.lepisAssoc[1] == 0:
			self.fill_hist('2lepMultitrk', 'nomuon_isAssociated', 1)
			tracks = helpers.Tracks(self.tree,self.fakeAOD)
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
			if self.MCEventType.isLNC: self.micro_ntuples["LNC_"+sel].fill()
			elif self.MCEventType.isLNV: self.micro_ntuples["LNV_"+sel].fill()
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
		
		# print truth_info.W_charge
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
			self.fill_hist(sel, 'DV_mass', DV_4vec.M())
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
			self.fill_hist(sel, 'lep1_trk_pt', self.MCEventType.p_2.Pt())
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

			self.fill_hist(sel, 'dlep1_pt', disp_lep[0].Pt())
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


			for itrk in range(2):
				self.fill_hist(sel, 'DV_trk_pt', truth_info.trkVec[itrk].Pt(), fill_ntuple=False)
				self.fill_hist(sel, 'DV_trk_eta', truth_info.trkVec[itrk].Eta(), fill_ntuple=False)
				self.fill_hist(sel, 'DV_trk_phi', truth_info.trkVec[itrk].Phi(), fill_ntuple=False)
				self.fill_hist(sel, 'DV_trk_d0',truth_info.dTrk_d0[itrk], fill_ntuple=False)

			for iel in range(len(truth_info.dEl)):
				self.fill_hist(sel, 'DV_El_pt', truth_info.dEl[iel].Pt(), fill_ntuple=False)
				self.fill_hist(sel, 'DV_El_eta', truth_info.dEl[iel].Eta(), fill_ntuple=False)
				self.fill_hist(sel, 'DV_El_phi', truth_info.dEl[iel].Phi(), fill_ntuple=False)
			
			for imu in range(len(truth_info.dMu)):
				self.fill_hist(sel, 'DV_Mu_pt', truth_info.dMu[imu].Pt(), fill_ntuple=False)
				self.fill_hist(sel, 'DV_Mu_eta', truth_info.dMu[imu].Eta(), fill_ntuple=False)
				self.fill_hist(sel, 'DV_Mu_phi', truth_info.dMu[imu].Phi(), fill_ntuple=False)

			# TODO: figure out a ntuple scheme that can store these variables as well
		if sel == self.saveNtuples or self.saveNtuples == 'allcuts': 
			if self.MCEventType.isLNC: self.micro_ntuples["LNC_"+sel].fill()
			elif self.MCEventType.isLNV: self.micro_ntuples["LNV_"+sel].fill()
			else: self.micro_ntuples[sel].fill()


	def _fill_selected_dv_histos(self, sel, do_lock=True):

		if self._locked < FILL_LOCKED and do_lock:
			# these are the histograms you only want to fill ONCE per DV
			# sel refers to the last selection that was applied

			# fill event weight. storing this per dv as weights include dv scale factor.
			self.fill_hist(sel, 'DV_weight', self.weight)


			tracks = helpers.Tracks(self.tree,self.fakeAOD)
			tracks.getTracks()
			trkVec = tracks.lepVec

			muons = helpers.Tracks(self.tree,self.fakeAOD)
			muons.getMuons()
			muVec = muons.lepVec

			electrons = helpers.Tracks(self.tree,self.fakeAOD)
			electrons.getElectrons()
			elVec = electrons.lepVec


			# fill histograms that require a prompt lepton to be identified
			if self.do_prompt_lepton_cut:
				plep_vec = self.plep_sel.plepVec
				plepd0 = self.plep_sel.plepd0
				plepz0 = self.plep_sel.plepz0
				plepcharge = self.plep_sel.plepcharge

				self.fill_hist(sel, 'plep_pt', plep_vec.Pt())
				self.fill_hist(sel, 'plep_eta', plep_vec.Eta())
				self.fill_hist(sel, 'plep_phi', plep_vec.Phi())
				self.fill_hist(sel, 'plep_d0', plepd0)
				self.fill_hist(sel, 'plep_z0', plepz0)
				self.fill_hist(sel, 'plep_charge', plepcharge)

				if tracks.ntracks == 2:
					Mlll = selections.Mlll(dv_type=self.dv_type, plep=plep_vec, dMu=muVec, dEl=elVec)
					Mhnl = selections.Mhnl(self.tree, self.dv_type, plep=plep_vec, dMu=muVec,dEl=elVec)
					# Mhnl_old = selections.Mhnl_old(self.tree, plep=plep_vec, trks=trkVec)
					Mhnl_fixWmass = selections.Mhnl(self.tree, self.dv_type, plep=plep_vec, dMu=muVec,dEl=elVec,fixWMass=True)

					self.fill_hist(sel, 'mvis', Mlll.mlll )
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
						self.fill_hist(sel, 'DV_redmass', self.tree.dv('mass')/dR)
						self.fill_hist(sel, 'DV_redmassvis', Mlll.mlll/dR)
						self.fill_hist(sel, 'DV_redmassHNL', Mhnl.mhnl/dR)

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

				self.fill_hist(sel, 'DV_trk_max_chi2_toSV', max(self.tree.dv('trk_chi2_toSV')[0],self.tree.dv('trk_chi2_toSV')[1] ) )
				self.fill_hist(sel, 'DV_trk_min_chi2_toSV', min(self.tree.dv('trk_chi2_toSV')[0],self.tree.dv('trk_chi2_toSV')[1] ) )
				self.fill_hist(sel, 'DV_trk_max_d0_wrtSV', max(self.tree.dv('trk_d0_wrtSV')[0],self.tree.dv('trk_d0_wrtSV')[1] ) )
				self.fill_hist(sel, 'DV_trk_min_d0_wrtSV', min(self.tree.dv('trk_d0_wrtSV')[0],self.tree.dv('trk_d0_wrtSV')[1] ) )
				self.fill_hist(sel, 'DV_trk_max_errd0_wrtSV', max(self.tree.dv('trk_errd0_wrtSV')[0],self.tree.dv('trk_errd0_wrtSV')[1] ) )
				self.fill_hist(sel, 'DV_trk_min_errd0_wrtSV', min(self.tree.dv('trk_errd0_wrtSV')[0],self.tree.dv('trk_errd0_wrtSV')[1] ) )
				self.fill_hist(sel, 'DV_trk_max_z0_wrtSV', max(self.tree.dv('trk_z0_wrtSV')[0],self.tree.dv('trk_z0_wrtSV')[1] ) )
				self.fill_hist(sel, 'DV_trk_min_z0_wrtSV', min(self.tree.dv('trk_z0_wrtSV')[0],self.tree.dv('trk_z0_wrtSV')[1] ) )
				self.fill_hist(sel, 'DV_trk_max_errz0_wrtSV', max(self.tree.dv('trk_errz0_wrtSV')[0],self.tree.dv('trk_errz0_wrtSV')[1] ) )
				self.fill_hist(sel, 'DV_trk_min_errz0_wrtSV', min(self.tree.dv('trk_errz0_wrtSV')[0],self.tree.dv('trk_errz0_wrtSV')[1] ) )

				DV_mumu = selections.DVtype(self.tree, dv_type="mumu",fakeAOD = self.fakeAOD).passes()
				DV_ee = selections.DVtype(self.tree, dv_type="ee",fakeAOD = self.fakeAOD).passes()
				DV_emu = selections.DVtype(self.tree, dv_type="emu",fakeAOD = self.fakeAOD).passes()
				DV_1lep = (len(muVec) ==  1 and len(elVec) == 0) or (len(muVec) ==  0 and len(elVec) == 1)

				self.fill_hist(sel, 'DV_mumu', DV_mumu)
				self.fill_hist(sel, 'DV_ee', DV_ee)
				self.fill_hist(sel, 'DV_emu', DV_emu)
				self.fill_hist(sel, 'DV_1lep', DV_1lep)
				

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

					self.fill_hist(sel, 'DV_trk_1_pt', self.tree.dv('trk_pt_wrtSV')[0])
					self.fill_hist(sel, 'DV_trk_1_eta', self.tree.dv('trk_eta_wrtSV')[0])
					self.fill_hist(sel, 'DV_trk_1_phi', self.tree.dv('trk_phi_wrtSV')[0])
					self.fill_hist(sel, 'DV_trk_1_d0', self.tree.dv('trk_d0')[0])
					self.fill_hist(sel, 'DV_trk_1_z0', self.tree.dv('trk_z0')[0])
					self.fill_hist(sel, 'DV_trk_1_charge', self.tree.dv('trk_charge')[0])
					self.fill_hist(sel, 'DV_trk_1_chi2', self.tree.dv('trk_chi2')[0])
					self.fill_hist(sel, 'DV_trk_1_isSelected', self.tree.dv('trk_isSelected')[0])
					self.fill_hist(sel, 'DV_trk_1_isAssociated', self.tree.dv('trk_isAssociated')[0])
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

					self.fill_hist(sel, 'DV_trk_1_pt', self.tree.dv('trk_pt_wrtSV')[1])
					self.fill_hist(sel, 'DV_trk_1_eta', self.tree.dv('trk_eta_wrtSV')[1])
					self.fill_hist(sel, 'DV_trk_1_phi', self.tree.dv('trk_phi_wrtSV')[1])
					self.fill_hist(sel, 'DV_trk_1_d0', self.tree.dv('trk_d0')[1])
					self.fill_hist(sel, 'DV_trk_1_z0', self.tree.dv('trk_z0')[1])
					self.fill_hist(sel, 'DV_trk_1_charge', self.tree.dv('trk_charge')[1])
					self.fill_hist(sel, 'DV_trk_1_chi2', self.tree.dv('trk_chi2')[1])
					self.fill_hist(sel, 'DV_trk_1_isSelected', self.tree.dv('trk_isSelected')[1])
					self.fill_hist(sel, 'DV_trk_1_isAssociated', self.tree.dv('trk_isAssociated')[1])



			# fill standard track variable histograms
			for i in xrange(tracks.ntracks):
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
				self.fill_hist(sel, 'DV_trk_nSiHits', self.tree.dv('trk_nSCTHits')[i]+self.tree.dv('trk_nPixelHits')[i])
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

			#is truth matched: 
			if not self.tree.is_data:
				maxlinkTruth_score = self.tree.dv('maxlinkTruth_score')
				maxlinkTruth_parent_pdgId = abs(self.tree.dv('maxlinkTruth_parent_pdgId'))
				is_truth_matched =  self.tree.dv('maxlinkTruth_score') > 0.75 and abs(self.tree.dv('maxlinkTruth_parent_pdgId')) == 50
				self.fill_hist(sel, 'DV_truth_matched', is_truth_matched)

			# compute alpha (3D angle between DV 3-momentum and rDV)
			dv = ROOT.TVector3( self.tree.dv('x'), self.tree.dv('y'),  self.tree.dv('z') )
			pv = ROOT.TVector3( self.tree['vertex_x'], self.tree['vertex_y'],  self.tree['vertex_z'])
			decayV = dv-pv

			dv_4vec = ROOT.TLorentzVector()
			dv_4vec.SetPtEtaPhiM(self.tree.dv('pt'), self.tree.dv('eta'),self.tree.dv('phi'), self.tree.dv('mass'))
			dv_mom_vec =  ROOT.TVector3( dv_4vec.Px(),  dv_4vec.Py(),  dv_4vec.Pz() )
			alpha = decayV.Angle(dv_mom_vec)

			self.fill_hist(sel, 'DV_alpha', alpha)

			if self.fakeAOD == False: 
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

			
			# fill TTree with ntuple information. Already set by fill_hist
			if sel == self.saveNtuples or self.saveNtuples == 'allcuts':  
				if self.MCEventType.isLNC: self.micro_ntuples["LNC_"+sel].fill()
				elif self.MCEventType.isLNV: self.micro_ntuples["LNV_"+sel].fill()
				else: self.micro_ntuples[sel].fill()

			if sel == "sel":
				self._locked = FILL_LOCKED  # this only becomes unlocked after the event loop finishes in makeHistograms so you can only fill one DV from each event.


class run2Analysis(Analysis):
	def __init__(self, name, tree, vtx_container, selections, outputFile, saveNtuples, debug_level,weight_override=None):
		
		Analysis.__init__(self, name, tree, vtx_container, selections, outputFile, saveNtuples, debug_level,weight_override)
		self.logger.info('Running  Full Run 2 Analysis cuts')

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
		if self.do_dv_type_cut:
			self.CutFlow.GetXaxis().SetBinLabel(11, "%s DV" % self.dv_type)
		if self.do_mat_veto_cut:
			self.CutFlow.GetXaxis().SetBinLabel(12, "mat. veto")
		if self.do_cosmic_veto_cut:
			self.CutFlow.GetXaxis().SetBinLabel(13, "cosmic veto")
		if self.do_dv_mass_cut:
			self.CutFlow.GetXaxis().SetBinLabel(14, "m_{DV}")
		if self.do_track_quality_cut:
			self.CutFlow.GetXaxis().SetBinLabel(15, "{}-lepton DV".format(self.track_quality))
		if self.do_trilepton_mass_cut:
			self.CutFlow.GetXaxis().SetBinLabel(16, "m_{lll}")
		self.CutFlow.GetXaxis().SetBinLabel(17, "truth matched")

		# Store LNC and LNV cutflows in the observables collection
		if not self.tree.is_data and not self.tree.notHNLmc: 
			self.CutFlow_LNV = self.CutFlow.Clone()
			self.CutFlow_LNC = self.CutFlow.Clone()
			self.CutFlow_LNV.SetName("CutFlow_LNV"+"_"+self.ch)
			self.CutFlow_LNC.SetName("CutFlow_LNC"+"_"+self.ch)
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_LNV'] = self.CutFlow_LNV
			self.observables.histogram_dict[self.cutflow_dir+'CutFlow_LNC'] = self.CutFlow_LNC

	def DVSelection(self):
		######################################################################################################
		# DV Selection is any cuts that are done per DV
		# Current cuts include: fiducial vol, ntrack, OS, DVtype, track quality, cosmic veto, mlll, mDV
		######################################################################################################

		# Fill all the histograms with ALL DVs (this could be more that 1 per event). Useful for vertexing efficiency studies.
		self._fill_all_dv_histos()

		# only do the DV selection if the preselction was passed for the event.
		if not self.passed_preselection_cuts:
			return

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
		if self.dv_type == "mumu":

			if self._multitrk_2lep_cut(): # no return becuase this is not an analysis cut, only used for studying S & B, only worked for uuu samples -DT
				self._fill_multitrk_histos()

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
				self._fill_selected_dv_histos("OS") # save OS histograms 
				if self._dv_type_cut(): 
					self._fill_selected_dv_histos("OS_DVtype") # save lepton OS histograms 
			elif SS_sel:
				self._fill_selected_dv_histos("SS") # save SS histograms 
				if self._dv_type_cut(): 
					self._fill_selected_dv_histos("SS_DVtype") # save lepton SS histograms 

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
			else:
				return
	
		self._fill_selected_dv_histos("DVtype")

		if self.do_mat_veto_cut:
			if self._mat_veto_cut():
				if not self.passed_mat_veto_cut:
					self._fill_cutflow(11)
					self.passed_mat_veto_cut = True
			else:
				return

		self._fill_selected_dv_histos("mat_veto")

		if self.do_cosmic_veto_cut:
			if self._cosmic_veto_cut():
				if not self.passed_cosmic_veto_cut:
					self._fill_cutflow(12)
					self.passed_cosmic_veto_cut = True
			else:
				return

		self._fill_selected_dv_histos("cosmic")

		if self.do_dv_mass_cut:
			if self._dv_mass_cut():
				if not self.passed_dv_mass_cut:
					self._fill_cutflow(13)
					self.passed_dv_mass_cut = True
			else:
				return
		self._fill_selected_dv_histos("mDV")

		
		if self.do_track_quality_cut:
			if self._track_quality_cut():
				if not self.passed_track_quality_cut:
					self._fill_cutflow(14)
					self.passed_track_quality_cut = True
			else:
				return

		self._fill_selected_dv_histos("trkqual")


		if self.do_trilepton_mass_cut:
			if self._trilepton_mass_cut():
				if not self.passed_trilepton_mass_cut:
					self._fill_cutflow(15)
					self.passed_trilepton_mass_cut = True
			else:
				return
		# self._fill_selected_dv_histos("mlll")
		 
		# Fill histos of truth-matched DVs
		if not self.tree.is_data and not self.tree.notHNLmc:
			if self._truth_match():
				self._fill_cutflow(16)
				# self.h['CutFlow'][self.ch].Fill(14)
				self._fill_selected_dv_histos("match")

		# Fill all the histograms with only selected DVs. (ie. the ones that pass the full selection)
		self._fill_selected_dv_histos("sel")



class KShort(Analysis):
	def __init__(self, name, tree, vtx_container, selections, outputFile, saveNtuples, debug_level,weight_override=None):
		Analysis.__init__(self, name, tree, vtx_container, selections, outputFile, saveNtuples, debug_level,weight_override)
		self.logger.info('Running KShort Analysis cuts', level=debug_level)

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
		if not self.tree.is_data and not self.tree.notHNLmc:
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
				if self.MCEventType.isLNC: self.micro_ntuples["LNC_"+sel].fill()
				elif self.MCEventType.isLNV: self.micro_ntuples["LNV_"+sel].fill()
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
		self.h['CutFlow'][self.ch].Fill(0)  # all

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
	def __init__(self, name, tree, vtx_container, selections, outputFile, saveNtuples, debug_level,weight_override=None):
		
		Analysis.__init__(self, name, tree, vtx_container, selections, outputFile, saveNtuples, debug_level,weight_override)
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

		# only do the DV selection if the preselction was passed for the event.
		if not self.passed_preselection_cuts:
			return

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

