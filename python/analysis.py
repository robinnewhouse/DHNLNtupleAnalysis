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
	def __init__(self, name, tree, vtx_container, selection_list, outputFile, saveNtuples, debug_level):
		self.logger = helpers.getLogger('dHNLAnalysis.analysis', level=debug_level)
		selections.set_debug_level(debug_level)
		self.name = name
		self.sel = selection_list
		self.outputFile = outputFile
		self.fi = ROOT.TFile.Open(outputFile, 'update')
		self.ch = vtx_container
		self.histSuffixes = [self.ch]
		self.h = {}
		self.micro_ntuples = {}
		self.tree = tree
		self.weight = 1
		self.saveNtuples = saveNtuples
		self._locked = UNLOCKED
		# create an instance of Observables to store histograms
		self.observables = observables.Observables()

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

		if 'CR' in self.sel:  # DO NOT CHANGE THESE CUTS OR YOU MIGHT UNBLIND DATA!!!
			self.do_CR = True
			self.do_trigger_cut = False  # do not apply trigger cut
			self.do_invert_trigger_cut = False  # do not apply inverted trigger cut
			self.do_filter_cut = False  # do not apply filter cut
			self.do_prompt_lepton_cut = False  # do not apply prompt lepton cut
			self.do_invert_prompt_lepton_cut = True  # invert prompt lepton cut
			self.logger.info('You are setup up to look in the inverted prompt lepton control region!')
		else:
			self.do_CR = False
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
		elif '2-tight' in self.sel:
			self.track_quality = '2-tight'
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
	def fill_hist(self, selection, hist_name, variable_1, variable_2=None, fill_ntuple=True):
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

		self.observables.fill_hist(directory, hist_name, variable_1, variable_2, self.tree.weight)

		# Unless suppressed, fill the corresponding micro-ntuple with the variable
		# Will not fill variables from 2D histograms to prevent double-counting
		# TODO Can we clean this up in some way?
		save_sel = self.saveNtuples == selection or 'truth_'+self.saveNtuples == selection or self.saveNtuples == 'allcuts'
		if fill_ntuple and variable_2 is None and save_sel:
			# Note: selection and hist_name will be overridden by full_name
			# Need selection to define ntuple tree
			# TODO redo this method to use the directory correctly
			self.fill_ntuple(selection, hist_name, variable_1)

	def fill_ntuple(self, selection, ntuple_name, variable, full_name=""):
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
				self.logger.error("You cannot calculate mlll, HNLpt or HNLm without selecting a prompt lepton!")
				sys.exit(1)  # abort because of error

		if self.do_opposite_sign_cut and self.do_same_sign_cut:
			self.logger.error("These cuts are mutually exclusive. You will get zero events!")
			sys.exit(1)  # abort because of error

	def unlock(self):
		self._locked = UNLOCKED

	def write(self):
		# Move ROOT to base directory
		self.fi.cd()

		# Store saved ntuple values to file
		# self.fi.mkdir(self.ch+'/ntuples', "Micro Ntuples " + self.ch)
		# self.fi.cd(self.ch+'/ntuples')
		[ntuple.write(self.ch+'_ntuples_'+key) for key, ntuple in self.micro_ntuples.items()]

		self.observables.write_histograms(root_file=self.fi)

		# # Make a subdirectory for vertex type. May be other channels in the future.
		# if not self.fi.FindObject(self.ch):
		# 	self.fi.mkdir(self.ch, "Analysis Channel " + self.ch)
		# # Move ROOT to the channel subdirectory
		# self.fi.cd(self.ch)

		# # Store saved histograms to file
		# # TODO: this should be saved in a different way. e.g. another level of dictionary. Parsing strings for variable names is not good.
		# for h_name in self.h:
		#
		# 	if self.h[h_name][self.ch].GetEntries() != 0:
		#
		# 		if not self.tree.is_data:
		# 			selection = h_name.split('_')[0]  # get LNC or LNV
		# 			EventType = '_'.join(h_name.split('_')[1:]).split('_')[0]  # get selection
		# 			sel_dir = self.ch + '/' +  selection + '/' + EventType
		# 			base_name = '_'.join(h_name.split('_')[2:])  # get base name
		# 		else:
		# 			selection = h_name.split('_')[0]  # get selection
		# 			sel_dir = self.ch + '/' + selection
		# 			base_name = '_'.join(h_name.split('_')[2:])  # get base name
		# 		# if "CutFlow" in h_name:
		# 		# 	print EventType
		# 		# 	print selection
		# 		# 	print sel_dir
		# 		# 	print base_name
		#
		# 		if not self.fi.GetDirectory(sel_dir):  # make TDirectory if necessary
		# 			self.fi.mkdir(sel_dir, "Analysis Selection " + selection)
		# 		self.fi.cd(sel_dir)  # change to TDirectory
		# 		self.h[h_name][self.ch].Write(base_name)  # save only the base name
		#
		# self.logger.info("Histograms written to {}".format(self.outputFile))

		self.fi.Close()

	def end(self):
		self.h['CutFlow_all_acceptance'] = {}
		self.h['CutFlow_all_acceptance'][self.ch] = self.CutFlow.Clone()
		self.h['CutFlow_all_acceptance'][self.ch].SetName("CutFlow_all_acceptance"+"_"+self.ch)
		self.h['CutFlow_all_acceptance'][self.ch].SetDirectory(0)
		if self.CutFlow.GetBinContent(1) != 0: # Protect against zero-division
			self.h['CutFlow_all_acceptance'][self.ch].Scale(1.0/self.CutFlow.GetBinContent(1))
		self.logger.info('Done with Channel("{}")'.format(self.ch))
		meta = []
		# if self.region:
		# 	meta.append('Region: {}'.format(self.region))
		# if self.period:
		# 	meta.append('Period: {}'.format(self.period))
		# if self.bcategory != None:
		# 	meta.append('B-tagging Categroy: {}'.format(self.bcategory))
		# self.logger.info('\t' + ' | '.join(meta))

		# gives warning messages if histograms are unfilled
		for histName in self.h:
			if self.h[histName][self.ch].GetEntries() == 0:
				self.logger.debug('\tUnfilled HIST({}<{}>)!'.format(histName, self.ch))

		# y_err = ROOT.Double()
		# for s in self.histSuffixes:
		# 	h = self.h['finalNEvents'][s]
		# 	self.logger.info('\tFinal NEvents <{syst}>: {y}'.format(syst = s, y = h.Integral()))
		# 	h = self.h['finalYields'][s]
		# 	self.logger.info('\tFinal Yields <{syst}>: {y:.4g}+/-{y_err:.2g}'.format(syst = s, y = h.IntegralAndError(0, h.GetNbinsX()+1, y_err), y_err = y_err))
		self.write()

		# Clean up memory
		del self.h
		del self.micro_ntuples


		# head, sep, tail = self._outputFile.partition('file://')
		# f = tail if head == '' else self._outputFile
		# f = tail if head == '' else self._outputFile
		# try:
		# 	os.rename(f + '.part', f)
		# except OSError as e:
		# 	self.logger.error(e, exc_info=True)


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
		self.plep_sel = selections.PromptLepton(self.tree, lepton=self.plep)
		# Add to histogram all prompt leptons that pass selection.
		# If _prompt_lepton_cut() is run after trigger and filter cut then those cuts will also be applied.
		if self.plep_sel.passes():
			self.fill_hist('all', 'plep_pt', self.plep_sel.plepVec.Pt())
			self.fill_hist('all', 'plep_eta', self.plep_sel.plepVec.Eta())
			self.fill_hist('all', 'plep_phi', self.plep_sel.plepVec.Phi())
			self.fill_hist('all', 'plep_d0', self.plep_sel.plepd0)
			self.fill_hist('all', 'plep_z0', self.plep_sel.plepz0)
		return self.plep_sel.passes()

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

	def _dv_type_cut(self):
		dv_sel = selections.DVtype(self.tree, dv_type=self.dv_type)
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

		mlll_sel = selections.Mlll(dv_type=self.dv_type, plep=plep_vec, dMu=muVec, dEl=elVec)
		return mlll_sel.passes()

	def _dv_mass_cut(self):
		dv_mass_sel = selections.DVmass(self.tree, dvmasscut=4)
		return dv_mass_sel.passes()

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
		self.passed_track_quality_cut = False
		self.passed_cosmic_veto_cut = False
		self.passed_trilepton_mass_cut = False
		self.passed_dv_mass_cut = False
		self.passed_HNL_mass_cut = False
		self.passed_HNL_pt_cut = False


	def preSelection(self):

		self.calculate_event_weight()

		######################################################################################################
		# Preselection are all the cuts that are requied per event
		# Current cuts include: trigger, filter, plepton, DV cut
		######################################################################################################

		#initialize the cut bools for every event
		self.initialize_cut_bools()

		self._fill_leptons()

		if not self.tree.is_data:
			self._fill_truth_histos(sel='truth/all')
			if self.MCEventType.isLNC: 
				self.CutFlow_LNC.SetBinContent(1, self.tree.all_entries)  # all events
			if self.MCEventType.isLNV:
				self.CutFlow_LNV.SetBinContent(1, self.tree.all_entries)  # all events
			self.CutFlow.SetBinContent(1, self.tree.all_entries)  # all events
		

		######################################################################################################
		# Selection code is deisgned so that it will pass the selection only if the cut true or cut is unused
		# ex. passTrigger is true if the trigcut is true OR if trigcut is not used)
		######################################################################################################

	

		if self.do_trigger_cut:
			if self._trigger_cut():
				# Fill the plot at the specified bin
				self._fill_cutflow(1)
				# self.h['CutFlow'][self.ch].Fill(1)
			else:
				return
		
		if self._pv_cut(): #Check to make sure event has a PV otherwise throw event away (this happens very rarely with data).
			self._fill_cutflow(2)
			# self.h['CutFlow'][self.ch].Fill(2)
		else:
			return


		if self.do_invert_trigger_cut:
			if self._invert_trigger_cut():
				self._fill_cutflow(2)
				# self.h['CutFlow'][self.ch].Fill(2)
			else:
				return

		if self.do_filter_cut:
			if self._filter_cut():
				self._fill_cutflow(3)
				# self.h['CutFlow'][self.ch].Fill(3)
			else:
				return

		if self.do_prompt_lepton_cut:
			if self._prompt_lepton_cut():
				self._fill_cutflow(4)
				# self.h['CutFlow'][self.ch].Fill(4)
			else:
				return

		if self.do_invert_prompt_lepton_cut:
			if self._invert_prompt_lepton_cut():
				self._fill_cutflow(4)
				# self.h['CutFlow'][self.ch].Fill(4)
			else:
				return

		if self.do_ndv_cut:
			if self._ndv_cut():
				self._fill_cutflow(5)
				# self.h['CutFlow'][self.ch].Fill(5)
			else:
				return

		# If you've made it here, preselection is passed
		self.passed_preselection_cuts = True
		if not self.tree.is_data:
			self._fill_truth_histos(sel='truth/presel')

	def calculate_event_weight(self):
		######################################################################################################
		# MC re-weighting to include spin correlations
		######################################################################################################
		if not self.tree.is_data:
			official_samples = True
			self.MCEventType = selections.MCEventType(self.tree, wrong_lep_order=official_samples)
			# self.weight = self.tree.mass_lt_weight*self.MCEventType.weight  #if not weight_override else weight_override
			self.weight = self.tree.mass_lt_weight  # dont apply the weighting
		else:
			self.weight = self.tree.mass_lt_weight  # for data, mass_lt_weight = 1

		# if self.MCEventType.weight < 0:
		# 	print "--------"
		# 	print "isLNC ", self.MCEventType.isLNC
		# 	print "isLNV ", self.MCEventType.isLNV
		# 	print "MC weight: ", self.MCEventType.weight
		# 	print "M2 spin corr: ", self.MCEventType.M2_spin_corr
		# 	print "M2 no corr: ", self.MCEventType.M2_nocorr
		# 	print "s13: ", self.MCEventType.s13
		# 	print "s24: ", self.MCEventType.s24
		# 	print "p1 (px,py,pz,m): ", self.MCEventType.p_1.Px(),self.MCEventType.p_1.Py(),self.MCEventType.p_1.Pz(),self.MCEventType.p_1.M()
		# 	print "p2 (px,py,pz,m): ", self.MCEventType.p_2.Px(),self.MCEventType.p_2.Py(),self.MCEventType.p_2.Pz(),self.MCEventType.p_2.M()
		# 	print "p3 (px,py,pz,m): ", self.MCEventType.p_3.Px(),self.MCEventType.p_3.Py(),self.MCEventType.p_3.Pz(),self.MCEventType.p_3.M()
		# 	print "p4 (px,py,pz,m): ", self.MCEventType.p_4.Px(),self.MCEventType.p_4.Py(),self.MCEventType.p_4.Pz(),self.MCEventType.p_4.M()

	def DVSelection(self):
		raise NotImplementedError("Please implement this method in your own Analysis subclass")

	def _fill_cutflow(self, nbin):
		if not self.tree.is_data:
			if self.MCEventType.isLNC:
				self.CutFlow_LNC.Fill(nbin)
			if self.MCEventType.isLNV:
				self.CutFlow_LNV.Fill(nbin)
			self.CutFlow.Fill(nbin)
		else:
			self.CutFlow.Fill(nbin)

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
		
		for i in xrange(self.tree.ntrk): 
			#if self.tree.dv('trk_electronIndex')[i] >= 0: 
			# if self.tree.dv('trk_muonIndex')[i] >= 0:
			self.fill_hist(sel, 'DV_trk_{}_pt'.format(i), self.tree.dv('trk_pt_wrtSV')[i])
			self.fill_hist(sel, 'DV_trk_{}_eta'.format(i), self.tree.dv('trk_eta_wrtSV')[i])
			self.fill_hist(sel, 'DV_trk_{}_phi'.format(i), self.tree.dv('trk_phi_wrtSV')[i])
			self.fill_hist(sel, 'DV_trk_{}_d0'.format(i), self.tree.dv('trk_d0')[i])
			self.fill_hist(sel, 'DV_trk_{}_z0'.format(i), self.tree.dv('trk_z0')[i])
			self.fill_hist(sel, 'DV_trk_{}_charge'.format(i), self.tree.dv('trk_charge')[i])
			self.fill_hist(sel, 'DV_trk_{}_chi2'.format(i), self.tree.dv('trk_chi2_toSV')[i])
				

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
			self.micro_ntuples[sel].fill()

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
		self.fill_hist(sel, 'HNL_pt', truth_info.HNL_vec.Pt())
		self.fill_hist(sel, 'HNL_eta', truth_info.HNL_vec.Eta())
		self.fill_hist(sel, 'HNL_phi', truth_info.HNL_vec.Phi())
		self.fill_hist(sel, 'HNL_mass', truth_info.HNL_vec.M())

		self.fill_hist(sel, 'mHNLcalc', truth_info.mhnl)

		self.fill_hist(sel, 'DV_r', truth_info.truth_dvr)
		
		self.fill_hist(sel, 'DV_x', truth_info.truth_dvx)
		self.fill_hist(sel, 'DV_y', truth_info.truth_dvy)
		self.fill_hist(sel, 'DV_z', truth_info.truth_dvz)
		self.fill_hist(sel, 'plep_pt', truth_info.plep_vec.Pt())
		self.fill_hist(sel, 'plep_eta', truth_info.plep_vec.Eta())
		self.fill_hist(sel, 'plep_phi', truth_info.plep_vec.Phi())
		self.fill_hist(sel, 'plep_mass', truth_info.plep_vec.M())

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
			self.fill_hist(sel, 'nu_trk_pt', truth_info.dNu_vec.Pt())
			self.fill_hist(sel, 'nu_trk_eta', truth_info.dNu_vec.Eta())
			self.fill_hist(sel, 'nu_trk_phi', truth_info.dNu_vec.Phi())

			for itrk in range(2):
				self.fill_hist(sel, 'DV_trk_pt', truth_info.trkVec[itrk].Pt(), fill_ntuple=False)
				self.fill_hist(sel, 'DV_trk_eta', truth_info.trkVec[itrk].Eta(), fill_ntuple=False)
				self.fill_hist(sel, 'DV_trk_phi', truth_info.trkVec[itrk].Phi(), fill_ntuple=False)
			# TODO: figure out a ntuple scheme that can store these variables as well
		if sel == self.saveNtuples or self.saveNtuples == 'allcuts': 
			self.micro_ntuples[sel].fill()


	def _fill_selected_dv_histos(self, sel, do_lock=True):

		if self._locked < FILL_LOCKED and do_lock:
			# these are the histograms you only want to fill ONCE per DV
			# sel refers to the last selection that was applied

			# fill event weight. storing this per dv as weights include dv scale factor.
			self.fill_hist(sel, 'DV_weight', self.weight)

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

				tracks = helpers.Tracks(self.tree)
				tracks.getTracks()
				trkVec = tracks.lepVec

				muons = helpers.Tracks(self.tree)
				muons.getMuons()
				muVec = muons.lepVec

				electrons = helpers.Tracks(self.tree)
				electrons.getElectrons()
				elVec = electrons.lepVec

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

					deta = abs(tracks.eta[0] - tracks.eta[1])
					dphi = abs(tracks.lepVec[0].DeltaPhi(tracks.lepVec[1]))
					dpt = abs(tracks.pt[0] - tracks.pt[1])
					dR = tracks.lepVec[0].DeltaR(tracks.lepVec[1])

					if dR == 0.0:
						self.fill_hist(sel, 'DV_redmass', -1)
						self.fill_hist(sel, 'DV_redmassvis', -1)
						self.fill_hist(sel, 'DV_redmassHNL', -1)
					else:
						self.fill_hist(sel, 'DV_redmass', self.tree.dv('mass')/dR)
						self.fill_hist(sel, 'DV_redmassvis', Mlll.mlll/dR)
						self.fill_hist(sel, 'DV_redmassHNL', Mhnl.mhnl/dR)

					self.fill_hist(sel, 'DV_trk_deta', deta)
					self.fill_hist(sel, 'DV_trk_dphi', dphi)
					self.fill_hist(sel, 'DV_trk_dpt', dpt)
					self.fill_hist(sel, 'DV_trk_dR', dR)

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

			trk_quality = selections.Trackqual(self.tree)
			self.fill_hist(sel, 'DV_2tight', trk_quality.DV_2tight)
			self.fill_hist(sel, 'DV_2medium', trk_quality.DV_2medium)
			self.fill_hist(sel, 'DV_2loose', trk_quality.DV_2loose)
			self.fill_hist(sel, 'DV_1tight', trk_quality.DV_1tight)
			self.fill_hist(sel, 'DV_1medium', trk_quality.DV_1medium)
			self.fill_hist(sel, 'DV_1loose', trk_quality.DV_1loose)

			# better to fill truth matched DVs... need to fix this -DT
			# if not self.tree.is_data:
			# 	truthInfo = helpers.Truth()
			# 	truthInfo.getTruthParticles(self.tree)
			# 	self.fill_hist('truth_'+sel, 'W_pt', truthInfo.W_vec.Pt())
			# 	self.fill_hist('truth_'+sel, 'W_eta', truthInfo.W_vec.Eta())
			# 	self.fill_hist('truth_'+sel, 'W_phi', truthInfo.W_vec.Phi())
			# 	self.fill_hist('truth_'+sel, 'W_mass', truthInfo.W_vec.M())
			# 	self.fill_hist('truth_'+sel, 'HNL_pt', truthInfo.HNL_vec.Pt())
			# 	self.fill_hist('truth_'+sel, 'HNL_eta', truthInfo.HNL_vec.Eta())
			# 	self.fill_hist('truth_'+sel, 'HNL_phi', truthInfo.HNL_vec.Phi())
			# 	self.fill_hist('truth_'+sel, 'HNL_mass', truthInfo.HNL_vec.M())

			# 	self.fill_hist('truth_'+sel, 'mHNLcalc', truthInfo.mhnl)

			# 	self.fill_hist('truth_'+sel, 'DV_r', truthInfo.truth_dvr)
			# 	self.fill_hist('truth_'+sel, 'DV_x', truthInfo.truth_dvx)
			# 	self.fill_hist('truth_'+sel, 'DV_y', truthInfo.truth_dvy)
			# 	self.fill_hist('truth_'+sel, 'DV_z', truthInfo.truth_dvz)
			# 	self.fill_hist('truth_'+sel, 'plep_pt', truthInfo.plep_vec.Pt())
			# 	self.fill_hist('truth_'+sel, 'plep_eta', truthInfo.plep_vec.Eta())
			# 	self.fill_hist('truth_'+sel, 'plep_phi', truthInfo.plep_vec.Phi())
			# 	self.fill_hist('truth_'+sel, 'plep_mass', truthInfo.plep_vec.M())

			# 	self.fill_hist('truth_'+sel, 'maxlinkTruth_score', self.tree.dv('maxlinkTruth_score') )
			# 	self.fill_hist('truth_'+sel, 'maxlinkTruth_parent_pdgId', self.tree.dv('maxlinkTruth_parent_pdgId') )


			# 	if len(truthInfo.trkVec) == 2: 

			# 		charge_1 = truthInfo.plep_charge # charge of prompt lepton
			# 		p_1 = truthInfo.plep_vec # prompt lepton 

			# 		# if self.get_LNC: 
			# 		# 	if charge_1 != truthInfo.dLepCharge[0]: 
			# 		# 		p_2 = truthInfo.dLepVec[0]
			# 		# 		p_3 = truthInfo.dLepVec[1]
			# 		# 		p_4 = truthInfo.dLepVec[2]
			# 		# 	else: 
			# 		# 		p_2 = truthInfo.dLepVec[1]
			# 		# 		p_3 = truthInfo.dLepVec[0]
			# 		# 		p_4 = truthInfo.dLepVec[2]
			# 		# if self.get_LNV: 
			# 		# 	if charge_1 == truthInfo.dLepCharge[0]: 
			# 		# 		p_2 = truthInfo.dLepVec[0]
			# 		# 		p_3 = truthInfo.dLepVec[1]
			# 		# 		p_4 = truthInfo.dLepVec[2]
			# 		# 	else: 
			# 		# 		p_2 = truthInfo.dLepVec[1]
			# 		# 		p_3 = truthInfo.dLepVec[0]
			# 		# 		p_4 = truthInfo.dLepVec[2]


			# 		self.fill_hist('truth_'+sel, 'lep1_trk_pt', truthInfo.dLepVec[0].Pt())
			# 		self.fill_hist('truth_'+sel, 'lep1_trk_eta', truthInfo.dLepVec[0].Eta())
			# 		self.fill_hist('truth_'+sel, 'lep1_trk_phi', truthInfo.dLepVec[0].Phi())

			# 		self.fill_hist('truth_'+sel, 'lep2_trk_pt', truthInfo.dLepVec[1].Pt())
			# 		self.fill_hist('truth_'+sel, 'lep2_trk_eta', truthInfo.dLepVec[1].Eta())
			# 		self.fill_hist('truth_'+sel, 'lep2_trk_phi', truthInfo.dLepVec[1].Phi())
			# 		for itrk in range(2):
			# 			self.fill_hist('truth_'+sel, 'DV_trk_pt', truthInfo.trkVec[itrk].Pt(), fill_ntuple=False)
			# 			self.fill_hist('truth_'+sel, 'DV_trk_eta', truthInfo.trkVec[itrk].Eta(), fill_ntuple=False)
			# 			self.fill_hist('truth_'+sel, 'DV_trk_phi', truthInfo.trkVec[itrk].Phi(), fill_ntuple=False)
				
			# 	#BUG here to add truth micro ntuples
			# 	if sel == self.saveNtuples or self.saveNtuples == 'allcuts': 
			# 		self.micro_ntuples['truth_'+sel].fill()

			# fill TTree with ntuple information. Already set by fill_hist
			if sel == self.saveNtuples or self.saveNtuples == 'allcuts':  
				self.micro_ntuples[sel].fill()

			if sel == "sel":
				self._locked = FILL_LOCKED  # this only becomes unlocked after the event loop finishes in makeHistograms so you can only fill one DV from each event.


class oldAnalysis(Analysis):
	def __init__(self, name, tree, vtx_container, selections, outputFile, saveNtuples, debug_level):
		
		Analysis.__init__(self, name, tree, vtx_container, selections, outputFile, saveNtuples, debug_level)
		self.logger.info('Running  Old Analysis cuts')

		# Define cutflow histogram "by hand"
		# TODO Maybe the directory here needs to change, or the selection needs to be set
		self.observables.histogram_dict['CutFlow'] = ROOT.TH1D('CutFlow', 'CutFlow', 15, -0.5, 14.5)
		self.CutFlow = self.observables.histogram_dict['CutFlow']
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
			self.CutFlow.GetXaxis().SetBinLabel(5, "tight prompt %s" % self.plep)
		if self.do_invert_prompt_lepton_cut:
			self.CutFlow.GetXaxis().SetBinLabel(5, "invert prompt lepton")
		if self.do_ndv_cut:
			self.CutFlow.GetXaxis().SetBinLabel(6, "DV")
		if self.do_fidvol_cut:
			self.CutFlow.GetXaxis().SetBinLabel(7, "fiducial")
		if self.do_ntrk_cut:
			self.CutFlow.GetXaxis().SetBinLabel(8, "%s-track DV" % self.ntrk)
		if self.do_opposite_sign_cut:
			self.CutFlow.GetXaxis().SetBinLabel(9, "OS DV")
		if self.do_same_sign_cut:
			self.CutFlow.GetXaxis().SetBinLabel(9, "SS DV")
		if self.do_dv_type_cut:
			self.CutFlow.GetXaxis().SetBinLabel(10, "%s DV" % self.dv_type)
		if self.do_track_quality_cut:
			self.CutFlow.GetXaxis().SetBinLabel(11, "{}-lepton DV".format(self.track_quality))
		if self.do_cosmic_veto_cut:
			self.CutFlow.GetXaxis().SetBinLabel(12, "cosmic veto")
		if self.do_trilepton_mass_cut:
			self.CutFlow.GetXaxis().SetBinLabel(13, "m_{lll}")
		if self.do_dv_mass_cut:
			self.CutFlow.GetXaxis().SetBinLabel(14, "m_{DV}")
		self.CutFlow.GetXaxis().SetBinLabel(15, "truth matched")

		# Store LNC and LNV cutflows in the observables collection
		self.CutFlow_LNV = self.CutFlow.Clone()
		self.CutFlow_LNC = self.CutFlow.Clone()
		self.CutFlow_LNV.SetName("CutFlow_LNV"+"_"+self.ch)
		self.CutFlow_LNC.SetName("CutFlow_LNC"+"_"+self.ch)
		self.observables.histogram_dict['CutFlow_LNV'] = self.CutFlow_LNV
		self.observables.histogram_dict['CutFlow_LNC'] = self.CutFlow_LNC

	
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
					self._fill_cutflow(6)
					# self.h['CutFlow'][self.ch].Fill(6)
					self.passed_fidvol_cut = True
			# If this cut doesn't pass, don't continue to check other cuts
			else:
				return

		if self.do_ntrk_cut:
			if self._ntrk_cut():
				if not self.passed_ntrk_cut:
					self._fill_cutflow(7)
					# self.h['CutFlow'][self.ch].Fill(7)
					self.passed_ntrk_cut = True
			else:
				return

		if self.do_opposite_sign_cut or self.do_same_sign_cut:
			if self._charge_cut():
				if not self.passed_charge_cut:
					self._fill_cutflow(8)
					# self.h['CutFlow'][self.ch].Fill(8)
					self.passed_charge_cut = True
			else:
				return

		# self._fill_selected_dv_histos("charge")

		if self.do_dv_type_cut:
			if self._dv_type_cut():
				if not self.passed_dv_type_cut:
					self._fill_cutflow(9)
					# self.h['CutFlow'][self.ch].Fill(9)
					self.passed_dv_type_cut = True
			else:
				return
	
		self._fill_selected_dv_histos("DVtype")

		if self.do_track_quality_cut:
			if self._track_quality_cut():
				if not self.passed_track_quality_cut:
					self._fill_cutflow(10)
					# self.h['CutFlow'][self.ch].Fill(10)
					self.passed_track_quality_cut = True
			else:
				return
		
		self._fill_selected_dv_histos("trkqual")

		if self.do_cosmic_veto_cut:
			if self._cosmic_veto_cut():
				if not self.passed_cosmic_veto_cut:
					self._fill_cutflow(11)
					# self.h['CutFlow'][self.ch].Fill(11)
					self.passed_cosmic_veto_cut = True
			else:
				return

		self._fill_selected_dv_histos("cosmic")


		if self.do_trilepton_mass_cut:
			if self._trilepton_mass_cut():
				if not self.passed_trilepton_mass_cut:
					self._fill_cutflow(12)
					# self.h['CutFlow'][self.ch].Fill(12)
					self.passed_trilepton_mass_cut = True
			else:
				return
	
		self._fill_selected_dv_histos("mlll")


		if self.do_dv_mass_cut:
			if self._dv_mass_cut():
				if not self.passed_dv_mass_cut:
					self._fill_cutflow(13)
					# self.h['CutFlow'][self.ch].Fill(13)
					self.passed_dv_mass_cut = True
			else:
				return

		# Fill histos of truth-matched DVs
		if not self.tree.is_data:
			if self._truth_match():
				self._fill_cutflow(14)
				# self.h['CutFlow'][self.ch].Fill(14)
				self._fill_selected_dv_histos("match")

		# Fill all the histograms with only selected DVs. (ie. the ones that pass the full selection)
		self._fill_selected_dv_histos("sel")


class ToyAnalysis(Analysis):
	def __init__(self, name, tree, vtx_container, selections, outputFile, saveNtuples, debug_level):
		Analysis.__init__(self, name, tree, vtx_container, selections, outputFile, saveNtuples, debug_level)
		self.logger.info('Running  Toy Analysis cuts')

		self.add('CutFlow_all', 14, -0.5, 13.5)
		# Bin labels are 1 greater than histogram bins
		self.CutFlow.GetXaxis().SetBinLabel(1, "all")
		if self.do_trigger_cut:
			if not self.do_CR:
				self.CutFlow.GetXaxis().SetBinLabel(2, "trigger")
			else:
				self.CutFlow.GetXaxis().SetBinLabel(2, "DAOD_RPVLL triggers")
		self.CutFlow.GetXaxis().SetBinLabel(3, "PV")
		if self.do_filter_cut:
			self.CutFlow.GetXaxis().SetBinLabel(4, "%s" % self.filter_type)
		if self.do_prompt_lepton_cut:
			self.CutFlow.GetXaxis().SetBinLabel(5, "tight prompt %s" % self.plep)
		if self.do_invert_prompt_lepton_cut:
			self.CutFlow.GetXaxis().SetBinLabel(5, "invert prompt lepton")
		if self.do_ndv_cut:
			self.CutFlow.GetXaxis().SetBinLabel(6, "DV")
		if self.do_fidvol_cut:
			self.CutFlow.GetXaxis().SetBinLabel(7, "fiducial")
		if self.do_ntrk_cut:
			self.CutFlow.GetXaxis().SetBinLabel(8, "%s-track DV" % self.ntrk)
		# if self.do_HNL_mass_cut:
		self.CutFlow.GetXaxis().SetBinLabel(9, "dR")
		if self.do_opposite_sign_cut:
			self.CutFlow.GetXaxis().SetBinLabel(10, "OS DV")
		if self.do_same_sign_cut:
			self.CutFlow.GetXaxis().SetBinLabel(10, "SS DV")
		# self.CutFlow.GetXaxis().SetBinLabel(11, "++ DV")
		# self.CutFlow.GetXaxis().SetBinLabel(12, "-- DV")
		# self.CutFlow.GetXaxis().SetBinLabel(13, "+++ lll")
		# self.CutFlow.GetXaxis().SetBinLabel(14, "+-- lll")
		# self.CutFlow.GetXaxis().SetBinLabel(15, "-++ lll")
		# self.CutFlow.GetXaxis().SetBinLabel(16, "--- lll")
		if self.do_dv_type_cut:
			self.CutFlow.GetXaxis().SetBinLabel(11, "%s DV" % self.dv_type)
		# self.CutFlow.GetXaxis().SetBinLabel(18, "++ DV")
		# self.CutFlow.GetXaxis().SetBinLabel(19, "-- DV")
		# self.CutFlow.GetXaxis().SetBinLabel(20, "+++ lll")
		# self.CutFlow.GetXaxis().SetBinLabel(21, "+-- lll")
		# self.CutFlow.GetXaxis().SetBinLabel(22, "-++ lll")
		# self.CutFlow.GetXaxis().SetBinLabel(23, "--- lll")
		if self.do_dv_mass_cut:
			self.CutFlow.GetXaxis().SetBinLabel(12, "m_{DV}")
		if self.do_trilepton_mass_cut:
			self.CutFlow.GetXaxis().SetBinLabel(13, "m_{lll}")
		if self.do_cosmic_veto_cut:
			self.CutFlow.GetXaxis().SetBinLabel(14, "cosmic veto")

		self.h['CutFlow_LNV'] = {}
		self.h['CutFlow_LNC'] = {}
		self.CutFlow_LNV = self.CutFlow.Clone()
		self.CutFlow_LNC = self.CutFlow.Clone()

	def _fill_correlation_histos(self, sel):
		w = self.weight
		# sel refers to the last selection that was applied

		if self.do_prompt_lepton_cut:

			tracks = helpers.Tracks(self.tree)
			tracks.getTracks()
			trkVec = tracks.lepVec

			muons = helpers.Tracks(self.tree)
			muons.getMuons()
			muVec = muons.lepVec

			electrons = helpers.Tracks(self.tree)
			electrons.getElectrons()
			elVec = electrons.lepVec

			if tracks.ntracks == 2:
				plep_vec = self.plep_sel.plepVec
				Mlll = selections.Mlll(dv_type=self.dv_type, plep=plep_vec, dMu=muVec, dEl=elVec)
				Mhnl = selections.Mhnl(self.tree, self.dv_type, plep=plep_vec, dMu=muVec,dEl=elVec)

			# if tracks.ntracks == 2:
			# 	Mltt = selections.Mltt(plep=self.plep_sel.plepVec, trks=trkVec)
			# 	Mhnl = selections.Mhnl(self.tree, plep=self.plep_sel.plepVec, trks=trkVec)
			# 	Mtrans = selections.Mtrans(plep=self.plep_sel.plepVec, trks=trkVec)

				# fill 2D mass correlation plots here
				self.fill_hist(sel, 'DVmass_mvis', self.tree.dv('mass'), Mlll.mlll)
				self.fill_hist(sel, 'DVmass_mhnl', self.tree.dv('mass'), Mhnl.mhnl)
				self.fill_hist(sel, 'DVmass_mtrans', self.tree.dv('mass'), Mlll.mtrans)
				self.fill_hist(sel, 'DVmass_hnlpt', self.tree.dv('mass'), Mhnl.hnlpt)
				self.fill_hist(sel, 'mvis_mhnl', Mlll.mlll, Mhnl.mhnl)
				self.fill_hist(sel, 'mvis_mtrans', Mlll.mlll, Mlll.mtrans)
				self.fill_hist(sel, 'mvis_hnlpt', Mlll.mlll, Mhnl.hnlpt)
				self.fill_hist(sel, 'mhnl_mtrans', Mhnl.mhnl, Mlll.mtrans)
				self.fill_hist(sel, 'mhnl_hnlpt', Mhnl.mhnl, Mhnl.hnlpt)
				self.fill_hist(sel, 'mhnl2D', Mhnl.mhnl, Mhnl.alt_mhnl)

	#########################################################################################################################
	# Define new cuts you want to apply here. This will overwrite whatever cuts are defined in the parent analysis class.
	#########################################################################################################################
	def _track_quality_cut_1tight(self):
		track_quality_sel = selections.Trackqual(self.tree, quality="1-tight")
		return track_quality_sel.passes()

	def _track_quality_cut_2tight(self):
		track_quality_sel = selections.Trackqual(self.tree, quality="2-tight")
		return track_quality_sel.passes()

	def _dR_cut(self):
		tracks = helpers.Tracks(self.tree)
		tracks.getTracks()
		trkVec = tracks.lepVec
		dR = tracks.lepVec[0].DeltaR(tracks.lepVec[1])
		return dR > 0.0

	def _trilepton_mass_cut(self):
		plep_vec = self.plep_sel.plepVec

		muons = helpers.Tracks(self.tree)
		muons.getMuons()
		muVec = muons.lepVec

		electrons = helpers.Tracks(self.tree)
		electrons.getElectrons()
		elVec = electrons.lepVec

		mlll_sel = selections.Mlll(dv_type=self.dv_type, plep=plep_vec, dMu=muVec, dEl=elVec)
		return mlll_sel.passes()

	def _dv_mass_cut(self):
		dv_mass_sel = selections.DVmass(self.tree, dvmasscut=2)  # changed the dvmass cut to 2 GeV
		return dv_mass_sel.passes()

	def _multitrk_2lep_cut(self):
		if self.tree.dv('ntrk') >= 2:  # 2+ trk vertex
			dv_type_sel = selections.DVtype(self.tree, dv_type=self.dv_type)
			if dv_type_sel.passes():  # 2 leptons in the DV
				sign_pair = "SS" if self.do_same_sign_cut else "OS"
				charge_sel = selections.ChargeDV(self.tree, sel=sign_pair, trk_charge=dv_type_sel.lepton_charge)
				return charge_sel.passes()

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


	def _HNL_mass_cut(self):  # not used
		tracks = helpers.Tracks(self.tree)
		tracks.getTracks()

		plep_vec = self.plep_sel.plepVec
		tracks_vec = tracks.lepVec

		mHNL_sel = selections.Mhnl(self.tree, plep=plep_vec, trks=tracks_vec, hnlmasscut=3)

		return mHNL_sel.passes()


	def _HNL_pt_cut(self):
		tracks = helpers.Tracks(self.tree)
		tracks.getTracks()

		plep_vec = self.plep_sel.plepVec
		tracks_vec = tracks.lepVec

		mHNL_sel = selections.Mhnl(plep=plep_vec,trks=tracks_vec)

		return (mHNL_sel.hnlpt > 20 and mHNL_sel.hnlpt < 60)

	def initialize_cut_bools(self):
		###########################################################################################################################
		# Initialize the cut bools every event. These bools tell the code if the cutflow has already been filled for this event.
		# Default is to select the first event that passes the selection
		###########################################################################################################################
		self.passed_preselection_cuts = False
		self.passed_fidvol_cut = False
		self.passed_ntrk_cut = False
		self.passed_charge_cut = False
		self.passed_dv_type_cut = False
		self.passed_dR_cut = False
		# self.passed_track_quality_cut = False
		self.passed_track_2tight_cut = False
		self.passed_track_1tight_cut = False
		self.passed_cosmic_veto_cut = False
		self.passed_trilepton_mass_cut = False
		self.passed_dv_mass_cut = False
		self.passed_HNL_mass_cut = False
		self.passed_HNL_pt_cut = False


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
					self._fill_cutflow(6)
					# self.h['CutFlow'][self.ch].Fill(6)
					self.passed_fidvol_cut = True
			# If this cut doesn't pass, don't continue to check other cuts
			else:
				return

		self.fill_hist('presel', 'num_trks', self.tree.dv('ntrk')) # fill ntrk cut after pre-selection but before ntrk cut
		
		if self._multitrk_2lep_cut(): # no return becuase this is not an analysis cut, only used for studying S & B 
			self._fill_multitrk_histos()

		if self.do_ntrk_cut:
			if self._ntrk_cut():
				if not self.passed_ntrk_cut:
					self._fill_cutflow(7)
					# self.h['CutFlow'][self.ch].Fill(7)
					self.passed_ntrk_cut = True
			else:
				return

		if self._dR_cut():
			if not self.passed_dR_cut:
				self._fill_cutflow(8)
				# self.h['CutFlow'][self.ch].Fill(8)
				self.passed_dR_cut = True
		else:
			return

		if self.do_opposite_sign_cut or self.do_same_sign_cut:
			if self._charge_cut():
				if not self.passed_charge_cut:
					self._fill_cutflow(9)
					# self.h['CutFlow'][self.ch].Fill(9)
					self.passed_charge_cut = True
			else:
				return

		if self.do_dv_type_cut:
			if self._dv_type_cut():
				if not self.passed_dv_type_cut:
					self._fill_cutflow(10)
					# self.h['CutFlow'][self.ch].Fill(10)
					self.passed_dv_type_cut = True
			else:
				return

		if not self.tree.is_data:
			if self.MCEventType.isLNC: 
				self._fill_selected_dv_histos("DVtype_LNC")
				self._fill_correlation_histos("DVtype_LNC")
			if self.MCEventType.isLNV:
				self._fill_selected_dv_histos("DVtype_LNV")
				self._fill_correlation_histos("DVtype_LNV")
		else:
			self._fill_selected_dv_histos("DVtype")
			self._fill_correlation_histos("DVtype")


		if self.do_dv_mass_cut:
			if self._dv_mass_cut():
				if not self.passed_dv_mass_cut:
					self._fill_cutflow(11)
					# self.h['CutFlow'][self.ch].Fill(11)
					self.passed_dv_mass_cut = True
			else:
				return

		if not self.tree.is_data:
			if self.MCEventType.isLNC: 
				self._fill_selected_dv_histos("mDV_LNC")
				self._fill_correlation_histos("mDV_LNC")
			if self.MCEventType.isLNV:
				self._fill_selected_dv_histos("mDV_LNV")
				self._fill_correlation_histos("mDV_LNV")
		else:
			self._fill_selected_dv_histos("mDV")
			self._fill_correlation_histos("mDV")

		if self.do_trilepton_mass_cut:
			if self._trilepton_mass_cut():
				if not self.passed_trilepton_mass_cut:
					self._fill_cutflow(12)
					# self.h['CutFlow'][self.ch].Fill(12)
					self.passed_trilepton_mass_cut = True
			else:
				return

		if not self.tree.is_data:
			if self.MCEventType.isLNC: 
				self._fill_selected_dv_histos("mlll_LNC")
				self._fill_correlation_histos("mlll_LNC")
			if self.MCEventType.isLNV:
				self._fill_selected_dv_histos("mlll_LNV")
				self._fill_correlation_histos("mlll_LNV")
		else:
			self._fill_selected_dv_histos("mlll")
			self._fill_correlation_histos("mlll")


		if self.do_cosmic_veto_cut:
			if self._cosmic_veto_cut():
				if not self.passed_cosmic_veto_cut:
					self._fill_cutflow(13)
					# self.h['CutFlow'][self.ch].Fill(14)
					self.passed_cosmic_veto_cut = True
			else:
				return
		
		if not self.tree.is_data:
			if self.MCEventType.isLNC: 
				self._fill_selected_dv_histos("cosmic_LNC")
			if self.MCEventType.isLNV:
				self._fill_selected_dv_histos("cosmic_LNV")

		else:
			self._fill_selected_dv_histos("cosmic")
			


		# Fill histos of truth-matched DVs
		if not self.tree.is_data:
			if self._truth_match():
				self._fill_cutflow(14)
				# self.h['CutFlow'][self.ch].Fill(14)
				if self.MCEventType.isLNC: 
					self._fill_selected_dv_histos("match_LNC")
				if self.MCEventType.isLNV:
					self._fill_selected_dv_histos("match_LNV")

		# Fill all the histograms with only selected DVs. (ie. the ones that pass the full selection)
		if not self.tree.is_data:
			if self.MCEventType.isLNC: 
				self._fill_selected_dv_histos("sel_LNC")
			if self.MCEventType.isLNV:
				self._fill_selected_dv_histos("sel_LNV")
		else:
			self._fill_selected_dv_histos("sel")
			


class KShort(Analysis):
	def __init__(self, name, tree, vtx_container, selections, outputFile, saveNtuples, debug_level):
		Analysis.__init__(self, name, tree, vtx_container, selections, outputFile, saveNtuples, debug_level)
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
			elif self.tree['muon_isMedium'][imu] == 1: self.fill_hist(sel, 'muon_quality', 2)
			elif self.tree['muon_isLoose'][imu] == 1:  self.fill_hist(sel, 'muon_quality', 1)
			else: self.fill_hist(sel, 'muon_quality', 0)

		for iel in range(len(self.tree['el_pt'])):
			self.fill_hist(sel, 'el_pt', self.tree['el_pt'][iel])
			self.fill_hist(sel, 'el_eta', self.tree['el_eta'][iel])
			self.fill_hist(sel, 'el_phi', self.tree['el_phi'][iel])
			if self.tree['el_LHTight'][iel] == 1:  self.fill_hist(sel, 'el_quality', 3)
			elif self.tree['el_LHMedium'][iel] == 1: self.fill_hist(sel, 'el_quality', 2)
			elif self.tree['el_LHLoose'][iel] == 1:  self.fill_hist(sel, 'el_quality', 1)
			else: self.fill_hist(sel, 'el_quality', 0)


	def _fill_selected_dv_histos(self, sel, do_lock=True):
		if not self.tree.is_data:
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
				self.micro_ntuples[sel].fill()


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
