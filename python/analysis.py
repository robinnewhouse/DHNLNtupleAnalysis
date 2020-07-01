# from ROOT import*
import ROOT
import numpy as np
import os
import sys
import helpers
import selections
import observables
import logging
import ntuples

logger = helpers.getLogger('dHNLAnalysis.analysis', level=logging.INFO)

UNLOCKED = 0
SELECTION_LOCKED = 1
FILL_LOCKED = 2


class Analysis(object):
	def __init__(self, tree, vtx_container, selections, outputFile, saveNtuples):
		self.sel = selections
		self.outputFile = outputFile
		self.fi = ROOT.TFile.Open(outputFile, 'update')
		self.ch = vtx_container
		self.histSuffixes = [self.ch]
		self.h = {}
		self.micro_ntuples = {}
		self.tree = tree
		self.saveNtuples = saveNtuples
		self._locked = UNLOCKED

		self.observables = [observable.registered(self) for observable in observables.ObservableList if ((observable.only is None) or any(only in self.sel for only in observable.only))]
		for observable in self.observables:
			if 'hist' in observable.do:
				if self.tree.is_data and observable.need_truth:
					continue
				elif type(observable.binning) == tuple:
					self.add(observable.name, *observable.binning)
				else:
					self.addVar(observable.name, observable.binning)

		# setting all the relevant variables for the cuts based on the input selections
		# trigger cut
		if 'alltriggers' in self.sel:
			self.trigger = 'all'
			self.do_trigger_cut = True
		else:
			if 'CR' not in self.sel:
				logger.warn('You did not specify a trigger configuration for this channel. Skipping trigger selection.')
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
				logger.warn('You did not specify a filter configuration for this channel. Skipping filter selection.')
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
				logger.warn('You did not specify a prompt lepton for this channel. Skipping prompt lepton selection.')
			self.do_prompt_lepton_cut = False

		if 'CR' in self.sel:  # DO NOT CHANGE THESE CUTS OR YOU MIGHT UNBLIND DATA!!!
			self.do_CR = True
			self.do_trigger_cut = False  # do not apply trigger cut
			self.do_invert_trigger_cut = False  # do not apply inverted trigger cut
			self.do_filter_cut = False  # do not apply filter cut
			self.do_prompt_lepton_cut = False  # do not apply prompt lepton cut
			self.do_invert_prompt_lepton_cut = True  # invert prompt lepton cut
			logger.info('You are setup up to look in the inverted prompt lepton control region!')
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
		if not self.do_ndv_cut: logger.warn('You did not add nDV cut. Skipping nDV selection.')

		# fiducial volume
		self.do_fidvol_cut = 'fidvol' in self.sel
		if not self.do_fidvol_cut:
			logger.warn('You did not add DV in fiducial volume cut. Skipping DV in fiducial volume selection.')

		# 2 (or more) track cut
		self.do_ntrk_cut = True
		if '2track' in self.sel:
			self.ntrk = 2
		elif '3track' in self.sel:
			self.ntrk = 3
		elif '4track' in self.sel:
			self.ntrk = 4
		else:
			logger.warn('You did not add an ntrack cut. Skipping ntrack selection.')
			self.do_ntrk_cut = False

		# Opposite sign children vertex cut
		self.do_opposite_sign_cut = 'OS' in self.sel
		# Same sign children vertex cut
		self.do_same_sign_cut = 'SS' in self.sel
		if not (self.do_opposite_sign_cut or self.do_same_sign_cut):
			logger.warn('You did not add an SS or OS track cut. Skipping SS/OS track selection.')

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
			logger.warn('You did not specify a DV type for this channel. Skipping DV type selection.')
			self.do_dv_type_cut = False

		# Track quality
		self.do_track_quality_cut = True
		if '1-tight' in self.sel:
			self.track_quality = '1-tight'
		elif '2-tight' in self.sel:
			self.track_quality = '2-tight'
		else:
			if "CR" not in self.sel:
				logger.warn('You did not specify a DV track quality for this channel. Skipping DV track quality selection.')
			self.do_track_quality_cut = False

		# cosmic veto cut
		self.do_cosmic_veto_cut = 'cosmicveto' in self.sel
		if not self.do_cosmic_veto_cut and 'CR' not in self.sel:
			logger.warn('You did not add a cosmic veto cut for this channel. Skipping cosmic veto selection.')

		# tri-lepton mass cut
		self.do_trilepton_mass_cut = 'mlll' in self.sel
		if not self.do_trilepton_mass_cut  and "CR" not in self.sel:
			logger.warn('You did not add a mlll cut for this channel. Skipping tri-lepton mass selection.')

		# DV mass cut
		self.do_dv_mass_cut = 'DVmass' in self.sel
		if not self.do_dv_mass_cut and "CR" not in self.sel:
			logger.warn('You did not add a DVmass cut for this channel. Skipping displaced vertex mass selection.')

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
		:return:
		"""
		if selection: full_name = selection + '_' + hist_name
		else: full_name = hist_name
		try:
			if variable_2 is None:
				self.h[full_name][self.ch].Fill(variable_1, self.tree.weight)
			else:
				self.h[full_name][self.ch].Fill(variable_1, variable_2, self.tree.weight)
		except KeyError as e:
			logger.error("Histogram {} not registered. Automatically adding with default binning.".format(full_name))
			observable = observables.Observable(full_name)
			observable.queue()
			self.add(observable.name, *observable.binning)
			self.fill_hist(selection, hist_name, variable_1, variable_2=variable_2, fill_ntuple=fill_ntuple)

		# Unless suppressed, fill the corresponding micro-ntuple with the variable
		# Will not fill variables from 2D histograms to prevent double-counting
		save_sel = self.saveNtuples == selection or 'truth_'+self.saveNtuples == selection or self.saveNtuples == 'allcuts'
		if fill_ntuple and variable_2 is None and save_sel:
			# Note: selection and hist_name will be overridden by full_name
			# Need selection to define ntuple tree
			self.fill_ntuple(selection, hist_name, variable_1)

	def add(self, hName, nBins, xLow, xHigh):
		self.h[hName] = {}
		self.h[hName][self.ch] = ROOT.TH1D(hName+"_"+self.ch, "", nBins, xLow, xHigh)
		self.h[hName][self.ch].Sumw2()
		self.h[hName][self.ch].SetDirectory(0)

	def add2D(self, hName, nBins, xLow, xHigh, nBinsY, yLow, yHigh):
		self.h[hName] = {}
		self.h[hName][self.ch] = ROOT.TH2D(hName + "_" + self.ch, "", nBins, xLow, xHigh, nBinsY, yLow, yHigh)
		self.h[hName][self.ch].Sumw2()
		self.h[hName][self.ch].SetDirectory(0)

	def fill_ntuple(self, selection, ntuple_name, variable, full_name=""):
		"""
		A helper function for filling micro-ntuples. Often called from the fill_hist function.
		If you are using this in you analysis,
		please check that it is not also being called by fill_hist to prevent double-counting.
		:param selection: the step of selection the analysis it at. May be "None" in which case there will be no prefix.
		:param ntuple_name: base name of the ntuple. When saved, a prefix and suffix will be appended.
		:param variable: variable you want to fill the histogram with.
		:param full_name: override the automatic naming of the ntuple.
		:return:
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
				logger.error("You are looking in the CR without prompt leptons so you cannot cut on mll, HNLpt or HNLm!!")
				sys.exit(1)  # abort because of error

			if not self.do_prompt_lepton_cut:
				logger.error("You cannot calculate mlll, HNLpt or HNLm without selecting a prompt lepton!")
				sys.exit(1)  # abort because of error

		if self.do_opposite_sign_cut and self.do_same_sign_cut:
			logger.error("These cuts are mutually exclusive. You will get zero events!")
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

		# Make a subdirectory for vertex type. May be other channels in the future.
		if not self.fi.FindObject(self.ch):
			self.fi.mkdir(self.ch, "Analysis Channel " + self.ch)
		# Move ROOT to the channel subdirectory
		self.fi.cd(self.ch)

		# Store saved histograms to file
		# TODO: this should be saved in a different way. e.g. another level of dictionary. Parsing strings for variable names is not good.
		for h_name in self.h:
			selection = h_name.split('_')[0]  # get selection
			sel_dir = self.ch + '/' + selection
			base_name = '_'.join(h_name.split('_')[1:])  # get base name
			if not self.fi.GetDirectory(sel_dir):  # make TDirectory if necessary
				self.fi.mkdir(sel_dir, "Analysis Selection " + selection)
			self.fi.cd(sel_dir)  # change to TDirectory
			self.h[h_name][self.ch].Write(base_name)  # save only the base name

		logger.info("Histograms written to {}".format(self.outputFile))

		self.fi.Close()

	def end(self):
		logger.info('Done with Channel("{}")'.format(self.ch))
		meta = []
		# if self.region:
		# 	meta.append('Region: {}'.format(self.region))
		# if self.period:
		# 	meta.append('Period: {}'.format(self.period))
		# if self.bcategory != None:
		# 	meta.append('B-tagging Categroy: {}'.format(self.bcategory))
		# logger.info('\t' + ' | '.join(meta))

		# gives warning messages if histograms are unfilled
		for histName in self.h:
			if self.h[histName][self.ch].GetEntries() == 0:
				logger.debug('\tUnfilled HIST({}<{}>)!'.format(histName, self.ch))

		# y_err = ROOT.Double()
		# for s in self.histSuffixes:
		# 	h = self.h['finalNEvents'][s]
		# 	logger.info('\tFinal NEvents <{syst}>: {y}'.format(syst = s, y = h.Integral()))
		# 	h = self.h['finalYields'][s]
		# 	logger.info('\tFinal Yields <{syst}>: {y:.4g}+/-{y_err:.2g}'.format(syst = s, y = h.IntegralAndError(0, h.GetNbinsX()+1, y_err), y_err = y_err))
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
		# 	logger.error(e, exc_info=True)


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

	def _prompt_lepton_cut(self, quality="tight", min_dR=0.05):
		self.plep_sel = selections.PromptLepton(self.tree, lepton=self.plep, quality=quality, min_dR=min_dR)
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
		return self.tree.dv('maxlinkTruth_score') > 0.75 and \
			   abs(self.tree.dv('maxlinkTruth_parent_pdgId')) == 50

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
		######################################################################################################
		# Preselection are all the cuts that are requied per event
		# Current cuts include: trigger, filter, plepton, DV cut
		######################################################################################################

		#initialize the cut bools for every event
		self.initialize_cut_bools()

		self._fill_leptons()

		if not self.tree.is_data:
			self._fill_truth_histos(sel='truth_all')

		self.h['CutFlow'][self.ch].SetBinContent(1, self.tree.cutflow[1])  # all events

		######################################################################################################
		# Selection code is deisgned so that it will pass the selection only if the cut true or cut is unused
		# ex. passTrigger is true if the trigcut is true OR if trigcut is not used)
		######################################################################################################

		if self._pv_cut(): #Check to make sure event has a PV otherwise throw event away (this happens very rarely with data).
			self.h['CutFlow'][self.ch].Fill(1)
		else:
			return

		if self.do_trigger_cut:
			if self._trigger_cut():
				# Fill the plot at the specified bin
				self.h['CutFlow'][self.ch].Fill(2)
			else:
				return

		if self.do_invert_trigger_cut:
			if self._invert_trigger_cut():
				self.h['CutFlow'][self.ch].Fill(2)
			else:
				return

		if self.do_filter_cut:
			if self._filter_cut():
				self.h['CutFlow'][self.ch].Fill(3)
			else:
				return

		if self.do_prompt_lepton_cut:
			if self._prompt_lepton_cut():
				self.h['CutFlow'][self.ch].Fill(4)
			else:
				return

		if self.do_invert_prompt_lepton_cut:
			if self._invert_prompt_lepton_cut():
				self.h['CutFlow'][self.ch].Fill(4)
			else:
				return

		if self.do_ndv_cut:
			if self._ndv_cut():
				self.h['CutFlow'][self.ch].Fill(5)
			else:
				return

		# If you've made it here, preselection is passed
		self.passed_preselection_cuts = True
		if not self.tree.is_data:
			self._fill_truth_histos(sel='truth_presel')

	# def preSelection(self):
	# 	raise NotImplementedError("Please implement this method in your own Analysis subclass")

	def DVSelection(self):
		raise NotImplementedError("Please implement this method in your own Analysis subclass")

	# Common histograms to fill
	def _fill_leptons(self):
		sel = 'all'
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

	def _fill_all_dv_histos(self):
		sel = 'all'
		# self.fill_hist(sel, 'charge_ntrk', self.tree.dv('charge'), self.tree.dv('ntrk'))
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
		
		if sel == self.saveNtuples or self.saveNtuples == 'allcuts': 
			self.micro_ntuples[sel].fill()

	def _fill_truth_histos(self, sel='truth_all'):
		truth_info = helpers.Truth()
		truth_info.getTruthParticles(self.tree)
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
			self.fill_hist(sel, 'lep1_trk_pt', truth_info.trkVec[0].Pt())
			self.fill_hist(sel, 'lep1_trk_eta', truth_info.trkVec[0].Eta())
			self.fill_hist(sel, 'lep1_trk_phi', truth_info.trkVec[0].Phi())
			self.fill_hist(sel, 'lep2_trk_pt', truth_info.trkVec[1].Pt())
			self.fill_hist(sel, 'lep2_trk_eta', truth_info.trkVec[1].Eta())
			self.fill_hist(sel, 'lep2_trk_phi', truth_info.trkVec[1].Phi())

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
			self.fill_hist(sel, 'DV_weight', self.tree.weight)

			# fill histograms that require a prompt lepton to be identified
			if self.do_prompt_lepton_cut:
				plep_vec = self.plep_sel.plepVec
				plepd0 = self.plep_sel.plepd0
				plepz0 = self.plep_sel.plepz0

				self.fill_hist(sel, 'plep_pt', plep_vec.Pt())
				self.fill_hist(sel, 'plep_eta', plep_vec.Eta())
				self.fill_hist(sel, 'plep_phi', plep_vec.Phi())
				self.fill_hist(sel, 'plep_d0', plepd0)
				self.fill_hist(sel, 'plep_z0', plepz0)

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

			# fill histograms for track values
			ntracks = self.tree.ntrk
			for itrk in range(ntracks):  # loop over tracks
				self.fill_hist(sel, 'DV_trk_pt', self.tree.dv('trk_pt_wrtSV')[itrk], fill_ntuple=False)
				self.fill_hist(sel, 'DV_trk_eta', self.tree.dv('trk_eta_wrtSV')[itrk], fill_ntuple=False)
				self.fill_hist(sel, 'DV_trk_phi', self.tree.dv('trk_phi_wrtSV')[itrk], fill_ntuple=False)
				self.fill_hist(sel, 'DV_trk_d0', self.tree.dv('trk_d0')[itrk], fill_ntuple=False)
				self.fill_hist(sel, 'DV_trk_z0', self.tree.dv('trk_z0')[itrk], fill_ntuple=False)
				self.fill_hist(sel, 'DV_trk_charge', self.tree.dv('trk_charge')[itrk], fill_ntuple=False)
				self.fill_hist(sel, 'DV_trk_chi2', self.tree.dv('trk_chi2_toSV')[itrk], fill_ntuple=False)

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

			# better to fill truth matched DVs... need to fix this -DT
			if not self.tree.is_data:
				truthInfo = helpers.Truth()
				truthInfo.getTruthParticles(self.tree)
				self.fill_hist('truth_'+sel, 'W_pt', truthInfo.W_vec.Pt())
				self.fill_hist('truth_'+sel, 'W_eta', truthInfo.W_vec.Eta())
				self.fill_hist('truth_'+sel, 'W_phi', truthInfo.W_vec.Phi())
				self.fill_hist('truth_'+sel, 'W_mass', truthInfo.W_vec.M())
				self.fill_hist('truth_'+sel, 'HNL_pt', truthInfo.HNL_vec.Pt())
				self.fill_hist('truth_'+sel, 'HNL_eta', truthInfo.HNL_vec.Eta())
				self.fill_hist('truth_'+sel, 'HNL_phi', truthInfo.HNL_vec.Phi())
				self.fill_hist('truth_'+sel, 'HNL_mass', truthInfo.HNL_vec.M())

				self.fill_hist('truth_'+sel, 'mHNLcalc', truthInfo.mhnl)

				self.fill_hist('truth_'+sel, 'DV_r', truthInfo.truth_dvr)
				self.fill_hist('truth_'+sel, 'DV_x', truthInfo.truth_dvx)
				self.fill_hist('truth_'+sel, 'DV_y', truthInfo.truth_dvy)
				self.fill_hist('truth_'+sel, 'DV_z', truthInfo.truth_dvz)
				self.fill_hist('truth_'+sel, 'plep_pt', truthInfo.plep_vec.Pt())
				self.fill_hist('truth_'+sel, 'plep_eta', truthInfo.plep_vec.Eta())
				self.fill_hist('truth_'+sel, 'plep_phi', truthInfo.plep_vec.Phi())
				self.fill_hist('truth_'+sel, 'plep_mass', truthInfo.plep_vec.M())
				if len(truthInfo.trkVec) == 2: 
					self.fill_hist('truth_'+sel, 'lep1_trk_pt', truthInfo.trkVec[0].Pt())
					self.fill_hist('truth_'+sel, 'lep1_trk_eta', truthInfo.trkVec[0].Eta())
					self.fill_hist('truth_'+sel, 'lep1_trk_phi', truthInfo.trkVec[0].Phi())

					self.fill_hist('truth_'+sel, 'lep2_trk_pt', truthInfo.trkVec[1].Pt())
					self.fill_hist('truth_'+sel, 'lep2_trk_eta', truthInfo.trkVec[1].Eta())
					self.fill_hist('truth_'+sel, 'lep2_trk_phi', truthInfo.trkVec[1].Phi())
					for itrk in range(2):
						self.fill_hist('truth_'+sel, 'DV_trk_pt', truthInfo.trkVec[itrk].Pt(), fill_ntuple=False)
						self.fill_hist('truth_'+sel, 'DV_trk_eta', truthInfo.trkVec[itrk].Eta(), fill_ntuple=False)
						self.fill_hist('truth_'+sel, 'DV_trk_phi', truthInfo.trkVec[itrk].Phi(), fill_ntuple=False)
				
				#BUG here to add truth micro ntuples
				if sel == self.saveNtuples or self.saveNtuples == 'allcuts': 
					self.micro_ntuples['truth_'+sel].fill()

			# fill TTree with ntuple information. Already set by fill_hist
			if sel == self.saveNtuples or self.saveNtuples == 'allcuts':  
				self.micro_ntuples[sel].fill()

			if sel == "sel":
				self._locked = FILL_LOCKED  # this only becomes unlocked after the event loop finishes in makeHistograms so you can only fill one DV from each event.


class oldAnalysis(Analysis):
	def __init__(self, tree, vtx_container, selections, outputFile, saveNtuples):
		logger.info('Running  Old Analysis cuts')
		Analysis.__init__(self, tree, vtx_container, selections, outputFile, saveNtuples)

		self.add2D('charge_ntrk', 11, -5.5, 5.5, 9, -0.5, 8.5)

		self.add('CutFlow', 14, -0.5, 13.5)
		# Bin labels are 1 greater than histogram bins
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(1, "all")
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(2, "PV")
		if self.do_trigger_cut:
			if self.do_CR == False:
				self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(3, "trigger")
			else:
				self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(3, "DAOD_RPVLL triggers")
		if self.do_invert_trigger_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(3, "invert trigger")
		if self.do_filter_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(4, "%s" % self.filter_type)
		if self.do_prompt_lepton_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(5, "tight prompt %s" % self.plep)
		if self.do_invert_prompt_lepton_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(5, "invert prompt lepton")
		if self.do_ndv_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(6, "DV")
		if self.do_fidvol_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(7, "fiducial")
		if self.do_ntrk_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(8, "%s-track DV" % self.ntrk)
		if self.do_opposite_sign_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(9, "OS DV")
		if self.do_same_sign_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(9, "SS DV")
		if self.do_dv_type_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(10, "%s DV" % self.dv_type)
		if self.do_track_quality_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(11, "{}-lepton DV".format(self.track_quality))
		if self.do_cosmic_veto_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(12, "cosmic veto")
		if self.do_trilepton_mass_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(13, "m_{lll}")
		if self.do_dv_mass_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(14, "m_{DV}")

	
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
					self.h['CutFlow'][self.ch].Fill(6)
					self.passed_fidvol_cut = True
			# If this cut doesn't pass, don't continue to check other cuts
			else:
				return

		if self.do_ntrk_cut:
			if self._ntrk_cut():
				if not self.passed_ntrk_cut:
					self.h['CutFlow'][self.ch].Fill(7)
					self.passed_ntrk_cut = True
			else:
				return

		if self.do_opposite_sign_cut or self.do_same_sign_cut:
			if self._charge_cut():
				if not self.passed_charge_cut:
					self.h['CutFlow'][self.ch].Fill(8)
					self.passed_charge_cut = True
			else:
				return

		# self._fill_selected_dv_histos("charge")

		if self.do_dv_type_cut:
			if self._dv_type_cut():
				if not self.passed_dv_type_cut:
					self.h['CutFlow'][self.ch].Fill(9)
					self.passed_dv_type_cut = True
			else:
				return
		self._fill_selected_dv_histos("DVtype")

		if self.do_track_quality_cut:
			if self._track_quality_cut():
				if not self.passed_track_quality_cut:
					self.h['CutFlow'][self.ch].Fill(10)
					self.passed_track_quality_cut = True
			else:
				return

		self._fill_selected_dv_histos("trkqual")

		if self.do_cosmic_veto_cut:
			if self._cosmic_veto_cut():
				if not self.passed_cosmic_veto_cut:
					self.h['CutFlow'][self.ch].Fill(11)
					self.passed_cosmic_veto_cut = True
			else:
				return


		self._fill_selected_dv_histos("cosmic")

		if self.do_trilepton_mass_cut:
			if self._trilepton_mass_cut():
				if not self.passed_trilepton_mass_cut:
					self.h['CutFlow'][self.ch].Fill(12)
					self.passed_trilepton_mass_cut = True
			else:
				return

		self._fill_selected_dv_histos("mlll")

		if self.do_dv_mass_cut:
			if self._dv_mass_cut():
				if not self.passed_dv_mass_cut:
					self.h['CutFlow'][self.ch].Fill(13)
					self.passed_dv_mass_cut = True
			else:
				return

		# Fill histos of truth-matched DVs
		if not self.tree.is_data:
			if self._truth_match():
				self.h['CutFlow'][self.ch].Fill(14)
				self._fill_selected_dv_histos("match")

		# Fill all the histograms with only selected DVs. (ie. the ones that pass the full selection)
		self._fill_selected_dv_histos("sel")


class ToyAnalysis(Analysis):
	def __init__(self, tree, vtx_container, selections, outputFile, saveNtuples):
		logger.info('Running  Toy Analysis cuts')
		Analysis.__init__(self, tree, vtx_container, selections, outputFile, saveNtuples)

		self.add('CutFlow', 21, -0.5, 20.5)
		# Bin labels are 1 greater than histogram bins
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(1, "all")
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(2, "PV")
		if self.do_trigger_cut:
			if not self.do_CR:
				self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(3, "trigger")
			else:
				self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(3, "DAOD_RPVLL triggers")
		if self.do_filter_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(4, "%s" % self.filter_type)
		if self.do_prompt_lepton_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(5, "tight prompt %s" % self.plep)
		if self.do_invert_prompt_lepton_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(5, "invert prompt lepton")
		if self.do_ndv_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(6, "DV")
		if self.do_fidvol_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(7, "fiducial")
		if self.do_ntrk_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(8, "%s-track DV" % self.ntrk)
		# if self.do_HNL_mass_cut:
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(9, "dR")
		if self.do_opposite_sign_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(10, "OS DV")
		if self.do_same_sign_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(10, "SS DV")
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(11, "++ DV")
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(12, "-- DV")
		if self.do_dv_type_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(13, "%s DV" % self.dv_type)
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(14, "++ DV")
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(15, "-- DV")
		if self.do_dv_mass_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(16, "m_{DV}")
		if self.do_trilepton_mass_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(17, "m_{lll}")
		if self.do_HNL_pt_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(18, "HNL p_{T}")
		if self.do_cosmic_veto_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(19, "cosmic veto")
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(20, "1-tight")
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(21, "2-tight")

		# 2D correlation plots after each cut in the DV
		sel_list = ["charge", "DVtype", "mDV", "mlll", "HNLpt", "sel"]
		for sel in sel_list:
			self.add2D( sel + '_ntrk', 11, -5.5, 5.5, 9, -0.5, 8.5)
			self.add2D( sel + '_DVmass_mvis', 1000, 0, 500, 1000, 0, 500)
			self.add2D( sel + '_DVmass_mhnl', 1000, 0, 500, 1010, -5, 500)
			self.add2D( sel + '_DVmass_mtrans', 1000, 0, 500, 1000, 0, 500)
			self.add2D( sel + '_DVmass_hnlpt', 1000, 0, 500, 1010, -5, 500)
			self.add2D( sel + '_mvis_mhnl', 1000, 0, 500, 1010, -5, 500)
			self.add2D( sel + '_mvis_mtrans', 1000, 0, 500, 1000, 0, 500)
			self.add2D( sel + '_mvis_hnlpt', 1000, 0, 500, 1010, -5, 500)
			self.add2D( sel + '_mhnl_hnlpt', 1010, -5, 500, 1010, -5, 500)
			self.add2D( sel + '_mhnl_mtrans', 1010, -5, 500, 1000, 0, 500)
			self.add2D( sel + '_mhnl2D', 1010, -5, 500, 1000, 0, 500)

	def _fill_correlation_histos(self, sel):
		w = self.tree.weight
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
					self.h['CutFlow'][self.ch].Fill(6)
					self.passed_fidvol_cut = True
			# If this cut doesn't pass, don't continue to check other cuts
			else:
				return

		if self.do_ntrk_cut:
			if self._ntrk_cut():
				if not self.passed_ntrk_cut:
					self.h['CutFlow'][self.ch].Fill(7)
					self.passed_ntrk_cut = True
			else:
				return

		if self._dR_cut():
			if not self.passed_dR_cut:
				self.h['CutFlow'][self.ch].Fill(8)
				self.passed_dR_cut = True
		else:
			return

		if self.do_opposite_sign_cut or self.do_same_sign_cut:
			if self._charge_cut():
				if not self.passed_charge_cut:
					self.h['CutFlow'][self.ch].Fill(9)
					sign_pair = "SS" if self.do_same_sign_cut else "OS"
					charge_sel = selections.ChargeDV(self.tree, sel=sign_pair)
					if charge_sel.two_plus: 
						self.h['CutFlow'][self.ch].Fill(10)
					if charge_sel.two_minus: 
						self.h['CutFlow'][self.ch].Fill(11)
					self.passed_charge_cut = True
			else:
				return

		# self._fill_selected_dv_histos("charge")
		# self._fill_correlation_histos("charge")



		if self.do_dv_type_cut:
			if self._dv_type_cut():
				if not self.passed_dv_type_cut:
					self.h['CutFlow'][self.ch].Fill(12)
					sign_pair_1 = "SS" if self.do_same_sign_cut else "OS"
					charge_sel_1 = selections.ChargeDV(self.tree, sel=sign_pair_1)
					# print charge_sel.charge_trk1
					# print charge_sel.charge_trk1
					# print "++ ", charge_sel.two_plus
					# print "-- ", charge_sel.two_minus
					if charge_sel_1.two_plus: 
						self.h['CutFlow'][self.ch].Fill(13)
					if charge_sel_1.two_minus: 
						self.h['CutFlow'][self.ch].Fill(14)
					self.passed_dv_type_cut = True
			else:
				return

		self._fill_selected_dv_histos("DVtype")
		self._fill_correlation_histos("DVtype")

		if self.do_dv_mass_cut:
			if self._dv_mass_cut():
				if not self.passed_dv_mass_cut:
					self.h['CutFlow'][self.ch].Fill(15)
					self.passed_dv_mass_cut = True
			else:
				return

		self._fill_selected_dv_histos("mDV")
		self._fill_correlation_histos("mDV")

		if self.do_trilepton_mass_cut:
			if self._trilepton_mass_cut():
				if not self.passed_trilepton_mass_cut:
					self.h['CutFlow'][self.ch].Fill(16)
					self.passed_trilepton_mass_cut = True
			else:
				return
		self._fill_selected_dv_histos("mlll")
		self._fill_correlation_histos("mlll")

		if self.do_HNL_pt_cut:
			if self._HNL_pt_cut():
				if not self.passed_HNL_pt_cut:
					self.h['CutFlow'][self.ch].Fill(17)
					self.passed_HNL_pt_cut = True

			else:
				return

		self._fill_selected_dv_histos("HNLpt")
		self._fill_correlation_histos("HNLpt")

		if self.do_cosmic_veto_cut:
			if self._cosmic_veto_cut():
				if not self.passed_cosmic_veto_cut:
					self.h['CutFlow'][self.ch].Fill(18)
					self.passed_cosmic_veto_cut = True
			else:
				return
		self._fill_selected_dv_histos("cosmic")

		if self._track_quality_cut_1tight():
			if not self.passed_track_1tight_cut:
				self.h['CutFlow'][self.ch].Fill(19)
				self.passed_track_1tight_cut = True
		else:
			return

		self._fill_selected_dv_histos("tight1")

		if self._track_quality_cut_2tight():
			if not self.passed_track_2tight_cut:
				self.h['CutFlow'][self.ch].Fill(20)
				self.passed_track_2tight_cut = True
		else:
			return

		# Fill histos of truth-matched DVs
		if not self.tree.is_data:
			if self._truth_match():
				self.h['CutFlow'][self.ch].Fill(21)
				self._fill_selected_dv_histos("match")

		# Fill all the histograms with only selected DVs. (ie. the ones that pass the full selection)
		self._fill_selected_dv_histos("sel")
		self._fill_correlation_histos("sel")


class FilterCompareData(Analysis):
	def __init__(self, tree, vtx_container, selections, outputFile, saveNtuples):
		Analysis.__init__(self, tree, vtx_container, selections, outputFile, saveNtuples)
		logger.info('Running FilterCompareData Analysis cuts')

		# check extra options
		# prompt lepton cut
		if 'prompt-muon-loose' in self.sel:
			self.plep = 'muon'
			self.do_prompt_lepton_cut = True
			self.plep_quality = 'loose'
		if 'prompt-muon-medium' in self.sel:
			self.plep = 'muon'
			self.do_prompt_lepton_cut = True
			self.plep_quality = 'medium'
		if 'prompt-electron-loose' in self.sel:
			self.plep = 'electron'
			self.do_prompt_lepton_cut = True
			self.plep_quality = 'loose'
		if 'prompt-electron-medium' in self.sel:
			self.plep = 'electron'
			self.do_prompt_lepton_cut = True
			self.plep_quality = 'medium'

		# displaced lepton cut
		if 'displaced-muon' in self.sel:
			self.dlep = 'muon'
			self.do_displaced_lepton_cut = True
		if 'displaced-electron' in self.sel:
			self.dlep = 'electron'
			self.do_displaced_lepton_cut = True


		# Make cutflow
		self.add('CutFlow', 10, -0.5, 9.5)
		# Bin labels are 1 greater than histogram bins
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(1, "all")
		if self.do_trigger_cut:
			if not self.do_CR:
				self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(2, "trigger")
			else:
				self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(2, "DAOD_RPVLL triggers")
		if self.do_filter_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(3, "%s" % self.filter_type)
		if self.do_prompt_lepton_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(4, "prompt {} {}".format(self.plep, self.plep_quality) )
		if self.do_displaced_lepton_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(5, "displaced {} loose".format(self.dlep))
		if self.do_displaced_lepton_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(6, "displaced {} medium".format(self.dlep))
		if self.do_ndv_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(7, "nDV")
		# if self.do_fidvol_cut:
		# 	self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(7, "fiducial")
		# if self.do_ntrk_cut:
		# 	self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(8, "%s-track DV" % self.ntrk)


	def initialize_cut_bools(self):
		###########################################################################################################################
		# Cut bools that will be intialized in the pre selection for every event. These bools tell the code if the cutflow has already been filled for this event.
		# Default is to select the first event that passes the selection.
		###########################################################################################################################
		self.passed_preselection_cuts = False
		self.passed_displaced_lepton_cut_loose = False
		self.passed_displaced_lepton_cut_medium = False


	def _displaced_lepton_cut(self, quality):
		displaced_lepton = selections.DisplacedLepton(self.tree, lepton=self.dlep, quality=quality)
		return displaced_lepton.passes()

	def preSelection(self):
		######################################################################################################
		# Preselection are all the cuts that are requied per event
		# Current cuts include: trigger, filter, plepton, DV cut
		######################################################################################################

		#initialize the cut bools for every event
		self.initialize_cut_bools()

		self._fill_leptons()

		if not self.tree.is_data:
			self._fill_truth_histos(sel='truth_all')

		self.h['CutFlow'][self.ch].SetBinContent(1, self.tree.max_entries)  # all events

		######################################################################################################
		# Selection code is deisgned so that it will pass the selection only if the cut true or cut is unused
		# ex. passTrigger is true if the trigcut is true OR if trigcut is not used)
		######################################################################################################

		if self.do_trigger_cut:
			if self._trigger_cut():
				# Fill the plot at the specified bin
				self.h['CutFlow'][self.ch].Fill(1)
			else:
				return

		if self.do_filter_cut:
			if self._filter_cut():
				self.h['CutFlow'][self.ch].Fill(2)
			else:
				return

		if self.do_prompt_lepton_cut:
			if self._prompt_lepton_cut(quality=self.plep_quality):
				self.h['CutFlow'][self.ch].Fill(3)
			else:
				return

		if self.do_displaced_lepton_cut:
			if self._displaced_lepton_cut(quality='loose'):
				print("passes loose displaced_lepton")
				self.h['CutFlow'][self.ch].Fill(4)
			else:
				return

		if self.do_displaced_lepton_cut:
			if self._displaced_lepton_cut(quality='medium'):
				self.h['CutFlow'][self.ch].Fill(5)
			else:
				return

		if self.do_ndv_cut:
			if self._ndv_cut():
				self.h['CutFlow'][self.ch].Fill(6)
			else:
				return

		self.passed_preselection_cuts = True

	def DVSelection(self):
		# only do the DV selection if the preselction was passed for the event.
		if not self.passed_preselection_cuts:
			return

		# There is an extra bit of logic here since we look at several DVs
		# but only want to fill the cutflow once per event

class KShort(Analysis):
	def __init__(self, tree, vtx_container, selections, outputFile, saveNtuples):
		Analysis.__init__(self, tree, vtx_container, selections, outputFile, saveNtuples)
		logger.info('Running KShort Analysis cuts')

		self.add('CutFlow', 17, -0.5, 16.5)
		# Bin labels are 1 greater than histogram bins
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(1, "all")
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(2, "PV")
		# if self.do_trigger_cut:
		# 	self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(3, "trigger")
		# if self.do_filter_cut:
		# 	self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(4, "%s" % self.filter_type)
		# if self.do_prompt_lepton_cut:
		# 	self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(5, "tight prompt %s" % self.plep)
		if self.do_invert_prompt_lepton_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(3, "invert prompt lepton")
		if self.do_alpha_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(4, "alpha")
		if self.do_mass_window_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(5, "K0 mass")
		# add other histograms
		self.add2D('charge_ntrk', 11, -5.5, 5.5, 9, -0.5, 8.5)

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

	# def _fill_truth_dv_histos(self,sel):
	# 	if not self.tree.is_data:
	# 		truthInfo = helpers.Truth()
	# 		truthInfo.getTruthParticles(self.tree)
	# 		self.fill_hist('truth_all', 'DV_r', truthInfo.truth_dvr)
	# 		self.fill_hist('truth_all', 'DV_x', truthInfo.truth_dvx)
	# 		self.fill_hist('truth_all', 'DV_y', truthInfo.truth_dvy)
	# 		self.fill_hist('truth_all', 'DV_z', truthInfo.truth_dvz)
			
	# 		if sel == self.saveNtuples or self.saveNtuples == 'allcuts': 
	# 			self.micro_ntuples['truth_all'].fill()

	def _fill_selected_dv_histos(self, sel, do_lock=True):
		if self._locked < FILL_LOCKED or not do_lock:
			# these are the histograms you only want to fill ONCE per DV

			# sel refers to the last selection that was applied
			# self.fill_hist(sel, 'charge_ntrk', self.tree.dv('charge'), self.tree.dv('ntrk'))
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
			self.fill_hist(sel, 'DV_alpha', selections.Alpha(self.tree).alpha)

			# kshort stuff
			track_sum = selections.SumTrack(self.tree)
			self.fill_hist(sel, 'DV_sum_track_pt', track_sum.sum_track_pt)
			self.fill_hist(sel, 'DV_sum_track_pt_wrt_pv', track_sum.sum_track_pt_wrt_pv)
			self.fill_hist(sel, 'DV_sum_track_pt_diff', track_sum.sum_track_pt_wrt_pv - track_sum.sum_track_pt)
			self.fill_hist(sel, 'DV_sum_track_charge', track_sum.sum_track_charge)

			if sel == self.saveNtuples or self.saveNtuples == 'allcuts': 
				self.micro_ntuples[sel].fill()

	# def _fill_selected_dv_ntuples(self, sel, do_lock=True):
	# 	if self._locked < FILL_LOCKED or not do_lock:
	#
	# 		# Fill micro ntuples
	# 		self.fill_ntuple(sel, 'DV_num_trks', self.tree.dv('ntrk'))
	# 		self.fill_ntuple(sel, 'DV_x', self.tree.dv('x'))
	# 		self.fill_ntuple(sel, 'DV_y', self.tree.dv('y'))
	# 		self.fill_ntuple(sel, 'DV_z', self.tree.dv('z'))
	# 		self.fill_ntuple(sel, 'DV_r', self.tree.dv('r'))
	# 		self.fill_ntuple(sel, 'DV_distFromPV', self.tree.dv('distFromPV'))
	# 		self.fill_ntuple(sel, 'DV_mass', self.tree.dv('mass'))
	# 		self.fill_ntuple(sel, 'DV_pt', self.tree.dv('pt'))
	# 		self.fill_ntuple(sel, 'DV_eta', self.tree.dv('eta'))
	# 		self.fill_ntuple(sel, 'DV_phi', self.tree.dv('phi'))
	# 		self.fill_ntuple(sel, 'DV_alpha', selections.Alpha(self.tree).alpha)
	#
	# 		track_sum = selections.SumTrack(self.tree)
	# 		self.fill_ntuple(sel, 'sum_track_pt', track_sum.sum_track_pt)
	# 		self.fill_ntuple(sel, 'DV_sum_track_pt', track_sum.sum_track_pt)
	# 		self.fill_ntuple(sel, 'DV_sum_track_pt_wrt_pv', track_sum.sum_track_pt_wrt_pv)
	# 		self.fill_ntuple(sel, 'DV_sum_track_pt_diff', track_sum.sum_track_pt_wrt_pv - track_sum.sum_track_pt)
	# 		self.fill_ntuple(sel, 'DV_sum_track_charge', track_sum.sum_track_charge)


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

		# if self.do_trigger_cut:
		# 	if self._trigger_cut():
		# 		# Fill the plot at the specified bin
		# 		self.h['CutFlow'][self.ch].Fill(2)
		# 	else:
		# 		return

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
