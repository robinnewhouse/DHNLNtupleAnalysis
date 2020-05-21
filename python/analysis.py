# from ROOT import*
import ROOT
import numpy as np
import os
import sys
import helpers
import selections
import treenames
import observables
import logging
logger = helpers.getLogger('dHNLAnalysis.analysis',level = logging.INFO)





UNLOCKED = 0
SELECTION_LOCKED = 1
FILL_LOCKED = 2


class Analysis(object):
	ch = ''
	h = {}  # map of histogram names to map of systematics to histograms
	blinded = True
	mapSel = {}

	def __init__(self, channel, selections, outputFile, isdata):
		self.sel = selections
		self._outputFile = outputFile
		self.fi = ROOT.TFile.Open(outputFile, 'update')
		self.ch = channel
		self.histSuffixes = [self.ch]
		self.h = {}
		self.DAOD_RPVLL_triggers = []
		self.inverted_triggers = []
		# make histograms (common for all channels)
		# self.add('CutFlow', 16, -0.5, 15.5)

		self._locked = UNLOCKED
		self._dv_truth_locked = UNLOCKED

		self.observables = [observable.registered(self) for observable in observables.ObservableList if ((observable.only is None) or any(only in self.sel for only in observable.only))]
		for observable in self.observables:
			if 'hist' in observable.do:
				if isdata and observable.need_truth: 
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
			logger.warn('You did not specify a prompt lepton for this channel. Skipping prompt lepton selection.')
			self.do_prompt_lepton_cut = False

		if 'CR' in self.sel: # DO NOT CHANGE THESE CUTS OR YOU MIGHT UNBLIND DATA!!!
			self.do_CR = True
			self.do_trigger_cut = False #do not apply trigger cut
			self.do_invert_trigger_cut = False #do not apply inverted trigger cut
			self.do_filter_cut = False #do not apply filter cut
			self.do_prompt_lepton_cut = False #do not apply prompt lepton cut
			self.do_invert_prompt_lepton_cut = True # invert prompt lepton cut
			logger.info('You are setup up to look in the inverted prompt lepton control region!')
		else: 
			self.do_CR = False
			self.do_invert_prompt_lepton_cut = False
			self.do_invert_trigger_cut = False

		# nDV cut
		self.do_ndv_cut = ('nDV' in self.sel)
		if not self.do_ndv_cut: logger.warn('You did not add nDV cut. Skipping nDV selection.')

		# fiducial volume
		self.do_fidvol_cut = 'fidvol' in self.sel
		if not self.do_fidvol_cut: logger.warn('You did not add DV in fiducial volume cut. Skipping DV in fiducial volume selection.')

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
		if not (self.do_opposite_sign_cut or self.do_same_sign_cut): logger.warn('You did not add an SS or OS track cut. Skipping SS/OS track selection.')

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
		elif '2-medium' in self.sel:
			self.track_quality = '2-medium'
		elif '2-loose' in self.sel:
			self.track_quality = '2-loose'
		elif '2-veryloose' in self.sel:
			self.track_quality = '2-veryloose'
		elif 'loose-veryloose' in self.sel:
			self.track_quality = 'loose-veryloose'
		elif 'medium-veryloose' in self.sel:
			self.track_quality = 'medium-veryloose'
		elif 'tight-veryloose' in self.sel:
			self.track_quality = 'tight-veryloose'
		elif '2-any' in self.sel:
			self.track_quality = '2-any'
		else:
			logger.warn('You did not specify a DV track quality for this channel. Skipping DV track quality selection.')
			self.do_track_quality_cut = False

		# cosmic veto cut
		self.do_cosmic_veto_cut = 'cosmicveto' in self.sel
		if not self.do_cosmic_veto_cut: logger.warn('You did not add a cosmic veto cut for this channel. Skipping cosmic veto selection.')

		# tri-lepton mass cut
		self.do_trilepton_mass_cut = 'mlll' in self.sel
		if not self.do_trilepton_mass_cut and "CR" not in self.sel: logger.warn('You did not add a mlll cut for this channel. Skipping tri-lepton mass selection.')

		# DV mass cut
		self.do_dv_mass_cut = 'DVmass' in self.sel
		if not self.do_dv_mass_cut: logger.warn('You did not add a DVmass cut for this channel. Skipping displaced vertex mass selection.')

		#HNL mass cut
		self.do_HNL_mass_cut = 'HNLmass' in self.sel

		#HNL pT cut
		self.do_HNL_pt_cut = 'HNLpt' in self.sel

		self.check_input_consistency()


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
		self._dv_truth_locked = UNLOCKED

	def add(self, hName, nBins, xLow, xHigh):
		self.h[hName] = {}
		for s in self.histSuffixes:
			self.h[hName][s] = ROOT.TH1D(hName+"_"+s, "", nBins, xLow, xHigh)
			self.h[hName][s].Sumw2()
			self.h[hName][s].SetDirectory(0)

	def add2D(self, hName, nBins, xLow, xHigh, nBinsY, yLow, yHigh):
		self.h[hName] = {}
		for s in self.histSuffixes:
			self.h[hName][s] = ROOT.TH2D(hName + "_" + s, "", nBins, xLow, xHigh, nBinsY, yLow, yHigh)
			self.h[hName][s].Sumw2()
			self.h[hName][s].SetDirectory(0)

	def write(self):
		self.fi.cd()
		for hName in self.h:
			for s in self.histSuffixes:
				self.h[hName][s].Write(hName + '_' + s)
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
			for s in self.histSuffixes:
				if self.h[histName][s].GetEntries() == 0:
					logger.debug('\tUnfilled HIST({}<{}>)!'.format(histName, s))

		# y_err = ROOT.Double()
		# for s in self.histSuffixes:
		# 	h = self.h['finalNEvents'][s]
		# 	logger.info('\tFinal NEvents <{syst}>: {y}'.format(syst = s, y = h.Integral()))
		# 	h = self.h['finalYields'][s]
		# 	logger.info('\tFinal Yields <{syst}>: {y:.4g}+/-{y_err:.2g}'.format(syst = s, y = h.IntegralAndError(0, h.GetNbinsX()+1, y_err), y_err = y_err))
		self.write()

		# head, sep, tail = self._outputFile.partition('file://')
		# f = tail if head == '' else self._outputFile
		# f = tail if head == '' else self._outputFile
		# try:
		# 	os.rename(f + '.part', f)
		# except OSError as e:
		# 	logger.error(e, exc_info=True)

	def _preSelection(self, evt):
		raise NotImplementedError

	def preSelection(self, evt):
		presel = self._preSelection(evt)
		return presel

	def _DVSelection(self, evt):
		raise NotImplementedError

	def DVSelection(self, evt):
		self._DVSelection(evt)

	# Protected function to create the selection object and return its success
	# The selection object may be used to fill additional histograms (see _prompt_lepton_cut)
	
	def _pv_cut(self, evt):
		pv_sel = selections.PV(evt=evt)
		return pv_sel.passes()

	def _trigger_cut(self, evt):
		trigger_sel = selections.Trigger(evt=evt, trigger=self.trigger)
		return trigger_sel.passes()

	def _invert_trigger_cut(self, evt):
		trigger_sel = selections.Trigger(evt=evt, trigger=self.trigger, invert= True)
		return trigger_sel.passes()

	def _nmuon_cut(self, evt):
		n_muons = len(evt.tree.muonpt[evt.ievt])
		ncomb_muon = 0
		for imu in range(n_muons): 
			mutype = evt.tree.muontype[evt.ievt][imu]
			if mutype == 0:
				ncomb_muon +=1
		
		return ncomb_muon !=  3
		# return ncomb_muon <  3

	def _filter_cut(self, evt):
		filter_sel = selections.Filter(evt=evt, filter_type=self.filter_type)
		return filter_sel.passes()

	def _prompt_lepton_cut(self, evt):
		self.plep_sel = selections.PromptLepton(evt=evt, lepton=self.plep)

		# Add to histogram all prompt leptons that pass selection.
		# If _prompt_lepton_cut() is run after trigger and filter cut then those cuts will also be applied.
		if self.plep_sel.passes(): 
			self.h["all_plep_pt"][self.ch].Fill(self.plep_sel.plepVec.Pt(), evt.weight)
			self.h["all_plep_eta"][self.ch].Fill(self.plep_sel.plepVec.Eta(), evt.weight)
			self.h["all_plep_phi"][self.ch].Fill(self.plep_sel.plepVec.Phi(), evt.weight)
			self.h["all_plep_d0"][self.ch].Fill(self.plep_sel.plepd0, evt.weight)
			self.h["all_plep_z0"][self.ch].Fill(self.plep_sel.plepz0, evt.weight)

		return self.plep_sel.passes()

	def _invert_prompt_lepton_cut(self, evt):
		self.invt_lep = selections.InvertedPromptLepton(evt=evt) #invert prompt lepton selection
		return self.invt_lep.passes()


	def _ndv_cut(self, evt):
		dv_sel = selections.nDV(evt=evt)
		return dv_sel.passes()

	def _fidvol_cut(self, evt):
		fidvol_sel = selections.DVradius(evt=evt)
		return fidvol_sel.passes()

	def _ntrk_cut(self, evt):
		ntracks_sel = selections.DVntracks(evt=evt, ntrk=self.ntrk)
		return ntracks_sel.passes()

	def _charge_cut(self, evt):
		sign_pair = "SS" if self.do_same_sign_cut else "OS"
		charge_sel = selections.ChargeDV(evt=evt, sel=sign_pair)
		return charge_sel.passes()

	def _dv_type_cut(self, evt):
		dv_sel = selections.DVtype(evt=evt, dv_type=self.dv_type)
		return dv_sel.passes()

	def _track_quality_cut(self, evt):
		track_quality_sel = selections.Trackqual(evt=evt, quality=self.track_quality)
		return track_quality_sel.passes()

	def _cosmic_veto_cut(self, evt):
		cosmic_veto_sel = selections.Cosmicveto(evt=evt)
		# self.h["DV_trk_sep"][self.ch].Fill(cosmic_veto_sel.separation)
		return cosmic_veto_sel.passes()

	def _trilepton_mass_cut(self, evt):
		plep_vec = self.plep_sel.plepVec

		muons = helpers.Tracks()
		muons.getMuons(evt=evt)
		muVec = muons.lepVec

		electrons = helpers.Tracks()
		electrons.getElectrons(evt=evt)
		elVec = electrons.lepVec

		mlll_sel = selections.Mlll(dv_type=self.dv_type, plep=plep_vec, dMu=muVec, dEl=elVec)
		return mlll_sel.passes()

	def _dv_mass_cut(self, evt):
		dv_mass_sel = selections.DVmass(evt=evt,dvmasscut=4)
		return dv_mass_sel.passes()


	def _fill_histos(self, evt):
		w = evt.weight
		for imu in range(len(evt.tree.muontype[evt.ievt])):
			self.h["muon_type"][self.ch].Fill(evt.tree.muontype[evt.ievt][imu], w)
			self.h["muon_pt"][self.ch].Fill(evt.tree.muonpt[evt.ievt][imu], w)
			self.h["muon_eta"][self.ch].Fill(evt.tree.muoneta[evt.ievt][imu], w)
			self.h["muon_phi"][self.ch].Fill(evt.tree.muonphi[evt.ievt][imu], w)

		for iel in range(len(evt.tree.elpt[evt.ievt])):
			self.h["el_pt"][self.ch].Fill(evt.tree.elpt[evt.ievt][iel], w)
			self.h["el_eta"][self.ch].Fill(evt.tree.eleta[evt.ievt][iel], w)
			self.h["el_phi"][self.ch].Fill(evt.tree.elphi[evt.ievt][iel], w)

		if evt.tree.isData == True: 
			pass
		else: 
			truthInfo = helpers.Truth()
			truthInfo.getTruthParticles(evt)
			self.h["truth_all_W_pt"][self.ch].Fill(truthInfo.W_vec.Pt(), w)
			self.h["truth_all_W_eta"][self.ch].Fill(truthInfo.W_vec.Eta(), w)
			self.h["truth_all_W_phi"][self.ch].Fill(truthInfo.W_vec.Phi(), w)
			self.h["truth_all_W_mass"][self.ch].Fill(truthInfo.W_vec.M(), w)
			self.h["truth_all_HNL_pt"][self.ch].Fill(truthInfo.HNL_vec.Pt(), w)
			self.h["truth_all_HNL_eta"][self.ch].Fill(truthInfo.HNL_vec.Eta(), w)
			self.h["truth_all_HNL_phi"][self.ch].Fill(truthInfo.HNL_vec.Phi(), w)
			self.h["truth_all_HNL_mass"][self.ch].Fill(truthInfo.HNL_vec.M(), w)

			self.h["truth_all_mHNLcalc"][self.ch].Fill(truthInfo.mhnl, w)

			self.h["truth_all_DV_r"][self.ch].Fill(truthInfo.truth_dvr, w)
			self.h["truth_all_DV_x"][self.ch].Fill(truthInfo.truth_dvx, w)
			self.h["truth_all_DV_y"][self.ch].Fill(truthInfo.truth_dvy, w)
			self.h["truth_all_DV_z"][self.ch].Fill(truthInfo.truth_dvz, w)
			self.h["truth_all_plep_pt"][self.ch].Fill(truthInfo.plep_vec.Pt(), w)
			self.h["truth_all_plep_eta"][self.ch].Fill(truthInfo.plep_vec.Eta(), w)
			self.h["truth_all_plep_phi"][self.ch].Fill(truthInfo.plep_vec.Phi(), w)
			self.h["truth_all_plep_mass"][self.ch].Fill(truthInfo.plep_vec.M(), w)
			self.h["truth_all_lep1_trk_pt"][self.ch].Fill(truthInfo.trkVec[0].Pt(), w)
			self.h["truth_all_lep1_trk_eta"][self.ch].Fill(truthInfo.trkVec[0].Eta(), w)
			self.h["truth_all_lep1_trk_phi"][self.ch].Fill(truthInfo.trkVec[0].Phi(), w)
			self.h["truth_all_lep2_trk_pt"][self.ch].Fill(truthInfo.trkVec[1].Pt(), w)
			self.h["truth_all_lep2_trk_eta"][self.ch].Fill(truthInfo.trkVec[1].Eta(), w)
			self.h["truth_all_lep2_trk_phi"][self.ch].Fill(truthInfo.trkVec[1].Phi(), w)

			for itrk in xrange(2): 
				self.h["truth_all_DV_trk_pt"][self.ch].Fill(truthInfo.trkVec[itrk].Pt(), w)
				self.h["truth_all_DV_trk_eta"][self.ch].Fill(truthInfo.trkVec[itrk].Eta(), w)
				self.h["truth_all_DV_trk_phi"][self.ch].Fill(truthInfo.trkVec[itrk].Phi(), w)


	def _fill_all_dv_histos(self, evt):
		w = evt.weight
		self.h["charge_ntrk"][self.ch].Fill(evt.tree.dvcharge[evt.ievt][evt.idv], evt.tree.dvntrk[evt.ievt][evt.idv], w)
		ntracks = len(evt.tree.trackd0[evt.ievt][evt.idv])
		for itrk in range(ntracks):  # loop over tracks
			self.h["all_DV_trk_pt"][self.ch].Fill(evt.tree.trackpt[evt.ievt][evt.idv][itrk], w )
			self.h["all_DV_trk_eta"][self.ch].Fill(evt.tree.tracketa[evt.ievt][evt.idv][itrk], w )
			self.h["all_DV_trk_phi"][self.ch].Fill(evt.tree.trackphi[evt.ievt][evt.idv][itrk], w )
			self.h["all_DV_trk_d0"][self.ch].Fill(evt.tree.trackd0[evt.ievt][evt.idv][itrk], w )
			self.h["all_DV_trk_z0"][self.ch].Fill(evt.tree.trackz0[evt.ievt][evt.idv][itrk], w )
			self.h["all_DV_trk_charge"][self.ch].Fill(evt.tree.trackcharge[evt.ievt][evt.idv][itrk], w )
			self.h["all_DV_trk_chi2"][self.ch].Fill(evt.tree.trackchi2[evt.ievt][evt.idv][itrk], w )

		self.h["all_DV_num_trks"][self.ch].Fill(evt.tree.dvntrk[evt.ievt][evt.idv], w)
		self.h["all_DV_x"][self.ch].Fill(evt.tree.dvx[evt.ievt][evt.idv], w)
		self.h["all_DV_y"][self.ch].Fill(evt.tree.dvy[evt.ievt][evt.idv], w)
		self.h["all_DV_z"][self.ch].Fill(evt.tree.dvz[evt.ievt][evt.idv], w)
		self.h["all_DV_r"][self.ch].Fill(evt.tree.dvr[evt.ievt][evt.idv], w)
		self.h["all_DV_distFromPV"][self.ch].Fill(evt.tree.dvdistFromPV[evt.ievt][evt.idv], w)
		self.h["all_DV_mass"][self.ch].Fill(evt.tree.dvmass[evt.ievt][evt.idv], w)
		self.h["all_DV_pt"][self.ch].Fill(evt.tree.dvpt[evt.ievt][evt.idv], w)
		self.h["all_DV_eta"][self.ch].Fill(evt.tree.dveta[evt.ievt][evt.idv], w)
		self.h["all_DV_phi"][self.ch].Fill(evt.tree.dvphi[evt.ievt][evt.idv], w)
		self.h["all_DV_minOpAng"][self.ch].Fill(evt.tree.dvminOpAng[evt.ievt][evt.idv], w)
		self.h["all_DV_maxOpAng"][self.ch].Fill(evt.tree.dvmaxOpAng[evt.ievt][evt.idv], w)
		self.h["all_DV_charge"][self.ch].Fill(evt.tree.dvcharge[evt.ievt][evt.idv], w)
		self.h["all_DV_chi2"][self.ch].Fill(evt.tree.dvchi2[evt.ievt][evt.idv], w)

	def _fill_selected_dv_histos(self, evt, sel):
		w = evt.weight
		if self._locked < FILL_LOCKED:
			# these are the histograms you only want to fill ONCE per DV
			# sel refers to the last selection that was applied 
				
			if self.do_prompt_lepton_cut:
				plep_vec = self.plep_sel.plepVec
				plepd0 = self.plep_sel.plepd0
				plepz0 = self.plep_sel.plepz0

				self.h[sel + "_plep_pt"][self.ch].Fill(plep_vec.Pt(), w)
				self.h[sel + "_plep_eta"][self.ch].Fill(plep_vec.Eta(), w)
				self.h[sel + "_plep_phi"][self.ch].Fill(plep_vec.Phi(), w)
				self.h[sel + "_plep_d0"][self.ch].Fill(plepd0, w)
				self.h[sel + "_plep_z0"][self.ch].Fill(plepz0, w)
				
				tracks = helpers.Tracks()
				tracks.getTracks(evt=evt)
				trkVec = tracks.lepVec	

				if tracks.ntracks == 2: 
					Mltt = selections.Mltt(plep=plep_vec, trks=trkVec)
					Mhnl = selections.Mhnl(evt=evt, plep=plep_vec, trks =trkVec )
					Mtrans = selections.Mtrans(plep=plep_vec, trks =trkVec )

					self.h[sel + "_mvis"][self.ch].Fill(Mltt.mltt, w)
					self.h[sel + "_HNLm"][self.ch].Fill(Mhnl.mhnl, w)
					self.h[sel + "_HNLm2"][self.ch].Fill(Mhnl.mhnl2, w)
					self.h[sel + "_HNLpt"][self.ch].Fill(Mhnl.hnlpt, w)
					self.h[sel + "_HNLeta"][self.ch].Fill(Mhnl.hnleta, w)
					self.h[sel + "_HNLphi"][self.ch].Fill(Mhnl.hnlphi, w)
					self.h[sel + "_mtrans"][self.ch].Fill(Mtrans.mtrans, w)
					self.h[sel + "_mtrans_rot"][self.ch].Fill(Mhnl.mtrans_rot, w)
					

					deta = abs(tracks.eta[0] - tracks.eta[1])
					dphi = abs(tracks.lepVec[0].DeltaPhi(tracks.lepVec[1]))
					dpt = abs(tracks.pt[0] - tracks.pt[1])
					dR = tracks.lepVec[0].DeltaR(tracks.lepVec[1])

					if dR == 0.0:
						self.h[sel + "_DV_redmass"][self.ch].Fill(-1, w)
						self.h[sel + "_DV_redmassvis"][self.ch].Fill(-1, w)
						self.h[sel + "_DV_redmassHNL"][self.ch].Fill(-1, w)
					else:
						self.h[sel + "_DV_redmass"][self.ch].Fill(evt.tree.dvmass[evt.ievt][evt.idv]/dR, w)
						self.h[sel + "_DV_redmassvis"][self.ch].Fill(Mltt.mltt/dR, w)
						self.h[sel + "_DV_redmassHNL"][self.ch].Fill(Mltt.mltt/dR, w)
					
					self.h[sel + "_DV_trk_deta"][self.ch].Fill(deta, w)
					self.h[sel + "_DV_trk_dphi"][self.ch].Fill(dphi, w)
					self.h[sel + "_DV_trk_dpt"][self.ch].Fill(dpt, w)
					self.h[sel + "_DV_trk_dR"][self.ch].Fill(dR, w)

			ntracks = len(evt.tree.trackd0[evt.ievt][evt.idv])
			for itrk in range(ntracks):  # loop over tracks
				self.h[sel + "_DV_trk_pt"][self.ch].Fill(evt.tree.trackpt[evt.ievt][evt.idv][itrk], w)
				self.h[sel + "_DV_trk_eta"][self.ch].Fill(evt.tree.tracketa[evt.ievt][evt.idv][itrk], w)
				self.h[sel + "_DV_trk_phi"][self.ch].Fill(evt.tree.trackphi[evt.ievt][evt.idv][itrk], w)
				self.h[sel + "_DV_trk_d0"][self.ch].Fill(evt.tree.trackd0[evt.ievt][evt.idv][itrk], w)
				self.h[sel + "_DV_trk_z0"][self.ch].Fill(evt.tree.trackz0[evt.ievt][evt.idv][itrk], w)
				self.h[sel + "_DV_trk_charge"][self.ch].Fill(evt.tree.trackcharge[evt.ievt][evt.idv][itrk], w)
				self.h[sel + "_DV_trk_chi2"][self.ch].Fill(evt.tree.trackchi2[evt.ievt][evt.idv][itrk], w)

			self.h[sel + "_DV_num_trks"][self.ch].Fill(evt.tree.dvntrk[evt.ievt][evt.idv], w)
			self.h[sel + "_DV_x"][self.ch].Fill(evt.tree.dvx[evt.ievt][evt.idv], w)
			self.h[sel + "_DV_y"][self.ch].Fill(evt.tree.dvy[evt.ievt][evt.idv], w)
			self.h[sel + "_DV_z"][self.ch].Fill(evt.tree.dvz[evt.ievt][evt.idv], w)
			self.h[sel + "_DV_r"][self.ch].Fill(evt.tree.dvr[evt.ievt][evt.idv], w)
			self.h[sel + "_DV_distFromPV"][self.ch].Fill(evt.tree.dvdistFromPV[evt.ievt][evt.idv], w)
			self.h[sel + "_DV_mass"][self.ch].Fill(evt.tree.dvmass[evt.ievt][evt.idv], w)
			self.h[sel + "_DV_pt"][self.ch].Fill(evt.tree.dvpt[evt.ievt][evt.idv], w)
			self.h[sel + "_DV_eta"][self.ch].Fill(evt.tree.dveta[evt.ievt][evt.idv], w)
			self.h[sel + "_DV_phi"][self.ch].Fill(evt.tree.dvphi[evt.ievt][evt.idv], w)
			self.h[sel + "_DV_minOpAng"][self.ch].Fill(evt.tree.dvminOpAng[evt.ievt][evt.idv], w)
			self.h[sel + "_DV_maxOpAng"][self.ch].Fill(evt.tree.dvmaxOpAng[evt.ievt][evt.idv], w)
			self.h[sel + "_DV_charge"][self.ch].Fill(evt.tree.dvcharge[evt.ievt][evt.idv], w)
			self.h[sel + "_DV_chi2"][self.ch].Fill(evt.tree.dvchi2[evt.ievt][evt.idv], w)

			 # better to fill truth matched DVs... need to fix this -DT
			if evt.tree.isData:
				pass
			else: 
				truthInfo = helpers.Truth()
				truthInfo.getTruthParticles(evt)
				self.h["truth_" + sel + "_W_pt"][self.ch].Fill(truthInfo.W_vec.Pt(), w)
				self.h["truth_" + sel + "_W_eta"][self.ch].Fill(truthInfo.W_vec.Eta(), w)
				self.h["truth_" + sel + "_W_phi"][self.ch].Fill(truthInfo.W_vec.Phi(), w)
				self.h["truth_" + sel + "_W_mass"][self.ch].Fill(truthInfo.W_vec.M(), w)
				self.h["truth_" + sel + "_HNL_pt"][self.ch].Fill(truthInfo.HNL_vec.Pt(), w)
				self.h["truth_" + sel + "_HNL_eta"][self.ch].Fill(truthInfo.HNL_vec.Eta(), w)
				self.h["truth_" + sel + "_HNL_phi"][self.ch].Fill(truthInfo.HNL_vec.Phi(), w)
				self.h["truth_" + sel + "_HNL_mass"][self.ch].Fill(truthInfo.HNL_vec.M(), w)

				self.h["truth_" + sel + "_mHNLcalc"][self.ch].Fill(truthInfo.mhnl, w)

				self.h["truth_" + sel + "_DV_r"][self.ch].Fill(truthInfo.truth_dvr, w)
				self.h["truth_" + sel + "_DV_x"][self.ch].Fill(truthInfo.truth_dvx, w)
				self.h["truth_" + sel + "_DV_y"][self.ch].Fill(truthInfo.truth_dvy, w)
				self.h["truth_" + sel + "_DV_z"][self.ch].Fill(truthInfo.truth_dvz, w)
				self.h["truth_" + sel + "_DV_r"][self.ch].Fill(truthInfo.truth_dvr, w)
				self.h["truth_" + sel + "_plep_pt"][self.ch].Fill(truthInfo.plep_vec.Pt(), w)
				self.h["truth_" + sel + "_plep_eta"][self.ch].Fill(truthInfo.plep_vec.Eta(), w)
				self.h["truth_" + sel + "_plep_phi"][self.ch].Fill(truthInfo.plep_vec.Phi(), w)
				self.h["truth_" + sel + "_plep_mass"][self.ch].Fill(truthInfo.plep_vec.M(), w)
				self.h["truth_" + sel + "_lep1_trk_pt"][self.ch].Fill(truthInfo.trkVec[0].Pt(), w)
				self.h["truth_" + sel + "_lep1_trk_eta"][self.ch].Fill(truthInfo.trkVec[0].Eta(), w)
				self.h["truth_" + sel + "_lep1_trk_phi"][self.ch].Fill(truthInfo.trkVec[0].Phi(), w)

				self.h["truth_" + sel + "_lep2_trk_pt"][self.ch].Fill(truthInfo.trkVec[1].Pt(), w)
				self.h["truth_" + sel + "_lep2_trk_eta"][self.ch].Fill(truthInfo.trkVec[1].Eta(), w)
				self.h["truth_" + sel + "_lep2_trk_phi"][self.ch].Fill(truthInfo.trkVec[1].Phi(), w)
				for itrk in xrange(2): 
					self.h["truth_" + sel + "_DV_trk_pt"][self.ch].Fill(truthInfo.trkVec[itrk].Pt(), w)
					self.h["truth_" + sel + "_DV_trk_eta"][self.ch].Fill(truthInfo.trkVec[itrk].Eta(), w)
					self.h["truth_" + sel + "_DV_trk_phi"][self.ch].Fill(truthInfo.trkVec[itrk].Phi(), w)
				

			if sel == "sel": 
				self._locked = FILL_LOCKED  # this only becomes unlocked after the event loop finishes in makeHistograms so you can only fill one DV from each event.



class oldAnalysis(Analysis):

	def __init__(self, channel, selections, outputFile, isdata):
		logger.info('Running old analysis cuts')
		Analysis.__init__(self, channel, selections, outputFile, isdata)

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

	def _preSelection(self, evt):
		######################################################################################################
		# Preselection are all the cuts that are requied per event
		# Current cuts include: trigger, filter, plepton, DV cut
		######################################################################################################


		###########################################################################################################################
		# Initialize the cut bools every event. These bools tell the code if the cutflow has already been filled for this event.
		# Default is to select the first event that passes the selection
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

		self._fill_histos(evt)

		# triggerList = evt.tree.passedtriggers[evt.ievt]
		# # for triggerList in passedTriggers:
		# for trigger in triggerList:
		# 	if trigger not in self.DAOD_RPVLL_triggers:
		# 		self.DAOD_RPVLL_triggers.append(trigger)
		# n = len(evt.tree.passedtriggers) - 1 
		# if evt.ievt == n:
		# 	print(self.DAOD_RPVLL_triggers)

		self.h['CutFlow'][self.ch].SetBinContent(1, evt.tree.allEvt)


		######################################################################################################
		# Selection code is deisgned so that it will pass the selection only if the cut true or cut is unused
		# ex. passTrigger is true if the trigcut is true OR if trigcut is not used)
		######################################################################################################

		if self._pv_cut(evt): #Check to make sure event has a PV otherwise throw event away (this happens very rarely with data).
			self.h['CutFlow'][self.ch].Fill(1)
		else:
			return

		if self.do_trigger_cut:
			if self._trigger_cut(evt):
				# Fill the plot at the specified bin
				self.h['CutFlow'][self.ch].Fill(2)
			else:
				return

		if self.do_invert_trigger_cut:
			if self._invert_trigger_cut(evt):
				self.h['CutFlow'][self.ch].Fill(2)
			else:
				return

		if self.do_filter_cut:
			if self._filter_cut(evt):
				self.h['CutFlow'][self.ch].Fill(3)
			else:
				return

		if self.do_prompt_lepton_cut:
			if self._prompt_lepton_cut(evt):
				self.h['CutFlow'][self.ch].Fill(4)
			else:
				return

		if self.do_invert_prompt_lepton_cut:
			if self._invert_prompt_lepton_cut(evt):
				self.h['CutFlow'][self.ch].Fill(4)
			else:
				return

		if self.do_ndv_cut:
			if self._ndv_cut(evt):
				self.h['CutFlow'][self.ch].Fill(5)
			else:
				return

		# If you've made it here, preselection is passed
		self.passed_preselection_cuts = True

	def _DVSelection(self, evt):

		######################################################################################################
		# DV Selection is any cuts that are done per DV
		# Current cuts include: fiducial vol, ntrack, OS, DVtype, track quality, cosmic veto, mlll, mDV
		######################################################################################################

		# Fill all the histograms with ALL DVs (this could be more that 1 per event). Useful for vertexing efficiency studies.

		self._fill_all_dv_histos(evt)

		# only do the DV selection if the preselction was passed for the event.
		if not self.passed_preselection_cuts:
			return

		# There is an extra bit of logic here since we look at several DVs
		# but only want to fill the cutflow once per event

		# Do we want to use this cut?
		if self.do_fidvol_cut:
			# Does this cut pass?
			if self._fidvol_cut(evt):
				# Has the cutflow already been filled for this event?
				if not self.passed_fidvol_cut:
					self.h['CutFlow'][self.ch].Fill(6)
					self.passed_fidvol_cut = True
			# If this cut doesn't pass, don't continue to check other cuts
			else:
				return

		if self.do_ntrk_cut:
			if self._ntrk_cut(evt):
				if not self.passed_ntrk_cut:
					self.h['CutFlow'][self.ch].Fill(7)
					self.passed_ntrk_cut = True
			else:
				return

		if self.do_opposite_sign_cut or self.do_same_sign_cut:
			if self._charge_cut(evt):
				if not self.passed_charge_cut:
					self.h['CutFlow'][self.ch].Fill(8)
					self.passed_charge_cut = True
			else:
				return

		self._fill_selected_dv_histos(evt, "charge")

		if self.do_dv_type_cut:
			if self._dv_type_cut(evt):	
				if not self.passed_dv_type_cut:
					self.h['CutFlow'][self.ch].Fill(9)
					self.passed_dv_type_cut = True
			else:
				return

		self._fill_selected_dv_histos(evt, "DVtype")

		if self.do_track_quality_cut:
			if self._track_quality_cut(evt):
				if not self.passed_track_quality_cut:
					self.h['CutFlow'][self.ch].Fill(10)
					self.passed_track_quality_cut = True
			else:
				return

		self._fill_selected_dv_histos(evt, "trkqual")

		if self.do_cosmic_veto_cut:
			if self._cosmic_veto_cut(evt):
				if not self.passed_cosmic_veto_cut:
					self.h['CutFlow'][self.ch].Fill(11)
					self.passed_cosmic_veto_cut = True
			else:
				return


		self._fill_selected_dv_histos(evt, "cosmic")

		if self.do_trilepton_mass_cut:
			if self._trilepton_mass_cut(evt):
				if not self.passed_trilepton_mass_cut:
					self.h['CutFlow'][self.ch].Fill(12)
					self.passed_trilepton_mass_cut = True
			else:
				return

		self._fill_selected_dv_histos(evt, "mlll")

		if self.do_dv_mass_cut:
			if self._dv_mass_cut(evt):
				if not self.passed_dv_mass_cut:
					self.h['CutFlow'][self.ch].Fill(13)
					self.passed_dv_mass_cut = True
			else:
				return


		self._fill_selected_dv_histos(evt,"sel")  # Fill all the histograms with only selected DVs. (ie. the ones that pass the full selection)






class ToyAnalysis(Analysis):

	def __init__(self, channel, selections, outputFile, isdata):
		logger.info('Running  Toy Analysis cuts')
		Analysis.__init__(self, channel, selections, outputFile, isdata)

	
		self.add('CutFlow', 17, -0.5, 16.5)
		# Bin labels are 1 greater than histogram bins
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(1, "all")
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(2, "PV")
		if self.do_trigger_cut:
			if self.do_CR == False:
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
		if self.do_dv_type_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(11, "%s DV" % self.dv_type)
		if self.do_dv_mass_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(12, "m_{DV}")
		if self.do_trilepton_mass_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(13, "m_{lll}")
		if self.do_HNL_pt_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(14, "HNL p_{T}")
		if self.do_cosmic_veto_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(15, "cosmic veto")
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(16, "1-tight")
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(17, "2-tight")
		
		# 2D correlation plots after each cut in the DV
		sel_list = ["charge", "DVtype", "mDV","mlll", "HNLpt","sel"]
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
			self.add2D( sel + '_pos_mhnl12_13', 1000, 0, 500, 1000, 0, 500)
			self.add2D( sel + '_pos_mhnl23_12', 1000, 0, 500, 1000, 0, 500)
			self.add2D( sel + '_pos_mhnl13_23', 1000, 0, 500, 1000, 0, 500)
			self.add2D( sel + '_neg_mhnl12_13', 1000, 0, 500, 1000, 0, 500)
			self.add2D( sel + '_neg_mhnl23_12', 1000, 0, 500, 1000, 0, 500)
			self.add2D( sel + '_neg_mhnl13_23', 1000, 0, 500, 1000, 0, 500)

	
	def _fill_histos(self, evt):
		w = evt.weight
		for imu in range(len(evt.tree.muontype[evt.ievt])):
			# print evt.tree.muonpt_wrtSV[evt.ievt][imu]
			self.h["muon_type"][self.ch].Fill(evt.tree.muontype[evt.ievt][imu], w)
			self.h["muon_pt"][self.ch].Fill(evt.tree.muonpt[evt.ievt][imu], w)
			self.h["muon_eta"][self.ch].Fill(evt.tree.muoneta[evt.ievt][imu], w)
			self.h["muon_phi"][self.ch].Fill(evt.tree.muonphi[evt.ievt][imu], w)

		for iel in range(len(evt.tree.elpt[evt.ievt])):
			self.h["el_pt"][self.ch].Fill(evt.tree.elpt[evt.ievt][iel], w)
			self.h["el_eta"][self.ch].Fill(evt.tree.eleta[evt.ievt][iel], w)
			self.h["el_phi"][self.ch].Fill(evt.tree.elphi[evt.ievt][iel], w)


	def _fill_correlation_histos(self, evt, sel):
		w = evt.weight
		# sel refers to the last selection that was applied 

		if self.do_prompt_lepton_cut:

			tracks = helpers.Tracks()
			tracks.getTracks(evt=evt)
			trkVec = tracks.lepVec	

			if tracks.ntracks == 2: 
				Mltt = selections.Mltt(plep=self.plep_sel.plepVec, trks=trkVec)
				Mhnl = selections.Mhnl(evt=evt, plep=self.plep_sel.plepVec, trks =trkVec )
				Mtrans = selections.Mtrans(plep=self.plep_sel.plepVec, trks =trkVec )
				
				# fill 2D mass correlation plots here 
				self.h[sel +'_DVmass_mvis'][self.ch].Fill(evt.tree.dvmass[evt.ievt][evt.idv],Mltt.mltt, w)
				self.h[sel +'_DVmass_mhnl'][self.ch].Fill(evt.tree.dvmass[evt.ievt][evt.idv],Mhnl.mhnl, w)
				self.h[sel +'_DVmass_mtrans'][self.ch].Fill(evt.tree.dvmass[evt.ievt][evt.idv],Mtrans.mtrans, w)
				self.h[sel +'_DVmass_hnlpt'][self.ch].Fill(evt.tree.dvmass[evt.ievt][evt.idv], Mhnl.hnlpt, w)
				self.h[sel +'_mvis_mhnl'][self.ch].Fill(Mltt.mltt,Mhnl.mhnl, w)
				self.h[sel +'_mvis_mtrans'][self.ch].Fill(Mltt.mltt,Mtrans.mtrans, w)
				self.h[sel +'_mvis_hnlpt'][self.ch].Fill(Mltt.mltt,Mhnl.hnlpt, w)
				self.h[sel +'_mhnl_mtrans'][self.ch].Fill(Mhnl.mhnl,Mtrans.mtrans, w)
				self.h[sel +'_mhnl_hnlpt'][self.ch].Fill(Mhnl.mhnl,Mhnl.hnlpt, w)
				self.h[sel +'_mhnl2D'][self.ch].Fill(Mhnl.mhnl,Mhnl.mhnl2, w)
				self.h[sel +'_neg_mhnl12_13'][self.ch].Fill(Mhnl.neg_mhnl12,Mhnl.neg_mhnl13, w)
				self.h[sel +'_neg_mhnl23_12'][self.ch].Fill(Mhnl.neg_mhnl23,Mhnl.neg_mhnl12, w)
				self.h[sel +'_neg_mhnl13_23'][self.ch].Fill(Mhnl.neg_mhnl13,Mhnl.neg_mhnl23, w)
				self.h[sel +'_pos_mhnl12_13'][self.ch].Fill(Mhnl.pos_mhnl12,Mhnl.pos_mhnl13, w)
				self.h[sel +'_pos_mhnl23_12'][self.ch].Fill(Mhnl.pos_mhnl23,Mhnl.pos_mhnl12, w)
				self.h[sel +'_pos_mhnl13_23'][self.ch].Fill(Mhnl.pos_mhnl13,Mhnl.pos_mhnl23, w)
	
	#########################################################################################################################
	# Define new cuts you want to apply here. This will overwrite whatever cuts are defined in the parent analysis class.
	#########################################################################################################################
	def _track_quality_cut_1tight(self, evt):
		track_quality_sel = selections.Trackqual(evt=evt, quality="1-tight")
		return track_quality_sel.passes()

	def _track_quality_cut_2tight(self, evt):
		track_quality_sel = selections.Trackqual(evt=evt, quality="2-tight")
		return track_quality_sel.passes()

	def _dR_cut(self, evt):
		tracks = helpers.Tracks()
		tracks.getTracks(evt=evt)
		trkVec = tracks.lepVec	
		dR = tracks.lepVec[0].DeltaR(tracks.lepVec[1])
		return dR > 0.0

	def _trilepton_mass_cut(self, evt):
		plep_vec = self.plep_sel.plepVec

		muons = helpers.Tracks()
		muons.getMuons(evt=evt)
		muVec = muons.lepVec

		electrons = helpers.Tracks()
		electrons.getElectrons(evt=evt)
		elVec = electrons.lepVec

		mlll_sel = selections.Mlll(dv_type=self.dv_type, plep=plep_vec, dMu=muVec, dEl=elVec)
		return mlll_sel.passes()

	def _dv_mass_cut(self, evt):
		dv_mass_sel = selections.DVmass(evt=evt,dvmasscut=2) # changed the dvmass cut to 2 GeV
		return dv_mass_sel.passes()

	def _HNL_mass_cut(self, evt): # not used 
		tracks = helpers.Tracks()
		tracks.getTracks(evt=evt)

		plep_vec = self.plep_sel.plepVec
		tracks_vec = tracks.lepVec

		mHNL_sel = selections.Mhnl(evt=evt,plep=plep_vec,trks=tracks_vec,hnlmasscut=3)

		return mHNL_sel.passes()


	def _HNL_pt_cut(self, evt):
		tracks = helpers.Tracks()
		tracks.getTracks(evt=evt)

		plep_vec = self.plep_sel.plepVec
		tracks_vec = tracks.lepVec

		mHNL_sel = selections.Mhnl(evt=evt,plep=plep_vec,trks=tracks_vec)

		return (mHNL_sel.hnlpt > 20 and mHNL_sel.hnlpt < 60)


	def _preSelection(self, evt):
		######################################################################################################
		# Preselection are all the cuts that are requied per event
		# Current cuts include: trigger, filter, plepton, DV cut
		######################################################################################################


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

		self._fill_histos(evt)

		self.h['CutFlow'][self.ch].SetBinContent(1, evt.tree.allEvt)


		######################################################################################################
		# Selection code is deisgned so that it will pass the selection only if the cut true or cut is unused
		# ex. passTrigger is true if the trigcut is true OR if trigcut is not used)
		######################################################################################################

		if self._pv_cut(evt): #Check to make sure event has a PV otherwise throw event away (this happens very rarely with data).
			self.h['CutFlow'][self.ch].Fill(1)
		else:
			return

		if self.do_trigger_cut:
			if self._trigger_cut(evt):
				# Fill the plot at the specified bin
				self.h['CutFlow'][self.ch].Fill(2)
			else:
				return

		if self.do_filter_cut:
			if self._filter_cut(evt):
				self.h['CutFlow'][self.ch].Fill(3)
			else:
				return

		if self.do_prompt_lepton_cut:
			if self._prompt_lepton_cut(evt):
				self.h['CutFlow'][self.ch].Fill(4)
			else:
				return

		if self.do_invert_prompt_lepton_cut:
			if self._invert_prompt_lepton_cut(evt):
				self.h['CutFlow'][self.ch].Fill(4)
			else:
				return

		if self.do_ndv_cut:
			if self._ndv_cut(evt):
				self.h['CutFlow'][self.ch].Fill(5)
			else:
				return

		# If you've made it here, preselection is passed
		self.passed_preselection_cuts = True

	def _DVSelection(self, evt):

		######################################################################################################
		# DV Selection is any cuts that are done per DV
		# Current cuts include: fiducial vol, ntrack, OS, DVtype, track quality, cosmic veto, mlll, mDV
		######################################################################################################

		# Fill all the histograms with ALL DVs (this could be more that 1 per event). Useful for vertexing efficiency studies.
		self._fill_all_dv_histos(evt)

		# only do the DV selection if the preselction was passed for the event.
		if not self.passed_preselection_cuts:
			return

		# There is an extra bit of logic here since we look at several DVs
		# but only want to fill the cutflow once per event

		# Do we want to use this cut?
		if self.do_fidvol_cut:
			# Does this cut pass?
			if self._fidvol_cut(evt):
				# Has the cutflow already been filled for this event?
				if not self.passed_fidvol_cut:
					self.h['CutFlow'][self.ch].Fill(6)
					self.passed_fidvol_cut = True
			# If this cut doesn't pass, don't continue to check other cuts
			else:
				return

		if self.do_ntrk_cut:
			if self._ntrk_cut(evt):
				if not self.passed_ntrk_cut:
					self.h['CutFlow'][self.ch].Fill(7)
					self.passed_ntrk_cut = True
			else:
				return

		if self._dR_cut(evt):
			if not self.passed_dR_cut:
				self.h['CutFlow'][self.ch].Fill(8)
				self.passed_dR_cut = True
		else:
			return

		if self.do_opposite_sign_cut or self.do_same_sign_cut:
			if self._charge_cut(evt):
				if not self.passed_charge_cut:
					self.h['CutFlow'][self.ch].Fill(9)
					self.passed_charge_cut = True
			else:
				return

		self._fill_selected_dv_histos(evt, "charge")
		self._fill_correlation_histos(evt, "charge")

		if self.do_dv_type_cut:
			if self._dv_type_cut(evt):
				if not self.passed_dv_type_cut:
					self.h['CutFlow'][self.ch].Fill(10)
					self.passed_dv_type_cut = True

			else:
				return
		self._fill_selected_dv_histos(evt, "DVtype")
		self._fill_correlation_histos(evt, "DVtype")

		if self.do_dv_mass_cut:
			if self._dv_mass_cut(evt):
				if not self.passed_dv_mass_cut:
					self.h['CutFlow'][self.ch].Fill(11)
					self.passed_dv_mass_cut = True
			else:
				return

		self._fill_selected_dv_histos(evt, "mDV")
		self._fill_correlation_histos(evt, "mDV")

		if self.do_trilepton_mass_cut:
			if self._trilepton_mass_cut(evt):
				if not self.passed_trilepton_mass_cut:
					self.h['CutFlow'][self.ch].Fill(12)
					self.passed_trilepton_mass_cut = True
			else:
				return
		self._fill_selected_dv_histos(evt, "mlll")
		self._fill_correlation_histos(evt, "mlll")

		if self.do_HNL_pt_cut:
			if self._HNL_pt_cut(evt):
				if not self.passed_HNL_pt_cut:
					self.h['CutFlow'][self.ch].Fill(13)
					self.passed_HNL_pt_cut = True

			else:
				return

		self._fill_selected_dv_histos(evt, "HNLpt")
		self._fill_correlation_histos(evt, "HNLpt")


		if self.do_cosmic_veto_cut:
			if self._cosmic_veto_cut(evt):
				if not self.passed_cosmic_veto_cut:
					self.h['CutFlow'][self.ch].Fill(14)
					self.passed_cosmic_veto_cut = True
			else:
				return
		self._fill_selected_dv_histos(evt, "cosmic")

		if self._track_quality_cut_1tight(evt):
			if not self.passed_track_1tight_cut:
				self.h['CutFlow'][self.ch].Fill(15)
				self.passed_track_1tight_cut = True
		else:
			return

		self._fill_selected_dv_histos(evt, "tight1")

		if self._track_quality_cut_2tight(evt):
			if not self.passed_track_2tight_cut:
				self.h['CutFlow'][self.ch].Fill(16)
				self.passed_track_2tight_cut = True
		else:
			return


		self._fill_selected_dv_histos(evt,"sel")  # Fill all the histograms with only selected DVs. (ie. the ones that pass the full selection)
		self._fill_correlation_histos(evt, "sel")





class TestAnalysis(Analysis):

	def __init__(self, channel, selections, outputFile, isdata):
		logger.info('Running  Test Analysis cuts (testing lepton quality)')
		Analysis.__init__(self, channel, selections, outputFile, isdata)

	
		self.add('CutFlow', 17, -0.5, 16.5)
		# Bin labels are 1 greater than histogram bins
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(1, "all")
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(2, "PV")
		if self.do_trigger_cut:
			if self.do_CR == False:
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
		if self.do_dv_type_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(11, "%s DV" % self.dv_type)
		if self.do_track_quality_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(12, "{}-lepton DV".format(self.track_quality))
		if self.do_dv_mass_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(13, "m_{DV}")
		if self.do_trilepton_mass_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(14, "m_{lll}")
		if self.do_HNL_pt_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(15, "HNL p_{T}")
		if self.do_cosmic_veto_cut:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(16, "cosmic veto")
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(17, "match")
		
		# 2D correlation plots after each cut in the DV
		sel_list = ["charge", "DVtype", "mDV","mlll", "HNLpt","sel"]
		for sel in sel_list:
			self.add2D( sel + '_ntrk', 11, -5.5, 5.5, 9, -0.5, 8.5)

	
	def _fill_histos(self, evt):
		w = evt.weight
		for imu in range(len(evt.tree.muontype[evt.ievt])):
			# print evt.tree.muonpt_wrtSV[evt.ievt][imu]
			self.h["muon_type"][self.ch].Fill(evt.tree.muontype[evt.ievt][imu], w)
			self.h["muon_pt"][self.ch].Fill(evt.tree.muonpt[evt.ievt][imu], w)
			self.h["muon_eta"][self.ch].Fill(evt.tree.muoneta[evt.ievt][imu], w)
			self.h["muon_phi"][self.ch].Fill(evt.tree.muonphi[evt.ievt][imu], w)

		for iel in range(len(evt.tree.elpt[evt.ievt])):
			self.h["el_pt"][self.ch].Fill(evt.tree.elpt[evt.ievt][iel], w)
			self.h["el_eta"][self.ch].Fill(evt.tree.eleta[evt.ievt][iel], w)
			self.h["el_phi"][self.ch].Fill(evt.tree.elphi[evt.ievt][iel], w)

		if evt.tree.isData == True: 
			pass
		else: 
			truthInfo = helpers.Truth()
			truthInfo.getTruthParticles(evt)
			self.h["truth_all_W_pt"][self.ch].Fill(truthInfo.W_vec.Pt(), w)
			self.h["truth_all_W_eta"][self.ch].Fill(truthInfo.W_vec.Eta(), w)
			self.h["truth_all_W_phi"][self.ch].Fill(truthInfo.W_vec.Phi(), w)
			self.h["truth_all_W_mass"][self.ch].Fill(truthInfo.W_vec.M(), w)
			self.h["truth_all_HNL_pt"][self.ch].Fill(truthInfo.HNL_vec.Pt(), w)
			self.h["truth_all_HNL_eta"][self.ch].Fill(truthInfo.HNL_vec.Eta(), w)
			self.h["truth_all_HNL_phi"][self.ch].Fill(truthInfo.HNL_vec.Phi(), w)
			self.h["truth_all_HNL_mass"][self.ch].Fill(truthInfo.HNL_vec.M(), w)

			self.h["truth_all_mHNLcalc"][self.ch].Fill(truthInfo.mhnl, w)

			self.h["truth_all_DV_r"][self.ch].Fill(truthInfo.truth_dvr, w)
			self.h["truth_all_DV_x"][self.ch].Fill(truthInfo.truth_dvx, w)
			self.h["truth_all_DV_y"][self.ch].Fill(truthInfo.truth_dvy, w)
			self.h["truth_all_DV_z"][self.ch].Fill(truthInfo.truth_dvz, w)
			self.h["truth_all_plep_pt"][self.ch].Fill(truthInfo.plep_vec.Pt(), w)
			self.h["truth_all_plep_eta"][self.ch].Fill(truthInfo.plep_vec.Eta(), w)
			self.h["truth_all_plep_phi"][self.ch].Fill(truthInfo.plep_vec.Phi(), w)
			self.h["truth_all_plep_mass"][self.ch].Fill(truthInfo.plep_vec.M(), w)
			self.h["truth_all_lep1_trk_pt"][self.ch].Fill(truthInfo.trkVec[0].Pt(), w)
			self.h["truth_all_lep1_trk_eta"][self.ch].Fill(truthInfo.trkVec[0].Eta(), w)
			self.h["truth_all_lep1_trk_phi"][self.ch].Fill(truthInfo.trkVec[0].Phi(), w)
			self.h["truth_all_lep2_trk_pt"][self.ch].Fill(truthInfo.trkVec[1].Pt(), w)
			self.h["truth_all_lep2_trk_eta"][self.ch].Fill(truthInfo.trkVec[1].Eta(), w)
			self.h["truth_all_lep2_trk_phi"][self.ch].Fill(truthInfo.trkVec[1].Phi(), w)

			for itrk in xrange(2): 
				self.h["truth_all_DV_trk_pt"][self.ch].Fill(truthInfo.trkVec[itrk].Pt(), w)
				self.h["truth_all_DV_trk_eta"][self.ch].Fill(truthInfo.trkVec[itrk].Eta(), w)
				self.h["truth_all_DV_trk_phi"][self.ch].Fill(truthInfo.trkVec[itrk].Phi(), w)

	
	#########################################################################################################################
	# Define new cuts you want to apply here. This will overwrite whatever cuts are defined in the parent analysis class.
	#########################################################################################################################

	def _dR_cut(self, evt):
		tracks = helpers.Tracks()
		tracks.getTracks(evt=evt)
		trkVec = tracks.lepVec	
		dR = tracks.lepVec[0].DeltaR(tracks.lepVec[1])
		return dR > 0.0

	def _trilepton_mass_cut(self, evt):
		plep_vec = self.plep_sel.plepVec

		muons = helpers.Tracks()
		muons.getMuons(evt=evt)
		muVec = muons.lepVec

		electrons = helpers.Tracks()
		electrons.getElectrons(evt=evt)
		elVec = electrons.lepVec

		mlll_sel = selections.Mlll(dv_type=self.dv_type, plep=plep_vec, dMu=muVec, dEl=elVec)
		return mlll_sel.passes()

	def _dv_mass_cut(self, evt):
		dv_mass_sel = selections.DVmass(evt=evt,dvmasscut=2) # changed the dvmass cut to 2 GeV
		return dv_mass_sel.passes()

	def _HNL_mass_cut(self, evt): # not used 
		tracks = helpers.Tracks()
		tracks.getTracks(evt=evt)

		plep_vec = self.plep_sel.plepVec
		tracks_vec = tracks.lepVec

		mHNL_sel = selections.Mhnl(evt=evt,plep=plep_vec,trks=tracks_vec,hnlmasscut=3)

		return mHNL_sel.passes()


	def _HNL_pt_cut(self, evt):
		tracks = helpers.Tracks()
		tracks.getTracks(evt=evt)

		plep_vec = self.plep_sel.plepVec
		tracks_vec = tracks.lepVec

		mHNL_sel = selections.Mhnl(evt=evt,plep=plep_vec,trks=tracks_vec)

		return (mHNL_sel.hnlpt > 20 and mHNL_sel.hnlpt < 60)


	def _preSelection(self, evt):
		######################################################################################################
		# Preselection are all the cuts that are requied per event
		# Current cuts include: trigger, filter, plepton, DV cut
		######################################################################################################


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
                self.passed_track_quality_cut = False
		self.passed_cosmic_veto_cut = False
		self.passed_trilepton_mass_cut = False
		self.passed_dv_mass_cut = False
		self.passed_HNL_mass_cut = False
		self.passed_HNL_pt_cut = False

		self._fill_histos(evt)

		self.h['CutFlow'][self.ch].SetBinContent(1, evt.tree.allEvt)


		######################################################################################################
		# Selection code is deisgned so that it will pass the selection only if the cut true or cut is unused
		# ex. passTrigger is true if the trigcut is true OR if trigcut is not used)
		######################################################################################################

		if self._pv_cut(evt): #Check to make sure event has a PV otherwise throw event away (this happens very rarely with data).
			self.h['CutFlow'][self.ch].Fill(1)
		else:
			return

		if self.do_trigger_cut:
			if self._trigger_cut(evt):
				# Fill the plot at the specified bin
				self.h['CutFlow'][self.ch].Fill(2)
			else:
				return

		if self.do_filter_cut:
			if self._filter_cut(evt):
				self.h['CutFlow'][self.ch].Fill(3)
			else:
				return

		if self.do_prompt_lepton_cut:
			if self._prompt_lepton_cut(evt):
				self.h['CutFlow'][self.ch].Fill(4)
			else:
				return

		if self.do_invert_prompt_lepton_cut:
			if self._invert_prompt_lepton_cut(evt):
				self.h['CutFlow'][self.ch].Fill(4)
			else:
				return

		if self.do_ndv_cut:
			if self._ndv_cut(evt):
				self.h['CutFlow'][self.ch].Fill(5)
			else:
				return

		# If you've made it here, preselection is passed
		self.passed_preselection_cuts = True

	def _DVSelection(self, evt):

		######################################################################################################
		# DV Selection is any cuts that are done per DV
		# Current cuts include: fiducial vol, ntrack, OS, DVtype, track quality, cosmic veto, mlll, mDV
		######################################################################################################

		# Fill all the histograms with ALL DVs (this could be more that 1 per event). Useful for vertexing efficiency studies.
		self._fill_all_dv_histos(evt)

		# only do the DV selection if the preselction was passed for the event.
		if not self.passed_preselection_cuts:
			return

		# There is an extra bit of logic here since we look at several DVs
		# but only want to fill the cutflow once per event

		# Do we want to use this cut?
		if self.do_fidvol_cut:
			# Does this cut pass?
			if self._fidvol_cut(evt):
				# Has the cutflow already been filled for this event?
				if not self.passed_fidvol_cut:
					self.h['CutFlow'][self.ch].Fill(6)
					self.passed_fidvol_cut = True
			# If this cut doesn't pass, don't continue to check other cuts
			else:
				return

		if self.do_ntrk_cut:
			if self._ntrk_cut(evt):
				if not self.passed_ntrk_cut:
					self.h['CutFlow'][self.ch].Fill(7)
					self.passed_ntrk_cut = True
			else:
				return

		if self._dR_cut(evt):
			if not self.passed_dR_cut:
				self.h['CutFlow'][self.ch].Fill(8)
				self.passed_dR_cut = True
		else:
			return

		if self.do_opposite_sign_cut or self.do_same_sign_cut:
			if self._charge_cut(evt):
				if not self.passed_charge_cut:
					self.h['CutFlow'][self.ch].Fill(9)
					self.passed_charge_cut = True
			else:
				return

		self._fill_selected_dv_histos(evt, "charge")

		if self.do_dv_type_cut:
			if self._dv_type_cut(evt):
				if not self.passed_dv_type_cut:
					self.h['CutFlow'][self.ch].Fill(10)
					self.passed_dv_type_cut = True

			else:
				return
		self._fill_selected_dv_histos(evt, "DVtype")

		if self.do_track_quality_cut:
			if self._track_quality_cut(evt):
				if not self.passed_track_quality_cut:
					self.h['CutFlow'][self.ch].Fill(11)
					self.passed_track_quality_cut = True
			else:
				return

		self._fill_selected_dv_histos(evt, "trkqual")


		if self.do_dv_mass_cut:
			if self._dv_mass_cut(evt):
				if not self.passed_dv_mass_cut:
					self.h['CutFlow'][self.ch].Fill(12)
					self.passed_dv_mass_cut = True
			else:
				return

		self._fill_selected_dv_histos(evt, "mDV")

		if self.do_trilepton_mass_cut:
			if self._trilepton_mass_cut(evt):
				if not self.passed_trilepton_mass_cut:
					self.h['CutFlow'][self.ch].Fill(13)
					self.passed_trilepton_mass_cut = True
			else:
				return
		self._fill_selected_dv_histos(evt, "mlll")

		if self.do_HNL_pt_cut:
			if self._HNL_pt_cut(evt):
				if not self.passed_HNL_pt_cut:
					self.h['CutFlow'][self.ch].Fill(14)
					self.passed_HNL_pt_cut = True

			else:
				return

		self._fill_selected_dv_histos(evt, "HNLpt")


		if self.do_cosmic_veto_cut:
			if self._cosmic_veto_cut(evt):
				if not self.passed_cosmic_veto_cut:
					self.h['CutFlow'][self.ch].Fill(15)
					self.passed_cosmic_veto_cut = True
			else:
				return
		self._fill_selected_dv_histos(evt, "cosmic")


#                print "hnl? ",  evt.tree.linkedtruthParentPdgId[evt.ievt][evt.idv]
 #               print "hnl? ",  evt.tree.linkedtruthscore[evt.ievt][evt.idv]
                if not evt.tree.isData:
                        if evt.tree.linkedtruthscore[evt.ievt][evt.idv] > 0.75 and (evt.tree.linkedtruthParentPdgId[evt.ievt][evt.idv] == 50 or evt.tree.linkedtruthParentPdgId[evt.ievt][evt.idv] == -50):
                                self.h['CutFlow'][self.ch].Fill(16)
                                self._fill_selected_dv_histos(evt,"match")  # Fill all the histograms with only selected DVs. (ie. the ones that pass the full selection)



		self._fill_selected_dv_histos(evt,"sel")  # Fill all the histograms with only selected DVs. (ie. the ones that pass the full selection)






