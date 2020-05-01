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
logger = helpers.getLogger('dHNLAnalysis.analysis',level = logging.WARNING)





UNLOCKED = 0
SELECTION_LOCKED = 1
FILL_LOCKED = 2


class Analysis(object):
	ch = ''
	h = {}  # map of histogram names to map of systematics to histograms
	blinded = True
	mapSel = {}

	def __init__(self, channel, selections, outputFile):
		# print channel
		# print outputFile

		self.sel = selections
		self._outputFile = outputFile
		self.fi = ROOT.TFile.Open(outputFile, 'update')
		self.ch = channel
		self.histSuffixes = [self.ch]
		self.h = {}
		# make histograms (common for all channels)
		self.add('CutFlow', 15, -0.5, 14.5)
		self.add2D('charge_ntrk', 11, -5.5, 5.5, 9, -0.5, 8.5)

		self.add2D( 'charge_DVmass_mvis', 1000, 0, 500, 1000, 0, 500)
		self.add2D( 'charge_DVmass_mhnl', 1000, 0, 500, 1005, -5, 500)
		self.add2D( 'charge_DVmass_mtrans', 1000, 0, 500, 1000, 0, 500)
		self.add2D( 'charge_mvis_mhnl', 1000, 0, 500, 1005, -5, 500)
		self.add2D( 'charge_mvis_mtrans', 1000, 0, 500, 1000, 0, 500)
		self.add2D( 'charge_mhnl_mtrans', 1005, -5, 500, 1000, 0, 500)
		self._locked = UNLOCKED

		self.observables = [observable.registered(self) for observable in observables.ObservableList if ((observable.only is None) or any(only in self.sel for only in observable.only))]
		for observable in self.observables:
			if 'hist' in observable.do:
				if type(observable.binning) == tuple:
					self.add(observable.name, *observable.binning)
				else:
					self.addVar(observable.name, observable.binning)
		# 	 if isdata and observable.need_truth: # need to add something like this for truth variables
		#            continue

		# Parse input cuts

		# trigger cut
		if 'alltriggers' in self.sel:
			self.trigger = 'all'
			self.do_trigger_cut = True
		if 'invertedtriggers' in self.sel:
			self.trigger = 'inverted'
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

		# nDV cut
		self.do_ndv_cut = ('nDV' in self.sel)
		if not self.do_ndv_cut: logger.warn('You did not add nDV cut. Skipping nDV selection.')

		# alpha cut
		self.do_alpha_cut = 'alpha' in self.sel

		# mass window cut
		self.do_mass_window_cut = 'mass_window' in self.sel

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
		else:
			logger.warn('You did not specify a DV track quality for this channel. Skipping DV track quality selection.')
			self.do_track_quality_cut = False

		# cosmic veto cut
		self.do_cosmic_veto_cut = 'cosmicveto' in self.sel
		if not self.do_cosmic_veto_cut: logger.warn('You did not add a cosmic veto cut for this channel. Skipping cosmic veto selection.')

		# tri-lepton mass cut
		self.do_trilepton_mass_cut = 'mlll' in self.sel
		if not self.do_trilepton_mass_cut: logger.warn('You did not add a mlll cut for this channel. Skipping tri-lepton mass selection.')

		# DV mass cut
		self.do_dv_mass_cut = 'DVmass' in self.sel
		if not self.do_dv_mass_cut: logger.warn('You did not add a DVmass cut for this channel. Skipping displaced vertex mass selection.')

		#HNL mass cut
		self.do_HNL_mass_cut = 'HNLmass' in self.sel
		# if not self.do_HNL_mass_cut: logger.warn('You did not add an HNLmass cut for this channel. Skipping HNL mass selection.')


		self.check_input_consistency()

		self.set_cutflow_labels()

	def check_input_consistency(self):
		if self.do_trilepton_mass_cut and not self.do_prompt_lepton_cut:
			logger.error("You cannot calculate mlll without selecting a prompt lepton!")
			sys.exit(1)  # abort because of error

		if self.do_opposite_sign_cut and self.do_same_sign_cut:
			logger.error("These cuts are mutually exclusive. You will get zero events!")
			sys.exit(1)  # abort because of error

	def set_cutflow_labels(self):
		# Bin labels are 1 greater than histogram bins
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(1, "all")
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(2, "PV")
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(3, "Kshort mass")
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(4, "alpha")

		# if self.do_trigger_cut:
		# 	self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(3, "trigger")
		# if self.do_filter_cut:
		# 	self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(4, "%s" % self.filter_type)
		# if self.do_prompt_lepton_cut:
		# 	self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(5, "tight prompt %s" % self.plep)
		# if self.do_ndv_cut:
		# 	self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(6, "DV")
		# if self.do_fidvol_cut:
		# 	self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(7, "fiducial")
		# if self.do_ntrk_cut:
		# 	self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(8, "%s-track DV" % self.ntrk)
		# if self.do_opposite_sign_cut:
		# 	self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(9, "OS DV")
		# if self.do_same_sign_cut:
		# 	self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(9, "SS DV")
		# if self.do_dv_type_cut:
		# 	self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(10, "%s DV" % self.dv_type)
		# if self.do_track_quality_cut:
		# 	self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(11, "{}-lepton DV".format(self.track_quality))
		# if self.do_cosmic_veto_cut:
		# 	self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(12, "cosmic veto")
		# if self.do_trilepton_mass_cut:
		# 	self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(13, "m_{lll}")
		# if self.do_dv_mass_cut:
		# 	self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(14, "m_{DV}")
		# if self.do_HNL_mass_cut:
		# 	self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(15, "m_{HNL}")

	def unlock(self):
		self._locked = UNLOCKED

	def add(self, hName, nBins, xLow, xHigh):
		self.h[hName] = {}
		# self.fi.cd()
		# self.h[hName] = ROOT.TH1D(hName+'_'+self.ch, "", nBins, xLow, xHigh)
		# self.h[hName].Sumw2()
		# self.h[hName].SetDirectory(0)
		for s in self.histSuffixes:
			# print 'adding histogram with name ', hName+self.ch+s
			self.h[hName][s] = ROOT.TH1D(hName+"_"+s, "", nBins, xLow, xHigh)
			self.h[hName][s].Sumw2()
			self.h[hName][s].SetDirectory(0)

	def addVar(self, hName, nBinsList):
		ar = np.array("d", nBinsList)
		# self.fi.cd()
		self.h[hName] = {}
		self.h[hName] = ROOT.TH1D(hName+'_'+self.ch, "", len(nBinsList) - 1, ar)
		self.h[hName].Sumw2()
		self.h[hName].SetDirectory(0)
		# for s in self.histSuffixes:
			# self.h[hName][s] = ROOT.TH1D(hName+self.ch+s, "", len(nBinsList) - 1, ar)
			# self.h[hName][s].Sumw2()
			# self.h[hName][s].SetDirectory(0)

	def add2D(self, hName, nBins, xLow, xHigh, nBinsY, yLow, yHigh):
		self.h[hName] = {}
		# self.h[hName] = ROOT.TH2D(hName+'_'+self.ch, "", nBins, xLow, xHigh, nBinsY, yLow, yHigh)
		# self.h[hName].Sumw2()
		# self.h[hName].SetDirectory(0)
		# self.fi.cd()
		for s in self.histSuffixes:
			# print 'adding histogram with name ', hName+self.ch+s
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

	def _filter_cut(self, evt):
		filter_sel = selections.Filter(evt=evt, filter_type=self.filter_type)
		return filter_sel.passes()

	def _prompt_lepton_cut(self, evt):
		self.plep_sel = selections.Plepton(evt=evt, lepton=self.plep)

		# Add to histogram all prompt leptons that pass selection.
		# If _prompt_lepton_cut() is run after trigger and filter cut then those cuts will also be applied.
		self.h["all_plep_pt"][self.ch].Fill(self.plep_sel.plepVec.Pt())
		self.h["all_plep_eta"][self.ch].Fill(self.plep_sel.plepVec.Eta())
		self.h["all_plep_phi"][self.ch].Fill(self.plep_sel.plepVec.Phi())
		self.h["all_plep_d0"][self.ch].Fill(self.plep_sel.plepd0)
		self.h["all_plep_z0"][self.ch].Fill(self.plep_sel.plepz0)
		return self.plep_sel.passes()

	def _ndv_cut(self, evt):
		dv_sel = selections.nDV(evt=evt)
		return dv_sel.passes()

	def _fidvol_cut(self, evt):
		fidvol_sel = selections.DVradius(evt=evt)
		return fidvol_sel.passes()
		
	def _alpha_cut(self, evt):
		alpha_sel = selections.Alpha(evt=evt)
		return alpha_sel.passes()

	def _ntrk_cut(self, evt):
		ntracks_sel = selections.DVntracks(evt=evt, ntrk=self.ntrk)
		return ntracks_sel.passes()

	def _charge_cut(self, evt):
		sign_pair = "SS" if self.do_same_sign_cut else "OS"
		charge_sel = selections.ChargeDV(evt=evt, sel=sign_pair)
		return charge_sel.passes()

	def _dv_type_cut(self, evt, dv_type):
		dv_sel = selections.DVtype(evt=evt, dv_type=dv_type)
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
		dv_mass_sel = selections.DVmass(evt=evt)
		return dv_mass_sel.passes()

	def _mass_window_cut(self, evt):
		dv_mass_window_sel = selections.DV_mass_window(evt=evt)
		return dv_mass_window_sel.passes()

	def _HNL_mass_cut(self, evt):
		tracks = helpers.Tracks()
		tracks.getTracks(evt=evt)

		plep_vec = self.plep_sel.plepVec
		tracks_vec = tracks.lepVec

		mHNL_sel = selections.Mhnl(evt=evt,plep=plep_vec,trks=tracks_vec)

		return mHNL_sel.passes()



	def _fill_histos(self, evt):
		for imu in range(len(evt.tree.muontype[evt.ievt])):
			self.h["muon_type"][self.ch].Fill(evt.tree.muontype[evt.ievt][imu])
			self.h["muon_pt"][self.ch].Fill(evt.tree.muonpt[evt.ievt][imu])
			self.h["muon_eta"][self.ch].Fill(evt.tree.muoneta[evt.ievt][imu])
			self.h["muon_phi"][self.ch].Fill(evt.tree.muonphi[evt.ievt][imu])

		for iel in range(len(evt.tree.elpt[evt.ievt])):
			self.h["el_pt"][self.ch].Fill(evt.tree.elpt[evt.ievt][iel])
			self.h["el_eta"][self.ch].Fill(evt.tree.eleta[evt.ievt][iel])
			self.h["el_phi"][self.ch].Fill(evt.tree.elphi[evt.ievt][iel])

		# fill truth
		# should add a flag like this to prevent filling
		# if isdata and observable.need_truth: # need to get the "needs truth variable"
			# return
 		#else:
			# BUG with the truth values in ntuple need to investigate this -DT
			# self.h["truth_DV_x"].Fill(evt.tree.truth_dvx[evt.ievt])
			# self.h["truth_DV_y"].Fill(evt.tree.truth_dvy[evt.ievt])
			# self.h["truth_DV_z"].Fill(evt.tree.truth_dvz[evt.ievt])
			# self.h["truth_DV_r"].Fill(evt.tree.truth_dvr[evt.ievt])
			# self.h["truth_DV_mass"].Fill(evt.tree.truth_dvmass[evt.ievt])
			# self.h["truth_DV_pt"].Fill(evt.tree.truth_dvpt[evt.ievt])
			# self.h["truth_DV_eta"].Fill(evt.tree.truth_dveta[evt.ievt])
			# self.h["truth_DV_phi"].Fill(evt.tree.truth_dvphi[evt.ievt])

	def _fill_all_dv_histos(self, evt):

		self.h["charge_ntrk"][self.ch].Fill(evt.tree.dvcharge[evt.ievt][evt.idv], evt.tree.dvntrk[evt.ievt][evt.idv])
		ntracks = len(evt.tree.trackd0[evt.ievt][evt.idv])
		for itrk in range(ntracks):  # loop over tracks
			self.h["all_DV_trk_pt"][self.ch].Fill(evt.tree.trackpt[evt.ievt][evt.idv][itrk])
			self.h["all_DV_trk_eta"][self.ch].Fill(evt.tree.tracketa[evt.ievt][evt.idv][itrk])
			self.h["all_DV_trk_phi"][self.ch].Fill(evt.tree.trackphi[evt.ievt][evt.idv][itrk])
			self.h["all_DV_trk_d0"][self.ch].Fill(evt.tree.trackd0[evt.ievt][evt.idv][itrk])
			self.h["all_DV_trk_z0"][self.ch].Fill(evt.tree.trackz0[evt.ievt][evt.idv][itrk])
			self.h["all_DV_trk_charge"][self.ch].Fill(evt.tree.trackcharge[evt.ievt][evt.idv][itrk])
			self.h["all_DV_trk_chi2"][self.ch].Fill(evt.tree.trackchi2[evt.ievt][evt.idv][itrk])

		self.h["all_DV_num_trks"][self.ch].Fill(evt.tree.dvntrk[evt.ievt][evt.idv])
		self.h["all_DV_x"][self.ch].Fill(evt.tree.dvx[evt.ievt][evt.idv])
		self.h["all_DV_y"][self.ch].Fill(evt.tree.dvy[evt.ievt][evt.idv])
		self.h["all_DV_z"][self.ch].Fill(evt.tree.dvz[evt.ievt][evt.idv])
		self.h["all_DV_r"][self.ch].Fill(evt.tree.dvr[evt.ievt][evt.idv])
		self.h["all_DV_distFromPV"][self.ch].Fill(evt.tree.dvdistFromPV[evt.ievt][evt.idv])
		self.h["all_DV_mass"][self.ch].Fill(evt.tree.dvmass[evt.ievt][evt.idv])
		self.h["all_DV_pt"][self.ch].Fill(evt.tree.dvpt[evt.ievt][evt.idv])
		self.h["all_DV_eta"][self.ch].Fill(evt.tree.dveta[evt.ievt][evt.idv])
		self.h["all_DV_phi"][self.ch].Fill(evt.tree.dvphi[evt.ievt][evt.idv])
		self.h["all_DV_minOpAng"][self.ch].Fill(evt.tree.dvminOpAng[evt.ievt][evt.idv])
		self.h["all_DV_maxOpAng"][self.ch].Fill(evt.tree.dvmaxOpAng[evt.ievt][evt.idv])
		self.h["all_DV_charge"][self.ch].Fill(evt.tree.dvcharge[evt.ievt][evt.idv])
		self.h["all_DV_chi2"][self.ch].Fill(evt.tree.dvchi2[evt.ievt][evt.idv])

		self.h["all_DV_alpha"][self.ch].Fill(selections.Alpha(evt=evt).alpha)


	def _fill_selected_dv_histos(self, evt, sel):
		if self._locked < FILL_LOCKED:
			# these are the histograms you only want to fill ONCE per DV
			# sel refers to the last selection that was applied 

			if self.do_prompt_lepton_cut:
				plep_vec = self.plep_sel.plepVec
				plepd0 = self.plep_sel.plepd0
				plepz0 = self.plep_sel.plepz0

				self.h[sel + "_plep_pt"][self.ch].Fill(plep_vec.Pt())
				self.h[sel + "_plep_eta"][self.ch].Fill(plep_vec.Eta())
				self.h[sel + "_plep_phi"][self.ch].Fill(plep_vec.Phi())
				self.h[sel + "_plep_d0"][self.ch].Fill(plepd0)
				self.h[sel + "_plep_z0"][self.ch].Fill(plepz0)
				
				tracks = helpers.Tracks()
				tracks.getTracks(evt=evt)
				trkVec = tracks.lepVec	

	

				if tracks.ntracks == 2: 
					# mltt = selections.Mltt(plep=plep_vec, trks=trkVec)
					Mhnl = selections.Mhnl(evt=evt, plep=plep_vec, trks =trkVec )
					Mtrans = selections.Mtrans(plep=plep_vec, trks =trkVec )

					self.h[sel + "_mvis"][self.ch].Fill(Mhnl.mvis)
					self.h[sel + "_HNLm"][self.ch].Fill(Mhnl.mhnl)
					self.h[sel + "_HNLpt"][self.ch].Fill(Mhnl.hnlpt)
					self.h[sel + "_HNLeta"][self.ch].Fill(Mhnl.hnleta)
					self.h[sel + "_HNLphi"][self.ch].Fill(Mhnl.hnlphi)
					self.h[sel + "_mtrans"][self.ch].Fill(Mtrans.mtrans)
					self.h[sel + "_mtrans_rot"][self.ch].Fill(Mhnl.mtrans_rot)
					
					if sel == "charge": # fill 2D mass correlation plots here 
						self.h[sel +'_DVmass_mvis'][self.ch].Fill(evt.tree.dvmass[evt.ievt][evt.idv],Mhnl.mvis)
						self.h[sel +'_DVmass_mhnl'][self.ch].Fill(evt.tree.dvmass[evt.ievt][evt.idv],Mhnl.mhnl)
						self.h[sel +'_DVmass_mtrans'][self.ch].Fill(evt.tree.dvmass[evt.ievt][evt.idv],Mtrans.mtrans)
						self.h[sel +'_mvis_mhnl'][self.ch].Fill(Mhnl.mvis,Mhnl.mhnl)
						self.h[sel +'_mvis_mtrans'][self.ch].Fill(Mhnl.mvis,Mtrans.mtrans)
						self.h[sel +'_mhnl_mtrans'][self.ch].Fill(Mhnl.mhnl,Mtrans.mtrans)
					
					deta = abs(tracks.eta[0] - tracks.eta[1])
					dphi = abs(tracks.phi[0] - tracks.phi[1])
					dpt = abs(tracks.pt[0] - tracks.pt[1])
					dR =  np.sqrt(deta**2 + dphi**2)

					self.h[sel + "_DV_trk_deta"][self.ch].Fill(deta)
					self.h[sel + "_DV_trk_dphi"][self.ch].Fill(dphi)
					self.h[sel + "_DV_trk_dpt"][self.ch].Fill(dpt)
					self.h[sel + "_DV_trk_dR"][self.ch].Fill(dR)

					# print evt.tree.dvmass[evt.ievt][evt.idv]/dR

					self.h[sel + "_DV_redmass"][self.ch].Fill(evt.tree.dvmass[evt.ievt][evt.idv]/dR)
					self.h[sel + "_DV_redmassvis"][self.ch].Fill(Mhnl.mvis/dR)
					self.h[sel + "_DV_redmassHNL"][self.ch].Fill(Mhnl.mhnl/dR)


			ntracks = len(evt.tree.trackd0[evt.ievt][evt.idv])
			for itrk in range(ntracks):  # loop over tracks
				self.h[sel + "_DV_trk_pt"][self.ch].Fill(evt.tree.trackpt[evt.ievt][evt.idv][itrk])
				self.h[sel + "_DV_trk_eta"][self.ch].Fill(evt.tree.tracketa[evt.ievt][evt.idv][itrk])
				self.h[sel + "_DV_trk_phi"][self.ch].Fill(evt.tree.trackphi[evt.ievt][evt.idv][itrk])
				self.h[sel + "_DV_trk_d0"][self.ch].Fill(evt.tree.trackd0[evt.ievt][evt.idv][itrk])
				self.h[sel + "_DV_trk_z0"][self.ch].Fill(evt.tree.trackz0[evt.ievt][evt.idv][itrk])
				self.h[sel + "_DV_trk_charge"][self.ch].Fill(evt.tree.trackcharge[evt.ievt][evt.idv][itrk])
				self.h[sel + "_DV_trk_chi2"][self.ch].Fill(evt.tree.trackchi2[evt.ievt][evt.idv][itrk])

			self.h[sel + "_DV_num_trks"][self.ch].Fill(evt.tree.dvntrk[evt.ievt][evt.idv])
			self.h[sel + "_DV_x"][self.ch].Fill(evt.tree.dvx[evt.ievt][evt.idv])
			self.h[sel + "_DV_y"][self.ch].Fill(evt.tree.dvy[evt.ievt][evt.idv])
			self.h[sel + "_DV_z"][self.ch].Fill(evt.tree.dvz[evt.ievt][evt.idv])
			self.h[sel + "_DV_r"][self.ch].Fill(evt.tree.dvr[evt.ievt][evt.idv])
			self.h[sel + "_DV_distFromPV"][self.ch].Fill(evt.tree.dvdistFromPV[evt.ievt][evt.idv])
			self.h[sel + "_DV_mass"][self.ch].Fill(evt.tree.dvmass[evt.ievt][evt.idv])
			self.h[sel + "_DV_pt"][self.ch].Fill(evt.tree.dvpt[evt.ievt][evt.idv])
			self.h[sel + "_DV_eta"][self.ch].Fill(evt.tree.dveta[evt.ievt][evt.idv])
			self.h[sel + "_DV_phi"][self.ch].Fill(evt.tree.dvphi[evt.ievt][evt.idv])
			self.h[sel + "_DV_minOpAng"][self.ch].Fill(evt.tree.dvminOpAng[evt.ievt][evt.idv])
			self.h[sel + "_DV_maxOpAng"][self.ch].Fill(evt.tree.dvmaxOpAng[evt.ievt][evt.idv])
			self.h[sel + "_DV_charge"][self.ch].Fill(evt.tree.dvcharge[evt.ievt][evt.idv])
			self.h[sel + "_DV_chi2"][self.ch].Fill(evt.tree.dvchi2[evt.ievt][evt.idv])
		
			self.h[sel + "_DV_alpha"][self.ch].Fill(selections.Alpha(evt=evt).alpha)

			

			# if sel == "sel": 
			# 	self._locked = FILL_LOCKED  # this only becomes unlocked after the event loop finishes in makeHistograms so you can only fill one DV from each event.

class KShort(Analysis):
	def _init__(self, channel, selections, outputFile):
		Analysis.__init__(self, channel, outputFile)
		self.passed_alpha_cut = False
		self.passed_mass_window_cut = False
		# bin labels are one greater than histogram bins
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(3, "Kshort mass")
		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(4, "alpha")

	def _preSelection(self, evt):
		self.passed_preselection_cuts = False
		self.passed_alpha_cut = False
		self.passed_mass_window_cut = False
		self._fill_histos(evt)
		self.h['CutFlow'][self.ch].Fill(0)

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

		if self.do_mass_window_cut:
			if self._mass_window_cut(evt):
				self.h['CutFlow'][self.ch].Fill(2)
				self.passed_mass_window_cut = True
			else:
				return
		self._fill_selected_dv_histos(evt, "mass")

		if self.do_alpha_cut:
			if self._alpha_cut(evt):
				self.h['CutFlow'][self.ch].Fill(3)
				self.passed_alpha_cut = True
			else:
				return
		self._fill_selected_dv_histos(evt, "alpha")

		self._fill_selected_dv_histos(evt, "sel")


class WmuHNL(Analysis):
	isdata = False  # not currently used...

	def _init__(self, channel, selections, outputFile):
		Analysis.__init__(self, channel, outputFile)

		# make histograms specific to this channel
		# self.add2D()
		# self.add()


	def _preSelection(self, evt):
		######################################################################################################
		# Preselection are all the cuts that are requied per event
		# Current cuts include: trigger, filter, plepton, DV cut
		######################################################################################################

		# if self.ch not in self.mapSel:
		# 	logger.warn('The selected channel '{}' is not registered. The events will be processed anyway without any further constraint.'.format(self.ch))
		# 	self.mapSel[self.ch] = [self.ch]

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
			if self._dv_type_cut(evt, self.dv_type):
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

		if self.do_HNL_mass_cut:
			if self._HNL_mass_cut(evt):
				if not self.passed_HNL_mass_cut:
					self.h['CutFlow'][self.ch].Fill(14)
					self.passed_HNL_mass_cut = True

			else:
				return

		self._fill_selected_dv_histos(evt,"sel")  # Fill all the histograms with only selected DVs. (ie. the ones that pass the full selection)

######################################################################################################################
# An example of a new class. Here you could add any new cuts you want without disturbing the main analysis cuts.
# To use your new class update the class name in analysis.py (e.g. anaClass = getattr(analysis, "new_class") )
######################################################################################################################

# class new_class(Analysis):

# 	isdata = False # not currently used...

# 	def _init__(self, channel, selections, outputFile):
# 		Analysis.__init__(self, channel, outputFile)
# 		#make histograms specfic to this channel
# 		# self.add2D()
# 		# self.add()

# 	def _DVtypeCut(self, evt):

# 		DV_sel = selections.DVtype(evt= evt,dv_type="ee")
# 		print "doing the ee DV type selection"
# 		return DV_sel.passes()


# 	def _preSelection(self, evt):
# 		###########################################################################################################################\
# 		# Preselection are all the cuts that are requied per event
# 		#Initialize the cut bools every event. These bools tell the code if the cutflow has already been filled for this event.
# 		#Default is to select the first event that passes the selection
# 		###########################################################################################################################

# 		# initiazlie booleans here for each cut you apply to be False
# 		self.passDVtype = False

# 		self.h['CutFlow'][self.ch].Fill(0)

# 		self._fill_histos(evt) # filling all


# 		######################################################################################################
# 		# Selection code is deisgned so that it will pass the selection only if the cut true or cut is unused
# 		# ex. passTrigger is true if the trigcut is true OR if trigcut is not used)
# 		######################################################################################################


# 		# add your pre-selection cuts here


# 	def _DVSelection(self, evt):

# 		######################################################################################################
# 		# DV Selection is any cuts that are done per DV
# 		######################################################################################################

# 		self._fill_all_dv_histos(evt) # Fill all the histograms with ALL DVs (this could be more that 1 per event). Useful for vertexing efficiency studies.



		# # Do we want to use this cut?
		# if self.do_fidvol_cut:
		# 	# Does this cut pass?
		# 	if self._fidvol_cut(evt):
		# 		# Has the cutflow already been filled for this event?
		# 		if not self.passed_fidvol_cut:
		# 			self.h['CutFlow'][self.ch].Fill(5)
		# 			self.passed_fidvol_cut = True
		# 	# If this cut doesn't pass, don't continue to check other cuts
		# 	else:
		# 		return

# 		self._fill_selected_dv_histos(evt) # Fill all the histograms with only selected DVs.
