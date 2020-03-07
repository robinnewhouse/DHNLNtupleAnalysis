# from ROOT import*
import ROOT
import numpy as np
import os
import helpers
import selections
import treenames
import observables
logger = helpers.getLogger('dHNLAnalysis.analysis')





UNLOCKED = 0
SELECTION_LOCKED = 1
FILL_LOCKED = 2


class Analysis(object):

	ch =''
	h = {} #map of hisotgram names to map of systematics to histograms
	blinded = True
	mapSel = {}

	def __init__(self, channel, selections, outputFile):
		# print channel
		# print outputFile

		self.sel = selections
		self._outputFile = outputFile
		# self.fi = ROOT.TFile.Open(outputFile + '.part', 'recreate')
		self.fi = ROOT.TFile.Open(outputFile , 'update')
		self.ch = channel
		logger.info('Initializing Channel("{}")'.format(self.ch))
		self.histSuffixes = [self.ch]
		self.h = {}
		# self.histSuffixes = self.histSuffixes + [self.ch]
		# print self.histSuffixes
		# make histograms (common for all channels)
		self.add('CutFlow', 13, -0.5, 12.5)
		self.add2D('charge_ntrk', 11, -5.5, 5.5,9,-0.5,8.5)
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

		# trigger cut
		if ('alltriggers' in self.sel): 
			self.triggger = 'all'
			self.do_trigger_cut = True
		else: 
			logger.warn('You did not specify a trigger configuration for this channel. Skipping trigger selection.')
			self.dotriggger = False

		# filter cut
		if ('4-filter' in self.sel):
			self.filter_type = '4-filter'
			self.do_filter_cut = True
		elif ('3-filter' in self.sel):
			self.filter_type = '3-fitler'
			self.do_filter_cut = True
		elif ('2-filter'in self.sel): 
			self.filter_type = '2-filter'
			self.do_filter_cut = True
		elif('1-filter'in self.sel): 
			self.filter_type = '1-fitler'
			self.do_filter_cut = True
		elif('mumu-filter' in self.sel):
			self.filter_type = 'mu-mu'
			self.do_filter_cut = True
		elif('elmu-filter' in self.sel):
			self.filter_type = 'el-mu'
			self.do_filter_cut = True
		elif('elel-filter' in self.sel):
			self.filter_type = 'el-el'
			self.do_filter_cut = True
		elif('muel-filter' in self.sel):
			self.filter_type = 'mu-el'
			self.do_filter_cut = True
		else: 
			logger.warn('You did not specify a filter configuration for this channel. Skipping filter selection.')
			self.do_filter_cut = False

		#prompt lepton cut
		if ('pmuon' in self.sel):
			self.plep = 'muon'
			self.do_prompt_lepton_cut = True
		elif ('pelectron' in self.sel):
			self.plep = 'electron'
			self.do_prompt_lepton_cut = True
		else: 
			logger.warn('You did not specify a prompt lepton for this channel. Skipping prompt lepton selection.')
			self.do_prompt_lepton_cut = False

		# nDV cut 
		self.do_ndv_cut = ('nDV' in self.sel)
		if not self.do_ndv_cut: logger.warn('You did not add nDV cut. Skipping nDV selection.')

		# fiducial volume
		if ('fidvol' in self.sel):
			self.dofivol = True
		else:
			logger.warn('You did not add DV in fiducial volume cut. Skipping DV in fiducial volume selection.')
			self.dofivol = False

		#2 track cut
		if ('2track' in self.sel): 
			self.donTrack = True
			self.ntrk = 2
		elif ('3track' in self.sel): 
			self.donTrack = True
			self.ntrk = 3
		elif ('4track' in self.sel): 
			self.donTrack = True
			self.ntrk = 4
		else: 
			logger.warn('You did not add an ntrack cut. Skipping ntrack selection.')
			self.donTrack = False


		# OS cut
		if ('OS' in self.sel):
			self.doOS = True
		else: 
			logger.warn('You did not add an OS track cut. Skipping OS track selection.')
			self.doOS = False

		if ('SS' in self.sel): 
			self.doSS = True
		else: 
			self.doSS = False

		# DV type
		if ('mumu' in self.sel): 
			self.DVtype = "mumu"
			self.doDVtype = True
		elif('emu' in self.sel):
			self.DVtype = "emu"
			self.doDVtype = True
		elif('ee' in self.sel):
			self.DVtype = "ee"
			self.doDVtype = True
		elif('1-lep' in self.sel):
			self.DVtype = "1-lep"
			self.doDVtype = True
		elif('2-lep' in self.sel):
			self.DVtype = "2-lep"
			self.doDVtype = True
		else: 
			logger.warn('You did not specify a DV type for this channel. Skipping DV type selection.')
			self.doDVtype = False

		#track quality 
		if ('1-tight' in self.sel):
			self.trackqul = '1-tight'
			self.dotrackqual = True
		elif ('2-tight' in self.sel): 
			self.trackqual = '2-tight'
			self.dotrackqual = True
		else: 
			logger.warn('You did not specify a DV track quality for this channel. Skipping DV track quality selection.')
			self.dotrackqual = False

		# cosmic veto cut
		if ('cosmicveto' in self.sel): 
			self.docosmicveto = True
		else: 
			logger.warn('You did not add a cosmic veto cut for this channel. Skipping cosmic veto selection.')
			self.docosmicveto = False


		# mlll cut
		if ('mlll' in self.sel): 
			self.domlll = True
		else: 
			self.domlll = False

		# DV mass cut
		if ('DVmass' in self.sel): 
			self.doDVmass = True
		else: 
			self.doDVmass = False

		if self.domlll == True and self.do_prompt_lepton_cut == False: 
			logger.error("You cannot calculate mlll without selecting a prompt lepton!")
			quit()

		self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(1, "all")
		if self.do_trigger_cut == True:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(2, "trigger")
		if self.do_filter_cut == True:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(3, "%s"%self.filter_type)
		if self.do_prompt_lepton_cut == True:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(4, "tight prompt %s"%self.plep)
		if self.do_ndv_cut == True:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(5, "DV")
		if self.dofivol == True:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(6, "fiducial")
		if self.donTrack == True:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(7, "%s-track DV"%self.ntrk)
		if self.doOS == True:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(8, "OS DV")
		if self.doSS == True:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(8, "SS DV")
		if self.doDVtype == True:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(9, "%s DV"%self.DVtype)
		if self.dotrackqual == True:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(10, "2-tight-lepton DV")
		if self. docosmicveto == True:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(11, "cosmic veto")
		if self.domlll == True:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(12, "m_{lll}")
		if self.doDVmass == True:
			self.h['CutFlow'][self.ch].GetXaxis().SetBinLabel(13, "mDV")




	def unlock(self):
		self._locked = UNLOCKED

	def add(self, hName, nBins, xLow, xHigh):
		self.h[hName] = {}
		#self.fi.cd()
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
		#self.fi.cd()
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
		#self.fi.cd()
		for s in self.histSuffixes:
			# print 'adding histogram with name ', hName+self.ch+s
			self.h[hName][s] = ROOT.TH2D(hName+"_"+s, "", nBins, xLow, xHigh, nBinsY, yLow, yHigh)
			self.h[hName][s].Sumw2()
			self.h[hName][s].SetDirectory(0)


	def write(self):
		self.fi.cd()
		for hName in self.h:
			for s in self.histSuffixes:
				self.h[hName][s].Write(hName+'_'+s)
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
	def _trigger_cut(self, evt):
		trigger_sel = selections.Trigger(evt=evt, trigger=self.triggger)
		return trigger_sel.passes()

	def _filter_cut(self, evt):
		filter_sel = selections.Filter(evt=evt, filter_type=self.filter_type)
		return filter_sel.passes()

	def _prompt_lepton_cut(self, evt):
		self.plep_sel = selections.Plepton(evt=evt, lepton=self.plep)

		# Add to histogram all prompt leptons that pass selection.
		# If _prompt_lepton_cut() is run after trigger and filter cut then those cuts will also be applied.
		self.h["plep_pt"][self.ch].Fill(self.plep_sel.plepVec.Pt())
		self.h["plep_eta"][self.ch].Fill(self.plep_sel.plepVec.Eta())
		self.h["plep_phi"][self.ch].Fill(self.plep_sel.plepVec.Phi())
		self.h["plep_d0"][self.ch].Fill(self.plep_sel.plepd0)
		self.h["plep_z0"][self.ch].Fill(self.plep_sel.plepz0)
		return self.plep_sel.passes()

	def _ndv_cut(self, evt):
		dv_sel = selections.nDV(evt=evt)
		return dv_sel.passes()

	def _fidvolCut(self, evt):
		if self.dofivol == True:
			fidvol_sel = selections.DVradius(evt= evt)
			return fidvol_sel.passes() 
		else: 
			return "unused cut"

	def _ntrackCut(self, evt): 
		if self.donTrack == True:
			ntracks_sel = selections.DVntracks(evt= evt,ntrk=self.ntrk)
			return ntracks_sel.passes()
		else: 
			return "unused cut"

	def _chargeCut(self, evt): 
		if self.doOS: 
			os_sel = selections.ChargeDV(evt= evt,sel='OS')
			return os_sel.passes()
		elif self.doSS: 
			ss_sel = selections.ChargeDV(evt= evt,sel='SS')
			return ss_sel.passes()
		else: 
			return "unused cut"

	def _DVtypeCut(self, evt):  
		if self.doDVtype: 
			DV_sel = selections.DVtype(evt= evt,decayprod=self.DVtype)
			return DV_sel.passes()
		else: 
			return "unused cut"

	def _trackqualCut(self, evt): 
		if self.dotrackqual: 
			trackqual_sel = selections.Trackqual(evt=evt, quality=self.trackqual)
			return trackqual_sel.passes()
		else: 
			return "unused cut"

	def _cosmicvetoCut(self, evt): 
		if self.docosmicveto: 
			cosmicveto_sel = selections.Cosmicveto(evt= evt)
			# self.h["DV_trk_sep"][self.ch].Fill(cosmicveto_sel.separation)
			return cosmicveto_sel.passes()
		else: 
			return "unused cut"

	def _mlllCut(self, evt): 
		if self.domlll: 
			plep_vec = self.plep_sel.plepVec

			muons = helpers.Tracks()
			muons.getMuons(evt= evt)
			muVec = muons.lepVec

			electrons = helpers.Tracks()
			electrons.getElectrons(evt= evt)
			elVec = electrons.lepVec
		
			mlll_sel = selections.Mlll(decayprod=self.DVtype,plep=plep_vec,dMu=muVec,dEl=elVec)
			return mlll_sel.passes()
		else:
			return "unused cut"

	def _DVmassCut(self, evt): 
		if self.doDVmass: 
			DVmass_sel = selections.DVmasscut(evt= evt)
			return DVmass_sel.passes()
		else: 
			return "unused cut"

	def _doCut(self, cut, passCut, nbin):
		##################################################################################################################################
		# DoCut is done for every cut. It needs to know the'cut' that should be applied, and if the cut was already applied ('passCut')
		# This ensures that the cutflow is only filled the first time the cut is passed and does not double count,
		# This should be modified eventually to decide which DV is the "best" as opposed to just counting the first DV found.
		###################################################################################################################################
		if cut == True:  # select events that pass the trigger 
			if passCut == False: 
				self.h['CutFlow'][self.ch].Fill(nbin)
			return True
		elif cut == False:
			return False
		elif cut == "unused cut": 
			return True

	def _fillHistos(self, evt):
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

	def _fillallDVHistos(self, evt): 

		self.h["charge_ntrk"][self.ch].Fill(evt.tree.dvcharge[evt.ievt][evt.idv], evt.tree.dvntrk[evt.ievt][evt.idv])
		ntracks = len(evt.tree.trackd0[evt.ievt][evt.idv])
		for itrk in range(ntracks): # loop over tracks 
			self.h["DV_trk_pt"][self.ch].Fill(evt.tree.trackpt[evt.ievt][evt.idv][itrk])
			self.h["DV_trk_eta"][self.ch].Fill(evt.tree.tracketa[evt.ievt][evt.idv][itrk])
			self.h["DV_trk_phi"][self.ch].Fill(evt.tree.trackphi[evt.ievt][evt.idv][itrk])
			self.h["DV_trk_d0"][self.ch].Fill(evt.tree.trackd0[evt.ievt][evt.idv][itrk])
			self.h["DV_trk_z0"][self.ch].Fill(evt.tree.trackz0[evt.ievt][evt.idv][itrk])
			self.h["DV_trk_charge"][self.ch].Fill(evt.tree.trackcharge[evt.ievt][evt.idv][itrk])
			self.h["DV_trk_chi2"][self.ch].Fill(evt.tree.trackchi2[evt.ievt][evt.idv][itrk])

		self.h["DV_num_trks"][self.ch].Fill(evt.tree.dvntrk[evt.ievt][evt.idv])
		self.h["DV_x"][self.ch].Fill(evt.tree.dvx[evt.ievt][evt.idv]) 
		self.h["DV_y"][self.ch].Fill(evt.tree.dvy[evt.ievt][evt.idv]) 
		self.h["DV_z"][self.ch].Fill(evt.tree.dvz[evt.ievt][evt.idv]) 
		self.h["DV_r"][self.ch].Fill(evt.tree.dvr[evt.ievt][evt.idv]) 
		self.h["DV_distFromPV"][self.ch].Fill(evt.tree.dvdistFromPV[evt.ievt][evt.idv]) 
		self.h["DV_mass"][self.ch].Fill(evt.tree.dvmass[evt.ievt][evt.idv]) 
		self.h["DV_pt"][self.ch].Fill(evt.tree.dvpt[evt.ievt][evt.idv]) 
		self.h["DV_eta"][self.ch].Fill(evt.tree.dveta[evt.ievt][evt.idv]) 
		self.h["DV_phi"][self.ch].Fill(evt.tree.dvphi[evt.ievt][evt.idv]) 
		self.h["DV_minOpAng"][self.ch].Fill(evt.tree.dvminOpAng[evt.ievt][evt.idv]) 
		self.h["DV_maxOpAng"][self.ch].Fill(evt.tree.dvmaxOpAng[evt.ievt][evt.idv]) 
		self.h["DV_charge"][self.ch].Fill(evt.tree.dvcharge[evt.ievt][evt.idv])
		self.h["DV_chi2"][self.ch].Fill(evt.tree.dvchi2[evt.ievt][evt.idv])

	def _fillselectedDVHistos(self, evt): 
		if self._locked < FILL_LOCKED:
			#these are the histograms you only want to fill ONCE per DV
			if self.do_prompt_lepton_cut:
				plep_vec = self.plep_sel.plepVec
				plepd0 = self.plep_sel.plepd0
				plepz0 = self.plep_sel.plepz0

				self.h["selplep_pt"][self.ch].Fill(plep_vec.Pt())
				self.h["selplep_eta"][self.ch].Fill(plep_vec.Eta())
				self.h["selplep_phi"][self.ch].Fill(plep_vec.Phi())
				self.h["selplep_d0"][self.ch].Fill(plepd0)
				self.h["selplep_z0"][self.ch].Fill(plepz0)


			ntracks = len(evt.tree.trackd0[evt.ievt][evt.idv])
			for itrk in range(ntracks): #loop over tracks
				self.h["selDV_trk_pt"][self.ch].Fill(evt.tree.trackpt[evt.ievt][evt.idv][itrk])
				self.h["selDV_trk_eta"][self.ch].Fill(evt.tree.tracketa[evt.ievt][evt.idv][itrk])
				self.h["selDV_trk_phi"][self.ch].Fill(evt.tree.trackphi[evt.ievt][evt.idv][itrk])
				self.h["selDV_trk_d0"][self.ch].Fill(evt.tree.trackd0[evt.ievt][evt.idv][itrk])
				self.h["selDV_trk_z0"][self.ch].Fill(evt.tree.trackz0[evt.ievt][evt.idv][itrk])
				self.h["selDV_trk_charge"][self.ch].Fill(evt.tree.trackcharge[evt.ievt][evt.idv][itrk])
				self.h["selDV_trk_chi2"][self.ch].Fill(evt.tree.trackchi2[evt.ievt][evt.idv][itrk])

			self.h["selDV_num_trks"][self.ch].Fill(evt.tree.dvntrk[evt.ievt][evt.idv])
			self.h["selDV_x"][self.ch].Fill(evt.tree.dvx[evt.ievt][evt.idv]) 
			self.h["selDV_y"][self.ch].Fill(evt.tree.dvy[evt.ievt][evt.idv]) 
			self.h["selDV_z"][self.ch].Fill(evt.tree.dvz[evt.ievt][evt.idv]) 
			self.h["selDV_r"][self.ch].Fill(evt.tree.dvr[evt.ievt][evt.idv]) 
			self.h["selDV_distFromPV"][self.ch].Fill(evt.tree.dvdistFromPV[evt.ievt][evt.idv]) 
			self.h["selDV_mass"][self.ch].Fill(evt.tree.dvmass[evt.ievt][evt.idv]) 
			self.h["selDV_pt"][self.ch].Fill(evt.tree.dvpt[evt.ievt][evt.idv]) 
			self.h["selDV_eta"][self.ch].Fill(evt.tree.dveta[evt.ievt][evt.idv]) 
			self.h["selDV_phi"][self.ch].Fill(evt.tree.dvphi[evt.ievt][evt.idv]) 
			self.h["selDV_minOpAng"][self.ch].Fill(evt.tree.dvminOpAng[evt.ievt][evt.idv]) 
			self.h["selDV_maxOpAng"][self.ch].Fill(evt.tree.dvmaxOpAng[evt.ievt][evt.idv]) 
			self.h["selDV_charge"][self.ch].Fill(evt.tree.dvcharge[evt.ievt][evt.idv])
			self.h["selDV_chi2"][self.ch].Fill(evt.tree.dvchi2[evt.ievt][evt.idv])
			

			self._locked = FILL_LOCKED # this only becomes unlocked after the event loop finishes in makeHistograms so you can only fill one DV from each event. 




class WmuHNL(Analysis):

	isdata = False # not currently used...

	def _init__(self, channel, selections, outputFile):
		Analysis.__init__(self, channel, outputFile)
		#make histograms specfic to this channel
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
		#Initialize the cut bools every event. These bools tell the code if the cutflow has already been filled for this event. 
		#Default is to select the first event that passes the selection
		###########################################################################################################################
		self.passed_preselection_cuts = False
		self.passed_trigger_cut = False
		self.passed_filter_cut = False
		self.passed_prompt_lepton_cut = False
		self.passed_ndv_cut = False
		self.passFid = False
		self.passntracks = False
		self.passChargeDV = False
		self.passDVtype = False
		self.passTrackqual = False
		self.passTrackqual_2 = False
		self.passCosmicveto = False
		self.passMlll = False
		self.passDVmass = False

		self.h['CutFlow'][self.ch].Fill(0)

		self._fillHistos(evt)


		######################################################################################################
		# Selection code is deisgned so that it will pass the selection only if the cut true or cut is unused 
		# ex. passTrigger is true if the trigcut is true OR if trigcut is not used) 
		######################################################################################################

		if self.do_trigger_cut and not self.passed_trigger_cut:
			if self._trigger_cut(evt):
				# Fill the plot at the specified bin
				self.h['CutFlow'][self.ch].Fill(1)
				# Record that this cut has been done so we don't accidentally reject it in later tests
				self.passed_trigger_cut = True
			else:
				return

		if self.do_filter_cut and not self.passed_filter_cut:
			if self._filter_cut(evt):
				self.h['CutFlow'][self.ch].Fill(2)
				self.passed_filter_cut = True
			else:
				return

		if self.do_prompt_lepton_cut and not self.passed_prompt_lepton_cut:
			if self._prompt_lepton_cut(evt):
				self.h['CutFlow'][self.ch].Fill(3)
				self.passed_prompt_lepton_cut = True
			else:
				return

		if self.do_ndv_cut and not self.passed_ndv_cut:
			if self._ndv_cut(evt):
				self.h['CutFlow'][self.ch].Fill(4)
				self.passed_ndv_cut = True
			else:
				return

		# If you've made it here, preselection is passed
		self.passed_preselection_cuts = True

	def _DVSelection(self, evt):

		######################################################################################################
		# DV Selection is any cuts that are done per DV 
		# Current cuts include: fiducial vol, ntrack, OS, DVtype, track quality, cosmic veto, mlll, mDV
		######################################################################################################
		
		self._fillallDVHistos(evt) # Fill all the histograms with ALL DVs (this could be more that 1 per event). Useful for vertexing efficiency studies.

		if self.passed_preselection_cuts: # only do the DV selection if the preselction was passed for the event. 
			
			fidvolCut = self._doCut(self._fidvolCut(evt), self.passFid, 5)
			if fidvolCut == True: 
				self.passFid = True
			else: 
				return

			ntrackCut = self._doCut(self._ntrackCut(evt), self.passntracks, 6)
			if ntrackCut == True: 
				self.passntracks = True
			else: 
				return

			chargeCut = self._doCut(self._chargeCut(evt), self.passChargeDV, 7)
			if chargeCut == True: 
				self.passChargeDV = True
			else: 
				return

			DVtypeCut = self._doCut(self._DVtypeCut(evt), self.passDVtype, 8)
			if DVtypeCut == True: 
				self.passDVtype = True
			else: 
				return

			trackqualCut = self._doCut(self._trackqualCut(evt), self.passTrackqual, 9)
			if trackqualCut == True: 
				self.passTrackqual = True
			else: 
				return

			cosmicvetoCut = self._doCut(self._cosmicvetoCut(evt), self.passCosmicveto, 10)
			if cosmicvetoCut == True: 
				self.passCosmicveto = True
			else: 
				return

			mlllCut = self._doCut(self._mlllCut(evt), self.passMlll, 11)
			if mlllCut == True: 
				self.passMlll = True
			else: 
				return

			DVmassCut = self._doCut(self._DVmassCut(evt), self.passDVmass, 12)
			if DVmassCut == True: 
				self.passDVmass = True	
			else: 
				return
			
			self._fillselectedDVHistos(evt) # Fill all the histrograms with only selected DVs.

	

	
######################################################################################################################
# An example of a new class. Here you could add any new cuts you want without distubing the main analysis cuts.
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
	
# 		DV_sel = selections.DVtype(evt= evt,decayprod="ee")
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

# 		self._fillHistos(evt) # filling all


# 		######################################################################################################
# 		# Selection code is deisgned so that it will pass the selection only if the cut true or cut is unused 
# 		# ex. passTrigger is true if the trigcut is true OR if trigcut is not used) 
# 		######################################################################################################


# 		# add your pre-selection cuts here



# 	def _DVSelection(self, evt):

# 		######################################################################################################
# 		# DV Selection is any cuts that are done per DV 
# 		######################################################################################################
		
# 		self._fillallDVHistos(evt) # Fill all the histograms with ALL DVs (this could be more that 1 per event). Useful for vertexing efficiency studies.

		
# 		# add your DV cuts here
# 		DVtypeCut = self._doCut(self._DVtypeCut(evt), self.passDVtype, 8)
# 		if DVtypeCut == True: 
# 			self.passDVtype = True
# 		else: 
# 			return

		
			
# 		self._fillselectedDVHistos(evt) # Fill all the histrograms with only selected DVs.		

		

		

		

		

		

































