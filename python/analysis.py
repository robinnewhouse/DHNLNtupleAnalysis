# from ROOT import*
import ROOT
import numpy as np
import os
import helpers
import selections
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
		self.fi = ROOT.TFile.Open(outputFile + '.part', 'recreate')
		self.ch = channel
		self.histSuffixes = []
		self.h = {}
		# make histograms (common for all channels)
		self.add('CutFlow', 13, -0.5, 12.5)
		self.add('trackd0', 200, -10, 10)

		# trigger cut
		if ('alltriggers' in self.sel): 
			self.triggger = 'alltriggers'
			self.dotrigger = True
		else: 
			logger.warn('You did not specify a trigger configuration for this channel. Skipping trigger selection.')
			self.dotriggger = False

		# filter cut
		if ('4-filter' in self.sel):
			self.filter_type = '4-filter'
			self.dofilter = True
		elif ('3-filter' in self.sel):
			self.filter_type = '3-fitler'
			self.dofilter = True
		elif ('2-filter'in self.sel): 
			self.filter_type = '2-filter'
			self.dofilter = True
		elif('1-filter'in self.sel): 
			self.filter_type = '1-fitler'
			self.dofilter = True
		elif('mumu' in self.sel):
			self.filter_type = 'mu-mu'
			self.dofilter = True
		elif('elmu' in self.sel):
			self.filter_type = 'el-mu'
			self.dofilter = True
		elif('elel' in self.sel):
			self.filter_type = 'el-el'
			self.dofilter = True
		elif('muel' in self.sel):
			self.filter_type = 'mu-el'
			self.dofilter = True
		else: 
			logger.warn('You did not specify a filter configuration for this channel. Skipping filter selection.')
			self.dofilter = False

		#prompt lepton cut
		if ('pmuon' in self.sel):
			self.plep = 'muon'
			self.doplep = True
		elif ('pelectron' in self.sel):
			self.plep = 'electron'
			self.doplep = True
		else: 
			logger.warn('You did not specify a prompt lepton for this channel. Skipping prompt lepton selection.')
			self.doplep = False

		#nDV cut 
		if ('nDV' in self.sel): 
			self.donDV = True
		else: 
			logger.warn('You did not add nDV cut. Skipping nDV selection.')
			self.donDV = False

		# OS cut
		if ('OS' in self.sel):
			self.doOS = True
		else: 
			logger.warn('You did not add an OS track cut. Skipping OS track selection.')
			self.doOS = False

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




	def unlock(self):
		self._locked = UNLOCKED

	def add(self, hName, nBins, xLow, xHigh):
		self.h[hName] = {}
		#self.fi.cd()
		self.h[hName] = ROOT.TH1D(hName+self.ch, "", nBins, xLow, xHigh)
		self.h[hName].Sumw2()
		self.h[hName].SetDirectory(0)
		# for s in self.histSuffixes:
			#print 'adding histogram with name ', hName+self.ch+s
			# self.h[hName][s] = ROOT.TH1D(hName+self.ch+s, ', nBins, xLow, xHigh)
			# self.h[hName][s].Sumw2()
			# self.h[hName][s].SetDirectory(0)


	def add2D(self, hName, nBins, xLow, xHigh, nBinsY, yLow, yHigh):
		self.h[hName] = {}
		#self.fi.cd()
		for s in self.histSuffixes:
			#print 'adding histogram with name ', hName+self.ch+s
			self.h[hName][s] = ROOT.TH2D(hName+self.ch+s, "", nBins, xLow, xHigh, nBinsY, yLow, yHigh)
			self.h[hName][s].Sumw2()
			self.h[hName][s].SetDirectory(0)

	def write(self):
		self.fi.cd()
		for hName in self.h:
			# for s in self.histSuffixes:
			self.h[hName].Write(hName)
		self.fi.Close()

	def end(self):
		logger.info('Channel("{}")'.format(self.ch))
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
		head, sep, tail = self._outputFile.partition('file://')
		f = tail if head == '' else self._outputFile
		try:
			os.rename(f + '.part', f)
		except OSError as e:
			logger.error(e, exc_info=True)


	# def _selectChannel(self, sel, syst):
	# 	raise NotImplementedError

	# def selectChannel(self, sel, syst):
	# 	if self._locked != UNLOCKED: # this locks the selectors to avoid rerun all the selections for the same events with different systematics
	# 		# logger.info('<{}, {}> cached'.format(sel, syst))
	# 		return self._passed
	# 	# logger.info('<{}, {}> new one'.format(sel, syst))
	# 	self._passed = self._selectChannel(sel, syst)
	# 	self._locked = SELECTION_LOCKED
	# 	return self._passed

	# def _run(self, sel, syst, wo, wTruth):
	# 	raise NotImplementedError

	# def run(self, sel, syst, wo, wTruth):
	# 	assert self._locked != UNLOCKED, '`run` can be only executed after `selectChannel`'
	# 	self._run(sel, syst, wo, wTruth)
	# 	self._locked = FILL_LOCKED

	def triggerSel(self, evt): 
		self.trigger = selections.Trigger(evt = evt, plepton='muon',trigger = 'HLT_mu26_ivarmedium')
		# print self.trigger
		if self.trigger.passes(): 
			return True
		else: 
			return False

	# def passSel(self, sel): 
	# 	# self.trigger = selections.Trigger(evt = evt, plepton='muon',trigger = 'HLT_mu26_ivarmedium')
	# 	# print self.trigger
	# 	if sel.passes(): 
	# 		return True
	# 	else: 
	# 		return False

	def _preSelection(self, evt):
		raise NotImplementedError

	def preSelection(self, evt): 
		# assert self._locked != UNLOCKED, '`run` can be only executed after `selectChannel`'
		presel = self._preSelection(evt)
		return presel
		# self._locked = FILL_LOCKED

	def _DVSelection(self, evt):
		raise NotImplementedError

	def DVSelection(self, evt): 
		# assert self._locked != UNLOCKED, '`run` can be only executed after `selectChannel`'
		self._DVSelection(evt)
		# self._locked = FILL_LOCKED




class WmuHNL(Analysis):
	mapSel = { 
			   'mumu' : ['alltriggers','pmuon', '4-filter','DV' 'mumu'],   # put a map for a 1 one word key to a list of inputs for the selections
			   'emu'  : ['alltriggers','pmuon', '4-filter' 'emu']}

	def _init__(self, channel, selections, outputFile):
		Analysis.__init__(self, channel, outputFile)
		#make histograms (specfic to this channel)
		self.add('mupt', 100, 0, 100)

	def _trigCut(self, evt): 

		if ('alltriggers' in self.sel):	
			trigger_sel = selections.Trigger(evt = evt,trigger = 'all') 
			return trigger_sel.passes()
		else: 
			return "unused cut"
		
	def _filterCut(self, evt): 
		if self.dofilter: 
			filter_sel =  selections.Filter(evt= evt, _filter=self.filter_type)
			return filter_sel.passes()
		else: 
			return "unused cut"

	def _plepCut(self, evt): 
		if self.doplep: 
			self.plep_sel = selections.Plepton(evt = evt, lepton=self.plep)
			return self.plep_sel.passes()
		else: 
			return "unused cut"

	def _nDVCut(self, evt): 
		if self.donDV: 
			DV_sel = selections.nDV(evt=evt)
			return DV_sel.passes()
		else: 
			return "unused cut"

	def _fidvolCut(self, evt): 
		if ('fidvol' in self.sel): 
			fidvol_sel = selections.DVradius(evt= evt)
			return fidvol_sel.passes()
		else: 
			return "unused cut"

	def _ntrackCut(self, evt): 
		if ('2track' in self.sel):
			ntracks_sel = selections.DVntracks(evt= evt,ntrk=2)
			return ntracks_sel.passes()
		else: 
			return "unused cut"

	def _OSCut(self, evt): 
		if self.doOS: 
			# print "OS in channel sel"
			os_sel = selections.OSDV(evt= evt)
			return os_sel.passes()
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
			return cosmicveto_sel.passes()
		else: 
			return "unused cut"

	def _mlllCut(self, evt): 

		plep_vec = self.plep_sel.plepVec

		muons = helpers.Tracks()
		muons.getMuons(evt= evt)
		muVec = muons.lepVec

		electrons = helpers.Tracks()
		electrons.getElectrons(evt= evt)
		elVec = electrons.lepVec

		if self.domlll: 
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
		if cut == True:  # select events that pass the trigger 
			if passCut == False: 
				self.h['CutFlow'].Fill(nbin)
			return True
		elif cut == False:
			return False
		elif cut == "unused cut": 
			return True


	

	def _preSelection(self, evt):
		# if self.ch not in self.mapSel:
		# 	logger.warn('The selected channel '{}' is not registered. The events will be processed anyway without any further constraint.'.format(self.ch))
		# 	self.mapSel[self.ch] = [self.ch]

		#initialize the DV cuts every event
		self.passTrigger = False
		self.passHNLfilter = False
		self.passPlep = False
		self.passnDV = False
		self.passFid = False
		self.passntracks = False
		self.passOSDV = False
		self.passDVtype = False
		self.passTrackqual = False
		self.passTrackqual_2 = False
		self.passCosmicveto = False
		self.passMlll = False
		self.passDVmass = False

		self.h['CutFlow'].Fill(0)
		
		trigCut = self._doCut(self._trigCut(evt), self.passTrigger, 1)
		if trigCut == True: 
			self.passTrigger = True
		else: 
			return


		filterCut = self._doCut(self._filterCut(evt), self.passHNLfilter, 2)
		if filterCut == True: 
			self.passHNLfilter = True
		else: 
			return

		plepCut = self._doCut(self._plepCut(evt), self.passPlep, 3)
		if plepCut == True: 
			self.passPlep = True
		else: 
			return

		DVCut = self._doCut(self._nDVCut(evt), self.passnDV, 4)
		if DVCut == True: 
			self.passnDV = True
		else: 
			return


		# self._doCut(self._nDVCut(evt), self.passDV, 4)


		# if self._trigCut(evt) == True: 
		# 	if self.passTrigger ==False: 	
		# 		self.h['CutFlow'].Fill(1)
		# 	self.passTrigger = True
		# elif self._trigCut(evt) == False:
		# 	return 
		# elif self._trigCut(evt) == "unused cut": 
		# 	self.passTrigger = True


		# if self._filterCut(evt) == True: 
		# 	if self.passHNLfilter ==False: 	
		# 		self.h['CutFlow'].Fill(2)
		# 	self.passHNLfilter = True
		# elif self._filterCut(evt) == False:
		# 	return 
		# elif self._trigCut(evt) == "unused cut": 
		# 	self.passHNLfilter = True



		# print self._plepCut(evt)

		# if self._plepCut(evt) == True: 
		# 	if self.passPlep ==False: 
		# 		self.h['CutFlow'].Fill(3)
		# 	self.passPlep = True
		# elif self._plepCut(evt) == False:
		# 	return 
		# elif self._trigCut(evt) == "unused cut": 
		# 	self.passPlep = True

		# if self._nDVCut(evt) == True: 
		# 	if self.passDV == False:	
		# 		self.h['CutFlow'].Fill(4)
		# 	self.passDV = True
		# elif self._nDVCut(evt) == False:
		# 	return 
		# elif self._trigCut(evt) == "unused cut": 
		# 	self.passDV = True

		if self.passTrigger and self.passHNLfilter and self.passPlep and self.passnDV: 
			return True
		else: 
			return False


	def _DVSelection(self, evt):
		
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

		OSCut = self._doCut(self._OSCut(evt), self.passOSDV, 7)
		if OSCut == True: 
			self.passOSDV = True
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



	

	

		

		

		

		

		

		




		# if self._fidvolCut(evt) == True: 
		# 	if self.passFid == False: 
		# 		self.h['CutFlow'].Fill(5)
		# 	self.passFid = True
		# elif self._fidvolCut(evt) == False:
		# 	return 
		# elif self._fidvolCut(evt) == "unused cut": 
		# 	self.passFid = True 

		# if self._ntrackCut(evt) == True: 
		# 	if self.passntracks == False: 
		# 		self.h['CutFlow'].Fill(6)
		# 	self.passntracks = True
		# elif self._ntrackCut(evt) == False:
		# 	return 
		# elif self._ntrackCut(evt) == "unused cut": 
		# 	self.passntracks = True 


		# if self._OSCut(evt) == True:
		# 	if self.passOSDV == False: 
		# 		self.h['CutFlow'].Fill(7)
		# 	self.passOSDV = True
		# elif self._OSCut(evt) == False:
		# 	return 
		# elif self._OSCut(evt) == "unused cut": 
		# 	self.passOSDV = True 

		# if self._DVtypeCut(evt) == True:
		# 	if self.passDVtype == False: 
		# 		self.h['CutFlow'].Fill(8)
		# 	self.passDVtype = True
		# elif self._DVtypeCut(evt) == False:
		# 	return 
		# elif self._DVtypeCut(evt) == "unused cut": 
		# 	self.passDVtype = True 
		
		# if self._trackqualCut(evt) == True:
		# 	if self.passTrackqual == False: 
		# 		self.h['CutFlow'].Fill(9)
		# 	self.passTrackqual = True
		# elif self._trackqualCut(evt) == False:
		# 	return 
		# elif self._trackqualCut(evt) == "unused cut": 
		# 	self.passTrackqual = True 

		# if self._cosmicvetoCut(evt) == True: 
		# 	if self.passCosmicveto == False: 
		# 		self.h['CutFlow'].Fill(10)
		# 	self.passCosmicveto = True
		# elif self._cosmicvetoCut(evt) == False:
		# 	return 
		# elif self._cosmicvetoCut(evt) == "unused cut": 
		# 	self.passCosmicveto = True 


		# if self._mlllCut(evt) == True: 
		# 	if self.passMlll == False: 
		# 		self.h['CutFlow'].Fill(11)
		# 	self.passMlll = True
		# elif self._mlllCut(evt) == False:
		# 	return 
		# elif self._mlllCut(evt) == "unused cut": 
		# 	self.passMlll = True 

		# if self._DVmassCut(evt) == True: 
		# 	if self.passDVmass == False: 
		# 		self.h['CutFlow'].Fill(12)
		# 	self.passDVmass = True
		# elif self._DVmassCut(evt) == False:
		# 	return 
		# elif self._DVmassCut(evt) == "unused cut": 
		# 	self.passDVmass = True 


































