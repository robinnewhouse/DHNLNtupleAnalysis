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

	# def passSel(self, sel, cutflowbin, DVpass):
	# 	# if DVpass == None: 
	# 	# 	if sel.passes(): 
	# 	# 		self.h['CutFlow'].Fill(cutflowbin)
	# 	# 	# else: 
	# 	# 	# 	continue 

	# 	# else: 
	# 	if sel.passes(): 
	# 		if DVpass == False: 
	# 			self.h['CutFlow'].Fill(cutflowbin)
	# 		DVpass = True
			# print DVpass
			# else: 
			# 	continue



class WmuHNL(Analysis):
	mapSel = { 
			   'mumu' : ['alltriggers','pmuon', '4-filter','DV' 'mumu'],   # put a map for a 1 one word key to a list of inputs for the selections
			   'emu'  : ['alltriggers','pmuon', '4-filter' 'emu']}

	def _init__(self, channel, selections, outputFile):
		Analysis.__init__(self, channel, outputFile)
		#make histograms (specfic to this channel)
		self.add('mupt', 100, 0, 100)

	def _trigSelection(self, evt): 

		if ('alltriggers' in self.sel):	
			trigger_sel = selections.Trigger(evt = evt,trigger = 'all') 
			return trigger_sel.passes()
		else: 
			return "unused cut"
		
	def _filterSelection(self, evt): 
		if self.dofilter: 
			filter_sel =  selections.Filter(evt= evt, _filter=self.filter_type)
			return filter_sel.passes()
		else: 
			return "unused cut"

	def _plepSelection(self, evt): 
		if self.doplep: 
			self.plep_sel = selections.Plepton(evt = evt, lepton=self.plep)
			return self.plep_sel.passes()
		else: 
			return "unused cut"

	def _nDVSelection(self, evt): 
		if self.donDV: 
			DV_sel = selections.nDV(evt=evt)
			return DV_sel.passes()
		else: 
			return "unused cut"

	def _fidvolSelection(self, evt): 
		if ('fidvol' in self.sel): 
			fidvol_sel = selections.DVradius(evt= evt)
			return fidvol_sel.passes()
		else: 
			return "unused cut"

	def _ntrackSelection(self, evt): 
		if ('2track' in self.sel):
			ntracks_sel = selections.DVntracks(evt= evt,ntrk=2)
			return ntracks_sel.passes()
		else: 
			return "unused cut"

	def _OSSelection(self, evt): 
		if self.doOS: 
			# print "OS in channel sel"
			os_sel = selections.OSDV(evt= evt)
			return os_sel.passes()
		else: 
			return "unused cut"

	def _DVtypeSelection(self, evt):  
		if self.doDVtype: 
			DV_sel = selections.DVtype(evt= evt,decayprod=self.DVtype)
			return DV_sel.passes()
		else: 
			return "unused cut"

	def _trackqualSelection(self, evt): 
		if self.dotrackqual: 
			trackqual_sel = selections.Trackqual(evt=evt, quality=self.trackqual)
			return trackqual_sel.passes()
		else: 
			return "unused cut"

	def _cosmicvetoSelection(self, evt): 
		if self.docosmicveto: 
			cosmicveto_sel = selections.Cosmicveto(evt= evt)
			return cosmicveto_sel.passes()
		else: 
			return "unused cut"

	def _mlllSelection(self, evt): 

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

	def _DVmassSelection(self, evt): 
		if self.doDVmass: 
			DVmass_sel = selections.DVmasscut(evt= evt)
			return DVmass_sel.passes()
		else: 
			return "unused cut"


	

	def _preSelection(self, evt):
		# if self.ch not in self.mapSel:
		# 	logger.warn('The selected channel '{}' is not registered. The events will be processed anyway without any further constraint.'.format(self.ch))
		# 	self.mapSel[self.ch] = [self.ch]

		#initialize the DV cuts every event
		self.passTrigger = False
		self.passHNLfilter = False
		self.passPlep = False
		self.passDV = False
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
		
		if self._trigSelection(evt) == True:  # select events that pass the trigger 
			if self.passTrigger == False: 
				self.h['CutFlow'].Fill(1)
			self.passTrigger = True
		elif self._trigSelection(evt) == False:
			return 
		elif self._trigSelection(evt) == "unused cut": 
			self.passTrigger = True

		if self._filterSelection(evt) == True: 
			cutflow_pass = self.passTrigger
			if self.passHNLfilter ==False and cutflow_pass: 
				self.h['CutFlow'].Fill(2)
			self.passHNLfilter = True
		elif self._filterSelection(evt) == False:
			return 
		elif self._trigSelection(evt) == "unused cut": 
			self.passHNLfilter = True

		# print self._plepSelection(evt)

		if self._plepSelection(evt) == True: 
			cutflow_pass = self.passTrigger and self.passHNLfilter
			if self.passPlep ==False and cutflow_pass: 
				self.h['CutFlow'].Fill(3)
			self.passPlep = True
		elif self._plepSelection(evt) == False:
			return 
		elif self._trigSelection(evt) == "unused cut": 
			self.passPlep = True

		if self._nDVSelection(evt) == True: 
			cutflow_pass = self.passTrigger and self.passHNLfilter and self.passPlep
			if self.passDV == False and cutflow_pass: 
				self.h['CutFlow'].Fill(4)
			self.passDV = True
		elif self._nDVSelection(evt) == False:
			return 
		elif self._trigSelection(evt) == "unused cut": 
			self.passDV = True

		if self.passTrigger and self.passHNLfilter and self.passPlep and self.passDV: 
			return True
		else: 
			return False


	def _DVSelection(self, evt):

		if self._fidvolSelection(evt) == True: 
			if self.passFid == False: 
				self.h['CutFlow'].Fill(5)
			self.passFid = True
		elif self._fidvolSelection(evt) == False:
			return 
		elif self._fidvolSelection(evt) == "unused cut": 
			self.passFid = True 

		if self._ntrackSelection(evt) == True: 
			cutflow_pass = self.passFid
			if self.passntracks == False and cutflow_pass: 
				self.h['CutFlow'].Fill(6)
			self.passntracks = True
		elif self._ntrackSelection(evt) == False:
			return 
		elif self._ntrackSelection(evt) == "unused cut": 
			self.passntracks = True 


		if self._OSSelection(evt) == True:
			cutflow_pass = self.passFid and self.passntracks
			if self.passOSDV == False and cutflow_pass: 
				self.h['CutFlow'].Fill(7)
			self.passOSDV = True
		elif self._OSSelection(evt) == False:
			return 
		elif self._OSSelection(evt) == "unused cut": 
			self.passOSDV = True 

		if self._DVtypeSelection(evt) == True:
			cutflow_pass = self.passFid and self.passntracks and self.passOSDV
			# print cutflow_xpass
			if self.passDVtype == False and cutflow_pass: 
				self.h['CutFlow'].Fill(8)
			self.passDVtype = True
		elif self._DVtypeSelection(evt) == False:
			return 
		elif self._DVtypeSelection(evt) == "unused cut": 
			self.passDVtype = True 
		
		if self._trackqualSelection(evt) == True:
			cutflow_pass = self.passFid and self.passntracks and self.passOSDV and self.passDVtype
			if self.passTrackqual == False and cutflow_pass: 
				self.h['CutFlow'].Fill(9)
			self.passTrackqual = True
		elif self._trackqualSelection(evt) == False:
			return 
		elif self._trackqualSelection(evt) == "unused cut": 
			self.passTrackqual = True 

		if self._cosmicvetoSelection(evt) == True: 
			cutflow_pass = self.passFid and self.passntracks and self.passOSDV and self.passDVtype and self.passTrackqual
			if self.passCosmicveto == False and cutflow_pass: 
				self.h['CutFlow'].Fill(10)
			self.passCosmicveto = True
		elif self._cosmicvetoSelection(evt) == False:
			return 
		elif self._cosmicvetoSelection(evt) == "unused cut": 
			self.passCosmicveto = True 


		if self._mlllSelection(evt) == True: 
			cutflow_pass = self.passFid and self.passntracks and self.passOSDV and self.passDVtype and self.passTrackqual and self.passCosmicveto
			if self.passMlll == False and cutflow_pass: 
				self.h['CutFlow'].Fill(11)
			self.passMlll = True
		elif self._mlllSelection(evt) == False:
			return 
		elif self._mlllSelection(evt) == "unused cut": 
			self.passMlll = True 

		if self._DVmassSelection(evt) == True: 
			cutflow_pass = self.passFid and self.passntracks and self.passOSDV and self.passDVtype and self.passTrackqual and self.passCosmicveto and self.passMlll
			if self.passDVmass == False and cutflow_pass: 
				self.h['CutFlow'].Fill(12)
			self.passDVmass = True
		elif self._DVmassSelection(evt) == False:
			return 
		elif self._DVmassSelection(evt) == "unused cut": 
			self.passDVmass = True 


































