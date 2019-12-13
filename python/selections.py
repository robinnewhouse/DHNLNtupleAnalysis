# from ROOT import*
import ROOT
import numpy as np
import helpers



class Trigger():
	def __init__(self, evt,  plepton, trigger):
		self.evt = evt
		self.plepton = plepton 
		self.trigger = trigger

	def passes(self):
		if self.plepton == "muon":
			ntriggers = len(self.evt.tree.passedtriggers[self.evt.ievt])

			self.passtrig = False
			# print ntriggers
			for itrig in range(ntriggers): 
				# print "trigger ", tree.passedtriggers[self.evt.ievt][itrig]
				# print self.trigger
				if self.evt.tree.passedtriggers[self.evt.ievt][itrig] == self.trigger:
					self.passtrig = True
		
			if self.passtrig == True: 
				return True
			else:
				return False




class Filter():
	def __init__(self, evt, _filter):
		self.evt = evt
		self.filter = _filter 
		
	def passes(self):
		if self.filter == "mu-mu":
			if self.evt.tree.mumufilter[self.evt.ievt] == True:
				return True
			else:
				return False

		if self.filter == "mu-el":
			if self.evt.tree.mumufilter[self.evt.ievt] == True:
				return True
			else:
				return False

		if self.filter == "el-mu":
			if self.evt.tree.mumufilter[self.evt.ievt] == True:
				return True
			else:
				return False

		if self.filter == "el-el":
			if self.evt.tree.mumufilter[self.evt.ievt] == True:
				return True
			else:
				return False

class Plepton():
	def __init__(self, evt, lepton, quality="tight", _mindR=0.05):
		self.evt = evt
		self.lepton = lepton
		self.quality = quality 
		self._mindR = _mindR
		
		self.plepVec = ROOT.TLorentzVector(0,0,0,0)
		ndv = len(self.evt.tree.dvx[self.evt.ievt])	

		if self.lepton == "muon":
			if self.quality == "tight": #tight muon is requested
				lepquality = self.evt.tree.tightmu
			if self.quality == "medium":
				lepquality = self.evt.tree.mediummu
			if self.quality == "loose":
				lepquality = self.evt.tree.loosemu

			nleps = len(self.evt.tree.muonpt[self.evt.ievt])
			passPfilter = self.evt.tree.muonpassPfilter

		if self.lepton == "electron":
			if self.quality == "tight": #tight muon is requested
				lepquality = self.evt.tree.tightel
			# if self.quality == "medium":
			# 	lepquality = self.evt.tree.mediumel
			# if self.quality == "loose":
			# 	lepquality = self.evt.tree.looseel

			nleps = len(self.evt.tree.elpt[self.evt.ievt])
			passPfilter = self.evt.tree.electronpassPfilter

						
		self.highestpt_plep = ROOT.TLorentzVector(0,0,0,0)

		for ilep in xrange(nleps): 
			overlap = False
			plepVec_i = ROOT.TLorentzVector()
			if self.lepton == "muon": 
				plepVec_i.SetPtEtaPhiM(self.evt.tree.muonpt[self.evt.ievt][ilep],self.evt.tree.muoneta[self.evt.ievt][ilep],self.evt.tree.muonphi[self.evt.ievt][ilep],self.evt.tree.muonmass[self.evt.ievt][ilep])
			if self.lepton == "electron":
				plepVec_i.SetPtEtaPhiM(self.evt.tree.elpt[self.evt.ievt][ilep],self.evt.tree.eleta[self.evt.ievt][ilep],self.evt.tree.elphi[self.evt.ievt][ilep],self.evt.tree.elmass[self.evt.ievt][ilep])

			if passPfilter[self.evt.ievt][ilep]:
				for idv in xrange(ndv):
					leptracks = helpers.Tracks()
					trackevt = helpers.Event(self.evt.tree, self.evt.ievt, idv)
					leptracks.getTracks(trackevt)
					dlepVec = leptracks.lepVec
					ndtracks = len(dlepVec)
						
					for itr in xrange(ndtracks): 
						dR = np.sqrt((dlepVec[itr].Eta() - plepVec_i.Eta())**2 + (dlepVec[itr].Phi() - plepVec_i.Phi())**2)

						if dR < self._mindR:  # set overlap to true if muon overlaps with displaced track
							overlap = True
			

				if overlap == False:
					# if self.evt.ievt == 424:
					# 	print self.evt.ievt
					# 	print plepVec_i.Pt(),plepVec_i.Eta(),plepVec_i.Phi()
					if lepquality[self.evt.ievt][ilep] == True or self.quality =="None": # if lepton qulaity requirement is met or no lepton quality is required 
						if (plepVec_i.Pt() > self.highestpt_plep.Pt()): # update highestpt_plep vector to find the largest pt prompt lepton
							self.highestpt_plep= plepVec_i 
						

	def passes(self):
		if self.highestpt_plep.Pt() != 0: 
			self.plepVec = self.highestpt_plep
			return True
		else: 
			return False


class nDV():
	def __init__(self, evt):
		self.evt = evt

	def passes(self):
		if len(self.evt.tree.dvx[self.evt.ievt]) > 0:
			return True
		else:
			return False



class DVradius():
	def __init__(self, evt):
		self.evt = evt
		# self._min = _min 
		# self._max = _max 

		self.rdv = -1
		self.ntracks = len(self.evt.tree.dvx[self.evt.ievt])

		if self.ntracks  > 0:
			dx = self.evt.tree.dvx[self.evt.ievt][self.evt.idv]
			dy = self.evt.tree.dvy[self.evt.ievt][self.evt.idv]
			self.rdv = np.sqrt(dx**2 + dy**2)

		# print self.rdv



	def passes(self, _min = 4,_max = 300):
		if (self.rdv > _min and self.rdv < _max):
			return True
		else: 
			return False




class DVntracks():
	def __init__(self, evt, ntrk=2,  decaymode="leptonic"):
		self.evt = evt
		self.ntrk = ntrk
		self.decaymode = decaymode

		self.ntracks = -1 

		if self.decaymode == "leptonic":
			self.ntracks = len(self.evt.tree.trackpt[self.evt.ievt][self.evt.idv])

	def passes(self):
		if self.ntracks == self.ntrk: 			
			return True
		else: 
			return False 




class OSDV(): 
	def __init__(self, evt, decaymode="leptonic"): 
		self.evt = evt
		self.decaymode = decaymode

		self.ntracks = -1 
		self.charge_trk1 = -2 # dont make default -1 since that's a valid charge! :)
		self.charge_trk2 = -2 # dont make default -1 since that's a valid charge! :)

		if self.decaymode == "leptonic":
			self.ntracks = len(self.evt.tree.trackpt[self.evt.ievt][self.evt.idv])
			# if self.evt.ievt == 217: 
			# 	print "-----"
			# 	print self.evt.ievt
			# 	print tree.trackcharge[self.evt.ievt][self.evt.idv][0]
			# 	print tree.trackcharge[self.evt.ievt][self.evt.idv][1]
			# 	print tree.trackpt[self.evt.ievt][self.evt.idv][0], tree.tracketa[self.evt.ievt][self.evt.idv][0],tree.trackphi[self.evt.ievt][self.evt.idv][0]
			# 	print tree.trackpt[self.evt.ievt][self.evt.idv][1], tree.tracketa[self.evt.ievt][self.evt.idv][1],tree.trackphi[self.evt.ievt][self.evt.idv][1]
			# if self.evt.ievt ==289 or self.evt.ievt ==708 or self.evt.ievt ==805 or self.evt.ievt ==2038 or self.evt.ievt ==2300:

			# 	print "----------"
			# 	print self.evt.ievt

			if self.ntracks == 2: 
				self.charge_trk1 = self.evt.tree.trackcharge[self.evt.ievt][self.evt.idv][0]
				self.charge_trk2 = self.evt.tree.trackcharge[self.evt.ievt][self.evt.idv][1]

	def passes(self): 
		if self.charge_trk1 != self.charge_trk2: 
			return True
		else:
			return False



class DVtype():
	def __init__(self, evt, decayprod, decaymode="leptonic"):
		self.evt = evt
		self.decaymode = decaymode
		self.decayprod = decayprod
		
		if self.decaymode == "leptonic": 
			self.ntracks = len(self.evt.tree.trackpt[self.evt.ievt][self.evt.idv])
			self.nel = -1
			self.nmu = -1


			self.muons = helpers.Tracks()
			self.muons.getMuons(self.evt)
			self.nmu = len(self.muons.lepVec)
			

			self.electrons = helpers.Tracks()
			self.electrons.getElectrons(self.evt)
			self.nel = len(self.electrons.lepVec)	



			# if self.nel == 1:
			# 	print "----------"
			# 	print self.evt.ievt
			# 	print self.nel
			# 	print
			# 	print "track 1: ", self.electrons.lepVec[0].Pt(), self.electrons.lepVec[0].Eta(), self.electrons.lepVec[0].Phi()
			# 	print "el 1: ", self.evt.tree.elpt[self.evt.ievt][self.electrons.lepIndex[0]], self.evt.tree.eleta[self.evt.ievt][self.electrons.lepIndex[0]], self.evt.tree.elphi[self.evt.ievt][self.electrons.lepIndex[0]]
			# 	print "el 1 quality: ", self.evt.tree.tightmu[self.evt.ievt][self.electrons.lepIndex[0]]
			# print "nmu ", self.nmu
			# print "nel", self.nel
			
			# if self.nmu == 2: 
			# 	print "----------"
			# 	print self.evt.ievt
			# 	print "track 1: ", self.muons.lepVec[0].Pt(), self.muons.lepVec[0].Eta(), self.muons.lepVec[0].Phi()
			# 	print "mu: ", self.evt.tree.muonpt[self.evt.ievt][self.muons.lepIndex[0]], self.evt.tree.muoneta[self.evt.ievt][self.muons.lepIndex[0]], self.evt.tree.muonphi[self.evt.ievt][self.muons.lepIndex[0]]
			# 	print "mu quality: ", self.evt.tree.tightmu[self.evt.ievt][self.muons.lepIndex[0]]
			# 	print ""

			# 	print "track 2: ", self.muons.lepVec[1].Pt(), self.muons.lepVec[1].Eta(), self.muons.lepVec[1].Phi()
			# 	print "mu: ", self.evt.tree.muonpt[self.evt.ievt][self.muons.lepIndex[1]], self.evt.tree.muoneta[self.evt.ievt][self.muons.lepIndex[1]], self.evt.tree.muonphi[self.evt.ievt][self.muons.lepIndex[1]]
			# 	print "mu  quality: ", self.evt.tree.tightmu[self.evt.ievt][self.muons.lepIndex[1]]
			# 	print ""


			# if self.nel ==1 and self.nmu == 1:

			# if self.evt.ievt ==289 or self.evt.ievt ==708 or self.evt.ievt ==805 or self.evt.ievt ==2038 or self.evt.ievt ==2300:
			# if (self.evt.ievt == 8) or (self.evt.ievt == 58) or (self.evt.ievt == 60) or (self.evt.ievt == 82) or (self.evt.ievt == 96):
			# 	print "----------"
			# 	print self.evt.ievt
				# print "n el: ", self.nel 
				# print "n mu: ", self.nmu
				# print "DV: ", self.evt.tree.dvx[self.evt.ievt][self.evt.idv], self.evt.tree.dvy[self.evt.ievt][self.evt.idv], self.evt.tree.dvz[self.evt.ievt][self.evt.idv]
				# print "track 1: ", self.electrons.lepVec[0].Pt(), self.electrons.lepVec[0].Eta(), self.electrons.lepVec[0].Phi()
				# print "el: ", self.evt.tree.elpt[self.evt.ievt][self.electrons.lepIndex[0]], self.evt.tree.eleta[self.evt.ievt][self.electrons.lepIndex[0]], self.evt.tree.elphi[self.evt.ievt][self.electrons.lepIndex[0]]
				# print "el quality: ", self.evt.tree.tightmu[self.evt.ievt][self.electrons.lepIndex[0]]
				# print ""
				# print "DV: ", self.evt.tree.dvx[self.evt.ievt][self.evt.idv], self.evt.tree.dvy[self.evt.ievt][self.evt.idv], self.evt.tree.dvz[self.evt.ievt][self.evt.idv]
				# print "track 1: ", self.muons.lepVec[0].Pt(), self.muons.lepVec[0].Eta(), self.muons.lepVec[0].Phi()
				# print "mu: ", self.evt.tree.muonpt[self.evt.ievt][self.muons.lepIndex[0]], self.evt.tree.muoneta[self.evt.ievt][self.muons.lepIndex[0]], self.evt.tree.muonphi[self.evt.ievt][self.muons.lepIndex[0]]
				# print "mu quality: ", self.evt.tree.tightmu[self.evt.ievt][self.muons.lepIndex[0]]
				# print ""
				# print "track 2: ", self.muons.lepVec[1].Pt(), self.muons.lepVec[1].Eta(), self.muons.lepVec[1].Phi()
				# print "mu: ", self.evt.tree.muonpt[self.evt.ievt][self.muons.lepIndex[1]], self.evt.tree.muoneta[self.evt.ievt][self.muons.lepIndex[1]], self.evt.tree.muonphi[self.evt.ievt][self.muons.lepIndex[1]]
				# print "mu  quality: ", self.evt.tree.tightmu[self.evt.ievt][self.muons.lepIndex[1]]

		


	def passes(self): 
		combined = 0 

		if self.decayprod == "emu": 
			if self.nel == 1 and self.nmu == 1: 
				mu1_type = self.evt.tree.muontype[self.evt.ievt][self.muons.lepIndex[0]]

				if mu1_type == combined:  # Only count combined muons 
					return True
				else:
					return False
			else:
				return False


		elif self.decayprod == "mumu":
			if self.nmu == 2: 
				mu1_type = self.evt.tree.muontype[self.evt.ievt][self.muons.lepIndex[0]]
				mu2_type = self.evt.tree.muontype[self.evt.ievt][self.muons.lepIndex[1]]

				# print self.muons.lepIndex[0], self.muons.lepIndex[1]

				# if self.evt.ievt == 8:
				# 	print "------"
				# 	print self.evt.ievt
				# 	print mu1_type, mu2_type
				# 	print "track 1: ", self.muons.lepVec[0].Pt(), self.muons.lepVec[0].Eta(), self.muons.lepVec[0].Phi()
				# 	print "track 2: ", self.muons.lepVec[1].Pt(), self.muons.lepVec[1].Eta(), self.muons.lepVec[1].Phi()
				# 	print "muon 1: ", self.evt.tree.muonpt[self.evt.ievt][self.muons.lepIndex[0]], self.evt.tree.muoneta[self.evt.ievt][self.muons.lepIndex[0]], self.evt.tree.muonphi[self.evt.ievt][self.muons.lepIndex[0]]
				# 	print "muon 2: ", self.evt.tree.muonpt[self.evt.ievt][self.muons.lepIndex[1]], self.evt.tree.muoneta[self.evt.ievt][self.muons.lepIndex[1]], self.evt.tree.muonphi[self.evt.ievt][self.muons.lepIndex[1]]
					
				# 	print "muons pt:",  self.evt.tree.muonpt[self.evt.ievt]
				# 	print "muons eta:",  self.evt.tree.muoneta[self.evt.ievt]
				# 	print "muons phi:",  self.evt.tree.muonphi[self.evt.ievt]

				if (mu1_type == combined and mu2_type == combined) :  # Only count combined muons 
					return True
				else:
					return False
			else:
				return False

		elif self.decayprod == "ee":
			if self.nel == 2: 
				return True
			else:
				return False



class Trackqual():
	def __init__(self,evt, decaymode="leptonic", quality="2-tight"):
		self.evt = evt
		self.decaymode = decaymode
		self.quality = quality 

		if self.decaymode == "leptonic": 
			self.ntracks = len(self.evt.tree.trackpt[self.evt.ievt][self.evt.idv])
			self.ntight = 0
			self.DVmuons = []
			self.DVelectrons = []

			muons = helpers.Tracks()
			muons.getMuons(self.evt)
		
			

			electrons = helpers.Tracks()
			electrons.getElectrons(self.evt)
			self.nel = len(electrons.lepVec)

			self.nmu_tight = 0
			self.nel_tight = 0

			self.ndvmu = len(muons.lepVec)
			self.ndvel = len(electrons.lepVec)
		
			for imu in range(self.ndvmu):
				muindex = muons.lepIndex[imu]
				muisTight = self.evt.tree.tightmu[self.evt.ievt][muindex]
				if muisTight: 
					self.nmu_tight = self.nmu_tight + 1

			for iel in range(self.ndvel):
				elindex = electrons.lepIndex[iel]
				# print elindex
				# print len(self.evt.tree.tightel[self.evt.ievt])
				elisTight = self.evt.tree.tightel[self.evt.ievt][elindex]
				if elisTight: 
					self.nel_tight = self.nel_tight + 1

			# if (self.evt.ievt == 875) or (self.evt.ievt == 2115) or (self.evt.ievt == 2995) or (self.evt.ievt == 44464) or (self.evt.ievt == 339):
			# print "----------"
			# print self.evt.ievt
			# print "track 1: ", electrons.lepVec[0].Pt(), electrons.lepVec[0].Eta(), electrons.lepVec[0].Phi()
			# print "el 1: ", self.evt.tree.elpt[self.evt.ievt][electrons.lepIndex[0]], self.evt.tree.eleta[self.evt.ievt][electrons.lepIndex[0]], self.evt.tree.elphi[self.evt.ievt][electrons.lepIndex[0]]
			# print "el 1 quality: ", self.evt.tree.tightmu[self.evt.ievt][electrons.lepIndex[0]]
			# print ""
			# print "track 2: ", electrons.lepVec[1].Pt(), electrons.lepVec[1].Eta(), electrons.lepVec[1].Phi()	
			# print "el 2: ", self.evt.tree.elpt[self.evt.ievt][electrons.lepIndex[1]], self.evt.tree.eleta[self.evt.ievt][electrons.lepIndex[1]], self.evt.tree.elphi[self.evt.ievt][electrons.lepIndex[1]]
			# print "el 2 quality: ", self.evt.tree.tightmu[self.evt.ievt][electrons.lepIndex[1]]

			# print "number of tight electrons",self.nel_tight


			# for iel in range(self.ndvel): 
			# 	elindex = electrons.lepIndex[iel]
			# 	elisTight = self.evt.tree.tightel[self.evt.ievt][elindex]
			# 	if elisTight:
			# 		self.nel_tight = self.nel_tight + 1


	def passes(self):
			if self.quality == "2-tight": 
				if (self.nmu_tight == 2 or self.nel_tight == 2 or (self.nmu_tight == 1 and self.nel_tight == 1) ):
					return True
				else:
					return False

			if self.quality == "1-tight": 
				if (self.nmu_tight == 1 or self.nel_tight == 1):
					return True
				else: 
					return False

class Cosmicveto():
	def __init__(self, evt,  decaymode="leptonic",cosmicvetocut=0.05 ):
		self.evt = evt
		self.decaymode = decaymode
		self.cosmicvetocut = cosmicvetocut

		self.seperation = -1 

		if self.decaymode == "leptonic": 
			ntracks = len(self.evt.tree.trackpt[self.evt.ievt][self.evt.idv])
			if ntracks == 2: 

				sumeta = self.evt.tree.tracketa[self.evt.ievt][self.evt.idv][0] + self.evt.tree.tracketa[self.evt.ievt][self.evt.idv][1]
				dphi = abs(self.evt.tree.trackphi[self.evt.ievt][self.evt.idv][0] - self.evt.tree.trackphi[self.evt.ievt][self.evt.idv][1])

				self.separation = np.sqrt( sumeta**2 + (np.pi -dphi)**2 )


	def passes(self):		
		if (self.separation > self.cosmicvetocut):
			return True
		else:
			return False
		


class Mlll():
	def __init__(self, decayprod, plep, dMu, dEl, decaymode="leptonic", 
				_minmlll= 50 , _maxmlll = 84):
		self.decaymode = decaymode
		self.decayprod = decayprod
		self.plep = plep
		self.dMu = dMu
		self.dEl = dEl
		self._minmlll = _minmlll
		self._maxmlll = _maxmlll

		self.mlll = -1
		self.plll = ROOT.TLorentzVector(0,0,0,0)

		if self.decaymode == "leptonic":	
		
			if self.decayprod == "emu": 
				self.plll = self.plep + self.dEl[0] + self.dMu[0]
				self.mlll = self.plll.M()

			if self.decayprod == "mumu": 
				self.plll = self.plep + self.dMu[0] + self.dMu[1]
				self.mlll = self.plll.M()

			if self.decayprod == "ee": 
				self.plll = self.plep + self.dEl[0] + self.dEl[1]
				self.mlll = self.plll.M()

	def passes(self):
		
		if (self.mlll> self._minmlll and self.mlll < self._maxmlll):
			return True
		else: 
			return False



class DVmasscut():
	def __init__(self, evt,  decaymode="leptonic",dvmasscut=4):
		self.evt = evt
		self.decaymode = decaymode
		self.dvmasscut = dvmasscut

	def passes(self):
		dvmass = self.evt.tree.dvmass[self.evt.ievt][self.evt.idv]
		if (dvmass > self.dvmasscut):
			return True
		else: 
			return False







