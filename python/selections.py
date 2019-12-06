# from ROOT import*
import ROOT
import numpy as np
import helpers



class Trigger():
	def __init__(self, plepton, trigger):
		self.plepton = plepton 
		self.trigger = trigger

	def passes(self, tree, ievt):
		if self.plepton == "muon":
			ntriggers = len(tree.passedtriggers[ievt])

			self.passtrig = False
			# print ntriggers
			for itrig in range(ntriggers): 
				# print "trigger ", tree.passedtriggers[ievt][itrig]
				# print self.trigger
				if tree.passedtriggers[ievt][itrig] == self.trigger:
					self.passtrig = True
		
			if self.passtrig == True: 
				return True
			else:
				return False




class Filter():
	def __init__(self, _filter):
		self.filter = _filter 
		
	def passes(self,tree,ievt):
		if self.filter == "mu-mu":
			if tree.mumufilter[ievt] == True:
				return True
			else:
				return False

		if self.filter == "mu-el":
			if tree.mumufilter[ievt] == True:
				return True
			else:
				return False

		if self.filter == "el-mu":
			if tree.mumufilter[ievt] == True:
				return True
			else:
				return False

		if self.filter == "el-el":
			if tree.mumufilter[ievt] == True:
				return True
			else:
				return False

class Plepton():
	def __init__(self, lepton, quality="tight", _mindR=0.05):
		self.quality = quality 
		self.lepton = lepton
		self._mindR = _mindR
		self.plepVec = ROOT.TLorentzVector(0,0,0,0)


	def passes(self, tree, ievt):
		ndv = len(tree.dvx[ievt])	

		if self.lepton == "muon":

			if self.quality == "tight": #tight muon is requested
				# muonquality = tree.tightmu
				muonquality = True  # NEED TO FIX THIS TO BE TIGHT MUON
			if self.quality == "medium":
				muonquality = tree.mediummu
			if self.quality == "loose":
				muonquality = tree.loosemu

			nmuons = len(tree.muonpt[ievt])			
			highestpt_pmuon = ROOT.TLorentzVector(0,0,0,0)

			for imu in xrange(nmuons): 
				overlap = False
				pmuVec_i = ROOT.TLorentzVector()
				pmuVec_i.SetPtEtaPhiM(tree.muonpt[ievt][imu],tree.muoneta[ievt][imu],tree.muonphi[ievt][imu],tree.muonmass[ievt][imu])
			
				if tree.muonpassPfilter[ievt][imu]:
					for idv in xrange(ndv):
						leptracks = helpers.Leptons()
						leptracks.getTracks(tree, ievt, idv)
						dlepVec = leptracks.lepVec
						# dlepVec = helpers.Leptons().getTracks(tree, ievt, idv)
						ndtracks = len(dlepVec)
						for idvmu in xrange(ndtracks): 
							dR = np.sqrt((dlepVec[idvmu].Eta() - pmuVec_i.Eta())**2 + (dlepVec[idvmu].Phi() - pmuVec_i.Phi())**2)
							if dR < self._mindR:  # set overlap to true if muon overlaps with displaced track
								overlap = True
			

					if overlap == False:
						# if muonquality[ievt][imu] == True or self.quality =="None": # if muon qulaity requirement is met or no muon quality is required # WANT THIS LINE NOT THE ONE BELOW
						if muonquality == True or self.quality =="None": # if muon qulaity requirement is met or no muon quality is required
							if (pmuVec_i.Pt() > highestpt_pmuon.Pt()): # update highestpt_pmuon vector to find the largest pt prompt muon
								highestpt_pmuon = pmuVec_i 
						


			if highestpt_pmuon.Pt() != 0: 
				self.plepVec = highestpt_pmuon
				return True
			else: 
				return False


class DV():

	def passes(self, tree, ievt	):
		if len(tree.dvx[ievt]) > 0:
			return True
		else:
			return False

class Fiducial():
	def __init__(self, _min = 4, _max = 300):
		self._min = _min 
		self._max = _max 
		self.rdv = -1

	def passes(self,tree,ievt,idv):
		if len(tree.dvx[ievt]) > 0:
			self.rdv = np.sqrt(tree.dvx[ievt][idv]**2 + tree.dvy[ievt][idv]**2)
		if (self.rdv > self._min and self.rdv < self._max):
			return True
		else: 
			return False


class DVntracks():
	def __init__(self, decaymode="leptonic"):
		self.decaymode = decaymode

	def passes(self, tree, ievt, idv):
		self.ntracks = len(tree.trackpt[ievt][idv])

		if self.decaymode == "leptonic": 
			if self.ntracks == 2: 			
				return True
			else: 
				return False 

class OSDV(): 
	def __init__(self, decaymode="leptonic"): 
		self.decaymode = decaymode
		self.ntracks = -1 

	def passes(self, tree, ievt, idv): 

		if self.decaymode == "leptonic":
			self.ntracks = len(tree.trackpt[ievt][idv])
			if self.ntracks == 2: 
				if tree.trackcharge[ievt][idv][0] != tree.trackcharge[ievt][idv][1]: 
					return True
				else:
					return False

class DVtype():
	def __init__(self, decayprod, decaymode="leptonic"):
		self.decaymode = decaymode
		self.decayprod = decayprod
		

	def passes(self, tree, ievt, idv): 

		if self.decaymode == "leptonic": 
			self.ntracks = len(tree.trackpt[ievt][idv])
			self.nel = 0
			self.nmu = 0


			muons = helpers.Leptons()
			muons.getMuons(tree, ievt, idv)
			self.nmu = len(muons.lepVec)
			

			electrons = helpers.Leptons()
			electrons.getElectrons(tree, ievt, idv)
			self.nel = len(electrons.lepVec)

			if self.decayprod == "emu": 
				if self.nel == 1 and self.nmu == 1: 
					index0 = muons.lepIndex[0]
					if tree.muontype[ievt][index0] == 0:  # Only count combined muons 
						return True
					else:
						return False
				else:
					return False

			elif self.decayprod == "mumu":
				if self.nmu == 2: 
					index0 = muons.lepIndex[0]
					index1 = muons.lepIndex[1]
					if (tree.muontype[ievt][index0] == 0 and tree.muontype[ievt][index1] == 0) :  # Only count combined muons 
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
	def __init__(self, decaymode="leptonic", quality="2-tight"):
		# self.decayprod = decayprod
		self.decaymode = decaymode
		self.quality = quality 

	def passes(self, tree, ievt, idv):

		if self.decaymode == "leptonic": 
			self.ntracks = len(tree.trackpt[ievt][idv])
			self.ntight = 0
			self.DVmuons = []
			self.DVelectrons = []

			muons = helpers.Leptons()
			muons.getMuons(tree, ievt, idv)
			self.nmu = len(muons.lepVec)
			

			electrons = helpers.Leptons()
			electrons.getElectrons(tree, ievt, idv)
			self.nel = len(electrons.lepVec)

			self.nmu_tight = 0
			self.nel_tight = 0
			ndvmuons = len(muons.lepVec)
			ndvelectrons = len(electrons.lepVec)
		
			for imu in range(ndvmuons):
				index = muons.lepIndex[imu]
				if tree.tightmu[ievt][index] : 
					self.nmu_tight = self.nmu_tight + 1

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

class DVmasscut():
	def __init__(self, decaymode="leptonic",dvmasscut=4):
		self.decaymode = decaymode
		self.dvmasscut = dvmasscut

	def passes(self, tree, ievt, idv):
		if (tree.dvmass[ievt][idv] > self.dvmasscut):
			return True
		else: 
			return False

class Mlllcut():
	def __init__(self, decayprod, plep, dMu, dEl, decaymode="leptonic", _minmlll= 50 , _maxmlll = 84):
		self.decaymode = decaymode
		self.decayprod = decayprod
		self.plep = plep
		self.dMu = dMu
		self.dEl = dEl
		self._minmlll = _minmlll
		self._maxmlll = _maxmlll

		self.mlll = -1
		self.plll = ROOT.TLorentzVector(0,0,0,0)

	def passes(self):
		
		if self.decaymode == "leptonic":	
		
			if self.decayprod == "emu": 
				self.plll = self.plep + self.dEl[0] + self.dMu[0]
				self.mlll = self.plll.M()
				if (self.mlll> self._minmlll and self.mlll < self._maxmlll):
					return True
				else: 
					return False

			if self.decayprod == "mumu": 
				self.plll = self.plep + self.dMu[0] + self.dMu[1]
				self.mlll = self.plll.M()
				if (self.mlll> self._minmlll and self.mlll < self._maxmlll):
					return True
				else: 
					return False

			if self.decayprod == "ee": 
				self.plll = self.plep + self.dEl[0] + self.dEl[1]
				self.mlll = self.plll.M()
				if (self.mlll> self._minmlll and self.mlll < self._maxmlll):
					return True
				else: 
					return False

class Cosmicveto():
	def __init__(self, decaymode="leptonic",cosmicvetocut=0.05 ):
		self.decaymode = decaymode
		self.cosmicvetocut = cosmicvetocut

		self.seperation = -1 

	def passes(self, tree, ievt, idv):	

		if self.decaymode == "leptonic": 
			ntracks = len(tree.trackpt[ievt][idv])
			if ntracks == 2: 

				sumeta = tree.tracketa[ievt][idv][0]+tree.tracketa[ievt][idv][1]
				dphi = abs(tree.trackphi[ievt][idv][0]-tree.trackphi[ievt][idv][1])
				self.separation = np.sqrt( sumeta**2 + (np.pi -dphi)**2 )
		
				if (self.separation > self.cosmicvetocut):
					return True
				else:
					return False
			else: 
				return False







