# from ROOT import*
import ROOT
import numpy as np
import helpers

muIDdict =	{
  2: "tight"
}


class Trigger():
	def __init__(self, plepton, trigger = None):
		self.plepton = plepton 
		self.trigger = trigger

	def passes(self, tree):
		if self.plepton == "muon":
			if self.trigger == "triggerpassmu26":
				if tree.triggerpassmu26 == True:
					return True 
		elif self.plepton == "el":
			return False
		else:
			return False


class Filter():
	def __init__(self, _filter):
		self.filter = _filter 
		
	def passes(self,tree):
		if self.filter == "mu-mu":
			if tree.mumufilter == True:
				return True
		if self.filter == "mu-el":
			if tree.mumufilter == True:
				return True
		if self.filter == "el-mu":
			if tree.mumufilter == True:
				return True
		if self.filter == "el-el":
			if tree.mumufilter == True:
				return True

class Plepton():
	def __init__(self, lepton, quality="tight", _mindR=0.05):
		self.quality = quality 
		self.lepton = lepton
		self._mindR = _mindR
		self.plepVec = None

	def passes(self, tree, ievt):
		ndv = len(tree.dvx[ievt])	
		dmuVec = []
		delVec = []
	
		for idv in xrange(ndv): 
			dmuVec.append(helpers.Leptons().getMuons(tree, ievt, idv))
			delVec.append(helpers.Leptons().getElectrons(tree, ievt, idv))

		print dmuVec

		if self.lepton == "muon":
			nmuons = len(tree.muonpt[ievt])
			ndmuons = len(dmuVec)
			overlap = False

			for imu in xrange(nmuons): 
				pmuVec = ROOT.TLorentzVector()
				pmuVec.SetPtEtaPhiM(tree.muonpt[ievt][imu],tree.muoneta[ievt][imu],tree.muonphi[ievt][imu],tree.muonmass[ievt][imu])
				# if tree.muonpassPfilter[ievt]:
				if imu == 0: #temporary until we have this variable!!!!
					for idvmu in xrange(ndmuons): 
						dR = np.sqrt((dmuVec[idvmu].Eta() - pmuVec.Eta())**2 + (dmuVec[idvmu].Phi() - pmuVec.Phi())**2)
						if dR > self._mindR: 
							plepVec = pmuVec
							return True
						else: 
							return False
				else: 
					return False
	def getPlep():
		return self.plepVec









			# if quality == "tight":
			# 	if tree.muTight == True:
			# 		return True
			# if quality == "medium":
			# 	if tree.muMedium == True:
			# 		return True
			# if quality == "loose":	
			# 	if tree.muLoose == True:
			# 		return True

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

	def passes(self,tree,ievt):
		if len(tree.dvx[ievt]) > 0:
			self.rdv = np.sqrt(tree.dvx[ievt][0]**2 + tree.dvy[ievt][0]**2)

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

			for itr in xrange(self.ntracks): #loop over tracks

				if (tree.muonindex[ievt][idv][itr] >= 0): # count muons 
					imu = tree.muonindex[ievt][idv][itr]
					if tree.muontype[ievt][imu]==0: # only combined muons!! (want to check this cut I think)
						self.nmu = self.nmu + 1

				if (tree.elindex[ievt][idv][itr] >= 0): # count electrons 
					self.nel = self.nel + 1		
			# end track loop


			if self.decayprod == "emu": 
				if self.nel == 1 and self.nmu == 1: 
					return True
				else:
					return False

			elif self.decayprod == "mumu":
				if self.nmu == 2: 
					return True
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

			for itr in xrange(self.ntracks):
				if (tree.muonindex[ievt][idv][itr] >= 0):
					self.DVmuons.append(tree.muonindex[ievt][idv][itr])

				if (tree.elindex[ievt][idv][itr] >= 0):
					self.DVelectrons.append(tree.muonindex[ievt][idv][itr])

			self.nmu_tight = 0
			self.nel_tight = 0

			for imu in self.DVmuons:
				if tree.muonistight[ievt][imu] == True: 
					self.nmu_tight = self.nmu_tight + 1

			# for iel in self.DVelectrons: 
			# 	if tree.eltight[ievt][imu] == True: 
			# 		self.nel_tight = self.nel_tight + 1

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

class Masscut():
	def __init__(self, decaymode="leptonic"):
		self.decaymode = decaymode

	# def mlllpassses(self, plepton, dlepton1, dlepton2):


	def mDVpassses(self, tree, ievt, idv, masscut=4):
		if (tree.dvmass[ievt][idv] > masscut):
			return True
		else: 
			return False



		




		










