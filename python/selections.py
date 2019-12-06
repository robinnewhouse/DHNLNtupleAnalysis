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
			


			# for itr in xrange(self.ntracks): #loop over tracks

			# 	if (tree.muonindex[ievt][idv][itr] >= 0): # count muons 
			# 		imu = tree.muonindex[ievt][idv][itr]
			# 		if tree.muontype[ievt][imu]==0: # only combined muons!! (want to check this cut I think)
			# 			self.nmu = self.nmu + 1

			# 	if (tree.elindex[ievt][idv][itr] >= 0): # count electrons 
			# 		self.nel = self.nel + 1		
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

			# muVec = helpers.Leptons().getMuons(tree, ievt, idv)
			# elVec = helpers.Leptons().getElectrons(tree, ievt, idv)

			# for itr in xrange(self.ntracks):
			# 	if (tree.muonindex[ievt][idv][itr] >= 0):
			# 		self.DVmuons.append(tree.muonindex[ievt][idv][itr])

			# 	if (tree.elindex[ievt][idv][itr] >= 0):
			# 		self.DVelectrons.append(tree.muonindex[ievt][idv][itr])

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
				# index = self.DVmuons[imu]
				index = muons.lepIndex[imu]
				# print index
				# print tree.tightmu[ievt][index]
				# if tree.tightmu[ievt][index] == True: 
				if tree.tightmu[ievt][index] : 
					self.nmu_tight = self.nmu_tight + 1

			# for iel in range(ndvelectrons): 
			# 	index = electrons.lepIndex[imu]

			# 	if tree.eltight[ievt][index] == True: 
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


	def passes(self, tree, ievt, idv, masscut=4):
		if (tree.dvmass[ievt][idv] > masscut):
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

# class Cosmicveto():
# 	def __init__(self, decayprod, decaymode="leptonic"):
# 		self.decaymode = decaymode
# 		self.decayprod = decayprod
# 		self.mlll = -1

# 	def passes(self, tree, ievt, idv):	
	
		




# 			sumeta = tracketa[ievt][idv][0]+tracketa[ievt][idv][1]
# 					dphi = abs(trackphi[ievt][idv][0]-trackphi[ievt][idv][1])
# 					separation = np.sqrt( sumeta**2 + (np.pi -dphi)**2 )







