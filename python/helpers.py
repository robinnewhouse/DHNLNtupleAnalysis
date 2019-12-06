import ROOT

class Leptons(): 
	def __init__(self, decaymode="leptonic"):
		self.decaymode = decaymode
		self.lepVec = []
		self.lepIndex = []

	def getMuons(self, tree, ievt, idv):
		
		if self.decaymode == "leptonic": 
			self.ntracks = len(tree.trackpt[ievt][idv])
			# self.muVec = []
			# self.muIndex = []

			for itr in xrange(self.ntracks):
				lepVec = ROOT.TLorentzVector()
				if (tree.muonindex[ievt][idv][itr] >= 0):
					lepVec.SetPtEtaPhiE(tree.trackpt[ievt][idv][itr],tree.tracketa[ievt][idv][itr],tree.trackphi[ievt][idv][itr],tree.tracke[ievt][idv][itr])
					self.lepVec.append(lepVec)
					self.lepIndex.append(tree.muonindex[ievt][idv][itr])
			# 	else: 
			# 		lepVec.SetPtEtaPhiE(0,0,0,0)
			# 		self.muVec.append(lepVec)

			# if self.ntracks == 0: 
			# 	self.muVec.append(ROOT.TLorentzVector(0,0,0,0))

			# return self.muVec
	
	def getElectrons(self, tree, ievt, idv):
		
		if self.decaymode == "leptonic": 
			self.ntracks = len(tree.trackpt[ievt][idv])
			# self.elVec = []
			# self.elIndex = []

			for itr in xrange(self.ntracks):
				lepVec = ROOT.TLorentzVector()

				if (tree.elindex[ievt][idv][itr] >= 0):
					lepVec.SetPtEtaPhiE(tree.trackpt[ievt][idv][itr],tree.tracketa[ievt][idv][itr],tree.trackphi[ievt][idv][itr],tree.tracke[ievt][idv][itr])
					self.lepVec.append(lepVec)
					self.lepIndex.append(tree.elindex[ievt][idv][itr])
				# else: 
				# 	lepVec.SetPtEtaPhiE(0,0,0,0)
				# 	self.elVec.append(lepVec)

			# return self.elVec

	def getTracks(self, tree, ievt, idv):
		
		if self.decaymode == "leptonic": 
			self.ntracks = len(tree.trackpt[ievt][idv])
			# self.lepVec = []

			for itr in xrange(self.ntracks):
				vec = ROOT.TLorentzVector()
				vec.SetPtEtaPhiE(tree.trackpt[ievt][idv][itr],tree.tracketa[ievt][idv][itr],tree.trackphi[ievt][idv][itr],tree.tracke[ievt][idv][itr])
				self.lepVec.append(vec)
				# else: 
				# 	lepVec.SetPtEtaPhiE(0,0,0,0)
				# 	self.elVec.append(lepVec)

			# return self.lepVec