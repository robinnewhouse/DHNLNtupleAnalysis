import ROOT

class Leptons(): 
	def __init__(self, decaymode="leptonic"):
		self.decaymode = decaymode


	def getMuons(self, tree, ievt, idv):
		
		if self.decaymode == "leptonic": 
			self.ntracks = len(tree.trackpt[ievt][idv])
			self.muVec = []

			for itr in xrange(self.ntracks):
				lepVec = ROOT.TLorentzVector()
				
				if (tree.muonindex[ievt][idv][itr] >= 0):
					lepVec.SetPtEtaPhiE(tree.trackpt[ievt][idv][itr],tree.tracketa[ievt][idv][itr],tree.trackphi[ievt][idv][itr],tree.tracke[ievt][idv][itr])
					self.muVec.append(lepVec)

			return self.muVec
	
	def getElectrons(self, tree, ievt, idv):
		
		if self.decaymode == "leptonic": 
			self.ntracks = len(tree.trackpt[ievt][idv])
			self.elVec = []

			for itr in xrange(self.ntracks):
				lepVec = ROOT.TLorentzVector()

				if (tree.elindex[ievt][idv][itr] >= 0):
					lepVec.SetPtEtaPhiE(tree.trackpt[ievt][idv][itr],tree.tracketa[ievt][idv][itr],tree.trackphi[ievt][idv][itr],tree.tracke[ievt][idv][itr])
					self.elVec.appned(lepVec)

			return self.elVec