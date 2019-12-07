import ROOT
from ROOT import * 
gROOT.LoadMacro("AtlasStyle.C")
gROOT.LoadMacro("AtlasUtils.C")
gROOT.LoadMacro("AtlasLabels.C")


#get note
def getNote(size=14):
	n = ROOT.TLatex()
	n.SetNDC()
	n.SetTextFont(43)
	n.SetTextColor(1)
	n.SetTextSize(size)
	return n

	
def drawNotes(MC_campaign,DV_type,mass,lifetime):
	a = getNote()
	b = getNote()
	c = getNote()
	d = getNote()
	e = getNote()
	ax = 0.50
	ay = 0.87
	if MC_campaign == "merged": 
		a.DrawLatex(ax,ay,'all MC campaigns')
	else:
		a.DrawLatex(ax,ay,'%s'%MC_campaign) 
	b.DrawLatex(ax,ay-0.05,'mass: %s GeV'%mass)
	c.DrawLatex(ax,ay-0.10,'lifetime: %s mm'%lifetime)
	if DV_type == "0":
		d.DrawLatex(ax,ay-0.15,'DV type: e\mu')
	else:
		d.DrawLatex(ax,ay-0.15,'DV type: \mu\mu')
	# if DV_Default == True:
	e.DrawLatex(ax,ay-0.20,'VSI')
	# else:
	# 	e.DrawLatex(ax,ay-0.20,'VSI Leptons')
	ATLASLabel(0.25,0.87,"Internal")




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