import ROOT
from ROOT import * 
import numpy as np
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
	e.DrawLatex(ax,ay-0.20,'VSI Lep')
	# else:
	# 	e.DrawLatex(ax,ay-0.20,'VSI Leptons')
	ATLASLabel(0.25,0.87,"Internal")

class Event(): 
	def __init__(self, tree, ievt, idv):
		self.tree = tree
		self.ievt = ievt
		self.idv = idv 



class Tracks(): 
	def __init__(self):
		self.lepVec = []
		self.lepIndex = []

	def getMuons(self, evt):
		self.evt = evt 
		self.ntracks = len(self.evt.tree.trackpt[self.evt.ievt][self.evt.idv])

		for itr in xrange(self.ntracks):
			lepVec = ROOT.TLorentzVector()
			if (self.evt.tree.trk_muonindex[self.evt.ievt][self.evt.idv][itr] >= 0): #matched muon!
				
				# find position of muon in the muon container that is matched to the sec vtx track (works for calibrated and uncalibrated containers) 
				muon_index = np.where(self.evt.tree.muonindex[self.evt.ievt] == self.evt.tree.trk_muonindex[self.evt.ievt][self.evt.idv][itr])[0][0]
				

				# use track quantities
				# pt = self.evt.tree.trackpt[self.evt.ievt][self.evt.idv][itr]
				# eta = self.evt.tree.tracketa[self.evt.ievt][self.evt.idv][itr]
				# phi = self.evt.tree.trackphi[self.evt.ievt][self.evt.idv][itr]
				# E = self.evt.tree.tracke[self.evt.ievt][self.evt.idv][itr]
				# lepVec.SetPtEtaPhiE(pt,eta, phi, E)

				# use calibrated muon quantities
				pt = self.evt.tree.muonpt[self.evt.ievt][muon_index]
				eta = self.evt.tree.muoneta[self.evt.ievt][muon_index]
				phi = self.evt.tree.muonphi[self.evt.ievt][muon_index]
				M = self.evt.tree.muonmass[self.evt.ievt][muon_index]
				lepVec.SetPtEtaPhiM(pt,eta, phi, M)
			

					
				self.lepVec.append(lepVec)
				self.lepIndex.append(muon_index)
	



	
	def getElectrons(self, evt):
		self.evt = evt 
		self.ntracks = len(self.evt.tree.trackpt[self.evt.ievt][self.evt.idv])

		for itr in xrange(self.ntracks):
			lepVec = ROOT.TLorentzVector()
			
			if (self.evt.tree.trk_elindex[self.evt.ievt][self.evt.idv][itr] >= 0): #matched electron!
				# find position of electron in the electron container that is matched to the sec vtx track (works for calibrated and uncalibrated containers)
				el_index = np.where(self.evt.tree.elindex[self.evt.ievt] == self.evt.tree.trk_elindex[self.evt.ievt][self.evt.idv][itr])[0][0]

				# use track quantities
				# pt = self.evt.tree.trackpt[self.evt.ievt][self.evt.idv][itr]
				# eta = self.evt.tree.tracketa[self.evt.ievt][self.evt.idv][itr]
				# phi = self.evt.tree.trackphi[self.evt.ievt][self.evt.idv][itr]
				# E = self.evt.tree.tracke[self.evt.ievt][self.evt.idv][itr]
				# lepVec.SetPtEtaPhiE(pt, eta, phi, E)

				# use calibrated electron quantities
				pt = self.evt.tree.elpt[self.evt.ievt][el_index]
				eta = self.evt.tree.eleta[self.evt.ievt][el_index]
				phi = self.evt.tree.elphi[self.evt.ievt][el_index]
				M = self.evt.tree.elmass[self.evt.ievt][el_index]
				lepVec.SetPtEtaPhiM(pt,eta, phi, M)


				self.lepVec.append(lepVec)
				self.lepIndex.append(el_index)
	



	def getTracks(self, evt):
		self.evt = evt
		self.ntracks = len(self.evt.tree.trackpt[self.evt.ievt][self.evt.idv])

		for itr in xrange(self.ntracks):
			trkvec = ROOT.TLorentzVector()
			pt = self.evt.tree.trackpt[self.evt.ievt][self.evt.idv][itr]
			eta = self.evt.tree.tracketa[self.evt.ievt][self.evt.idv][itr]
			phi = self.evt.tree.trackphi[self.evt.ievt][self.evt.idv][itr]
			E = self.evt.tree.tracke[self.evt.ievt][self.evt.idv][itr]
			trkvec.SetPtEtaPhiE(pt, eta, phi, E)

			self.lepVec.append(trkvec)
			

	













