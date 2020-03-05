import ROOT
from ROOT import * 
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import sys
import ast
import os
import re
import subprocess
import urlparse



import logging
# logging.captureWarnings(True)
msgfmt = '%(asctime)s %(levelname)-7s %(name)-35s %(message)s'
datefmt = '%H:%M:%S'

def getLogger(name = None, level = logging.INFO):
    logger = logging.getLogger(name)
    try:
        import coloredlogs
        coloredlogs.install(logger = logger, level = level, fmt = msgfmt, datefmt = datefmt)
    except ImportError:
        logging.basicConfig(format = msgfmt, datefmt = datefmt)
        logger.setLevel(level)
    return logger
logger = getLogger('dHNLAnalysis')



#get note
def getNote(size=14):
	n = ROOT.TLatex()
	n.SetNDC()
	n.SetTextFont(43)
	n.SetTextColor(1)
	n.SetTextSize(size)
	return n

	
def drawNotes(MC_campaign,DV_type,mass,lifetime,plepton,VtxConfig):
	a = getNote()
	b = getNote()
	c = getNote()
	d = getNote()
	e = getNote()
	f = getNote()
	ax = 0.50
	ay = 0.87

	if mass == "-1" or lifetime == "-1": 
		a.DrawLatex(ax,ay,'%s'%MC_campaign) 
		if plepton == "muon":
			b.DrawLatex(ax,ay-0.05,'Prompt muon')
		if plepton == "electron":
			b.DrawLatex(ax,ay-0.05,'Prompt electron')
		c.DrawLatex(ax,ay-0.10,'DV type: 1-lepton')
		d.DrawLatex(ax,ay-0.15,'%s'%(VtxConfig))

	else: 
		a.DrawLatex(ax,ay,'%s'%MC_campaign) 
		b.DrawLatex(ax,ay-0.05,'mass: %s GeV'%mass)
		c.DrawLatex(ax,ay-0.10,'lifetime: %s mm'%lifetime)
		if plepton == "muon":
			d.DrawLatex(ax,ay-0.15,'Prompt muon')
		if plepton == "electron":
			d.DrawLatex(ax,ay-0.15,'Prompt electron')
		if DV_type == "mumu":
			e.DrawLatex(ax,ay-0.20,'DV type: \mu\mu\\nu')
		if DV_type == "emu":
			e.DrawLatex(ax,ay-0.20,'DV type: e\mu\\nu')
		if DV_type == "ee":
			e.DrawLatex(ax,ay-0.20,'DV type: ee\\nu')
		f.DrawLatex(ax,ay-0.25,'%s'%(VtxConfig))
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
		# print "number of tracks: ", self.ntracks
		for itr in xrange(self.ntracks):
			lepVec = ROOT.TLorentzVector()
			if (self.evt.tree.trk_muonindex[self.evt.ievt][self.evt.idv][itr] >= 0): #matched muon!
				# find position of muon in the muon container that is matched to the sec vtx track (works for calibrated and uncalibrated containers)
				if len(self.evt.tree.muonindex[self.evt.ievt]) > 0: 
					muon_index = np.where(self.evt.tree.muonindex[self.evt.ievt] == self.evt.tree.trk_muonindex[self.evt.ievt][self.evt.idv][itr])[0][0]
					# print "muon index: ", muon_index
					# print  "track index", self.evt.tree.trk_muonindex[self.evt.ievt][self.evt.idv][itr]

					# use track quantities
					# pt = self.evt.tree.trackpt[self.evt.ievt][self.evt.idv][itr]
					# eta = self.evt.tree.tracketa[self.evt.ievt][self.evt.idv][itr]
					# phi = self.evt.tree.trackphi[self.evt.ievt][self.evt.idv][itr]
					# E = self.evt.tree.tracke[self.evt.ievt][self.evt.idv][itr]
					# lepVec.SetPtEtaPhiE(pt,eta, phi, E)

					# use calibrated muon quantities
					pt = self.evt.tree.muonpt[self.evt.ievt][muon_index]
					# print "mu pt", pt
					eta = self.evt.tree.muoneta[self.evt.ievt][muon_index]
					phi = self.evt.tree.muonphi[self.evt.ievt][muon_index]
					M = self.evt.tree.muonmass[self.evt.ievt][muon_index]
					lepVec.SetPtEtaPhiM(pt,eta, phi, M)
				

						
					self.lepVec.append(lepVec)
					self.lepIndex.append(muon_index)
				else:
					continue
	



	
	def getElectrons(self, evt):
		self.evt = evt 
		self.ntracks = len(self.evt.tree.trackpt[self.evt.ievt][self.evt.idv])

		for itr in xrange(self.ntracks):
			lepVec = ROOT.TLorentzVector()
			
			if (self.evt.tree.trk_elindex[self.evt.ievt][self.evt.idv][itr] >= 0): #matched electron!
				# find position of electron in the electron container that is matched to the sec vtx track (works for calibrated and uncalibrated containers)
				if len(self.evt.tree.elindex[self.evt.ievt]) > 0: 
					el_index = np.where(self.evt.tree.elindex[self.evt.ievt] == self.evt.tree.trk_elindex[self.evt.ievt][self.evt.idv][itr])[0][0]
					# print "el_index", el_index
					# print "track index", self.evt.tree.trk_elindex[self.evt.ievt][self.evt.idv][itr]
					# use track quantities
					# pt = self.evt.tree.trackpt[self.evt.ievt][self.evt.idv][itr]
					# eta = self.evt.tree.tracketa[self.evt.ievt][self.evt.idv][itr]
					# phi = self.evt.tree.trackphi[self.evt.ievt][self.evt.idv][itr]
					# E = self.evt.tree.tracke[self.evt.ievt][self.evt.idv][itr]
					# lepVec.SetPtEtaPhiE(pt, eta, phi, E)

					# use calibrated electron quantities
					pt = self.evt.tree.elpt[self.evt.ievt][el_index]
					# print "el pt", pt
					eta = self.evt.tree.eleta[self.evt.ievt][el_index]
					phi = self.evt.tree.elphi[self.evt.ievt][el_index]
					M = self.evt.tree.elmass[self.evt.ievt][el_index]
					lepVec.SetPtEtaPhiM(pt,eta, phi, M)


					self.lepVec.append(lepVec)
					self.lepIndex.append(el_index)
				else: 
					continue 
	



	def getTracks(self, evt):
		self.evt = evt
		self.ntracks = len(self.evt.tree.trackpt[self.evt.ievt][self.evt.idv])
	

		for itr in xrange(self.ntracks):
			trkvec = ROOT.TLorentzVector()
			pt = self.evt.tree.trackpt[self.evt.ievt][self.evt.idv][itr]
			eta = self.evt.tree.tracketa[self.evt.ievt][self.evt.idv][itr]
			phi = self.evt.tree.trackphi[self.evt.ievt][self.evt.idv][itr]
			M = self.evt.tree.trackmass[self.evt.ievt][self.evt.idv][itr]
			
			trkvec.SetPtEtaPhiM(pt, eta, phi, M)

			self.lepVec.append(trkvec)

			

	













