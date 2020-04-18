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
import atlas_style, selections

import logging
# logging.captureWarnings(True)
msgfmt = '%(asctime)s %(levelname)-7s %(name)-35s %(message)s'
datefmt = '%H:%M:%S'

def getLogger(name = None, level = logging.DEBUG):
    logger = logging.getLogger(name)
    try:
        import coloredlogs
        coloredlogs.install(logger = logger, level = level, fmt = msgfmt, datefmt = datefmt)
    except ImportError:
        logging.basicConfig(format = msgfmt, datefmt = datefmt)
        logger.setLevel(level)
    return logger
logger = getLogger('dHNLAnalysis.helpers')



class Event(): 
	def __init__(self, tree, ievt, idv=None, mass=1.0, ctau=1.0):
		self.tree = tree
		self.ievt = ievt
		self.idv = idv 
		LNV = False

		if tree.isData: # you are running on data
			self.weight = 1
		else: # you are running on MC file
			if mass == -1 or ctau == -1: # MC weighting error
				logger.warning("Can't determine the mass and lifetime of signal sample. MC weight will be set to 1!!")
				self.weight = 1
			else: 
				mW = 80.379 # mass of W boson in GeV
				U2Gronau=4.49e-12*3e8*mass**(-5.19)/(ctau/1000) #LNC prediction	
				if(LNV):  U2=0.5*U2
				else: U2 =  U2Gronau

				xsec = 20.6e6*U2*((1-(mass/mW)**2)**2)*(1+(mass**2)/(2*mW**2))#in fb
				self.weight = 1*xsec/tree.allEvt #scale to 1 fb^-1  of luminosity



class Truth(): 
	def __init__(self):
		self.HNL_vec = []
		self.dMuon_vec = []
		self.dEl_vec = []
		self.dNu_vec = []
		self.HNL_pdgID = 50

	def getTruthParticles(self, evt): 
		print "------------------------------------------------------------------------------"
		print "event number: ", evt.tree.evtNum[evt.ievt]
		print ""
		print ""
		print "#########"
		print "TRUTH "
		print "#########"
		print ""
		print ""
		ntruthDV = len(evt.tree.truth_parent_pdgId[evt.ievt])
		trkVec = []
		for idvtru in xrange(ntruthDV):
			if abs(evt.tree.truth_parent_pdgId[evt.ievt][idvtru]) == 50: 
				if len(evt.tree.truth_outP_pdgId[evt.ievt][idvtru]) == 3:
					# print evt.tree.truth_parent_pdgId[evt.ievt][idvtru]
					print "DV:     ", evt.tree.truth_x[evt.ievt][idvtru], "     ",evt.tree.truth_y[evt.ievt][idvtru], "     ", evt.tree.truth_z[evt.ievt][idvtru] 
					print "HNL:     ", evt.tree.truth_parent_pt[evt.ievt][idvtru], "     ", evt.tree.truth_parent_eta[evt.ievt][idvtru], "     ", evt.tree.truth_parent_phi[evt.ievt][idvtru], "     ",evt.tree.truth_parent_m[evt.ievt][idvtru]
					print "d lep 0: ", evt.tree.truth_outP_pt[evt.ievt][idvtru][0], "     ", evt.tree.truth_outP_eta[evt.ievt][idvtru][0], "     ", evt.tree.truth_outP_phi[evt.ievt][idvtru][0], "     ",evt.tree.truth_outP_m[evt.ievt][idvtru][0]
					print "d lep 1: ", evt.tree.truth_outP_pt[evt.ievt][idvtru][1], "     ", evt.tree.truth_outP_eta[evt.ievt][idvtru][1], "     ", evt.tree.truth_outP_phi[evt.ievt][idvtru][1], "     ",evt.tree.truth_outP_m[evt.ievt][idvtru][1]
					print "nu:      ", evt.tree.truth_outP_pt[evt.ievt][idvtru][2], "      ", evt.tree.truth_outP_eta[evt.ievt][idvtru][2], "     ", evt.tree.truth_outP_phi[evt.ievt][idvtru][2], "     ",evt.tree.truth_outP_m[evt.ievt][idvtru][2]
					
					trkVec0 =  ROOT.TLorentzVector()
					trkVec1 =  ROOT.TLorentzVector()

					trkVec0.SetPtEtaPhiM(evt.tree.truth_outP_pt[evt.ievt][idvtru][0],
										evt.tree.truth_outP_eta[evt.ievt][idvtru][0],
										evt.tree.truth_outP_phi[evt.ievt][idvtru][0],
										evt.tree.truth_outP_m[evt.ievt][idvtru][0]
										)
					trkVec1.SetPtEtaPhiM(evt.tree.truth_outP_pt[evt.ievt][idvtru][1],
										evt.tree.truth_outP_eta[evt.ievt][idvtru][1],
										evt.tree.truth_outP_phi[evt.ievt][idvtru][1],
										evt.tree.truth_outP_m[evt.ievt][idvtru][1]
										)

					trkVec.append(trkVec0)
					trkVec.append(trkVec1)

			if abs(evt.tree.truth_parent_pdgId[evt.ievt][idvtru]) == 24:
				if len(evt.tree.truth_outP_pdgId[evt.ievt][idvtru]) == 2:
					plep_vec = ROOT.TLorentzVector()
					plep_vec.SetPtEtaPhiM(evt.tree.truth_outP_pt[evt.ievt][idvtru][0],
										evt.tree.truth_outP_eta[evt.ievt][idvtru][0],
										evt.tree.truth_outP_phi[evt.ievt][idvtru][0],
										evt.tree.truth_outP_m[evt.ievt][idvtru][0])

					print "PV:    ", evt.tree.truth_x[evt.ievt][idvtru], "     ",evt.tree.truth_y[evt.ievt][idvtru], "     ", evt.tree.truth_z[evt.ievt][idvtru] 
					print "W:     ", evt.tree.truth_parent_pt[evt.ievt][idvtru], "     ", evt.tree.truth_parent_eta[evt.ievt][idvtru], "     ", evt.tree.truth_parent_phi[evt.ievt][idvtru], "     ",evt.tree.truth_parent_m[evt.ievt][idvtru]
					print "p lep: ", evt.tree.truth_outP_pt[evt.ievt][idvtru][0], "     ", evt.tree.truth_outP_eta[evt.ievt][idvtru][0], "     ", evt.tree.truth_outP_phi[evt.ievt][idvtru][0], "     ",evt.tree.truth_outP_m[evt.ievt][idvtru][0]
					print ""
					print ""
				# 	# print tree.truth_parent_pdgId[ievt][idvtru]
				# 	print "DV:   ", tree.truth_x[ievt][idvtru], "     ",tree.truth_y[ievt][idvtru], "     ", tree.truth_z[ievt][idvtru] 
				# 	print "HNL:   ", tree.truth_parent_pt[ievt][idvtru], "     ", tree.truth_parent_eta[ievt][idvtru], "     ", tree.truth_parent_phi[ievt][idvtru], "     ",tree.truth_parent_m[ievt][idvtru]
				# 	print "lep 0: ", tree.truth_outP_pt[ievt][idvtru][0], "     ", tree.truth_outP_eta[ievt][idvtru][0], "     ", tree.truth_outP_phi[ievt][idvtru][0], "     ",tree.truth_outP_m[ievt][idvtru][0]
				# 	print "lep 1: ", tree.truth_outP_pt[ievt][idvtru][1], "     ", tree.truth_outP_eta[ievt][idvtru][1], "     ", tree.truth_outP_phi[ievt][idvtru][1], "     ",tree.truth_outP_m[ievt][idvtru][1]
				# 	print "nu:    ", tree.truth_outP_pt[ievt][idvtru][2], "      ", tree.truth_outP_eta[ievt][idvtru][2], "     ", tree.truth_outP_phi[ievt][idvtru][2], "     ",tree.truth_outP_m[ievt][idvtru][2]

		Mhnl = selections.Mhnl(evt=evt, plep=plep_vec, trks =trkVec )
		print ""
		print ""
		print "negative HNL mass solution: ",  Mhnl.mhnl
		print "negative solution for neutrino 4-vector: ",  Mhnl.nu_vec.Pt(),"     ",Mhnl.nu_vec.Eta(),"     ",Mhnl.nu_vec.Phi(),"     ",Mhnl.nu_vec.M()
		print""
		print "postive HNL mass solution: ",  Mhnl.mhnl2
		print "positive solution for neutrino 4-vector: ",  Mhnl.nu_vec2.Pt(),"     ",Mhnl.nu_vec2.Eta(),"     ",Mhnl.nu_vec2.Phi(),"     ",Mhnl.nu_vec.M()

					# HNLVec = ROOT.TLorentzVector()
					# HNLVec.SetPtEtaPhiM(tree.truth_parent_pt[ievt][idvtru],
					# 					tree.truth_parent_eta[ievt][idvtru],
					# 					tree.truth_parent_phi[ievt][idvtru],
					# 					tree.truth_parent_m[ievt][idvtru])
					# self.HNL_vec.append(HNLVec)
					# print HNLVec.M()
				# self.num_HNL = self.num_HNL + 1
				# print self.num_HNL



class Tracks(): 
	def __init__(self):
		self.lepVec = []
		self.lepIndex = []
		self.eta = []
		self.phi = []
		self.pt = []
		self.ntracks = -1 

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
					self.pt.append(pt)
					self.eta.append(eta)
					self.phi.append(phi)

						
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
					self.pt.append(pt)
					self.eta.append(eta)
					self.phi.append(phi)


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
			self.eta.append(eta)
			self.phi.append(phi)
			self.pt.append(pt)

class File_info():

	def __init__(self,file, channel):
		self.Output_filename = "histograms.root"
		self.mass = -1 # signal mass of HNL in GeV
		self.ctau = -1 # in mm

		MC_campaign =""
		ctau_str = ""
		mass_str = ""

		if "lt1dd" in file:
			self.ctau = 1.0
			ctau_str = "1mm"
		elif "lt10dd" in file:
			self.ctau = 10.0
			ctau_str = "10mm"
		elif "lt100dd" in file:
			self.ctau = 100.0
			ctau_str = "100mm"

		if "3G" in file:
			self.mass = 3.0
			mass_str = "3G"
		elif "4G" in file:
			self.mass = 4.0
			mass_str = "4G"
		elif "4p5G" in file:
			self.mass = 4.5
			mass_str = "4p5G"
		elif "5G" in file:
			self.mass = 5.0
			mass_str = "5G"
		elif "7p5G" in file:
			self.mass = 7.5
			mass_str = "7p5G"
		elif "10G" in file:
			self.mass = 10.0
			mass_str = "10G" 
		elif "12p5G" in file:
			self.mass = 12.5
			mass_str = "12p5G"
		elif "15G" in file:
			self.mass = 15.0
			mass_str = "15G"
		elif "17p5G" in file:
			self.mass = 17.5
			mass_str = "17p5G"
		elif "20G" in file:
			self.mass = 20.0
			mass_str = "20G"
		

		if "r10740" in file:
			MC_campaign = "mc16a"
		if "r10739" in file:
			MC_campaign = "mc16d"
		if "r10790" in file:
			MC_campaign = "mc16e"
		


		self.Output_filename = "histograms_%s_%s_%s_%s.root"%(MC_campaign, mass_str, ctau_str, channel) 
		

#get note
def getNote(size=14):
	n = ROOT.TLatex()
	n.SetNDC()
	n.SetTextFont(43)
	n.SetTextColor(1)
	n.SetTextSize(size)
	return n

	
def drawNotes(DV_type,plepton,VtxConfig):
	a = getNote()
	b = getNote()
	c = getNote()
	d = getNote()
	e = getNote()
	f = getNote()
	ax = 0.25
	ay = 0.87

	if plepton == "muon":
		a.DrawLatex(ax,ay-0.05,'Prompt muon')
	if plepton == "electron":
		a.DrawLatex(ax,ay-0.05,'Prompt electron')
	if DV_type == "mumu":
		b.DrawLatex(ax,ay-0.10,'DV type: \mu\mu\\nu')
	if DV_type == "emu":
		b.DrawLatex(ax,ay-0.10,'DV type: e\mu\\nu')
	c.DrawLatex(ax,ay-0.15,'%s'%(VtxConfig))

	# else: 
	# 	a.DrawLatex(ax,ay,'%s'%MC_campaign) 
	# 	b.DrawLatex(ax,ay-0.05,'mass: %s GeV'%mass)
	# 	c.DrawLatex(ax,ay-0.10,'lifetime: %s mm'%lifetime)
	# 	if plepton == "muon":
	# 		d.DrawLatex(ax,ay-0.15,'Prompt muon')
	# 	if plepton == "electron":
	# 		d.DrawLatex(ax,ay-0.15,'Prompt electron')
	# 	if DV_type == "mumu":
	# 		e.DrawLatex(ax,ay-0.20,'DV type: \mu\mu\\nu')
	# 	if DV_type == "emu":
	# 		e.DrawLatex(ax,ay-0.20,'DV type: e\mu\\nu')
	# 	if DV_type == "ee":
	# 		e.DrawLatex(ax,ay-0.20,'DV type: ee\\nu')
	# 	f.DrawLatex(ax,ay-0.25,'%s'%(VtxConfig))
	# 	# else:
		# 	e.DrawLatex(ax,ay-0.20,'VSI Leptons')
	atlas_style.ATLASLabel(0.25,0.87,"Internal")


def drawNotesMC(MC_campaign,Vertextype, channel,mass,lifetime):
	a = getNote()
	b = getNote()
	c = getNote()
	d = getNote()
	e = getNote()
	ax = 0.25
	ay = 0.87
	if MC_campaign == "merged": 
		a.DrawLatex(ax,ay,'all MC campaigns')
	else:
		a.DrawLatex(ax,ay,'%s'%MC_campaign) 
	b.DrawLatex(ax,ay-0.05,'mass: %s GeV'%mass)
	c.DrawLatex(ax,ay-0.10,'lifetime: %s mm'%lifetime)
	if channel == "uuu":
		d.DrawLatex(ax,ay-0.15,'channel: \mu\mu\mu')
	elif channel == "ueu":
		d.DrawLatex(ax,ay-0.15,'channel: \mue\mu')
	elif channel == "uee":
		d.DrawLatex(ax,ay-0.15,'channel: \muee')
	elif channel == "eee":
		d.DrawLatex(ax,ay-0.15,'channel: eee')
	elif channel == "eeu":
		d.DrawLatex(ax,ay-0.15,'channel: ee\mu')
	elif channel == "euu":
		d.DrawLatex(ax,ay-0.15,'channel: e\mu\mu')
	# if DV_Default == True:
	# 	e.DrawLatex(ax,ay-0.20,'VSI')
	# else:
	e.DrawLatex(ax,ay-0.20,Vertextype)
	atlas_style.ATLASLabel(0.25,0.87,"Internal")

def drawNotesData(datarun,Vertextype):
	a = getNote()
	b = getNote()
	
	ax = 0.25
	ay = 0.82

	a.DrawLatex(ax,ay,Vertextype)
	b.DrawLatex(ax,ay-0.05,datarun)
	atlas_style.ATLASLabel(0.25,0.87,"Internal")


def drawNotesVertextype(Vertextype, lumi=1):
	a = getNote()
	b = getNote()
	
	ax = 0.25
	ay = 0.82

	a.DrawLatex(ax,ay,"%s fb^{-1}"%lumi)
	b.DrawLatex(ax,ay-0.05,Vertextype)
	


def xlabelhistograms(hist): 
	if "DV_r" in hist:
		if  ("redmassvis" in hist):
			return "reduced visible mass [GeV]"
		elif  ("redmass" in hist):
			if "redmassHNL" in hist:
				return "reduced HNL mass [GeV]"
			else:
				return "reduced DV mass [GeV]"
		else: 
			return "DV r [mm]"
	if "DV_mass" in hist: 
		return "DV mass [GeV]"
	if "trk_pt" in hist:
		return "track p_{T} [GeV]"
	if "trk_eta" in hist:
		return "track \eta"
	if "trk_phi" in hist:
		return "track \phi"
	if "trk_d0" in hist:
		return "track d_{0}"
	if "mvis" in hist:
		return "Visible mass (m_{lll}) [GeV]"
	if "dpt" in hist: 
		return "\Deltap_{T} between tracks in DV [GeV]"
	if "deta" in hist: 
		return "\Delta\eta between tracks in DV"
	if "dphi" in hist: 
		return "\Delta\phi between tracks in DV"
	if "dR" in hist: 
		return "\DeltaR between tracks in DV"
	if "mtrans" in hist:
		return "m_{T} [GeV]"
	if "HNLm" in hist: 
		return "HNL mass [GeV]"
	if "HNLpt" in hist: 
		return "HNL p_{T} [GeV]"
	if "HNLphi" in hist: 
		return "HNL \phi"
	if "HNLeta" in hist: 
		return "HNL \eta"
	else: 
		return ""


def histColours(nhist): 
	if nhist== 0:
		return kAzure+6
	if nhist== 1:
		return kViolet+8
	if nhist== 2:
		return kRed
	if nhist== 3:
		return kGreen+1
	if nhist== 4:
		return kOrange -3
	else: 
		return kBlack




			

	













