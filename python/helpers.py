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
		self.HNL_vec = ROOT.TLorentzVector()
		self.dNu_vec =  ROOT.TLorentzVector()
		self.trkVec = []
		self.truth_dvx = -1 
		self.truth_dvy = -1 
		self.truth_dvz = -1 
		self.truth_pvx = -1 
		self.truth_pvy = -1 
		self.truth_pvz = -1 
		self.truth_dvr = -1 
		self.W_vec =  ROOT.TLorentzVector()
		self.plep_vec =  ROOT.TLorentzVector()
		self.mhnl = -1
		self.HNL_pdgID = 50

	def getTruthParticles(self, evt): 
		ntruthDV = len(evt.tree.truth_parent_pdgId[evt.ievt])

		for idvtru in xrange(ntruthDV):
			if abs(evt.tree.truth_parent_pdgId[evt.ievt][idvtru]) == 50:  # get the DV!
				if len(evt.tree.truth_outP_pdgId[evt.ievt][idvtru]) == 3:

					self.truth_dvx = evt.tree.truth_x[evt.ievt][idvtru]
					self.truth_dvy = evt.tree.truth_y[evt.ievt][idvtru]
					self.truth_dvz = evt.tree.truth_z[evt.ievt][idvtru] 
					self.truth_dvr = np.sqrt(self.truth_dvx**2 + self.truth_dvy**2)
					
					trkVec0 =  ROOT.TLorentzVector()
					trkVec1 =  ROOT.TLorentzVector()
					nu_vec =  ROOT.TLorentzVector()


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

					self.trkVec.append(trkVec0)
					self.trkVec.append(trkVec1)


					nu_vec.SetPtEtaPhiM(evt.tree.truth_outP_pt[evt.ievt][idvtru][2],
										evt.tree.truth_outP_eta[evt.ievt][idvtru][2],
										evt.tree.truth_outP_phi[evt.ievt][idvtru][2],
										evt.tree.truth_outP_m[evt.ievt][idvtru][2]
										)

					self.HNL_vec.SetPtEtaPhiM(evt.tree.truth_parent_pt[evt.ievt][idvtru],
										evt.tree.truth_parent_eta[evt.ievt][idvtru],
										evt.tree.truth_parent_phi[evt.ievt][idvtru],
										evt.tree.truth_parent_m[evt.ievt][idvtru])




			if abs(evt.tree.truth_parent_pdgId[evt.ievt][idvtru]) == 24: # get the PV!
				if len(evt.tree.truth_outP_pdgId[evt.ievt][idvtru]) == 2:

					self.truth_pvx = evt.tree.truth_x[evt.ievt][idvtru]
					self.truth_pvy = evt.tree.truth_y[evt.ievt][idvtru]
					self.truth_pvz = evt.tree.truth_z[evt.ievt][idvtru] 

					self.plep_vec.SetPtEtaPhiM(  evt.tree.truth_outP_pt[evt.ievt][idvtru][0],
											evt.tree.truth_outP_eta[evt.ievt][idvtru][0],
											evt.tree.truth_outP_phi[evt.ievt][idvtru][0],
											evt.tree.truth_outP_m[evt.ievt][idvtru][0]
											)
					self.W_vec.SetPtEtaPhiM( evt.tree.truth_parent_pt[evt.ievt][idvtru],
										evt.tree.truth_parent_eta[evt.ievt][idvtru],
										evt.tree.truth_parent_phi[evt.ievt][idvtru],
										evt.tree.truth_parent_m[evt.ievt][idvtru]
										)

		# Mhnl = selections.Mhnl(evt=evt, plep=self.plep_vec, trks =self.trkVec ) # bug need to adapt mhnl to run with truth vertices 
		# self.mhnl = Mhnl.mhnl




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
					pt = self.evt.tree.trackpt[self.evt.ievt][self.evt.idv][itr]
					eta = self.evt.tree.tracketa[self.evt.ievt][self.evt.idv][itr]
					phi = self.evt.tree.trackphi[self.evt.ievt][self.evt.idv][itr]
					M = self.evt.tree.trackmass[self.evt.ievt][self.evt.idv][itr]
					# E = self.evt.tree.tracke[self.evt.ievt][self.evt.idv][itr]
					lepVec.SetPtEtaPhiM(pt,eta, phi, M)

					# use calibrated muon quantities
					# pt = self.evt.tree.muonpt[self.evt.ievt][muon_index]
					# eta = self.evt.tree.muoneta[self.evt.ievt][muon_index]
					# phi = self.evt.tree.muonphi[self.evt.ievt][muon_index]
					# M = self.evt.tree.muonmass[self.evt.ievt][muon_index]
					# lepVec.SetPtEtaPhiM(pt,eta, phi, M)


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
					

					if (self.evt.tree.trk_muonindex[self.evt.ievt][self.evt.idv][itr] >= 0): # remove electrons that are also matched to muons!
						if len(self.evt.tree.muonindex[self.evt.ievt]) > 0: 
							muon_index = np.where(self.evt.tree.muonindex[self.evt.ievt] == self.evt.tree.trk_muonindex[self.evt.ievt][self.evt.idv][itr])[0][0]
							# print muon_index
							# print "track is matched to both muon and electron!"
							continue

					#use track quantities
					pt = self.evt.tree.trackpt[self.evt.ievt][self.evt.idv][itr]
					eta = self.evt.tree.tracketa[self.evt.ievt][self.evt.idv][itr]
					phi = self.evt.tree.trackphi[self.evt.ievt][self.evt.idv][itr]
					M = self.evt.tree.trackmass[self.evt.ievt][self.evt.idv][itr]
					# E = self.evt.tree.tracke[self.evt.ievt][self.evt.idv][itr]
					lepVec.SetPtEtaPhiM(pt, eta, phi, M)

					# use calibrated electron quantities
					# pt = self.evt.tree.elpt[self.evt.ievt][el_index]
					# eta = self.evt.tree.eleta[self.evt.ievt][el_index]
					# phi = self.evt.tree.elphi[self.evt.ievt][el_index]
					# M = self.evt.tree.elmass[self.evt.ievt][el_index]
					# lepVec.SetPtEtaPhiM(pt,eta, phi, M)
					
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

		self.MC_campaign =""
		self.ctau_str = ""
		self.mass_str = ""

		if "lt1dd" in file or "1mm" in file:
			self.ctau = 1.0
			self.ctau_str = "1mm"
		elif "lt10dd" in file  or "10mm" in file:
			self.ctau = 10.0
			self.ctau_str = "10mm"
		elif "lt100dd" in file  or "100mm" in file:
			self.ctau = 100.0
			self.ctau_str = "100mm"

		if "3G" in file:
			self.mass = 3.0
			self.mass_str = "3G"
		elif "4G" in file:
			self.mass = 4.0
			self.mass_str = "4G"
		elif "4p5G" in file:
			self.mass = 4.5
			self.mass_str = "4p5G"
		elif "5G" in file:
			self.mass = 5.0
			self.mass_str = "5G"
		elif "7p5G" in file:
			self.mass = 7.5
			self.mass_str = "7p5G"
		elif "10G" in file:
			self.mass = 10.0
			self.mass_str = "10G" 
		elif "12p5G" in file:
			self.mass = 12.5
			self.mass_str = "12p5G"
		elif "15G" in file:
			self.mass = 15.0
			self.mass_str = "15G"
		elif "17p5G" in file:
			self.mass = 17.5
			self.mass_str = "17p5G"
		elif "20G" in file:
			self.mass = 20.0
			self.mass_str = "20G"
		

		if "r10740" in file or "mc16a" in file:
			self.MC_campaign = "mc16a"
		if "r10739" in file or "mc16d" in file:
			self.MC_campaign = "mc16d"
		if "r10790" in file or "mc16e" in file:
			self.MC_campaign = "mc16e"
		


		self.Output_filename = "histograms%s_%s_%s_%s.root"%(self.MC_campaign, self.mass_str, self.ctau_str, channel) 
		




			

	













