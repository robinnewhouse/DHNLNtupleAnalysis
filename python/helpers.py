import ROOT
import numpy as np
import sys
import os
import time

ROOT.PyConfig.IgnoreCommandLineOptions = True
import atlas_style

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
# get a logger for own module
logger = getLogger('dHNLAnalysis.helpers')
logger_debug_level = logging.INFO # default

def get_debug_level(level):
	import logging
	# supported levels are: CRITICAL, ERROR, WARNING, INFO, DEBUG
	debug_level = None
	if level == "INFO":
		debug_level = logging.INFO
	elif level == "DEBUG":
		debug_level = logging.DEBUG
	elif level == "CRITICAL":
		debug_level = logging.CRITICAL
	elif level == "WARNING":
		debug_level = logging.WARNING
	elif level == "ERROR":
		debug_level = logging.ERROR
	return debug_level


class ReadBRdat:
    def __init__(self, filename):
        f = open(filename, 'r')
        content = f.read()

        list_rows = []
        row = ''
        for c in content:
            if c == '\n':
                list_rows.append(row)
                row = ''
            else:
                row = row + c

        list_content = list_rows[18:]

        mass = []
        U_e = []
        U_mu = []
        U_tau = []
        BR_uu = []
        BR_mue = []
        BR_ee = []

        number_t = []
        for row in list_content:
            count = 0
            entry = ''
            for c in row:
                if c != '\t':
                    entry = entry + c
                if c == '\t':
                    if count == 0:
                        mass.append(entry)
                    if count == 1:
                        U_e.append(entry)
                    if count == 2:
                        U_mu.append(entry)
                    if count == 3:
                        U_tau.append(entry)
                    if count == 14:
                        BR_ee.append(entry)
                    if count == 15:
                        BR_mue.append(entry)

                    count = count + 1
                    entry = ''
            BR_uu.append(entry)
            number_t.append(count)

        self.M = []
        self.BRuuu = []
        self.BRuue = []
        for i in range(0, len(BR_uu)):
            if U_e[i] == '0.' and U_tau[i] == '0.' and U_mu[i] != '0.':
                self.M.append(float(mass[i]))
                self.BRuuu.append(float(BR_uu[i]))
                self.BRuue.append(float(BR_mue[i]))

        self.M_el = []
        self.BReee = []
        self.BReeu = []
        for i in range(0, len(BR_uu)):
            if U_e[i] != '0.' and U_tau[i] == '0.' and U_mu[i] == '0.':
                self.M_el.append(float(mass[i]))
                self.BReee.append(float(BR_ee[i]))
                self.BReeu.append(float(BR_mue[i]))

    def get_BR(self, channel, mass):

        i = -1
        for m in self.M:
            i = i + 1
            if m == mass:
                if channel == 'uuu':
                    return self.BRuuu[i]
                if channel == 'uue':
                    return self.BRuue[i]
                if channel == 'eee':
                    return self.BReee[i]
                if channel == 'eeu':
                    return self.BReeu[i]


def hnl_xsec(br, U2, mass):
	mW = 80.379  # mass of W boson in GeV
	xsec = br * 20.6e6 * U2 * ((1 - (mass / mW) ** 2) ** 2) * (1 + (mass ** 2) / (2 * mW ** 2))  # in fb
	return xsec


def get_mass_lt_weight(tree, lnc_plus_lnv=False):
	"""
	Calculates the weight of the event based on the Gronau parametrization
	https://journals.aps.org/prd/abstract/10.1103/PhysRevD.29.2539
	Sets the weight of events for this tree
	:param tree: Tree object with mass and lifetime info
	:param lnc_plus_lnv: if both lnc and lnv decays are possible then lifetime is reduced by a factor of 2
	:return: calculated weight.
	"""
	mass = tree.mass # GeV
	ctau = tree.ctau # mm
	mc_campaign = tree.mc_campaign
	channel = tree.channel

	if mass == -1 or ctau == -1:  # MC weighting error
		logger.debug("Can't determine the mass and lifetime of signal sample. MC mass-lifetime weight will be set to 1!!")
		return 1

	# define luminosity for the different mc campaigns
	lumi = {'mc16a': 36.20766, 'mc16d': 44.30740, 'mc16e': 58.45010, None: 1.0}
	# lumi_tot = lumi['mc16a'] + lumi['mc16d'] + lumi['mc16e']

	if tree.is_data or tree.not_hnl_mc:  # you are running on data non non-hnl MC
		return 1
	else:  # you are running on a signal MC file
		# calculate Gronau coupling; parametrization depends on coupling flavour you are probing
		if channel == 'uuu' or channel == 'uue':
			U2Gronau = 4.49e-12 * 3e8 * mass ** (-5.19) / (ctau / 1000)  # LNC prediction
		if channel == 'eee' or channel == 'eeu':
			U2Gronau = 4.15e-12 * 3e8 * mass ** (-5.17) / (ctau / 1000)  # LNC prediction

		# HNL decays in only a lepton-number conserving way
		U2_LNC_only = U2Gronau
		xsec_LNC_only = hnl_xsec(br=tree.br, U2=U2_LNC_only, mass=mass)  # in fb
		# if HNL decays to LNC and LNV, then lifetime is reduced by a factor of 2 (more decay channels available)
		U2_LNC_plus_LNV = 0.5 * U2Gronau
		xsec_LNC_plus_LNV = hnl_xsec(br=tree.br, U2=U2_LNC_plus_LNV, mass=mass)  # in fb
		# mass-lifetime weight = BR(N->llv) * L * xsec / total num. of MC events
		# LNC and LNV branches are split into into separate LNC and LNV branches
		# total num. of MC events = (tree.all_entries / 2) because pythia samples have a 50% mix of LNC+ LNV
		weight_LNC_only = lumi[mc_campaign] * xsec_LNC_only / (tree.all_entries / 2)
		# for LNC plus LNV then all MC events are added to the output tree so total number of events == tree.all_entries
		weight_LNC_plus_LNV = lumi[mc_campaign] * xsec_LNC_plus_LNV / tree.all_entries

	if lnc_plus_lnv:
		return weight_LNC_plus_LNV
	else:
		return weight_LNC_only


class Truth:
	def __init__(self):
		self.HNL_vec = ROOT.TLorentzVector()
		self.dNu_vec =  ROOT.TLorentzVector()
		self.trkVec = []
		self.dLepVec = []
		self.dLepCharge = []
		self.dEl = []
		self.dEl_charge = []
		self.dEl_d0 = []
		self.dMu = []
		self.dMu_charge = []
		self.dMu_d0 = []
		self.dTrk_d0 = []
		self.truth_dvx = -1
		self.truth_dvy = -1
		self.truth_dvz = -1
		self.truth_dv = ROOT.TLorentzVector()
		self.truth_dvr = -1
		self.truth_pvx = -1
		self.truth_pvy = -1
		self.truth_pvz = -1
		self.truth_pv = ROOT.TLorentzVector()
		self.W_vec = ROOT.TLorentzVector()
		self.W_charge =  -2
		self.plep_vec = ROOT.TLorentzVector()
		self.plep_charge = -99
		self.mhnl = -1
		self.dvmass = -1
		self.HNL_pdgID = 50
		self.gamma = 1
		self.beta = 1
		self.properLifetime = -1


	def getTruthParticles(self, tree):
		for ivx in range(len(tree['truthVtx_parent_pdgId'])):
			# get the DV!
			if abs(tree['truthVtx_parent_pdgId'][ivx]) == 50:  # PDGID 50: Heavy Neutral Lepton
				if len(tree['truthVtx_outP_pdgId'][ivx]) == 3:  # Has three children (two leptons and neutrino)

					self.truth_dvx = tree['truthVtx_x'][ivx]
					self.truth_dvy = tree['truthVtx_y'][ivx]
					self.truth_dvz = tree['truthVtx_z'][ivx]
					# mod by Christian for proper lifetime calculation
					self.gamma = tree['truthVtx_parent_E'][ivx] / tree['truthVtx_parent_M'][ivx]
					self.beta = np.sqrt(1 - 1 / self.gamma ** 2)

					self.truth_dv = ROOT.TVector3( self.truth_dvx, self.truth_dvy, self.truth_dvz )
					self.truth_dvr = np.sqrt(self.truth_dvx**2 + self.truth_dvy**2)
					visTrkVec =  ROOT.TLorentzVector()
					truthVec  =  ROOT.TLorentzVector()
					nu_vec =  ROOT.TLorentzVector()

					for i in range(len(tree['truthVtx_outP_pt'][ivx])):
						trk_pdgId = abs(tree['truthVtx_outP_pdgId'][ivx][i])
						if trk_pdgId == 13:
							TrkVec  =  ROOT.TLorentzVector()
							TrkVec.SetPtEtaPhiM(tree['truthVtx_outP_pt'][ivx][i],
												tree['truthVtx_outP_eta'][ivx][i],
												tree['truthVtx_outP_phi'][ivx][i],
												tree['truthVtx_outP_M'][ivx][i]
												)
							mu_charge = tree['truthVtx_outP_charge'][ivx][i]
							mu_d0 = tree['truthVtx_outP_d0'][ivx][i]
							self.dMu.append(TrkVec)
							self.dMu_charge.append(mu_charge)
							self.dMu_d0.append(mu_d0)
						if trk_pdgId == 11:
							TrkVec  =  ROOT.TLorentzVector()
							TrkVec.SetPtEtaPhiM(tree['truthVtx_outP_pt'][ivx][i],
												tree['truthVtx_outP_eta'][ivx][i],
												tree['truthVtx_outP_phi'][ivx][i],
												tree['truthVtx_outP_M'][ivx][i]
												)
							el_charge = tree['truthVtx_outP_charge'][ivx][i]
							el_d0 = tree['truthVtx_outP_d0'][ivx][i]
							self.dEl.append(TrkVec)
							self.dEl_charge.append(el_charge)
							self.dEl_d0.append(el_d0)

						if trk_pdgId == 13 or trk_pdgId == 11: # is track a muon of electron? Then these are our visible (charged) truth tracks
							visTrkVec =  ROOT.TLorentzVector()
							visTrkVec.SetPtEtaPhiM(tree['truthVtx_outP_pt'][ivx][i],
												tree['truthVtx_outP_eta'][ivx][i],
												tree['truthVtx_outP_phi'][ivx][i],
												tree['truthVtx_outP_M'][ivx][i]
												)
							self.trkVec.append(visTrkVec) #only add visible leptons to trkVec list
						else: # remaining child is the neutrino
							nu_vec.SetPtEtaPhiM(tree['truthVtx_outP_pt'][ivx][i],
												tree['truthVtx_outP_eta'][ivx][i],
												tree['truthVtx_outP_phi'][ivx][i],
												tree['truthVtx_outP_M'][ivx][i]
												)
							self.dNu_vec = nu_vec

					for i in range(len(tree['truthVtx_outP_pt'][ivx])):
						dLepVec  =  ROOT.TLorentzVector()
						dLepVec.SetPtEtaPhiM(tree['truthVtx_outP_pt'][ivx][i],
												tree['truthVtx_outP_eta'][ivx][i],
												tree['truthVtx_outP_phi'][ivx][i],
												tree['truthVtx_outP_M'][ivx][i]
												)
						self.dLepVec.append(dLepVec) #add all the displaced leptons to one list in the order they are in pythia
						self.dLepCharge.append(tree['truthVtx_outP_charge'][ivx][i])
						# self.dTrk_d0.append(tree['truthVtx_outP_d0'][ivx][i])
						self.dTrk_d0.append(-1) # fill with -1 for now, default DHNLalg does not have truth d0

					self.HNL_vec.SetPtEtaPhiM(tree['truthVtx_parent_pt'][ivx],
											tree['truthVtx_parent_eta'][ivx],
											tree['truthVtx_parent_phi'][ivx],
											tree['truthVtx_parent_M'][ivx]
											)


			# get the primary vertex
			if abs(tree['truthVtx_parent_pdgId'][ivx]) == 24:  # PDGID 24: W Boson
				if len(tree['truthVtx_outP_pdgId'][ivx]) == 2:  # Has two children (HNL and lepton)
					# TODO: Should we be checking if one of the children is an HNL?
					self.truth_pvx = tree['truthVtx_x'][ivx]
					self.truth_pvy = tree['truthVtx_y'][ivx]
					self.truth_pvz = tree['truthVtx_z'][ivx]
					self.truth_pv = ROOT.TVector3( self.truth_pvx, self.truth_pvy, self.truth_pvz )

					self.plep_vec.SetPtEtaPhiM(tree['truthVtx_outP_pt'][ivx][0],
											   tree['truthVtx_outP_eta'][ivx][0],
											   tree['truthVtx_outP_phi'][ivx][0],
											   tree['truthVtx_outP_M'][ivx][0]
											   )
					self.plep_charge = tree['truthVtx_outP_charge'][ivx][0]
					self.W_vec.SetPtEtaPhiM(tree['truthVtx_parent_pt'][ivx],
											tree['truthVtx_parent_eta'][ivx],
											tree['truthVtx_parent_phi'][ivx],
											tree['truthVtx_parent_M'][ivx]
											)
					self.W_charge = tree['truthVtx_parent_charge'][ivx]

			# calculate proper lifetime
			dx = np.abs(self.truth_pvx - self.truth_dvx)
			dy = np.abs(self.truth_pvy - self.truth_dvy)
			dz = np.abs(self.truth_pvz - self.truth_dvz)
			dr = np.sqrt(dx**2 + dy**2 + dz**2)
			self.properLifetime = dr/(self.gamma * self.beta)

		# TO DO: bug with truth mHNL calculation
		# try:
		# 	import selections
		# 	if len(dMu) == 2: dv_type = "mumu"
		# 	if len(dEl) == 2: dv_type = "ee"
		# 	if len(dEl) == 1 and len(dMu) == 1: dv_type = "emu"
		# 	Mhnl = selections.Mhnl(tree=tree, dv_type = dv_type, plep=self.plep_vec, dMu=dMu, dEl=dEl,use_truth=True,truth_pv=self.truth_pv,truth_dv=self.truth_dv)
		# 	self.mhnl = Mhnl.alt_mhnl
		# 	dv_vec = self.trkVec[0] + self.trkVec[1]
		# 	self.dvmass = dv_vec.M()
		# except:
		# 	pass

class Tracks:
	def __init__(self, tree):
		self.tree = tree
		# track vector with trk wrtSV quantities
		self.lepVec = []
		# track vector with standard trk quantities (usual tracking def.)
		self.std_lepVec = []
		# track vector with lepton calibrated trk quantities
		self.lepmatched_lepVec = []
		truthlepCharge = []
		self.lepIndex = []
		self.lepCharge = []
		self.lepisAssoc = []
		self.muonType = []
		self.eta = []
		self.phi = []
		self.pt = []
		self.ntracks = -1
		self.muon_isTight = []
		self.muon_isLoose = []
		self.muon_isMedium = []
		self.el_isTight = []
		self.el_isLoose = []
		self.el_isveryLoose = []
		self.el_isveryveryLoose = []
		self.el_isveryveryLooseSi = []
		self.el_isMedium = []

	def getMuons(self):
		self.ntracks = self.tree.ntrk
		# print "number of tracks: ", self.ntracks
		for itrk in range(self.ntracks):
			lepVec = ROOT.TLorentzVector()
			std_lepVec = ROOT.TLorentzVector()
			lepmatched_lepVec =  ROOT.TLorentzVector()
			if self.tree.dv('trk_muonIndex')[itrk] >= 0:  # matched muon!
				# find position of muon in the muon container that is matched to the sec vtx track
				# (works for calibrated and uncalibrated containers)
				# Default: use track quantities wrt SV
				pt = self.tree.dv('trk_pt_wrtSV')[itrk]
				eta = self.tree.dv('trk_eta_wrtSV')[itrk]
				phi = self.tree.dv('trk_phi_wrtSV')[itrk]
				# get standard track qualtities
				std_pt = self.tree.dv('trk_pt')[itrk]
				std_eta = self.tree.dv('trk_eta')[itrk]
				std_phi = self.tree.dv('trk_phi')[itrk]
				M = self.tree.dv('trk_M')[itrk]
				lepVec.SetPtEtaPhiM(pt, eta, phi, M)
				std_lepVec.SetPtEtaPhiM(std_pt, std_eta, std_phi, M)

				if not self.tree.fake_aod:
					if len(self.tree['muon_index']) > 0:
						muon_index = np.where(self.tree['muon_index'] == self.tree.dv('trk_muonIndex')[itrk])[0][0]
						self.lepIndex.append(muon_index)
						# get calibrated muon quantities (not calculated wrt DV!)
						lep_pt = self.tree['muon_pt'][muon_index]
						lep_eta = self.tree['muon_eta'][muon_index]
						lep_phi = self.tree['muon_phi'][muon_index]
						lepmatched_lepVec.SetPtEtaPhiM(lep_pt, lep_eta, lep_phi, M)
				else:
					self.lepIndex.append(-1)

				if self.tree.fake_aod:
					self.muonType.append(self.tree.dv('trk_muonType')[itrk]) # add muon type to the track class if running on fakeAODs
					self.muon_isTight.append(self.tree.dv('trk_isTight')[itrk]) # add muon quality info
					self.muon_isMedium.append(self.tree.dv('trk_isMedium')[itrk])
					self.muon_isLoose.append(self.tree.dv('trk_isLoose')[itrk])

				self.pt.append(pt)
				self.eta.append(eta)
				self.phi.append(phi)

				self.lepVec.append(lepVec)
				self.lepmatched_lepVec.append(lepmatched_lepVec)
				self.std_lepVec.append(std_lepVec)

				self.lepCharge.append(self.tree.dv('trk_charge')[itrk])
				self.lepisAssoc.append(self.tree.dv('trk_isAssociated')[itrk])
			else:
				continue

	def getElectrons(self):
		self.ntracks = self.tree.ntrk

		for itrk in range(self.ntracks):
			lepVec = ROOT.TLorentzVector()
			std_lepVec = ROOT.TLorentzVector()
			lepmatched_lepVec =  ROOT.TLorentzVector()
			if self.tree.dv('trk_electronIndex')[itrk] >= 0:  # matched electron!
				# remove electrons that are also matched to muons!
				if self.tree.dv('trk_muonIndex')[itrk] >= 0:
					if len(self.tree['muon_index']) > 0 and self.tree.fake_aod == False: # dont think we need this unless debugging overlapping muons -DT
						muon_index = np.where(self.tree['muon_index'] == self.tree.dv('trk_muonIndex')[itrk])[0][0]
						# print muon_index
						# print "track is matched to both muon and electron!"
					continue

				# Default: use track quantities wrt SV
				pt = self.tree.dv('trk_pt_wrtSV')[itrk]
				eta = self.tree.dv('trk_eta_wrtSV')[itrk]
				phi = self.tree.dv('trk_phi_wrtSV')[itrk]
				# get standard track qualtities
				std_pt = self.tree.dv('trk_pt')[itrk]
				std_eta = self.tree.dv('trk_eta')[itrk]
				std_phi = self.tree.dv('trk_phi')[itrk]
				M = self.tree.dv('trk_M')[itrk]
				lepVec.SetPtEtaPhiM(pt, eta, phi, M)
				std_lepVec.SetPtEtaPhiM(std_pt, std_eta, std_phi, M)

				# find position of electron in the electron container that is matched to the sec vtx track
				# (works for calibrated and uncalibrated containers)
				if not self.tree.fake_aod: 
					if len(self.tree['el_index']) > 0:
						el_index = np.where(self.tree['el_index'] == self.tree.dv('trk_electronIndex')[itrk])[0][0]
						# use calibrated muon quantities (not calculated wrt DV!)
						lep_pt = self.tree['el_pt'][el_index]
						lep_eta = self.tree['el_eta'][el_index]
						lep_phi = self.tree['el_phi'][el_index]
						lepmatched_lepVec.SetPtEtaPhiM(lep_pt, lep_eta, lep_phi, M)

						self.lepIndex.append(el_index)
				else:
					self.lepIndex.append(-1)

				if self.tree.fake_aod:
					self.el_isTight.append(self.tree.dv('trk_isTight')[itrk]) # add muon quality info
					self.el_isMedium.append(self.tree.dv('trk_isMedium')[itrk])
					self.el_isLoose.append(self.tree.dv('trk_isLoose')[itrk])
					self.el_isveryLoose.append(self.tree.dv('trk_isVeryLoose')[itrk])
					self.el_isveryveryLoose.append(self.tree.dv('trk_isVeryVeryLoose')[itrk])
					self.el_isveryveryLooseSi.append(self.tree.dv('trk_isVeryVeryLooseSi')[itrk])

				self.pt.append(pt)
				self.eta.append(eta)
				self.phi.append(phi)

				self.lepVec.append(lepVec)
				self.lepmatched_lepVec.append(lepmatched_lepVec)
				self.std_lepVec.append(std_lepVec)

				self.lepCharge.append(self.tree.dv('trk_charge')[itrk])
				self.lepisAssoc.append(self.tree.dv('trk_isAssociated')[itrk])
			else:
				continue

	def getTracks(self, idv=-1):
		"""Fills the Track object with a collection of track vectors.
		No lepton matching enforced."""
		at_idv = self.tree.idv if idv < 0 else idv
		at_ievt = self.tree.ievt
		prefix = self.tree.dv_prefix+'_'
		self.ntracks = self.tree.ntrk
		for itrk in range(self.tree.ntrk):
			trkvec = ROOT.TLorentzVector()
			pt = self.tree.get_at(prefix+'trk_pt_wrtSV', at_ievt, at_idv, itrk)
			eta = self.tree.get_at(prefix+'trk_eta_wrtSV', at_ievt, at_idv, itrk)
			phi = self.tree.get_at(prefix+'trk_phi_wrtSV', at_ievt, at_idv, itrk)
			M = self.tree.get_at(prefix+'trk_M', at_ievt, at_idv, itrk)

			trkvec.SetPtEtaPhiM(pt, eta, phi, M)

			self.lepVec.append(trkvec)
			self.eta.append(eta)
			self.phi.append(phi)
			self.pt.append(pt)
			self.lepCharge.append(self.tree.dv('trk_charge')[itrk])
			self.lepisAssoc.append(self.tree.dv('trk_isAssociated')[itrk])


class Muons(Tracks, object):
	def __init__(self, tree):
		super(Muons, self).__init__(tree)
		self.getMuons()


class Electrons(Tracks, object):
	def __init__(self, tree):
		super(Electrons, self).__init__(tree)
		self.getElectrons()


class FileInfo:
	def __init__(self, infile, channel=""):
		self.mass = -1  # signal mass of HNL in GeV
		self.ctau = -1  # in mm
		self.dsid = None
		# read dsids from offical samples by splitting according to '.' in sample name
		if len([int(s) for s in infile.split(".") if s.isdigit()]) != 0:
			self.dsid = [int(s) for s in infile.split(".") if s.isdigit()][0]

		self.MC_campaign = None
		self.ctau_str = ""
		self.mass_str = ""
		
		sig_info = mc_info(self.dsid)
		self.file_ch = sig_info.ch_str # used in job submission do not delete!

		if "lt1dd" in infile or "1mm" in infile or sig_info.ctau_str == "lt1dd":
			self.ctau = 1.0
			self.ctau_str = "1mm"
		elif "lt10dd" in infile or "10mm" in infile or sig_info.ctau_str == "lt10dd":
			self.ctau = 10.0
			self.ctau_str = "10mm"
		elif "lt100dd" in infile or "100mm" in infile or sig_info.ctau_str == "lt100dd":
			self.ctau = 100.0
			self.ctau_str = "100mm"

		if "_3G" in infile or sig_info.mass_str == "3G":
			self.mass = 3.0
			self.mass_str = "3G"
		elif "_4G" in infile or sig_info.mass_str == "4G":
			self.mass = 4.0
			self.mass_str = "4G"
		elif "_4p5G" in infile or sig_info.mass_str == "4p5G":
			self.mass = 4.5
			self.mass_str = "4p5G"
		elif "_5G" in infile or sig_info.mass_str == "5G":
			self.mass = 5.0
			self.mass_str = "5G"
		elif "_7p5G" in infile or sig_info.mass_str == "7p5G":
			self.mass = 7.5
			self.mass_str = "7p5G"
		elif "_10G" in infile or sig_info.mass_str == "10G":
			self.mass = 10.0
			self.mass_str = "10G"
		elif "_12p5G" in infile or sig_info.mass_str == "12p5G":
			self.mass = 12.5
			self.mass_str = "12p5G"
		elif "_15G" in infile or sig_info.mass_str == "15G":
			self.mass = 15.0
			self.mass_str = "15G"
		elif "_17p5G" in infile or sig_info.mass_str == "17p5G":
			self.mass = 17.5
			self.mass_str = "17p5G"
		elif "_20G" in infile or sig_info.mass_str == "20G":
			self.mass = 20.0
			self.mass_str = "20G"

		# two rtags for different reconstruction of our signal samples
		# r11915,r11916,r11891 are the latest ones
		# r10740,r10740,r10740 are the original ones

		if "r11915" in infile or "mc16a" in infile or "r10740" in infile:
			self.MC_campaign = "mc16a"
		if "r11916" in infile or "mc16d" in infile or "r10739" in infile:
			self.MC_campaign = "mc16d"
		if "r11891" in infile or "mc16e" in infile or "r10790" in infile:
			self.MC_campaign = "mc16e"

		# More flexibility for non-signal samples
		self.output_filename = "histograms"
		if self.MC_campaign:
			self.output_filename += "_" + self.MC_campaign
		else:
			self.output_filename += "_mc"
		if self.mass_str: self.output_filename += "_" + self.mass_str
		if self.ctau_str: self.output_filename += "_" + self.ctau_str
		self.output_filename += "_" + channel + ".root"

		f_br = ReadBRdat(os.path.dirname(os.path.abspath(__file__))+'/../data/BR/HNL_branching_20GeV.dat')
		self.br = f_br.get_BR(channel, self.mass)


class mc_info:
	def __init__(self, dsid):
		mc_info = {}
		mc_info[311602] = ["uuu", "3G", "lt1dd"]
		mc_info[311603] = ["uuu", "3G", "lt10dd"]
		mc_info[311604] = ["uuu", "3G", "lt100dd"]
		mc_info[311605] = ["uue", "3G", "lt1dd"]
		mc_info[311606] = ["uue", "3G", "lt10dd"]
		mc_info[311607] = ["uue", "3G", "lt100dd"]
		mc_info[311608] = ["uuu", "4G", "lt1dd"]
		mc_info[311609] = ["uuu", "4G", "lt10dd"]
		mc_info[311610] = ["uuu", "4G", "lt100dd"]
		mc_info[311611] = ["uue", "4G", "lt1dd"]
		mc_info[311612] = ["uue", "4G", "lt10dd"]
		mc_info[311613] = ["uue", "4G", "lt100dd"]
		mc_info[311614] = ["uuu", "4p5G", "lt1dd"]
		mc_info[311615] = ["uuu", "4p5G", "lt10dd"]
		mc_info[311616] = ["uuu", "4p5G", "lt100dd"]
		mc_info[311617] = ["uue", "4p5G", "lt1dd"]
		mc_info[311618] = ["uue", "4p5G", "lt10dd"]
		mc_info[311619] = ["uue", "4p5G", "lt100dd"]
		mc_info[311620] = ["uuu", "5G", "lt1dd"]
		mc_info[311621] = ["uuu", "5G", "lt10dd"]
		mc_info[311622] = ["uuu", "5G", "lt100dd"]
		mc_info[311623] = ["uue", "5G", "lt1dd"]
		mc_info[311624] = ["uue", "5G", "lt10dd"]
		mc_info[311625] = ["uue", "5G", "lt100dd"]
		mc_info[311626] = ["uuu", "7p5G", "lt1dd"]
		mc_info[311627] = ["uuu", "7p5G", "lt10dd"]
		mc_info[311628] = ["uuu", "7p5G", "lt100dd"]
		mc_info[311629] = ["uue", "7p5G", "lt1dd"]
		mc_info[311630] = ["uue", "7p5G", "lt10dd"]
		mc_info[311631] = ["uue", "7p5G", "lt100dd"]
		mc_info[311632] = ["uuu", "10G", "lt1dd"]
		mc_info[311633] = ["uuu", "10G", "lt10dd"]
		mc_info[311634] = ["uuu", "10G", "lt100dd"]
		mc_info[311635] = ["uue", "10G", "lt1dd"]
		mc_info[311636] = ["uue", "10G", "lt10dd"]
		mc_info[311637] = ["uue", "10G", "lt100dd"]
		mc_info[311638] = ["uuu", "12p5G", "lt1dd"]
		mc_info[311639] = ["uuu", "12p5G", "lt10dd"]
		mc_info[311640] = ["uuu", "12p5G", "lt100dd"]
		mc_info[311641] = ["uue", "12p5G", "lt1dd"]
		mc_info[311642] = ["uue", "12p5G", "lt10dd"]
		mc_info[311643] = ["uue", "12p5G", "lt100dd"]
		mc_info[311644] = ["uuu", "15G", "lt1dd"]
		mc_info[311645] = ["uuu", "15G", "lt10dd"]
		mc_info[311646] = ["uuu", "15G", "lt100dd"]
		mc_info[311647] = ["uue", "15G", "lt1dd"]
		mc_info[311648] = ["uue", "15G", "lt10dd"]
		mc_info[311649] = ["uue", "15G", "lt100dd"]
		mc_info[311650] = ["uuu", "17p5G", "lt1dd"]
		mc_info[311651] = ["uuu", "17p5G", "lt10dd"]
		mc_info[311652] = ["uuu", "17p5G", "lt100dd"]
		mc_info[311653] = ["uue", "17p5G", "lt1dd"]
		mc_info[311654] = ["uue", "17p5G", "lt10dd"]
		mc_info[311655] = ["uue", "17p5G", "lt100dd"]
		mc_info[311656] = ["uuu", "20G", "lt1dd"]
		mc_info[311657] = ["uuu", "20G", "lt10dd"]
		mc_info[311658] = ["uuu", "20G", "lt100dd"]
		mc_info[311659] = ["uue", "20G", "lt1dd"]
		mc_info[311660] = ["uue", "20G", "lt10dd"]
		mc_info[311661] = ["uue", "20G", "lt100dd"]
		mc_info[312956] = ["eee", "3G", "lt1dd"]
		mc_info[312957] = ["eee", "3G", "lt10dd"]
		mc_info[312958] = ["eee", "3G", "lt100dd"]
		mc_info[312959] = ["eeu", "3G", "lt1dd"]
		mc_info[312960] = ["eeu", "3G", "lt10dd"]
		mc_info[312961] = ["eeu", "3G", "lt100dd"]
		mc_info[312962] = ["eee", "4G", "lt1dd"]
		mc_info[312963] = ["eee", "4G", "lt10dd"]
		mc_info[312964] = ["eee", "4G", "lt100dd"]
		mc_info[312965] = ["eeu", "4G", "lt1dd"]
		mc_info[312966] = ["eeu", "4G", "lt10dd"]
		mc_info[312967] = ["eeu", "4G", "lt100dd"]
		mc_info[312968] = ["eee", "4p5G", "lt1dd"]
		mc_info[312969] = ["eee", "4p5G", "lt10dd"]
		mc_info[312970] = ["eee", "4p5G", "lt100dd"]
		mc_info[312971] = ["eeu", "4p5G", "lt1dd"]
		mc_info[312972] = ["eeu", "4p5G", "lt10dd"]
		mc_info[312973] = ["eeu", "4p5G", "lt100dd"]
		mc_info[312974] = ["eee", "5G", "lt1dd"]
		mc_info[312975] = ["eee", "5G", "lt10dd"]
		mc_info[312976] = ["eee", "5G", "lt100dd"]
		mc_info[312977] = ["eeu", "5G", "lt1dd"]
		mc_info[312978] = ["eeu", "5G", "lt10dd"]
		mc_info[312979] = ["eeu", "5G", "lt100dd"]
		mc_info[312980] = ["eee", "7p5G", "lt1dd"]
		mc_info[312981] = ["eee", "7p5G", "lt10dd"]
		mc_info[312982] = ["eee", "7p5G", "lt100dd"]
		mc_info[312983] = ["eeu", "7p5G", "lt1dd"]
		mc_info[312984] = ["eeu", "7p5G", "lt10dd"]
		mc_info[312985] = ["eeu", "7p5G", "lt100dd"]
		mc_info[312986] = ["eee", "10G", "lt1dd"]
		mc_info[312987] = ["eee", "10G", "lt10dd"]
		mc_info[312988] = ["eee", "10G", "lt100dd"]
		mc_info[312989] = ["eeu", "10G", "lt1dd"]
		mc_info[312990] = ["eeu", "10G", "lt10dd"]
		mc_info[312991] = ["eeu", "10G", "lt100dd"]
		mc_info[312992] = ["eee", "12p5G", "lt1dd"]
		mc_info[312993] = ["eee", "12p5G", "lt10dd"]
		mc_info[312994] = ["eee", "12p5G", "lt100dd"]
		mc_info[312995] = ["eeu", "12p5G", "lt1dd"]
		mc_info[312996] = ["eeu", "12p5G", "lt10dd"]
		mc_info[312997] = ["eeu", "12p5G", "lt100dd"]
		mc_info[312998] = ["eee", "15G", "lt1dd"]
		mc_info[312999] = ["eee", "15G", "lt10dd"]
		mc_info[313000] = ["eee", "15G", "lt100dd"]
		mc_info[313001] = ["eeu", "15G", "lt1dd"]
		mc_info[313002] = ["eeu", "15G", "lt10dd"]
		mc_info[313003] = ["eeu", "15G", "lt100dd"]
		mc_info[313004] = ["eee", "17p5G", "lt1dd"]
		mc_info[313005] = ["eee", "17p5G", "lt10dd"]
		mc_info[313006] = ["eee", "17p5G", "lt100dd"]
		mc_info[313007] = ["eeu", "17p5G", "lt1dd"]
		mc_info[313008] = ["eeu", "17p5G", "lt10dd"]
		mc_info[313009] = ["eeu", "17p5G", "lt100dd"]
		mc_info[313010] = ["eee", "20G", "lt1dd"]
		mc_info[313011] = ["eee", "20G", "lt10dd"]
		mc_info[313012] = ["eee", "20G", "lt100dd"]
		mc_info[313013] = ["eeu", "20G", "lt1dd"]
		mc_info[313014] = ["eeu", "20G", "lt10dd"]
		mc_info[313015] = ["eeu", "20G", "lt100dd"]

		if dsid is None:
			logger.warning("No dsid")
			self.mass_str = None
			self.ctau_str = None
			self.ch_str = None
		else:
			pmuon_dsid = (311602 <= dsid) and (dsid <= 311661)
			pel_dsid = (312956 <= dsid) and (dsid <= 313015)

			if pmuon_dsid or pel_dsid:
				self.mass_str = mc_info[dsid][1]
				self.ctau_str = mc_info[dsid][2]
				self.ch_str = mc_info[dsid][0]
				logger.info("This sample is type: {}, mass: {}, lifetime: {}".format(self.ch_str, self.mass_str, self.ctau_str, ))
			else:
				logger.warning("dsid {} is not registered. If running on HNL signal, please check your signal sample".format(dsid))
				self.mass_str = None
				self.ctau_str = None
				self.ch_str = None


# Define trigger lists here
# trigger lists taken from https://acode-browser1.usatlas.bnl.gov/lxr/source/athena/PhysicsAnalysis/SUSYPhys/LongLivedParticleDPDMaker/share/PhysDESDM_HNL.py?v=21.0#0008
# seperated by year using comments from the above link and cross checking with this twiki: https://twiki.cern.ch/twiki/bin/view/Atlas/LowestUnprescaled

# Single muon triggers used in DHNL analysis
SingleMuonTriggerlist = [
	"HLT_mu26_ivarmedium", # 2017-2018
	# "HLT_mu60_0eta105_msonly", # all
	"HLT_mu20_iloose_L1MU15", # 2015
	# "HLT_mu24_iloose", 
	# "HLT_mu24_ivarloose", 
	# "HLT_mu24_imedium",
	# "HLT_mu24_ivarmedium",
	# "HLT_mu26_imedium", 
	]

SingleMuonTriggerlist_2018 = [
	"HLT_mu26_ivarmedium", 
	# "HLT_mu60_0eta105_msonly",
	]

SingleMuonTriggerlist_2017 = [
	"HLT_mu26_ivarmedium", 
	# "HLT_mu60_0eta105_msonly",
	]

SingleMuonTriggerlist_2015_2016 = [
	"HLT_mu26_ivarmedium", 
	# "HLT_mu60_0eta105_msonly",
	"HLT_mu20_iloose_L1MU15", 
	# "HLT_mu24_iloose", 
	# "HLT_mu24_ivarloose", 
	# "HLT_mu24_ivarmedium",
	# "HLT_mu24_imedium", 
	# "HLT_mu26_imedium", 
	]

# Single electron triggers used in DHNL analysis
SingleElectronTriggerlist = [
	"HLT_e24_lhmedium_L1EM20VH", # 2015
	"HLT_e26_lhtight_nod0_ivarloose", # 2016-2018
	# "HLT_e60_lhmedium_nod0", # 2016-2018 # Kick this out?
	# "HLT_e140_lhloose_nod0", # 2016-2018 # Kick this out?
	# "HLT_e24_lhtight_nod0_ivarloose",  # not in recommendations
	# "HLT_e26_lhtight_nod0", # not in recommendations
	]

SingleElectronTriggerlist_2018 = [
	"HLT_e26_lhtight_nod0_ivarloose",
	# "HLT_e60_lhmedium_nod0",
	# "HLT_e140_lhloose_nod0",
	# "HLT_e26_lhtight_nod0", # not in recommendations
	]

SingleElectronTriggerlist_2017 = [
	"HLT_e26_lhtight_nod0_ivarloose",
	# "HLT_e60_lhmedium_nod0",
	# "HLT_e140_lhloose_nod0",
	 ]

SingleElectronTriggerlist_2015_2016 = [
	"HLT_e24_lhmedium_L1EM20VH",
	"HLT_e26_lhtight_nod0_ivarloose",
	# "HLT_e60_lhmedium_nod0",
	# "HLT_e140_lhloose_nod0",
	# "HLT_e24_lhtight_nod0_ivarloose", # not in recommendations
	]




# full list of triggers from DHNL filter. Some high-threshold triggers removed from the list as they did not provide addtional signal efficiency
# muons triggers
apiSingleMuonTriggerlist = ["HLT_mu20_iloose_L1MU15", "HLT_mu24_iloose", "HLT_mu24_ivarloose", "HLT_mu24_imedium","HLT_mu24_ivarmedium",
							"HLT_mu26_imedium", "HLT_mu26_ivarmedium", "HLT_mu40", "HLT_mu50",
							"HLT_mu60_0eta105_msonly"]
apiSingleMuonTriggerlist_2018 = ["HLT_mu26_ivarmedium", "HLT_mu50", "HLT_mu60_0eta105_msonly" ]
apiSingleMuonTriggerlist_2017 = ["HLT_mu26_ivarmedium", "HLT_mu50", "HLT_mu60_0eta105_msonly" ]
apiSingleMuonTriggerlist_2015_2016 = ["HLT_mu20_iloose_L1MU15", "HLT_mu24_iloose", "HLT_mu24_ivarloose", "HLT_mu24_ivarmedium",
									       "HLT_mu24_imedium", "HLT_mu26_imedium", "HLT_mu26_ivarmedium", "HLT_mu40", "HLT_mu50", "HLT_mu60_0eta105_msonly"]
#electron triggers
apiSingleElectronTriggerlist = ["HLT_e24_lhmedium_L1EM20VH", "HLT_e24_lhtight_nod0_ivarloose", "HLT_e26_lhtight_nod0",
                                "HLT_e26_lhtight_nod0_ivarloose", "HLT_e60_lhmedium_nod0", "HLT_e60_lhmedium",
                                "HLT_e60_medium", "HLT_e120_lhloose", "HLT_e140_lhloose_nod0", "HLT_e300_etcut"]

apiSingleElectronTriggerlist_2018 = ["HLT_e26_lhtight_nod0_ivarloose", "HLT_e26_lhtight_nod0", "HLT_e60_lhmedium_nod0",
									 "HLT_e140_lhloose_nod0", "HLT_e300_etcut"]

apiSingleElectronTriggerlist_2017 = ["HLT_e26_lhtight_nod0_ivarloose", "HLT_e60_lhmedium_nod0", "HLT_e140_lhloose_nod0", "HLT_e300_etcut"]

apiSingleElectronTriggerlist_2015_2016 = ["HLT_e24_lhmedium_L1EM20VH", "HLT_e24_lhtight_nod0_ivarloose", "HLT_e26_lhtight_nod0_ivarloose", "HLT_e60_lhmedium_nod0",
										  "HLT_e60_lhmedium", "HLT_e60_medium", "HLT_e120_lhloose", "HLT_e140_lhloose_nod0", "HLT_e300_etcut"]


apiMultiMuonTriggerlist = ["HLT_2mu10", "HLT_2mu10_nomucomb", "HLT_2mu14", "HLT_2mu14_nomucomb",
                           "HLT_mu20_nomucomb_mu6noL1_nscan03", "HLT_mu20_mu8noL1", "HLT_mu22_mu8noL1",
                           "HLT_mu22_mu8noL1_calotag_0eta010", "HLT_3mu4", "HLT_mu6_2mu4", "HLT_3mu6",
                           "HLT_3mu6_msonly", "HLT_mu20_2mu4noL1", "HLT_4mu4",
                           "HLT_mu11_nomucomb_2mu4noL1_nscan03_L1MU11_2MU6",
                           "HLT_mu11_nomucomb_2mu4noL1_nscan03_L1MU11_2MU6_bTau",
                           "HLT_mu20_msonly_mu10noL1_msonly_nscan05_noComb",
                           "HLT_mu11_nomucomb_mu6noL1_nscan03_L1MU11_2MU6_bTau",
                           "HLT_mu6_nomucomb_2mu4_nomucomb_bTau_L1MU6_3MU4",
                           "HLT_2mu6_nomucomb_mu4_nomucomb_bTau_L12MU6_3MU4"]
apiMultiElectronTriggerlist = ["HLT_2e12_lhloose_L12EM10VH", "HLT_2e15_lhvloose_nod0_L12EM13VH",
                               "HLT_e17_lhloose_2e9_lhloose", "HLT_e17_lhloose_nod0_2e9_lhloose_nod0",
                               "HLT_e17_lhloose_nod0_2e10_lhloose_nod0_L1EM15VH_3EM8VH", "HLT_2e17_lhvloose_nod0",
                               "HLT_2e17_lhvloose_nod0_L12EM15VHI", "HLT_2e24_lhvloose_nod0",
                               "HLT_e24_lhvloose_nod0_2e12_lhvloose_nod0_L1EM20VH_3EM10VH"]
apiElectronMuonTriggerlist = ["HLT_e17_lhloose_mu14", "HLT_e17_lhloose_nod0_mu14",
                              "HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1", "HLT_e26_lhmedium_nod0_mu8noL1",
                              "HLT_e7_lhmedium_nod0_mu24", "HLT_e12_lhloose_nod0_2mu10", "HLT_2e12_lhloose_mu10",
                              "HLT_2e12_lhloose_nod0_mu10", "HLT_e7_lhmedium_mu24", "HLT_e12_lhloose_2mu10"]


DAOD_RPVLLTriggerlist = ['HLT_mu26_ivarmedium', 'HLT_xe110_pufit_xe65_L1XE50', 'HLT_e26_lhtight_nod0_ivarloose',
                         'HLT_mu50', 'HLT_mu60_0eta105_msonly', 'HLT_mu80_msonly_3layersEC', 'HLT_mu22_mu8noL1',
                         'HLT_xe120_pufit_L1XE50', 'HLT_xe110_pufit_xe70_L1XE50', 'HLT_j0_perf_ds1_L1J100',
                         'HLT_g25_medium_L1EM22VHI_4j35_0eta490_invm1000', 'HLT_5j70_0eta240_L14J15',
                         'HLT_5j50_gsc70_boffperf_split_0eta240_L14J15', 'HLT_6j45_gsc55_boffperf_split_0eta240_L14J15',
                         'HLT_7j45_L14J15', 'HLT_7j35_gsc45_boffperf_split_L14J15', 'HLT_2mu14',
                         'HLT_mu6_dRl1_mu20_msonly_iloosems_mu6noL1_dRl1_msonly', 'HLT_3j200', 'HLT_4j120',
                         'HLT_4j85_gsc115_boffperf_split', 'HLT_ht1000_L1J100', 'HLT_2g20_tight_icalovloose_L12EM15VHI',
                         'HLT_2g22_tight_L12EM15VHI', 'HLT_2g25_tight_L12EM20VH',
                         'HLT_e26_lhtight_nod0_e15_etcut_L1EM7_Zee', 'HLT_2e17_lhvloose_nod0_L12EM15VHI',
                         'HLT_2e24_lhvloose_nod0', 'HLT_e24_lhmedium_nod0_L1EM20VH_g25_medium',
                         'HLT_g35_medium_g25_medium_L12EM20VH', 'HLT_j30_jes_cleanLLP_PS_llp_L1LLP-NOMATCH',
                         'HLT_j30_jes_cleanLLP_PS_llp_noiso_L1LLP-NOMATCH', 'HLT_6j50_gsc70_boffperf_split_L14J15',
                         'HLT_6j55_0eta240_L14J15', 'HLT_2j35_bmv2c1060_split_2j35_L14J15.0ETA25',
                         'HLT_2j35_bmv2c1070_split_2j35_bmv2c1085_split_L14J15.0ETA25', 'HLT_g140_loose', 'HLT_j420',
                         'HLT_j360_a10t_lcw_jes_60smcINF_j360_a10t_lcw_jes_L1SC111',
                         'HLT_j370_a10t_lcw_jes_35smcINF_j370_a10t_lcw_jes_L1SC111', 'HLT_e140_lhloose_nod0',
                         'HLT_tau160_medium1_tracktwoEF_L1TAU100', 'HLT_j30_jes_cleanLLP_PS_llp_L1LLP-RO',
                         'HLT_j30_jes_cleanLLP_PS_llp_noiso_L1LLP-RO', 'HLT_2g50_loose_L12EM20VH',
                         'HLT_j460_a10_lcw_subjes_L1J100', 'HLT_j460_a10_lcw_subjes_L1SC111', 'HLT_j460_a10r_L1J100',
                         'HLT_j460_a10r_L1SC111', 'HLT_j420_a10t_lcw_jes_35smcINF_L1J100', 'HLT_5j85_L14J15',
                         'HLT_5j60_gsc85_boffperf_split_L14J15', 'HLT_6j70_L14J15',
                         'HLT_mu4_j90_xe90_pufit_2dphi10_L1MU4_J50_XE50_DPHI-J20s2XE30',
                         'HLT_mu4_j90_xe90_pufit_2dphi10_L1MU4_XE60',
                         'HLT_e5_lhloose_nod0_j50_xe70_pufit_2dphi10_L1J40_XE50_DPHI-J20s2XE50',
                         'HLT_g80_loose_xe80noL1', 'HLT_mu20_msonly_iloosems_mu6noL1_msonly_nscan05_L1MU20_XE30',
                         'HLT_mu20_msonly_iloosems_mu6noL1_msonly_nscan05_L1MU20_J40', 'HLT_e60_lhmedium_nod0',
                         'HLT_tau200_medium1_tracktwoEF_L1TAU100', 'HLT_mu20_mu6noL1_bTau',
                         'HLT_e17_lhmedium_nod0_ivarloose_tau25_medium1_tracktwoEF', 'HLT_mu20_2mu4noL1',
                         'HLT_mu11_2mu4noL1_bNocut_L1MU11_2MU6', 'HLT_g85_tight_L1EM22VHI_3j50noL1',
                         'HLT_e5_lhvloose_nod0_mu4_j30_xe40_pufit_2dphi10_L1MU4_J30_XE40_DPHI-J20s2XE30',
                         'HLT_g300_etcut', 'HLT_e300_etcut', 'HLT_j225_gsc420_boffperf_split',
                         'HLT_j460_a10t_lcw_jes_L1J100', 'HLT_j460_a10t_lcw_jes_L1SC111',
                         'HLT_2j330_a10t_lcw_jes_35smcINF_L1J100', 'HLT_2j330_a10t_lcw_jes_35smcINF_L1SC111',
                         'HLT_g45_loose_6j45_0eta240', 'HLT_tau35_medium1_tracktwoEF_xe70_L1XE45',
                         'HLT_j55_gsc75_bmv2c1040_split_3j55_gsc75_boffperf_split',
                         'HLT_j60_gsc85_bmv2c1050_split_3j60_gsc85_boffperf_split',
                         'HLT_j45_gsc55_bmv2c1050_split_ht700_L1HT190-J15s5.ETA21',
                         'HLT_mu14_ivarloose_tau25_medium1_tracktwoEF',
                         'HLT_mu14_ivarloose_tau25_medium1_tracktwoEF_L1DR-MU10TAU12I_TAU12I-J25',
                         'HLT_mu14_tau25_medium1_tracktwoEF_xe50', 'HLT_mu14_ivarloose_tau25_medium1_tracktwoEF_xe50',
                         'HLT_j30_muvtx', 'HLT_j30_muvtx_noiso', 'HLT_g35_loose_L1EM24VHI_mu18',
                         'HLT_2j35_gsc45_bmv2c1050_split_2j35_gsc45_boffperf_split_L14J15.0ETA25',
                         'HLT_2j45_gsc55_bmv2c1060_split_2j45_gsc55_boffperf_split_L14J15.0ETA25',
                         'HLT_2j35_bmv2c1060_split_3j35', 'HLT_2j35_gsc45_bmv2c1060_split_3j35_gsc45_boffperf_split',
                         'HLT_3j50_gsc65_bmv2c1077_split_L13J35.0ETA23', 'HLT_3j35_bmv2c1070_split_j35_L14J15.0ETA25',
                         'HLT_4j35_bmv2c1077_split_L14J15.0ETA25',
                         'HLT_3j35_bmv2c1077_split_xe70_pufit_xe50_L13J15.0ETA25_XE40', 'HLT_e17_lhloose_nod0_mu14',
                         'HLT_e26_lhmedium_nod0_mu8noL1', 'HLT_e7_lhmedium_nod0_mu24', 'HLT_g25_medium_mu24',
                         'HLT_2e5_lhvloose_nod0_j40_xe70_pufit_2dphi10_L1J40_XE50_DPHI-J20s2XE50',
                         'HLT_g35_tight_icalotight_L1EM24VHI_mu18noL1',
                         'HLT_e17_lhmedium_nod0_ivarloose_tau25_medium1_tracktwoEF_xe50',
                         'HLT_e17_lhmedium_nod0_tau25_medium1_tracktwoEF_xe50',
                         'HLT_e17_lhmedium_nod0_ivarloose_tau25_medium1_tracktwoEF_L1DR-EM15TAU12I-J25',
                         'HLT_e24_lhmedium_nod0_ivarloose_tau35_medium1_tracktwoEF',
                         'HLT_e17_lhmedium_nod0_ivarloose_tau25_medium1_tracktwoEF_L1EM15VHI_2TAU12IM_4J12',
                         'HLT_g25_medium_L1EM22VHI_j35_0eta490_bmv2c1077_split_3j35_0eta490_invm700',
                         'HLT_g25_medium_L1EM22VHI_2j35_0eta490_bmv2c1077_split_2j35_0eta490',
                         'HLT_g20_tight_icaloloose_j35_bmv2c1077_split_3j35_0eta490_invm500',
                         'HLT_j85_gsc100_bmv2c1050_split_xe85_pufit_xe50_L1XE55', 'HLT_g0_hiptrt_L1EM22VHI',
                         'HLT_j30_jes_cleanLLP_PS_llp_L1TAU100', 'HLT_j30_jes_cleanLLP_PS_llp_noiso_L1TAU100',
                         'HLT_e5_lhloose_nod0_j40_xe70_pufit_2dphi10_L1XE60',
                         'HLT_j110_gsc150_boffperf_split_2j45_gsc55_bmv2c1070_split_L1J85_3J30',
                         'HLT_tau40_medium1_tracktwoEF_tau35_medium1_tracktwoEF',
                         'HLT_tau35_medium1_tracktwoEF_tau25_medium1_tracktwoEF_L1DR-TAU20ITAU12I-J25',
                         'HLT_tau35_medium1_tracktwoEF_tau25_medium1_tracktwoEF_03dR30_L1DR-TAU20ITAU12I-J25',
                         'HLT_tau80_medium1_tracktwoEF_L1TAU60_tau35_medium1_tracktwoEF_L1TAU12IM_L1TAU60_DR-TAU20ITAU12I',
                         'HLT_tau80_medium1_tracktwoEF_L1TAU60_tau60_medium1_tracktwoEF_L1TAU40',
                         'HLT_2j35_gsc45_bmv2c1070_split_xe85_pufit_xe50_L12J15_XE55',
                         'HLT_j175_gsc225_bmv2c1040_split', 'HLT_j225_gsc275_bmv2c1060_split',
                         'HLT_j225_gsc300_bmv2c1070_split', 'HLT_j225_gsc360_bmv2c1077_split',
                         'HLT_tau60_medium1_tracktwoEF_tau25_medium1_tracktwoEF_xe50',
                         'HLT_2mu4_invm1_j20_xe40_pufit_2dphi10_L12MU4_J20_XE30_DPHI-J20s2XE30',
                         'HLT_2mu4_invm1_j20_xe60_pufit_2dphi10_L12MU4_J40_XE50',
                         'HLT_mu14_ivarloose_tau35_medium1_tracktwoEF_L1MU10_TAU20IM_J25_2J20',
                         'HLT_mu14_ivarloose_tau35_medium1_tracktwoEF', 'HLT_mu11_mu6_bDimu2700_L1LFV-MU11',
                         'HLT_mu11_mu6_bTau_L1LFV-MU11', 'HLT_mu11_mu6_bBmumuxv2_L1LFV-MU11', 'HLT_mu11_mu6_bDimu2700',
                         'HLT_mu11_mu6_bBmumuxv2', 'HLT_mu11_mu6_bTau', 'HLT_e12_lhloose_nod0_2mu10',
                         'HLT_g15_loose_2mu10_msonly', 'HLT_e5_lhmedium_nod0_mu4_j30_xe65_pufit_2dphi10_L1MU4_XE60',
                         'HLT_g35_tight_icalotight_L1EM24VHI_mu15noL1_mu2noL1', 'HLT_g35_loose_L1EM24VHI_mu15_mu2noL1',
                         'HLT_3mu6_msonly', 'HLT_e70_lhloose_nod0_xe70noL1',
                         'HLT_2e5_lhvloose_nod0_j40_xe70_pufit_2dphi10_L1XE60',
                         'HLT_2j45_gsc55_bmv2c1050_split_ht300_L1HT190-J15s5.ETA21', 'HLT_mu11_mu6_bDimu_L1LFV-MU11',
                         'HLT_mu11_mu6_bDimu2700_Lxy0_L1LFV-MU11', 'HLT_mu11_mu6_bDimu_Lxy0_L1LFV-MU11',
                         'HLT_mu11_mu6_bDimu_Lxy0', 'HLT_mu11_mu6_bDimu2700_Lxy0', 'HLT_mu11_mu6_bDimu',
                         'HLT_ht300_2j40_0eta490_invm700_L1HT150-J20s5.ETA31_MJJ-400-CF_AND_2j35_gsc45_bmv2c1070_split',
                         'HLT_j55_gsc80_bmv2c1070_split_j45_gsc60_bmv2c1085_split_j45_320eta490',
                         'HLT_j80_0eta240_j60_j45_320eta490_AND_2j35_gsc45_bmv2c1070_split',
                         'HLT_j150_gsc175_bmv2c1060_split_j45_gsc60_bmv2c1060_split', 'HLT_mu6_2mu4_bDimu2700',
                         'HLT_mu10_mgonly_L1LATE-MU10_J50', 'HLT_mu10_mgonly_L1LATE-MU10_XE40',
                         'HLT_2g25_loose_g15_loose', 'HLT_mu11_mu6_bJpsimumu_L1LFV-MU11',
                         'HLT_mu11_mu6_bBmumux_BpmumuKp_L1LFV-MU11', 'HLT_mu11_mu6_bJpsimumu_Lxy0_L1LFV-MU11',
                         'HLT_mu11_mu6_bBmumux_BpmumuKp', 'HLT_mu11_mu6_bJpsimumu_Lxy0', 'HLT_mu11_mu6_bJpsimumu',
                         'HLT_mu11_mu6_bUpsimumu_L1LFV-MU11', 'HLT_mu11_mu6_bUpsimumu', 'HLT_3mu4_bTau',
                         'HLT_mu11_mu6_bBmumu_L1LFV-MU11', 'HLT_3mu6_bDimu',
                         'HLT_2mu6_bBmumu_Lxy0_L1BPH-2M9-2MU6_BPH-2DR15-2MU6', 'HLT_mu11_mu6_bBmumu', 'HLT_3mu6']

def decode_list(in_list, encoding='utf8'):
	"""Helper fucntion to make bytes python2/3 compatible"""
	out_list = []
	for item in in_list:
		if isinstance(item, bytes):
			out_list.append(item.decode(encoding))
		else:
			out_list.append(item)
	return out_list


def mom_perp(pvec, decay_vector):
	mom_perp_vec = pvec.Cross(decay_vector)
	if mom_perp_vec.Theta() < np.pi / 2.0:
		return mom_perp_vec.Mag() / decay_vector.Mag()
	elif mom_perp_vec.Theta() > np.pi / 2.0:
		return -1 * mom_perp_vec.Mag() / decay_vector.Mag()
	else:
		return -1


def mom_parall(pvec, decay_vector):
	return pvec.Dot(decay_vector) / decay_vector.Mag()


def mom_frac_parall(pvec, decay_vector):
	if pvec.Mag() == 0.0:  # protect against div by 0...
		return -1
	else:
		return mom_parall(pvec, decay_vector) / pvec.Mag()

def get_time():
	# https://bugs.python.org/issue36895#msg342267
	if sys.version_info >= (3, 3):
		return time.perf_counter()
	else: # python 2
		return time.clock()
