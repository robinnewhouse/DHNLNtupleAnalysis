import ROOT
import numpy as np
import sys
import os
import time

ROOT.PyConfig.IgnoreCommandLineOptions = True

import logging

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
				# update the br with something so the code runs -DT
				if channel == 'uee':
					return self.BRuuu[i]
				if channel == 'euu':
					return self.BReee[i]


def hnl_xsec_single_flavour_mixing(channel, br, mass, ctau, LNC_only = True):

	# calculate Gronau coupling; parametrization depends on coupling flavour you are probing
	if channel == 'uuu' or channel == 'uue' or channel == 'uee':
		U2Gronau = 4.49e-12 * 3e8 * mass ** (-5.19) / (ctau / 1000)  # LNC prediction
	if channel == 'eee' or channel == 'eeu' or channel == 'euu':
		U2Gronau = 4.15e-12 * 3e8 * mass ** (-5.17) / (ctau / 1000)  # LNC prediction

	# HNL decays in only a lepton-number conserving way
	if LNC_only: k = 1
	# if HNL decays to LNC and LNV, then lifetime is reduced by a factor of 2 (more decay channels available)
	else: k = 0.5

	mW = 80.379  # mass of W boson in GeV
	xsec = br * 20.6e6 * k * U2Gronau * ((1 - (mass / mW) ** 2) ** 2) * (1 + (mass ** 2) / (2 * mW ** 2))  # in fb
	return xsec

def hnl_xsec_generic_model(channel, br_single_flavour_mixing, mass, ctau, LNC_only = True, x_e = 1, x_mu = 0, x_tau = 0):
	# HNL decays in only a lepton-number conserving way
	if LNC_only: k = 1
	# if HNL decays to LNC and LNV, then lifetime is reduced by a factor of 2 (more decay channels available)
	else: k = 0.5

	#Gronau parametrization relates mass, coupling and lifetime

	#constants for Gronau approx
	tau_0_e = 4.15e-12 # in seconds
	tau_0_mu = 4.49e-12  # in seconds
	tau_0_tau = 1.08e-11  # in seconds	
	b_e = 5.17
	b_mu = 5.19
	b_tau = 5.44

	e_term   = x_e   * mass ** b_e / (tau_0_e * 3e8) # in 1/m 
	mu_term  = x_mu  * mass ** b_mu / (tau_0_mu * 3e8) # in 1/m 
	tau_term = x_tau * mass ** b_tau / (tau_0_tau * 3e8) # in 1/m 

	U2Gronau = (1 / (ctau / 1000))*(1 / (mu_term + e_term + tau_term )) # unitless (ctau is in mm )
	
	# compute branching ratios with single flavour mixing. 
	# Get partial widths by multiplying br total width from single flavour mixing
	if channel == "uuu" or channel == "uue" or channel == "euu" or channel == "eue":
		tau_0_single_flavour = tau_0_mu	
		b_single_flavour = b_mu
	elif channel == "eee" or channel == "eeu" or channel == "uee" or channel == "ueu": 
		tau_0_single_flavour = tau_0_e
		b_single_flavour = b_e
	else:
		logger.error("Can't determine the what constants to use to compute the partial width. Please check your sample!")
		sys.exit(1) # abort becuase of an error
	# total width with unit coupling
	total_width_single_flavour_mixing =  mass**b_single_flavour / (tau_0_single_flavour *3e8) # in 1/m
	
	# partial width should be independent of the interpretation
	# Used single flavour mixing br numbers to compute the partial widths
	partial_width = br_single_flavour_mixing * total_width_single_flavour_mixing

	# define the prod and decay ratios depending on the channel
	if channel == 'uuu' or channel == 'uue': 
		x_prod  = x_mu
		x_decay = x_mu
	if channel == 'eee' or channel == 'eeu': 
		x_prod  = x_e
		x_decay = x_e
	if channel == 'ueu' or channel == 'uee': 
		x_prod  = x_mu
		x_decay = x_e
	if channel == 'eue' or channel == 'euu':
		x_prod  = x_e
		x_decay = x_mu

	mW = 80.379  # mass of W boson in GeV
	xsec = 20.6e6 * ((1 - (mass / mW) ** 2) ** 2) * (1 + (mass ** 2) / (2 * mW ** 2)) * partial_width * k *  U2Gronau ** 2 * (ctau/1000) * x_prod * x_decay  # in fb
	return xsec


def get_mass_lt_weight(tree, lnc_plus_lnv=False, single_flavour_mixing = True, ih_mixing = False, nh_mixing = False):
	"""
	Calculates the weight of the event based on the Gronau parametrization
	https://journals.aps.org/prd/abstract/10.1103/PhysRevD.29.2539
	Sets the weight of events for this tree
	:param tree: Tree object with mass and lifetime info
	:param lnc_plus_lnv: if both lnc and lnv decays are possible then lifetime is reduced by a factor of 2
	:return: calculated weight.
	"""
	mass = tree.mass  # GeV
	ctau = tree.ctau  # mm
	mc_campaign = tree.mc_campaign
	channel = tree.channel

	if single_flavour_mixing: 
		if channel == "uuu" or channel == "uue" or channel == "uee":
			x_e = 0
			x_mu = 1
			x_tau = 0
		if channel == "eee" or channel == "eeu" or channel == "euu":
			x_e = 1
			x_mu = 0
			x_tau = 0

	if ih_mixing:
		x_e = 1.0/3.0 
		x_mu = 1.0/3.0
		x_tau = 1.0/3.0
	
	if nh_mixing:
		x_e   = 0.06
    	x_mu  = 0.48
        x_tau = 0.46

		
	if mass == -1 or ctau == -1:  # MC weighting error
		logger.debug("Can't determine the mass and lifetime of signal sample. MC mass-lifetime weight will be set to 1!!")
		return 1

	# define luminosity for the different mc campaigns
	lumi = {'mc16a': 36.20766, 'mc16d': 44.30740, 'mc16e': 58.45010, None: 1.0}
	# lumi_tot = lumi['mc16a'] + lumi['mc16d'] + lumi['mc16e']

	if tree.is_data or tree.not_hnl_mc:  # you are running on data non non-hnl MC
		return 1
	else:  # you are running on a signal MC file


		xsec_LNC_only = hnl_xsec_generic_model(channel = channel, x_e = x_e, x_mu = x_mu, x_tau = x_tau, br_single_flavour_mixing=tree.br, 
											   mass=mass, ctau = ctau, LNC_only = True)  # in fb
		xsec_LNC_plus_LNV = hnl_xsec_generic_model(channel= channel, x_e = x_e, x_mu = x_mu, x_tau = x_tau, br_single_flavour_mixing=tree.br, 
											       mass=mass, ctau = ctau, LNC_only = False)  # in fb

		# mass-lifetime weight = L * HNL_xsec / total num. of MC events
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
		self.dNu_vec = ROOT.TLorentzVector()
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
		self.W_charge = -2
		self.plep_vec = ROOT.TLorentzVector()
		self.plep_charge = -99
		self.mhnl = -1
		self.dvmass = -1
		self.HNL_pdgID = 50
		self.gamma = 1
		self.beta = 1
		self.properLifetime = -1

	def get_truth_particles(self, tree):
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

					self.truth_dv = ROOT.TVector3(self.truth_dvx, self.truth_dvy, self.truth_dvz)
					self.truth_dvr = np.sqrt(self.truth_dvx ** 2 + self.truth_dvy ** 2)
					visTrkVec = ROOT.TLorentzVector()
					truthVec = ROOT.TLorentzVector()
					nu_vec = ROOT.TLorentzVector()

					for i in range(len(tree['truthVtx_outP_pt'][ivx])):
						trk_pdgId = abs(tree['truthVtx_outP_pdgId'][ivx][i])
						if trk_pdgId == 13:
							TrkVec = ROOT.TLorentzVector()
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
							TrkVec = ROOT.TLorentzVector()
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

						if trk_pdgId == 13 or trk_pdgId == 11:  # is track a muon of electron? Then these are our visible (charged) truth tracks
							visTrkVec = ROOT.TLorentzVector()
							visTrkVec.SetPtEtaPhiM(tree['truthVtx_outP_pt'][ivx][i],
												   tree['truthVtx_outP_eta'][ivx][i],
												   tree['truthVtx_outP_phi'][ivx][i],
												   tree['truthVtx_outP_M'][ivx][i]
												   )
							self.trkVec.append(visTrkVec)  # only add visible leptons to trkVec list
						else:  # remaining child is the neutrino
							nu_vec.SetPtEtaPhiM(tree['truthVtx_outP_pt'][ivx][i],
												tree['truthVtx_outP_eta'][ivx][i],
												tree['truthVtx_outP_phi'][ivx][i],
												tree['truthVtx_outP_M'][ivx][i]
												)
							self.dNu_vec = nu_vec

					for i in range(len(tree['truthVtx_outP_pt'][ivx])):
						dLepVec = ROOT.TLorentzVector()
						dLepVec.SetPtEtaPhiM(tree['truthVtx_outP_pt'][ivx][i],
											 tree['truthVtx_outP_eta'][ivx][i],
											 tree['truthVtx_outP_phi'][ivx][i],
											 tree['truthVtx_outP_M'][ivx][i]
											 )
						self.dLepVec.append(dLepVec)  # add all the displaced leptons to one list in the order they are in pythia
						self.dLepCharge.append(tree['truthVtx_outP_charge'][ivx][i])
						# self.dTrk_d0.append(tree['truthVtx_outP_d0'][ivx][i])
						self.dTrk_d0.append(-1)  # fill with -1 for now, default DHNLalg does not have truth d0

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
					self.truth_pv = ROOT.TVector3(self.truth_pvx, self.truth_pvy, self.truth_pvz)

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

	def get_muons(self):
		self.ntracks = self.tree.ntrk
		# print "number of tracks: ", self.ntracks
		for itrk in range(self.ntracks):
			lepVec = ROOT.TLorentzVector()
			std_lepVec = ROOT.TLorentzVector()
			lepmatched_lepVec =  ROOT.TLorentzVector()
			if self.tree.dv('trk_muonIndex')[itrk] >= 0:  # matched muon!
				if self.tree.dv('trk_electronIndex')[itrk] >= 0:  # also matched to an electron!
					# get the muon index
					if len(self.tree['muon_index']) > 0 and self.tree.fake_aod == False:
						muon_index = np.where(self.tree['muon_index'] == self.tree.dv('trk_muonIndex')[itrk])[0][0]
					pass_muon_loose = self.tree['muon_isLoose'][muon_index]
					#If track is NOT matched to a loose muon --> no muon match!
					if not pass_muon_loose == 1:
						# skip tracks matched to muons with no quality!
						continue

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

	def get_electrons(self):
		self.ntracks = self.tree.ntrk

		for itrk in range(self.ntracks):
			lepVec = ROOT.TLorentzVector()
			std_lepVec = ROOT.TLorentzVector()
			lepmatched_lepVec = ROOT.TLorentzVector()
			if self.tree.dv('trk_electronIndex')[itrk] >= 0:  # matched electron!
				# Check if track is also matched to a muon!
				if self.tree.dv('trk_muonIndex')[itrk] >= 0:
					# get the muon and electron indicies
					if len(self.tree['muon_index']) > 0 and self.tree.fake_aod == False:
						muon_index = np.where(self.tree['muon_index'] == self.tree.dv('trk_muonIndex')[itrk])[0][0]
					if len(self.tree['el_index']) > 0:
						el_index = np.where(self.tree['el_index'] == self.tree.dv('trk_electronIndex')[itrk])[0][0]
					pass_muon_loose = self.tree['muon_isLoose'][muon_index]
					pass_electron_vvl = self.tree['el_isLHVeryLoose_mod1'][el_index]
					# If track is matched to a loose muon --> no electron match!
					if pass_muon_loose == 1:
						# skip tracks matched to loose muons
						continue
					# else if track is NOT vvl electron --> no electron match!
					elif not pass_electron_vvl == 1:
						# skip track with no electron quality!
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
					self.el_isTight.append(self.tree.dv('trk_isTight')[itrk])  # add muon quality info
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

	def get_tracks(self, idv=-1):
		"""Fills the Track object with a collection of track vectors.
		No lepton matching enforced."""
		at_idv = self.tree.idv if idv < 0 else idv
		at_ievt = self.tree.ievt
		prefix = self.tree.dv_prefix + '_'
		self.ntracks = self.tree.ntrk
		for itrk in range(self.tree.ntrk):
			trkvec = ROOT.TLorentzVector()
			pt = self.tree.get_at(prefix + 'trk_pt_wrtSV', at_ievt, at_idv, itrk)
			eta = self.tree.get_at(prefix + 'trk_eta_wrtSV', at_ievt, at_idv, itrk)
			phi = self.tree.get_at(prefix + 'trk_phi_wrtSV', at_ievt, at_idv, itrk)
			M = self.tree.get_at(prefix + 'trk_M', at_ievt, at_idv, itrk)

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
		self.get_muons()


class Electrons(Tracks, object):
	def __init__(self, tree):
		super(Electrons, self).__init__(tree)
		self.get_electrons()


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

		sig_info = MCInfo(self.dsid)
		self.file_ch = sig_info.ch_str  # used in job submission do not delete!

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

		f_br = ReadBRdat(os.path.dirname(os.path.abspath(__file__)) + '/../data/BR/HNL_branching_20GeV.dat')
		self.br = f_br.get_BR(channel, self.mass)


class MCInfo:
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
# separated by year using comments from the above link and cross checking with this twiki: https://twiki.cern.ch/twiki/bin/view/Atlas/LowestUnprescaled

# Single muon triggers used in DHNL analysis
SingleMuonTriggerlist = [
	"HLT_mu26_ivarmedium",  # 2017-2018
	"HLT_mu20_iloose_L1MU15",  # 2015
]

SingleMuonTriggerlist_2018 = [
	"HLT_mu26_ivarmedium",
]

SingleMuonTriggerlist_2017 = [
	"HLT_mu26_ivarmedium",
]

SingleMuonTriggerlist_2015_2016 = [
	"HLT_mu26_ivarmedium",
	"HLT_mu20_iloose_L1MU15",
]

# Single electron triggers used in DHNL analysis
SingleElectronTriggerlist = [
	"HLT_e24_lhmedium_L1EM20VH",  # 2015
	"HLT_e26_lhtight_nod0_ivarloose",  # 2016-2018
]

SingleElectronTriggerlist_2018 = [
	"HLT_e26_lhtight_nod0_ivarloose",
]

SingleElectronTriggerlist_2017 = [
	"HLT_e26_lhtight_nod0_ivarloose",
]

SingleElectronTriggerlist_2015_2016 = [
	"HLT_e24_lhmedium_L1EM20VH",
	"HLT_e26_lhtight_nod0_ivarloose",
]


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

def pT_diff(track_pT_1 , track_pT_2):
	abs_diff = abs(track_pT_1- track_pT_2)
	return abs_diff/track_pT_1
