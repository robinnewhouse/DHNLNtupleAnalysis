import ROOT
import numpy as np
import sys
import os
import time
import json

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


class ReadJsonFiles:
	def __init__(self, filename):
		with open(filename) as fp:
			self.json_string = json.load(fp)

	def get_coupling(self, mass, ctau, model, use_str=False):
		if use_str:
			mass_str = mass
			ctau_str = ctau
		else:
			mass_str = str(mass)
			ctau_str = '{ctau} mm'.format(ctau = int(ctau))
		coupling = self.json_string[mass_str][ctau_str][model]
		return coupling

	def get_br(self, channel, mass, model, use_str=False):
		if use_str: mass_str = mass
		else: mass_str = str(mass)

		if   channel == "uuu": decay_str = "mmv"
		elif channel == "uue": decay_str = "mev"
		elif channel == "ueu": decay_str = "emv"
		elif channel == "uee": decay_str = "eev"
		elif channel == "eee": decay_str = "eev"
		elif channel == "eeu": decay_str = "emv"
		elif channel == "eue": decay_str = "mev"
		elif channel == "euu": decay_str = "mmv"

		br = self.json_string[mass_str][decay_str][model]
		return br


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
				if channel == 'uee':
					return self.BReee[i]
				if channel == 'euu':
					return self.BRuuu[i]


class MCEventWeight:
	def __init__(self, tree, mixing_type, dirac_limit=False, flip_e_and_mu=False, use_gronau=False):
		"""
		Class use for computing MC event weights and HNL cross sections.
		@parm tree: ntuple analysis tree
		@parm mixing_type: HNL model to be used in the computation: single-flavour, IH or NH
		@parm  dirac_limit: Dirac limit means that only LNC decays are allowed, LNV processes are suppressed to zero. 
						    Default: False means that Majorana limit is used with 50% mix of LNC and LNV decays.
		@parm flip_e_and_mu: flip order of e and mu in channel name
		@parm use_gronau: use Gronau parametrization (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.29.2539) in the computation of U2 and br. Otherwise Oleg's parametrization will be used
		"""
		self.tree = tree
		self.use_gronau = use_gronau
		self.flip_e_and_mu = flip_e_and_mu
		self.mixing_type = mixing_type
		self.dirac_limit = dirac_limit
		self.channel = tree.channel
		# Overwrite channel name if flip_e_and_mu is true!
		if self.channel == "uue" and self.flip_e_and_mu: self.channel = "ueu"
		if self.channel == "eeu" and self.flip_e_and_mu: self.channel = "eue"
		# Name the single flavour mixing model depending on the channel
		if self.mixing_type == "single-flavour":
			if self.channel == "uuu" or self.channel == "uue" or self.channel == "uee":
				self.mixing_type = "mu_only"
			if self.channel == "eee" or self.channel == "eeu" or self.channel == "euu":
				self.mixing_type = "e_only"

		# ###########################################################
		# Define the model dependent coupling fractions
		# x_alpha = U_alpha^2 / U_tot^2 such that alpha = e, mu, tau
		# ###########################################################
		if self.mixing_type == "mu_only":
			self.x_e = 0
			self.x_mu = 1
			self.x_tau = 0
		if self.mixing_type == "e_only":
			self.x_e = 1
			self.x_mu = 0
			self.x_tau = 0
		if self.mixing_type == "IH":
			self.x_e = 1.0 / 3.0
			self.x_mu = 1.0 / 3.0
			self.x_tau = 1.0 / 3.0
		if self.mixing_type == "NH":
			self.x_e = 0.06
			self.x_mu = 0.48
			self.x_tau = 0.46

		# Define the prod and decay coupling ratios depending on the channel,
		# where x_alpha = U_alpha^2 / U_tot^2 such that alpha = e, mu, tau
		if self.channel == 'uuu' or self.channel == 'uue':
			self.x_prod = self.x_mu
			self.x_decay = self.x_mu
		if self.channel == 'eee' or self.channel == 'eeu':
			self.x_prod = self.x_e
			self.x_decay = self.x_e
		if self.channel == 'ueu' or self.channel == 'uee':
			self.x_prod = self.x_mu
			self.x_decay = self.x_e
		if self.channel == 'eue' or self.channel == 'euu':
			self.x_prod = self.x_e
			self.x_decay = self.x_mu

		if not self.tree.is_data:
			# ######################################################################################################################################################
			# Get branching ratios (BR). BR depend on the mass, decay mode and model
			# BR json files contain a dictionaries of BR[mass][decay][model] for e_only, mu_only, NH and IH models.
			# ######################################################################################################################################################
			br_root_path = os.path.dirname(os.path.abspath(__file__)) + '/../data/BR/'
			#  Computes lifetime via Gronau formulas and thus has 10-20% difference
			if self.use_gronau: f_br = ReadJsonFiles(br_root_path + 'BranchingRatios_DifferentMixings_Gronau_lifetime.json')
			#  BranchingRatios_DifferentMixings_Olegs_lifetime.json computes lifetime via Oleg's parametrization (better agreement with MG)
			else: f_br = ReadJsonFiles(br_root_path + 'BranchingRatios_DifferentMixings_Olegs_lifetime.json')
			self.br = f_br.get_br(self.channel, self.tree.mass, self.mixing_type)

			# ######################################################################################################################################################
			# Get coupling squared. Coupling depends on the mass, lifetime and model
			# Theta2 json files contain a dictionaries of BR[mass][lifetime][model] for e_only, mu_only, NH and IH models.
			# ######################################################################################################################################################
			#  Theta2 files contain a dictionaries of coupling2[mass][ctau][model] for signal-flavour mixing, NH and IH models.
			if self.use_gronau: f_coupling = ReadJsonFiles(br_root_path + '/Theta2_DifferentMixings_Gronau_lifetime.json')
			else: f_coupling = ReadJsonFiles(br_root_path + '/Theta2_DifferentMixings_Olegs_lifetime.json')
			self.U2 = f_coupling.get_coupling(self.tree.mass, self.tree.ctau, self.mixing_type)
		else:
			self.U2 = -1
			self.br = -1


	def get_mc_event_weight(self):
		"""
		Calculates the MC event weight as luminosity x cross section / total number of MC events
		Cross section is model dependent an therefore the mc event weight is also model dependent
		@return: calculated weight.
		"""
		if self.tree.mass == -1 or self.tree.ctau == -1:  # MC weighting error
			logger.debug("Can't determine the mass and lifetime of signal sample. MC mass-lifetime weight will be set to 1!!")
			return 1

		# #############################################################################
		# Get the luminosity for the different mc campaigns
		# lumi_tot = lumi['mc16a'] + lumi['mc16d'] + lumi['mc16e']
		# #############################################################################
		lumi = {'mc16a': 36.10416, 'mc16d': 44.30740, 'mc16e': 58.45010, None: 1.0}

		if self.tree.is_data or self.tree.not_hnl_mc:  # Running on data non non-hnl MC
			return 1 # mc event weight equals 1
		else:  # Running on an HNL signal file
			# #############################################################################
			# Get the cross sections for the different HNL models
			# #############################################################################
			# One HNL Majorana model
			one_hnl_majorana_hnl_xsec = self.hnl_xsec_generic_model(channel = self.channel,mass=self.tree.mass, ctau = self.tree.ctau,  
																	x_e = self.x_e, x_mu = self.x_mu, x_tau = self.x_tau )  # in fb
			# HNL Dirac models only have LNC decays. LNC rates are coherently enhanced by a factor of 2 compared to Majorana HNLs
			xsec_one_hnl_dirac = one_hnl_majorana_hnl_xsec * 2
			# Quasi-Dirac with "Majorana limit" model with 50/50 LNC/LNV decays has two Majorana particles mediating the process, enhances cross sections by a factor of 2 compared to on Majorana rates
			xsec_quasi_dirac_pair_majorana_limit = one_hnl_majorana_hnl_xsec * 2
			# Quasi-Dirac with "Dirac limit" model has only LNC decays and two Majorana particles mediating the process, enhances cross sections by a factor of 4 compared to on Majorana rates
			# One factor of 2 is from the quasi-dirac pair and one factor of two becuase of the change in lifetime/ coupling
			xsec_quasi_dirac_pair_dirac_limit = one_hnl_majorana_hnl_xsec * 4
			if self.mixing_type == "NH" or self.mixing_type == "IH":
				if self.dirac_limit: hnl_xsec = xsec_quasi_dirac_pair_dirac_limit
				else: hnl_xsec = xsec_quasi_dirac_pair_majorana_limit

			if self.mixing_type == "e_only" or self.mixing_type == "mu_only":
				if self.dirac_limit: hnl_xsec = xsec_one_hnl_dirac
				else: hnl_xsec = one_hnl_majorana_hnl_xsec

			# Compute the cross section for LNC or LNV decay process. Pythia8 samples have a 50% mix of LNC+ LNV of the number of LNC or LNV events
			# Thus, the number of mc generated LNC or LNV events is equal to "all_entries / 2"
			n_mc_events = self.tree.all_entries / 2

			# Compute MC event weight as "L * xsec / total num. of MC events"
			weight = lumi[self.tree.mc_campaign] * hnl_xsec / n_mc_events

		return weight
	
	def hnl_xsec_generic_model(self, channel, mass, ctau , x_e = 1, x_mu = 0, x_tau = 0, use_U2_input = True, majorana_particle = True ):
		"""
		Calculates the HNL cross section given a br, mass, ctau
		Cross section is model dependent an therefore the mc event weight is also model dependent

		@param channel: HNL channel
		@param mass: HNL mass in GeV
		@param mass: HNL lifetime in mm
		@param x_e: U_e^2 / U_tot ^2
		@param x_mu: U_mu^2 / U_tot ^2
		@param x_tau: U_tau^2 / U_tot ^2
		@param use_U2_input use self.U2 read from the json file as the coupling squared
		
		@parm  dirac_particle: Dirac particle nature, otherwise assume Majorana
		@parm flip_e_and_mu: flip order of e and mu in channel name
		@return: calculated weight.
		"""

		# Get lifetime either from class attribute or compute it by hand using the Gronau parametrization.
		# If you are using the Gronau json file to get the couplings and majorana_particle= True, then these two methods are identical
		if use_U2_input: U2 = self.U2
		else: U2_for_xsec = U2 = self.get_Gronau_U2(mass, ctau, majorana_particle)

		# ########################################################################################################################
		# Calculate the HNL cross section
		# ########################################################################################################################
		# constants used in cross section computation
		mW = 80.379  # mass of W boson in GeV from PDG (https://pdg.lbl.gov/2019/listings/rpp2019-list-w-boson.pdf)
		xsec_pp_times_br_W_to_l_nu = 20.6e6 # in fb from ATLAS W cross section measurement (https://arxiv.org/abs/1603.09222)
		hnl_xsec =  xsec_pp_times_br_W_to_l_nu * ((1 - (mass / mW) ** 2) ** 2) * (1 + (mass ** 2) / (2 * mW ** 2))  * self.x_prod *  U2  * self.br  # in fb

		return hnl_xsec
	
	def get_Gronau_U2(self,mass, ctau, x_e, x_mu, x_tau, majorana_particle):
		''''
		Computes U2 for a given mass and lifetime and model parameters (x_e, x_mu, x_tau) using Gronau parametrization (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.29.2539)
		'''
		# Compute coupling using the Gronau parametrization
		# Model constants from Gronau
		tau_0_e = 4.15e-12 # in seconds
		tau_0_mu = 4.49e-12  # in seconds
		tau_0_tau = 1.08e-11  # in seconds
		b_e = 5.17
		b_mu = 5.19
		b_tau = 5.44
		# For Majorana model, both LNC and LNV decays are allowed, which means twice as many decay channels are opened and, as a result, the lifetime is reduced by a factor of 2.
		if majorana_particle: 
			tau_0_e = tau_0_e * 0.5
			tau_0_mu = tau_0_mu * 0.5
			tau_0_tau = tau_0_tau * 0.5
		e_term   = x_e   * mass ** b_e / (tau_0_e * 3e8) # in 1/m 
		mu_term  = x_mu  * mass ** b_mu / (tau_0_mu * 3e8) # in 1/m 
		tau_term = x_tau * mass ** b_tau / (tau_0_tau * 3e8) # in 1/m 
		U2Gronau = (1 / (ctau / 1000))*(1 / (mu_term + e_term + tau_term )) # unitless (ctau is in mm )
		
		return U2Gronau


class Truth:
	def __init__(self):
		self.HNL_vec = ROOT.TLorentzVector()
		self.dNu_vec = ROOT.TLorentzVector()
		self.trkVec = []
		self.dLepVec = []
		self.dLep_pdgID = []
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
						dLep_pdgID =  abs(tree['truthVtx_outP_pdgId'][ivx][i])
						self.dLep_pdgID.append(dLep_pdgID)
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


def get_lepton_index(tree, itrk, lepton_type):
	"""
	Retrieves the lepton index based on the associated track
	@param tree: input TTree
	@param itrk: index of track to be matched
	@param lepton_type: 'muon' or 'electron'
	@return: lepton index or None if index not found
	"""
	lepton_index = None
	if lepton_type == 'muon':
		lepton_index = np.where(tree['muon_index'] == tree.dv('trk_muonIndex')[itrk])
	if lepton_type == 'electron':
		lepton_index = np.where(tree['el_index'] == tree.dv('trk_electronIndex')[itrk])
	if len(lepton_index[0]) > 0:
		if tree['muon_phi' if lepton_type == 'muon' else 'el_phi'][lepton_index[0][0]] - tree.dv('trk_phi')[itrk] > 0.02:
			logger.error("Lepton and track phi to not match. Check index counting. phi_lep: {}, phi_track: {}".format(
				tree['muon_phi' if lepton_type == 'muon' else 'el_phi'][lepton_index[0][0]], tree.dv('trk_phi')[itrk]))
		return lepton_index[0][0]
	else:
		return None


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
			lepmatched_lepVec = ROOT.TLorentzVector()
			if self.tree.dv('trk_muonIndex')[itrk] >= 0:  # matched muon!
				if self.tree.dv('trk_electronIndex')[itrk] >= 0:  # also matched to an electron!
					# get the muon index
					if len(self.tree['muon_index']) > 0 and not self.tree.fake_aod:
						muon_index = get_lepton_index(self.tree, itrk, 'muon')
						if muon_index is None: continue
					pass_muon_loose = self.tree['muon_isLoose'][muon_index]
					# If track is NOT matched to a loose muon --> no muon match!
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
						muon_index = get_lepton_index(self.tree, itrk, 'muon')
						if muon_index is None: continue
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
					if len(self.tree['muon_index']) > 0 and not self.tree.fake_aod:
						muon_index = get_lepton_index(self.tree, itrk, 'muon')
						if muon_index is None: continue
					if len(self.tree['el_index']) > 0:
						el_index = get_lepton_index(self.tree, itrk, 'electron')
						if el_index is None: continue
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
						el_index = get_lepton_index(self.tree, itrk, 'electron')
						if el_index is None: continue
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

		self.mc_campaign = None
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

		if "_2p5G" in infile or sig_info.mass_str == "2p5G":
			self.mass = 2.5
			self.mass_str = "2p5G"
		elif "_3G" in infile or sig_info.mass_str == "3G":
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
			self.mc_campaign = "mc16a"
		if "r11916" in infile or "mc16d" in infile or "r10739" in infile:
			self.mc_campaign = "mc16d"
		if "r11891" in infile or "mc16e" in infile or "r10790" in infile:
			self.mc_campaign = "mc16e"
		logger.info("This mc campaign is: {}".format(self.mc_campaign))

		# More flexibility for non-signal samples
		self.output_filename = "histograms"
		if self.mc_campaign:
			self.output_filename += "_" + self.mc_campaign
		else:
			self.output_filename += "_mc"
		if self.mass_str: self.output_filename += "_" + self.mass_str
		if self.ctau_str: self.output_filename += "_" + self.ctau_str
		self.output_filename += "_" + channel + ".root"


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
		mc_info[313419] = ["uee", "3G", "lt1dd"]
		mc_info[313420] = ["uee", "3G", "lt10dd"]
		mc_info[313421] = ["uee", "3G", "lt100dd"]
		mc_info[313422] = ["euu", "3G", "lt1dd"]
		mc_info[313423] = ["euu", "3G", "lt10dd"]
		mc_info[313424] = ["euu", "3G", "lt100dd"]
		mc_info[313425] = ["uee", "4G", "lt1dd"]
		mc_info[313426] = ["uee", "4G", "lt10dd"]
		mc_info[313427] = ["uee", "4G", "lt100dd"]
		mc_info[313428] = ["euu", "4G", "lt1dd"]
		mc_info[313429] = ["euu", "4G", "lt10dd"]
		mc_info[313430] = ["euu", "4G", "lt100dd"]
		mc_info[313431] = ["uee", "4p5G", "lt1dd"]
		mc_info[313432] = ["uee", "4p5G", "lt10dd"]
		mc_info[313433] = ["uee", "4p5G", "lt100dd"]
		mc_info[313434] = ["euu", "4p5G", "lt1dd"]
		mc_info[313435] = ["euu", "4p5G", "lt10dd"]
		mc_info[313436] = ["euu", "4p5G", "lt100dd"]
		mc_info[313437] = ["uee", "5G", "lt1dd"]
		mc_info[313438] = ["uee", "5G", "lt10dd"]
		mc_info[313439] = ["uee", "5G", "lt100dd"]
		mc_info[313440] = ["euu", "5G", "lt1dd"]
		mc_info[313441] = ["euu", "5G", "lt10dd"]
		mc_info[313442] = ["euu", "5G", "lt100dd"]
		mc_info[313443] = ["uee", "7p5G", "lt1dd"]
		mc_info[313444] = ["uee", "7p5G", "lt10dd"]
		mc_info[313445] = ["uee", "7p5G", "lt100dd"]
		mc_info[313446] = ["euu", "7p5G", "lt1dd"]
		mc_info[313447] = ["euu", "7p5G", "lt10dd"]
		mc_info[313448] = ["euu", "7p5G", "lt100dd"]
		mc_info[313449] = ["uee", "10G", "lt1dd"]
		mc_info[313450] = ["uee", "10G", "lt10dd"]
		mc_info[313451] = ["uee", "10G", "lt100dd"]
		mc_info[313452] = ["euu", "10G", "lt1dd"]
		mc_info[313453] = ["euu", "10G", "lt10dd"]
		mc_info[313454] = ["euu", "10G", "lt100dd"]
		mc_info[313455] = ["uee", "12p5G", "lt1dd"]
		mc_info[313456] = ["uee", "12p5G", "lt10dd"]
		mc_info[313457] = ["uee", "12p5G", "lt100dd"]
		mc_info[313458] = ["euu", "12p5G", "lt1dd"]
		mc_info[313459] = ["euu", "12p5G", "lt10dd"]
		mc_info[313460] = ["euu", "12p5G", "lt100dd"]
		mc_info[313461] = ["uee", "15G", "lt1dd"]
		mc_info[313462] = ["uee", "15G", "lt10dd"]
		mc_info[313463] = ["uee", "15G", "lt100dd"]
		mc_info[313464] = ["euu", "15G", "lt1dd"]
		mc_info[313465] = ["euu", "15G", "lt10dd"]
		mc_info[313466] = ["euu", "15G", "lt100dd"]
		mc_info[313467] = ["uee", "17p5G", "lt1dd"]
		mc_info[313468] = ["uee", "17p5G", "lt10dd"]
		mc_info[313469] = ["uee", "17p5G", "lt100dd"]
		mc_info[313470] = ["euu", "17p5G", "lt1dd"]
		mc_info[313471] = ["euu", "17p5G", "lt10dd"]
		mc_info[313472] = ["euu", "17p5G", "lt100dd"]
		mc_info[313473] = ["uee", "20G", "lt1dd"]
		mc_info[313474] = ["uee", "20G", "lt10dd"]
		mc_info[313475] = ["uee", "20G", "lt100dd"]
		mc_info[313476] = ["euu", "20G", "lt1dd"]
		mc_info[313477] = ["euu", "20G", "lt10dd"]
		mc_info[313478] = ["euu", "20G", "lt100dd"]
		mc_info[313479] = ["uue", "2p5G", "lt1dd"]
		mc_info[313480] = ["uue", "2p5G", "lt10dd"]
		mc_info[313481] = ["uue", "2p5G", "lt100dd"]
		mc_info[313482] = ["eeu", "2p5G", "lt1dd"]
		mc_info[313483] = ["eeu", "2p5G", "lt10dd"]
		mc_info[313484] = ["eeu", "2p5G", "lt100dd"]
		mc_info[313485] = ["utt", "10G", "lt1dd"]
		mc_info[313486] = ["utt", "10G", "lt10dd"]
		mc_info[313487] = ["utt", "10G", "lt100dd"]
		mc_info[313488] = ["ett", "10G", "lt1dd"]
		mc_info[313489] = ["ett", "10G", "lt10dd"]
		mc_info[313490] = ["ett", "10G", "lt100dd"]

		if dsid is None:
			logger.warning("No dsid")
			self.mass_str = None
			self.ctau_str = None
			self.ch_str = None
		else:
			pmuon_dsid = dsid in range(311602, 311661 + 1) or dsid in range(313482, 313484 + 1)
			pel_dsid = dsid in range(312956, 313015 + 1) or dsid in range(313479, 313481 + 1)
			mixed_coupling_dsid = dsid in range(313419, 313490 + 1) and not dsid in range(313479, 313484 + 1)

			if pmuon_dsid or pel_dsid or mixed_coupling_dsid:
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
