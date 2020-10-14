# from ROOT import*
import ROOT
import numpy as np
import helpers
import logging

#make a global logger variable for all selection classes
logger = helpers.getLogger('dHNLAnalysis.selections', level=logging.WARNING)

def set_debug_level(level): 
	#update debug level to match debug level set at run time
	logger = helpers.getLogger('dHNLAnalysis.selections', level=level)


class Trigger():
	def __init__(self, tree, trigger, invert=False):
		self.event_triggers = tree['passedTriggers']
		self.invert = invert
		if trigger == "muononly":
			if tree.mc_campaign == "mc16a":
				self.allowed_trigger_list = helpers.apiSingleMuonTriggerlist_2015_2016
			if tree.mc_campaign == "mc16d":
				self.allowed_trigger_list = helpers.apiSingleMuonTriggerlist_2017
			if tree.mc_campaign == "mc16e":
				self.allowed_trigger_list = helpers.apiSingleMuonTriggerlist_2018
			else:
				self.allowed_trigger_list = helpers.apiSingleMuonTriggerlist
		elif trigger == "electrononly":
			if tree.mc_campaign == "mc16a":
				self.allowed_trigger_list = helpers.apiSingleElectronTriggerlist_2015_2016
			if tree.mc_campaign == "mc16d":
				self.allowed_trigger_list = helpers.apiSingleElectronTriggerlist_2017
			if tree.mc_campaign == "mc16e":
				self.allowed_trigger_list = helpers.apiSingleElectronTriggerlist_2018
			else:
				self.allowed_trigger_list = helpers.apiSingleElectronTriggerlist
		elif trigger == "all":
			if tree.mc_campaign == "mc16a":
				self.allowed_trigger_list = helpers.apiSingleMuonTriggerlist_2015_2016 + helpers.apiSingleElectronTriggerlist_2015_2016
			if tree.mc_campaign == "mc16d":
				self.allowed_trigger_list = helpers.apiSingleMuonTriggerlist_2017 + helpers.apiSingleElectronTriggerlist_2017
			if tree.mc_campaign == "mc16e":
				self.allowed_trigger_list = helpers.apiSingleMuonTriggerlist_2018 + helpers.apiSingleElectronTriggerlist_2018
			else:
				self.allowed_trigger_list = helpers.apiSingleMuonTriggerlist + helpers.apiSingleElectronTriggerlist

		elif trigger == "DAOD_RPVLL":
			self.allowed_trigger_list = helpers.DAOD_RPVLLTriggerlist
		else:
			self.allowed_trigger_list = list(trigger)

	def overlap(self, event_triggers, trigger_list):
		"""
		Evaluates whether the event triggers are found in the given trigger list
		https://stackoverflow.com/questions/3170055/test-if-lists-share-any-items-in-python
		overlap([a,b,c], [x,y,z]) # False
		overlap([a,b,c], [b,c,d]) # True
		"""
		return not set(helpers.decode_list(event_triggers)).isdisjoint(trigger_list)

	def passes(self):
		if self.invert:  # invert trigger requirement
			# check if event_triggers includes a trigger in the DAOD_RPVLL trigger list
			if self.overlap(self.event_triggers, helpers.DAOD_RPVLLTriggerlist):
				# If the event passes a DAOD_RPVLL trigger, return the inversion of the standard trigger requirement
				return not self.overlap(self.event_triggers, self.allowed_trigger_list)
			else:
				return False
		else:
			# default; check if event_triggers includes a trigger on the allowed trigger list
			return self.overlap(self.event_triggers, self.allowed_trigger_list)


class Filter():
	def __init__(self, tree, filter_type):
		self.passes_filter = False

		if filter_type == "mu-mu":
			self.passes_filter = tree['passesHnlMuMuFilter']

		if filter_type == "mu-el":
			self.passes_filter = tree['passesHnlMuElFilter']

		if filter_type == "el-mu":
			self.passes_filter = tree['passesHnlElMuFilter']

		if filter_type == "el-el":
			self.passes_filter = tree['passesHnlElElFilter']

		if filter_type == "4-filter":
			self.passes_filter = (tree['passesHnlMuMuFilter']
					or tree['passesHnlElMuFilter']
					or tree['passesHnlElElFilter']
					or tree['passesHnlMuElFilter'])

		if filter_type == "3-filter":
			self.passes_filter = (tree['passesHnlMuMuFilter']
					or tree['passesHnlElMuFilter']
					or tree['passesHnlElElFilter'])

		if filter_type == "2-filter":
			self.passes_filter = (tree['passesHnlMuMuFilter']
					or tree['passesHnlElMuFilter'])

		if filter_type == "1-filter":
			self.passes_filter = tree['passesHnlMuMuFilter']

	def passes(self):
		return self.passes_filter

class InvertedPromptLepton():
	def __init__(self, tree, d0_cut=3.0, z0_sin_theta_cut=0.5):
		self.n_prompt_leptons = 0
		self.n_prompt_muons = 0 
		self.n_prompt_electrons = 0 
		n_muons = len(tree['muon_pt'])
		n_electrons = len(tree['el_pt'])

		for imu in range(n_muons):
			# check the muon has some quality (at least loose)
			if not ((tree['muon_isTight'][imu] == 1) or
					(tree['muon_isMedium'][imu] == 1) or
					(tree['muon_isLoose'][imu]) == 1):
				# muon doesn't satisfy any quality, ignore it
				continue

			mu_d0 = tree['muon_trkd0'][imu]
			mu_z0_sin_theta = tree['muon_trkz0sintheta'][imu]
			# check that the muon satisfies prompt lepton requirements
			if abs(mu_d0) < d0_cut and abs(mu_z0_sin_theta) < z0_sin_theta_cut:
				self.n_prompt_muons += 1

		for iel in range(n_electrons):
			# make sure the electron is at least loose
			if not ((tree['el_LHTight'][iel] == 1) or
					(tree['el_LHMedium'][iel] == 1) or
					(tree['el_LHLoose'][iel]) == 1):
				# electron doesn't satisfy any quality, ignore it
				continue

			el_d0 = tree['el_trkd0'][iel]
			el_z0_sin_theta = tree['el_trkz0sintheta'][iel]
			if (abs(el_d0) < d0_cut) and (abs(el_z0_sin_theta) < z0_sin_theta_cut):
				# electron satisfies prompt lepton requirements
				self.n_prompt_electrons += 1

		self.n_prompt_leptons = self.n_prompt_electrons + self.n_prompt_muons

	def passes(self):
		return self.n_prompt_leptons == 0


class PromptLepton():
	def __init__(self, tree, lepton="any", quality="tight", min_dR=0.05):
		self.plepVec = ROOT.TLorentzVector(0,0,0,0)
		self.plepcharge = 0
		self.plepd0 = -2000
		self.plepz0 = -2000
		self.nPlep = 0


		lepquality = ""
		passPfilter = False
		if lepton == "muon":
			if quality == "tight":  # tight muon is requested
				lepquality = 'muon_isTight'
			if quality == "medium":
				lepquality = 'muon_isMedium'
			if quality == "loose":
				lepquality = 'muon_isLoose'

			nleps = len(tree['muon_pt'])
			passPfilter = tree['muon_passesPromptCuts']

		if lepton == "electron":
			if quality == "tight":  # tight electron is requested
				lepquality = 'el_LHTight'
			if quality == "medium":
				lepquality = 'el_LHMedium'
			if quality == "loose":
				lepquality = 'el_LHLoose'


			nleps = len(tree['el_pt'])
			passPfilter = tree['el_passesPromptCuts']

		# variable for the highest pt lepton				
		self.highestpt_lep = ROOT.TLorentzVector(0, 0, 0, 0)
		self.highestpt_lep_d0 = -2000
		self.highestpt_lep_z0 = -2000

		for ilep in range(nleps): 
			overlap = False
			plepVec_i = ROOT.TLorentzVector()

			if lepton == "muon":
				pt = tree['muon_pt'][ilep]
				eta = tree['muon_eta'][ilep]
				phi = tree['muon_phi'][ilep]
				mass = tree['muon_m'][ilep]
				charge = tree['muon_charge'][ilep]
				plepVec_i.SetPtEtaPhiM(pt, eta, phi, mass)

				lepd0 = tree['muon_trkd0'][ilep]
				lepz0 = tree['muon_trkz0'][ilep]
				lepz0sintheta = tree['muon_trkz0sintheta'][ilep]

			if lepton == "electron":
				pt = tree['el_pt'][ilep]
				eta = tree['el_eta'][ilep]
				phi = tree['el_phi'][ilep]
				mass = tree['el_m'][ilep]
				charge = tree['el_charge'][ilep]
				plepVec_i.SetPtEtaPhiM(pt, eta, phi, mass)

				lepd0 = tree['el_trkd0'][ilep]
				lepz0 = tree['el_trkz0'][ilep]
				lepz0sintheta = tree['el_trkz0sintheta'][ilep]

			# check if the plep passes the DRAW filter and passes quality before looping over tracks
			# changed to be careful with negative electron quality values # RN
			passes_lep_quality = lepquality == "" or tree[lepquality][ilep] > 0
			
			#apply filter cut seperately do not need to addtionally require our lepton passes the prompt filter cuts -DT
			# if passPfilter[ilep] and passes_lep_quality and abs(lepd0) < 3 and abs(lepz0sintheta) < 0.5:
			if passes_lep_quality and abs(lepd0) < 3 and abs(lepz0sintheta) < 0.5:
				# Check the overlap between the prompt lepton and every displaced vertex track
				for idv in range(tree.ndv):
					prefix = tree.dv_prefix + '_'
					ntrks = tree.get_at(prefix+'ntrk', tree.ievt, idv)
					for itrk in range(ntrks):
						# Currently the only live example of get_at() which gives full control over the tree access.
						pt = tree.get_at(prefix + 'trk_pt', tree.ievt, idv, itrk)
						eta = tree.get_at(prefix + 'trk_eta', tree.ievt, idv, itrk)
						phi = tree.get_at(prefix + 'trk_phi', tree.ievt, idv, itrk)
						M = tree.get_at(prefix + 'trk_M', tree.ievt, idv, itrk)
						track_vector = ROOT.TLorentzVector()
						track_vector.SetPtEtaPhiM(pt, eta, phi, M)

						dR = track_vector.DeltaR(plepVec_i)
						if dR < min_dR:  # set overlap to true if muon overlaps with displaced track
							overlap = True
		
				# if lepton doesnt overlap with and DV tracks
				if not overlap:
					self.nPlep = self.nPlep + 1
					# if pt is larger then the previous prompt lepton found
					if plepVec_i.Pt() > self.highestpt_lep.Pt():
						self.highestpt_lep = plepVec_i  # get highest pt prompt lepton!
						self.highestpt_lep_d0 = lepd0
						self.highestpt_lep_z0 = lepz0
						self.highestpt_lep_charge = charge
						self.highestpt_lep_z0sintheta = lepz0sintheta

				#for trigger matching
						# if tree["muon_isTrigMatched"][ilep] == 0:
						# 	print "is muon trig matched?", tree["muon_isTrigMatched"][ilep]
						# 	print "pt of the highest pt TIGHT muon: ", self.highestpt_plep.Pt() 
						# 	print "muon pt: ", self.evt.tree.muonpt[self.evt.ievt]
						# 	print "trigger matched: ", tree["muon_isTrigMatched"]
						# 	print "lepton quality: ", lepquality[self.evt.ievt]
						# 	print "muon type: ", self.evt.tree.muontype[self.evt.ievt]

	def passes(self):
		# check if you found a prompt lepton
		if self.nPlep > 0 and self.highestpt_lep.Pt() > 0:
			self.plepVec = self.highestpt_lep
			self.plepd0 = self.highestpt_lep_d0
			self.plepz0 = self.highestpt_lep_z0
			self.plepcharge = self.highestpt_lep_charge
			return True
		else: 
			return False
		
class Alpha():
	def __init__(self, tree, max_alpha=0.01):
		self.max_alpha = max_alpha
		# calculate alpha
		# secondary vertex position vector
		sv_vector = ROOT.TVector3(tree.dv('x'),
								  tree.dv('y'),
								  tree.dv('z'))
		# primary vertex position vector
		pv_vector = ROOT.TVector3(tree['vertex_x'],
								  tree['vertex_y'],
								  tree['vertex_z'])
		# vector from pv to sv
		pv_sv_vector = sv_vector - pv_vector

		# vector difference between momentum vector and position vector		
		alpha = pv_sv_vector.Phi() - tree.dv('phi')
		# put in -pi to pi range
		self.alpha = (alpha + np.pi/2) % np.pi*2 - np.pi

	def passes(self):
		# pass if alpha is less than the sent cut
		return abs(self.alpha) < self.max_alpha


class DVradius():
	def __init__(self, tree):
		self.rdv = -1
		if tree.ntrk > 0:
			dx = tree.dv('x')
			dy = tree.dv('y')
			self.rdv = np.sqrt(dx**2 + dy**2)

	def passes(self, rdv_min=4, rdv_max=300):
		if self.rdv > rdv_min and self.rdv < rdv_max:
			return True
		else:
			return False


class DVntracks():
	def __init__(self, tree, ntrk=2, decaymode="leptonic"):
		self.ntrk = ntrk
		self.decaymode = decaymode

		self.ntracks = -1 

		if self.decaymode == "leptonic":
			self.ntracks = tree.ntrk

	def passes(self):
		if self.ntracks == self.ntrk:
			return True
		else:
			return False


class ChargeDV():
	def __init__(self, tree, sel="OS", decaymode="leptonic",trk_charge=[]):
		self.decaymode = decaymode
		self.sel = sel
		self.ntracks = -1
		self.charge_trk1 = -111  # dont make default -1 since that's a valid charge! :)
		self.charge_trk2 = -222  # dont make default -1 since that's a valid charge! :)
		# also, don't make the same since the equality is checked later
		self.two_plus = False
		self.two_minus = False


		if self.decaymode == "leptonic":
			self.ntracks = tree.ntrk

			if self.ntracks == 2: 
				self.charge_trk1 = tree.dv('trk_charge')[0]
				self.charge_trk2 = tree.dv('trk_charge')[1]	
			else: 
				if len(trk_charge) == 2: 
					self.charge_trk1 = trk_charge[0]
					self.charge_trk2 = trk_charge[1]

			if self.charge_trk1 == 1 and self.charge_trk2 == 1: 
					self.two_plus = True
			if self.charge_trk1 == -1 and self.charge_trk2 == -1: 
				self.two_minus = True


	def passes(self): 
		if self.sel == 'OS':
			if self.charge_trk1 != self.charge_trk2: 
				return True
				
		elif self.sel == 'SS':
			if self.charge_trk1 == self.charge_trk2: 
				return True
		else:
			return False



class DVtype():
	def __init__(self, tree, dv_type, decaymode="leptonic"):
		self.tree = tree
		self.decaymode = decaymode
		self.dv_type = dv_type
		self.lepton_charge = []

		if self.decaymode == "leptonic":
			self.ntracks = self.tree.ntrk
			self.nel = -1
			self.nmu = -1

			self.muons = helpers.Tracks(self.tree)
			self.muons.getMuons()
			self.nmu = len(self.muons.lepVec)

			self.electrons = helpers.Tracks(self.tree)
			self.electrons.getElectrons()
			self.nel = len(self.electrons.lepVec)

	def passes(self):
		combined = 0

		if self.dv_type == "emu": 
			if self.nel == 1 and self.nmu == 1: 
				mu1_type = self.tree['muon_type'][self.muons.lepIndex[0]]

				if mu1_type == combined:  # Only count combined muons 
					self.lepton_charge.append(self.electrons.lepCharge[0])
					self.lepton_charge.append(self.muons.lepCharge[0])
					return True
				else:
					return False
			else:
				return False

		elif self.dv_type == "mumu":
			if self.nmu == 2: 
				mu1_type = self.tree['muon_type'][self.muons.lepIndex[0]]
				mu2_type = self.tree['muon_type'][self.muons.lepIndex[1]]

				if mu1_type == combined and mu2_type == combined:  # Only count combined muons
					self.lepton_charge.append(self.muons.lepCharge[0])
					self.lepton_charge.append(self.muons.lepCharge[1])
					return True
				else:
					return False
			else:
				return False

		elif self.dv_type == "ee":
			if self.nel == 2: 
				return True
				lepton_charge.append(self.electrons.lepCharge[0])
				lepton_charge.append(self.electrons.lepCharge[1])
			else:
				return False

		elif self.dv_type == "mumu-notcomb":
			if self.nmu == 2:
				lepton_charge.append(self.muons.lepCharge[0])
				lepton_charge.append(self.muons.lepCharge[1])
				return True
			else: 
				return False

		elif self.dv_type == "1-lep":
			if self.nmu > 0 or self.nel> 0: 			
				return True
			else: 
				return False
		elif self.dv_type == "2-lep":
			if self.nmu == 2 or (self.nmu == 1 and self.nel == 1) or self.nel ==2:
				return True
			else:
				return False


class Trackqual():
	def __init__(self, tree, decaymode="leptonic", quality="2-tight"):

		self.decaymode = decaymode
		self.quality = quality 
		self.DV_2tight = False
		self.DV_1tight = False
		self.DV_2medium = False
		self.DV_1medium = False
		self.DV_2loose = False
		self.DV_1loose = False
		self.DV_med_vl = False
		self.DV_med_vvl = False
		self.DV_med_vvlSi = False
		self.DV_loose_vl = False
		self.DV_loose_vvl = False
		self.DV_loose_vvlSi = False
		self.DV_2vl = False
		self.DV_2vvl = False
		self.DV_2vvlSi = False



		if self.decaymode == "leptonic": 
			muons = helpers.Tracks(tree)
			muons.getMuons()

			electrons = helpers.Tracks(tree)
			electrons.getElectrons()

			self.nmu_tight = 0
			self.nmu_medium = 0
			self.nmu_loose = 0
			self.nel_tight = 0
			self.nel_medium = 0
			self.nel_loose = 0
			self.nel_vl = 0
			self.nel_vvl = 0
			self.nel_vvlSi = 0

			self.ndvmu = len(muons.lepVec)
			self.ndvel = len(electrons.lepVec)
		
			for imu in range(self.ndvmu):
				muindex = muons.lepIndex[imu]
				muisTight = tree['muon_isTight'][muindex]
				muisMedium = tree['muon_isMedium'][muindex]
				muisLoose = tree['muon_isLoose'][muindex]
				# check if Tight == 1 to incase safeFill was used and isTight == -1 (which is also not Tight!) -DT
				if muisTight == 1:
					self.nmu_tight = self.nmu_tight + 1
				if muisMedium == 1:
					self.nmu_medium = self.nmu_medium + 1
				if muisLoose == 1:
					self.nmu_loose = self.nmu_loose + 1

			for iel in range(self.ndvel):
				elindex = electrons.lepIndex[iel]
				elisTight = tree['el_LHTight'][elindex]
				elisMedium = tree['el_LHMedium'][elindex]
				elisLoose = tree['el_LHLoose'][elindex]
				elisvl = tree['el_isLHVeryLoose'][elindex]
				elisvvl = tree['el_isLHVeryLoose_mod1'][elindex]
				elisvvlSi = tree['el_isLHVeryLoose_modSi'][elindex]
				if elisTight == 1:
					self.nel_tight = self.nel_tight + 1
				if elisMedium == 1:
					self.nel_medium = self.nel_medium + 1
				if elisLoose == 1:
					self.nel_loose = self.nel_loose + 1
				if elisvl ==1: 
					self.nel_vl = self.nel_vl + 1
				if elisvvl == 1: 
					self.nel_vvl = self.nel_vvl + 1
				if elisvvlSi == 1: 
					self.nel_vvlSi = self.nel_vvlSi + 1

			self.DV_2tight = self.nmu_tight == 2 or self.nel_tight == 2 or (self.nmu_tight == 1 and self.nel_tight == 1) 
			self.DV_2medium = self.nmu_medium == 2 or self.nel_medium == 2 or (self.nmu_medium == 1 and self.nel_medium == 1)
			self.DV_2loose = self.nmu_loose == 2 or self.nel_loose == 2 or (self.nmu_loose == 1 and self.nel_loose == 1)
			self.DV_1tight = self.nmu_tight > 0 or self.nel_tight > 0
			self.DV_1medium = self.nmu_medium > 0 or self.nel_medium > 0
			self.DV_1loose = self.nmu_loose > 0 or self.nel_loose > 0
			self.DV_2vvlSi =  self.nel_vvlSi  == 2
			self.DV_2vvl =  self.nel_vvl  == 2
			self.DV_2vl =  self.nel_vl  == 2
			self.DV_loose_vl = self.nmu_loose == 1 and self.nel_vl == 1
			self.DV_loose_vvl = self.nmu_loose == 1 and self.nel_vvl == 1
			self.DV_loose_vvlSi = self.nmu_loose == 1 and self.nel_vvlSi == 1
			self.DV_med_vl = self.nmu_medium == 1 and self.nel_vl == 1
			self.DV_med_vvl = self.nmu_medium == 1 and self.nel_vvl == 1
			self.DV_med_vvlSi = self.nmu_medium == 1 and self.nel_vvlSi == 1

	def passes(self):
		if self.quality == "2-tight":
			return self.DV_2tight
			
		if self.quality == "2-medium":
			return self.DV_2medium
		
		if self.quality == "2-loose":
			return self.DV_2loose 

		if self.quality == "1-tight":
			return self.DV_1tight

		if self.quality == "1-medium":
			return self.DV_1medium

		if self.quality == "1-loose":
			return self.DV_1loose

class Cosmicveto():
	def __init__(self, tree, decaymode="leptonic", cosmicvetocut=0.05):

		self.decaymode = decaymode
		self.cosmicvetocut = cosmicvetocut

		self.separation = -1

		if self.decaymode == "leptonic":
			ntracks = tree.ntrk
			if ntracks == 2:
				sumeta = tree.dv('trk_eta_wrtSV')[0] + tree.dv('trk_eta_wrtSV')[1]
				dphi = abs(tree.dv('trk_phi_wrtSV')[0] - tree.dv('trk_phi_wrtSV')[1])

				self.separation = np.sqrt(sumeta ** 2 + (np.pi - dphi) ** 2)

	def passes(self):
		if self.separation > self.cosmicvetocut:
			return True
		else:
			return False


class Mlll():
	"""
	Trilepton mass calculation.
	The invariant mass of the prompt lepton and both diplaced leptons.
	"""

	def __init__(self, dv_type, plep, dMu, dEl, decaymode="leptonic", minmlll=50, maxmlll=84):
		self.decaymode = decaymode
		self.dv_type = dv_type
		self.plep = plep
		self.dMu = dMu
		self.dEl = dEl
		self.minmlll = minmlll
		self.maxmlll = maxmlll

		self.mlll = -1
		self.mtrans = -1
		self.plll = ROOT.TLorentzVector(0, 0, 0, 0)

		if self.decaymode == "leptonic":

			if self.dv_type == "emu":
				self.plll = self.plep + self.dEl[0] + self.dMu[0]

			if self.dv_type == "mumu":
				self.plll = self.plep + self.dMu[0] + self.dMu[1]

			if self.dv_type == "ee":
				self.plll = self.plep + self.dEl[0] + self.dEl[1]

			self.mlll = self.plll.M()
			self.mtrans = self.plll.Perp()

	def passes(self):

		if self.mlll > self.minmlll and self.mlll < self.maxmlll:
			return True
		else:
			return False


class Mltt():

	def __init__(self, plep, trks, decaymode="leptonic", minmltt=50, maxmltt=84):
		self.decaymode = decaymode
		self.plep = plep
		self.trks = trks
		self.minmltt = minmltt
		self.maxmltt = maxmltt

		self.mltt = -1
		self.plll = ROOT.TLorentzVector(0, 0, 0, 0)

		if self.decaymode == "leptonic":
			self.plll = self.plep + self.trks[0] + self.trks[1]
			self.mltt = self.plll.M()

	def passes(self):
		if (self.mltt > self.minmltt and self.mltt < self.maxmltt):
			return True
		else:
			return False


class Mtrans():
	def __init__(self, plep, trks, decaymode="leptonic", minmtrans=50, maxmtrans=84):
		self.decaymode = decaymode
		self.plep = plep
		self.trks = trks
		self.minmtrans = minmtrans
		self.maxmtrans = maxmtrans

		self.mtrans = -1
		self.mvis = -1
		self.plll = ROOT.TLorentzVector(0, 0, 0, 0)

		if self.decaymode == "leptonic":
			self.plll = self.plep + self.trks[0] + self.trks[1]
			self.mvis = self.plll.M()
			self.mtrans = self.plll.Perp()

	def passes(self):

		if self.mtrans > self.minmtrans and self.mtrans < self.maxmtrans:
			return True
		else:
			return False


class DVmass():
	def __init__(self, tree, decaymode="leptonic", dvmasscut=4):

		self.decaymode = decaymode
		self.dvmasscut = dvmasscut
		self.dvmass = tree.dv('mass')

	def passes(self):
		if self.dvmass > self.dvmasscut:
			return True
		else:
			return False


class Mhnl_old():
	def __init__(self, tree, plep, trks, hnlmasscut=4):

		self.plep = plep
		self.trks = trks
		self.hnlmasscut = hnlmasscut

		self.mvis = -1 
		self.mhnl = -1
		self.mhnl2 = -1
		self.hnlpt = -1
		self.hnleta = -99
		self.hnlphi = -99
		self.mtrans_rot = -1 
		self.neg_mhnl12 = -1
		self.neg_mhnl13 = -1
		self.neg_mhnl23 = -1
		self.pos_mhnl12 = -1
		self.pos_mhnl13 = -1
		self.pos_mhnl23 = -1

		self.nu_vec = ROOT.TLorentzVector()
		self.nu_vec2 = ROOT.TLorentzVector()

		def rotate_vector(r,v):
			r_new = ROOT.TVector3(v)
			rotation_axis = ROOT.TVector3(-1*r.Y(),r.X(),0.0)
			rotation_angle = -1*r.Theta() 
			r_new.Rotate(rotation_angle,rotation_axis)
			
			if (r== v and r_new.X() > 0.001): 
				#if r=v then you should end up with a vector all in the z component
				logger.error("Roatating vectors did not work!! Check HNL mass calculation.")
			return r_new

		def unrotate_vector(r,v):
			r_new = v

			rotation_axis = ROOT.TVector3(-1*r.Y(),r.X(),0.0)
			rotation_angle = r.Theta() 
			r_new.Rotate(rotation_angle,rotation_axis)
		
			return r_new

		# primary vertex vector
		dv_vec = ROOT.TVector3(tree.dv('x'),
							   tree.dv('y'),
							   tree.dv('z'))
		# primary vertex position vector
		pv_vec = ROOT.TVector3(tree['vertex_x'],
							   tree['vertex_y'],
							   tree['vertex_z'])

		# vector defining direction hnl trajectory
		hnl_vec =  dv_vec- pv_vec
		
		lepp_vec = ROOT.TVector3(plep.Px(),plep.Py(),plep.Pz()) 
		trkp_vec = []
		ntrk = len(self.trks)
		for i in range(ntrk):
			trkp_vec.append( ROOT.TVector3(self.trks[i].Px(),self.trks[i].Py(),self.trks[i].Pz()) )

		# print "Original DV vector: (", hnl_vec.X(),",", hnl_vec.Y(),",",hnl_vec.Z() , ")"

		#rotate coordinate system so hnl vector = z-axis
		lepp_vec_rot = rotate_vector(hnl_vec,lepp_vec)
		hnl_vec_rot = rotate_vector(hnl_vec,hnl_vec)
		trkp_vec_rot = []
		for i in range(ntrk):
			trkp_vec_rot.append(rotate_vector(hnl_vec,trkp_vec[i]))

		# check hnl_vec_rot is all in z component (x & y maybe be very small (~10^-15) due to precision when rotating vector)
		# print "New vector: (", hnl_vec_rot.X(),",", hnl_vec_rot.Y(),",",hnl_vec_rot.Z() , ")"

		# neutrino x-y momentum in rotated plane equals -1* sum of visible x-y momentum 
		pnu_rot = ROOT.TLorentzVector() 
		pnu_rot_x = -1*trkp_vec_rot[0].Px()-trkp_vec_rot[1].Px()
		pnu_rot_y = -1*trkp_vec_rot[0].Py()-trkp_vec_rot[1].Py()

		pnu_rot.SetPx(pnu_rot_x)
		pnu_rot.SetPy(pnu_rot_y)

		#make 4 vectors for visible particles in the rotated coordinate system (assume m = 0)
		plep_rot = ROOT.TLorentzVector(lepp_vec_rot,lepp_vec_rot.Mag())
		ptk1_rot = ROOT.TLorentzVector(trkp_vec_rot[0],trkp_vec_rot[0].Mag())
		ptk2_rot = ROOT.TLorentzVector(trkp_vec_rot[1],trkp_vec_rot[1].Mag())

		# make 4-vector for the visible particles in rotated coordinate system
		pvis_rot = plep_rot + ptk1_rot + ptk2_rot
		self.mvis = pvis_rot.M() # visible mass of the system
		self.mtrans_rot = pvis_rot.Perp() 

		# solve pw = pvis + pnu to find pnu_z
		m_w = 80.379 # mass of W boson in GeV
		K = (m_w**2 - self.mvis**2)/2 + pvis_rot.Px()*pnu_rot.Px() + pvis_rot.Py()*pnu_rot.Py()
		A = pvis_rot.Pz()**2 - pvis_rot.E()**2
		B = 2*K*pvis_rot.Pz()
		C = K**2 - pvis_rot.E()**2*(pnu_rot.Px()**2+ pnu_rot.Py()**2)

		# two solutions to the quadratic equation
		if (B**2 - 4*A*C) < 0: 
			#pw != pvis+ pnu so we can't solve quadratic
			return
		else: 
			pnu_rot_z1= (-B + np.sqrt(B**2 - 4*A*C) )/(2*A) 
			pnu_rot_z2= (-B - np.sqrt(B**2 - 4*A*C) )/(2*A)

			# two solutions for neutrino momentum vector
			pnu_rot1_vec = ROOT.TVector3(pnu_rot_x,pnu_rot_y,pnu_rot_z1) 
			pnu_rot2_vec = ROOT.TVector3(pnu_rot_x,pnu_rot_y,pnu_rot_z2) 

			# make a 4-vector for the neutrino
			pnu_rot1 = ROOT.TLorentzVector(pnu_rot1_vec,pnu_rot1_vec.Mag())
			pnu_rot2 = ROOT.TLorentzVector(pnu_rot2_vec,pnu_rot2_vec.Mag())

			# 2 solutions for HNL 4-vector
			pHNL_1 = ptk1_rot + ptk2_rot + pnu_rot1 # postive root
			pHNL_2 = ptk1_rot + ptk2_rot + pnu_rot2 # negative root
						#1        #2          #3
			
			neg_mhnl12_vec = ptk1_rot + ptk2_rot
			neg_mhnl13_vec = ptk1_rot + pnu_rot2
			neg_mhnl23_vec = ptk2_rot + pnu_rot2


			pos_mhnl12_vec = ptk1_rot + ptk2_rot
			pos_mhnl13_vec = ptk1_rot + pnu_rot1
			pos_mhnl23_vec = ptk2_rot + pnu_rot1



			# truth studies show pHNL_2 is the solution that gets us the HNL mass
			self.mhnl2 = pHNL_1.M()
			self.mhnl = pHNL_2.M()

			pHNL_2_lab = unrotate_vector(hnl_vec,pHNL_2)
			hnl_vec_unrot = unrotate_vector(hnl_vec,hnl_vec_rot)

			nu_vec_unrot2 = unrotate_vector(hnl_vec,pnu_rot1)

			nu_vec_unrot = unrotate_vector(hnl_vec,pnu_rot2)

			self.nu_vec = nu_vec_unrot
			self.nu_vec2 = nu_vec_unrot2

			# check that you rotated back to the original reference frame
			#print "Unrotated DV vector: (", hnl_vec_unrot.X(),",", hnl_vec_unrot.Y(),",",hnl_vec_unrot.Z() , ")"

			self.hnlpt = pHNL_2_lab.Pt()
			self.hnleta = pHNL_2_lab.Eta()
			self.hnlphi = pHNL_2_lab.Phi()

			# print self.mhnl, " ", self.hnlpt, " ", self.hnleta, " ",self.hnlphi

			self.neg_mhnl12 = neg_mhnl12_vec.M()
			self.neg_mhnl13 = neg_mhnl13_vec.M()
			self.neg_mhnl23 = neg_mhnl23_vec.M()

			self.pos_mhnl12 = pos_mhnl12_vec.M()
			self.pos_mhnl13 = pos_mhnl13_vec.M()
			self.pos_mhnl23 = pos_mhnl23_vec.M()

	def passes(self):
		
		if (self.mhnl > self.hnlmasscut ):
			return True
		else: 
			return False

class Mhnl():
	def __init__(self, tree, dv_type, plep, dMu, dEl, MW = 80.379,fixWMass=False, hnlmasscut = 4,use_truth=False,truth_pv=ROOT.TVector3(), truth_dv=ROOT.TVector3() ):
		# Global W pole mass
		MW = 80.379
		MW2 = MW**2
		WGamma = 2.085

		self.mhnl = -1
		self.alt_mhnl = -1
		self.hnlpt = -1
		self.hnleta = -99
		self.hnlphi = -99
		self.mlll = -1

		dtrks = []
		if dv_type == "emu":
			dtrks.append(dEl[0])
			dtrks.append(dMu[0])

		if dv_type == "mumu":
			dtrks.append(dMu[0])
			dtrks.append(dMu[1])

		if dv_type == "ee":
			dtrks.append(dEl[0])
			dtrks.append(dEl[1])

		# Get 3 vectors
		if not use_truth: 	
			pv = ROOT.TVector3( tree.dv('x'), tree.dv('y'),  tree.dv('z') )
			dv = ROOT.TVector3( tree['vertex_x'], tree['vertex_y'],  tree['vertex_z'])
		else: 
			pv = truth_pv
			dv = truth_dv

		p0 = ROOT.TVector3( plep.Px(), plep.Py(), plep.Pz() )

		d0 = ROOT.TVector3( dtrks[0].Px(), dtrks[0].Py(), dtrks[0].Pz() )
		d1 = ROOT.TVector3( dtrks[1].Px(), dtrks[1].Py(), dtrks[1].Pz() )

		def findMass(pv, dv, p0, d0, d1, MW2, fixWMass):
			# Choose z direction to be along decay
			decayV = dv-pv
			z = decayV*(1.0/decayV.Mag())
			
			# Visible (2 decay lepton) system
			dvis = d0+d1
			mvis2 = 2*(d0.Mag()*d1.Mag() - d0.Dot(d1))
			
			# Define plane perpendicular to direction of decay
			x = dvis.Cross(z)
			qperp = x.Mag()
			qperp2 = qperp*qperp
			
			# x and y unit vectors
			x = x*(1.0/qperp)
			y = z.Cross(x)

			# Visible momentum in the (x,y,z) coordinates
			qv = ROOT.TVector3( 0., qperp, dvis.Dot(z) )
			Ev = np.sqrt(mvis2 + qv.Dot(qv))
			qperp3 = ROOT.TVector3( 0., qperp, 0 )  #not needed, just for checking
			# Prompt lepton in new (x,y,z) coordinates
			pp = ROOT.TVector3(  p0.Dot(x), p0.Dot(y), p0.Dot(z) ) 
			Ep = np.sqrt(pp.Dot(pp))

			# Terms from conservation of 4 momentum that involve various powers of the neutrino z momentum
		#    A = (MW2 - mvis2)/2. - Ep*Ev - qperp2 + pp[2]*qv[2] 
			B = (pp[2] + qv[2])
			E = Ep + Ev
		#    A = A/E
			B = B/E
			
			# Min W mass possible
			alpha = qperp * B / np.sqrt(1-B**2)
			qn1 = ROOT.TVector3(  0, -qperp, alpha )
			En1 = qn1.Mag()
			qtot1 = pp + qv + qn1
			Etot1 = Ep + Ev + En1
			mWMin2 = Etot1**2 - qtot1.Dot(qtot1)
			mWMin = np.sqrt(mWMin2)

			cdVal = 0.5 + np.arctan((mWMin2-MW2)/MW/WGamma)/np.pi
			cdMed = (1+cdVal)/2
		#   invert to get mass of median allowed range
			mMed2 = MW2 + MW*WGamma*np.tan(np.pi*(cdMed-0.5))
		#    mMed = MW + WGamma*np.tan(np.pi*(cdMed-0.5))

			mMed = np.sqrt(mMed2)
			
		#    MW2fit = max(MW2, mWMin2+1)
			if fixWMass:
				MW2fit = MW2
			else:
				MW2fit = mMed2
			
			#    MW2fit = MW2
			A = (MW2fit - mvis2)/2. - Ep*Ev - qperp2 + pp[2]*qv[2] 
			A = A/E

			# Coefficients of the quadratic to solve
			b = 2*A*B/(B*B-1)
			c = (A*A - qperp2)/(B*B-1)
			
			arg = b*b - 4*c
			# protect against imaginary solutions (investigate any of these)
			noSol = 0
			if arg>0:
				rad = np.sqrt(arg)
			else:
				rad = -1000
				noSol = 1
				
			# These are the possible z momenta for the neutrino
			sol1 = (-b + rad)/2
			sol2 = (-b - rad)/2
			
			# Make vectors of the two z-momentum solutions for the neutrino
			qn1 = ROOT.TVector3(0, -qperp, sol1 )
			En1 = qn1.Mag()
			qn2 = ROOT.TVector3(0, -qperp, sol2 )
			En2 = qn2.Mag()

			# Total momentum of the decaying system
			qtot1 = qv + qn1
			qtot2 = qv + qn2
			
			# neutrino momentum in original coordinates 
			dn1 = y*(-1*qperp) + z*sol1
			dn2 = y*(-1*qperp) + z*sol2
			
			# mass of nuetrino+d0 lepton
			qn0 = dn1+d0
			mn01 = np.sqrt( (d0.Mag()+En1)**2 - qn0.Dot(qn0))
			qn0 = dn2+d0
			mn02 = np.sqrt( (d0.Mag()+En2)**2 - qn0.Dot(qn0))

			# And the mass of the N
			mN21 = (Ev+En1)**2 - qtot1.Dot(qtot1)
			mN22 = (Ev+En2)**2 - qtot2.Dot(qtot2)
			
			# Calculate the dependence of the mass solutions on the W mass
			dmN1dMW = (Ev*sol1/np.sqrt(qperp2+sol1**2) - qv[2])/np.sqrt(mN21) * (A + B*sol1) * np.sqrt(MW2) / \
			(((B**2-1)*sol1 + A*B)*E)
			
			dmN2dMW = (Ev*sol2/np.sqrt(qperp2+sol2**2) - qv[2])/np.sqrt(mN22) * (A + B*sol2) * np.sqrt(MW2) / \
			(((B**2-1)*sol2 + A*B)*E)


			# make 4-vectors in original coordinates 
			pnu2 = ROOT.TLorentzVector(dn2,dn2.Mag())
			pnu1 = ROOT.TLorentzVector(dn1,dn1.Mag())
			pion_mass = 0.139 # GeV
			plep0 = ROOT.TLorentzVector()
			ptrk0 = ROOT.TLorentzVector()
			ptrk1 = ROOT.TLorentzVector()
			# using tracking assumption of pion mass for consistency with trk quantities  
			plep0.SetPxPyPzE( p0.X(), p0.Y() , p0.Z(), np.sqrt(p0.Mag()**2 +pion_mass**2) )
			ptrk0.SetPxPyPzE( d0.X(), d0.Y() , d0.Z(), np.sqrt(d0.Mag()**2 +pion_mass**2) )
			ptrk1.SetPxPyPzE( d1.X(), d1.Y() , d1.Z(), np.sqrt(d1.Mag()**2 +pion_mass**2) )

			# assume massless tracks 
			# plep0 = ROOT.TLorentzVector( p0 ,p0.Mag())
			# ptrk0 = ROOT.TLorentzVector( d0 ,d0.Mag())
			# ptrk0 = ROOT.TLorentzVector( d1,d1.Mag())

			pHNL1 = pnu1 + ptrk0 + ptrk1
			pHNL2 = pnu2 + ptrk0 + ptrk1

			plll = plep0 + ptrk0 + ptrk1
		   
			# set the attributes of the class
			self.mhnl =pHNL2.M()
			self.hnlpt =pHNL2.Pt()
			self.hnleta =pHNL2.Eta()
			self.hnlphi =pHNL2.Phi()
			self.mlll = plll.M()
		  
			self.alt_mhnl = pHNL1.M()

		findMass(pv, dv, p0, d0, d1, MW2,fixWMass)

	def passes(self):
		return self.mhnl > self.hnlmasscut
		
class PV():
	def __init__(self, tree):

		self.pv_x = tree['vertex_x']
		self.pv_y = tree['vertex_y']
		self.pv_z = tree['vertex_z']

	def passes(self):
		
		if (self.pv_x != -999.0 and self.pv_y != -999.0 and self.pv_z != -999.0 ): 
			return True
		else: 
			return False # no primary vertex in the event

class MCEventType:
	def __init__(self, tree, wrong_lep_order = True):
		self.weight = 1 # if not data weight is default, event is neither LNC or LNV
		self.M2_spin_corr = -1
		self.M2_nocorr = -1
		self.isLNC = False
		self.isLNV = False
		# Matrix elements for the trilepton process, when only the charged-current contribution is present.

		# Input:
		# - All masses are in GeV, and all Mandelstam variables in GeV^2.
		# - pW2 = p(W)^2 = mW^2
		# - The Mandelstam variables are defined as follows:
		#       s13 = (p(l1)+p(l3))^2
		#       s24 = (p(l2)+p(nu))^2
		#   with:
		#     - l1 the charged lepton produced along with the HNL;
		#     - l2 the charged lepton produced in the HNL decay, on the same fermion line as the HNL;
		#     - l3 the charged lepton produced in the HNL decay, on the same fermion line as the light neutrino;
		#     - nu the light neutrino produced in the HNL decay.
		#   There are two additional, independent Mandelstam variables which do not enter the matrix elements.

		# Assumptions:
		# - Light lepton (e, mu) masses are neglected.
		# - Some off-shell effects coming from the finite width of the W are neglected.
		# - Matrix elements are given up to a dimensionful constant, but should be consistent among themselves.
		#   In particular: M2_nocorr = M2_LNC + M2_LNV.


		MW = 80.379 # Change the W mass to match your particle data. This is the latest PDG value.

		# Lepton number conserving charged-current trilepton process.
		def M2_LNC(MN, pW2, s13, s24):
		    return s24 * ( 2*MW**2*pW2*(MN**2-s24) + MN**4*(s13-2*MW**2) + 2*(s24-s13)*MN**2*MW**2 ) / ( 6*MW**6 )

		# Lepton number violating charged-current trilepton process.
		def M2_LNV(MN, pW2, s13, s24):
		    return - s24 * MN**2 * ( pW2*(s24-MN**2) + (s13-s24)*MN**2 + MN**4 - 2*s13*MW**2 ) / ( 6*MW**6 )

		# Charged-current trilepton process, ignoring spin correlations between the HNL production and its decay.
		def M2_nocorr(MN, pW2, s24):
		    return ( s24*(MN**2-s24)*(MN**2+2*MW**2)*(pW2-MN**2) ) / ( 6*MW**6 )

		if not tree.is_data:
			truth_info = helpers.Truth()
			truth_info.getTruthParticles(tree)
			pW2 = truth_info.W_vec.Mag2()
			MN = truth_info.HNL_vec.M()
			# print "HNL mass ", MN
			charge_1 = truth_info.plep_charge # charge of prompt lepton
			self.p_1 = truth_info.plep_vec # prompt lepton 
			if wrong_lep_order: 
				self.p_2 = truth_info.dLepVec[0]
				self.p_3 = truth_info.dLepVec[1]
				self.p_4 = truth_info.dLepVec[2]
			else: 
				self.p_2 = truth_info.dLepVec[2]
				self.p_3 = truth_info.dLepVec[0]
				self.p_4 = truth_info.dLepVec[1]

			if charge_1 != truth_info.dLepCharge[0]: 
				self.isLNC = True
			else: 
				self.isLNV = True
				
			if self.isLNC == self.isLNV: 
				logger.error("MCEventType selection found that this event is both LNC and LNV. Check this event!")


			p12 = self.p_1 + self.p_2
			p13 = self.p_1 + self.p_3
			p14 = self.p_1 + self.p_4

			p23 = self.p_2 + self.p_3
			p24 = self.p_2 + self.p_4
			p34 = self.p_3 + self.p_4

			self.s12 = p12.Mag2()
			self.s13 = p13.Mag2()
			self.s14 = p14.Mag2()
			self.s23 = p23.Mag2()
			self.s24 = p24.Mag2()
			self.s34 = p34.Mag2()
			
			# calculate the correct matrix that takes into account LNC or LNV decay
			if self.isLNC: 
				self.M2_spin_corr = M2_LNC(MN=MN, pW2=pW2, s13=self.s13, s24=self.s24)
			if self.isLNV: 
				self.M2_spin_corr = M2_LNV(MN=MN, pW2=pW2, s13=self.s13, s24=self.s24)

			if wrong_lep_order: 
			# N.B Official samples have wrong lepton ordering where lepton 2 and lepton 4 are swapped i.e instead of 1234 we have 1423.
			# For official samples, swap s24 -> s34 
				self.M2_nocorr = M2_nocorr(MN=MN, pW2=pW2,s24= self.s34) # wrong matrix used when generating pythia samples, includes lepton permutation
			else: 
				self.M2_nocorr = M2_nocorr(MN=MN, pW2=pW2, s24= self.s24) # wrong matrix used when generating pythia samples, NO lepton permutation

			self.weight = 2*self.M2_spin_corr/self.M2_nocorr  # factor of 2 here is becuase M2_nocorr as calculated includes LNC + LNV decays

		


class SumTrack:
	def __init__(self, tree):

		self.sum_track_pt = 0
		n_tracks = tree.ntrk
		self.sum_track_pt_wrt_pv = 0
		self.sum_track_charge = 0

		for k in range(n_tracks):
			self.sum_track_pt += tree.dv('trk_pt_wrtSV')[k]
			self.sum_track_pt_wrt_pv += tree.dv('trk_pt')[k]
			self.sum_track_charge += tree.dv('trk_charge')[k]





