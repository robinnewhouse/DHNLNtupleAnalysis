import ROOT
import numpy as np
import helpers
import logging

# make a global logger variable for all selection classes
logger = helpers.getLogger('dHNLAnalysis.selections', level=logging.INFO)


class Trigger:
	def __init__(self, tree, trigger, invert=False):
		self.event_triggers = tree['passedTriggers']
		self.invert = invert
		if trigger == "muononly":
			if tree.mc_campaign == "mc16a":
				self.allowed_trigger_list = helpers.SingleMuonTriggerlist_2015_2016
			if tree.mc_campaign == "mc16d":
				self.allowed_trigger_list = helpers.SingleMuonTriggerlist_2017
			if tree.mc_campaign == "mc16e":
				self.allowed_trigger_list = helpers.SingleMuonTriggerlist_2018
			else:
				self.allowed_trigger_list = helpers.SingleMuonTriggerlist
		elif trigger == "electrononly":
			if tree.mc_campaign == "mc16a":
				self.allowed_trigger_list = helpers.SingleElectronTriggerlist_2015_2016
			if tree.mc_campaign == "mc16d":
				self.allowed_trigger_list = helpers.SingleElectronTriggerlist_2017
			if tree.mc_campaign == "mc16e":
				self.allowed_trigger_list = helpers.SingleElectronTriggerlist_2018
			else:
				self.allowed_trigger_list = helpers.SingleElectronTriggerlist
		elif trigger == "alltriggers":
			if tree.mc_campaign == "mc16a":
				self.allowed_trigger_list = helpers.SingleMuonTriggerlist_2015_2016 + helpers.SingleElectronTriggerlist_2015_2016
			if tree.mc_campaign == "mc16d":
				self.allowed_trigger_list = helpers.SingleMuonTriggerlist_2017 + helpers.SingleElectronTriggerlist_2017
			if tree.mc_campaign == "mc16e":
				self.allowed_trigger_list = helpers.SingleMuonTriggerlist_2018 + helpers.SingleElectronTriggerlist_2018
			else:
				self.allowed_trigger_list = helpers.SingleMuonTriggerlist + helpers.SingleElectronTriggerlist

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


class Filter:
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


class InvertedPromptLepton:
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


class PromptTrack:
	def __init__(self, tree, lepton="any", quality="tight", min_dR=0.05):
		self.trkVec = ROOT.TLorentzVector(0, 0, 0, 0)
		self.trkd0 = -2000
		self.trkz0 = -2000
		self.nPtrk = 0
		self.found_trk = False
		self.trk_Index = -1

		ntrk = len(tree['tracks_pt'])
		# variable for the highest pt lepton				
		self.highestpt_trk = ROOT.TLorentzVector(0, 0, 0, 0)
		self.highestpt_trk_d0 = -2000
		self.highestpt_trk_z0 = -2000
		self.highestpt_trk_Index = -1

		for itrk in range(ntrk):
			overlap = False
			trkVec_i = ROOT.TLorentzVector()
			pt = tree['tracks_pt'][itrk]
			eta = tree['tracks_eta'][itrk]
			phi = tree['tracks_phi'][itrk]
			mass = 0.139  # pion mass in GeV
			trkVec_i.SetPtEtaPhiM(pt, eta, phi, mass)

			trkd0 = tree['tracks_d0'][itrk]
			trkz0 = tree['tracks_z0'][itrk]
			trktheta = tree['tracks_theta'][itrk]
			trkz0sintheta = trkz0 * np.sin(trktheta)
			# print trktheta,np.sin(trktheta)
			# print trkz0sintheta
			if abs(trkd0) < 3 and abs(trkz0sintheta) < 0.5:
				# Check the overlap between the prompt lepton and every displaced vertex track
				self.found_trk = True
				for idv in range(tree.ndv):
					prefix = tree.dv_prefix + '_'
					ntrks = tree.get_at(prefix + 'ntrk', tree.ievt, idv)
					for itrk in range(ntrks):
						# Currently the only live example of get_at() which gives full control over the tree access.
						pt = tree.get_at(prefix + 'trk_pt', tree.ievt, idv, itrk)
						eta = tree.get_at(prefix + 'trk_eta', tree.ievt, idv, itrk)
						phi = tree.get_at(prefix + 'trk_phi', tree.ievt, idv, itrk)
						M = tree.get_at(prefix + 'trk_M', tree.ievt, idv, itrk)
						track_vector = ROOT.TLorentzVector()
						track_vector.SetPtEtaPhiM(pt, eta, phi, M)

						dR = track_vector.DeltaR(trkVec_i)
						if dR < min_dR:  # set overlap to true if muon overlaps with displaced track
							overlap = True

				# if lepton doesnt overlap with and DV tracks
				if not overlap:
					self.nPtrk = self.nPtrk + 1
					# if pt is larger then the previous prompt lepton found
					if trkVec_i.Pt() > self.highestpt_trk.Pt():
						self.highestpt_trk = trkVec_i  # get highest pt prompt lepton!
						self.highestpt_trk_d0 = trkd0
						self.highestpt_trk_z0 = trkz0
						self.highestpt_trk_z0sintheta = trkz0sintheta
						self.highestpt_trk_Index = itrk

	def passes(self):
		# check if you found a prompt lepton
		if self.nPtrk > 0 and self.highestpt_trk.Pt() > 0:
			self.trkVec = self.highestpt_trk
			self.trkd0 = self.highestpt_trk_d0
			self.trkz0 = self.highestpt_trk_z0
			self.trk_Index = self.highestpt_trk_Index
			return True
		else:
			return False


class PromptLepton:
	def __init__(self, tree, lepton="any", quality="tight", min_dR=0.05):
		self.plepVec = ROOT.TLorentzVector(0, 0, 0, 0)
		self.plepcharge = 0
		self.plep_isTight = False
		self.plepd0 = -2000
		self.plepz0 = -2000
		self.nPlep = 0
		self.found_plep = False
		self.plep_Index = -1

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
		self.highestpt_lep_Index = -1
		self.highestpt_lep_isTight = False

		for ilep in range(nleps):
			overlap = False
			plepVec_i = ROOT.TLorentzVector()

			if lepton == "muon":
				min_pt = 3  # GeV
				pt = tree['muon_pt'][ilep]
				eta = tree['muon_eta'][ilep]
				phi = tree['muon_phi'][ilep]
				mass = tree['muon_m'][ilep]
				charge = tree['muon_charge'][ilep]
				plepVec_i.SetPtEtaPhiM(pt, eta, phi, mass)

				lepd0 = tree['muon_trkd0'][ilep]
				lepz0 = tree['muon_trkz0'][ilep]
				lepz0sintheta = tree['muon_trkz0sintheta'][ilep]
				muon_type = tree['muon_type'][ilep]
				# skip any muons that are not combined muons!
				if muon_type != 0: continue

			if lepton == "electron":
				min_pt = 4.5  # GeV
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

			if lepton == "electron": plep_isTight = tree["el_LHTight"][ilep] > 0
			if lepton == "muon": plep_isTight = tree["muon_isTight"][ilep] > 0

			# apply filter cut seperately do not need to addtionally require our lepton passes the prompt filter cuts -DT
			# if passPfilter[ilep] and passes_lep_quality and abs(lepd0) < 3 and abs(lepz0sintheta) < 0.5:
			if passes_lep_quality and abs(lepd0) < 3 and abs(lepz0sintheta) < 0.5 and pt > min_pt:
				# Check the overlap between the prompt lepton and every displaced vertex track
				self.found_plep = True
				for idv in range(tree.ndv):
					prefix = tree.dv_prefix + '_'
					ntrks = tree.get_at(prefix + 'ntrk', tree.ievt, idv)
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
						self.highestpt_lep_Index = ilep
						self.highestpt_lep_isTight = plep_isTight

	def passes(self):
		# check if you found a prompt lepton
		if self.nPlep > 0 and self.highestpt_lep.Pt() > 0:
			self.plepVec = self.highestpt_lep
			self.plepd0 = self.highestpt_lep_d0
			self.plepz0 = self.highestpt_lep_z0
			self.plepcharge = self.highestpt_lep_charge
			self.plep_Index = self.highestpt_lep_Index
			self.plep_isTight = self.highestpt_lep_isTight
			return True
		else:
			return False


class Alpha:
	def __init__(self, tree, max_alpha=0.5):
		self.max_alpha = max_alpha
		# compute alpha (3D angle between DV 3-momentum and rDV)
		dv = ROOT.TVector3(tree.dv('x'), tree.dv('y'), tree.dv('z'))
		pv = ROOT.TVector3(tree['vertex_x'], tree['vertex_y'], tree['vertex_z'])
		decay_vector = dv - pv

		dv_4vec = ROOT.TLorentzVector()
		dv_4vec.SetPtEtaPhiM(tree.dv('pt'), tree.dv('eta'), tree.dv('phi'), tree.dv('mass'))
		dv_mom_vec = ROOT.TVector3(dv_4vec.Px(), dv_4vec.Py(), dv_4vec.Pz())
		self.alpha = decay_vector.Angle(dv_mom_vec)

	def passes(self):
		# pass if alpha is less than the sent cut
		return self.alpha < self.max_alpha


class DVRadius:
	def __init__(self, tree):
		self.rdv = -1
		if tree.ntrk > 0:
			dx = tree.dv('x')
			dy = tree.dv('y')
			self.rdv = np.sqrt(dx ** 2 + dy ** 2)

	def passes(self, rdv_min=4, rdv_max=300):
		if rdv_min < self.rdv < rdv_max:
			return True
		else:
			return False


class DVNTracks:
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


class ChargeDV:
	def __init__(self, tree, sel="OS", decaymode="leptonic", trk_charge=[]):
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


class DVtype:
	def __init__(self, tree, dv_type, decaymode="leptonic"):
		self.tree = tree
		self.decaymode = decaymode
		self.dv_type = dv_type
		self.lepton_charge = []
		self.dEl_Index = []
		self.dMu_Index = []

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
				if self.tree.fake_aod:  # skip muon type cut for now with fakeAOD
					# return True
					mu1_type = self.muons.muonType[0]
				else:
					mu1_type = self.tree['muon_type'][self.muons.lepIndex[0]]

				if mu1_type == combined:  # Only count combined muons 
					self.lepton_charge.append(self.electrons.lepCharge[0])
					self.lepton_charge.append(self.muons.lepCharge[0])
					self.dEl_Index.append(self.electrons.lepIndex[0])
					self.dMu_Index.append(self.muons.lepIndex[0])
					return True
				else:
					return False
			else:
				return False

		elif self.dv_type == "mumu":

			if self.nmu == 2:
				if self.tree.fake_aod:  # skip muon type cut for now with fakeAOD
					# return True
					mu1_type = self.muons.muonType[0]
					mu2_type = self.muons.muonType[1]
				else:
					mu1_type = self.tree['muon_type'][self.muons.lepIndex[0]]
					mu2_type = self.tree['muon_type'][self.muons.lepIndex[1]]

				if mu1_type == combined and mu2_type == combined:  # Only count combined muons
					self.lepton_charge.append(self.muons.lepCharge[0])
					self.lepton_charge.append(self.muons.lepCharge[1])
					self.dMu_Index.append(self.muons.lepIndex[0])
					self.dMu_Index.append(self.muons.lepIndex[1])
					return True
				else:
					return False
			else:
				return False

		elif self.dv_type == "ee":
			if self.nel == 2:
				self.lepton_charge.append(self.electrons.lepCharge[0])
				self.lepton_charge.append(self.electrons.lepCharge[1])
				self.dEl_Index.append(self.electrons.lepIndex[0])
				self.dEl_Index.append(self.electrons.lepIndex[1])
				return True
			else:
				return False

		elif self.dv_type == "mumu-notcomb":
			if self.nmu == 2:
				self.lepton_charge.append(self.muons.lepCharge[0])
				self.lepton_charge.append(self.muons.lepCharge[1])
				return True
			else:
				return False

		elif self.dv_type == "1-lep":
			if self.nmu > 0 or self.nel > 0:
				return True
			else:
				return False
		elif self.dv_type == "2-lep":
			if self.nmu == 2 or (self.nmu == 1 and self.nel == 1) or self.nel == 2:
				return True
			else:
				return False


class TrackQuality:
	def __init__(self, tree, decaymode="leptonic", quality="2-tight"):

		self.tree = tree
		self.decaymode = decaymode
		self.quality = quality
		# self.DV_2tight = False
		# self.DV_1tight = False
		# self.DV_2medium = False
		# self.DV_1medium = False
		# self.DV_2loose = False
		# self.DV_1loose = False
		# self.DV_med_vl = False
		# self.DV_med_vvl = False
		# self.DV_med_vvlSi = False
		# self.DV_loose_vl = False
		# self.DV_loose_vvl = False
		# self.DV_loose_vvlSi = False
		# self.DV_2vl = False
		# self.DV_2vvl = False
		# self.DV_2vvlSi = False

		if self.decaymode == "leptonic":
			muons = helpers.Tracks(self.tree)
			muons.getMuons()

			electrons = helpers.Tracks(self.tree)
			electrons.getElectrons()

			self.nmu_tight = 0
			self.nmu_medium = 0
			self.nmu_loose = 0
			self.nel_tight = 0
			self.nel_medium = 0
			self.nel_loose = 0
			self.nel_veryloose = 0
			self.nel_veryveryloose = 0
			self.nel_veryveryloosesi = 0

			self.ndvmu = len(muons.lepVec)
			self.ndvel = len(electrons.lepVec)

			for imu in range(self.ndvmu):
				if self.tree.fake_aod:  # get quality information from info decorated on tracks
					muisTight = muons.muon_isTight[imu]
					muisMedium = muons.muon_isMedium[imu]
					muisLoose = muons.muon_isLoose[imu]
				else:  # get quality information from matching muons
					muindex = muons.lepIndex[imu]
					muisTight = self.tree['muon_isTight'][muindex]
					muisMedium = self.tree['muon_isMedium'][muindex]
					muisLoose = self.tree['muon_isLoose'][muindex]
				# check if Tight == 1 to in case safeFill was used and isTight == -1 (which is also not Tight!) -DT
				if muisTight == 1:
					self.nmu_tight = self.nmu_tight + 1
				if muisMedium == 1:
					self.nmu_medium = self.nmu_medium + 1
				if muisLoose == 1:
					self.nmu_loose = self.nmu_loose + 1

			for iel in range(self.ndvel):
				if self.tree.fake_aod:  # get quality infomation from info decorated on tracks
					elisTight = electrons.el_isTight[iel]
					elisMedium = electrons.el_isMedium[iel]
					elisLoose = electrons.el_isLoose[iel]
					elisVeryLoose = electrons.el_isveryLoose[iel]
					elisVeryVeryLoose = electrons.el_isveryveryLoose[iel]
					elisVeryVeryLooseSi = electrons.el_isveryveryLooseSi[iel]
				else:  # get quality infomation from matching electrons
					elindex = electrons.lepIndex[iel]
					elisTight = self.tree['el_LHTight'][elindex]
					elisMedium = self.tree['el_LHMedium'][elindex]
					elisLoose = self.tree['el_LHLoose'][elindex]
					elisVeryLoose = self.tree['el_isLHVeryLoose'][elindex]
					elisVeryVeryLoose = self.tree['el_isLHVeryLoose_mod1'][elindex]
					elisVeryVeryLooseSi = self.tree['el_isLHVeryLoose_modSi'][elindex]

				if elisTight == 1:
					self.nel_tight = self.nel_tight + 1
				if elisMedium == 1:
					self.nel_medium = self.nel_medium + 1
				if elisLoose == 1:
					self.nel_loose = self.nel_loose + 1
				if elisVeryLoose == 1:
					self.nel_veryloose = self.nel_veryloose + 1
				if elisVeryVeryLoose == 1:
					self.nel_veryveryloose = self.nel_veryveryloose + 1
				if elisVeryVeryLooseSi == 1:
					self.nel_veryveryloosesi = self.nel_veryveryloosesi + 1

			self.DV_2tight = self.nmu_tight == 2 or self.nel_tight == 2 or (self.nmu_tight == 1 and self.nel_tight == 1)
			self.DV_2medium = self.nmu_medium == 2 or self.nel_medium == 2 or (self.nmu_medium == 1 and self.nel_medium == 1)
			self.DV_2loose = self.nmu_loose == 2 or self.nel_loose == 2 or (self.nmu_loose == 1 and self.nel_loose == 1)
			self.DV_1tight = self.nmu_tight > 0 or self.nel_tight > 0
			self.DV_1medium = self.nmu_medium > 0 or self.nel_medium > 0
			self.DV_1loose = self.nmu_loose > 0 or self.nel_loose > 0
			self.DV_tight_medium = (self.nmu_tight == 1 and self.nel_medium == 1) or (self.nel_tight == 1 and self.nmu_medium == 1) or (self.nmu_tight >= 1 and self.nmu_medium == 2)
			self.DV_tight_loose = (self.nmu_tight == 1 and self.nel_loose == 1) or (self.nel_tight == 1 and self.nmu_loose == 1) or (self.nmu_tight >= 1 and self.nmu_loose == 2)
			self.DV_medium_loose = (self.nmu_medium == 1 and self.nel_loose == 1) or (self.nel_medium == 1 and self.nmu_loose == 1) or (self.nmu_medium >= 1 and self.nmu_loose == 2)
			self.DV_tight_veryloose = self.nmu_tight == 1 and self.nel_veryloose == 1
			self.DV_medium_veryloose = self.nmu_medium == 1 and self.nel_veryloose == 1
			self.DV_loose_veryloose = self.nmu_loose == 1 and self.nel_veryloose == 1
			self.DV_tight_veryveryloose = self.nmu_tight == 1 and self.nel_veryveryloose == 1
			self.DV_medium_veryveryloose = self.nmu_medium == 1 and self.nel_veryveryloose == 1
			self.DV_loose_veryveryloose = self.nmu_loose == 1 and self.nel_veryveryloose == 1
			self.DV_2veryveryloose = self.nel_veryveryloose == 2
			self.DV_1veryveryloose = self.nel_veryveryloose == 1

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
		if self.quality == "tight-loose":
			return self.DV_tight_loose

		if self.quality == "tight-medium":
			return self.DV_tight_medium

		if self.quality == "medium-loose":
			return self.DV_medium_loose

		if self.quality == "tight-veryloose":
			return self.DV_tight_veryloose

		if self.quality == "medium-veryloose":
			return self.DV_medium_veryloose

		if self.quality == "loose-veryloose":
			return self.DV_loose_veryloose

		if self.quality == "tight-veryveryloose":
			return self.DV_tight_veryveryloose

		if self.quality == "medium-veryveryloose":
			return self.DV_medium_veryveryloose

		if self.quality == "2-veryveryloose":
			return self.DV_2veryveryloose

		if self.quality == "loose-veryveryloose":
			return self.DV_loose_veryveryloose

		if self.quality == "any-loose":
			return self.DV_any_loose

		if self.quality == "any-veryveryloose":
			return self.DV_loose_veryveryloose

		if self.quality == "2-any":
			return self.DV_2any


class CosmicVeto:
	def __init__(self, tree, decaymode="leptonic", cosmicvetocut=0.05):

		self.decaymode = decaymode
		self.cosmicvetocut = cosmicvetocut

		self.separation = -1

		if self.decaymode == "leptonic":
			ntracks = tree.ntrk
			if ntracks == 2:
				tracks = helpers.Tracks(tree)
				tracks.getTracks()
				trkVec = tracks.lepVec

				sumeta = tracks.lepVec[0].Eta() + tracks.lepVec[1].Eta()
				dphi = abs(tracks.lepVec[0].DeltaPhi(tracks.lepVec[1]))
				self.separation = np.sqrt(sumeta ** 2 + (np.pi - dphi) ** 2)

	def passes(self):
		return self.separation > self.cosmicvetocut


class Mlll:
	"""
	Trilepton mass calculation.
	The invariant mass of the prompt lepton and both diplaced leptons.
	"""

	def __init__(self, dv_type, plep, dMu, dEl, decaymode="leptonic", minmlll=50, maxmlll=84, invert=False):
		self.decaymode = decaymode
		self.dv_type = dv_type
		self.plep = plep
		self.dMu = dMu
		self.dEl = dEl
		self.minmlll = minmlll
		self.maxmlll = maxmlll
		self.invert = invert

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
		if not self.invert:
			return self.minmlll < self.mlll < self.maxmlll
		else:
			return self.mlll < self.minmlll or self.mlll > self.maxmlll


class MaterialVeto:
	"""
	Material Veto.
	This cut rejects any vertices whose (r, z, phi) position coincides with the location of known detector elements. 
	"""

	def __init__(self, tree):
		self.pass_mat = tree.dv('pass_mat')

	# print self.pass_mat

	def passes(self):
		return self.pass_mat


class Mltt:
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
		if self.minmltt < self.mltt < self.maxmltt:
			return True
		else:
			return False


class Mtrans:
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

		if self.minmtrans < self.mtrans < self.maxmtrans:
			return True
		else:
			return False


class DVmass:
	def __init__(self, tree, decaymode="leptonic", dvmasscut=5.5):

		self.decaymode = decaymode
		self.dvmasscut = dvmasscut
		self.dvmass = tree.dv('mass')

	def passes(self):
		if self.dvmass > self.dvmasscut:
			return True
		else:
			return False

class BHadronVeto:
	def __init__(self, tree, dv_type, dv_mass_cut=5.5):

		self.dv_mass_cut = dv_mass_cut
		self.dv_mass = tree.dv('mass')
		self.dv_type = dv_type
		self.rdv = -1

		if tree.ntrk > 0:
			dx = tree.dv('x')
			dy = tree.dv('y')
			self.rdv = np.sqrt(dx**2 + dy**2)

		self.pass_diagonal_cut = self.dv_mass > (-(7.0 / 150.0) * self.rdv + 7.0)
		# (((DV_mass >2 && DV_mass <5.5) && DV_mass > -(7/150)*DV_r + 7 ) || DV_mass > 5.5)


	def passes(self):
		if self.dv_type == "mumu" or self.dv_type == "ee":
			return self.dv_mass > self.dv_mass_cut
		elif self.dv_type == "emu":
			return ((2 < self.dv_mass < self.dv_mass_cut) and self.pass_diagonal_cut) or self.dv_mass > self.dv_mass_cut
		else:
			return False

class ZMassVeto:
	def __init__(self, tree, plep_vec, plep, plep_charge, dv_type):
		self.mll_dMu_plep_is_SS = None
		self.mll_dMu_plep_is_OS = None
		self.mll_dMu_plep = -1
		self.mll_dEl_plep_is_SS = None
		self.mll_dEl_plep_is_OS = None
		self.mll_dEl_plep = -1
		
		self.plep = plep
		self.dv_type = dv_type

		muons = helpers.Tracks(tree)
		muons.getMuons()
		muVec = muons.lepVec
		muVec_lepmatched = muons.lepmatched_lepVec
		mu_Index = muons.lepIndex

		electrons = helpers.Tracks(tree)
		electrons.getElectrons()
		elVec = electrons.lepVec
		elVec_lepmatched = electrons.lepmatched_lepVec
		el_Index = electrons.lepIndex

		# get invariant mass between prompt lepton + same flavour displaced lepton for Z->ll veto
		if plep == 'muon' and len(muVec) > 0:
			# only 1 muon in DV then choose this lepton for dlep in veto
			if len(muVec) == 1: mu_index = 0
			# choose highest pt muon in DV for dlep in veto
			elif len(muVec) == 2:
				if muVec_lepmatched[0].Pt() > muVec_lepmatched[1].Pt(): mu_index = 0
				else: mu_index = 1
			# get the 4-vector for highest pT same flavour lepton as plep + prompt lepton
			mu_plep_vec = muVec_lepmatched[mu_index] + plep_vec
			self.mll_dMu_plep_is_SS = int(plep_charge) == int(muons.lepCharge[mu_index])
			self.mll_dMu_plep_is_OS = not self.mll_dMu_plep_is_SS
			self.mll_dMu_plep = mu_plep_vec.M()

		if plep == 'electron' and len(elVec) > 0:
			# only 1 electron in DV then choose this lepton for dlep in veto
			if len(elVec) == 1: el_index = 0
			# choose highest pt electron in DV for dlep in veto
			elif len(elVec) == 2:
				if elVec_lepmatched[0].Pt() >  elVec_lepmatched[1].Pt(): el_index = 0
				else: el_index = 1
			# get the 4-vector for highest pT same flavour lepton as plep + prompt lepton
			el_plep_vec = elVec_lepmatched[el_index] + plep_vec
			self.mll_dEl_plep_is_SS = int(plep_charge) == int(electrons.lepCharge[el_index])
			self.mll_dEl_plep_is_OS = not self.mll_dEl_plep_is_SS
			self.mll_dEl_plep = el_plep_vec.M()

	def passes(self):
		if self.plep == "muon":
			if self.dv_type == "mumu" or self.dv_type == "emu":
				return (self.mll_dMu_plep_is_OS == True and (self.mll_dMu_plep < 80 or self.mll_dMu_plep > 100)) or self.mll_dMu_plep_is_SS == True
			else:
				return True
		elif self.plep == "electron":
			if self.dv_type == "emu" or self.dv_type == "ee":
				return (self.mll_dEl_plep_is_OS == True and (self.mll_dEl_plep < 80 or self.mll_dEl_plep > 100)) or self.mll_dEl_plep_is_SS == True
			else:
				return True
		else:
			return False


class DVLepPt:
	def __init__(self, tree, dv_type, el_pt_cut=4.5, mu_pt_cut=3):
		self.pass_pt_cut = False

		# get muons
		muons = helpers.Tracks(tree)
		muons.getMuons()
		# get lep matched 4 vector
		muVec = muons.lepmatched_lepVec

		# get electrons
		electrons = helpers.Tracks(tree)
		electrons.getElectrons()
		# get lep matched 4 vector
		elVec = electrons.lepmatched_lepVec

		if dv_type == "mumu":
			self.pass_pt_cut = muVec[0].Pt() > mu_pt_cut and muVec[1].Pt() > mu_pt_cut
		if dv_type == "emu":
			self.pass_pt_cut = elVec[0].Pt() > el_pt_cut and muVec[0].Pt() > mu_pt_cut
		if dv_type == "ee":
			self.pass_pt_cut = elVec[0].Pt() > el_pt_cut and elVec[1].Pt() > el_pt_cut

	def passes(self):
		return self.pass_pt_cut

class Mhnl:
	def __init__(self, tree, dv_type, plep, dMu, dEl, MW=80.379, fixWMass=False, hnlmasscut=20, use_truth=False, truth_pv=ROOT.TVector3(), truth_dv=ROOT.TVector3(), trks=[], use_tracks=False, invert=False):
		# Global W pole mass
		MW = 80.379
		MW2 = MW**2
		WGamma = 2.085
		self.hnlmasscut = hnlmasscut
		self.mhnl = -1
		self.alt_mhnl = -1
		self.hnlpt = -1
		self.hnleta = -99
		self.hnlphi = -99
		self.mlll = -1
		self.invert = invert

		dtrks = []
		if use_tracks:
			dtrks.append(trks[0])
			dtrks.append(trks[1])
		else:
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
			dv = ROOT.TVector3(tree.dv('x'), tree.dv('y'), tree.dv('z'))
			pv = ROOT.TVector3(tree['vertex_x'], tree['vertex_y'], tree['vertex_z'])
		else:
			pv = truth_pv
			dv = truth_dv

		p0 = ROOT.TVector3(plep.Px(), plep.Py(), plep.Pz())

		d0 = ROOT.TVector3(dtrks[0].Px(), dtrks[0].Py(), dtrks[0].Pz())
		d1 = ROOT.TVector3(dtrks[1].Px(), dtrks[1].Py(), dtrks[1].Pz())

		def findMass(pv, dv, p0, d0, d1, MW2, fixWMass):
			# Choose z direction to be along decay
			decayV = dv - pv
			z = decayV * (1.0 / decayV.Mag())

			# Visible (2 decay lepton) system
			dvis = d0 + d1
			mvis2 = 2 * (d0.Mag() * d1.Mag() - d0.Dot(d1))

			# Define plane perpendicular to direction of decay
			x = dvis.Cross(z)
			qperp = x.Mag()
			qperp2 = qperp * qperp

			# x and y unit vectors
			x = x * (1.0 / qperp)
			y = z.Cross(x)

			# Visible momentum in the (x,y,z) coordinates
			qv = ROOT.TVector3(0., qperp, dvis.Dot(z))
			Ev = np.sqrt(mvis2 + qv.Dot(qv))
			qperp3 = ROOT.TVector3(0., qperp, 0)  # not needed, just for checking
			# Prompt lepton in new (x,y,z) coordinates
			pp = ROOT.TVector3(p0.Dot(x), p0.Dot(y), p0.Dot(z))
			Ep = np.sqrt(pp.Dot(pp))

			# Terms from conservation of 4 momentum that involve various powers of the neutrino z momentum
			# A = (MW2 - mvis2)/2. - Ep*Ev - qperp2 + pp[2]*qv[2]
			B = (pp[2] + qv[2])
			E = Ep + Ev
			# A = A/E
			B = B / E

			# Min W mass possible
			alpha = qperp * B / np.sqrt(1 - B ** 2)
			qn1 = ROOT.TVector3(0, -qperp, alpha)
			En1 = qn1.Mag()
			qtot1 = pp + qv + qn1
			Etot1 = Ep + Ev + En1
			mWMin2 = Etot1 ** 2 - qtot1.Dot(qtot1)
			mWMin = np.sqrt(mWMin2)

			cdVal = 0.5 + np.arctan((mWMin2 - MW2) / MW / WGamma) / np.pi
			cdMed = (1 + cdVal) / 2
			# invert to get mass of median allowed range
			mMed2 = MW2 + MW * WGamma * np.tan(np.pi * (cdMed - 0.5))
			# mMed = MW + WGamma*np.tan(np.pi*(cdMed-0.5))

			mMed = np.sqrt(mMed2)

			# MW2fit = max(MW2, mWMin2+1)
			if fixWMass:
				MW2fit = MW2
			else:
				MW2fit = mMed2

			# MW2fit = MW2
			A = (MW2fit - mvis2) / 2. - Ep * Ev - qperp2 + pp[2] * qv[2]
			A = A / E

			# Coefficients of the quadratic to solve
			b = 2 * A * B / (B * B - 1)
			c = (A * A - qperp2) / (B * B - 1)

			arg = b * b - 4 * c
			# protect against imaginary solutions (investigate any of these)
			noSol = 0
			if arg > 0:
				rad = np.sqrt(arg)
			else:
				rad = -1000
				noSol = 1

			# These are the possible z momenta for the neutrino
			sol1 = (-b + rad) / 2
			sol2 = (-b - rad) / 2

			# Make vectors of the two z-momentum solutions for the neutrino
			qn1 = ROOT.TVector3(0, -qperp, sol1)
			En1 = qn1.Mag()
			qn2 = ROOT.TVector3(0, -qperp, sol2)
			En2 = qn2.Mag()

			# Total momentum of the decaying system
			qtot1 = qv + qn1
			qtot2 = qv + qn2

			# neutrino momentum in original coordinates 
			dn1 = y * (-1 * qperp) + z * sol1
			dn2 = y * (-1 * qperp) + z * sol2

			# mass of nuetrino+d0 lepton
			qn0 = dn1 + d0
			mn01 = np.sqrt((d0.Mag() + En1) ** 2 - qn0.Dot(qn0))
			qn0 = dn2 + d0
			mn02 = np.sqrt((d0.Mag() + En2) ** 2 - qn0.Dot(qn0))

			# And the mass of the N
			mN21 = (Ev + En1) ** 2 - qtot1.Dot(qtot1)
			mN22 = (Ev + En2) ** 2 - qtot2.Dot(qtot2)

			# Calculate the dependence of the mass solutions on the W mass
			dmN1dMW = (Ev * sol1 / np.sqrt(qperp2 + sol1 ** 2) - qv[2]) / np.sqrt(mN21) * (A + B * sol1) * np.sqrt(MW2) / \
					  (((B ** 2 - 1) * sol1 + A * B) * E)

			dmN2dMW = (Ev * sol2 / np.sqrt(qperp2 + sol2 ** 2) - qv[2]) / np.sqrt(mN22) * (A + B * sol2) * np.sqrt(MW2) / \
					  (((B ** 2 - 1) * sol2 + A * B) * E)

			# make 4-vectors in original coordinates
			pnu2 = ROOT.TLorentzVector(dn2, dn2.Mag())
			pnu1 = ROOT.TLorentzVector(dn1, dn1.Mag())
			pion_mass = 0.139  # GeV
			plep0 = ROOT.TLorentzVector()
			ptrk0 = ROOT.TLorentzVector()
			ptrk1 = ROOT.TLorentzVector()
			# using tracking assumption of pion mass for consistency with trk quantities  
			plep0.SetPxPyPzE(p0.X(), p0.Y(), p0.Z(), np.sqrt(p0.Mag() ** 2 + pion_mass ** 2))
			ptrk0.SetPxPyPzE(d0.X(), d0.Y(), d0.Z(), np.sqrt(d0.Mag() ** 2 + pion_mass ** 2))
			ptrk1.SetPxPyPzE(d1.X(), d1.Y(), d1.Z(), np.sqrt(d1.Mag() ** 2 + pion_mass ** 2))

			# assume massless tracks 
			# plep0 = ROOT.TLorentzVector( p0 ,p0.Mag())
			# ptrk0 = ROOT.TLorentzVector( d0 ,d0.Mag())
			# ptrk0 = ROOT.TLorentzVector( d1,d1.Mag())

			pHNL1 = pnu1 + ptrk0 + ptrk1
			pHNL2 = pnu2 + ptrk0 + ptrk1

			plll = plep0 + ptrk0 + ptrk1

			# set the attributes of the class
			self.mhnl = pHNL1.M()
			self.hnlpt = pHNL1.Pt()
			self.hnleta = pHNL1.Eta()
			self.hnlphi = pHNL1.Phi()
			self.mlll = plll.M()

			self.alt_mhnl = pHNL1.M()

		findMass(pv, dv, p0, d0, d1, MW2, fixWMass)

	def passes(self):
		if not self.invert:
			return self.mhnl < self.hnlmasscut
		else:
			return self.mhnl > self.hnlmasscut


class PV:
	def __init__(self, tree):

		self.pv_x = tree['vertex_x']
		self.pv_y = tree['vertex_y']
		self.pv_z = tree['vertex_z']

	def passes(self):

		if self.pv_x != -999.0 and self.pv_y != -999.0 and self.pv_z != -999.0:
			return True
		else:
			return False  # no primary vertex in the event


class MCEventType:
	"""
	Matrix elements for the trilepton process, when only the charged-current contribution is present.

	Input:
	- All masses are in GeV, and all Mandelstam variables in GeV^2.
	- pW2 = p(W)^2 = mW^2
	- The Mandelstam variables are defined as follows:
	      s13 = (p(l1)+p(l3))^2
	      s24 = (p(l2)+p(nu))^2
	  with:
	    - l1 the charged lepton produced along with the HNL;
	    - l2 the charged lepton produced in the HNL decay, on the same fermion line as the HNL;
	    - l3 the charged lepton produced in the HNL decay, on the same fermion line as the light neutrino;
	    - nu the light neutrino produced in the HNL decay.
	  There are two additional, independent Mandelstam variables which do not enter the matrix elements.

	Assumptions:
	- Light lepton (e, mu) masses are neglected.
	- Some off-shell effects coming from the finite width of the W are neglected.
	- Matrix elements are given up to a dimensionful constant, but should be consistent among themselves.
	  In particular: M2_nocorr = M2_LNC + M2_LNV.
	"""

	def __init__(self, tree, wrong_lep_order=True):
		self.weight = 1  # if not data weight is default, event is neither LNC or LNV
		self.M2_spin_corr = -1
		self.M2_nocorr = -1
		self.isLNC = False
		self.isLNV = False

		MW = 80.379  # Change the W mass to match your particle data. This is the latest PDG value.

		# Lepton number conserving charged-current trilepton process.
		def M2_LNC(MN, pW2, s13, s24):
			return s24 * (2 * MW ** 2 * pW2 * (MN ** 2 - s24) + MN ** 4 * (s13 - 2 * MW ** 2) + 2 * (s24 - s13) * MN ** 2 * MW ** 2) / (6 * MW ** 6)

		# Lepton number violating charged-current trilepton process.
		def M2_LNV(MN, pW2, s13, s24):
			return - s24 * MN ** 2 * (pW2 * (s24 - MN ** 2) + (s13 - s24) * MN ** 2 + MN ** 4 - 2 * s13 * MW ** 2) / (6 * MW ** 6)

		# Charged-current trilepton process, ignoring spin correlations between the HNL production and its decay.
		def M2_nocorr(MN, pW2, s24):
			return (s24 * (MN ** 2 - s24) * (MN ** 2 + 2 * MW ** 2) * (pW2 - MN ** 2)) / (6 * MW ** 6)

		if not tree.is_data and not tree.not_hnl_mc:
			truth_info = helpers.Truth()
			truth_info.getTruthParticles(tree)
			pW2 = truth_info.W_vec.Mag2()
			MN = truth_info.HNL_vec.M()
			# print "HNL mass ", MN
			charge_1 = truth_info.plep_charge  # charge of prompt lepton
			self.p_1 = truth_info.plep_vec  # prompt lepton
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
				self.M2_nocorr = M2_nocorr(MN=MN, pW2=pW2, s24=self.s34)  # wrong matrix used when generating pythia samples, includes lepton permutation
			else:
				self.M2_nocorr = M2_nocorr(MN=MN, pW2=pW2, s24=self.s24)  # wrong matrix used when generating pythia samples, NO lepton permutation

			self.weight = 2 * self.M2_spin_corr / self.M2_nocorr  # factor of 2 here is because M2_nocorr as calculated includes LNC + LNV decays


class LeptonTriggerMatching:
	def __init__(self, tree, lep, lepton_index):
		self.is_trigger_matched = False
		if lep == "muon":
			lep_matched = tree["muon_isTrigMatched"]
		if lep == "electron":
			lep_matched = tree["el_isTrigMatched"]
		if lep_matched[lepton_index] == 1:
			self.is_trigger_matched = True


class RequireMediumTriggerMatching:
	def __init__(self, tree, prompt_lepton_index, prompt_lepton_type, muons, electrons, dv_type):

		# check if prompt lepton is trigger matched
		prompt_is_trigger_matched = LeptonTriggerMatching(tree, prompt_lepton_type, prompt_lepton_index).is_trigger_matched
		prompt_is_medium = False
		if prompt_lepton_type == 'muon':
			prompt_is_medium = tree.get('muon_isMedium')[prompt_lepton_index] == 1
		if prompt_lepton_type == 'electron':
			prompt_is_medium = tree.get('el_LHMedium')[prompt_lepton_index] == 1

		# check if displaced is trigger matched
		displaced_0_is_trigger_matched = False
		displaced_0_is_medium = False
		displaced_1_is_trigger_matched = False
		displaced_1_is_medium = False
		if dv_type == 'mumu':
			displaced_0_is_trigger_matched = LeptonTriggerMatching(tree, 'muon', muons.lepIndex[0]).is_trigger_matched
			displaced_0_is_medium = tree.get('muon_isMedium')[muons.lepIndex[0]] == 1
			displaced_1_is_trigger_matched = LeptonTriggerMatching(tree, 'muon', muons.lepIndex[1]).is_trigger_matched
			displaced_1_is_medium = tree.get('muon_isMedium')[muons.lepIndex[1]] == 1
		if dv_type == 'ee':
			displaced_0_is_trigger_matched = LeptonTriggerMatching(tree, 'electron', electrons.lepIndex[0]).is_trigger_matched
			displaced_0_is_medium = tree.get('el_LHMedium')[electrons.lepIndex[0]] == 1
			displaced_1_is_trigger_matched = LeptonTriggerMatching(tree, 'electron', electrons.lepIndex[1]).is_trigger_matched
			displaced_1_is_medium = tree.get('el_LHMedium')[electrons.lepIndex[1]] == 1
		if dv_type == 'emu':
			displaced_0_is_trigger_matched = LeptonTriggerMatching(tree, 'muon', muons.lepIndex[0]).is_trigger_matched
			displaced_0_is_medium = tree.get('muon_isMedium')[muons.lepIndex[0]] == 1
			displaced_1_is_trigger_matched = LeptonTriggerMatching(tree, 'electron', electrons.lepIndex[0]).is_trigger_matched
			displaced_1_is_medium = tree.get('el_LHMedium')[electrons.lepIndex[0]] == 1

		# check for at least one lepton passing criteria
		self.prompt_trigger_matched_medium = prompt_is_trigger_matched and prompt_is_medium
		self.displaced_0_trigger_matched_medium = displaced_0_is_trigger_matched and displaced_0_is_medium
		self.displaced_1_trigger_matched_medium = displaced_1_is_trigger_matched and displaced_1_is_medium
		self.n_trigger_matched_medium = sum([
			self.prompt_trigger_matched_medium,
			self.displaced_0_trigger_matched_medium,
			self.displaced_1_trigger_matched_medium,
		])

	def passes(self):
		return self.n_trigger_matched_medium > 0

# class TriggerMatching_disp:
# 	def __init__(self, tree, dv_type, dMu_Index, dEl_Index):
# 		self.dlep_isTrigMatched = False
# 		if dv_type == "emu":
# 			self.dlep_isTrigMatched = tree["muon_isTrigMatched"][dMu_Index[0]] == 1 or tree["el_isTrigMatched"][dEl_Index[0]] == 1
# 		if dv_type == "ee":
# 			self.dlep_isTrigMatched = tree["el_isTrigMatched"][dEl_Index[0]] == 1 or tree["el_isTrigMatched"][dEl_Index[1]] == 1
# 		if dv_type == "mumu":
# 			self.dlep_isTrigMatched = tree["muon_isTrigMatched"][dMu_Index[0]] == 1 or tree["muon_isTrigMatched"][dMu_Index[1]] == 1


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
