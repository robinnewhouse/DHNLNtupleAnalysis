# from ROOT import*
import ROOT
import numpy as np
import helpers
import logging

logger = helpers.getLogger('dHNLAnalysis.selections', level=logging.WARNING)


class Trigger():
	def __init__(self, tree, trigger):
		self.tree = tree

		# trigger lists taken from https://acode-browser1.usatlas.bnl.gov/lxr/source/athena/PhysicsAnalysis/SUSYPhys/LongLivedParticleDPDMaker/share/PhysDESDM_HNL.py?v=21.0#0008
		apiSingleMuonTriggerlist = ["HLT_mu20_iloose_L1MU15", "HLT_mu24_iloose", "HLT_mu24_ivarloose", "HLT_mu24_ivarmedium", "HLT_mu26_imedium", "HLT_mu26_ivarmedium", "HLT_mu40", "HLT_mu50", "HLT_mu60_0eta105_msonly"]
		apiSingleElectronTriggerlist = ["HLT_e24_lhmedium_L1EM20VH", "HLT_e24_lhtight_nod0_ivarloose", "HLT_e26_lhtight_nod0", "HLT_e26_lhtight_nod0_ivarloose", "HLT_e60_lhmedium_nod0", "HLT_e60_lhmedium", "HLT_e60_medium",
										"HLT_e120_lhloose", "HLT_e140_lhloose_nod0", "HLT_e300_etcut"]
		apiMultiMuonTriggerlist = ["HLT_2mu10",  "HLT_2mu10_nomucomb",  "HLT_2mu14",  "HLT_2mu14_nomucomb",  "HLT_mu20_nomucomb_mu6noL1_nscan03",  "HLT_mu20_mu8noL1",  "HLT_mu22_mu8noL1",  "HLT_mu22_mu8noL1_calotag_0eta010",  "HLT_3mu4",  "HLT_mu6_2mu4",  "HLT_3mu6",  "HLT_3mu6_msonly",  "HLT_mu20_2mu4noL1",  "HLT_4mu4",  "HLT_mu11_nomucomb_2mu4noL1_nscan03_L1MU11_2MU6",  "HLT_mu11_nomucomb_2mu4noL1_nscan03_L1MU11_2MU6_bTau",  "HLT_mu20_msonly_mu10noL1_msonly_nscan05_noComb",  "HLT_mu11_nomucomb_mu6noL1_nscan03_L1MU11_2MU6_bTau",  "HLT_mu6_nomucomb_2mu4_nomucomb_bTau_L1MU6_3MU4",  "HLT_2mu6_nomucomb_mu4_nomucomb_bTau_L12MU6_3MU4"]
		apiMultiElectronTriggerlist = ["HLT_2e12_lhloose_L12EM10VH", "HLT_2e15_lhvloose_nod0_L12EM13VH", "HLT_e17_lhloose_2e9_lhloose", "HLT_e17_lhloose_nod0_2e9_lhloose_nod0", "HLT_e17_lhloose_nod0_2e10_lhloose_nod0_L1EM15VH_3EM8VH", "HLT_2e17_lhvloose_nod0", "HLT_2e17_lhvloose_nod0_L12EM15VHI", "HLT_2e24_lhvloose_nod0", "HLT_e24_lhvloose_nod0_2e12_lhvloose_nod0_L1EM20VH_3EM10VH"]
		apiElectronMuonTriggerlist = ["HLT_e17_lhloose_mu14", "HLT_e17_lhloose_nod0_mu14", "HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1", "HLT_e26_lhmedium_nod0_mu8noL1", "HLT_e7_lhmedium_nod0_mu24", "HLT_e12_lhloose_nod0_2mu10", "HLT_2e12_lhloose_mu10", "HLT_2e12_lhloose_nod0_mu10", "HLT_e7_lhmedium_mu24", "HLT_e12_lhloose_2mu10"]
		
		if trigger == "muononly":
			self.allowed_trigger_list = apiSingleMuonTriggerlist
		elif trigger == "electrononly":
			self.allowed_trigger_list = apiSingleElectronTriggerlist
		elif trigger == "all":
			self.allowed_trigger_list = apiSingleMuonTriggerlist + apiSingleElectronTriggerlist
		else:
			self.allowed_trigger_list = list(trigger)

	def passes(self):
		# evaluate whether the event passed the trigger
		# This method checks if there is any overlap between sets a and b
		# https://stackoverflow.com/questions/3170055/test-if-lists-share-any-items-in-python
		event_triggers = self.tree['passedTriggers']
		return not set(event_triggers).isdisjoint(self.allowed_trigger_list)


class Filter():
	def __init__(self, tree, filter_type):
		self.tree = tree
		self.filter_type = filter_type

	def passes(self):
		if self.filter_type == "mu-mu":
			return self.tree['passesHnlMuMuFilter']

		if self.filter_type == "mu-el":
			return self.tree['passesHnlMuElFilter']

		if self.filter_type == "el-mu":
			return self.tree['passesHnlElMuFilter']

		if self.filter_type == "el-el":
			return self.tree['passesHnlElElFilter']

		if self.filter_type == "4-filter":
			return (self.tree['passesHnlMuMuFilter']
					or self.tree['passesHnlElMuFilter']
					or self.tree['passesHnlElElFilter']
					or self.tree['passesHnlMuElFilter'])

		if self.filter_type == "3-filter":
			return (self.tree['passesHnlMuMuFilter']
					or self.tree['passesHnlElMuFilter']
					or self.tree['passesHnlElElFilter'])

		if self.filter_type == "2-filter":
			return (self.tree['passesHnlMuMuFilter']
					or self.tree['passesHnlElMuFilter'])

		if self.filter_type == "1-filter":
			return self.tree['passesHnlMuMuFilter']

class InvertedPromptLepton():
	def __init__(self, tree, d0_cut=3.0, z0_sin_theta_cut=0.5, pt_cut=25.0):
		self.tree = tree
		self.n_prompt_leptons = 0 
		self.n_prompt_muons = 0 
		self.n_prompt_electrons = 0 
		n_muons = len(self.tree['muon_pt'])
		n_electrons = len(self.tree['el_pt'])

		for imu in range(n_muons):
			# make sure the muon is at least loose
			if not ((self.tree['muon_isTight'][imu] == 1) or
					(self.tree['muon_isMedium'][imu] == 1) or
					(self.tree['muon_isLoose'][imu]) == 1):
				# muon doesn't satisfy any quality, ignore it
				continue

			# check muon pt
			mupt = self.tree['muon_pt'][imu]
			# if mupt > pt_cut:
			# 	# muon satisfies "fast" lepton requirements
			# 	self.n_prompt_muons += 1
			# 	continue
			
			# check muon promptness
			mumass = self.tree['muon_m'][imu]
			mueta = self.tree['muon_eta'][imu]
			muphi = self.tree['muon_phi'][imu]
			muz0 = self.tree['muon_trkz0'][imu]
			mud0 = self.tree['muon_trkd0'][imu]

			muVec_i = ROOT.TLorentzVector()
			muVec_i.SetPtEtaPhiM(mupt, mueta, muphi, mumass)
			sintheta = np.sin(muVec_i.Theta())
			# muon satisfies prompt lepton requirements
			if (abs(mud0) < d0_cut) and (abs(muz0*sintheta) < z0_sin_theta_cut):
				self.n_prompt_muons += 1
				continue

		for iel in range(n_electrons):
			# make sure the electron is at least loose
			if not ((self.tree['el_LHTight'][iel] == 1) or
					(self.tree['el_LHMedium'][iel] == 1) or
					(self.tree['el_LHLoose'][iel]) == 1):
				# electron doesn't satisfy any quality, ignore it
				continue
			
			# check electron pt
			elpt = self.tree['el_pt'][iel]
			# if elpt > pt_cut:
			# 	# electron satisfies "fast" lepton requirements
			# 	self.n_prompt_electrons += 1
			# 	continue
			
			# check electron promptness
			elmass = self.tree['el_pt'][iel]
			eleta = self.tree['el_eta'][iel]
			elphi = self.tree['el_phi'][iel]
			elz0 = self.tree['el_trkz0'][iel]
			eld0 = self.tree['el_trkd0'][iel]

			elVec_i = ROOT.TLorentzVector()
			elVec_i.SetPtEtaPhiM(elpt, eleta, elphi, elmass)
			sintheta = np.sin(elVec_i.Theta())
			if (abs(eld0) < d0_cut) and (abs(elz0 * sintheta) < z0_sin_theta_cut):
				# electron satisfies prompt lepton requirements
				self.n_prompt_electrons += 1
				continue

		self.n_prompt_leptons = self.n_prompt_electrons + self.n_prompt_muons

	def passes(self):
		return self.n_prompt_leptons == 0


class Plepton():
	def __init__(self, tree, lepton, quality="tight", mindR=0.05):
		self.tree = tree
		self.lepton = lepton
		self.quality = quality 
		self.mindR = mindR

		self.plepVec = ROOT.TLorentzVector(0, 0, 0, 0)
		self.plepd0 = -2000
		self.plepz0 = -2000
		ndv = tree.ndv
		nleps = 0
		self.nPlep = 0

		lepquality = ""
		passPfilter = False
		if self.lepton == "muon":
			if self.quality == "tight":  # tight muon is requested
				lepquality = 'muon_isTight'
			if self.quality == "medium":
				lepquality = 'muon_isMedium'
			if self.quality == "loose":
				lepquality = 'muon_isLoose'

			nleps = len(self.tree['muon_pt'])
			passPfilter = self.tree['muon_passesPromptCuts']

		if self.lepton == "electron":
			if self.quality == "tight":  # tight electron is requested
				lepquality = 'el_LHTight'
			if self.quality == "medium":
				lepquality = 'el_LHMedium'
			if self.quality == "loose":
				lepquality = 'el_LHLoose'


			nleps = len(self.tree['el_pt'])
			passPfilter = self.tree['el_passesPromptCuts']

		# variable for the highest pt lepton				
		self.highestpt_lep = ROOT.TLorentzVector(0, 0, 0, 0)
		self.highestpt_lep_d0 = -2000
		self.highestpt_lep_z0 = -2000

		for ilep in range(nleps): 
			overlap = False
			plepVec_i = ROOT.TLorentzVector()

			if self.lepton == "muon": 
				pt = self.tree['muon_pt'][ilep]
				eta = self.tree['muon_eta'][ilep]
				phi = self.tree['muon_phi'][ilep]
				mass = self.tree['muon_m'][ilep]
				plepVec_i.SetPtEtaPhiM(pt, eta, phi, mass)

				lepd0 = self.tree['muon_trkd0'][ilep]
				lepz0 = self.tree['muon_trkz0'][ilep]

			if self.lepton == "electron":
				pt = self.tree['el_pt'][ilep]
				eta = self.tree['el_eta'][ilep]
				phi = self.tree['el_phi'][ilep]
				mass = self.tree['el_m'][ilep]
				plepVec_i.SetPtEtaPhiM(pt, eta, phi, mass)

				lepd0 = self.tree['el_trkd0'][ilep]
				lepz0 = self.tree['el_trkz0'][ilep]

			# check if the plep passes the DRAW filter and passes quality before looping over tracks
			# changed to be careful with negative electron quality values # RN
			passes_lep_quality = lepquality == "" or self.tree[lepquality][ilep] > 0
			if passPfilter[ilep] and passes_lep_quality:
				for idv in range(ndv):
					leptracks = helpers.Tracks(self.tree)
					# trackevt = helpers.Event(self.evt.tree, self.evt.ievt, idv)
					leptracks.getTracks()
					dlepVec = leptracks.lepVec
					ndtracks = len(dlepVec)
						
					for itr in range(ndtracks): # check overlap with DVs
						dR = dlepVec[itr].DeltaR(plepVec_i)
						if dR < self.mindR:  # set overlap to true if muon overlaps with displaced track
							overlap = True
		
				if overlap == False: # if lepton doesnt overlap with and DV tracks
					sintheta = np.sin(plepVec_i.Theta())
					if lepd0 < 3 and lepz0*sintheta < 0.5: # if lepton pass the track significance cuts 
						self.nPlep = self.nPlep + 1 
						if (plepVec_i.Pt() > self.highestpt_lep.Pt()): # if pt is larger then the previous prompt lepton found 
							self.highestpt_lep = plepVec_i  #get highest pt prompt lepton!
							self.highestpt_lep_d0 = lepd0
							self.highestpt_lep_z0 = lepz0

						#for trigger matching
						# if self.evt.tree.muontrigmatched[self.evt.ievt][ilep] == 0:
						# 	print "is muon trig matched?", self.evt.tree.muontrigmatched[self.evt.ievt][ilep]
						# 	print "pt of the highest pt TIGHT muon: ", self.highestpt_plep.Pt() 
						# 	print "muon pt: ", self.evt.tree.muonpt[self.evt.ievt]
						# 	print "trigger matched: ", self.evt.tree.muontrigmatched[self.evt.ievt]
						# 	print "lepton quality: ", lepquality[self.evt.ievt]
						# 	print "muon type: ", self.evt.tree.muontype[self.evt.ievt]

	def passes(self):
		if self.nPlep > 0 and self.highestpt_lep.Pt() > 0: 
			self.plepVec = self.highestpt_lep
			self.plepd0 = self.highestpt_lep_d0
			self.plepz0 = self.highestpt_lep_z0
			return True
		else: 
			return False
		
class Alpha():
	def __init__(self, tree, max_alpha=0.01):
		self.tree = tree
		self.max_alpha = max_alpha
		# calculate alpha
		# secondary vertex position vector
		sv_vector = ROOT.TVector3(self.tree.dv('x'),
								  self.tree.dv('y'),
								  self.tree.dv('z'))
		# primary vertex position vector
		pv_vector = ROOT.TVector3(self.tree['vertex_x'],
								  self.tree['vertex_y'],
								  self.tree['vertex_z'])
		# vector from pv to sv
		pv_sv_vector = sv_vector - pv_vector

		# vector difference between momentum vector and position vector		
		alpha = pv_sv_vector.Phi() - self.tree.dv('phi')
		# put in -pi to pi range
		self.alpha = (alpha + np.pi/2) % np.pi*2 - np.pi

	def passes(self):
		# pass if alpha is less than the sent cut
		return abs(self.alpha) < self.max_alpha

# class nDV():
# 	def __init__(self, evt):
# 		self.evt = evt
#
# 	def passes(self):
# 		return len(self.evt.tree.dvx[self.evt.ievt]) > 0
#


class DVradius():
	def __init__(self, tree):
		self.tree = tree
		self.rdv = -1
		if self.tree.ntrk > 0:
			dx = self.tree.dv('x')
			dy = self.tree.dv('y')
			self.rdv = np.sqrt(dx**2 + dy**2)

	def passes(self, rdv_min=4, rdv_max=300):
		if self.rdv > rdv_min and self.rdv < rdv_max:
			return True
		else:
			return False




class DVntracks():
	def __init__(self, tree, ntrk=2, decaymode="leptonic"):
		self.tree = tree
		self.ntrk = ntrk
		self.decaymode = decaymode

		self.ntracks = -1 

		if self.decaymode == "leptonic":
			self.ntracks = self.tree.ntrk

	def passes(self):
		if self.ntracks == self.ntrk:
			return True
		else:
			return False


class ChargeDV():
	def __init__(self, tree, sel="OS", decaymode="leptonic"):
		self.tree = tree
		self.decaymode = decaymode
		self.sel = sel
		self.ntracks = -1
		self.charge_trk1 = -111  # dont make default -1 since that's a valid charge! :)
		self.charge_trk2 = -222  # dont make default -1 since that's a valid charge! :)
		# also, don't make the same since the equality is checked later

		if self.decaymode == "leptonic":
			self.ntracks = self.tree.ntrk

			if self.ntracks == 2: 
				self.charge_trk1 = self.tree.dv('trk_charge')[0]
				self.charge_trk2 = self.tree.dv('trk_charge')[1]

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
					return True
				else:
					return False
			else:
				return False

		elif self.dv_type == "ee":
			if self.nel == 2: 
				return True
			else:
				return False

		elif self.dv_type == "mumu-notcomb":
			if self.nmu == 2:
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
		self.tree = tree
		self.decaymode = decaymode
		self.quality = quality 

		if self.decaymode == "leptonic": 
			self.ntracks = self.tree.ntrk
			self.ntight = 0
			self.DVmuons = []
			self.DVelectrons = []

			muons = helpers.Tracks(self.tree)
			muons.getMuons()

			electrons = helpers.Tracks(self.tree)
			electrons.getElectrons()
			self.nel = len(electrons.lepVec)

			self.nmu_tight = 0
			self.nel_tight = 0

			self.ndvmu = len(muons.lepVec)
			self.ndvel = len(electrons.lepVec)
		
			for imu in range(self.ndvmu):
				muindex = muons.lepIndex[imu]
				muisTight = self.tree['muon_isTight'][muindex]
				if muisTight: 
					self.nmu_tight = self.nmu_tight + 1

			for iel in range(self.ndvel):
				elindex = electrons.lepIndex[iel]
				elisTight = self.tree['el_LHTight'][elindex]
				if elisTight:
					self.nel_tight = self.nel_tight + 1

	def passes(self):
		if self.quality == "2-tight":
			# print self.nmu_tight
			if self.nmu_tight == 2 or self.nel_tight == 2 or (self.nmu_tight == 1 and self.nel_tight == 1):
				return True
			else:
				return False

		if self.quality == "1-tight":
			if self.nmu_tight > 0 or self.nel_tight > 0:
				return True
			else:
				return False


class Cosmicveto():
	def __init__(self, tree, decaymode="leptonic", cosmicvetocut=0.05):
		self.tree = tree
		self.decaymode = decaymode
		self.cosmicvetocut = cosmicvetocut

		self.separation = -1

		if self.decaymode == "leptonic":
			ntracks = self.tree.ntrk
			if ntracks == 2:
				sumeta = self.tree.dv('trk_eta_wrtSV')[0] + self.tree.dv('trk_eta_wrtSV')[1]
				dphi = abs(self.tree.dv('trk_phi_wrtSV')[0] - self.tree.dv('trk_phi_wrtSV')[1])

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
		self.plll = ROOT.TLorentzVector(0, 0, 0, 0)

		if self.decaymode == "leptonic":

			if self.dv_type == "emu":
				self.plll = self.plep + self.dEl[0] + self.dMu[0]

			if self.dv_type == "mumu":
				self.plll = self.plep + self.dMu[0] + self.dMu[1]

			if self.dv_type == "ee":
				self.plll = self.plep + self.dEl[0] + self.dEl[1]

			self.mlll = self.plll.M()

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
		self.tree = tree
		self.decaymode = decaymode
		self.dvmasscut = dvmasscut
		self.dvmass = self.tree.dv('mass')

	def passes(self):
		if self.dvmass > self.dvmasscut:
			return True
		else:
			return False


class Mhnl():
	def __init__(self, tree, plep, trks, hnlmasscut=4):
		self.tree = tree
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
				logger.ERROR("Roatating vectors did not work!! Check HNL mass calculation.")
			return r_new

		def unrotate_vector(r,v):
			r_new = v

			rotation_axis = ROOT.TVector3(-1*r.Y(),r.X(),0.0)
			rotation_angle = r.Theta() 
			r_new.Rotate(rotation_angle,rotation_axis)
		
			return r_new

		# primary vertex vector
		dv_vec = ROOT.TVector3(self.tree.dv('x'),
							   self.tree.dv('y'),
							   self.tree.dv('z'))
		# primary vertex position vector
		pv_vec = ROOT.TVector3(self.tree['vertex_x'],
							   self.tree['vertex_y'],
							   self.tree['vertex_z'])

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



class PV():
	def __init__(self, tree):
		self.tree = tree
		self.pv_x = self.tree['vertex_x']
		self.pv_y = self.tree['vertex_y']
		self.pv_z = self.tree['vertex_z']

	def passes(self):
		
		if (self.pv_x != -999.0 and self.pv_y != -999.0 and self.pv_z != -999.0 ): 
			return True
		else: 
			return False # no primary vertex in the event


class SumTrack:
	def __init__(self, tree):
		self.tree = tree
		self.sum_track_pt = 0
		n_tracks = self.tree.ntrk
		self.sum_track_pt_wrt_pv = 0
		self.sum_track_charge = 0

		for k in range(n_tracks):
			self.sum_track_pt += self.tree.dv('trk_pt_wrtSV')[k]
			self.sum_track_pt_wrt_pv += self.tree.dv('trk_pt')[k]
			self.sum_track_charge += self.tree.dv('trk_charge')[k]