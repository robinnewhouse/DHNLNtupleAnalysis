import ROOT
import uproot
import helpers

logger = helpers.getLogger('dHNLAnalysis.treenames')


class Event:
	def __init__(self, ievt, idv=None, mass=1.0, ctau=1.0):
		self.ievt = ievt
		self.idv = idv



class NewTree:
	def __init__(self, file_name, tree_name, max_entries, mass=1.0, ctau=1.0):
		"""
		Tree is the primary class that stores all information about the variables in a loaded ntuple
		and the information about wh
		:param file_name:
		:param tree_name:
		:param vtx_container:
		:param nentries:
		:param mass:
		:param ctau:
		"""
		# Set class attributes
		self.ievt = 0
		self.idv = 0
		self.max_entries = max_entries
		self.mass = mass
		self.ctau = ctau
		self.arrays = {}
		# Calculated class attributes
		# Open and load uproot trees
		self.file = uproot.open(file_name)
		self.tree = self.file[tree_name]
		self.cutflow = self.file["cutflow"]
		self.vtx_container = ""
		self.weight = self.get_weight(mass, ctau)

	def increment_event(self):
		self.ievt += 1

	def reset_event(self):
		"""Sets the event index to zero."""
		self.ievt = 0

	def increment_dv(self):
		self.idv += 1
		if self.idv > self.ndv:
			raise IndexError("DV index is too large. Is your loop resetting it?")

	def reset_dv(self):
		"""Sets the displaced vertex index to zero."""
		self.idv = 0

	def add(self, key):
		"""Loads a tree from uproot into a numpy array.
		key is the parameter to be loaded."""
		logger.debug("Accessing tree {} for the first time. Loading from file.".format(key))
		try:
			self.arrays[key] = self.tree[key].array(entrystop=self.max_entries)
		except KeyError as e:
			print("Key not found", e)

	def get(self, key):
		return self[key]

	def get_dv(self, key):
		"""
		A helper function for accessing variables specific to the displaced vertex.
		For example if you wish to access the variable secVtx_VSI_trk_pt_wrtSV and you are
		currently analyzing the VSI channel, call tree.get_dv('trk_pt_wrtSV') and this method
		will return the variable at the appropriate event and dv indices.
		For full control over key and index use the `get_at()` method.
		:param key: key of variable in tree without the secVtx_<vertex type> prefix
		:return: variable specified or awkward array of variables
		"""
		return self['{}_{}'.format(self.dv_prefix, key)]

	def get_at(self, key, ievt=None, idv=None, itrk=None):
		"""
		A method for getting value with full control over the index and variable name.
		This method will also not attempt to prefix the variable with the displaced vertex container name.
		:param key: name of the variable to be accessed
		:param ievt: specific index of event
		:param idv: specific index of displaced vertex
		:param itrk: specific index of track
		:return: variable specified or awkward array of variables
		"""
		if key not in self.arrays:
			# Load the tree from uproot if it hasn't been loaded yet.
			self.add(key)
		val = self.arrays[key][ievt]
		if idv is not None:
			val = val[idv]
		if itrk is not None:
			val = val[itrk]
		return val

	def __getitem__(self, key):
		"""
		__getitem__ defines the behaviour of the [] operator.
		This method will try to return an event-level variable at the current event index.
		If the event level variable is an awkward array, it will attempt to return the value at the dv index.
		:param key: The variable to be accessed
		:return: The variable for the event and displaced vertex
		"""
		if key.startswith(b'secVtx'):
			# This variable is a likely dv type. Return variable at dv index.
			return self.get_at(key, self.ievt, self.idv)
		else:
			# This variable is likely event type (e.g. ndv)
			return self.get_at(key, self.ievt)

	def __setitem__(self, key, value):
		raise AttributeError("Can't set attribute")

	def get_weight(self, mass, ctau, lnv=False):
		"""
		Calculates the weight of the event based on the Gronau parametrization
		https://journals.aps.org/prd/abstract/10.1103/PhysRevD.29.2539
		Sets the weight of events for this tree
		:param mass: HNL sample mass
		:param ctau: HNL sample lifetime
		:param lnv: Use Lepton Number Violating calculation
		:return:
		"""
		if self.is_data:  # you are running on data
			self.weight = 1
		else:  # you are running on MC file
			if self.mass == -1 or ctau == -1:  # MC weighting error
				logger.debug("Can't determine the mass and lifetime of signal sample. MC weight will be set to 1!!")
				self.weight = 1
			else:
				mW = 80.379  # mass of W boson in GeV
				U2Gronau = 4.49e-12 * 3e8 * mass ** (-5.19) / (ctau / 1000)  # LNC prediction
				if (lnv):
					U2 = 0.5 * U2Gronau
				else:
					U2 = U2Gronau

				xsec = 20.6e6 * U2 * ((1 - (mass / mW) ** 2) ** 2) * (1 + (mass ** 2) / (2 * mW ** 2))  # in fb
				self.weight = 1 * xsec / self.numentries  # scale to 1 fb^-1  of luminosity

		return self.weight

	@property
	def numentries(self):
		return self.tree.numentries

	@property
	def is_data(self):
		"""Checks if this tree represents real data"""
		return "truthVtx_x" not in self.tree.keys()

	@property
	def ndv(self):
		return self['n'+self.dv_prefix]

	@property
	def ntrk(self):
		return self.get_dv('ntrk')

	@property
	def dv_prefix(self):
		if not self.vtx_container:
			raise AttributeError("vtx_container is not set. please specify a vtx_container (e.g. VSI or VSI_Leptons")
		return "secVtx_{}".format(self.vtx_container)

class Tree():
	def __init__(self, fileName, treeName, vtx_container, nentries):
		logger.info("Importing trees from ntuple for %s."%vtx_container)
		self.file = uproot.open(fileName)	
		self.tree = self.file[treeName]
		
		if "truthVtx_x" in self.file[treeName].keys(): 
			self.isData = False
		else: 
			self.isData = True

		DVprefix = "secVtx_" + vtx_container

		self.cutflow = self.file["cutflow"]
		
		self.allEvt =  self.cutflow[1]
		self.npassTrig = int(self.cutflow[4])

		# -----------------------------------------------------------------------
		self.evtNum = self.tree["eventNumber"].array(entrystop=nentries)
		self.passedtriggers = self.tree["passedTriggers"].array(entrystop=nentries)
		self.mumufilter = self.tree["passesHnlMuMuFilter"].array(entrystop=nentries)
		self.muelfilter = self.tree["passesHnlMuElFilter"].array(entrystop=nentries)
		self.elelfilter = self.tree["passesHnlElElFilter"].array(entrystop=nentries)
		self.elmufilter = self.tree["passesHnlElMuFilter"].array(entrystop=nentries)

		if "vertex_x" in self.file[treeName]: 
			self.pvx = self.tree["vertex_x"].array(entrystop=nentries)  
			self.pvy = self.tree["vertex_y"].array(entrystop=nentries)  
			self.pvz = self.tree["vertex_z"].array(entrystop=nentries) 
		
		# -----------------------------------------------------------------------
		# DV track variables
		# -----------------------------------------------------------------------
		self.trk_muonindex = self.tree[DVprefix + "_trk_muonIndex"].array(entrystop=nentries)
		self.trk_elindex = self.tree[DVprefix + "_trk_electronIndex"].array(entrystop=nentries)
		self.trackpt = self.tree[DVprefix + "_trk_pt_wrtSV"].array(entrystop=nentries) 
		self.trackpt_wrt_pv = self.tree[DVprefix + "_trk_pt"].array(entrystop=nentries) 
		self.tracketa = self.tree[DVprefix + "_trk_eta_wrtSV"].array(entrystop=nentries) 
		self.trackphi = self.tree[DVprefix + "_trk_phi_wrtSV"].array(entrystop=nentries)
		self.trackmass = self.tree[DVprefix + "_trk_M"].array(entrystop=nentries)
		self.trackcharge = self.tree[DVprefix + "_trk_charge"].array(entrystop=nentries)
		self.trackchi2 = self.tree[DVprefix + "_trk_chi2_toSV"].array(entrystop=nentries)
		
		self.trackd0 = self.tree[DVprefix + "_trk_d0"].array(entrystop=nentries) #d0 and z0 are track quantites calculated wrt IP we dont want to look at the _wrtSV variable!
		self.trackz0 = self.tree[DVprefix + "_trk_z0"].array(entrystop=nentries) 

		# self.trackz0 = self.tree[DVprefix + "_trk_z0"].array(entrystop=nentries) 
		# self.trackd0 = self.tree[DVprefix + "_trk_d0"].array(entrystop=nentries) 

		# -----------------------------------------------------------------------
		# DV reco variables
		# -----------------------------------------------------------------------
		self.dvx = self.tree[DVprefix + "_x"].array(entrystop=nentries)  
		self.dvy = self.tree[DVprefix + "_y"].array(entrystop=nentries)  
		self.dvz = self.tree[DVprefix + "_z"].array(entrystop=nentries) 
		self.dvr = self.tree[DVprefix + "_r"].array(entrystop=nentries) 
		self.dvmass = self.tree[DVprefix + "_mass"].array(entrystop=nentries)  
		self.dvpt = self.tree[DVprefix + "_pt"].array(entrystop=nentries)
		self.dveta = self.tree[DVprefix + "_eta"].array(entrystop=nentries)  
		self.dvphi = self.tree[DVprefix + "_phi"].array(entrystop=nentries) 
		self.dvminOpAng = self.tree[DVprefix + "_minOpAng"].array(entrystop=nentries) 
		self.dvmaxOpAng = self.tree[DVprefix + "_maxOpAng"].array(entrystop=nentries)   
		self.dvntrk = self.tree[DVprefix + "_ntrk"].array(entrystop=nentries) 
		self.dvdistFromPV = self.tree[DVprefix + "_distFromPV"].array(entrystop=nentries) 
		self.dvcharge = self.tree[DVprefix + "_charge"].array(entrystop=nentries) 
		self.dvchi2 = self.tree[DVprefix + "_chi2"].array(entrystop=nentries) 
		# -----------------------------------------------------------------------
		
		# -----------------------------------------------------------------------
		# muon variables 
		# -----------------------------------------------------------------------
		self.muonindex = self.tree["muon_index"].array(entrystop=nentries)
		self.muonpt = self.tree["muon_pt"].array(entrystop=nentries) 
		self.muoneta = self.tree["muon_eta"].array(entrystop=nentries)  
		self.muonphi = self.tree["muon_phi"].array(entrystop=nentries) 
		self.muond0 = self.tree["muon_trkd0"].array(entrystop=nentries) 
		self.muonz0 = self.tree["muon_trkz0"].array(entrystop=nentries)
		self.muonpx = self.tree["muon_px"].array(entrystop=nentries) 
		self.muonpy = self.tree["muon_py"].array(entrystop=nentries) 
		self.muonpz = self.tree["muon_pz"].array(entrystop=nentries) 
		self.muonmass = self.tree["muon_m"].array(entrystop=nentries)  
		self.tightmu = self.tree["muon_isTight"].array(entrystop=nentries)
		self.mediummu = self.tree["muon_isMedium"].array(entrystop=nentries)
		self.loosemu = self.tree["muon_isLoose"].array(entrystop=nentries)
		self.muontype = self.tree["muon_type"].array(entrystop=nentries)
		self.muonpassPfilter = self.tree["muon_passesPromptCuts"].array(entrystop=nentries) 
		self.muontrigmatched = self.tree["muon_isTrigMatched"].array(entrystop=nentries)
		# -----------------------------------------------------------------------
		
		# -----------------------------------------------------------------------
		# electron variables
		# -----------------------------------------------------------------------
		self.elindex = self.tree["el_index"].array(entrystop=nentries)
		self.elpt = self.tree["el_pt"].array(entrystop=nentries)  
		self.eleta = self.tree["el_eta"].array(entrystop=nentries)  
		self.elphi = self.tree["el_phi"].array(entrystop=nentries)  
		self.eld0 = self.tree["el_trkd0"].array(entrystop=nentries) 
		self.elz0 = self.tree["el_trkz0"].array(entrystop=nentries)
		self.elmass = self.tree["el_m"].array(entrystop=nentries)  
		self.tightel = self.tree["el_LHTight"].array(entrystop=nentries)
		self.mediumel = self.tree["el_LHMedium"].array(entrystop=nentries)
		self.looseel = self.tree["el_LHLoose"].array(entrystop=nentries)
		self.elpassPfilter = self.tree["el_passesPromptCuts"].array(entrystop=nentries) 
		# -----------------------------------------------------------------------



		# -----------------------------------------------------------------------
		# DV truth variables
		# -----------------------------------------------------------------------
		if self.isData == False:
			self.truth_x = self.tree["truthVtx_x"].array(entrystop=nentries)  
			self.truth_y = self.tree["truthVtx_y"].array(entrystop=nentries)  
			self.truth_z = self.tree["truthVtx_z"].array(entrystop=nentries) 
			self.truth_dvr = self.tree["truthVtx_r"].array(entrystop=nentries) 
			self.truth_dvmass = self.tree["truthVtx_mass"].array(entrystop=nentries)  
			self.truth_dvpt = self.tree["truthVtx_pt"].array(entrystop=nentries)
			self.truth_dveta = self.tree["truthVtx_eta"].array(entrystop=nentries)  
			self.truth_dvphi = self.tree["truthVtx_phi"].array(entrystop=nentries) 
			self.truth_outP_pdgId = self.tree["truthVtx_outP_pdgId"].array(entrystop=nentries)
			self.truth_outP_pt = self.tree["truthVtx_outP_pt"].array(entrystop=nentries)
			self.truth_outP_eta = self.tree["truthVtx_outP_eta"].array(entrystop=nentries)
			self.truth_outP_phi = self.tree["truthVtx_outP_phi"].array(entrystop=nentries)
			self.truth_outP_m = self.tree["truthVtx_outP_M"].array(entrystop=nentries)
			self.truth_parent_pdgId = self.tree["truthVtx_parent_pdgId"].array(entrystop=nentries)
			self.truth_parent_pt = self.tree["truthVtx_parent_pt"].array(entrystop=nentries)
			self.truth_parent_eta = self.tree["truthVtx_parent_eta"].array(entrystop=nentries)
			self.truth_parent_phi = self.tree["truthVtx_parent_phi"].array(entrystop=nentries)
			self.truth_parent_m = self.tree["truthVtx_parent_M"].array(entrystop=nentries)
			# self.dvntrk = self.tree["truthVtx_ntrk"].array(entrystop=nentries) 
			# self.dvdistFromPV = self.tree["truthVtx_distFromPV"].array(entrystop=nentries) 
			# self.dvcharge = self.tree["truthVtx_charge"].array(entrystop=nentries) 
			# self.dvcharge = self.tree["truthVtx_chi2"].array(entrystop=nentries) 
		# -----------------------------------------------------------------------
		logger.info("Done importing trees!")
