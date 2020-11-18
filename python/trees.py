import uproot
import helpers
import logging




class Tree:
	def __init__(self, file_name, tree_name, max_entries, debug_level, mc_campaign=None, mass=1.0, ctau=1.0,notHNLmc=False):
		"""
		Tree is the primary class that stores all information about the variables in a loaded ntuple
		and the information about the indices of the current event (ievt) and displaced vertex (idv).
		The variables can be accessed with simple getter functions without specifying the index.
		The index incrementation will be handled in the analysis loop.

		A tree is loaded from the ntuple file only when it is needed for the first time, hopefully
		reducing overhead time in case an analysis does not need every tree.

		The Tree class also contains information about the mass, lifetime, and event weight.

		:param file_name: Name of input ntuple file to be read.
		:param tree_name: The name of the tree in the root file. often outTree or outTree+<systematic variation>
		:param max_entries: The maximum number of entries to load when reading a tree from the file.
		:param mass: The HNL mass of the sample
		:param ctau: The mean lifetime of the sample
		"""
		# Set class attributes
		self.logger = helpers.getLogger('dHNLAnalysis.trees',level=debug_level)
		self.ievt = 0
		self.idv = 0
		self.max_entries = max_entries
		self.mass = mass
		self.ctau = ctau
		self.mc_campaign = mc_campaign
		self.arrays = {}
		# Calculated class attributes
		# Open and load uproot trees
		self.file = uproot.open(file_name)
		self.tree = self.file[tree_name]
		self.cutflow = self.file["cutflow"] 
		self.all_entries = self.cutflow[2]  
		self.init_entries = self.cutflow[2]  # init entries
		self.vtx_container = ""
		self.notHNLmc = notHNLmc

	
		


	def increment_event(self):
		self.ievt += 1
		if self.ievt > self.max_entries:
			raise IndexError("Event index is too large. Is your loop resetting it?")

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
		"""Loads a tree from uproot into a numpy array."""
		self.logger.debug("Accessing tree {} for the first time. Loading from file.".format(key))
		try:
			self.arrays[key] = self.tree[key].array(entrystop=self.max_entries)
		except KeyError as e:
			raise KeyError("Key not found: {} Make sure this key is in your ntuple.".format(key))

	def get(self, key):
		"""A helper function if you want to avoid the [] operator."""
		return self[key]

	def dv(self, key):
		"""
		A helper function for accessing variables specific to the displaced vertex.
		For example if you wish to access the variable secVtx_VSI_trk_pt_wrtSV and you are
		currently analyzing the VSI channel, call tree.dv('trk_pt_wrtSV') and this method
		will return the variable at the appropriate event and dv indices with the `sexVtx_VSI_` prefix.
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
		if key.startswith('secVtx'):
			# This variable is a likely dv type. Return variable at dv index.
			return self.get_at(key, self.ievt, self.idv)
		else:
			# This variable is likely event type (e.g. ndv)
			return self.get_at(key, self.ievt)

	def __setitem__(self, key, value):
		raise AttributeError("Can't set attribute")

	@property
	def numentries(self):
		return self.tree.numentries

	@property
	def is_data(self):
		"""Checks if this tree represents real data"""
		return b'truthVtx_x' not in self.tree.keys()

	@property
	def ndv(self):
		return self['n' + self.dv_prefix]

	@property
	def ntrk(self):
		return self.dv('ntrk')

	@property
	def dv_prefix(self):
		if not self.vtx_container:
			raise AttributeError("vtx_container is not set. please specify a vtx_container (e.g. VSI or VSI_Leptons")
		return "secVtx_{}".format(self.vtx_container)
