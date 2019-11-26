import ROOT
import uproot

class Tree():
	def __init__(self, fileName,treeName):
		self.file = uproot.open(fileName)
		self.tree = self.file[treeName]
		
		#define variables
		self.muonindex = self.tree["secVtx_trk_muonIndex"].array()
		self.elindex = self.tree["secVtx_trk_electronIndex"].array()
		self.dvmass = self.tree["secVtx_mass"].array() 
		self.trackpt = self.tree["secVtx_trk_pt"].array() 
		self.tracketa = self.tree["secVtx_trk_eta"].array() 
		self.trackphi = self.tree["secVtx_trk_phi"].array()
		self.tracke = self.tree["secVtx_trk_E"].array() 
		self.trackd0 = self.tree["secVtx_trk_d0"].array() 
		self.trackz0 = self.tree["secVtx_trk_d0"].array() 
		self.trackcharge = self.tree["secVtx_trk_charge"].array()
		self.dvx = self.tree["secVtx_x"].array()  
		self.dvy = self.tree["secVtx_y"].array()  
		self.dvz = self.tree["secVtx_z"].array() 
		self.muonpt = self.tree["muon_pt"].array() 
		self.muoneta = self.tree["muon_pt"].array()  
 		self.muonphi = self.tree["muon_phi"].array() 
		self.muond0 = self.tree["muon_trkd0"].array() 
		self.muonz0 = self.tree["muon_trkz0"].array()
		self.muonmass = self.tree["muon_m"].array()   
		self.muonistight = self.tree["muon_isTight"].array()
		self.muontype = self.tree["muon_type"].array()
		# self.muonpassPfilter = self.tree["muon_passesPromptCuts"].array()  #need this!! 1 per event :) 

		self.elpt = self.tree["el_pt"].array()  
