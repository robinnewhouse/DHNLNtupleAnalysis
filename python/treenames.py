import ROOT
import uproot

class Tree():
	def __init__(self, fileName,treeName):
		self.file = uproot.open(fileName)
		self.tree = self.file[treeName]
		
		#define variables
		self.passedtriggers = self.tree["passedTriggers"].array()
		self.mumufilter = self.tree["passesHnlMuMuFilter"].array()
		self.muelfilter = self.tree["passesHnlMuElFilter"].array()
		self.elelfilter = self.tree["passesHnlElElFilter"].array()
		self.elmufilter = self.tree["passesHnlElMuFilter"].array()
		self.trk_muonindex = self.tree["secVtx_trk_muonIndex"].array()
		self.muonindex = self.tree["muon_index"].array()
		self.trk_elindex = self.tree["secVtx_trk_electronIndex"].array()
		self.elindex = self.tree["el_index"].array()
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
		self.muoneta = self.tree["muon_eta"].array()  
 		self.muonphi = self.tree["muon_phi"].array() 
		self.muond0 = self.tree["muon_trkd0"].array() 
		self.muonz0 = self.tree["muon_trkz0"].array()
		self.muonpx = self.tree["muon_px"].array() 
		self.muonpy = self.tree["muon_py"].array() 
		self.muonpz = self.tree["muon_pz"].array() 
		self.muonmass = self.tree["muon_m"].array() 
		self.muonmass = self.tree["muon_m"].array()   
		self.tightmu = self.tree["muon_isTight"].array()
		self.mediummu = self.tree["muon_isMedium"].array()
		self.loosemu = self.tree["muon_isLoose"].array()
		self.muontype = self.tree["muon_type"].array()
		self.muonpassPfilter = self.tree["muon_passesPromptCuts"].array() 
		self.tightel = self.tree["el_LHTight"].array()



		self.elpt = self.tree["el_pt"].array()  
