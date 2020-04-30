import ROOT
import uproot
import helpers

logger = helpers.getLogger('dHNLAnalysis.treenames')

class Tree():
	def __init__(self, fileName,treeName,vtx_container):
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
		self.evtNum = self.tree["eventNumber"].array()
		self.passedtriggers = self.tree["passedTriggers"].array()
		self.mumufilter = self.tree["passesHnlMuMuFilter"].array()
		self.muelfilter = self.tree["passesHnlMuElFilter"].array()
		self.elelfilter = self.tree["passesHnlElElFilter"].array()
		self.elmufilter = self.tree["passesHnlElMuFilter"].array()

		if "vertex_x" in self.file[treeName].keys(): 
			self.pvx = self.tree["vertex_x"].array()  
			self.pvy = self.tree["vertex_y"].array()  
			self.pvz = self.tree["vertex_z"].array() 
		
		# -----------------------------------------------------------------------
		# DV track variables
		# -----------------------------------------------------------------------
		self.trk_muonindex = self.tree[DVprefix + "_trk_muonIndex"].array()
		self.trk_elindex = self.tree[DVprefix + "_trk_electronIndex"].array()
		self.trackpt = self.tree[DVprefix + "_trk_pt_wrtSV"].array() 
		self.tracketa = self.tree[DVprefix + "_trk_eta_wrtSV"].array() 
		self.trackphi = self.tree[DVprefix + "_trk_phi_wrtSV"].array()
		self.trackmass = self.tree[DVprefix + "_trk_M"].array()
		self.trackcharge = self.tree[DVprefix + "_trk_charge"].array()
		self.trackchi2 = self.tree[DVprefix + "_trk_chi2_toSV"].array()
		
		self.trackd0 = self.tree[DVprefix + "_trk_d0"].array() #d0 and z0 are track quantites calculated wrt IP we dont want to look at the _wrtSV variable!
		self.trackz0 = self.tree[DVprefix + "_trk_z0"].array() 

		self.trackz0 = self.tree[DVprefix + "_trk_z0"].array() 
		self.trackd0 = self.tree[DVprefix + "_trk_d0"].array() 

		# -----------------------------------------------------------------------
		# DV reco variables
		# -----------------------------------------------------------------------
		self.dvx = self.tree[DVprefix + "_x"].array()  
		self.dvy = self.tree[DVprefix + "_y"].array()  
		self.dvz = self.tree[DVprefix + "_z"].array() 
		self.dvr = self.tree[DVprefix + "_r"].array() 
		self.dvmass = self.tree[DVprefix + "_mass"].array()  
		self.dvpt = self.tree[DVprefix + "_pt"].array()
		self.dveta = self.tree[DVprefix + "_eta"].array()  
		self.dvphi = self.tree[DVprefix + "_phi"].array() 
		self.dvminOpAng = self.tree[DVprefix + "_minOpAng"].array() 
		self.dvmaxOpAng = self.tree[DVprefix + "_maxOpAng"].array()   
		self.dvntrk = self.tree[DVprefix + "_ntrk"].array() 
		self.dvdistFromPV = self.tree[DVprefix + "_distFromPV"].array() 
		self.dvcharge = self.tree[DVprefix + "_charge"].array() 
		self.dvchi2 = self.tree[DVprefix + "_chi2"].array() 
		# -----------------------------------------------------------------------
		
		# -----------------------------------------------------------------------
		# muon variables 
		# -----------------------------------------------------------------------
		self.muonindex = self.tree["muon_index"].array()
		self.muonpt = self.tree["muon_pt"].array() 
		self.muoneta = self.tree["muon_eta"].array()  
 		self.muonphi = self.tree["muon_phi"].array() 
		self.muond0 = self.tree["muon_trkd0"].array() 
		self.muonz0 = self.tree["muon_trkz0"].array()
		self.muonpx = self.tree["muon_px"].array() 
		self.muonpy = self.tree["muon_py"].array() 
		self.muonpz = self.tree["muon_pz"].array() 
		self.muonmass = self.tree["muon_m"].array()  
		self.tightmu = self.tree["muon_isTight"].array()
		self.mediummu = self.tree["muon_isMedium"].array()
		self.loosemu = self.tree["muon_isLoose"].array()
		self.muontype = self.tree["muon_type"].array()
		self.muonpassPfilter = self.tree["muon_passesPromptCuts"].array() 
		self.muontrigmatched = self.tree["muon_isTrigMatched"].array()
		# -----------------------------------------------------------------------
		
		# -----------------------------------------------------------------------
		# electron variables
		# -----------------------------------------------------------------------
		self.elindex = self.tree["el_index"].array()
		self.elpt = self.tree["el_pt"].array()  
		self.eleta = self.tree["el_eta"].array()  
		self.elphi = self.tree["el_phi"].array()  
		self.eld0 = self.tree["el_trkd0"].array() 
		self.elz0 = self.tree["el_trkz0"].array()
		self.elmass = self.tree["el_m"].array()  
		self.tightel = self.tree["el_LHTight"].array()
		self.elpassPfilter = self.tree["el_passesPromptCuts"].array() 
		# -----------------------------------------------------------------------



		# -----------------------------------------------------------------------
		# DV truth variables
		# -----------------------------------------------------------------------
		if self.isData == False:
			self.truth_x = self.tree["truthVtx_x"].array()  
			self.truth_y = self.tree["truthVtx_y"].array()  
			self.truth_z = self.tree["truthVtx_z"].array() 
			self.truth_dvr = self.tree["truthVtx_r"].array() 
			self.truth_dvmass = self.tree["truthVtx_mass"].array()  
			self.truth_dvpt = self.tree["truthVtx_pt"].array()
			self.truth_dveta = self.tree["truthVtx_eta"].array()  
			self.truth_dvphi = self.tree["truthVtx_phi"].array() 
			self.truth_outP_pdgId = self.tree["truthVtx_outP_pdgId"].array()
			self.truth_outP_pt = self.tree["truthVtx_outP_pt"].array()
			self.truth_outP_eta = self.tree["truthVtx_outP_eta"].array()
			self.truth_outP_phi = self.tree["truthVtx_outP_phi"].array()
			self.truth_outP_m = self.tree["truthVtx_outP_M"].array()
			self.truth_parent_pdgId = self.tree["truthVtx_parent_pdgId"].array()
			self.truth_parent_pt = self.tree["truthVtx_parent_pt"].array()
			self.truth_parent_eta = self.tree["truthVtx_parent_eta"].array()
			self.truth_parent_phi = self.tree["truthVtx_parent_phi"].array()
			self.truth_parent_m = self.tree["truthVtx_parent_M"].array()
			# self.dvntrk = self.tree["truthVtx_ntrk"].array() 
			# self.dvdistFromPV = self.tree["truthVtx_distFromPV"].array() 
			# self.dvcharge = self.tree["truthVtx_charge"].array() 
			# self.dvcharge = self.tree["truthVtx_chi2"].array() 
		# -----------------------------------------------------------------------
		logger.info("Done importing trees!")
