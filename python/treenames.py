import ROOT
import uproot
import helpers

logger = helpers.getLogger('dHNLAnalysis.treenames')

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

		# -----------------------------------------------------------------------
		self.passedtriggers = self.tree["passedTriggers"].array(entrystop=nentries)
		self.mumufilter = self.tree["passesHnlMuMuFilter"].array(entrystop=nentries)
		self.muelfilter = self.tree["passesHnlMuElFilter"].array(entrystop=nentries)
		self.elelfilter = self.tree["passesHnlElElFilter"].array(entrystop=nentries)
		self.elmufilter = self.tree["passesHnlElMuFilter"].array(entrystop=nentries)

		if "vertex_x" in self.file[treeName].keys(): 
			self.pvx = self.tree["vertex_x"].array(entrystop=nentries)  
			self.pvy = self.tree["vertex_y"].array(entrystop=nentries)  
			self.pvz = self.tree["vertex_z"].array(entrystop=nentries) 
		
		# -----------------------------------------------------------------------
		# DV track variables
		# -----------------------------------------------------------------------
		self.trk_muonindex = self.tree[DVprefix + "_trk_muonIndex"].array(entrystop=nentries)
		self.trk_elindex = self.tree[DVprefix + "_trk_electronIndex"].array(entrystop=nentries)
		self.trackpt = self.tree[DVprefix + "_trk_pt_wrtSV"].array(entrystop=nentries) 
		self.tracketa = self.tree[DVprefix + "_trk_eta_wrtSV"].array(entrystop=nentries) 
		self.trackphi = self.tree[DVprefix + "_trk_phi_wrtSV"].array(entrystop=nentries)
		self.trackmass = self.tree[DVprefix + "_trk_M"].array(entrystop=nentries)
		self.trackd0 = self.tree[DVprefix + "_trk_d0_wrtSV"].array(entrystop=nentries) 
		self.trackz0 = self.tree[DVprefix + "_trk_z0_wrtSV"].array(entrystop=nentries) 
		self.trackcharge = self.tree[DVprefix + "_trk_charge"].array(entrystop=nentries)
		self.trackchi2 = self.tree[DVprefix + "_trk_chi2_toSV"].array(entrystop=nentries)

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
		self.elpassPfilter = self.tree["el_passesPromptCuts"].array(entrystop=nentries) 
		# -----------------------------------------------------------------------



		# -----------------------------------------------------------------------
		# DV truth variables
		# -----------------------------------------------------------------------
		if self.isData == False:
			self.truth_dvx = self.tree["truthVtx_x"].array(entrystop=nentries)  
			self.truth_dvy = self.tree["truthVtx_y"].array(entrystop=nentries)  
			self.truth_dvz = self.tree["truthVtx_z"].array(entrystop=nentries) 
			self.truth_dvr = self.tree["truthVtx_r"].array(entrystop=nentries) 
			self.truth_dvmass = self.tree["truthVtx_mass"].array(entrystop=nentries)  
			self.truth_dvpt = self.tree["truthVtx_pt"].array(entrystop=nentries)
			self.truth_dveta = self.tree["truthVtx_eta"].array(entrystop=nentries)  
			self.truth_dvphi = self.tree["truthVtx_phi"].array(entrystop=nentries) 
			# self.dvntrk = self.tree["truthVtx_ntrk"].array(entrystop=nentries) 
			# self.dvdistFromPV = self.tree["truthVtx_distFromPV"].array(entrystop=nentries) 
			# self.dvcharge = self.tree["truthVtx_charge"].array(entrystop=nentries) 
			# self.dvcharge = self.tree["truthVtx_chi2"].array(entrystop=nentries) 
		# -----------------------------------------------------------------------
		logger.info("Done importing trees!")
