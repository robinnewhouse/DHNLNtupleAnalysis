import ROOT
import uproot


class Tree():
    def __init__(self, fileName, treeName):
        self.file = uproot.open(fileName)
        self.tree = self.file[treeName]

        # -----------------------------------------------------------------------
        try:
            self.passedtriggers = self.tree["passedTriggers"].array()
            self.mumufilter = self.tree["passesHnlMuMuFilter"].array()
            self.muelfilter = self.tree["passesHnlMuElFilter"].array()
            self.elelfilter = self.tree["passesHnlElElFilter"].array()
            self.elmufilter = self.tree["passesHnlElMuFilter"].array()
        except:
            print "Trigger/Filter information not available in tree"

        # -----------------------------------------------------------------------
        # DV track variables
        # -----------------------------------------------------------------------
        try:
            self.trk_muonindex = self.tree["secVtx_trk_muonIndex"].array()
            self.trk_elindex = self.tree["secVtx_trk_electronIndex"].array()
            self.trackpt = self.tree["secVtx_trk_pt_wrtSV"].array()
            self.tracketa = self.tree["secVtx_trk_eta_wrtSV"].array()
            self.trackphi = self.tree["secVtx_trk_phi_wrtSV"].array()
            self.trackmass = self.tree["secVtx_trk_M"].array()
            self.trackd0 = self.tree["secVtx_trk_d0_wrtSV"].array()
            self.trackz0 = self.tree["secVtx_trk_z0_wrtSV"].array()
            self.trackcharge = self.tree["secVtx_trk_charge"].array()
            self.trackchi2 = self.tree["secVtx_trk_chi2_toSV"].array()
        except:
            print "Secondary vertex track information not available in tree"

        # -----------------------------------------------------------------------
        # DV reco variables
        # -----------------------------------------------------------------------
        try:
            self.dvx = self.tree["secVtx_x"].array()
            self.dvy = self.tree["secVtx_y"].array()
            self.dvz = self.tree["secVtx_z"].array()
            self.dvr = self.tree["secVtx_r"].array()
            self.dvmass = self.tree["secVtx_mass"].array()
            self.dvpt = self.tree["secVtx_pt"].array()
            self.dveta = self.tree["secVtx_eta"].array()
            self.dvphi = self.tree["secVtx_phi"].array()
            self.dvminOpAng = self.tree["secVtx_minOpAng"].array()
            self.dvmaxOpAng = self.tree["secVtx_maxOpAng"].array()
            self.dvntrk = self.tree["secVtx_ntrk"].array()
            self.dvdistFromPV = self.tree["secVtx_distFromPV"].array()
            self.dvcharge = self.tree["secVtx_charge"].array()
            self.dvchi2 = self.tree["secVtx_chi2"].array()
        except:
            print "Secondary vertex information not available in tree"
        # -----------------------------------------------------------------------

        # -----------------------------------------------------------------------
        # muon variables
        # -----------------------------------------------------------------------
        try:
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
        except:
            print "Muon information not available in tree"
        # -----------------------------------------------------------------------

        # -----------------------------------------------------------------------
        # electron variables
        # -----------------------------------------------------------------------
        try:
            self.elindex = self.tree["el_index"].array()
            self.elpt = self.tree["el_pt"].array()
            self.eleta = self.tree["el_eta"].array()
            self.elphi = self.tree["el_phi"].array()
            self.eld0 = self.tree["el_trkd0"].array()
            self.elz0 = self.tree["el_trkz0"].array()
            self.elmass = self.tree["el_m"].array()
            self.tightel = self.tree["el_LHTight"].array()
            self.elpassPfilter = self.tree["el_passesPromptCuts"].array()
        except:
            print "Muon information not available in tree"
        # -----------------------------------------------------------------------

        # -----------------------------------------------------------------------
        # DV truth variables
        # -----------------------------------------------------------------------
        try:
            self.truth_dvx = self.tree["truthVtx_x"].array()
            self.truth_dvy = self.tree["truthVtx_y"].array()
            self.truth_dvz = self.tree["truthVtx_z"].array()
            self.truth_dvr = self.tree["truthVtx_r"].array()
            self.truth_dvmass = self.tree["truthVtx_mass"].array()
            self.truth_dvpt = self.tree["truthVtx_pt"].array()
            self.truth_dveta = self.tree["truthVtx_eta"].array()
            self.truth_dvphi = self.tree["truthVtx_phi"].array()
        # self.dvntrk = self.tree["truthVtx_ntrk"].array()
        # self.dvdistFromPV = self.tree["truthVtx_distFromPV"].array()
        # self.dvcharge = self.tree["truthVtx_charge"].array()
        # self.dvcharge = self.tree["truthVtx_chi2"].array()
        except:
            print "Truth information not available in tree"
        # -----------------------------------------------------------------------
