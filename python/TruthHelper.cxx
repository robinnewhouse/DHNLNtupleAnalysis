#include "NtupleAnalysisUtils/Ntuple/NTAUNtupleIncludes.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "DHNLTree.h"

struct TruthOutput{

	TLorentzVector HNL_vec;
	TLorentzVector dNu_vec;
	std::vector<TLorentzVector> trkVec={}; 
	std::vector<TLorentzVector> dLepVec={}; 
	std::vector<int> dLep_pdgID={}; 
	std::vector<int> dLepCharge={}; 
	std::vector<TLorentzVector> dEl={}; 
	std::vector<int> dEl_charge={}; 
	std::vector<float> dEl_d0={}; 
	std::vector<TLorentzVector> dMu={}; 
	std::vector<int> dMu_charge={}; 
	std::vector<float> dMu_d0={}; 
	std::vector<float> dTrk_d0={}; 
	double truth_dvx = -1;
	double truth_dvy = -1;
	double truth_dvz = -1;
	TVector3 truth_dv;
	double truth_dvr = -1;
	double truth_pvx = -1;
	double truth_pvy = -1;
	double truth_pvz = -1;
	TVector3 truth_pv;
	TLorentzVector W_vec;
	int W_charge = -2;
	TLorentzVector plep_vec;
	int plep_charge = -99;
	double mhnl = -1;
	double dvmass = -1;
	int HNL_pdgID = 50;
	double gamma = 1;
	double beta = 1;
	double properLifetime = -1;

};

TruthOutput fillTruth(DHNLTree &t){
    TruthOutput self{}; 
    int nVx = t.ntruthVtx(); 
    bool have_DV = false;  
    bool have_PV = false; 
    for (size_t ivx = 0; ivx < nVx; ++ivx){
        if (have_DV and have_PV) break;
        // get the DV!
        int vx_pdg = std::abs(t.truthVtx_parent_pdgId(ivx));
        int vx_nChildren = t.truthVtx_nOutP(ivx);
        if (vx_pdg == 50 && vx_nChildren == 3){
            self.truth_dvx = t.truthVtx_x(ivx);
            self.truth_dvy = t.truthVtx_y(ivx);
            self.truth_dvz = t.truthVtx_z(ivx);
            // # for proper lifetime calculation
            self.gamma = t.truthVtx_parent_E(ivx) / t.truthVtx_parent_M(ivx); 
            self.beta = std::sqrt(1. - 1. / self.gamma /self.gamma);
            self.truth_dv = TVector3(self.truth_dvx, self.truth_dvy, self.truth_dvz);
            self.truth_dvr = std::sqrt(self.truth_dvx *self.truth_dvx + self.truth_dvy *self.truth_dvy);
            for (size_t i = 0; i < vx_nChildren; ++i){
                int trk_pdgId = std::abs(t.truthVtx_outP_pdgId(ivx)[i]); 
                TLorentzVector p4_out;
                p4_out.SetPtEtaPhiM(t.truthVtx_outP_pt(ivx)[i],
                                        t.truthVtx_outP_eta(ivx)[i],
                                        t.truthVtx_outP_phi(ivx)[i],
                                        t.truthVtx_outP_M(ivx)[i]
                                        );
                int q_out =  t.truthVtx_outP_charge(ivx)[i];
                float d0_out = t.truthVtx_outP_d0(ivx)[i];
                if (trk_pdgId == 13){
                    self.dMu.push_back(p4_out);
                    self.dMu_charge.push_back(q_out);
                    self.dMu_d0.push_back(d0_out);
                }
                if (trk_pdgId == 11){
                    self.dEl.push_back(p4_out);
                    self.dEl_charge.push_back(q_out);
                    self.dEl_d0.push_back(d0_out);
                }

                if (trk_pdgId == 13 or trk_pdgId == 11){
                    //   # is track a muon of electron? Then these are our visible (charged) truth tracks
                    self.trkVec.push_back(p4_out);  
                    // # only add visible leptons to trkVec list
                }
                else{//  # remaining child is the neutrino
                    self.dNu_vec = p4_out;
                }
                self.dLep_pdgID.push_back(trk_pdgId);
                self.dLepVec.push_back(p4_out);//   # add all the displaced leptons to one list in the order they are in pythia
                self.dLepCharge.push_back(q_out);
                // # self.dTrk_d0.push_back(d0_out);
                self.dTrk_d0.push_back(-1); //  # fill with -1 for now, default DHNLalg does not have truth d0
            }
            self.HNL_vec.SetPtEtaPhiM(t.truthVtx_parent_pt(ivx),
                                        t.truthVtx_parent_eta(ivx),
                                        t.truthVtx_parent_phi(ivx),
                                        t.truthVtx_parent_M(ivx)
                                        );
            have_DV = true;
        }

        // # get the primary vertex
        else if (vx_pdg == 24 && vx_nChildren == 2){
                // # TODO: Should we be checking if one of the children is an HNL?
            self.truth_pvx = t.truthVtx_x(ivx);
            self.truth_pvy = t.truthVtx_y(ivx);
            self.truth_pvz = t.truthVtx_z(ivx);
            self.truth_pv = TVector3(self.truth_pvx, self.truth_pvy, self.truth_pvz);

            self.plep_vec.SetPtEtaPhiM(t.truthVtx_outP_pt(ivx)[0],
                                        t.truthVtx_outP_eta(ivx)[0],
                                        t.truthVtx_outP_phi(ivx)[0],
                                        t.truthVtx_outP_M(ivx)[0]
                                        );
            self.plep_charge = t.truthVtx_outP_charge(ivx)[0];
            self.W_vec.SetPtEtaPhiM(t.truthVtx_parent_pt(ivx),
                                    t.truthVtx_parent_eta(ivx),
                                    t.truthVtx_parent_phi(ivx),
                                    t.truthVtx_parent_M(ivx)
                                    );
            self.W_charge = t.truthVtx_parent_charge(ivx);
            have_PV = true;
        }
    }

    // # calculate proper lifetime
    double dx = std::abs(self.truth_pvx - self.truth_dvx);
    double dy = std::abs(self.truth_pvy - self.truth_dvy);
    double dz = std::abs(self.truth_pvz - self.truth_dvz);
    double dr = std::sqrt(dx*dx + dy*dy + dz*dz);
    self.properLifetime = dr/(self.gamma * self.beta);
    return self;
}