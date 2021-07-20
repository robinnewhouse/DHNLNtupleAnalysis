import os, math, ROOT, json,sys
import numpy as np
from ROOT import *
from pylab import *
import atlas_style
ROOT.gROOT.SetBatch(True)
sys.path.append('/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/python/')
sys.path.append('/home/dtrischuk/HNLAnalysis/DHNLPlotting/util/')
from plot_classes import Hist1D,Hist1DRatio,Hist2D
import helpers, selections, ntuples
gROOT.SetStyle("ATLAS")


def main():
    # update all the flags from arg parser
    makeHist = options.makeHist
    plot = options.plot
    looser_plep = options.looser_plep
    use_loose_DVs = options.use_loose_DVs
    use_CR_for_DVs = options.use_CR_for_DVs
    mix_channels = options.mix_channels
    shuffle_PV = options.shuffle_PV
    CR_charge = options.CR_charge
    signal = options.signal
    version = options.version

    # Get SS background files
    if mix_channels:
        if signal == "uuu": f_SS_bkg =  ROOT.TFile('/data/hnl/v6_histograms/jun15_SSbkg_v6p2_histograms/fullrun2_histograms_SSbkg_uue.root')
        if signal == "uee": f_SS_bkg =  ROOT.TFile('/data/hnl/v6_histograms/jun15_SSbkg_v6p2_histograms/fullrun2_histograms_SSbkg_uue.root')
        if signal == "uue": f_SS_bkg = ROOT.TFile('/data/hnl/v6_histograms/jun15_SSbkg_v6p2_histograms/fullrun2_histograms_SSbkg_uee.root')
        if signal == "eeu": f_SS_bkg = ROOT.TFile('/data/hnl/v6_histograms/jun15_SSbkg_v6p2_histograms/fullrun2_histograms_SSbkg_eee.root')
        if signal == "euu": f_SS_bkg = ROOT.TFile('/data/hnl/v6_histograms/jun15_SSbkg_v6p2_histograms/fullrun2_histograms_SSbkg_eee.root')
        if signal == "eee": f_SS_bkg = ROOT.TFile('/data/hnl/v6_histograms/jun15_SSbkg_v6p2_histograms/fullrun2_histograms_SSbkg_eeu.root')
    else: 
        f_SS_bkg =  ROOT.TFile(f'/data/hnl/{version}_histograms/jun15_SSbkg_v6p2_histograms/fullrun2_histograms_SSbkg_{signal}.root')

    # Get inverted prompt lepton VR files used to select DVs
    if signal == "uuu" or signal == "euu":
        f_CR = ROOT.TFile(f'/data/hnl/v6_histograms/jun15_CR_v6p2_histograms/fullrun2_CR_histograms_data_{CR_charge}_uu.root')
    if signal == "uue" or signal == "eeu":
        f_CR = ROOT.TFile(f'/data/hnl/v6_histograms/jun15_CR_v6p2_histograms/fullrun2_CR_histograms_data_{CR_charge}_eu.root')
    if signal == "eee" or signal == "uee":
        f_CR = ROOT.TFile(f'/data/hnl/v6_histograms/jun15_CR_v6p2_histograms/fullrun2_CR_histograms_data_{CR_charge}_ee.root')

    # make output dir and output mini-tree root file
    if use_loose_DVs:  outputDir = f'../shuffled_bkg_output/shuffled_background_estimate_looseDVs_{version}_{signal}/'
    elif use_CR_for_DVs: outputDir =  f'../shuffled_bkg_output/shuffled_background_estimate_{version}_{CR_charge}_{signal}/'
    else: outputDir = f'../shuffled_bkg_output/shuffled_background_estimate_{version}_{signal}/'

    if use_CR_for_DVs:
        if CR_charge == "OS": 
            ouputf_name = f"shuffled_bkg_ntuple_nominal_{signal}.root"
        elif CR_charge == "SS": 
            ouputf_name = f"shuffled_bkg_ntuple_DV_systematic_{signal}.root"
        else: 
            ouputf_name = f"shuffled_bkg_ntuple_{signal}.root"
    elif mix_channels:
        ouputf_name = f"shuffled_bkg_ntuple_plep_systematic_{signal}.root"
    else: 
        ouputf_name = f"shuffled_bkg_ntuple_{signal}.root"
    if not os.path.exists(outputDir): os.mkdir(outputDir)
    if not os.path.exists(outputDir + "plots/"): os.mkdir(outputDir + "plots/")

    # set default vertexing algorithm and selection to get mini-trees from DHNL Ntuple Analysis output
    vtx_alg = "VSI_LeptonsMod"
    selection = "DVtype"

    # Get the mini-tree used to select DVs from the DHNLNtupleAnalysis output files
    if use_CR_for_DVs: tree = f_CR.Get('{}_ntuples_{}'.format(vtx_alg, selection))  # get TTree
    else: tree = f_SS_bkg.Get('{}_ntuples_{}'.format(vtx_alg, selection))  # get TTree

    # Get the mini-tree used to select prompt leptons from the DHNLNtupleAnalysis output files
    tree_shuffle = f_SS_bkg.Get('{}_ntuples_{}'.format(vtx_alg, selection))  # get TTree

    # Get the number of entries in the respective trees
    dv_nentries = tree.GetEntries()
    plep_nentries = tree_shuffle.GetEntries()

    output_text = f"""
    total entries for DV: {dv_nentries}
    total entries for plep: {plep_nentries}
    """
    print (f'total entries for DV: {dv_nentries}')
    print (f'total entries for plep: {plep_nentries}')
    offset = 1

    # Define histograms
    h ={}
    h["mDV"] = ROOT.TH1D('h_mDV', 'h_mDV', 50, 0, 50)
    h["DVpt"] = ROOT.TH1D('h_DVpt', 'h_DVpt', 25, 0, 100)
    h["DVeta"] = ROOT.TH1D('h_DVeta', 'h_DVeta', 100, -3, 3)
    h["DVphi"] = ROOT.TH1D('h_DVphi', 'h_DVphi', 100, -4, 4)
    h["plep_pt"] = ROOT.TH1D('h_plep_pt', 'h_plep_pt', 50, 0, 200)
    h["plep_eta"] = ROOT.TH1D('h_plep_eta', 'h_plep_eta', 100, -3, 3)
    h["plep_phi"] = ROOT.TH1D('h_plep_phi', 'h_plep_phi', 100, -4, 4)
    h["Lxy"] = ROOT.TH1D('h_Lxy', 'h_Lxy', 30, 0, 300)
    h["HNLm_orig"] = ROOT.TH1D('h_mHNL', 'h_mHNL', 25, 0, 50)
    h["mvis_orig"] = ROOT.TH1D('h_mvis_orig', 'h_mvis_orig', 25, 0, 200)
    h["mll_dEl_plep_orig"] = ROOT.TH1D('h_mll_dEl_plep_orig', 'h_mll_dEl_plep_orig', 100, 0, 200)

    # for shuffled branch (matching name convensions in original ntuple code.)
    h["HNLm"] = ROOT.TH1D('h_HNLm', 'h_HNLm', 25, 0, 50)
    h["mvis"]  = ROOT.TH1D('h_mvis', 'h_mvis', 25, 0, 200)
    h["mll_1"] = ROOT.TH1D('h_mll_1', 'h_mll_1', 100, 0, 200)
    h["mll_0"] = ROOT.TH1D('h_mll_0', 'h_mll_0', 100, 0, 200)
    h["DV_weight_LNC_only"] = ROOT.TH1D('weight_LNC_only', 'weight_LNC_only', 100, 0, 200)
    h["DV_weight_LNC_plus_LNV"] = ROOT.TH1D('weight_LNC_plus_LNV', 'weight_LNC_plus_LNV', 100, 0, 200)
    h["DV_mass"] = ROOT.TH1D('h_DV_mass','h_DV_mass', 100,0,50 ) 
    h["DV_r"] = ROOT.TH1D('h_DV_r','h_DV_r', 30,0,300)
    h["DV_index"] = ROOT.TH1D('h_DV_index','h_DV_index', 500,0,500)
    h["plep_index"] = ROOT.TH1D('h_plep_index','h_plep_index', 500,0,500)
    h["DV_z"] = ROOT.TH1D('h_DV_z','h_DV_z', 30,0,300)
    h["DV_charge"] = ROOT.TH1D('h_DV_charge','h_DV_charge', 2, -0.5, 1.5 ) 
    h["DV_2medium"] = ROOT.TH1D('h_2med','h_2med', 2, -0.5, 1.5 ) 
    h["DV_2loose"] = ROOT.TH1D('h_2loose','h_2loose', 2, -0.5, 1.5) 
    h["DV_2veryveryloose"] = ROOT.TH1D('h_2veryveryloose','h_2veryveryloose', 2, -0.5, 1.5) 
    h["DV_medium_veryveryloose"] = ROOT.TH1D('h_0','h_0',2, -0.5, 1.5)
    h["DV_alpha"] = ROOT.TH1D('h_med_loose','h_med_loose',10,0,3)
    h["DV_cosmic_sep"] = ROOT.TH1D('h_sep','h_sep', 100,0,1 )
    h["DV_pass_lep_pt"] = ROOT.TH1D('h_pass_lep_pt','h_pass_lep_pt',2, -0.5, 1.5 )
    h["mll_dMu_plep_is_OS"] = ROOT.TH1D('h_mu_is_OS','h_mu_is_OS',2, -0.5, 1.5 )
    h["mll_dMu_plep_is_SS"] = ROOT.TH1D('h_mu_is_SS','h_mu_is_SS',2, -0.5, 1.5 )
    h["mll_dMu_plep"] = ROOT.TH1D('h_mll_dMu','h_mll_dMu', 100,0,200)
    h["mll_dEl_plep_is_OS"] = ROOT.TH1D('h_el_is_OS','h_el_is_OS',2, -0.5, 1.5 )
    h["mll_dEl_plep_is_SS"] = ROOT.TH1D('h_el_is_SS','h_el_is_SS',2, -0.5, 1.5 )
    h["mll_dEl_plep"] = ROOT.TH1D('h_mll_dEl','h_mll_dEl',100,0,200 )
    h["DV_pass_mat_veto"] = ROOT.TH1D('h_mat_veto','h_mat_veto',2, -0.5, 1.5)

    micro_ntuples ={}

    def fill_ntuple(selection, ntuple_name, variable, full_name="", p=None):
        """
        A helper function for filling micro-ntuples. Often called from the fill_hist function.
        If you are using this in you analysis,
        please check that it is not also being called by fill_hist to prevent double-counting.
        :param selection: the step of selection the analysis it at. May be "None" in which case there will be no prefix.
        :param ntuple_name: base name of the ntuple. When saved, a prefix and suffix will be appended.
        :param variable: variable you want to fill the histogram with.
        :param full_name: override the automatic naming of the ntuple.
        """
        if not selection:
            raise ValueError("You must indicate a selection in order to store the ntuple. Use 'all' if no selection.")

        if selection not in micro_ntuples:
            micro_ntuples[selection] = ntuples.Ntuples('ntuples_{}_{}'.format(selection, vtx_alg))  # temp name. not written
        # The name of the ntuple
        if not full_name:
            full_name = ntuple_name
        micro_ntuples[selection][full_name] = variable

    def fill_hist(sel,name, variable ,weight=1): 
        h[name].Fill(variable,weight) 
        fill_ntuple(sel, name, variable)


    def pass_dv_cuts(channel, tree): 
        if mix_channels:
            if channel == "euu": pass_extra_cuts = tree.n_trigger_matched_medium > 0 and (((tree.DV_mass > 2 and tree.DV_mass < 5.5) and tree.DV_mass > -7/150*tree.DV_r + 7 ) or tree.DV_mass > 5.5 ) and tree.DV_2veryveryloose == 1 and tree.DV_pass_mat_veto == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dEl_plep_is_OS==1 and (tree.mll_dEl_plep<80 or tree.mll_dEl_plep>100)) or tree.mll_dEl_plep_is_SS==1 )
            if channel == "uee": pass_extra_cuts = tree.n_trigger_matched_medium > 0 and (((tree.DV_mass > 2 and tree.DV_mass < 5.5) and tree.DV_mass > -7/150*tree.DV_r + 7 ) or tree.DV_mass > 5.5 ) and tree.DV_medium_veryveryloose == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dMu_plep_is_OS==1 and (tree.mll_dMu_plep<80 or tree.mll_dMu_plep>100)) or tree.mll_dMu_plep_is_SS==1 )
            if channel == "uue": pass_extra_cuts = tree.n_trigger_matched_medium > 0 and (((tree.DV_mass > 2 and tree.DV_mass < 5.5) and tree.DV_mass > -7/150*tree.DV_r + 7 ) or tree.DV_mass > 5.5 ) and tree.DV_2veryveryloose == 1 and tree.DV_pass_mat_veto == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05
            if channel == "eeu": pass_extra_cuts = tree.n_trigger_matched_medium > 0 and (((tree.DV_mass > 2 and tree.DV_mass < 5.5) and tree.DV_mass > -7/150*tree.DV_r + 7 ) or tree.DV_mass > 5.5 ) and tree.DV_2veryveryloose == 1 and tree.DV_pass_mat_veto == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dEl_plep_is_OS==1 and (tree.mll_dEl_plep<80 or tree.mll_dEl_plep>100)) or tree.mll_dEl_plep_is_SS==1 )
            if channel == "eee": pass_extra_cuts = tree.n_trigger_matched_medium > 0 and (((tree.DV_mass > 2 and tree.DV_mass < 5.5) and tree.DV_mass > -7/150*tree.DV_r + 7 ) or tree.DV_mass > 5.5 ) and tree.DV_medium_veryveryloose == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dEl_plep_is_OS==1 and (tree.mll_dEl_plep<80 or tree.mll_dEl_plep>100)) or tree.mll_dEl_plep_is_SS==1 ) 
            if channel == "uuu": pass_extra_cuts = tree.n_trigger_matched_medium > 0 and (((tree.DV_mass > 2 and tree.DV_mass < 5.5) and tree.DV_mass > -7/150*tree.DV_r + 7 ) or tree.DV_mass > 5.5 ) and tree.DV_medium_veryveryloose == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dMu_plep_is_OS==1 and (tree.mll_dMu_plep<80 or tree.mll_dMu_plep>100)) or tree.mll_dMu_plep_is_SS==1 )
        else: 
            if channel == "euu": pass_extra_cuts = tree.n_trigger_matched_medium > 0 and tree.DV_mass > 5.5 and tree.DV_2medium== 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05
            if channel == "uee": pass_extra_cuts = tree.n_trigger_matched_medium > 0 and (((tree.DV_mass > 2 and tree.DV_mass < 5.5) and tree.DV_mass > -7/150*tree.DV_r + 7 ) or tree.DV_mass > 5.5 ) and tree.DV_2veryveryloose == 1 and tree.DV_pass_mat_veto == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05
            if channel == "uuu": pass_extra_cuts = tree.n_trigger_matched_medium > 0 and tree.DV_mass > 5.5 and tree.DV_2medium== 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dMu_plep_is_OS==1 and (tree.mll_dMu_plep<80 or tree.mll_dMu_plep>100)) or tree.mll_dMu_plep_is_SS==1)
            if channel == "eee": pass_extra_cuts = tree.n_trigger_matched_medium > 0 and (((tree.DV_mass > 2 and tree.DV_mass < 5.5) and tree.DV_mass > -7/150*tree.DV_r + 7 ) or tree.DV_mass > 5.5 ) and tree.DV_2veryveryloose == 1 and tree.DV_pass_mat_veto == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dEl_plep_is_OS==1 and (tree.mll_dEl_plep<80 or tree.mll_dEl_plep>100)) or tree.mll_dEl_plep_is_SS==1 )
            if channel == "eeu": pass_extra_cuts = tree.n_trigger_matched_medium > 0 and (((tree.DV_mass > 2 and tree.DV_mass < 5.5) and tree.DV_mass > -7/150*tree.DV_r + 7 ) or tree.DV_mass > 5.5 ) and tree.DV_medium_veryveryloose == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dEl_plep_is_OS==1 and (tree.mll_dEl_plep<80 or tree.mll_dEl_plep>100)) or tree.mll_dEl_plep_is_SS==1 ) 
            if channel == "uue": pass_extra_cuts = tree.n_trigger_matched_medium > 0 and (((tree.DV_mass > 2 and tree.DV_mass < 5.5) and tree.DV_mass > -7/150*tree.DV_r + 7 ) or tree.DV_mass > 5.5 ) and tree.DV_medium_veryveryloose == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dMu_plep_is_OS==1 and (tree.mll_dMu_plep<80 or tree.mll_dMu_plep>100)) or tree.mll_dMu_plep_is_SS==1 )
        return pass_extra_cuts

    def pass_loose_dv_cuts(channel, tree): 
        if channel == "euu": pass_extra_cuts = tree.DV_mass > 2 and tree.DV_2medium== 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 
        if channel == "uee": pass_extra_cuts = tree.DV_mass > 2 and tree.DV_2veryveryloose == 1 and tree.DV_pass_mat_veto == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05
        if channel == "uuu": pass_extra_cuts = tree.DV_mass > 2 and tree.DV_2medium== 1 and tree.DV_cosmic_sep > 0.05 and tree.DV_pass_lep_pt == 1  and ((tree.mll_dMu_plep_is_OS==1 and (tree.mll_dMu_plep<80 or tree.mll_dMu_plep>100)) or tree.mll_dMu_plep_is_SS==1)
        if channel == "eee": pass_extra_cuts = tree.DV_mass > 2 and tree.DV_2veryveryloose == 1 and tree.DV_pass_mat_veto == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dEl_plep_is_OS==1 and (tree.mll_dEl_plep<80 or tree.mll_dEl_plep>100)) or tree.mll_dEl_plep_is_SS==1 )
        if channel == "eeu": pass_extra_cuts = tree.DV_mass > 2 and tree.DV_medium_veryveryloose == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dEl_plep_is_OS==1 and (tree.mll_dEl_plep<80 or tree.mll_dEl_plep>100)) or tree.mll_dEl_plep_is_SS==1 ) 
        if channel == "uue": pass_extra_cuts = tree.DV_mass > 2 and tree.DV_medium_veryveryloose == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dMu_plep_is_OS==1 and (tree.mll_dMu_plep<80 or tree.mll_dMu_plep>100)) or tree.mll_dMu_plep_is_SS==1 )
        return pass_extra_cuts

    def pass_dv_cuts_SR(channel, tree): 
        if channel == "euu": pass_extra_cuts  = tree.DV_mass > 5.5 and tree.DV_2medium== 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and tree.HNLm < 50
        if channel == "uee": pass_extra_cuts = (((tree.DV_mass > 2 and tree.DV_mass < 5.5) and tree.DV_mass > -7/150*tree.DV_r + 7 ) or tree.DV_mass > 5.5 ) and tree.DV_2veryveryloose == 1 and tree.DV_pass_mat_veto == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and tree.HNLm < 50
        if channel == "uuu": pass_extra_cuts  = tree.DV_mass > 5.5 and tree.DV_2medium== 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dMu_plep_is_OS==1 and (tree.mll_dMu_plep<80 or tree.mll_dMu_plep>100)) or tree.mll_dMu_plep_is_SS == 1) and tree.HNLm < 50
        if channel == "eee": pass_extra_cuts = (((tree.DV_mass > 2 and tree.DV_mass < 5.5) and tree.DV_mass > -7/150*tree.DV_r + 7 ) or tree.DV_mass > 5.5 ) and tree.DV_2veryveryloose == 1 and tree.DV_pass_mat_veto == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dEl_plep_is_OS==1 and (tree.mll_dEl_plep<80 or tree.mll_dEl_plep>100)) or tree.mll_dEl_plep_is_SS==1 ) and tree.HNLm < 50
        if channel == "eeu": pass_extra_cuts = (((tree.DV_mass > 2 and tree.DV_mass < 5.5) and tree.DV_mass > -7/150*tree.DV_r + 7 ) or tree.DV_mass > 5.5 ) and tree.DV_medium_veryveryloose == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dEl_plep_is_OS==1 and (tree.mll_dEl_plep<80 or tree.mll_dEl_plep>100)) or tree.mll_dEl_plep_is_SS==1 ) and tree.HNLm < 50 
        if channel == "uue": pass_extra_cuts = (((tree.DV_mass > 2 and tree.DV_mass < 5.5) and tree.DV_mass > -7/150*tree.DV_r + 7 ) or tree.DV_mass > 5.5 ) and tree.DV_medium_veryveryloose == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dMu_plep_is_OS==1 and (tree.mll_dMu_plep<80 or tree.mll_dMu_plep>100)) or tree.mll_dMu_plep_is_SS==1 ) and tree.HNLm < 50
        return pass_extra_cuts

    def pass_dv_CR_cuts(channel, tree): 
        if channel == "euu": pass_extra_cuts  = tree.DV_2medium == 1 and tree.DV_mass > 5.5 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05
        if channel == "uee": pass_extra_cuts  = (((tree.DV_mass > 2 and tree.DV_mass < 5.5) and tree.DV_mass > -7/150*tree.DV_r + 7 ) or tree.DV_mass > 5.5 ) and tree.DV_2veryveryloose == 1 and tree.DV_pass_mat_veto == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05
        if channel == "uuu": pass_extra_cuts  = tree.DV_2medium == 1 and tree.DV_mass > 5.5 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05   
        if channel == "eee": pass_extra_cuts  = (((tree.DV_mass > 2 and tree.DV_mass < 5.5) and tree.DV_mass > -7/150*tree.DV_r + 7 ) or tree.DV_mass > 5.5 ) and tree.DV_2veryveryloose == 1 and tree.DV_pass_mat_veto == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05
        if channel == "eeu": pass_extra_cuts  = (((tree.DV_mass > 2 and tree.DV_mass < 5.5) and tree.DV_mass > -7/150*tree.DV_r + 7 ) or tree.DV_mass > 5.5 ) and tree.DV_medium_veryveryloose == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05
        if channel == "uue": pass_extra_cuts  = (((tree.DV_mass > 2 and tree.DV_mass < 5.5) and tree.DV_mass > -7/150*tree.DV_r + 7 ) or tree.DV_mass > 5.5 ) and tree.DV_medium_veryveryloose == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05
        return pass_extra_cuts

    def pass_looser_prompt_cuts(channel, tree): 
        if mix_channels:
            if channel == "euu": pass_extra_cuts =  tree.n_trigger_matched_medium > 0 and tree.DV_mass > 1 and tree.DV_pass_mat_veto == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dEl_plep_is_OS==1 and (tree.mll_dEl_plep<80 or tree.mll_dEl_plep>100)) or tree.mll_dEl_plep_is_SS==1 )
            if channel == "uee": pass_extra_cuts =  tree.n_trigger_matched_medium > 0 and tree.DV_mass > 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dMu_plep_is_OS==1 and (tree.mll_dMu_plep<80 or tree.mll_dMu_plep>100)) or tree.mll_dMu_plep_is_SS==1 )
            if channel == "uue": pass_extra_cuts =  tree.n_trigger_matched_medium > 0 and tree.DV_mass > 1 and tree.DV_pass_mat_veto == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05
            if channel == "eeu": pass_extra_cuts =  tree.n_trigger_matched_medium > 0 and tree.DV_mass > 1 and tree.DV_pass_mat_veto == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dEl_plep_is_OS==1 and (tree.mll_dEl_plep<80 or tree.mll_dEl_plep>100)) or tree.mll_dEl_plep_is_SS==1 )
            if channel == "eee": pass_extra_cuts =  tree.n_trigger_matched_medium > 0 and tree.DV_mass > 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dEl_plep_is_OS==1 and (tree.mll_dEl_plep<80 or tree.mll_dEl_plep>100)) or tree.mll_dEl_plep_is_SS==1 ) 
            if channel == "uuu": pass_extra_cuts =  tree.n_trigger_matched_medium > 0 and tree.DV_mass > 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dMu_plep_is_OS==1 and (tree.mll_dMu_plep<80 or tree.mll_dMu_plep>100)) or tree.mll_dMu_plep_is_SS==1 )
        else:
            if channel == "euu": pass_extra_cuts =  tree.n_trigger_matched_medium > 0 and tree.DV_mass > 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05
            if channel == "uee": pass_extra_cuts =  tree.n_trigger_matched_medium > 0 and tree.DV_mass > 1 and tree.DV_pass_mat_veto == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05
            if channel == "uuu": pass_extra_cuts =  tree.n_trigger_matched_medium > 0 and tree.DV_mass > 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dMu_plep_is_OS==1 and (tree.mll_dMu_plep<80 or tree.mll_dMu_plep>100)) or tree.mll_dMu_plep_is_SS == 1)
            if channel == "eee": pass_extra_cuts =  tree.n_trigger_matched_medium > 0 and tree.DV_mass > 1 and tree.DV_pass_mat_veto == 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dEl_plep_is_OS==1 and (tree.mll_dEl_plep<80 or tree.mll_dEl_plep>100)) or tree.mll_dEl_plep_is_SS==1 )
            if channel == "eeu": pass_extra_cuts =  tree.n_trigger_matched_medium > 0 and tree.DV_mass > 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dEl_plep_is_OS==1 and (tree.mll_dEl_plep<80 or tree.mll_dEl_plep>100)) or tree.mll_dEl_plep_is_SS==1 )
            if channel == "uue": pass_extra_cuts =  tree.n_trigger_matched_medium > 0 and tree.DV_mass > 1 and tree.DV_pass_lep_pt == 1 and tree.DV_cosmic_sep > 0.05 and ((tree.mll_dMu_plep_is_OS==1 and (tree.mll_dMu_plep<80 or tree.mll_dMu_plep>100)) or tree.mll_dMu_plep_is_SS==1 )
        
        return pass_extra_cuts

    if makeHist:
        file = ROOT.TFile.Open(outputDir+"/plots/" +  ouputf_name, 'recreate') 
        num_dv = 0.0
        num_plep = 0.0
        num_loose_dv = 0.0
        num_SR_dv= 0.0
        
        for i in range(dv_nentries):
            tree.GetEntry(i)
            if use_loose_DVs: pass_dv = pass_loose_dv_cuts(signal, tree)
            elif use_CR_for_DVs: pass_dv = pass_dv_CR_cuts(signal, tree)
            else: pass_dv = pass_dv_cuts(signal, tree)
            if pass_dv:
                num_dv = num_dv+1.0
        
        for j in range(plep_nentries):
            tree_shuffle.GetEntry(j)
            if use_loose_DVs: pass_plep = pass_loose_dv_cuts(signal, tree_shuffle)
            elif looser_plep: pass_plep = pass_looser_prompt_cuts(signal, tree_shuffle)
            else: pass_plep = pass_dv_cuts(signal, tree_shuffle)
            if pass_plep:
                num_plep = num_plep+1.0
        
        if use_CR_for_DVs:
            #get number of SR DVs from SS bkg file
            for k in range(plep_nentries):
                tree_shuffle.GetEntry(k)
                # if pass_dv_cuts(signal, tree_shuffle):
                #     num_loose_dv = num_loose_dv+1.0
                if pass_dv_cuts(signal, tree_shuffle):
                    num_SR_dv = num_SR_dv+1.0
        
        if mix_channels and (signal == "euu" or signal == "uee"): num_SR_dv = 1.0
        # if use_CR_for_DVs: shuffling_factor =num_dv*num_plep*num_SR_dv/num_loose_dv**2
        if use_CR_for_DVs: shuffling_factor =num_dv*num_plep/num_SR_dv
        elif mix_channels: shuffling_factor = 1 # calculate shuffling factor later
        else: shuffling_factor =num_dv*num_plep/num_dv

        print ("number of dvs: ", num_dv)
        print ("number of prompt leptons: ", num_plep)
        if use_CR_for_DVs: print ("number of SR events: ", num_SR_dv)
        print ("shuffling factor: ", shuffling_factor)

        output_text = output_text + f"""
                                    number of dvs: {num_dv}
                                    number of prompt leptons: {num_plep}
                                    """
        if use_CR_for_DVs: output_text = output_text + f"""
                                        number of SR events: {num_SR_dv}
                                        number of loose SR events (no mHNL cut): {num_loose_dv}
                                        shuffling factor: {num_plep*num_dv}
                                        1/w: {shuffling_factor}
                                        """
        else: output_text = output_text + f"""
                                        number of SR events: {num_dv}
                                        shuffling factor: {shuffling_factor}
                                        """

        with open(outputDir+f'info_about_shuffling_{signal}_{CR_charge}.txt', 'w') as f:
                print(output_text, file=f)
        

        for i in range(dv_nentries):
            if i % 10000 == 0:
                print( "Processing event {} / {}".format(i, dv_nentries))
            if signal == "uuu":
                if i % 10 == 0:
                    print("Processing event {} / {}".format(i, dv_nentries))

            tree.GetEntry(i)

            rDV = np.sqrt(tree.DV_x**2 + tree.DV_y**2)
            dv = ROOT.TVector3( tree.DV_x, tree.DV_y ,tree.DV_z)
            pv = ROOT.TVector3(tree.PV_x,tree.PV_y,tree.PV_z)
            
            pt_0 = tree.DV_trk_0_pt
            eta_0 = tree.DV_trk_0_eta
            phi_0 = tree.DV_trk_0_phi
            pt_1 = tree.DV_trk_1_pt
            eta_1 = tree.DV_trk_1_eta
            phi_1 = tree.DV_trk_1_phi
            M = 0.139 # pion mass assumption
            
            lepVec_0 = ROOT.TLorentzVector()
            lepVec_1 = ROOT.TLorentzVector()
            lepVec_0.SetPtEtaPhiM(pt_0, eta_0, phi_0, M)
            lepVec_1.SetPtEtaPhiM(pt_1, eta_1, phi_1, M)
            DV_4vec = lepVec_0 + lepVec_1
            elVec = [lepVec_0,lepVec_1]
            muVec = []

            if not use_CR_for_DVs:
                plep_vec = ROOT.TLorentzVector()
                plep_pt = tree.plep_pt
                plep_eta = tree.plep_eta
                plep_phi= tree.plep_phi
                plep_vec.SetPtEtaPhiM(plep_pt, plep_eta, plep_phi, M)

                Mhnl_original = selections.Mhnl(tree, "ee", plep=plep_vec, dMu=muVec,dEl=elVec,truth_pv=pv, truth_dv=dv,use_truth=True)
                Mlll_original = selections.Mlll(dv_type="ee", plep=plep_vec, dMu=muVec, dEl=elVec)
                mll_0 = lepVec_0 + plep_vec
                mll_1 = lepVec_1 + plep_vec


            # only fill histograms if DV passes some basic selections
            if use_loose_DVs: pass_dv = pass_loose_dv_cuts(signal, tree)
            elif use_CR_for_DVs: pass_dv = pass_dv_CR_cuts(signal, tree)
            else: pass_dv = pass_dv_cuts(signal, tree)
            if pass_dv: 
                # DV variables
                fill_hist("SS_bkg", "mDV", DV_4vec.M())
                fill_hist("SS_bkg", "DVpt", DV_4vec.Pt())
                fill_hist("SS_bkg", "DVeta", DV_4vec.Eta())
                fill_hist("SS_bkg", "DVphi", DV_4vec.Phi())
                fill_hist("SS_bkg", "Lxy", rDV)
                
                # plep variables
                if not use_CR_for_DVs:
                    fill_hist("SS_bkg", "plep_pt", plep_pt)
                    fill_hist("SS_bkg", "plep_eta", plep_eta)
                    fill_hist("SS_bkg", "plep_phi", plep_phi)
                
                    # HNL variables
                    fill_hist("SS_bkg", "HNLm_orig", Mhnl_original.mhnl)
                    fill_hist("SS_bkg", "mvis_orig", Mlll_original.mlll)
                    fill_hist("SS_bkg", "mll_1", mll_1.M())
                    fill_hist("SS_bkg", "mll_0", mll_0.M())                
                
                micro_ntuples["SS_bkg"].fill()
                
                for j in range(plep_nentries): 
                    tree_shuffle.GetEntry(j)
                    if use_loose_DVs: plep_pass = pass_loose_dv_cuts(signal, tree_shuffle)
                    elif looser_plep: plep_pass = pass_looser_prompt_cuts(signal,tree_shuffle)
                    else: plep_pass = pass_dv_cuts(signal,tree_shuffle)
                    if plep_pass:
                        plep_shuffle_vec = ROOT.TLorentzVector()
                        plep_shuffle_pt = tree_shuffle.plep_pt
                        plep_shuffle_eta = tree_shuffle.plep_eta
                        plep_shuffle_phi= tree_shuffle.plep_phi
                        plep_shuffle_vec.SetPtEtaPhiM(plep_shuffle_pt, plep_shuffle_eta, plep_shuffle_phi, M)
                        pv_shuffle = ROOT.TVector3(tree_shuffle.PV_x,tree_shuffle.PV_y,tree_shuffle.PV_z)
                        
                        if shuffle_PV: pv = pv_shuffle
                        else: pv=pv

                        Mhnl_shuffle = selections.Mhnl(tree, "ee", plep=plep_shuffle_vec, dMu=muVec,dEl=elVec,truth_pv=pv, truth_dv=dv,use_truth=True)
                        Mlll_shuffle = selections.Mlll("ee", plep=plep_shuffle_vec, dMu=muVec,dEl=elVec)
            
                        fill_hist(selection,"HNLm", Mhnl_shuffle.mhnl, weight = 1.0/shuffling_factor)
                        fill_hist(selection,"mvis", Mlll_shuffle.mlll, weight = 1.0/shuffling_factor)
                        fill_hist(selection, "DV_mass", tree.DV_mass)
                        fill_hist(selection, "DV_r", tree.DV_r)
                        fill_hist(selection, "DV_z", tree.DV_z)
                        fill_hist(selection, "DV_index", i)
                        fill_hist(selection, "plep_index", j)
                        fill_hist(selection, "DV_charge", tree.DV_charge)
                        fill_hist(selection, "DV_2medium", tree.DV_2medium)
                        fill_hist(selection, "DV_2loose", tree.DV_2loose)
                        fill_hist(selection, "DV_medium_veryveryloose",tree.DV_medium_veryveryloose)
                        fill_hist(selection, "DV_2veryveryloose",tree.DV_2veryveryloose)
                        fill_hist(selection, "DV_alpha",tree.DV_alpha)
                        fill_hist(selection, "DV_cosmic_sep",tree.DV_cosmic_sep)
                        fill_hist(selection, "DV_pass_lep_pt",tree.DV_pass_lep_pt)
                        # no prompt lepton in CR so there is no mlll variable
                        # always make shuffled pass Z veto by setting mll_is_SS = 1
                        if use_CR_for_DVs: 
                            fill_hist(selection, "mll_dMu_plep_is_OS",0)
                            fill_hist(selection, "mll_dMu_plep_is_SS",1)
                            fill_hist(selection, "mll_dMu_plep",0.0)
                            fill_hist(selection, "mll_dEl_plep_is_OS",0)
                            fill_hist(selection, "mll_dEl_plep_is_SS",1)
                            fill_hist(selection, "mll_dEl_plep",0.0)
                        else:
                            if signal == "uue" or signal == "uuu":
                                fill_hist(selection, "mll_dMu_plep_is_OS",tree.mll_dMu_plep_is_OS)
                                fill_hist(selection, "mll_dMu_plep_is_SS",tree.mll_dMu_plep_is_SS)
                                fill_hist(selection, "mll_dMu_plep",tree.mll_dMu_plep)
                            if signal == "eee" or signal == "eeu":
                                fill_hist(selection, "mll_dEl_plep_is_OS",tree.mll_dEl_plep_is_OS)
                                fill_hist(selection, "mll_dEl_plep_is_SS",tree.mll_dEl_plep_is_SS)
                                fill_hist(selection, "mll_dEl_plep",tree.mll_dEl_plep)
                        fill_hist(selection, "DV_pass_mat_veto",tree.DV_pass_mat_veto)
                        fill_hist(selection, "DV_weight_LNC_only", 1.0/shuffling_factor)
                        fill_hist(selection, "DV_weight_LNC_plus_LNV", 1.0/shuffling_factor)
                        micro_ntuples[selection].fill()
        print("Done looping through events. Writing ntuples to output file!")
        # write the output to ntuple
        file.cd()    
        [ntuple.write(vtx_alg+'_ntuples_'+key) for key, ntuple in micro_ntuples.items()]
        for key in h:
            h[key].Write(key)
        print("Done writing ntuples")
        # file = micro_ntuples[selection].ttree.GetCurrentFile()
        file.Close()
        del micro_ntuples # Clean up memory
        print("Close file")
        print("Job finished successfully")

    if plot: 
        # Define signal plotting colours 
        if signal == "uuu": 
            samples = ["uu_SS", "uuu_10G"]
            ratio_samples = [ "uuu_10G","uu_SS"]
        if signal == "eee": 
            samples = ["ee_SS", "eee_10G"]
            ratio_samples = [ "eee_10G","ee_SS"]
        if signal == "uue": 
            samples = ["eu_SS", "uue_10G"]
            ratio_samples = [ "uue_10G","eu_SS"]
        if signal == "eeu": 
            samples = ["eu_SS", "eeu_10G"]
            ratio_samples = [ "eeu_10G","eu_SS"]
        if signal == "uee": 
            samples = ["ee_SS", "uee_10G"]
            ratio_samples = [ "uee_10G","ee_SS"]
        if signal == "euu": 
            samples = ["uu_SS", "euu_10G"]
            ratio_samples = [ "euu_10G","uu_SS"]
        
        #defaults for plotting
        norm = False
        log_scale_y = True
        draw_markers= True

        #set extra info for plotting
        if signal == "uuu": 
            ch_str = "\\mu-\\mu\\mu"
            y_min = 0.0001
        if signal == "eee": 
            ch_str = "e-ee"
            y_min = 0.0001
        if signal == "uue": 
            ch_str = "\\mu-\\mue"
            y_min = 0.8
        if signal == "eeu": 
            ch_str = "e-e\\mu"
            y_min = 0.8
        if signal == "uee": 
            ch_str = "\\mu-ee"
            y_min = 0.8
        if signal == "euu": 
            ch_str = "e-\\mu\\mu"
            y_min = 0.8
        extra_info = "channel: {}".format(ch_str)

        pfile  = ROOT.TFile.Open(outputDir+"plots/" +  ouputf_name)
        h_to_plot ={}
        for key in h:
            h_to_plot[key] = pfile.Get(key) 
        
        
        # shuffling_factor = 32.0
        #KS test for HNL histograms
        # https://root.cern/doc/master/classTH1.html#TH1:KolmogorovTest
        h_to_plot["HNLm"].KolmogorovTest( h_to_plot["HNLm_orig"],option="D" ) # shuffled

        h_to_plot["mvis"].KolmogorovTest( h_to_plot["mvis_orig"],option="D" ) # shuffled

        shuffling_factor = 59


        Hist1D(hist_channels= [],
        hists = [h_to_plot["plep_pt"]],
        labels = ["ch. {}".format(ch_str) ],
        vtx_alg = vtx_alg,
        types = samples,
        outputDir = outputDir,
        name="plep_pt",
        x_title ="Prompt lepton p_{T}",
        x_units ="GeV",
        y_min = y_min,
        x_min = 0,
        x_max = 200,
        scaleLumi = 139.0,
        norm = norm,
        log_scale_y = log_scale_y,
        draw_markers = draw_markers,
        atlas_mod = "Internal",
        extra_legend_lines = [vtx_alg]
        )

        Hist1D(hist_channels= [],
        hists = [h_to_plot["plep_eta"]],
        labels = ["ch. {}".format(ch_str)],
        vtx_alg = vtx_alg,
        types = samples,
        outputDir = outputDir,
        name="plep_eta",
        x_title ="Prompt Lepton \\eta",
        x_units ="",
        y_min = y_min,
        x_min = -3,
        x_max = 3,
        ntup_nbins =100,
        scaleLumi = 139.0,
        norm = norm,
        log_scale_y = log_scale_y,
        draw_markers = draw_markers,
        atlas_mod = "Internal",
        extra_legend_lines = [vtx_alg ]
        )

        Hist1D(hist_channels= [],
        hists = [h_to_plot["plep_phi"]],
        labels = ["ch. {}".format(ch_str)],
        vtx_alg = vtx_alg,
        types = samples,
        outputDir = outputDir,
        name="plep_phi",
        x_title ="Prompt Lepton \\phi",
        x_units ="",
        y_min = y_min,
        x_min = -4,
        x_max = 4,
        ntup_nbins =100,
        scaleLumi = 139.0,
        norm = norm,
        log_scale_y = log_scale_y,
        draw_markers = draw_markers,
        atlas_mod = "Internal",
        extra_legend_lines = [vtx_alg ]
        )

        Hist1DRatio(hist_channels= [],
        hists = [h_to_plot["HNLm"],h_to_plot["HNLm_orig"]],
        labels = ["random plep","orig. plep"],
        vtx_alg = vtx_alg,
        types = ratio_samples,
        outputDir = outputDir,
        name="HNLm",
        x_title ="HNL mass",
        x_units ="GeV",
        # y_min = 0.1,
        y_min = 0.1,
        x_min = 0,
        x_max = 50,
        # rebin = 2,
        ratio_ymin = -0.5,
        ratio_ymax =2.5,
        ntup_nbins= 50,
        scaleLumi = 139.0,
        show_overflow = False,
        show_underflow = False,
        norm = norm,
        empty_scale = 3.8,
        log_scale_y = log_scale_y,
        draw_markers = draw_markers,
        atlas_mod = "Internal",
        extra_legend_lines = [vtx_alg,extra_info, f"shuffling factor f={shuffling_factor}"]
        )

        Hist1DRatio(hist_channels= [],
        hists = [h_to_plot["mvis"],h_to_plot["mvis_orig"]],
        labels = ["random plep","orig. plep"],
        vtx_alg = vtx_alg,
        types = ratio_samples,
        outputDir = outputDir,
        name="mvis",
        x_title ="Tri-lepton mass",
        x_units ="GeV",
        y_min = 0.1,
        x_min = 0,
        x_max = 200,
        ratio_ymin = -0.5,
        ratio_ymax = 2.5,
        # rebin =2,
        ntup_nbins =100,
        scaleLumi = 139.0,
        show_overflow = False,
        show_underflow = False,
        norm = norm,
        empty_scale = 3.8,
        log_scale_y = log_scale_y,
        draw_markers = draw_markers,
        atlas_mod = "Internal",
        extra_legend_lines = [vtx_alg,extra_info,f"shuffling factor f={shuffling_factor}"]
        )

        Hist1D(hist_channels= [],
        hists = [h_to_plot["Lxy"]],
        labels = ["ch. {}".format(ch_str)],
        vtx_alg = vtx_alg,
        types = samples,
        outputDir = outputDir,
        name="DV_r",
        x_title ="L_{xy}",
        x_units ="mm",
        y_min = y_min,
        x_min = 0,
        x_max = 300,
        ntup_1D= True,
        ntup_nbins =60,
        scaleLumi = 139.0,
        empty_scale = 1.8,
        norm = norm,
        log_scale_y = log_scale_y,
        draw_markers = draw_markers,
        atlas_mod = "Internal",
        extra_legend_lines = [vtx_alg ]
        )
            
        Hist1D(hist_channels= [],
        hists = [h_to_plot["mDV"]],
        labels = ["ch. {}".format(ch_str)],
        vtx_alg = vtx_alg,
        types = samples,
        outputDir = outputDir,
        name="DV_mass",
        x_title ="DV mass",
        x_units ="GeV",
        # y_min = y_min,
        x_min = 0,
        x_max = 50,
        ntup_nbins =0,
        scaleLumi = 139.0,
        ntup_1D = True,
        norm = norm,
        draw_cut = False,
        cut = 4,
        empty_scale = 1.8,
        hide_lumi = True,
        log_scale_y = log_scale_y,
        draw_markers = draw_markers,
        atlas_mod = "Internal",
        extra_legend_lines = [vtx_alg ]
        )

        Hist1D(hist_channels= [],
        hists = [h_to_plot["mll_0"]],
        labels = ["SS bkg. {}".format(ch_str)],
        vtx_alg = vtx_alg,
        types = samples,
        outputDir = outputDir,
        name="mll_0",
        x_title ="Mass of plep + Leading p_{T} Track",
        x_units ="GeV",
        # y_min = y_min,
        x_min = 0,
        x_max = 200,
        ntup_nbins =0,
        scaleLumi = 139.0,
        ntup_1D = True,
        norm = norm,
        draw_cut = False,
        cut = 4,
        rebin = 2,
        empty_scale = 1.8,
        hide_lumi = True,
        log_scale_y = log_scale_y,
        draw_markers = draw_markers,
        atlas_mod = "Internal",
        extra_legend_lines = [vtx_alg ]
        )

        Hist1D(hist_channels= [],
        hists = [h_to_plot["mll_1"]],
        labels = ["SS bkg. {}".format(ch_str)],
        vtx_alg = vtx_alg,
        types = samples,
        outputDir = outputDir,
        name="mll_1",
        x_title ="Mass of plep + Sub-leading p_{T} Track",
        x_units ="GeV",
        # y_min = y_min,
        x_min = 0,
        x_max = 200,
        ntup_nbins =0,
        scaleLumi = 139.0,
        ntup_1D = True,
        norm = norm,
        draw_cut = False,
        cut = 4,
        rebin = 2,
        empty_scale = 1.8,
        hide_lumi = True,
        log_scale_y = log_scale_y,
        draw_markers = draw_markers,
        atlas_mod = "Internal",
        extra_legend_lines = [vtx_alg ]
        )

        Hist1D(hist_channels= [],
        hists = [h_to_plot["DVpt"]],
        labels = ["ch. {}".format(ch_str)],
        vtx_alg = vtx_alg,
        types = samples,
        outputDir = outputDir,
        name="DV_pt",
        x_title ="DV p_{T}",
        x_units ="GeV",
        y_min = y_min,
        x_min = 0,
        x_max = 100,
        ntup_nbins =100,
        scaleLumi = 139.0,
        norm = norm,
        log_scale_y = log_scale_y,
        draw_markers = draw_markers,
        atlas_mod = "Internal",
        extra_legend_lines = [vtx_alg ]
        )


        Hist1D(hist_channels= [],
        hists = [h_to_plot["DVeta"]],
        labels = ["ch. {}".format(ch_str)],
        vtx_alg = vtx_alg,
        types = samples,
        outputDir = outputDir,
        name="DV_eta",
        x_title ="DV \\eta",
        x_units ="",
        y_min = y_min,
        x_min = -3,
        x_max = 3,
        ntup_nbins =100,
        scaleLumi = 139.0,
        norm = norm,
        log_scale_y = log_scale_y,
        draw_markers = draw_markers,
        atlas_mod = "Internal",
        extra_legend_lines = [vtx_alg ]
        )


        Hist1D(hist_channels= [],
        hists = [h_to_plot["DVphi"]],
        labels = ["ch. {}".format(ch_str)],
        vtx_alg = vtx_alg,
        types = samples,
        outputDir = outputDir,
        name="DV_phi",
        x_title ="DV \\phi",
        x_units ="",
        y_min = y_min,
        x_min = -4,
        x_max = 4,
        ntup_nbins =100,
        scaleLumi = 139.0,
        norm = norm,
        log_scale_y = log_scale_y,
        draw_markers = draw_markers,
        atlas_mod = "Internal",
        extra_legend_lines = [vtx_alg ]
        )   

if __name__ == "__main__":
    import argparse

    class AppendActionCleanDefault(argparse._AppendAction):
	    def __init__(self, *args, **kwargs):
		    super(argparse._AppendAction, self).__init__(*args,**kwargs)
		    self.index = 0
	    def __call__(self, parser, namespace, values, option_string = None):
		    items = argparse._copy.copy(argparse._ensure_value(namespace, self.dest, [])) if self.index else []
		    if values:
			    self.index += 1
			    items.append(values)
			    setattr(namespace, self.dest, items)

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

    parser.add_argument("--channel",
                        required = True,
						dest="signal",
						type = str,
						help="Choose what channel you want to run or plot shuffled background estimate."
				        )
    
    parser.add_argument("--version",
						dest="version",
						type = str,
                        default = "v6",
						help="Choose the version of the ntuples that you want to run on."
				        )

    parser.add_argument("-r", "--run",
						dest="makeHist",
						action="store_true",
						default = False,
						help="Run shuffled background estimate. This will process run through all the events and output the shuffled background mini-trees"
						)
    
    parser.add_argument("-p", "--plot",
						dest="plot",
						action="store_true",
						default = False,
						help="Plot the output of the shuffled background estimate."
				        )
    
    parser.add_argument("--VR_charge",
						dest="CR_charge",
						type = str, 
                        default = "OS",
						help="The charge selection you want to use if running background estimate with DVs selected from the VR."
				        )
    
    parser.add_argument("--select_plep_from_loose_regionB",
						dest="looser_plep",
						action="store_true",
						default = False,
						help="Add this flag to use looser prompt lepton (region B') selection."
				        )
    
    parser.add_argument("--use_loose_DVs",
						dest="use_loose_DVs",
						action="store_true",
						default = False,
						help="Add this flag to use loose DV selection. (This option is included for testing purposes only)."
				        )
    
    parser.add_argument("--select_DV_from_VR",
						dest="use_CR_for_DVs",
						action="store_true",
						default = False,
						help="Add this flag to select DVs from the inverted prompt lepton validation region."
				        )
    
    parser.add_argument("--select_plep_from_diff_channel",
						dest="mix_channels",
						action="store_true",
						default = False,
						help="Add this flag to use loose DV selection. This option is included to run prompt lepton systematic."
				        )

    parser.add_argument("--shuffle_PV",
						dest="shuffle_PV",
						action="store_true",
						default = False,
						help="Add this flag to shuffle PV location. (This option is included for testing purposes only)."
				        )


    parent_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter, parents = [parser])
    options = parent_parser.parse_args()
	
    main()


