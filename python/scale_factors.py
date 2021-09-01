import helpers

# trigger systematics order
# TODO: systematics are repeated in data tree (e.g. len(muon_TrigEff_SF_HLT_mu26_ivarmedium_RecoMedium)[0] = 15 when it should = 5 (nominal +stat up/down +syst up/down)
muon_trigger_systematics = ['nominal', 'MUON_EFF_TrigStatUncertainty__1down', 'MUON_EFF_TrigStatUncertainty__1up', 'MUON_EFF_TrigSystUncertainty__1down', 'MUON_EFF_TrigSystUncertainty__1up']

# indices for scale factor systematics as stored in trees
syst_index = {
    'nominal': 0,
    # MuonEfficiencyCorrector_RecoSyst_RecoMedium
    'MUON_EFF_RECO_SYS__1down': 7,
    'MUON_EFF_RECO_SYS__1up': 8,
    # MuonEfficiencyCorrector_TrigSyst_RecoMedium
    'MUON_EFF_TrigSystUncertainty__1down': 3,
    'MUON_EFF_TrigSystUncertainty__1up': 4,
    # EleEffCorr_RecoSyst_Reconstruction
    'EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR__1down': 1,
    'EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR__1up': 2,
    # EleEffCorr_PIDSyst_LooseBLayer & EleEffCorr_PIDSyst_Medium
    'EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR__1down': 1,
    'EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR__1up': 2,
    # EleEffCorr_TrigSyst_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_Medium_isolGradient
    'EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1down': 1,
    'EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR__1up': 2,
}


def get_trigger_scale_factor(analysis, systematic='nominal'):
    """
    Calculates the scale factor associated with the lepton trigger correction.
    @param analysis: Analysis class instance
    @param systematic: the systematic variation to use ('nominal', 'down', or 'up')
    @return: Scale factor for lepton trigger
    """
    # set default to zero
    trigger_scale_factor = 0

    # get index and type from instance of RequireMediumTriggerMatching class
    lep_index = analysis.trigger_matched_medium.trigger_matched_lepton_index
    lep_type = analysis.trigger_matched_medium.trigger_matched_lepton_type
    if lep_type == 'muon':
        trigger_scale_factor_mu26_ivarmedium = analysis.tree['muon_TrigEff_SF_HLT_mu26_ivarmedium_RecoMedium'][lep_index]
        if len(trigger_scale_factor_mu26_ivarmedium) > 1:  # get systematic
            trigger_scale_factor_mu26_ivarmedium = trigger_scale_factor_mu26_ivarmedium[
                syst_index[systematic if 'MUON_EFF_Trig' in systematic else 'nominal']
            ]
        else:
            trigger_scale_factor_mu26_ivarmedium = trigger_scale_factor_mu26_ivarmedium[syst_index['nominal']]
        trigger_scale_factor_mu20_iloose_L1MU15 = analysis.tree['muon_TrigEff_SF_HLT_mu20_iloose_L1MU15_RecoMedium'][lep_index]
        if len(trigger_scale_factor_mu20_iloose_L1MU15) > 1:  # get systematic
            trigger_scale_factor_mu20_iloose_L1MU15 = trigger_scale_factor_mu20_iloose_L1MU15[
                syst_index[systematic if 'MUON_EFF_Trig' in systematic else 'nominal']
            ]
        else:
            trigger_scale_factor_mu20_iloose_L1MU15 = trigger_scale_factor_mu20_iloose_L1MU15[syst_index['nominal']]
        # deal wih any weird circumstances
        if trigger_scale_factor_mu26_ivarmedium >= 0.0 and trigger_scale_factor_mu20_iloose_L1MU15 >= 0.0:
            analysis.logger.warning("Two trigger scale factors. What do we do here?")
        if trigger_scale_factor_mu26_ivarmedium == 0.0 and trigger_scale_factor_mu20_iloose_L1MU15 == 0.0:
            analysis.logger.warning("Trigger scale factor is zero")
        trigger_scale_factor = max(trigger_scale_factor_mu26_ivarmedium, trigger_scale_factor_mu20_iloose_L1MU15)  # safe for now

    if lep_type == 'electron':
        # only ever use medium
        trigger_scale_factor_el_medium = analysis.tree[
            'el_TrigEff_SF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_Medium_isolGradient'
        ][lep_index]
        if len(trigger_scale_factor_el_medium) >= 1:
            trigger_scale_factor_el_medium = trigger_scale_factor_el_medium[
                syst_index[systematic if 'EL_EFF_Trigger' in systematic else 'nominal']
            ]
        trigger_scale_factor = trigger_scale_factor_el_medium

    return trigger_scale_factor


def get_reco_scale_factor(analysis, systematic='nominal'):
    """
    Calculates the combination of scale factors for all three leptons.
    For muons the working point is Medium.
    For electrons the working point is Medium or LooseBLayer,
    depending on if the electron is prompt or displaced respectively.
    @param analysis: Analysis class instance
    @param systematic: the systematic variation to use ('nominal', 'down', or 'up')
    @return: Scale factor for lepton reconstruction/identification
    """
    # get muons and electrons
    muons = helpers.Muons(analysis.tree)
    electrons = helpers.Electrons(analysis.tree)

    # set default to 0
    prompt_lepton_scale_factor = displaced_lepton_0_scale_factor = displaced_lepton_1_scale_factor = 0

    # get prompt scale factors
    if analysis.plep == 'muon':
        prompt_lepton_scale_factor = analysis.tree['muon_RecoEff_SF_RecoMedium'][analysis.plep_sel.plep_index][
            syst_index[systematic if 'MUON_EFF_RECO_SYS' in systematic else 'nominal']  # only get systematic related to the quantity of interest
        ]
    elif analysis.plep == 'electron':
        # electron SF is split into ID and Reco
        prompt_lepton_scale_factor = analysis.tree['el_PIDEff_SF_Medium'][analysis.plep_sel.plep_index][
            syst_index[systematic if 'EL_EFF_ID' in systematic else 'nominal']
        ]
        prompt_lepton_scale_factor *= analysis.tree['el_RecoEff_SF'][analysis.plep_sel.plep_index][
            syst_index[systematic if 'EL_EFF_Reco' in systematic else 'nominal']
        ]

    # get displaced scale factors depending on type
    if analysis.dv_type == 'mumu':
        displaced_lepton_0_scale_factor = analysis.tree.get('muon_RecoEff_SF_RecoMedium')[muons.lepIndex[0]][
            syst_index[systematic if 'MUON_EFF_RECO_SYS' in systematic else 'nominal']  # only get systematic related to the quantity of interest
        ]
        displaced_lepton_1_scale_factor = analysis.tree.get('muon_RecoEff_SF_RecoMedium')[muons.lepIndex[1]][
            syst_index[systematic if 'MUON_EFF_RECO_SYS' in systematic else 'nominal']
        ]
    if analysis.dv_type == 'emu':
        displaced_lepton_0_scale_factor = analysis.tree.get('muon_RecoEff_SF_RecoMedium')[muons.lepIndex[0]][
            syst_index[systematic if 'MUON_EFF_RECO_SYS' in systematic else 'nominal']
        ]
        displaced_lepton_1_scale_factor = analysis.tree.get('el_PIDEff_SF_LooseBLayer')[electrons.lepIndex[0]][
            syst_index[systematic if 'EL_EFF_ID' in systematic else 'nominal']
        ]
        displaced_lepton_1_scale_factor *= analysis.tree.get('el_RecoEff_SF')[electrons.lepIndex[0]][
            syst_index[systematic if 'EL_EFF_Reco' in systematic else 'nominal']
        ]  # electron SF is split into ID and Reco
    if analysis.dv_type == 'ee':
        displaced_lepton_0_scale_factor = analysis.tree.get('el_PIDEff_SF_LooseBLayer')[electrons.lepIndex[0]][
            syst_index[systematic if 'EL_EFF_ID' in systematic else 'nominal']
        ]
        displaced_lepton_0_scale_factor *= analysis.tree.get('el_RecoEff_SF')[electrons.lepIndex[0]][
            syst_index[systematic if 'EL_EFF_Reco' in systematic else 'nominal']
        ]  # electron SF is split into ID and Reco
        displaced_lepton_1_scale_factor = analysis.tree.get('el_PIDEff_SF_LooseBLayer')[electrons.lepIndex[1]][
            syst_index[systematic if 'EL_EFF_ID' in systematic else 'nominal']
        ]
        displaced_lepton_1_scale_factor *= analysis.tree.get('el_RecoEff_SF')[electrons.lepIndex[1]][
            syst_index[systematic if 'EL_EFF_Reco' in systematic else 'nominal']
        ]  # electron SF is split into ID and Reco

    return prompt_lepton_scale_factor * displaced_lepton_0_scale_factor * displaced_lepton_1_scale_factor
