import helpers

# trigger systematics order
# TODO: systematics are repeated in data tree (e.g. len(muon_TrigEff_SF_HLT_mu26_ivarmedium_RecoMedium)[0] = 15 when it should = 5 (nominal +stat up/down +syst up/down)
muon_trigger_systematics = ['nominal', 'MUON_EFF_TrigStatUncertainty__1down', 'MUON_EFF_TrigStatUncertainty__1up', 'MUON_EFF_TrigSystUncertainty__1down', 'MUON_EFF_TrigSystUncertainty__1up']


def get_trigger_scale_factor(analysis):
    """
    Calculates the scale factor associated with the lepton trigger correction.
    @param analysis: Analysis class instance
    @return: Scale factor for lepton trigger
    """
    # getting nominal value from systematic array (update this later to deal with systematics)
    syst = 0

    # set default to zero
    trigger_scale_factor = 0

    # get index and type from instance of RequireMediumTriggerMatching class
    lep_index = analysis.trigger_matched_medium.trigger_matched_lepton_index
    lep_type = analysis.trigger_matched_medium.trigger_matched_lepton_type
    if lep_type == 'muon':
        trigger_scale_factor_mu26_ivarmedium = analysis.tree['muon_TrigEff_SF_HLT_mu26_ivarmedium_RecoMedium'][lep_index][syst]
        trigger_scale_factor_mu20_iloose_L1MU15 = analysis.tree['muon_TrigEff_SF_HLT_mu20_iloose_L1MU15_RecoMedium'][lep_index][syst]
        if trigger_scale_factor_mu26_ivarmedium >= 0.0 and trigger_scale_factor_mu20_iloose_L1MU15 >= 0.0:
            analysis.logger.warning("Two trigger scale factors. What do we do here?")
        if trigger_scale_factor_mu26_ivarmedium == 0.0 and trigger_scale_factor_mu20_iloose_L1MU15 == 0.0:
            analysis.logger.warning("Trigger scale factor is zero")
        trigger_scale_factor = max(trigger_scale_factor_mu26_ivarmedium, trigger_scale_factor_mu20_iloose_L1MU15)  # safe for now

    if lep_type == 'electron':
        # only ever use medium
        trigger_scale_factor_el_medium = analysis.tree[
            'el_TrigEff_SF_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_2018_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0_Medium_isolGradient'
        ][lep_index][syst]
        trigger_scale_factor = trigger_scale_factor_el_medium

    return trigger_scale_factor


def get_reco_scale_factor(analysis):
    """
    Calculates the combination of scale factors for all three leptons.
    For muons the working point is Medium.
    For electrons the working point is Medium or LooseBLayer,
    depending on if the electron is prompt or displaced respectively.
    @param analysis: Analysis class instance
    @param muons: Collection of displaced muons. Should be len 2 or None
    @param electrons: Collection if displaced electrons. Should be len 2 or None
    @return: Scale factor for lepton reconstruction/identification
    """
    # getting nominal value from systematic array (update this later to deal with systematics)
    syst = 0

    # get muons and electrons
    muons = helpers.Muons(analysis.tree)
    electrons = helpers.Electrons(analysis.tree)

    # set default to 0
    prompt_lepton_scale_factor = displaced_lepton_0_scale_factor = displaced_lepton_1_scale_factor = 0

    # get prompt scale factors
    if analysis.plep == 'muon':
        prompt_lepton_scale_factor = analysis.tree['muon_RecoEff_SF_RecoMedium'][analysis.plep_sel.plep_index][syst]
    elif analysis.plep == 'electron':
        prompt_lepton_scale_factor = analysis.tree['el_PIDEff_SF_Medium'][analysis.plep_sel.plep_index][syst]

    # get displaced scale factors depending on type
    if analysis.dv_type == 'mumu':
        displaced_lepton_0_scale_factor = analysis.tree.get('muon_RecoEff_SF_RecoMedium')[muons.lepIndex[0]][syst]
        displaced_lepton_1_scale_factor = analysis.tree.get('muon_RecoEff_SF_RecoMedium')[muons.lepIndex[1]][syst]
    if analysis.dv_type == 'emu':
        displaced_lepton_0_scale_factor = analysis.tree.get('muon_RecoEff_SF_RecoMedium')[muons.lepIndex[0]][syst]
        displaced_lepton_1_scale_factor = analysis.tree.get('el_PIDEff_SF_LooseBLayer')[electrons.lepIndex[0]][syst]
        displaced_lepton_1_scale_factor *= analysis.tree.get('el_RecoEff_SF')[electrons.lepIndex[0]][syst]  # electron SF is split into ID and Reco
    if analysis.dv_type == 'ee':
        displaced_lepton_0_scale_factor = analysis.tree.get('el_PIDEff_SF_LooseBLayer')[electrons.lepIndex[0]][syst]
        displaced_lepton_0_scale_factor *= analysis.tree.get('el_RecoEff_SF')[electrons.lepIndex[0]][syst]  # electron SF is split into ID and Reco
        displaced_lepton_1_scale_factor = analysis.tree.get('el_PIDEff_SF_LooseBLayer')[electrons.lepIndex[1]][syst]
        displaced_lepton_1_scale_factor *= analysis.tree.get('el_RecoEff_SF')[electrons.lepIndex[1]][syst]  # electron SF is split into ID and Reco

    return prompt_lepton_scale_factor * displaced_lepton_0_scale_factor * displaced_lepton_1_scale_factor
