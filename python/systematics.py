import numpy as np

# ____________________________________________________________
# Displaced vertexing + tracking systematic

# symmetrized map
vertexing_syst_map = np.array(
       [[1.        , 0.99765399, 0.95751971, 0.98843573, 0.96177865,
        0.84549048, 0.83384723, 0.79266834, 0.78091003, 0.78123461,
        0.79843703, 0.83767861, 0.80301152, 0.84257266, 0.79142634],
       [1.        , 0.9939737 , 0.92992138, 0.9465566 , 0.99875801,
        0.93629238, 0.94860912, 0.89891569, 0.87500332, 0.84247451,
        0.93407382, 0.88835452, 0.84731776, 0.83135761, 0.9060773 ],
       [1.        , 0.92518352, 0.92572888, 0.94112718, 0.92803603,
        0.94149876, 0.94417874, 0.92213096, 0.94610309, 0.91173704,
        0.9898899 , 0.9377975 , 0.98097787, 0.97702143, 0.923494  ],
       [1.        , 0.94301689, 0.99521855, 0.99160435, 0.99658204,
        0.97652272, 0.95793751, 0.94057302, 0.86645476, 0.94175798,
        0.87889807, 0.97827604, 0.91699002, 0.86055727, 0.79744354],
       [1.        , 0.94840858, 0.89280047, 0.89604313, 0.92735812,
        0.92614727, 0.98058096, 0.98708438, 0.95229023, 0.93131701,
        0.97014519, 0.9834055 , 0.96845412, 0.99608135, 0.90383356],
       [1.        , 1.        , 0.92956113, 0.90166105, 0.86396207,
        0.88974794, 0.93588828, 0.95051872, 0.98965435, 0.94993889,
        0.98308102, 0.98280788, 0.95714947, 0.99542689, 0.99101487],
       [1.        , 1.        , 0.89732126, 0.82665807, 0.76033377,
        0.81207448, 0.81805187, 0.96374076, 0.91099702, 0.90715136,
        0.95310178, 0.9039393 , 0.94287732, 0.95884645, 0.83209124],
       [1.        , 1.        , 1.        , 0.86304445, 0.81023136,
        0.87729683, 0.76804521, 0.85222382, 0.944624  , 0.9240908 ,
        0.84105131, 0.81887012, 0.94265321, 0.68529412, 0.98107417],
       [1.        , 1.        , 1.        , 1.        , 0.71755048,
        0.88658497, 0.52988744, 0.60121701, 0.94722048, 0.93331027,
        0.85251395, 0.65078236, 0.99052467, 0.75568572, 0.87441042]])

vertexing_syst_pt_bins = [2, 4, 6, 8, 10, 15, 20, 25, 35]
vertexing_syst_dvr_bins = [20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, np.inf]


def get_vertexing_systematic(r, pt):
    """retrieve the calculated vertexing + tracking uncertainty using numpy's digitize"""
    # no overflow was used in dvr bins so extending 300mm to infinity
    pt_bin = np.digitize(pt, vertexing_syst_pt_bins)
    dvr_bin = np.digitize(r, vertexing_syst_dvr_bins)
    # deal with overflow (underflow behaves as expected)
    if pt_bin >= len(vertexing_syst_pt_bins): pt_bin = len(vertexing_syst_pt_bins) - 1
    if dvr_bin >= len(vertexing_syst_dvr_bins): dvr_bin = len(vertexing_syst_dvr_bins) - 1
    return vertexing_syst_map[pt_bin, dvr_bin]


def plot_vertexing_systematic():
    """a utility for plotting the vertexing uncertainty if you want to see it"""
    import matplotlib.pyplot as plt
    vertexing_syst_dvr_bins[-1] = 300
    plt.figure(figsize=[10, 7])
    plt.pcolormesh(vertexing_syst_dvr_bins, vertexing_syst_pt_bins, vertexing_syst_map)
    plt.colorbar()
    plt.xlabel('DV Radius [mm]')
    plt.ylabel('DV pT [GeV]')
    plt.title('vertexing systematics map')


# ____________________________________________________________
# Displaced lepton d0 extrapolation systematic

d0_extrapolation_systematic_bins = [0, 3, 10, 20, 30, 40, 50, 60]
d0_extrapolation_electron_systematic = [1.0, 0.9797764420509338, 0.9509709477424622, 0.9440935254096985, 0.9337978959083557, 0.9336843490600586, 0.943777322769165]
d0_extrapolation_muon_systematic = [1.0, 0.983723521232605, 0.9758386015892029, 0.9704379439353943, 0.9621100425720215, 0.9447205662727356, 0.9487572908401489]


def get_d0_extrapolation_systematic(d0, lepton):
    ind = int(np.digitize(abs(d0), d0_extrapolation_systematic_bins))
    ind = min(ind, len(d0_extrapolation_systematic_bins) - 1)  # deal with overflow
    if lepton == 'electron':
        return d0_extrapolation_electron_systematic[ind - 1]
    if lepton == 'muon':
        return d0_extrapolation_muon_systematic[ind - 1]


def get_combined_d0_extrapolation_systematic(lepton_0_d0, lepton_0_type, lepton_1_d0, lepton_1_type):
    # get uncertainties
    lepton_0_systematic = get_d0_extrapolation_systematic(lepton_0_d0, lepton_0_type)
    lepton_1_systematic = get_d0_extrapolation_systematic(lepton_1_d0, lepton_1_type)
    # calculate in quadrature
    # total_systematic = 1 - np.sqrt(np.square(1 - lepton_0_systematic) + np.square(1 - lepton_1_systematic))
    # just kidding, calculate linearly since they are correlated systematics
    total_systematic = lepton_0_systematic * lepton_1_systematic
    return total_systematic


def get_d0_extrapolation_electron_systematic(d0):
    return get_d0_extrapolation_systematic(d0, lepton='electron')


def get_d0_extrapolation_muon_systematic(d0):
    return get_d0_extrapolation_systematic(d0, lepton='muon')
