import numpy as np


# ____________________________________________________________
# Displaced vertexing + tracking systematic

vertexing_syst_map = np.array(
	  [[0.96      , 0.95993126, 0.94165127, 0.95836189, 0.94467486,
        0.84039676, 0.8291002 , 0.78884503, 0.77728849, 0.77760779,
        0.79450637, 0.83282275, 0.79899139, 0.83757042, 0.78762539],
       [0.96      , 0.95954859, 0.91930915, 0.93324524, 0.95998072,
        0.92477593, 0.93487687, 0.8912892 , 0.86875911, 0.8374753 ,
        0.922888  , 0.88140525, 0.84216507, 0.82667875, 0.89791438],
       [0.96      , 0.91516189, 0.91564243, 0.9288241 , 0.91766646,
        0.92913114, 0.93132677, 0.91245808, 0.93288162, 0.90309618,
        0.9587421 , 0.92604629, 0.95570732, 0.95386959, 0.91366827],
       [0.96      , 0.93037907, 0.95971524, 0.95912841, 0.95985423,
        0.95361916, 0.94195473, 0.92836505, 0.86059293, 0.92934497,
        0.87246303, 0.95448154, 0.90785524, 0.85493355, 0.79353179],
       [0.96      , 0.9347185 , 0.88558086, 0.88861315, 0.91707327,
        0.91601057, 0.95553542, 0.95796652, 0.93774068, 0.92051822,
        0.95008698, 0.95669437, 0.94905746, 0.95980851, 0.89584635],
       [0.96      , 0.96      , 0.91899608, 0.89383716, 0.85820325,
        0.88271608, 0.92443339, 0.93637299, 0.95868375, 0.93592103,
        0.956569  , 0.95646187, 0.94138116, 0.95973943, 0.95900326],
       [0.96      , 0.96      , 0.88980506, 0.82210277, 0.75701872,
        0.80786463, 0.81370689, 0.94601174, 0.90242167, 0.89890168,
        0.93836038, 0.89594396, 0.93026478, 0.94260998, 0.82739249],
       [0.96      , 0.96      , 0.96      , 0.85732267, 0.80606151,
        0.87094161, 0.76462153, 0.84690591, 0.93168821, 0.9141967 ,
        0.8360955 , 0.81450598, 0.93008109, 0.68276225, 0.95574859],
       [0.96      , 0.96      , 0.96      , 0.96      , 0.71473217,
        0.87973791, 0.52818879, 0.59921592, 0.93377555, 0.92223419,
        0.84718595, 0.64849899, 0.95889304, 0.7524329 , 0.8681943 ]])

vertexing_syst_pt_bins = [2, 4, 6, 8, 10, 15, 20, 25, 35]
vertexing_syst_dvr_bins = [20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, np.inf]


def get_vertexing_uncertainty(r, pt):
    """retrieve the calculated vertexing + tracking uncertainty using numpy's digitize"""
    # no overflow was used in dvr bins so extending 300mm to infinity
    pt_bin = np.digitize(pt, vertexing_syst_pt_bins)
    dvr_bin = np.digitize(r, vertexing_syst_dvr_bins)
    # deal with overflow (underflow behaves as expected)
    if pt_bin >= len(vertexing_syst_pt_bins): pt_bin = len(vertexing_syst_pt_bins) - 1
    if dvr_bin >= len(vertexing_syst_dvr_bins): dvr_bin = len(vertexing_syst_dvr_bins) - 1
    return vertexing_syst_map[pt_bin, dvr_bin]


def plot_vertexing_uncertainty():
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
# Displaced lepton identification systematic

displaced_electron_uncertainty = [1.0, 0.9813603162765503, 0.9508885741233826, 0.9436681866645813, 0.9272117018699646, 0.9385538101196289, 0.9473791718482971]
displaced_muon_uncertainty = [1.0, 0.9852185249328613, 0.9768053293228149, 0.9713043570518494, 0.9561651945114136, 0.952542245388031, 0.9498283267021179]
displaced_lepton_uncertainty_d0_bins = [0, 3, 10, 20, 30, 40, 50, 60]

def get_displaced_lepton_uncertainty(d0, lepton):
    ind = int(np.digitize(abs(d0), displaced_lepton_uncertainty_d0_bins))
    ind = min(ind, len(displaced_lepton_uncertainty_d0_bins)-1) # deal with overflow
    if lepton == 'electron':
        return displaced_electron_uncertainty[ind-1]
    if lepton == 'muon':
        return displaced_muon_uncertainty[ind-1]


def get_combined_displaced_lepton_uncertainty(lepton_0_d0, lepton_0_type, lepton_1_d0, lepton_1_type):
    # get uncertainties
    lepton_0_uncertainty = get_displaced_lepton_uncertainty(lepton_0_d0, lepton_0_type)
    lepton_1_uncertainty = get_displaced_lepton_uncertainty(lepton_1_d0, lepton_1_type)
    # calculate in quadrature
    total_uncertainty = 1 - np.sqrt(np.square(1 - lepton_0_uncertainty) + np.square(1 - lepton_1_uncertainty))
    return total_uncertainty

def get_displaced_electron_uncertainty(d0):
    return get_displaced_lepton_uncertainty(d0, lepton='electron')

def get_displaced_muon_uncertainty(d0):
    return get_displaced_lepton_uncertainty(d0, lepton='muon')