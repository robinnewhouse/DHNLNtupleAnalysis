import ROOT
import helpers


class Observables:
	"""
	A class to contain all of the histograms in a directory
	"""

	def __init__(self):
		self.histogram_dict = {}  # stores a collection of Observable objects
		self.logger = helpers.getLogger('dHNLAnalysis.observables')

	def fill_hist(self, directory, hist_name, variable_1, variable_2, weight):
		"""
		This function will fill a histogram that is stored in the histogram_dict, and if it does not exist, will create it.
		:param directory: The path where the histogram will be stored in the root file.
		:param hist_name: The name of the histogram
		:param variable_1: Must be a number
		:param variable_2: Optional, will fill a 2d histogram
		:param weight: The weight assigned to this variable when filling
		"""
		if 'MCweight' in hist_name:
			weight = 1  # don't weight the weight

		# create histogram if it doesn't exist
		if not directory + hist_name in self.histogram_dict:
			self.add_histogram(directory, hist_name, is_2d=(variable_2 is not None))

		# retrieve the histogram from the collection
		histogram = self.histogram_dict[directory + hist_name]
		# fill 1d or 2d histogram
		if variable_2 is None:
			histogram.Fill(variable_1, weight)
		else:
			histogram.Fill(variable_1, variable_2, weight)

	def add_histogram(self, directory, hist_name, is_2d):
		"""
		Create a new 1D histogram. If a defined binning exists, use that. Otherwise use defaults.
		:param directory:
		:param hist_name:
		:param is_2d:
		"""
		if hist_name in binning_definitions:
			# get binning
			binning = binning_definitions[hist_name]
			# validate binning
			if not is_2d and len(binning) != 3:
				raise Exception("1D Hist binning need to be (nbins, xlow, xhigh)")
			if is_2d and len(binning) != 6:
				raise Exception("2D Hist binning need to be (nbinsx, xlow, xhigh, nbinsy, ylow, yhigh)")
		else:
			binning = (2000, -1000, 1000) if not is_2d else (2000, -1000, 1000, 2000, -1000, 1000)
			self.logger.debug("No bins defined for histogram {}. Using defaults {}.".format(hist_name, binning))
		# TODO if the name has to be unique for ROOT to understand, use this line to include the directory. Still testing
		if not is_2d:
			nbins, xlow, xhigh = binning
			hist = ROOT.TH1D(directory + hist_name, hist_name, nbins, xlow, xhigh)
		else:
			nbinsx, xlow, xhigh, nbinsy, ylow, yhigh = binning
			hist = ROOT.TH1D(directory + hist_name, hist_name, nbinsx, xlow, xhigh, nbinsy, ylow, yhigh)
		hist.Sumw2()
		hist.SetDirectory(0)
		self.histogram_dict[directory + hist_name] = hist

	def write_histograms(self, root_file):
		"""
		Write all histograms to the specified root file.
		Root file must be open.
		:param root_file: a ROOT.TFile() object
		"""
		root_file.cd()

		for hist_path, histogram in self.histogram_dict.items():
			# split the identifier into directory and name
			tokens = hist_path.split('/')
			hist_name = tokens[-1]
			directory = '/'.join(tokens[:-1])
			# make the directory if it doesn't exist and cd to it
			if not root_file.GetDirectory(directory):
					root_file.mkdir(directory)	
			root_file.cd(directory)
			# write the histogram
			histogram.Write(hist_name)


binning_definitions = {
	"""
	Binning for 1d histograms is:
		(nbins, xlow, xhigh)
	Binning for 2d histograms is:
		(nbinsx, xlow, xhigh, nbinsy, ylow, yhigh)
	"""
	'nmuon': (50, 0, 50),
	'nel': (50, 0, 50),
	'muon_type': (6, -0.5, 5.5),
	'muon_pt': (1000, 0, 1000),
	'muon_eta': (40, -10, 10),
	'muon_phi': (16, -4, 4),
	'muon_quality': (4, -.5, 3.5),
	'el_pt': (1000, 0, 1000),
	'el_eta': (40, -10, 10),
	'el_phi': (16, -4, 4),
	'el_quality': (4, -.5, 3.5),
	'prompt_muon': (50, 0, 50),
	'prompt_electron': (50, 0, 50),
	'prompt_lepton': (50, 0, 50),
	'presel_num_trks': (21, -0.5, 20.5),
	'nonlep_pt': (1000, 0, 1000),
	'nonlep_eta': (40, -10, 10),
	'nonlep_phi': (16, -4, 4),
	'nonlep_d0': (2000, -1000, 1000),
	'plep_pt': (1000, 0, 1000),
	'plep_eta': (40, -10, 10),
	'plep_phi': (16, -4, 4),
	'ptrk_pt': (1000, 0, 1000),
	'ptrk_eta': (40, -10, 10),
	'ptrk_phi': (16, -4, 4),
	'largew_plep_pt': (1000, 0, 1000),
	'largew_plep_eta': (40, -10, 10),
	'largew_plep_phi': (16, -4, 4),
	'plep_d0': (2000, -1000, 1000),
	'plep_z0': (2000, -1000, 1000),
	'ptrk_d0': (2000, -1000, 1000),
	'ptrk_z0': (2000, -1000, 1000),
	'plep_charge': (3, -1.5, 1.5),
	'DV_weight': (1, -10, 10),
	'DV_cosmic_sep': (1000, 0, 10),
	'DV_alpha': (350, 0, 3.5),
	'DV_pass_mat_veto' : (2, -0.5, 1.5),
	'DV_truth_matched' : (2, -0.5, 1.5),
	'DV_x': (2000, -500, 500),
	'DV_y': (2000, -500, 500),
	'DV_z': (2000, -500, 500),
	'DV_r': (500, 0, 500),
	'PV_x': (2000, -500, 500),
	'PV_y': (2000, -500, 500),
	'PV_z': (2000, -500, 500),
	'DV_num_trks': (6, -0.5, 5.5),
	'DV_distFromPV': (500, 0, 500),
	'DV_mass': (4000, 0, 100),
	'DV_pt': (1000, 0, 1000),
	'DV_sum_track_pt': (1000, 0, 1000),
	'DV_sum_track_pt_wrt_pv': (1000, 0, 1000),
	'DV_sum_track_pt_diff': (10000, 0, 100),
	'DV_sum_track_charge': (4, -1.5, 1.5),
	'DV_eta': (160, -10, 10),
	'DV_phi': (64, -4, 4),
	'DV_charge': (11, -5.5, 5.5),
	'DV_chi2': (600, 0, 30),
	'DV_chi2_assoc': (600, 0, 30),
	'DV_alpha': (400, -4, 4),
	'mvis': (10000, 0, 5000),
	'HNLm': (10005, -5, 5000),
	'HNLm_altbinning': (15, 0, 30),
	'alt_HNLm': (10005, -5, 5000),
	'HNLm_fixWmass': (10005, -5, 5000),
	'HNLpt': (1000, 0, 1000),
	'HNLeta': (160, -10, 10),
	'HNLphi': (64, -4, 4),
	'Wminus_HNLpt': (1000, 0, 1000),
	'Wminus_HNLeta': (160, -10, 10),
	'Wminus_HNLphi': (64, -4, 4),
	'Wminus_HNLE': (4000, 0, 2000),
	'Wplus_HNLpt': (1000, 0, 1000),
	'Wplus_HNLeta': (160, -10, 10),
	'Wplus_HNLphi': (64, -4, 4),
	'Wplus_HNLE': (4000, 0, 2000),
	'mtrans': (10000, 0, 5000),
	'DV_redmass': (10005, -5, 5000),
	'DV_redmassvis': (10005, -5, 5000),
	'DV_redmassHNL': (10005, -5, 5000),
	'DV_2tight': (2, -0.5, 1.5),
	'DV_2medium': (2, -0.5, 1.5),
	'DV_2loose': (2, -0.5, 1.5),
	'DV_1tight': (2, -0.5, 1.5),
	'DV_1medium': (2, -0.5, 1.5),
	'DV_1loose': (2, -0.5, 1.5),
	'DV_tight_loose' : (2, -0.5, 1.5), 
	'DV_tight_medium' : (2, -0.5, 1.5), 
	'DV_medium_loose' : (2, -0.5, 1.5), 	
	'DV_tight_veryloose' : (2, -0.5, 1.5), 
	'DV_medium_veryloose' : (2, -0.5, 1.5), 
	'DV_loose_veryloose' : (2, -0.5, 1.5), 
	'DV_tight_veryveryloose' : (2, -0.5, 1.5), 
	'DV_medium_veryveryloose' : (2, -0.5, 1.5), 
	'DV_loose_veryveryloose' : (2, -0.5, 1.5), 
	'DV_2veryveryloose' : (2, -0.5, 1.5),
	'DV_1veryveryloose' : (2, -0.5, 1.5),
	'DV_mumu': (2, -0.5, 1.5),
	'DV_ee': (2, -0.5, 1.5),
	'DV_emu': (2, -0.5, 1.5),
	'DV_1lep': (2, -0.5, 1.5),
	'DV_trk_max_chi2_toSV': (600, 0, 30),
	'DV_trk_min_chi2_toSV': (600, 0, 30),
	'DV_trk_max_d0_wrtSV': (1000, -50, 50),
	'DV_trk_min_d0_wrtSV': (1000, -50, 50),
	'DV_trk_max_errd0_wrtSV': (10000, 0, 50),
	'DV_trk_min_errd0_wrtSV': (10000, 0, 50),
	'DV_trk_max_z0_wrtSV': (3000, -30, 30),
	'DV_trk_min_z0_wrtSV': (3000, -30, 30),
	'DV_trk_max_errz0_wrtSV': (10000, 0, 50),
	'DV_trk_min_errz0_wrtSV': (10000, 0, 50),
	'DV_max_dR': (1000, 0, 10),
	'DV_max_dR_wrtSV': (1000, 0, 10),
	'DV_maxOpAng': (2000, -1, 1),
	'DV_minOpAng': (2000, -1, 1),
	'DV_maxd0': (1000, -250, 250),
	'DV_mind0': (1000, -250, 250),
	'DV_ntrk': (3, -0.5, 2.5),
	'DV_ntrk_lrt': (3, -0.5, 2.5),
	'DV_ntrk_sel': (3, -0.5, 2.5),
	'DV_ntrk_assoc': (3, -0.5, 2.5),
	'DV_trk_dpt': (2000, 0, 1000),
	'DV_trk_deta': (160, 0, 20),
	'DV_trk_dphi': (64, 0, 8),
	'DV_trk_dR': (1010, -5, 10),
	'DV_vertexing_uncertainty': (100, -1, 1),
	'event_type_MCweight': (1000, -5, 5),
	'M2_spin_corr_MCweight': (1000, -0.02, 0.02),
	'M2_nocorr_MCweight': (100000, -5, 5),
	'W_mass': (10000, 0, 5000),
	'W_pt': (1000, 0, 1000),
	'W_eta': (40, -10, 10),
	'W_phi': (16, -4, 4),
	'HNL_mass': (10000, 0, 5000),
	'mHNLcalc': (10000, 0, 5000),
	'HNL_pt': (1000, 0, 1000),
	'HNL_eta': (40, -10, 10),
	'HNL_phi': (16, -4, 4),
	'plep_mass': (10000, 0, 5000),
	'DV_trk_pt': (1000, 0, 1000),
	'DV_trk_eta': (160, -10, 10),
	'DV_trk_phi': (64, -4, 4),
	'DV_trk_d0': (2000, -200, 200),
	'DV_trk_0_pt': (1000, 0, 1000),
	'DV_trk_0_eta': (160, -10, 10),
	'DV_trk_0_phi': (64, -4, 4),
	'DV_trk_0_d0': (2000, -200, 200),
	'DV_trk_0_z0': (3000, -1500, 1500),
	'DV_trk_0_charge': (11, -5.5, 5.5),
	'DV_trk_0_chi2': (600, 0, 30),
	'DV_trk_0_isSelected': (2, -0.5, 1.5),
	'DV_trk_0_isAssociated': (2, -0.5, 1.5),
	'DV_trk_0_mom_parall': (1000, 0, 500),
	'DV_trk_0_mom_perp': (2000, -500, 500),
	'DV_trk_0_mom_mag': (2000, -500, 500),
	'DV_trk_0_mom_frac_parall': (1000, 0, 1),
	'DV_trk_1_pt': (1000, 0, 1000),
	'DV_trk_1_eta': (160, -10, 10),
	'DV_trk_1_phi': (64, -4, 4),
	'DV_trk_1_d0': (2000, -200, 200),
	'DV_trk_1_z0': (3000, -1500, 1500),
	'DV_trk_1_charge': (11, -5.5, 5.5),
	'DV_trk_1_chi2': (600, 0, 30),
	'DV_trk_1_isSelected': (2, -0.5, 1.5),
	'DV_trk_1_isAssociated': (2, -0.5, 1.5),
	'DV_trk_1_mom_parall': (1000, 0, 500),
	'DV_trk_1_mom_perp': (1000, 0, 500),
	'DV_trk_1_mom_mag': (1000, 0, 500),
	'DV_trk_1_mom_frac_parall': (1000, 0, 1),
	'DV_d0_cut': (2, -0.5, 1.5),
	'DV_El_pt': (1000, 0, 1000),
	'DV_El_eta': (160, -10, 10),
	'DV_El_phi': (64, -4, 4),
	'DV_Mu_pt': (1000, 0, 1000),
	'DV_Mu_eta': (160, -10, 10),
	'DV_Mu_phi': (64, -4, 4),
	'lep1_trk_pt': (1000, 0, 1000),
	'lep1_trk_eta': (160, -10, 10),
	'lep1_trk_phi': (64, -4, 4),
	'lep2_trk_pt': (1000, 0, 1000),
	'lep2_trk_eta': (160, -10, 10),
	'lep2_trk_phi': (64, -4, 4),
	'nu_trk_pt': (1000, 0, 1000),
	'nu_trk_eta': (160, -10, 10),
	'nu_trk_phi': (64, -4, 4),
	'largew_lep1_trk_pt': (1000, 0, 1000),
	'largew_lep1_trk_eta': (160, -10, 10),
	'largew_lep1_trk_phi': (64, -4, 4),
	'largew_lep2_trk_pt': (1000, 0, 1000),
	'largew_lep2_trk_eta': (160, -10, 10),
	'largew_lep2_trk_phi': (64, -4, 4),
	'largew_nu_trk_pt': (1000, 0, 1000),
	'largew_nu_trk_eta': (160, -10, 10),
	'largew_nu_trk_phi': (64, -4, 4),
	'dlep1_pt': (1000, 0, 1000),
	'dlep1_eta': (160, -10, 10),
	'dlep1_phi': (64, -4, 4),
	'dlep2_pt': (1000, 0, 1000),
	'dlep2_eta': (160, -10, 10),
	'dlep2_phi': (64, -4, 4),
	'dlep3_pt': (1000, 0, 1000),
	'dlep3_eta': (160, -10, 10),
	'dlep3_phi': (64, -4, 4),
	
	'm12': (200, 0, 100),
	'm23': (200, 0, 100),
	'm13': (200, 0, 100),
	'm12_sq': (200, 0, 100),
	'm23_sq': (200, 0, 100),
	'm13_sq': (200, 0, 100),
	's23': (200, 0, 100),
	's24': (200, 0, 100),
	's34': (200, 0, 100),
	's12': (12000, 0, 6000),
	's13': (12000, 0, 6000),
	's14': (12000, 0, 6000),
	'maxlinkTruth_score': (1000, 0, 1),
	'maxlinkTruth_parent_pdgId': (201, -0.5, 200),
	'num_trks': (21, -0.5, 20.5),
	'muon_isAssociated': (2, -0.5, 1.5),
	'bothmuon_isAssociated': (2, -0.5, 1.5),
	'nomuon_isAssociated': (2, -0.5, 1.5),
	'onemuon_isAssociated': (2, -0.5, 1.5),
	'onetrk_isAssociated': (2, -0.5, 1.5),
	'nonlep_trk_pt': (1000, 0, 1000),
	'nonlep_trk_eta': (160, -10, 10),
	'nonlep_trk_phi': (64, -4, 4),
	'nonlep_trk_d0': (1000, -250, 250),
	'nonlep_trk_z0': (3000, -1500, 1500),
	'nonlep_trk_absz0': (1500, 0, 1500),
	'nonlep_trk_charge': (11, -5.5, 5.5),
	'nonlep_trk_chi2': (600, 0, 30),
	'nonlep_trk_isLRT': (2, -0.5, 1.5),
	'nonlep_trk_isSelected': (2, -0.5, 1.5),
	'nonlep_trk_isAssociated': (2, -0.5, 1.5),
	'nonlep_trk_nPixelHits': (15, -0.5, 14.5),
	'nonlep_trk_nSCTHits': (15, -0.5, 14.5),
	'nonlep_trk_nSCTHoles': (4, -0.5, 3.5),
	'nonlep_trk_nSiHits': (24, -0.5, 23.5),
	'nonlep_trk_dTheta': (500, 0, 250),
	'nonlep_trk_chi2_toSV': (600, 0, 30),
	'nonlep_trk_d0_wrtSV': (1000, -50, 50),
	'nonlep_trk_errd0_wrtSV': (10000, 0, 50),
	'nonlep_trk_z0_wrtSV': (3000, -30, 30),
	'nonlep_trk_errz0_wrtSV': (10000, 0, 50),
	'lep_trk_pt': (1000, 0, 1000),
	'lep_trk_eta': (160, -10, 10),
	'lep_trk_phi': (64, -4, 4),
	'lep_trk_d0': (1000, -250, 250),
	'lep_trk_z0': (3000, -1500, 1500),
	'lep_trk_absz0': (1500, 0, 1500),
	'lep_trk_charge': (11, -5.5, 5.5),
	'lep_trk_chi2': (600, 0, 30),
	'lep_trk_isLRT': (2, -0.5, 1.5),
	'lep_trk_isSelected': (2, -0.5, 1.5),
	'lep_trk_isAssociated': (2, -0.5, 1.5),
	'lep_trk_nPixelHits': (15, -0.5, 14.5),
	'lep_trk_nSCTHits': (15, -0.5, 14.5),
	'lep_trk_nSCTHoles': (4, -0.5, 3.5),
	'lep_trk_nSiHits': (24, -0.5, 23.5),
	'lep_trk_dTheta': (500, 0, 250),
	'lep_trk_chi2_toSV': (600, 0, 30),
	'lep_trk_d0_wrtSV': (1000, -50, 50),
	'lep_trk_errd0_wrtSV': (10000, 0, 50),
	'lep_trk_z0_wrtSV': (3000, -30, 30),
	'lep_trk_errz0_wrtSV': (10000, 0, 50),
	'all_trk_pt': (1000, 0, 1000),
	'all_trk_eta': (160, -10, 10),
	'all_trk_phi': (64, -4, 4),
	'all_trk_d0': (1000, -250, 250),
	'all_trk_z0': (3000, -1500, 1500),
	'all_trk_absz0': (1500, 0, 1500),
	'all_trk_charge': (11, -5.5, 5.5),
	'all_trk_chi2': (600, 0, 30),
	'all_trk_isLRT': (2, -0.5, 1.5),
	'all_trk_isSelected': (2, -0.5, 1.5),
	'all_trk_isAssociated': (2, -0.5, 1.5),
	'all_trk_nPixelHits': (15, -0.5, 14.5),
	'all_trk_nSCTHits': (15, -0.5, 14.5),
	'all_trk_nSCTHoles': (4, -0.5, 3.5),
	'all_trk_nSiHits': (24, -0.5, 23.5),
	'all_trk_dTheta': (500, 0, 250),
	'all_trk_chi2_toSV': (600, 0, 30),
	'all_trk_d0_wrtSV': (1000, -50, 50),
	'all_trk_errd0_wrtSV': (10000, 0, 50),
	'all_trk_z0_wrtSV': (3000, -30, 30),
	'all_trk_errz0_wrtSV': (10000, 0, 50),

	# 2d histograms
	'charge_ntrk': (11, -5.5, 5.5, 9, -0.5, 8.5),
	'DVmass_mvis': (1000, 0, 500, 1000, 0, 500),
	'DVmass_mhnl': (1000, 0, 500, 1010, -5, 500),
	'DVmass_mtrans': (1000, 0, 500, 1000, 0, 500),
	'DVmass_hnlpt': (1000, 0, 500, 1010, -5, 500),
	'mvis_mhnl': (1000, 0, 500, 1010, -5, 500),
	'mvis_mtrans': (1000, 0, 500, 1000, 0, 500),
	'mvis_hnlpt': (1000, 0, 500, 1010, -5, 500),
	'mhnl_hnlpt': (1010, -5, 500, 1010, -5, 500),
	'mhnl_mtrans': (1010, -5, 500, 1000, 0, 500),
	'mhnl2D': (1010, -5, 500, 1000, 0, 500),
}

# add track binnings
for trk_str in ['', '_0', '_1']:
	binning_definitions['DV_trk{}_pt'.format(trk_str)] = (1000, 0, 1000)
	binning_definitions['DV_trk{}_eta'.format(trk_str)] = (160, -10, 10)
	binning_definitions['DV_trk{}_phi'.format(trk_str)] = (64, -4, 4)
	binning_definitions['DV_trk{}_d0'.format(trk_str)] = (1000, -250, 250)
	binning_definitions['DV_trk{}_z0'.format(trk_str)] = (3000, -1500, 1500)
	binning_definitions['DV_trk{}_absz0'.format(trk_str)] = (1500, 0, 1500)
	binning_definitions['DV_trk{}_charge'.format(trk_str)] = (11, -5.5, 5.5)
	binning_definitions['DV_trk{}_chi2'.format(trk_str)] = (600, 0, 30)
	binning_definitions['DV_trk{}_isLRT'.format(trk_str)] = (2, -0.5, 1.5)
	binning_definitions['DV_trk{}_isSelected'.format(trk_str)] = (2, -0.5, 1.5)
	binning_definitions['DV_trk{}_isAssociated'.format(trk_str)] = (2, -0.5, 1.5)
	binning_definitions['DV_trk{}_nPixelHits'.format(trk_str)] = (15, -0.5, 14.5)
	binning_definitions['DV_trk{}_nSCTHits'.format(trk_str)] = (15, -0.5, 14.5)
	binning_definitions['DV_trk{}_nSCTHoles'.format(trk_str)] = (4, -0.5, 3.5)
	binning_definitions['DV_trk{}_nSiHits'.format(trk_str)] = (24, -0.5, 23.5)
	binning_definitions['DV_trk{}_dTheta'.format(trk_str)] = (500, 0, 250)
	binning_definitions['DV_trk{}_chi2_toSV'.format(trk_str)] = (600, 0, 30)
	binning_definitions['DV_trk{}_d0_wrtSV'.format(trk_str)] = (1000, -50, 50)
	binning_definitions['DV_trk{}_errd0_wrtSV'.format(trk_str)] = (10000, 0, 50)
	binning_definitions['DV_trk{}_z0_wrtSV'.format(trk_str)] = (3000, -30, 30)
	binning_definitions['DV_trk{}_errz0_wrtSV'.format(trk_str)] = (10000, 0, 50)

for i in range(8):
	binning_definitions['num_assoc_{}trk'.format(i)] = (21, -0.5, 20.5)
