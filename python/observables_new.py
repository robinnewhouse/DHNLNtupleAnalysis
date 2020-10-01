import ROOT
import helpers
ObservableList = []


class Observables():
	"""
	A class to contain all of the histograms in a directory
	"""

	def __init__(self):
		self.histogram_dict = {}  # stores a collection of Observable objects
		self.logger = helpers.getLogger('dHNLAnalysis.observables')

	def fill_hist(self, directory, hist_name, variable_1, variable_2):
		"""
		This function will fill a histogram that is stored in the histogram_dict, and if it does not exist, will create it.
		:param directory: The path where the histogram will be stored in the root file.
		:param hist_name: The name of the histogram
		:param variable_1: Must be a number
		:param variable_2: Optional, will fill a 2d histogram
		"""
		if not directory + hist_name in self.histogram_dict:
			self.add_histogram(directory, hist_name)

	def add_histogram(self, directory, hist_name):
		"""
		Create a new histogram. If a defined binning exists, use that. Otherwise use defaults.
		:param directory:
		:param hist_name:
		:return:
		"""
		if hist_name in binning_definitions:
			binning = binning_definitions[hist_name]
			if len(binning) != 3:
				raise Exception("1D Hist binning need to be (nbins, xlow, xhigh).")
			nbins, xlow, xhigh = binning
		else:
			self.logger.warning("No bins defined for 1D histogram. Using defaults.")
			nbins, xlow, xhigh = (2000, -1000, 1000)
		# TODO if the name has to be unique for ROOT to understand, use this line to include the directory. Still testing
		# hist = ROOT.TH1D(directory + hist_name, hist_name, nBins, xLow, xHigh)
		hist = ROOT.TH1D(hist_name, hist_name, nbins, xlow, xhigh)


binning_definitions = {
	"""
	Binning for 1d histograms is:
		(nbins, xlow, xhigh)
	Binning for 2d histograms is:
		(nbinsx, xlow, xhigh, nbinsy, ylow, yhigh)
	"""
	"all_nmuon", (50, 0, 50),
	"all_nel", (50, 0, 50),
	"all_muon_type", (6, -0.5, 5.5),
	"all_muon_pt", (1000, 0, 1000),
	"all_muon_eta", (40, -10, 10),
	"all_muon_phi", (16, -4, 4),
	"all_muon_quality", (4, -.5, 3.5),
	"all_el_pt", (1000, 0, 1000),
	"all_el_eta", (40, -10, 10),
	"all_el_phi", (16, -4, 4),
	"all_el_quality", (4, -.5, 3.5),
	"all_prompt_muon", (50, 0, 50),
	"all_prompt_electron", (50, 0, 50),
	"all_prompt_lepton", (50, 0, 50),
	"presel_num_trks", (21, -0.5, 20.5),
}

def queue_all_observables(analysis_name,is_data):
	# Observable("sel_plep_z0", do=['hist']).queue() # what is this histogram? -DT



	reco_histograms("all")
	if is_data:
		EventType = [""]
	else:
		EventType = ["_LNC","_LNV"]
	for i in range(len(EventType)):
		if "oldAnalysis" or "ToyAnalysis" in analysis_name:
			reco_histograms("presel" + EventType[i])
			reco_histograms("DVtype" + EventType[i])
			reco_histograms("cosmic" + EventType[i])
			reco_histograms("mlll" + EventType[i])
			reco_histograms("sel" + EventType[i])
			reco_histograms("trkqual" + EventType[i])
			reco_histograms("match" + EventType[i])
		if analysis_name == "KShort":
			reco_histograms("mass" + EventType[i])
			reco_histograms("alpha" + EventType[i])
		if analysis_name == "ToyAnalysis":
			reco_histograms("mDV" + EventType[i])

		truth_histograms("truth_all")
		truth_histograms("truth_LNC")
		truth_histograms("truth_LNV")
		truth_histograms("truth_LNC_presel")
		# truth_histograms("truth_LNC_mass")
		# truth_histograms("truth_LNC_charge")
		# truth_histograms("truth_LNC_DVtype")
		# truth_histograms("truth_LNC_trkqual")
		# truth_histograms("truth_LNC_cosmic")
		# truth_histograms("truth_LNC_mlll")
		# truth_histograms("truth_LNC_mDV")
		# truth_histograms("truth_LNC_sel")
		# truth_histograms("truth_LNC_match")
		truth_histograms("truth_LNV_presel")
		# truth_histograms("truth_LNV_mass")
		# truth_histograms("truth_LNV_charge")
		# truth_histograms("truth_LNV_DVtype")
		# truth_histograms("truth_LNV_trkqual")
		# truth_histograms("truth_LNV_cosmic")
		# truth_histograms("truth_LNV_mlll")
		# truth_histograms("truth_LNV_mDV")
		# truth_histograms("truth_LNV_sel")
		# truth_histograms("truth_LNV_match")

	if analysis_name == "ToyAnalysis":
		twolep_multi()


def reco_histograms(selection):
	# does the histogram need truth?
	need_truth = "LNC" or "LNV" in selection

	Observable(selection + "_nonlep_pt", binning=(1000, 0, 1000), do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_nonlep_eta", binning=(40, -10, 10), do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_nonlep_phi", binning=(16, -4, 4), do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_nonlep_d0", binning=(2000, -1000, 1000), do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_plep_pt", binning=(1000, 0, 1000), do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_plep_eta", binning=(40, -10, 10), do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_plep_phi", binning=(16, -4, 4), do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_plep_d0", binning=(2000, -1000, 1000), do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_plep_z0", binning=(2000, -1000, 1000), do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_plep_charge", binning=(3,-1.5,1.5), do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_weight", binning=(1, -10, 10), do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()  # dummy histogram
	Observable(selection + "_DV_x",binning = (2000,-500,500), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_y",binning = (2000,-500,500), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_z",binning = (2000,-500,500), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_r",binning = (500,0,500), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_num_trks",binning = (6,-0.5,5.5), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_presel_num_trks",binning = (21,-0.5,20.5), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_distFromPV",binning = (500,0,500), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_mass", binning=(2000, 0, 1000), do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_pt",binning = (1000,0,1000), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_sum_track_pt",binning = (1000,0,1000), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_sum_track_pt_wrt_pv",binning = (1000,0,1000), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_sum_track_pt_diff",binning = (10000,0,100), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_sum_track_charge",binning = (4,-1.5,1.5), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_eta",binning = (160,-10,10),do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_phi", binning = (64,-4,4), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_charge",binning = (11,-5.5,5.5), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_chi2",binning = (600,0,30), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_chi2_assoc",binning = (600,0,30), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_alpha",binning = (400,-4,4), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_mvis",binning = (10000,0,5000), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_HNLm",binning = (10005,-5,5000), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_HNLm_altbinning",binning = (15,0,30), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_alt_HNLm",binning = (10005,-5,5000), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_HNLm_fixWmass",binning = (10005,-5,5000), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_HNLpt",binning = (1000,0,1000), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_HNLeta",binning = (160,-10,10),do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_HNLphi", binning = (64,-4,4), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_mtrans",binning = (10000,0,5000), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_redmass",binning = (10005,-5,5000), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_redmassvis",binning = (10005,-5,5000), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_redmassHNL",binning = (10005,-5,5000), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_2tight",binning = (2,-0.5,1.5), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_2medium",binning = (2,-0.5,1.5), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_2loose",binning = (2,-0.5,1.5), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_1tight",binning = (2,-0.5,1.5), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_1medium",binning = (2,-0.5,1.5), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_1loose",binning = (2,-0.5,1.5), do = ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_trk_max_chi2_toSV", binning = (600,0,30), do= ['ToyAnalysis'],need_truth = need_truth).queue()
	Observable(selection + "_DV_trk_min_chi2_toSV", binning = (600,0,30), do= ['ToyAnalysis'],need_truth = need_truth).queue()
	Observable(selection + "_DV_trk_max_d0_wrtSV", binning = (1000,-50,50), do= ['ToyAnalysis'],need_truth = need_truth).queue()
	Observable(selection + "_DV_trk_min_d0_wrtSV", binning = (1000,-50,50), do= ['ToyAnalysis'],need_truth = need_truth).queue()
	Observable(selection + "_DV_trk_max_errd0_wrtSV", binning = (10000,0,50), do= ['ToyAnalysis'],need_truth = need_truth).queue()
	Observable(selection + "_DV_trk_min_errd0_wrtSV", binning = (10000,0,50), do= ['ToyAnalysis'],need_truth = need_truth).queue()
	Observable(selection + "_DV_trk_max_z0_wrtSV", binning = (3000,-30,30), do= ['ToyAnalysis'],need_truth = need_truth).queue()
	Observable(selection + "_DV_trk_min_z0_wrtSV", binning = (3000,-30,30), do= ['ToyAnalysis'],need_truth = need_truth).queue()
	Observable(selection + "_DV_trk_max_errz0_wrtSV", binning = (10000,0,50), do= ['ToyAnalysis'],need_truth = need_truth).queue()
	Observable(selection + "_DV_trk_min_errz0_wrtSV", binning = (10000,0,50), do= ['ToyAnalysis'],need_truth = need_truth).queue()
	Observable(selection + "_DV_max_dR", binning = (1000,0,10),do=['ToyAnalysis'],need_truth = need_truth).queue()
	Observable(selection + "_DV_max_dR_wrtSV", binning = (1000,0,10),do=['ToyAnalysis'],need_truth = need_truth).queue()
	Observable(selection + "_DV_maxOpAng", binning = (2000,-1,1),do=['ToyAnalysis'],need_truth = need_truth).queue()
	Observable(selection + "_DV_minOpAng", binning = (2000,-1,1),do=['ToyAnalysis'],need_truth = need_truth).queue()
	Observable(selection + "_DV_maxd0", binning = (1000,-250,250),do=['ToyAnalysis'],need_truth = need_truth).queue()
	Observable(selection + "_DV_mind0", binning = (1000,-250,250),do=['ToyAnalysis'],need_truth = need_truth).queue()
	Observable(selection + "_DV_ntrk", binning = (3,-0.5,2.5),do=['ToyAnalysis'],need_truth = need_truth).queue()
	Observable(selection + "_DV_ntrk_lrt", binning = (3,-0.5,2.5),do=['ToyAnalysis'],need_truth = need_truth).queue()
	Observable(selection + "_DV_ntrk_sel", binning = (3,-0.5,2.5),do=['ToyAnalysis'],need_truth = need_truth).queue()
	Observable(selection + "_DV_ntrk_assoc", binning = (3,-0.5,2.5),do=['ToyAnalysis'],need_truth = need_truth).queue()

	trk_str = ["","_0","_1"]
	for i in range(len(trk_str)):
		Observable(selection + "_DV_trk{}_pt".format(trk_str[i]),binning = (1000,0,1000), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
		Observable(selection + "_DV_trk{}_eta".format(trk_str[i]),binning = (160,-10,10),do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
		Observable(selection + "_DV_trk{}_phi".format(trk_str[i]), binning =(64,-4,4), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
		Observable(selection + "_DV_trk{}_d0".format(trk_str[i]),binning = (1000,-250,250), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
		Observable(selection + "_DV_trk{}_z0".format(trk_str[i]),binning = (3000,-1500,1500), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
		Observable(selection + "_DV_trk{}_absz0".format(trk_str[i]),binning = (1500,0,1500), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
		Observable(selection + "_DV_trk{}_charge".format(trk_str[i]),binning = (11,-5.5,5.5), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
		Observable(selection + "_DV_trk{}_chi2".format(trk_str[i]),binning = (600,0,30), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
		Observable(selection + "_DV_trk{}_isLRT".format(trk_str[i]),binning = (2,-0.5,1.5), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
		Observable(selection + "_DV_trk{}_isSelected".format(trk_str[i]),binning = (2,-0.5,1.5), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
		Observable(selection + "_DV_trk{}_isAssociated".format(trk_str[i]),binning = (2,-0.5,1.5), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
		Observable(selection + "_DV_trk{}_nPixelHits".format(trk_str[i]),binning = (15,-0.5,14.5), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
		Observable(selection + "_DV_trk{}_nSCTHits".format(trk_str[i]),binning = (15,-0.5,14.5), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
		Observable(selection + "_DV_trk{}_nSCTHoles".format(trk_str[i]),binning = (4,-0.5,3.5), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
		Observable(selection + "_DV_trk{}_nSiHits".format(trk_str[i]),binning = (24,-0.5,23.5), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
		Observable(selection + "_DV_trk{}_dTheta".format(trk_str[i]),binning = (500,0,250), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
		Observable(selection + "_DV_trk{}_chi2_toSV".format(trk_str[i]),binning = (600,0,30), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
		Observable(selection + "_DV_trk{}_d0_wrtSV".format(trk_str[i]),binning = (1000,-50,50), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
		Observable(selection + "_DV_trk{}_errd0_wrtSV".format(trk_str[i]),binning = (10000,0,50), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
		Observable(selection + "_DV_trk{}_z0_wrtSV".format(trk_str[i]),binning = (3000,-30,30), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
		Observable(selection + "_DV_trk{}_errz0_wrtSV".format(trk_str[i]),binning = (10000,0,50), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()

	Observable(selection + "_DV_trk_dpt",binning = (2000,0,1000), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
	Observable(selection + "_DV_trk_deta",binning = (160,0,20),do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
	Observable(selection + "_DV_trk_dphi", binning = (64,0,8), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()
	Observable(selection + "_DV_trk_dR", binning = (1010,-5,10), do = ['oldAnalysis','ToyAnalysis','KShort']).queue()


def truth_histograms(selection):
	# histograms after selection selection
	Observable( selection + "_event_type_MCweight",binning = (1000,-5,5), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_M2_spin_corr_MCweight",binning = (1000,-0.02,0.02), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_M2_nocorr_MCweight",binning = (100000,-5,5), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_W_mass",binning = (10000,0,5000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_W_pt",binning = (1000,0,1000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_W_eta",binning = (40,-10,10),do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_W_phi", binning = (16,-4,4), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_HNL_mass",binning = (10000,0,5000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_mHNLcalc",binning = (10000,0,5000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_HNL_pt",binning = (1000,0,1000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_HNL_eta",binning = (40,-10,10),do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_HNL_phi", binning = (16,-4,4), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()

	Observable( selection + "_plep_pt",binning = (1000,0,1000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_plep_eta",binning = (40,-10,10),do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_plep_phi", binning = (16,-4,4), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_plep_mass",binning = (10000,0,5000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()

	# Observable( selection + "_plep_d0", do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue() # no d0 or z0 for truth particles
	# Observable( selection + "_plep_z0", do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()

	Observable( selection + "_DV_trk_pt",binning = (1000,0,1000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_DV_trk_eta",binning = (160,-10,10),do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_DV_trk_phi", binning =(64,-4,4), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_lep1_trk_pt", binning=(1000, 0, 1000), do=['oldAnalysis','ToyAnalysis'], need_truth=True).queue()
	Observable( selection + "_lep1_trk_eta", binning=(160, -10, 10), do=['oldAnalysis','ToyAnalysis'], need_truth=True).queue()
	Observable( selection + "_lep1_trk_phi", binning=(64, -4, 4), do=['oldAnalysis','ToyAnalysis'], need_truth=True).queue()
	Observable( selection + "_lep2_trk_pt", binning=(1000, 0, 1000), do=['oldAnalysis','ToyAnalysis'], need_truth=True).queue()
	Observable( selection + "_lep2_trk_eta", binning=(160, -10, 10), do=['oldAnalysis','ToyAnalysis'], need_truth=True).queue()
	Observable( selection + "_lep2_trk_phi", binning=(64, -4, 4), do=['oldAnalysis','ToyAnalysis'], need_truth=True).queue()
	Observable( selection + "_nu_trk_pt", binning=(1000, 0, 1000), do=['oldAnalysis','ToyAnalysis'], need_truth=True).queue()
	Observable( selection + "_nu_trk_eta", binning=(160, -10, 10), do=['oldAnalysis','ToyAnalysis'], need_truth=True).queue()
	Observable( selection + "_nu_trk_phi", binning=(64, -4, 4), do=['oldAnalysis','ToyAnalysis'], need_truth=True).queue()
	# Observable( selection + "_DV_trk_d0", do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue() # no d0 or z0 for truth particles
	# Observable( selection + "_DV_trk_z0", do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_DV_trk_charge",binning = (11,-5.5,5.5), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_DV_trk_chi2",binning = (40,0,20), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_DV_trk_dpt",binning = (2000,0,1000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_DV_trk_deta",binning = (160,0,20),do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_DV_trk_dphi", binning = (64,0,8), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_DV_trk_dR", binning = (1010,-5,10), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()

	Observable( selection + "_DV_x",binning = (2000,-500,500), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_DV_mass",binning = (200,0,100), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_m12",binning = (200,0,100), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_m23",binning = (200,0,100), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_m13",binning = (200,0,100), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_m12_sq",binning = (200,0,100), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_m23_sq",binning = (200,0,100), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_m13_sq",binning = (200,0,100), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_s23",binning = (200,0,100), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_s24",binning = (200,0,100), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_s34",binning = (200,0,100), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_s12",binning = (12000,0,6000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_s13",binning = (12000,0,6000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_s14",binning = (12000,0,6000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_m13",binning = (200,0,100), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_DV_y",binning = (2000,-500,500), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_DV_z",binning = (2000,-500,500), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_DV_r",binning = (500,0,500), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()

	Observable( selection + "_maxlinkTruth_score",binning = (1000,0,1), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	Observable( selection + "_maxlinkTruth_parent_pdgId",binning = (201,-0.5,200), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_DV_num_trks",binning = (6,-0.5,5.5), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_DV_distFromPV",binning = (500,0,500), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_DV_mass",binning = (10000,0,5000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_DV_pt",binning = (1000,0,1000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_DV_eta",binning = (160,-10,10),do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_DV_phi", binning = (64,-4,4), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_DV_minOpAng",binning = (100,0,1), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_DV_maxOpAng",binning = (100,0,1), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_DV_charge",binning = (11,-5.5,5.5), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_DV_chi2",binning = (40,0,20), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_mvis",binning = (10000,0,5000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_HNLm",binning = (10005,-5,5000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_HNLm2",binning = (10005,-5,5000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_HNLpt",binning = (1000,0,1000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_HNLeta",binning = (160,-10,10),do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_HNLphi", binning = (64,-4,4), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_mtrans",binning = (10000,0,5000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_mtrans_rot",binning = (10000,0,5000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_DV_redmass",binning = (10005,-5,5000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_DV_redmassvis",binning = (10005,-5,5000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()
	# Observable( selection + "_DV_redmassHNL",binning = (10005,-5,5000), do = ['oldAnalysis','ToyAnalysis'], need_truth = True).queue()

def twolep_multi():
	Observable("2lepMultitrk_num_trks",binning = (21,-0.5,20.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_muon_isAssociated",binning = (2,-0.5,1.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_bothmuon_isAssociated",binning = (2,-0.5,1.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nomuon_isAssociated",binning = (2,-0.5,1.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_onemuon_isAssociated",binning = (2,-0.5,1.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_onetrk_isAssociated",binning = (2,-0.5,1.5), do = ['ToyAnalysis']).queue()

	Observable("2lepMultitrk_nonlep_trk_pt", binning = (1000,0,1000), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nonlep_trk_eta", binning = (160,-10,10),do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nonlep_trk_phi",  binning =(64,-4,4), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nonlep_trk_d0", binning = (1000,-250,250), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nonlep_trk_z0", binning = (3000,-1500,1500), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nonlep_trk_absz0", binning = (1500,0,1500), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nonlep_trk_charge", binning = (11,-5.5,5.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nonlep_trk_chi2", binning = (600,0,30), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nonlep_trk_isLRT", binning = (2,-0.5,1.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nonlep_trk_isSelected", binning = (2,-0.5,1.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nonlep_trk_isAssociated", binning = (2,-0.5,1.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nonlep_trk_nPixelHits", binning = (15,-0.5,14.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nonlep_trk_nSCTHits", binning = (15,-0.5,14.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nonlep_trk_nSCTHoles", binning = (4,-0.5,3.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nonlep_trk_nSiHits", binning = (24,-0.5,23.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nonlep_trk_dTheta", binning = (500,0,250), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nonlep_trk_chi2_toSV", binning = (600,0,30), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nonlep_trk_d0_wrtSV", binning = (1000,-50,50), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nonlep_trk_errd0_wrtSV", binning = (10000,0,50), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nonlep_trk_z0_wrtSV", binning = (3000,-30,30), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_nonlep_trk_errz0_wrtSV", binning = (10000,0,50), do = ['ToyAnalysis']).queue()

	Observable("2lepMultitrk_lep_trk_pt", binning = (1000,0,1000), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_lep_trk_eta", binning = (160,-10,10),do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_lep_trk_phi",  binning =(64,-4,4), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_lep_trk_d0", binning = (1000,-250,250), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_lep_trk_z0", binning = (3000,-1500,1500), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_lep_trk_absz0", binning = (1500,0,1500), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_lep_trk_charge", binning = (11,-5.5,5.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_lep_trk_chi2", binning = (600,0,30), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_lep_trk_isLRT", binning = (2,-0.5,1.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_lep_trk_isSelected", binning = (2,-0.5,1.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_lep_trk_isAssociated", binning = (2,-0.5,1.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_lep_trk_nPixelHits", binning = (15,-0.5,14.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_lep_trk_nSCTHits", binning = (15,-0.5,14.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_lep_trk_nSCTHoles", binning = (4,-0.5,3.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_lep_trk_nSiHits", binning = (24,-0.5,23.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_lep_trk_dTheta", binning = (500,0,250), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_lep_trk_chi2_toSV", binning = (600,0,30), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_lep_trk_d0_wrtSV", binning = (1000,-50,50), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_lep_trk_errd0_wrtSV", binning = (10000,0,50), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_lep_trk_z0_wrtSV", binning = (3000,-30,30), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_lep_trk_errz0_wrtSV", binning = (10000,0,50), do = ['ToyAnalysis']).queue()

	Observable("2lepMultitrk_all_trk_pt", binning = (1000,0,1000), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_all_trk_eta", binning = (160,-10,10),do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_all_trk_phi",  binning =(64,-4,4), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_all_trk_d0", binning = (1000,-250,250), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_all_trk_z0", binning = (3000,-1500,1500), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_all_trk_absz0", binning = (1500,0,1500), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_all_trk_charge", binning = (11,-5.5,5.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_all_trk_chi2", binning = (600,0,30), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_all_trk_isLRT", binning = (2,-0.5,1.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_all_trk_isSelected", binning = (2,-0.5,1.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_all_trk_isAssociated", binning = (2,-0.5,1.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_all_trk_nPixelHits", binning = (15,-0.5,14.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_all_trk_nSCTHits", binning = (15,-0.5,14.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_all_trk_nSCTHoles", binning = (4,-0.5,3.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_all_trk_nSiHits", binning = (24,-0.5,23.5), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_all_trk_dTheta", binning = (500,0,250), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_all_trk_chi2_toSV", binning = (600,0,30), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_all_trk_d0_wrtSV", binning = (1000,-50,50), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_all_trk_errd0_wrtSV", binning = (10000,0,50), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_all_trk_z0_wrtSV", binning = (3000,-30,30), do = ['ToyAnalysis']).queue()
	Observable("2lepMultitrk_all_trk_errz0_wrtSV", binning = (10000,0,50), do = ['ToyAnalysis']).queue()

	ntrk = [3,4,5,6,7,8,9,10]
	for i in range(len(ntrk)):
		Observable("2lepMultitrk_num_assoc_{}trk".format(i), binning = (21,-0.5,20.5), do = ['hist']).queue()
