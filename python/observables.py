import math
import ROOT
ObservableList = []


class Observable(object):
	def __init__(self, name, binning = (2000, -1000, 1000), script = None, style = 'single', do = ['hist'], only = None, title = '{self.name}', dtype = float, default = -1, need_truth = False):
		"""Observable for filling histograms or trees

			Parameters
			----------
			name : {str}
				Name of the histogram/branch
			binning : {tuple, array-like}, optional
				If bining is an tuple, it defines equal-width bins in the given range. If bins is an array-like, it defines the bin edges, including the rightmost edge, allowing for non-uniform bin widths.
				(the default is (20, -1000, 1000))
			script : {str}, optional
				[description] (the default is None, which does nothing)
			style : {str}, optional
				This should be either "foreach" or "single". Decide how you fill the hist/branch (the default is 'single')
			do : {list[str]}, optional
				This should be either or both "hist" and "tree"
			only : {list[str]}, optional
				This should be either or both "be", "bmu", "re" or "rmu" (the default is None, which means fill it for all channels)
			title : {str}, optional
				Title of the histogram (the default is '{self.name}')
		"""
		self.name = name
		self.binning = binning
		self.title = title.format(self = self)
		# self.script = compile(script, '<TopNtupleAnalysis.Observable>', 'eval')
		self.script = script
		self.style = style
		self.only = only
		self.do = do
		self.dtype = dtype
		self.default = default
		self.need_truth = need_truth

	def registered(self, analysis):
		return  _Observable(self.name, analysis, binning = self.binning, script = self.script, style = self.style, only = self.only, title = self.title, do = self.do, dtype = self.dtype, default = self.default, need_truth = self.need_truth)

	def queue(self):
		ObservableList.append(self)

class _Observable(object):
    def __init__(self, name, analysis, binning = (20, -1000, 1000), script = None, style = 'single', do = ['hist'], only = None, title = '{self.name}', dtype = None, default = -1, need_truth = False):
        self.name = name
        self.binning = binning
        self.title = title.format(self = self)
        self.script = script
        self.style = style
        self.only = only
        self.analysis = analysis
        self.do = do
        self.dtype = dtype
        self.default = default
        self.need_truth = need_truth
        # self._globals = {'analysis': self.analysis}
        # self._globals.update(self.analysis.run.im_func.func_globals)
        # self._globals.update(globals())
    # def __call__(self, _type = None, _locals = None):
    #     if self.script != None:
    #         if _locals != None:
    #             self._globals.update(_locals)
    #         ret = eval(self.script, self._globals)
    #         if _type != None:
    #             return _type(ret)
    #         else:
    #             return ret

def queue_all_observables(analysis_name,is_data): 
	# Observable("sel_plep_z0", do=['hist']).queue() # what is this histogram? -DT
	Observable("all_nmuon", binning=(50, 0, 50), do=['oldAnalysis','ToyAnalysis']).queue()
	Observable("all_nel", binning=(50, 0, 50), do=['oldAnalysis','ToyAnalysis']).queue()
	Observable("all_muon_type", binning=(6, -0.5, 5.5), do=['oldAnalysis','ToyAnalysis']).queue()
	Observable("all_muon_pt", binning=(1000, 0, 1000), do=['oldAnalysis','ToyAnalysis']).queue()
	Observable("all_muon_eta", binning=(40, -10, 10), do=['oldAnalysis','ToyAnalysis']).queue()
	Observable("all_muon_phi", binning=(16, -4, 4), do=['oldAnalysis','ToyAnalysis']).queue()
	Observable("all_muon_quality", binning=(4, -.5, 3.5), do=['oldAnalysis','ToyAnalysis']).queue()
	Observable("all_el_pt", binning=(1000, 0, 1000), do=['oldAnalysis','ToyAnalysis']).queue()
	Observable("all_el_eta", binning=(40, -10, 10), do=['oldAnalysis','ToyAnalysis']).queue()
	Observable("all_el_phi", binning=(16, -4, 4), do=['oldAnalysis','ToyAnalysis']).queue()
	Observable("all_el_quality", binning=(4, -.5, 3.5), do=['oldAnalysis','ToyAnalysis']).queue()
	Observable("all_prompt_muon", binning=(50, 0, 50), do=['oldAnalysis','ToyAnalysis']).queue()
	Observable("all_prompt_electron", binning=(50, 0, 50), do=['oldAnalysis','ToyAnalysis']).queue()
	Observable("all_prompt_lepton", binning=(50, 0, 50), do=['oldAnalysis','ToyAnalysis']).queue()
	Observable("presel_num_trks",binning = (21,-0.5,20.5), do = ['oldAnalysis','ToyAnalysis']).queue()


	reco_histograms("all")
	if is_data: 
		EventType = [""]
	else: 
		EventType = ["","LNC_","LNV_"]
	for i in range(len(EventType)):
		reco_histograms("presel" + EventType[i])
		reco_histograms("DVtype" + EventType[i])
		reco_histograms("cosmic" + EventType[i])
		reco_histograms("mlll" + EventType[i])
		reco_histograms("mDV" + EventType[i])
		reco_histograms("sel" + EventType[i])
		reco_histograms("trkqual" + EventType[i])
		reco_histograms("match" + EventType[i])
		if analysis_name == "KShort": 
			reco_histograms("mass" + EventType[i])
			reco_histograms("alpha" + EventType[i])

		truth_histograms("truth_all")
		truth_histograms("truth_LNC")
		truth_histograms("truth_LNV")
		truth_histograms("truth_LNC_presel")
		truth_histograms("truth_LNC_mass")
		truth_histograms("truth_LNC_charge")
		truth_histograms("truth_LNC_DVtype")
		truth_histograms("truth_LNC_trkqual")
		truth_histograms("truth_LNC_cosmic")
		truth_histograms("truth_LNC_mlll")
		truth_histograms("truth_LNC_mDV")
		truth_histograms("truth_LNC_sel")
		truth_histograms("truth_LNC_match")
		truth_histograms("truth_LNV_presel")
		truth_histograms("truth_LNV_mass")
		truth_histograms("truth_LNV_charge")
		truth_histograms("truth_LNV_DVtype")
		truth_histograms("truth_LNV_trkqual")
		truth_histograms("truth_LNV_cosmic")
		truth_histograms("truth_LNV_mlll")
		truth_histograms("truth_LNV_mDV")
		truth_histograms("truth_LNV_sel")
		truth_histograms("truth_LNV_match")

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
	Observable(selection + "_DV_trk_max_chi2_toSV", binning = (600,0,30), do= ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_trk_min_chi2_toSV", binning = (600,0,30), do= ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_trk_max_d0_wrtSV", binning = (1000,-50,50), do= ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_trk_min_d0_wrtSV", binning = (1000,-50,50), do= ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_trk_max_errd0_wrtSV", binning = (10000,0,50), do= ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_trk_min_errd0_wrtSV", binning = (10000,0,50), do= ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_trk_max_z0_wrtSV", binning = (3000,-30,30), do= ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_trk_min_z0_wrtSV", binning = (3000,-30,30), do= ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_trk_max_errz0_wrtSV", binning = (10000,0,50), do= ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_trk_min_errz0_wrtSV", binning = (10000,0,50), do= ['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_max_dR", binning = (1000,0,10),do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_max_dR_wrtSV", binning = (1000,0,10),do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_maxOpAng", binning = (2000,-1,1),do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_minOpAng", binning = (2000,-1,1),do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_maxd0", binning = (1000,-250,250),do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_mind0", binning = (1000,-250,250),do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_ntrk", binning = (3,-0.5,2.5),do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_ntrk_lrt", binning = (3,-0.5,2.5),do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_ntrk_sel", binning = (3,-0.5,2.5),do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()
	Observable(selection + "_DV_ntrk_assoc", binning = (3,-0.5,2.5),do=['oldAnalysis','ToyAnalysis','KShort'],need_truth = need_truth).queue()

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
