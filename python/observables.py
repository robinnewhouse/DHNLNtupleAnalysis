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


Observable("sel_plep_z0", do=['hist']).queue()

Observable("all_nmuon", binning=(50, 0, 50), do=['hist']).queue()
Observable("all_nel", binning=(50, 0, 50), do=['hist']).queue()
Observable("all_muon_type", binning=(6, -0.5, 5.5), do=['hist']).queue()
Observable("all_muon_pt", binning=(1000, 0, 1000), do=['hist']).queue()
Observable("all_muon_eta", binning=(40, -10, 10), do=['hist']).queue()
Observable("all_muon_phi", binning=(16, -4, 4), do=['hist']).queue()
Observable("all_muon_quality", binning=(4, -.5, 3.5), do=['hist']).queue()
Observable("all_el_pt", binning=(1000, 0, 1000), do=['hist']).queue()
Observable("all_el_eta", binning=(40, -10, 10), do=['hist']).queue()
Observable("all_el_phi", binning=(16, -4, 4), do=['hist']).queue()
Observable("all_el_quality", binning=(4, -.5, 3.5), do=['hist']).queue()

Observable("all_prompt_muon", binning=(50, 0, 50), do=['hist']).queue()
Observable("all_prompt_electron", binning=(50, 0, 50), do=['hist']).queue()
Observable("all_prompt_lepton", binning=(50, 0, 50), do=['hist']).queue()


def reco_histograms(selection):
	# histograms after selection selection
	Observable(selection + "_nonlep_pt", binning=(1000, 0, 1000), do=['hist']).queue()
	Observable(selection + "_nonlep_eta", binning=(40, -10, 10), do=['hist']).queue()
	Observable(selection + "_nonlep_phi", binning=(16, -4, 4), do=['hist']).queue()
	Observable(selection + "_nonlep_d0", binning=(2000, -1000, 1000), do=['hist']).queue()
	Observable(selection + "_plep_pt", binning=(1000, 0, 1000), do=['hist']).queue()
	Observable(selection + "_plep_eta", binning=(40, -10, 10), do=['hist']).queue()
	Observable(selection + "_plep_phi", binning=(16, -4, 4), do=['hist']).queue()
	Observable(selection + "_plep_d0", binning=(2000, -1000, 1000), do=['hist']).queue()
	Observable(selection + "_plep_z0", binning=(2000, -1000, 1000), do=['hist']).queue()

	Observable(selection + "_DV_trk_pt",binning = (1000,0,1000), do = ['hist']).queue()
	Observable(selection + "_DV_trk_eta",binning = (160,-10,10),do = ['hist']).queue()
	Observable(selection + "_DV_trk_phi", binning =(64,-4,4), do = ['hist']).queue()
	Observable(selection + "_DV_trk_d0", do = ['hist']).queue()
	Observable(selection + "_DV_trk_z0", do = ['hist']).queue()
	Observable(selection + "_DV_trk_charge",binning = (11,-5.5,5.5), do = ['hist']).queue()
	Observable(selection + "_DV_trk_chi2",binning = (40,0,20), do = ['hist']).queue()
	Observable(selection + "_DV_trk_dpt",binning = (2000,0,1000), do = ['hist']).queue()
	Observable(selection + "_DV_trk_deta",binning = (160,0,20),do = ['hist']).queue()
	Observable(selection + "_DV_trk_dphi", binning = (64,0,8), do = ['hist']).queue()
	Observable(selection + "_DV_trk_dR", binning = (1010,-5,10), do = ['hist']).queue()

	Observable(selection + "_DV_x",binning = (2000,-500,500), do = ['hist']).queue()
	Observable(selection + "_DV_y",binning = (2000,-500,500), do = ['hist']).queue()
	Observable(selection + "_DV_z",binning = (2000,-500,500), do = ['hist']).queue()
	Observable(selection + "_DV_r",binning = (500,0,500), do = ['hist']).queue()
	Observable(selection + "_DV_num_trks",binning = (6,-0.5,5.5), do = ['hist']).queue()
	Observable(selection + "_DV_distFromPV",binning = (500,0,500), do = ['hist']).queue()
	Observable(selection + "_DV_mass", binning=(2000, 0, 1000), do=['hist']).queue()
	Observable(selection + "_DV_pt",binning = (1000,0,1000), do = ['hist']).queue()
	Observable(selection + "_DV_sum_track_pt",binning = (1000,0,1000), do = ['hist']).queue()
	Observable(selection + "_DV_sum_track_pt_wrt_pv",binning = (1000,0,1000), do = ['hist']).queue()
	Observable(selection + "_DV_sum_track_pt_diff",binning = (10000,0,100), do = ['hist']).queue()
	Observable(selection + "_DV_sum_track_charge",binning = (4,-1.5,1.5), do = ['hist']).queue()
	Observable(selection + "_DV_eta",binning = (160,-10,10),do = ['hist']).queue()
	Observable(selection + "_DV_phi", binning = (64,-4,4), do = ['hist']).queue()
	Observable(selection + "_DV_minOpAng",binning = (100,0,1), do = ['hist']).queue()
	Observable(selection + "_DV_maxOpAng",binning = (100,0,1), do = ['hist']).queue()
	Observable(selection + "_DV_charge",binning = (11,-5.5,5.5), do = ['hist']).queue()
	Observable(selection + "_DV_chi2",binning = (40,0,20), do = ['hist']).queue()
	Observable(selection + "_DV_alpha",binning = (400,-4,4), do = ['hist']).queue()
	Observable(selection + "_mvis",binning = (10000,0,5000), do = ['hist']).queue()
	Observable(selection + "_HNLm",binning = (10005,-5,5000), do = ['hist']).queue()
	Observable(selection + "_HNLm_altbinning",binning = (15,0,30), do = ['hist']).queue()
	Observable(selection + "_alt_HNLm",binning = (10005,-5,5000), do = ['hist']).queue()
	Observable(selection + "_HNLm_fixWmass",binning = (10005,-5,5000), do = ['hist']).queue()
	Observable(selection + "_HNLpt",binning = (1000,0,1000), do = ['hist']).queue()
	Observable(selection + "_HNLeta",binning = (160,-10,10),do = ['hist']).queue()
	Observable(selection + "_HNLphi", binning = (64,-4,4), do = ['hist']).queue()
	Observable(selection + "_mtrans",binning = (10000,0,5000), do = ['hist']).queue()
	Observable(selection + "_DV_redmass",binning = (10005,-5,5000), do = ['hist']).queue()
	Observable(selection + "_DV_redmassvis",binning = (10005,-5,5000), do = ['hist']).queue()
	Observable(selection + "_DV_redmassHNL",binning = (10005,-5,5000), do = ['hist']).queue()

reco_histograms("all")
reco_histograms("presel")
reco_histograms("alpha")
reco_histograms("mass")
reco_histograms("charge")
reco_histograms("DVtype")
reco_histograms("trkqual")
reco_histograms("cosmic")
reco_histograms("mlll")
reco_histograms("mDV")
reco_histograms("HNLpt")
reco_histograms("tight1")
reco_histograms("tight2")
reco_histograms("sel")
reco_histograms("match")


def truth_histograms(selection): 
	# histograms after selection selection
	Observable("truth_" + selection + "_W_mass",binning = (10000,0,5000), do = ['hist'], need_truth = True).queue()
	Observable("truth_" + selection + "_W_pt",binning = (1000,0,1000), do = ['hist'], need_truth = True).queue()
	Observable("truth_" + selection + "_W_eta",binning = (40,-10,10),do = ['hist'], need_truth = True).queue()
	Observable("truth_" + selection + "_W_phi", binning = (16,-4,4), do = ['hist'], need_truth = True).queue()
	Observable("truth_" + selection + "_HNL_mass",binning = (10000,0,5000), do = ['hist'], need_truth = True).queue()
	Observable("truth_" + selection + "_mHNLcalc",binning = (10000,0,5000), do = ['hist'], need_truth = True).queue()
	Observable("truth_" + selection + "_HNL_pt",binning = (1000,0,1000), do = ['hist'], need_truth = True).queue()
	Observable("truth_" + selection + "_HNL_eta",binning = (40,-10,10),do = ['hist'], need_truth = True).queue()
	Observable("truth_" + selection + "_HNL_phi", binning = (16,-4,4), do = ['hist'], need_truth = True).queue()

	Observable("truth_" + selection + "_plep_pt",binning = (1000,0,1000), do = ['hist'], need_truth = True).queue()
	Observable("truth_" + selection + "_plep_eta",binning = (40,-10,10),do = ['hist'], need_truth = True).queue()
	Observable("truth_" + selection + "_plep_phi", binning = (16,-4,4), do = ['hist'], need_truth = True).queue()
	Observable("truth_" + selection + "_plep_mass",binning = (10000,0,5000), do = ['hist'], need_truth = True).queue()

	# Observable("truth_" + selection + "_plep_d0", do = ['hist'], need_truth = True).queue() # no d0 or z0 for truth particles
	# Observable("truth_" + selection + "_plep_z0", do = ['hist'], need_truth = True).queue()

	Observable("truth_" + selection + "_DV_trk_pt",binning = (1000,0,1000), do = ['hist'], need_truth = True).queue()
	Observable("truth_" + selection + "_DV_trk_eta",binning = (160,-10,10),do = ['hist'], need_truth = True).queue()
	Observable("truth_" + selection + "_DV_trk_phi", binning =(64,-4,4), do = ['hist'], need_truth = True).queue()
	Observable("truth_" + selection + "_lep1_trk_pt", binning=(1000, 0, 1000), do=['hist'], need_truth=True).queue()
	Observable("truth_" + selection + "_lep1_trk_eta", binning=(160, -10, 10), do=['hist'], need_truth=True).queue()
	Observable("truth_" + selection + "_lep1_trk_phi", binning=(64, -4, 4), do=['hist'], need_truth=True).queue()
	Observable("truth_" + selection + "_lep2_trk_pt", binning=(1000, 0, 1000), do=['hist'], need_truth=True).queue()
	Observable("truth_" + selection + "_lep2_trk_eta", binning=(160, -10, 10), do=['hist'], need_truth=True).queue()
	Observable("truth_" + selection + "_lep2_trk_phi", binning=(64, -4, 4), do=['hist'], need_truth=True).queue()
	# Observable("truth_" + selection + "_DV_trk_d0", do = ['hist'], need_truth = True).queue() # no d0 or z0 for truth particles
	# Observable("truth_" + selection + "_DV_trk_z0", do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_DV_trk_charge",binning = (11,-5.5,5.5), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_DV_trk_chi2",binning = (40,0,20), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_DV_trk_dpt",binning = (2000,0,1000), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_DV_trk_deta",binning = (160,0,20),do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_DV_trk_dphi", binning = (64,0,8), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_DV_trk_dR", binning = (1010,-5,10), do = ['hist'], need_truth = True).queue()

	Observable("truth_" + selection + "_DV_x",binning = (2000,-500,500), do = ['hist'], need_truth = True).queue()
	Observable("truth_" + selection + "_DV_y",binning = (2000,-500,500), do = ['hist'], need_truth = True).queue()
	Observable("truth_" + selection + "_DV_z",binning = (2000,-500,500), do = ['hist'], need_truth = True).queue()
	Observable("truth_" + selection + "_DV_r",binning = (500,0,500), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_DV_num_trks",binning = (6,-0.5,5.5), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_DV_distFromPV",binning = (500,0,500), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_DV_mass",binning = (10000,0,5000), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_DV_pt",binning = (1000,0,1000), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_DV_eta",binning = (160,-10,10),do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_DV_phi", binning = (64,-4,4), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_DV_minOpAng",binning = (100,0,1), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_DV_maxOpAng",binning = (100,0,1), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_DV_charge",binning = (11,-5.5,5.5), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_DV_chi2",binning = (40,0,20), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_mvis",binning = (10000,0,5000), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_HNLm",binning = (10005,-5,5000), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_HNLm2",binning = (10005,-5,5000), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_HNLpt",binning = (1000,0,1000), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_HNLeta",binning = (160,-10,10),do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_HNLphi", binning = (64,-4,4), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_mtrans",binning = (10000,0,5000), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_mtrans_rot",binning = (10000,0,5000), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_DV_redmass",binning = (10005,-5,5000), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_DV_redmassvis",binning = (10005,-5,5000), do = ['hist'], need_truth = True).queue()
	# Observable("truth_" + selection + "_DV_redmassHNL",binning = (10005,-5,5000), do = ['hist'], need_truth = True).queue()

truth_histograms("all")
truth_histograms("presel")
truth_histograms("alpha")
truth_histograms("mass")
truth_histograms("charge")
truth_histograms("DVtype")
truth_histograms("trkqual")
truth_histograms("cosmic")
truth_histograms("mlll")
truth_histograms("mDV")
truth_histograms("HNLpt")
truth_histograms("tight1")
truth_histograms("tight2")
truth_histograms("sel")
truth_histograms("match")


