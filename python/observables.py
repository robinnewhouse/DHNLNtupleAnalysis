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

# all DV
Observable("DV_trk_pt",binning = (5000,0,1000), do = ['hist']).queue()
Observable("DV_trk_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("DV_trk_phi", binning =(64,-4,4), do = ['hist']).queue()
Observable("DV_trk_d0", do = ['hist']).queue()
Observable("DV_trk_z0", do = ['hist']).queue()
Observable("DV_trk_charge",binning = (11,-5.5,5.5), do = ['hist']).queue()
Observable("DV_trk_chi2",binning = (40,0,20), do = ['hist']).queue()

Observable("DV_x",binning = (2000,-500,500), do = ['hist']).queue()
Observable("DV_y",binning = (2000,-500,500), do = ['hist']).queue()
Observable("DV_z",binning = (2000,-500,500), do = ['hist']).queue()
Observable("DV_r",binning = (500,0,500), do = ['hist']).queue()
Observable("DV_num_trks",binning = (6,-0.5,5.5), do = ['hist']).queue()
Observable("DV_distFromPV",binning = (1000,0,1000), do = ['hist']).queue()
Observable("DV_mass",binning = (10000,0,5000), do = ['hist']).queue()
Observable("DV_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("DV_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("DV_phi", binning = (16,-4,4), do = ['hist']).queue()
Observable("DV_minOpAng",binning = (100,0,1), do = ['hist']).queue()
Observable("DV_maxOpAng",binning = (100,0,1), do = ['hist']).queue()
Observable("DV_charge",binning = (11,-5.5,5.5), do = ['hist']).queue()
Observable("DV_chi2",binning = (40,0,20), do = ['hist']).queue()
# Observable("DV_trk_sep", binning = (30,0,6), do = ['hist']).queue()

# selected DV
Observable("selDV_trk_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("selDV_trk_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("selDV_trk_phi", binning =(16,-4,4), do = ['hist']).queue()
Observable("selDV_trk_d0", do = ['hist']).queue()
Observable("selDV_trk_z0", do = ['hist']).queue()
Observable("selDV_trk_charge",binning = (11,-5.5,5.5), do = ['hist']).queue()
Observable("selDV_trk_chi2",binning = (40,0,20), do = ['hist']).queue()

Observable("selDV_x",binning = (2000,-500,500), do = ['hist']).queue()
Observable("selDV_y",binning = (2000,-500,500), do = ['hist']).queue()
Observable("selDV_z",binning = (2000,-500,500), do = ['hist']).queue()
Observable("selDV_r",binning = (500,0,500), do = ['hist']).queue()
Observable("selDV_num_trks",binning = (6,-0.5,5.5), do = ['hist']).queue()
Observable("selDV_distFromPV",binning = (500,0,500), do = ['hist']).queue()
Observable("selDV_mass",binning = (10000,0,5000), do = ['hist']).queue()
Observable("selDV_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("selDV_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("selDV_phi", binning = (16,-4,4), do = ['hist']).queue()
Observable("selDV_minOpAng",binning = (100,0,1), do = ['hist']).queue()
Observable("selDV_maxOpAng",binning = (100,0,1), do = ['hist']).queue()
Observable("selDV_charge",binning = (11,-5.5,5.5), do = ['hist']).queue()
Observable("selDV_chi2",binning = (40,0,20), do = ['hist']).queue()
# Observable("selDV_trk_sep", binning = (30,0,6), do = ['hist']).queue()




Observable("plep_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("plep_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("plep_phi", binning = (16,-4,4), do = ['hist']).queue()
Observable("plep_d0", do = ['hist']).queue()
Observable("plep_z0", do = ['hist']).queue()

Observable("selplep_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("selplep_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("selplep_phi", binning = (16,-4,4), do = ['hist']).queue()
Observable("selplep_d0", do = ['hist']).queue()
Observable("selplep_z0", do = ['hist']).queue()

Observable("nmuon", binning=(50, 0, 50), do=['hist']).queue()
Observable("nel", binning=(50, 0, 50), do=['hist']).queue()
Observable("muon_type", binning = (6 ,-0.5,5.5), do = ['hist']).queue()
Observable("muon_pt", binning = (1000 ,0,1000), do = ['hist']).queue()
Observable("muon_eta", binning = (40 ,-10,10), do = ['hist']).queue()
Observable("muon_phi", binning = (16 ,-4,4), do = ['hist']).queue()
Observable("el_pt", binning = (1000 ,0,1000), do = ['hist']).queue()
Observable("el_eta", binning = (40 ,-10,10), do = ['hist']).queue()
Observable("el_phi", binning = (16 ,-4,4), do = ['hist']).queue()


# SS
Observable("SSplep_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("SSplep_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("SSplep_phi", binning = (16,-4,4), do = ['hist']).queue()
Observable("SSplep_d0", do = ['hist']).queue()
Observable("SSplep_z0", do = ['hist']).queue()

Observable("SSDV_trk_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("SSDV_trk_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("SSDV_trk_phi", binning =(16,-4,4), do = ['hist']).queue()
Observable("SSDV_trk_d0", do = ['hist']).queue()
Observable("SSDV_trk_z0", do = ['hist']).queue()
Observable("SSDV_trk_charge",binning = (11,-5.5,5.5), do = ['hist']).queue()
Observable("SSDV_trk_chi2",binning = (40,0,20), do = ['hist']).queue()

Observable("SSDV_x",binning = (2000,-500,500), do = ['hist']).queue()
Observable("SSDV_y",binning = (2000,-500,500), do = ['hist']).queue()
Observable("SSDV_z",binning = (2000,-500,500), do = ['hist']).queue()
Observable("SSDV_r",binning = (500,0,500), do = ['hist']).queue()
Observable("SSDV_num_trks",binning = (6,-0.5,5.5), do = ['hist']).queue()
Observable("SSDV_distFromPV",binning = (500,0,500), do = ['hist']).queue()
Observable("SSDV_mass",binning = (10000,0,5000), do = ['hist']).queue()
Observable("SSDV_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("SSDV_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("SSDV_phi", binning = (16,-4,4), do = ['hist']).queue()
Observable("SSDV_minOpAng",binning = (100,0,1), do = ['hist']).queue()
Observable("SSDV_maxOpAng",binning = (100,0,1), do = ['hist']).queue()
Observable("SSDV_charge",binning = (11,-5.5,5.5), do = ['hist']).queue()
Observable("SSDV_chi2",binning = (40,0,20), do = ['hist']).queue()

#DV type
Observable("DVtypeplep_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("DVtypeplep_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("DVtypeplep_phi", binning = (16,-4,4), do = ['hist']).queue()
Observable("DVtypeplep_d0", do = ['hist']).queue()
Observable("DVtypeplep_z0", do = ['hist']).queue()

Observable("DVtypeDV_trk_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("DVtypeDV_trk_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("DVtypeDV_trk_phi", binning =(16,-4,4), do = ['hist']).queue()
Observable("DVtypeDV_trk_d0", do = ['hist']).queue()
Observable("DVtypeDV_trk_z0", do = ['hist']).queue()
Observable("DVtypeDV_trk_charge",binning = (11,-5.5,5.5), do = ['hist']).queue()
Observable("DVtypeDV_trk_chi2",binning = (40,0,20), do = ['hist']).queue()

Observable("DVtypeDV_x",binning = (2000,-500,500), do = ['hist']).queue()
Observable("DVtypeDV_y",binning = (2000,-500,500), do = ['hist']).queue()
Observable("DVtypeDV_z",binning = (2000,-500,500), do = ['hist']).queue()
Observable("DVtypeDV_r",binning = (500,0,500), do = ['hist']).queue()
Observable("DVtypeDV_num_trks",binning = (6,-0.5,5.5), do = ['hist']).queue()
Observable("DVtypeDV_distFromPV",binning = (500,0,500), do = ['hist']).queue()
Observable("DVtypeDV_mass",binning = (10000,0,5000), do = ['hist']).queue()
Observable("DVtypeDV_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("DVtypeDV_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("DVtypeDV_phi", binning = (16,-4,4), do = ['hist']).queue()
Observable("DVtypeDV_minOpAng",binning = (100,0,1), do = ['hist']).queue()
Observable("DVtypeDV_maxOpAng",binning = (100,0,1), do = ['hist']).queue()
Observable("DVtypeDV_charge",binning = (11,-5.5,5.5), do = ['hist']).queue()
Observable("DVtypeDV_chi2",binning = (40,0,20), do = ['hist']).queue()

## trkqual

Observable("trkqualplep_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("trkqualplep_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("trkqualplep_phi", binning = (16,-4,4), do = ['hist']).queue()
Observable("trkqualplep_d0", do = ['hist']).queue()
Observable("trkqualplep_z0", do = ['hist']).queue()

Observable("trkqualDV_trk_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("trkqualDV_trk_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("trkqualDV_trk_phi", binning =(16,-4,4), do = ['hist']).queue()
Observable("trkqualDV_trk_d0", do = ['hist']).queue()
Observable("trkqualDV_trk_z0", do = ['hist']).queue()
Observable("trkqualDV_trk_charge",binning = (11,-5.5,5.5), do = ['hist']).queue()
Observable("trkqualDV_trk_chi2",binning = (40,0,20), do = ['hist']).queue()

Observable("trkqualDV_x",binning = (2000,-500,500), do = ['hist']).queue()
Observable("trkqualDV_y",binning = (2000,-500,500), do = ['hist']).queue()
Observable("trkqualDV_z",binning = (2000,-500,500), do = ['hist']).queue()
Observable("trkqualDV_r",binning = (500,0,500), do = ['hist']).queue()
Observable("trkqualDV_num_trks",binning = (6,-0.5,5.5), do = ['hist']).queue()
Observable("trkqualDV_distFromPV",binning = (500,0,500), do = ['hist']).queue()
Observable("trkqualDV_mass",binning = (10000,0,5000), do = ['hist']).queue()
Observable("trkqualDV_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("trkqualDV_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("trkqualDV_phi", binning = (16,-4,4), do = ['hist']).queue()
Observable("trkqualDV_minOpAng",binning = (100,0,1), do = ['hist']).queue()
Observable("trkqualDV_maxOpAng",binning = (100,0,1), do = ['hist']).queue()
Observable("trkqualDV_charge",binning = (11,-5.5,5.5), do = ['hist']).queue()
Observable("trkqualDV_chi2",binning = (40,0,20), do = ['hist']).queue()

# cosmic

Observable("cosmicplep_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("cosmicplep_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("cosmicplep_phi", binning = (16,-4,4), do = ['hist']).queue()
Observable("cosmicplep_d0", do = ['hist']).queue()
Observable("cosmicplep_z0", do = ['hist']).queue()

Observable("cosmicDV_trk_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("cosmicDV_trk_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("cosmicDV_trk_phi", binning =(16,-4,4), do = ['hist']).queue()
Observable("cosmicDV_trk_d0", do = ['hist']).queue()
Observable("cosmicDV_trk_z0", do = ['hist']).queue()
Observable("cosmicDV_trk_charge",binning = (11,-5.5,5.5), do = ['hist']).queue()
Observable("cosmicDV_trk_chi2",binning = (40,0,20), do = ['hist']).queue()

Observable("cosmicDV_x",binning = (2000,-500,500), do = ['hist']).queue()
Observable("cosmicDV_y",binning = (2000,-500,500), do = ['hist']).queue()
Observable("cosmicDV_z",binning = (2000,-500,500), do = ['hist']).queue()
Observable("cosmicDV_r",binning = (500,0,500), do = ['hist']).queue()
Observable("cosmicDV_num_trks",binning = (6,-0.5,5.5), do = ['hist']).queue()
Observable("cosmicDV_distFromPV",binning = (500,0,500), do = ['hist']).queue()
Observable("cosmicDV_mass",binning = (10000,0,5000), do = ['hist']).queue()
Observable("cosmicDV_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("cosmicDV_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("cosmicDV_phi", binning = (16,-4,4), do = ['hist']).queue()
Observable("cosmicDV_minOpAng",binning = (100,0,1), do = ['hist']).queue()
Observable("cosmicDV_maxOpAng",binning = (100,0,1), do = ['hist']).queue()
Observable("cosmicDV_charge",binning = (11,-5.5,5.5), do = ['hist']).queue()
Observable("cosmicDV_chi2",binning = (40,0,20), do = ['hist']).queue()

# mlll

Observable("mlllplep_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("mlllplep_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("mlllplep_phi", binning = (16,-4,4), do = ['hist']).queue()
Observable("mlllplep_d0", do = ['hist']).queue()
Observable("mlllplep_z0", do = ['hist']).queue()

Observable("mlllDV_trk_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("mlllDV_trk_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("mlllDV_trk_phi", binning =(16,-4,4), do = ['hist']).queue()
Observable("mlllDV_trk_d0", do = ['hist']).queue()
Observable("mlllDV_trk_z0", do = ['hist']).queue()
Observable("mlllDV_trk_charge",binning = (11,-5.5,5.5), do = ['hist']).queue()
Observable("mlllDV_trk_chi2",binning = (40,0,20), do = ['hist']).queue()

Observable("mlllDV_x",binning = (2000,-500,500), do = ['hist']).queue()
Observable("mlllDV_y",binning = (2000,-500,500), do = ['hist']).queue()
Observable("mlllDV_z",binning = (2000,-500,500), do = ['hist']).queue()
Observable("mlllDV_r",binning = (500,0,500), do = ['hist']).queue()
Observable("mlllDV_num_trks",binning = (6,-0.5,5.5), do = ['hist']).queue()
Observable("mlllDV_distFromPV",binning = (500,0,500), do = ['hist']).queue()
Observable("mlllDV_mass",binning = (10000,0,5000), do = ['hist']).queue()
Observable("mlllDV_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("mlllDV_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("mlllDV_phi", binning = (16,-4,4), do = ['hist']).queue()
Observable("mlllDV_minOpAng",binning = (100,0,1), do = ['hist']).queue()
Observable("mlllDV_maxOpAng",binning = (100,0,1), do = ['hist']).queue()
Observable("mlllDV_charge",binning = (11,-5.5,5.5), do = ['hist']).queue()
Observable("mlllDV_chi2",binning = (40,0,20), do = ['hist']).queue()




# Bug with the truth variables in DHNL alg. Need to fix -DT
# Observable("truth_DV_x",binning = (2000,-500,500), do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_y",binning = (2000,-500,500), do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_z",binning = (2000,-500,500), do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_r",binning = (500,0,500), do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_mass",binning = (100,0,50), do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_pt",binning = (1000,0,1000), do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_eta",binning = (40,-10,10),do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_phi", binning = (16,-4,4), do = ['hist'], need_truth = True).queue()





# Observable("truth_DV_trk_pt",binning = (1000,0,1000), do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_trk_eta",binning = (40,-10,10),do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_trk_phi", binning = (16,-3,3), do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_trk_d0", do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_trk_z0", do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_trk_charge",binning = (12,-5.5,5.5), do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_num_trks",binning = (6,-0.5,5.5), do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_trk_chi2_toSV",binning = (40,0,20), do = ['hist'], need_truth = True).queue()

# Observable("truth_DV_distFromPV",binning = (500,0,500), do = ['hist']).queue()

# Observable("truth_DV_minOpAng",binning = (100,0,1), do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_maxOpAng",binning = (100,0,1), do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_charge",binning = (12,-5.5,5.5), do = ['hist']).queue()



