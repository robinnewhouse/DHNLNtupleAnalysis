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
Observable("DV_trk_charge",binning = (12,-5.5,5.5), do = ['hist']).queue()
Observable("DV_trk_chi2",binning = (40,0,20), do = ['hist']).queue()

Observable("DV_x",binning = (2000,-500,500), do = ['hist']).queue()
Observable("DV_y",binning = (2000,-500,500), do = ['hist']).queue()
Observable("DV_z",binning = (2000,-500,500), do = ['hist']).queue()
Observable("DV_r",binning = (500,0,500), do = ['hist']).queue()
Observable("DV_num_trks",binning = (6,-0.5,5.5), do = ['hist']).queue()
Observable("DV_distFromPV",binning = (1000,0,1000), do = ['hist']).queue()
Observable("DV_mass",binning = (100,0,50), do = ['hist']).queue()
Observable("DV_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("DV_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("DV_phi", binning = (16,-4,4), do = ['hist']).queue()
Observable("DV_minOpAng",binning = (100,0,1), do = ['hist']).queue()
Observable("DV_maxOpAng",binning = (100,0,1), do = ['hist']).queue()
Observable("DV_charge",binning = (12,-5.5,5.5), do = ['hist']).queue()
Observable("DV_chi2",binning = (40,0,20), do = ['hist']).queue()
Observable("DV_trk_sep", binning = (30,0,6), do = ['hist']).queue()

# selected DV
Observable("selDV_trk_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("selDV_trk_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("selDV_trk_phi", binning =(16,-4,4), do = ['hist']).queue()
Observable("selDV_trk_d0", do = ['hist']).queue()
Observable("selDV_trk_z0", do = ['hist']).queue()
Observable("selDV_trk_charge",binning = (12,-5.5,5.5), do = ['hist']).queue()
Observable("selDV_trk_chi2",binning = (40,0,20), do = ['hist']).queue()

Observable("selDV_x",binning = (2000,-500,500), do = ['hist']).queue()
Observable("selDV_y",binning = (2000,-500,500), do = ['hist']).queue()
Observable("selDV_z",binning = (2000,-500,500), do = ['hist']).queue()
Observable("selDV_r",binning = (500,0,500), do = ['hist']).queue()
Observable("selDV_num_trks",binning = (6,-0.5,5.5), do = ['hist']).queue()
Observable("selDV_distFromPV",binning = (500,0,500), do = ['hist']).queue()
Observable("selDV_mass",binning = (100,0,50), do = ['hist']).queue()
Observable("selDV_pt",binning = (1000,0,1000), do = ['hist']).queue()
Observable("selDV_eta",binning = (40,-10,10),do = ['hist']).queue()
Observable("selDV_phi", binning = (16,-4,4), do = ['hist']).queue()
Observable("selDV_minOpAng",binning = (100,0,1), do = ['hist']).queue()
Observable("selDV_maxOpAng",binning = (100,0,1), do = ['hist']).queue()
Observable("selDV_charge",binning = (12,-5.5,5.5), do = ['hist']).queue()
Observable("selDV_chi2",binning = (40,0,20), do = ['hist']).queue()
Observable("selDV_trk_sep", binning = (30,0,6), do = ['hist']).queue()




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


Observable("muon_type", binning = (6 ,-0.5,5.5), do = ['hist']).queue()
Observable("muon_pt", binning = (1000 ,0,1000), do = ['hist']).queue()
Observable("muon_eta", binning = (40 ,-10,10), do = ['hist']).queue()
Observable("muon_phi", binning = (16 ,-4,4), do = ['hist']).queue()
Observable("el_pt", binning = (1000 ,0,1000), do = ['hist']).queue()
Observable("el_eta", binning = (40 ,-10,10), do = ['hist']).queue()
Observable("el_phi", binning = (16 ,-4,4), do = ['hist']).queue()




# Observable("truth_DV_trk_pt",binning = (1000,0,1000), do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_trk_eta",binning = (40,-10,10),do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_trk_phi", binning = (16,-3,3), do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_trk_d0", do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_trk_z0", do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_trk_charge",binning = (12,-5.5,5.5), do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_num_trks",binning = (6,-0.5,5.5), do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_trk_chi2_toSV",binning = (40,0,20), do = ['hist'], need_truth = True).queue()
Observable("truth_DV_x",binning = (2000,-500,500), do = ['hist'], need_truth = True).queue()
Observable("truth_DV_y",binning = (2000,-500,500), do = ['hist'], need_truth = True).queue()
Observable("truth_DV_z",binning = (2000,-500,500), do = ['hist'], need_truth = True).queue()
Observable("truth_DV_r",binning = (500,0,500), do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_distFromPV",binning = (500,0,500), do = ['hist']).queue()
Observable("truth_DV_mass",binning = (100,0,50), do = ['hist'], need_truth = True).queue()
Observable("truth_DV_pt",binning = (1000,0,1000), do = ['hist'], need_truth = True).queue()
Observable("truth_DV_eta",binning = (40,-10,10),do = ['hist'], need_truth = True).queue()
Observable("truth_DV_phi", binning = (16,-4,4), do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_minOpAng",binning = (100,0,1), do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_maxOpAng",binning = (100,0,1), do = ['hist'], need_truth = True).queue()
# Observable("truth_DV_charge",binning = (12,-5.5,5.5), do = ['hist']).queue()



