import math
import ROOT
ObservableList = []



class Observable(object):
	def __init__(self, name, binning = (20, -1000, 1000), script = None, style = 'single', do = ['hist'], only = None, title = '{self.name}', dtype = float, default = -1, need_truth = False):