import treenames	
import selections 
import helpers
import ROOT, time
from ROOT import *
gROOT.LoadMacro("AtlasStyle.C")
gROOT.LoadMacro("AtlasUtils.C")
gROOT.LoadMacro("AtlasLabels.C")



if __name__ == '__main__':
	
	# emu file: 
	# file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/fixElmatching_mc16_13TeV.311660.Pythia8EvtGen_A14NNPDF23LO_WmuHNL50_20G_lt10dd_el.merge.DAOD_RPVLL.e7422_e5984_a875_r10739_r10706.root"
	# file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_20G_lt10dd_emu_wTightmuel.root"
	# mumu file:
	# file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_10G_lt10dd_mumu.root"
	file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_10G_lt10dd_mumu_wTightmuel.root"

	treename = "outTree"
	tree = treenames.Tree(file, treename)

	entries = len(tree.dvmass)
	# entries = 300

	hcutflow= ROOT.TH1D("hcutflow","cutflow",13,-0.5,12.5)
  	hcutflow.GetXaxis().SetBinLabel(1, "all")
  	hcutflow.GetXaxis().SetBinLabel(2, "trigger")
  	hcutflow.GetXaxis().SetBinLabel(3, "filter")
  	hcutflow.GetXaxis().SetBinLabel(4, "tight pmu")
  	hcutflow.GetXaxis().SetBinLabel(5, "DV")  
  	hcutflow.GetXaxis().SetBinLabel(6, "fiducial")  
  	hcutflow.GetXaxis().SetBinLabel(7, "2-track DV")  
  	hcutflow.GetXaxis().SetBinLabel(8, "OS DV")  
  	# if DV_type == '1':
  		# hcutflow.GetXaxis().SetBinLabel(9, "2-muon DV")  
  	# if DV_type == '0':
  	hcutflow.GetXaxis().SetBinLabel(9, "2-muon DV")  
  	# hcutflow.GetXaxis().SetBinLabel(9, "mu-el DV")  
  	hcutflow.GetXaxis().SetBinLabel(10, "2-tight-lepton DV")  
  	hcutflow.GetXaxis().SetBinLabel(11, "cosmic veto")    
  	hcutflow.GetXaxis().SetBinLabel(12, "m_{lll}")      
  	hcutflow.GetXaxis().SetBinLabel(13, "mDV") 

  	nmuons =0 

  	print "looping over ", entries, " events!"
  	start = time.clock()
	for ievt in xrange(entries):


		passTrigger = False
		passHNLfilter = False
		passPmuon = False
		passDV = False
		passFid = False
		passDVntracks = False
		passOSDV = False
		passDVtype = False
		passTrackqual = False
		passCosmicveto = False
		passMlll = False
		passDVmasscut = False



		evt = helpers.Event(tree=tree, ievt = ievt , idv = None)
		ndv = len(tree.dvx[ievt])





		hcutflow.Fill(0)

		#------------------------------------------------------------------------------------------
		# Check if event passes trigger
		#------------------------------------------------------------------------------------------
		Trigger = selections.Trigger(evt = evt, plepton="muon",trigger = "HLT_mu26_ivarmedium")
		Triggercut = Trigger.passes()
		if Triggercut: 
			hcutflow.Fill(1)
			passesTrigger = Triggercut
		else:
			continue 	
		#------------------------------------------------------------------------------------------


		#------------------------------------------------------------------------------------------
		# Check if event passes HNL filter 
		#------------------------------------------------------------------------------------------
		# HNLfilter = selections.Filter(evt= evt, _filter="mu-mu")
		# HNLfiltercut = HNLfilter.passes()
		# if HNLfiltercut: 
		# 	# hcutflow.Fill(2)
		# 	passHNLfilter = HNLfiltercut
		# else:
		# 	continue 	
		#------------------------------------------------------------------------------------------


		#------------------------------------------------------------------------------------------
		# Find prompt muon for event 
		#------------------------------------------------------------------------------------------
		Pmuon = selections.Plepton(evt = evt, lepton="muon")
		Pmuoncut = Pmuon.passes()

		if Pmuoncut: 
			hcutflow.Fill(3)
			passPmuon = Pmuoncut
		else:
			continue 	
		#------------------------------------------------------------------------------------------

		#------------------------------------------------------------------------------------------
		# Count number of DV in the event
		#------------------------------------------------------------------------------------------
		nDV = selections.nDV(evt=evt)
		nDVcut = nDV.passes()
		if nDVcut: 
			hcutflow.Fill(4)
			passDV = nDVcut	
		else:
			continue 
		#------------------------------------------------------------------------------------------

		#------------------------------------------------------------------------------------------
		# Pre-selection 
		preSel = Triggercut and Pmuoncut and nDVcut
		# preSel = passTrigger and passHNLfilter and  passPmuon and passDV
		#------------------------------------------------------------------------------------------

	

		if preSel: #only continue if pre-selection is True
			for idv in xrange(ndv):
				DVevt = helpers.Event(tree=tree, ievt = ievt , idv = idv)

				muons = helpers.Tracks()
				muons.getMuons(evt= DVevt)
				if muons.lepIndex >= 0:
					nmuons = nmuons + 1


				#------------------------------------------------------------------------------------------
				# Calcuate the radius of DV  
				#------------------------------------------------------------------------------------------
				DVradius = selections.DVradius(evt= DVevt)
				Fidvolcut = DVradius.passes(_min=4,_max=300)
				if Fidvolcut:
					if passFid == False:  #only fill cut flow once per DV!!
						hcutflow.Fill(5)
					passFid = Fidvolcut		
				else: 
					continue
				#------------------------------------------------------------------------------------------
				

				#------------------------------------------------------------------------------------------
				# Calculate number of tracks in the DV
				#------------------------------------------------------------------------------------------
				DVntracks = selections.DVntracks(evt= DVevt,ntrk=2)
				DVntrackscut = DVntracks.passes()
				if DVntrackscut:
					if passDVntracks == False:  #only fill cut flow once per DV!!
						hcutflow.Fill(6)
					passDVntracks = DVntrackscut
				else:
					continue
				#------------------------------------------------------------------------------------------

				
				#------------------------------------------------------------------------------------------
				# Check charge of the tracks in DV 
				#------------------------------------------------------------------------------------------
				OSDV = selections.OSDV(evt= DVevt)
				OSDVcut = OSDV.passes()
				if OSDVcut:
					if passOSDV == False:  #only fill cut flow once per DV!!
						hcutflow.Fill(7)		
					passOSDV = OSDVcut	
				else:
					continue
				#------------------------------------------------------------------------------------------		
				# Count if DV has the correct number of leptons 
				#------------------------------------------------------------------------------------------
				DVtype = selections.DVtype(evt= DVevt,decayprod="mumu")
				DVtypecut = DVtype.passes()
				if DVtypecut:
					if passDVtype == False:  #only fill cut flow once per DV!!
						hcutflow.Fill(8)
					passDVtype = DVtypecut
				else:
					continue
				#------------------------------------------------------------------------------------------
				# count the number of lepton tracks with a given quality 
				#------------------------------------------------------------------------------------------
				Trackqual = selections.Trackqual(evt=DVevt, quality="2-tight")
				Trackqualcut = Trackqual.passes()
				if Trackqualcut: 
					if passTrackqual == False:
						hcutflow.Fill(9)
					passTrackqual = Trackqualcut
				else: 
					continue 




				#------------------------------------------------------------------------------------------

		
				if Fidvolcut and DVntrackscut and OSDVcut and DVtypecut and Trackqual: # only do cosmicveto & mass cuts if you have a good 2 lep DV
					#------------------------------------------------------------------------------------------
					# Get prompt muon
					#------------------------------------------------------------------------------------------
					pMuonvec = Pmuon.plepVec
					#------------------------------------------------------------------------------------------

					#------------------------------------------------------------------------------------------
					# Get DV muons 
					#------------------------------------------------------------------------------------------
					muons = helpers.Tracks()
					muons.getMuons(evt= DVevt)
					muVec = muons.lepVec
					#------------------------------------------------------------------------------------------

					#------------------------------------------------------------------------------------------
					# Get DV electrons 
					#------------------------------------------------------------------------------------------
					electrons = helpers.Tracks()
					electrons.getElectrons(evt= DVevt)
					elVec = electrons.lepVec
					#------------------------------------------------------------------------------------------

					#------------------------------------------------------------------------------------------
					# Calculate the seperation of two tracks 
					#------------------------------------------------------------------------------------------
					Cosmicveto = selections.Cosmicveto(evt= DVevt)
					Cosmicvetocut = Cosmicveto.passes()

					if Cosmicvetocut: 
						if passCosmicveto == False: #only fill cut flow once per DV!!
							hcutflow.Fill(10)
						passCosmicveto = Cosmicvetocut
					else: 
						continue 
					#------------------------------------------------------------------------------------------
					
					#------------------------------------------------------------------------------------------
					# Calculate the tri-lepton mass (prompt + 2 displaced)
					#------------------------------------------------------------------------------------------
					Mlll = selections.Mlll(decayprod="mumu",plep=pMuonvec,dMu=muVec,dEl=elVec)
					Mlllcut = Mlll.passes()

					if Mlllcut: 
						if passMlll == False: #only fill cut flow once per DV!!
							hcutflow.Fill(11)
						passMlll = Mlllcut
					else: 
						continue 
					#------------------------------------------------------------------------------------------

					#------------------------------------------------------------------------------------------
					# Calculate the DV mass
					#------------------------------------------------------------------------------------------
					DVmass = selections.DVmasscut(evt= DVevt)
					DVmasscut = DVmass.passes()

					if DVmasscut:
						if passDVmasscut == False: #only fill cut flow once per DV!!
							hcutflow.Fill(12)
						passDVmasscut = DVmasscut
					else:
						continue
		


	print "Time elapsed after event selection: " , time.clock()-start
	ROOT.gStyle.SetOptStat(0)
	ROOT.gROOT.SetBatch(True)
	SetAtlasStyle()

	MyC01= ROOT.TCanvas("MyC01","cutflow",600,400)
  	MyC01.Divide(1,1)
 	MyC01.cd(1)
 	ymax_cutflow = hcutflow.GetMaximum()
 	hcutflow.GetYaxis().SetRangeUser(0,ymax_cutflow*1.05)
 	hcutflow.SetStats(0)
 	hcutflow.SetFillColor(kAzure-4)
 	hcutflow.SetLineWidth(0)
  	hcutflow.Draw("HIST TEXT0 SAME")
  	if file == "/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_20G_lt10dd_emu_wTightmuel.root":
  		helpers.drawNotes("Test without filter or mu & track quality cuts","0","20","10")	
  		MyC01.SaveAs("/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis" +'/plots/hcutflow_newEvtSel_20G_10mm_emu_wTightmuel'+'.pdf')
  	if file == "/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_10G_lt10dd_mumu_wTightmuel.root":
  		helpers.drawNotes("Test without filter","1","10","10")
  		MyC01.SaveAs("/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis" +'/plots/hcutflow_newEvtSel_10G_10mm_mumu_wTightmuel'+'.pdf')
  	


  	print "number of muons: ", nmuons 
				# print elVec[1].Pt(), elVec[1].Eta(), elVec[1].Phi()


				
				# print len(muVec)
			
			# print passTrackqual
			# print passDVtype
		# print passFid
		