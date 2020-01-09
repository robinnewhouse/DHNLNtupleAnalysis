import treenames	
import selections 
import helpers
import ROOT, time
from ROOT import *
gROOT.LoadMacro("AtlasStyle.C")
gROOT.LoadMacro("AtlasUtils.C")
gROOT.LoadMacro("AtlasLabels.C")


# emu file: 
# file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/fixElmatching_mc16_13TeV.311660.Pythia8EvtGen_A14NNPDF23LO_WmuHNL50_20G_lt10dd_el.merge.DAOD_RPVLL.e7422_e5984_a875_r10739_r10706.root"
# file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_20G_lt10dd_emu_wTightmuel.root"
# file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_20G_lt10dd_emu_VSILep.root"


# mumu file:
# file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_10G_lt10dd_mumu.root"
# file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_10G_lt10dd_mumu_wTightmuel.root"
# file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/rootfiles/newframework_Ntuple_WmuHNL_10G_lt10dd_mumu_VSILep.root"



def Event_Sel(fileName, dsid, plepton, DV_type, mass, lifetime, MC_campaign, VtxConfig):
	ROOT.gStyle.SetOptStat(0)
	ROOT.gROOT.SetBatch(True)
	SetAtlasStyle()

	file = fileName
	treename = "outTree"
	tree = treenames.Tree(file, treename)

	nentries = len(tree.dvmass)
	# nentries = 50

	hcutflow= ROOT.TH1D("hcutflow","cutflow",13,-0.5,12.5)
	

	hcutflow.GetXaxis().SetBinLabel(1, "all")
	# hcutflow.GetXaxis().SetBinLabel(2, "HNL mu26")
	hcutflow.GetXaxis().SetBinLabel(2, "all single lep triggers")
	hcutflow.GetXaxis().SetBinLabel(3, "4-filters")
	if plepton == "muon":
		hcutflow.GetXaxis().SetBinLabel(4, "tight pmu")
	if plepton == "electron":
		hcutflow.GetXaxis().SetBinLabel(4, "tight pel")	
	hcutflow.GetXaxis().SetBinLabel(5, "DV")  
	hcutflow.GetXaxis().SetBinLabel(6, "fiducial")  
	hcutflow.GetXaxis().SetBinLabel(7, "2-track DV")  
	hcutflow.GetXaxis().SetBinLabel(8, "2 mu DV") 
	hcutflow.GetXaxis().SetBinLabel(9, "1-tight")  
	hcutflow.GetXaxis().SetBinLabel(10, "2-tight")  

	# hcutflow.GetXaxis().SetBinLabel(8, "OS DV")  
	# if DV_type == 'mumu':
	# 	hcutflow.GetXaxis().SetBinLabel(9, "2-muon DV")  
	# if DV_type == 'emu':
	# 	hcutflow.GetXaxis().SetBinLabel(9, "mu-el DV")  
	# hcutflow.GetXaxis().SetBinLabel(10, "2-tight-lepton DV")  
	# hcutflow.GetXaxis().SetBinLabel(11, "cosmic veto")    
	# hcutflow.GetXaxis().SetBinLabel(12, "m_{lll}")      
	# hcutflow.GetXaxis().SetBinLabel(13, "mDV") 

	hd0= ROOT.TH1D("hd0_%s_%sG_%smm"%(dsid,mass,lifetime),"hd0_%s_%sG_%smm"%(dsid,mass,lifetime),200,-10,10)
	hd0.GetXaxis().SetTitle("d0 [mm]")
	hd0.GetYaxis().SetTitle("selected events")



	
	npasstrig = 0
	npassfilter = 0 
	npasssel = 0 

	for ievt in xrange(nentries):


		passTrigger = False
		passHNLfilter = False
		passPmuon = False
		passDV = False
		passFid = False
		passDVntracks = False
		passOSDV = False
		passDVtype = False
		passTrackqual = False
		passTrackqual_2 = False
		passCosmicveto = False
		passMlll = False
		passDVmasscut = False

		evt = helpers.Event(tree=tree, ievt = ievt , idv = None)
		ndv = len(tree.dvx[ievt])

		hcutflow.Fill(0)

		#------------------------------------------------------------------------------------------
		# Check if event passes trigger
		#------------------------------------------------------------------------------------------
		if plepton == "muon":
			Trigger = selections.Trigger(evt = evt, plepton="muon",trigger = "HLT_mu26_ivarmedium")
		if plepton == "electron": 
			Trigger = selections.Trigger(evt = evt, plepton="electron",trigger = "HLT_e26_lhtight_nod0_ivarloose")
		Triggercut = Trigger.passes()

		if Triggercut: 
			hcutflow.Fill(1)
			npasstrig = npasstrig + 1
			passesTrigger = Triggercut
		else:
			continue 	
		# #------------------------------------------------------------------------------------------


		#------------------------------------------------------------------------------------------
		# Check if event passes HNL filter 
		#------------------------------------------------------------------------------------------
		# if plepton == "muon": 
		# 	HNLfilter = selections.Filter(evt= evt, _filter="mu-el")
		# if plepton == "electron": 
		# 	HNLfilter = selections.Filter(evt= evt, _filter="el-mu")

		HNLfilter = selections.Filter(evt= evt, _filter="4-filter")

		HNLfiltercut = HNLfilter.passes()
		if HNLfiltercut: 
			npassfilter = npassfilter + 1
			hcutflow.Fill(2)
			passHNLfilter = HNLfiltercut
		else:
			continue 	
		# #------------------------------------------------------------------------------------------


		#------------------------------------------------------------------------------------------
		# Find prompt muon for event 
		#------------------------------------------------------------------------------------------
		if plepton == "muon": 
			Pmuon = selections.Plepton(evt = evt, lepton="muon")
		if plepton == "electron": 
			Pmuon = selections.Plepton(evt = evt, lepton="electron")

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
		# preSel = Triggercut and Pmuoncut and nDVcut
		# preSel = Triggercut and nDVcut
		preSel = nDVcut
		# preSel = passTrigger and passHNLfilter and  passPmuon and passDV
		#------------------------------------------------------------------------------------------


		if preSel: #only continue if pre-selection is True
			for idv in xrange(ndv):

				DVevt = helpers.Event(tree=tree, ievt = ievt , idv = idv)
				for itr in range(len(tree.trackpt[ievt][idv])): 
					hd0.Fill(tree.trackd0[ievt][idv][itr])


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
				# Count if DV has the correct number of leptons 
				#------------------------------------------------------------------------------------------
				# print "##############################################"
				# print "##############################################"
				# print "##############################################"
				# print "##############################################"
				DVtype = selections.DVtype(evt= DVevt,decayprod="mumu")

				DVtypecut = DVtype.passes()
				if DVtypecut:
					if passDVtype == False:  #only fill cut flow once per DV!!
						hcutflow.Fill(7)
					passDVtype = DVtypecut
				else:
					continue
				#------------------------------------------------------------------------------------------

			
				# #------------------------------------------------------------------------------------------
				# # Check charge of the tracks in DV 
				# #------------------------------------------------------------------------------------------
				# OSDV = selections.OSDV(evt= DVevt)
				# OSDVcut = OSDV.passes()
				# if OSDVcut:
				# 	if passOSDV == False:  #only fill cut flow once per DV!!
				# 		hcutflow.Fill(7)		
				# 	passOSDV = OSDVcut	
				# else:
				# 	continue
				# #------------------------------------------------------------------------------------------		
				# # Count if DV has the correct number of leptons 
				# #------------------------------------------------------------------------------------------
				# if DV_type == "mumu": 
				# 	DVtype = selections.DVtype(evt= DVevt,decayprod="mumu")
				# if DV_type == "emu": 
				# 	DVtype = selections.DVtype(evt= DVevt,decayprod="emu")
				# if DV_type == "ee":
				# 	DVtype = selections.DVtype(evt= DVevt,decayprod="ee")

				# DVtypecut = DVtype.passes()
				# if DVtypecut:
				# 	if passDVtype == False:  #only fill cut flow once per DV!!
				# 		hcutflow.Fill(8)
				# 	passDVtype = DVtypecut
				# else:
				# 	continue
				# #------------------------------------------------------------------------------------------
				# # count the number of lepton tracks with a given quality 
				# #------------------------------------------------------------------------------------------
				Trackqual = selections.Trackqual(evt=DVevt, quality="1-tight")
				Trackqualcut = Trackqual.passes()
				if Trackqualcut: 
					if passTrackqual == False:
						hcutflow.Fill(8)
					passTrackqual = Trackqualcut
				# else: 
				# 	continue 

				Trackqual_2 = selections.Trackqual(evt=DVevt, quality="2-tight")
				Trackqualcut_2 = Trackqual_2.passes()
				if Trackqualcut_2: 
					if passTrackqual_2 == False:
						hcutflow.Fill(9)
					passTrackqual_2 = Trackqualcut_2
				# else: 
				# 	continue 




				# #------------------------------------------------------------------------------------------

		
				# if Fidvolcut and DVntrackscut and OSDVcut and DVtypecut and Trackqual: # only do cosmicveto & mass cuts if you have a good 2 lep DV
				# 	#------------------------------------------------------------------------------------------
				# 	# Get prompt muon
				# 	#------------------------------------------------------------------------------------------
				# 	pMuonvec = Pmuon.plepVec
				# 	#------------------------------------------------------------------------------------------

				# 	#------------------------------------------------------------------------------------------
				# 	# Get DV muons 
				# 	#------------------------------------------------------------------------------------------
				# 	muons = helpers.Tracks()
				# 	muons.getMuons(evt= DVevt)
				# 	muVec = muons.lepVec
				# 	#------------------------------------------------------------------------------------------

				# 	#------------------------------------------------------------------------------------------
				# 	# Get DV electrons 
				# 	#------------------------------------------------------------------------------------------
				# 	electrons = helpers.Tracks()
				# 	electrons.getElectrons(evt= DVevt)
				# 	elVec = electrons.lepVec
				# 	#------------------------------------------------------------------------------------------

				# 	#------------------------------------------------------------------------------------------
				# 	# Calculate the seperation of two tracks 
				# 	#------------------------------------------------------------------------------------------
				# 	Cosmicveto = selections.Cosmicveto(evt= DVevt)
				# 	Cosmicvetocut = Cosmicveto.passes()

				# 	if Cosmicvetocut: 
				# 		if passCosmicveto == False: #only fill cut flow once per DV!!
				# 			hcutflow.Fill(10)
				# 		passCosmicveto = Cosmicvetocut
				# 	else: 
				# 		continue 
				# 	#------------------------------------------------------------------------------------------
					
				# 	#------------------------------------------------------------------------------------------
				# 	# Calculate the tri-lepton mass (prompt + 2 displaced)
				# 	#------------------------------------------------------------------------------------------
				# 	if DV_type == "mumu": 
				# 		Mlll = selections.Mlll(decayprod="mumu",plep=pMuonvec,dMu=muVec,dEl=elVec)
				# 	if DV_type == "emu": 
				# 		Mlll = selections.Mlll(decayprod="emu",plep=pMuonvec,dMu=muVec,dEl=elVec)
				# 	if DV_type == "ee":
				# 		Mlll = selections.Mlll(decayprod="ee",plep=pMuonvec,dMu=muVec,dEl=elVec)

				# 	Mlllcut = Mlll.passes()
				# 	if Mlllcut: 
				# 		if passMlll == False: #only fill cut flow once per DV!!
				# 			hcutflow.Fill(11)
				# 		passMlll = Mlllcut
				# 	else: 
				# 		continue 
				# 	#------------------------------------------------------------------------------------------

				# 	#------------------------------------------------------------------------------------------
				# 	# Calculate the DV mass
				# 	#------------------------------------------------------------------------------------------
				# 	DVmass = selections.DVmasscut(evt= DVevt)
				# 	DVmasscut = DVmass.passes()

				# 	if DVmasscut:
				# 		if passDVmasscut == False: #only fill cut flow once per DV!!
				# 			npasssel = npasssel + 1
				# 			hcutflow.Fill(12)
				# 		passDVmasscut = DVmasscut
				# 	else:
				# 		continue


	print npasstrig, npassfilter, nentries

	MyC01= ROOT.TCanvas("MyC01","cutflow",600,400)
  	MyC01.Divide(1,1)
 	MyC01.cd(1)
 	# ROOT.gPad.SetLogy()
 	ymax_cutflow = hcutflow.GetMaximum()
 	hcutflow.GetYaxis().SetRangeUser(0,ymax_cutflow*1.05)
 	hcutflow.SetStats(0)
 	hcutflow.SetFillColor(kAzure-4)
 	hcutflow.SetLineWidth(0)
  	hcutflow.Draw("HIST TEXT0 SAME")
  	helpers.drawNotes(MC_campaign,DV_type,mass,lifetime,plepton,VtxConfig)

  	if mass == "-1": 
  		if plepton == "muon": 
	  		# MyC01.SaveAs('/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis' +'/plots/hcutflow_WmuN_1-filter_%s_%s_%s'%(MC_campaign, DV_type,VtxConfig)+'.pdf')
	  		MyC01.SaveAs('/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis' +'/plots/hcutflow_WmuN_2muon_%s_%s_%s'%(MC_campaign, DV_type,VtxConfig)+'.pdf')

	  	if plepton == "electron":
	  		MyC01.SaveAs('/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis' +'/plots/hcutflow_WelN_1-filter_%s_%s_%s'%(MC_campaign, DV_type,VtxConfig)+'.pdf')
  	else: 
	  	if plepton == "muon": 
	  		MyC01.SaveAs('/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis' +'/plots/hcutflow_WmuN_%sG_%smm_%s_%s_allcuts'%(mass, lifetime, DV_type,VtxConfig)+'.pdf')
	  	if plepton == "electron":
	  		MyC01.SaveAs('/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis' +'/plots/hcutflow_WelN_%sG_%smm_%s_%s_allcuts'%(mass, lifetime, DV_type,VtxConfig)+'.pdf')


  	MyC02= ROOT.TCanvas("MyC02","cutflow",600,400)
  	MyC02.Divide(1,1)
 	MyC02.cd(1)
 	ymax_cutflow = hd0.GetMaximum()
 	hd0.GetYaxis().SetRangeUser(0,ymax_cutflow*1.05)
 	hd0.SetStats(0)
 	# hd0.SetFillColor(kAzure-4)
 	# hd0.SetLineWidth(0)
  	hd0.Draw()
  	helpers.drawNotes(MC_campaign,DV_type,mass,lifetime,plepton,VtxConfig)

  	if mass != "-1": 
	  	if plepton == "muon": 
	  		MyC02.SaveAs('/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis' +'/plots/d0_WmuN_%sG_%smm_%s_%s_alltracks'%(mass, lifetime, DV_type,VtxConfig)+'.pdf')
	  	if plepton == "electron":
	  		MyC02.SaveAs('/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis' +'/plots/hd0_WelN_%sG_%smm_%s_%s_alltracks'%(mass, lifetime, DV_type,VtxConfig)+'.pdf')




  
