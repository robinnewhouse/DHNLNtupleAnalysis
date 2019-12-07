import treenames
import selections 
import helpers
import ROOT
from ROOT import *
gROOT.LoadMacro("AtlasStyle.C")
gROOT.LoadMacro("AtlasUtils.C")
gROOT.LoadMacro("AtlasLabels.C")
from pprint import pprint



def dump(obj):
  for attr in dir(obj):
    print("obj.%s = %r" % (attr, getattr(obj, attr)))

if __name__ == '__main__':
	
	# file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/mc16_13TeV.311633.Pythia8EvtGen_A14NNPDF23LO_WmuHNL50_10G_lt10dd.merge.DAOD_RPVLL.e7422_e5984_a875_r10790_r10726_newVariables.root"
	file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/fixElmatching_mc16_13TeV.311660.Pythia8EvtGen_A14NNPDF23LO_WmuHNL50_20G_lt10dd_el.merge.DAOD_RPVLL.e7422_e5984_a875_r10739_r10706.root"
	
	treename = "outTree"
	tree = treenames.Tree(file, treename)

	entries = len(tree.dvmass)
	# entries = 300

	hcutflow= ROOT.TH1D("hcutflow","cutflow",13,-0.5,12.5)
  	hcutflow.GetXaxis().SetBinLabel(1, "all")
  	hcutflow.GetXaxis().SetBinLabel(2, "trigger")
  	hcutflow.GetXaxis().SetBinLabel(3, "filter")
  	hcutflow.GetXaxis().SetBinLabel(4, "pmu")
  	hcutflow.GetXaxis().SetBinLabel(5, "DV")  
  	hcutflow.GetXaxis().SetBinLabel(6, "fiducial")  
  	hcutflow.GetXaxis().SetBinLabel(7, "2-track DV")  
  	hcutflow.GetXaxis().SetBinLabel(8, "OS DV")  
  	# if DV_type == '1':
  		# hcutflow.GetXaxis().SetBinLabel(9, "2-muon DV")  
  	# if DV_type == '0':
  	hcutflow.GetXaxis().SetBinLabel(9, "mu-el DV")  
  	hcutflow.GetXaxis().SetBinLabel(10, "2-tight-lepton DV")  
  	hcutflow.GetXaxis().SetBinLabel(11, "cosmic veto")    
  	hcutflow.GetXaxis().SetBinLabel(12, "m_{lll}")      
  	hcutflow.GetXaxis().SetBinLabel(13, "mDV") 




	for ievt in xrange(entries):
		passTrigger = False
		passFilter = False
		passPmuon = False
		passDV = False
		passFid = False
		passDVntracks = False
		passOSDV = False
		passDVtype = False
		passCosmicveto = False
		passMlllcut = False
		passDVmasscut = False





		hcutflow.Fill(0)
		passTrigger = selections.Trigger(plepton="muon",trigger = "HLT_mu26_ivarmedium").passes(tree, ievt)
		if passTrigger: 
			hcutflow.Fill(1)
		else:
			continue 	

		passFilter = selections.Filter(_filter="mu-mu").passes(tree, ievt)

		# if passFilter: 
		# 	hcutflow.Fill(2)
		# else:
		# 	continue 	
		# print passFilter
		pMuon = selections.Plepton(lepton="muon")
		passPmuon = pMuon.passes(tree, ievt)

		if passPmuon: 
			hcutflow.Fill(3)
		else:
			continue 	

		# if passPmuon:
			# print "-------"
			# print ievt
			# print pMuon.plepVec.Pt()
			# print pMuon.plepVec.Eta()
			# print pMuon.plepVec.Phi()
			


		passDV = selections.DV().passes(tree, ievt)
		if passDV: 
			hcutflow.Fill(4)
		else:
			continue 
		
		preSel = passTrigger and passPmuon and passDV

		# print passTri
		# print 
		# print preSel

		ndv = len(tree.dvmass[ievt])
		if preSel: 
			for idv in xrange(ndv):
				
				Fid = selections.Fiducial(_min=4, _max=300).passes(tree, ievt,idv)
				if Fid:
					if passFid == False:  #only fill cut flow once per DV!!
						hcutflow.Fill(5)
					passFid = selections.Fiducial(_min=4, _max=300).passes(tree, ievt,idv)
					
				else: 
					continue

					# if ievt == 217: 
					# 	print "-----"
					# 	print ievt
					# 	print tree.trackcharge[ievt][idv][0]
					# 	print tree.trackcharge[ievt][idv][1]
					# 	print tree.trackpt[ievt][idv][0], tree.tracketa[ievt][idv][0],tree.trackphi[ievt][idv][0]
					# 	print tree.trackpt[ievt][idv][1], tree.tracketa[ievt][idv][1],tree.trackphi[ievt][idv][1]
				
				DVntracks = selections.DVntracks().passes(tree, ievt, idv)
				if DVntracks:
					if passDVntracks == False:  #only fill cut flow once per DV!!
						hcutflow.Fill(6)
					passDVntracks = selections.DVntracks().passes(tree, ievt, idv)
					
				else:
					continue
				# if passDVntracks: 
				# 	hcutflow.Fill(6)
				# else:
				# 	continue 
				OSDV = selections.OSDV().passes(tree, ievt, idv)
				if OSDV:
					if passOSDV == False:  #only fill cut flow once per DV!!
						hcutflow.Fill(7)
					passOSDV = selections.OSDV().passes(tree, ievt, idv)
					
				else:
					continue
				# if passOSDV: 
				# 	hcutflow.Fill(7)
				# else:
				# 	continue 
				DVtype = selections.DVtype(decayprod="emu").passes(tree,ievt,idv)
				if DVtype:
					if passDVtype == False:  #only fill cut flow once per DV!!
						hcutflow.Fill(8)
					passDVtype = selections.DVtype(decayprod="emu").passes(tree,ievt,idv)
				else:
					continue
				# if passDVtype: 
				# 	hcutflow.Fill(8)
				# else:
				# 	continue 
			
				# if selections.Trackqual(quality="2-tight").passes(tree,ievt,idv):
				# 	passTrackqual = selections.Trackqual(quality="2-tight").passes(tree,ievt,idv)
				# if passTrackqual: 
				# 	hcutflow.Fill(9)
				# else:
				# 	continue 
				# print passDVtype
				if Fid and DVntracks and OSDV and DVtype: 
					# print "HELLO"
			
					pMuonvec = pMuon.plepVec
					# print pMuonvec.Pt()
					# print pMuonvec.Eta()
					# print pMuonvec.Phi()

					muons = helpers.Leptons()
					muons.getMuons(tree, ievt, idv)
					muVec = muons.lepVec
					# print muVec
					muIndex = muons.lepIndex

					electrons = helpers.Leptons()
					electrons.getElectrons(tree, ievt, idv)
					elVec = electrons.lepVec
					# print elVec
					Cosmicveto = selections.Cosmicveto().passes(tree, ievt, idv)
					if Cosmicveto: 
						if passCosmicveto == False: #only fill cut flow once per DV!!
							hcutflow.Fill(10)
						passCosmicveto = Cosmicveto
					
					else: 
						continue 

					# if passDVtype: 
					# 	hcutflow.Fill(10)
					# else:
					# 	continue 
			
					Mlllcut = selections.Mlllcut(decayprod="emu",plep=pMuonvec,dMu=muVec,dEl=elVec).passes()
				
					if Mlllcut: 
						# print Mlllcut
						if passMlllcut == False: #only fill cut flow once per DV!!
							hcutflow.Fill(11)
						passMlllcut = Mlllcut
					else: 
						continue 

					# if passMllcut: 
					# 	hcutflow.Fill(11)
					# else:
					# 	continue 
			
					# print passMllcut
					DVmasscut = selections.DVmasscut().passes(tree, ievt, idv)
					if DVmasscut:
						if passDVmasscut == False: #only fill cut flow once per DV!!
							hcutflow.Fill(12)
						passDVmasscut = DVmasscut
						hcutflow.Fill(12)
					else:
						continue
					# if passDVmasscut: 
					# 	hcutflow.Fill(12)
					# else:
					# 	continue 
					# print passDVmasscut

					

					# if passMllcut and passCosmicVeto and passDVmasscut: 
						# print "-------"
						# print ievt
						# print "trigger pass: ", passTrigger
						# print "muon ", muVec[0].Pt(), muVec[0].Eta(), muVec[0].Phi()
						# # print muVec[1].Pt(), muVec[1].Eta(), muVec[1].Phi()
						
						# print "electron ",elVec[0].Pt(), elVec[0].Eta(), elVec[0].Phi()

		# if passTrigger: 
		# 	hcutflow.Fill(1)
		# 	if passPmuon: 
		# 		# print "--------"
		# 		# print ievt
		# 		hcutflow.Fill(3)
		# 		if passDV: 
		# 			hcutflow.Fill(4)
		# 			if passFid: 
		# 				# if ievt == 217: 
		# 				# 	print "-----"
		# 				# 	print ievt
		# 					# print tree.trackcharge[ievt][idv][0]
		# 					# print tree.trackcharge[ievt][idv][1]
		# 					# print tree.trackpt[ievt][idv][0], tree.tracketa[ievt][idv][0],tree.trackphi[ievt][idv][0]
		# 					# print tree.trackpt[ievt][idv][1], tree.tracketa[ievt][idv][1],tree.trackphi[ievt][idv][1]
		# 				hcutflow.Fill(5)
		# 				if passDVntracks: 
		# 					hcutflow.Fill(6)
		# 					if passOSDV: 
		# 						# print "--------"
		# 						# print ievt
		# 						hcutflow.Fill(7)
		# 						if passDVtype: 
		# 							hcutflow.Fill(8)
		# 							if passDVtype: 
		# 								hcutflow.Fill(10)
		# 								if passMlllcut: 
		# 									hcutflow.Fill(11)
		# 									if passDVmasscut: 
		# 										hcutflow.Fill(12)

	MyC01= ROOT.TCanvas("MyC01","cutflow",600,400)
  	MyC01.Divide(1,1)
 	MyC01.cd(1)
 	ymax_cutflow = hcutflow.GetMaximum()
 	hcutflow.GetYaxis().SetRangeUser(0,ymax_cutflow*1.05)
 	hcutflow.SetStats(0)
 	hcutflow.SetFillColor(kAzure-4)
 	hcutflow.SetLineWidth(0)
  	hcutflow.Draw("HIST TEXT0 SAME")
  	# drawNotes("Test without filter or mu/track ID","0","20","10")	
  	MyC01.SaveAs("/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis" +'/hcutflow_newEvtSel'+'.pdf')

				# print elVec[1].Pt(), elVec[1].Eta(), elVec[1].Phi()


				
				# print len(muVec)
			
			# print passTrackqual
			# print passDVtype
		# print passFid
		