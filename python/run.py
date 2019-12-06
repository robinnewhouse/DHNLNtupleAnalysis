import treenames
import selections 
import helpers
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
	entries = 300


	for ievt in xrange(entries):

		pMuon = selections.Plepton(lepton="muon")
		pMuonpass = pMuon.passes(tree, ievt)

		# if pMuonpass:
			# print "-------"
			# print ievt
			# print pMuon.plepVec.Pt()
			# print pMuon.plepVec.Eta()
			# print pMuon.plepVec.Phi()
			


		passDV = selections.DV().passes(tree, ievt)
		

		ndv = len(tree.dvmass[ievt])

		for idv in xrange(ndv):

			passfid = selections.Fiducial(_min=4, _max=300).passes(tree, ievt,idv)
	
			passDVntracks = selections.DVntracks().passes(tree, ievt, idv)
		
			passOSDV = selections.OSDV().passes(tree, ievt, idv)

			passDVtype = selections.DVtype(decayprod="emu").passes(tree,ievt,idv)
		
	
			passTrackqual = selections.Trackqual(quality="2-tight").passes(tree,ievt,idv)
		
			if pMuonpass and passDV and passfid and passDVntracks and passOSDV and passDVtype: 
				print "-------"
				print ievt
				pMuonvec = pMuon.plepVec
				# print pMuonvec.Pt()
				# print pMuonvec.Eta()
				# print pMuonvec.Phi()

				muons = helpers.Leptons()
				muons.getMuons(tree, ievt, idv)
				muVec = muons.lepVec
				muIndex = muons.lepIndex

				electrons = helpers.Leptons()
				electrons.getElectrons(tree, ievt, idv)
				elVec = electrons.lepVec
		
				passMllcut = selections.Mlllcut(decayprod="emu",plep=pMuonvec,dMu=muVec,dEl=elVec).passes()
				# print passMllcut
				passDVmasscut = selections.DVmasscut().passes(tree, ievt, idv)
				# print passDVmasscut

				passCosmicVeto = selections.Cosmicveto().passes(tree, ievt, idv)
				print passCosmicVeto

				if passMllcut and passCosmicVeto and passDVmasscut: 

					print "muon ", muVec[0].Pt(), muVec[0].Eta(), muVec[0].Phi()
					# print muVec[1].Pt(), muVec[1].Eta(), muVec[1].Phi()
					
					print "electron ",elVec[0].Pt(), elVec[0].Eta(), elVec[0].Phi()

				# print elVec[1].Pt(), elVec[1].Eta(), elVec[1].Phi()


				
				# print len(muVec)
			
			# print passTrackqual
			# print passDVtype
		# print passfid
		