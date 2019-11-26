import treenames
import selections 
import helpers
from pprint import pprint


def dump(obj):
  for attr in dir(obj):
    print("obj.%s = %r" % (attr, getattr(obj, attr)))

if __name__ == '__main__':
	
	file ="/home/dtrischuk/HNLAnalysis/DHNLNtupleAnalysis/mc16_13TeV.311633.Pythia8EvtGen_A14NNPDF23LO_WmuHNL50_10G_lt10dd.merge.DAOD_RPVLL.e7422_e5984_a875_r10790_r10726_newVariables.root"
	
	treename = "outTree"
	tree = treenames.Tree(file, treename)

	entries = len(tree.dvmass)

	for ievt in xrange(entries):
		passPlepton = selections.Plepton(lepton="muon").passes(tree, ievt)
		print passPlepton

		passDV = selections.DV().passes(tree, ievt)
		
		passfid = selections.Fiducial(_min=4, _max=300).passes(tree, ievt)

		ndv = len(tree.dvmass[ievt])

		for idv in xrange(ndv):
			passDVntracks = selections.DVntracks().passes(tree, ievt, idv)
			passOSDV = selections.OSDV().passes(tree, ievt, idv)
			passDVtype = selections.DVtype(decayprod="ee").passes(tree,ievt,idv)
			passTrackqual = selections.Trackqual(quality="2-tight").passes(tree,ievt,idv)

			if passTrackqual: 
				muVec = helpers.Leptons().getMuons(tree, ievt, idv)
				elVec = helpers.Leptons().getElectrons(tree, ievt, idv)


				
				print len(muVec)
			
			# print passTrackqual
			# print passDVtype
		# print passfid
		