makeHistograms.py \
--input root://eosuser.cern.ch//eos/user/d/dhnl/ntuple-ci/ntuples/v5p0_ntuples/data/group.phys-exotics.data18_ntuple_v5p0.v2_tree.root/group.phys-exotics.25411062._000215.tree.root \
--config $WorkDir/data/SSbkg/config_SSbkg_uuu.json \
--output_file $WorkDir/test/output/config_data_SSbkg_uuu.root \
--analysis run2Analysis \
--nevents 10000 \
--force 
