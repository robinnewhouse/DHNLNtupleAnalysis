makeHistograms.py \
--input root://eosuser.cern.ch//eos/user/d/dhnl/ntuple-ci/ntuples/v6p0_ntuples/data/group.phys-exotics.data15_ntuple_v6p0.v4_tree.root/group.phys-exotics.25819340._000068.tree.root \
--config $WorkDir/data/SSbkg/config_SSbkg_uuu.json \
--output_file $WorkDir/test/output/config_data_SSbkg_uuu.root \
--analysis run2Analysis \
--nevents 10000 \
--force 
