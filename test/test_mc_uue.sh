makeHistograms.py \
--input root://eosuser.cern.ch//eos/user/d/dhnl/ntuple-ci/ntuples/v6p0_ntuples/mc/group.phys-exotics.mc16_13TeV.311636.DAOD_EXOT29.e7422_e5984_a875_r11891_r11748_p4482.dhnl_ntuple_v6p0.v2_tree.root/group.phys-exotics.25819678._000001.tree.root \
--config $WorkDir/data/mc/config_mc_uue.json \
--output_file $WorkDir/test/output/config_mc_uue.root \
--analysis run2Analysis \
--nevents 1000 \
--force 
