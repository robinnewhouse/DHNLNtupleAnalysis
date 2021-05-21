makeHistograms.py \
--input root://eosuser.cern.ch//eos/user/d/dhnl/ntuple-ci/ntuples/v5p2_ntuples/mc/group.phys-exotics.mc16_13TeV.311636.DAOD_EXOT29.e7422_e5984_a875_r11891_r11748_p4482.dhnl_ntuple_v5p2.v1_tree.root/group.phys-exotics.25605424._000001.tree.root \
--config $WorkDir/data/mc/config_mc_uue.json \
--output_file $WorkDir/test/output/config_mc_uue.root \
--analysis run2Analysis \
--nevents 1000 \
--force 
