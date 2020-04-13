python makeHistograms.py -i ../../../NtupleMaker/Robin/DHNLAlgorithm/run/eee_10GeV_10mm.root --config ../data/config_mc_eee.json --force
mv ../output/histograms_mc_eee.root ../output/histograms_mc_eee_10GeV_10mm_SS.root
python makeHistograms.py -i ../../../NtupleMaker/Robin/DHNLAlgorithm/run/eee_5GeV_10mm.root --config ../data/config_mc_eee.json --force
mv ../output/histograms_mc_eee.root ../output/histograms_mc_eee_5GeV_10mm_SS.root
python makeHistograms.py -i ../../../NtupleMaker/Robin/DHNLAlgorithm/run/eee_3GeV_10mm.root --config ../data/config_mc_eee.json --force
mv ../output/histograms_mc_eee.root ../output/histograms_mc_eee_3GeV_10mm_SS.root

# python plotHistograms.py --config ../data/config_plotting_3GeV_10mm.json
# \mv ../output/plots/selMC_mass_VSI.pdf ../output/plots/selMC_mass_VSI_3GeV_10mm.pdf
# \mv ../output/plots/selMC_mass_VSI_Leptons.pdf ../output/plots/selMC_mass_VSI_Leptons_3GeV_10mm.pdf
# \mv ../output/plots/Cutflows/CutFlow_VSI_mc.pdf ../output/plots/Cutflows/CutFlow_VSI_mc_3GeV_10mm.pdf
# \mv ../output/plots/Cutflows/CutFlow_VSI_Leptons_mc.pdf ../output/plots/Cutflows/CutFlow_VSI_Leptons_mc_3GeV_10mm.pdf
# python plotHistograms.py --config ../data/config_plotting_5GeV_10mm.json
# \mv ../output/plots/selMC_mass_VSI.pdf ../output/plots/selMC_mass_VSI_5GeV_10mm.pdf
# \mv ../output/plots/selMC_mass_VSI_Leptons.pdf ../output/plots/selMC_mass_VSI_Leptons_5GeV_10mm.pdf
# \mv ../output/plots/Cutflows/CutFlow_VSI_mc.pdf ../output/plots/Cutflows/CutFlow_VSI_mc_5GeV_10mm.pdf
# \mv ../output/plots/Cutflows/CutFlow_VSI_Leptons_mc.pdf ../output/plots/Cutflows/CutFlow_VSI_Leptons_mc_5GeV_10mm.pdf
# python plotHistograms.py --config ../data/config_plotting_10GeV_10mm.json
# \mv ../output/plots/selMC_mass_VSI.pdf ../output/plots/selMC_mass_VSI_10GeV_10mm.pdf
# \mv ../output/plots/selMC_mass_VSI_Leptons.pdf ../output/plots/selMC_mass_VSI_Leptons_10GeV_10mm.pdf
# \mv ../output/plots/Cutflows/CutFlow_VSI_mc.pdf ../output/plots/Cutflows/CutFlow_VSI_mc_10GeV_10mm.pdf
# \mv ../output/plots/Cutflows/CutFlow_VSI_Leptons_mc.pdf ../output/plots/Cutflows/CutFlow_VSI_Leptons_mc_10GeV_10mm.pdf
