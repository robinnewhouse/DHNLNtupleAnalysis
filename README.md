# DHNLNtupleAnalysis

`DHNLNtupleAnalysis` is a framework for reading ntuples generated by [DHNLAlgorithm](https://gitlab.cern.ch/atlas-phys/exot/ueh/EXOT-2017-19/DHNLAlgorithm) 
that was developed for the displaced heavy neutral lepton (DHNL) analysis. For installation instructions, see the "Getting Started" section below.

## Getting Started

To clone the project: 

```
setupATLAS
lsetup git
git clone ssh://git@gitlab.cern.ch:7999/atlas-phys/exot/ueh/EXOT-2017-19/DHNLNtupleAnalysis.git
```

The analysis code uses the package `uproot` to load root files. For more details see [uproot documentation](https://pypi.org/project/uproot/). To setup python to include the `uproot` package on lxplus do the following: 

```
source /cvmfs/sft.cern.ch/lcg/views/LCG_latest/x86_64-centos7-gcc9-opt/setup.sh
```

If you are running on a local cluster the complier version might change (i.e. gcc9 -> gcc8). Alternatively you can setup python and install the uproot package locally yourself. 

### Running makeHistograms.py

To make histograms from DHNL ntuples do the following: 
```
cd python 
python makeHistograms.py -i path_to_dHNLntuple --config ../data/config_mc_all.json
```
One root file per channel will be saved in DHNLNtupleAnalysis/output/ folder. 


If you only want to run on a single channel (e.g. uuu channel) then run: 
```
cd python 
python makeHistograms.py -i path_to_dHNLntuple --config ../data/config_mc_uuu.json
```
For this example, the code will save the file DHNLNtupleAnalysis/output/histograms_uuu.root. To see a full list of channels see this [section](#list-of-configuration-files).


If you have previously created a histogram file for the channel you are trying to run, you can run the following to overwrite the output file: 
```
python makeHistograms.py --force -i path_to_dHNLntuple --config ../data/config_mc_all.json
```
Without the `--force` option, the code will not run to prevent you from accidentally overwriting your histogram files!


If you would like to run over a different set of analysis cuts (e.g using the ToyAnalysis cuts) you can run: 
```
python makeHistograms.py -i path_to_dHNLntuple --config ../data/config_mc_uuu.json --analysis ToyAnalysis
```

## A note about MC weighting

Event weighting for monte carlo samples is implmented in the framework. In order to properly do the weighting the mass and lifetime of the sample you are running is required to be in the name of the input string in the same format as the DAOD_RPVLL container names. If you see a warning such as:

```
"Can't determine the mass and lifetime of signal sample. MC weight will be set to 1!!"
```
this is becuase your file is not appropriately named.Either rename your ntuple file following the convention from the DAOD_RPVLL conatiner name you used to make the ntuple or make due without MC event weighting. See the list of DAOD_RPVLL samples [here](https://twiki.cern.ch/twiki/pub/AtlasProtected/ExoticLongLivedHeavyNeutralLeptonRel21/MC16a_MC16d_MC16e_dHNL_DAOD_RPVLLonly_corr_new.txt) for naming conventions.


### Running plotHistograms.py

Once you have run makeHistograms.py edit `../data/config_plotting.json` to include a path to the histogram files you want to plot, as well as an identifying label for each file.

Then run:
```
python plotHistograms.py --config ../data/config_plotting.json
```

The plotHisotograms.py code will be user specific depending on what histograms you want to plot. Please feel free to use the plotting functions in plotting.py in your own plotting scripts. Or use the example plotHistograms.py that is provided. 



## Quick Start Guide

This ntuple analysis code is designed to output histograms for analysis selection variables in the dHNL analysis. 

The `makeHistogram.py` file is the steering code for making histograms. To run the code, a configuration file is needed that provides the analysis code with a list of selections associated to a given channel name. For a list of the supported channel names see [List of configuration files](#list-of-configuration-files).

To add a new channel, edit the corresponding config file and make a new channel that includes the cuts you wish to apply. 

For example if you edit `../data/config_mc_uuu.json` and add "my_new_channel", when you run `makeHistograms.py -i inputFile --config ../data/config_mc_uuu.json`, a new histograms file will be created called `histograms_my_new_channel.root` in the output directory. This output file will include histograms with your new cuts. 

N.B If you wish to apply a new set of cuts that are not implemented you may wish to make a new analysis class. See Toy Analysis Class in `analysis.py` as an example.

For a list of the currently implemented cuts see the Analysis class definition in `analysis.py`.

To add new histograms, add a new observable to the list in `observables.py`. Then fill the histogram in the corresponding function in `analysis.py`.

To add a new selection, make a new class in `selections.py`.


## List of configuration files

The following are the default channels that have corresponding config files in DHNLNtupleAnalysis/data: 
- uuu: Prompt muon, 2-displaced muons
- ueu: Prompt muon, 1-displaced muon & 1-displaced electron 
- uee: Prompt muon, 2-displaced electrons
- euu: Prompt electron, 2-displaced muons
- eeu: Prompt electron, 1-displaced muon & 1-displaced electron 
- eee: Prompt electron, 2-displaced electrons



## Making Pretty Plots

Plotting functions are defined in `plotting.py`. The steering code `plotHistograms.py` defines what kinds of plots to save in the default output directory `DHNLNtupleAnalysis/output/plots/`. Running plotHistograms.py you can plot different distribtions on the same canvas and save them to the output directory.

Make sure you edit `DHNLNtupleAnalysis/data/config_plotting.json` to include the files you wish to plot.

For making cutflow plots: 

1. open `plotHistogram.py`
2. edit the function makeCutflows() following the example provided. 
3. run `plotHistograms.py` and the output will be in `DHNLNtupleAnalysis/output/plots/Cutflows`

For comparing data and MC distributions: 
1. open `plotHistogram.py`
2. edit compareMCdata() following the examples provided. 
3. run `plotHistograms.py` and the output will be in `DHNLNtupleAnalysis/output/plots/`

For comparing N histograms in the same file: 
1. open `plotHistogram.py`
2. edit compareHistograms() following the examples provided. 
3. run `plotHistograms.py` and the output will be in `DHNLNtupleAnalysis/output/plots/`

