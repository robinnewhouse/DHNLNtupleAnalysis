# DHNLNtupleAnalysis

`DHNLNtupleAnalysis` is a framework for reading Ntuples generated by [DHNLAlgorithm](https://gitlab.cern.ch/atlas-phys/exot/ueh/EXOT-2017-19/DHNLAlgorithm) 
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

### Running

To run the code: 

```
cd python 
python makeHistograms.py -f path_to_dHNLntuple
```

Histograms will be saved in output/histograms.root. 

By default the histograms.root file is recreated every time `makeHistograms.py` is run. If you wish to update the root file, instead of recreating it, run: 

```
python makeHistograms.py --update True -f path_to_dHNLntuple
```



## Quick Start Guide

This plotting framework is designed to output histograms for analysis selection variables in the dHNL analysis. 

The `makeHistogram.py` file is the steering code for this framework. In this file you will find a list of channels that map a 1 word key to a list of selections. To add a new channel, modify this list to include the cuts you wish to apply. When a new channel (e.g. "my_new_channel") is created and `makeHistograms.py`is run, histograms will be save in `histograms.root` with the following name format: histName_my_new_channel. 

N.B If you wish to apply a new set of cuts that are not implemented you may wish to make a new analysis class. An example can be found in the comments at the end of `analysis.py`

For a list of the currently implemented cuts see the Analysis class definition in `analysis.py`.

To add new histograms, add a new observable to the list in `observables.py`. Then fill the histogram in the corresponding function in `analysis.py`.

To add a new selection, make a new class in `selections.py`.



## Making Pretty Plots

The plan is to add `plotHistograms.py` that can format plots from `histograms.root` output file. More details to come soon!



