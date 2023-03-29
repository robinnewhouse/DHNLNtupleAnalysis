# DHNLNtupleAnalysis (R22 dev branch)

Work in progress version for R22 samples, only tested for uuu mc20e so far. We will update the readme once the output is validated.
Therefore, use with caution ;)



`DHNLNtupleAnalysis` is a framework for reading ntuples generated by [DHNLAlgorithm](https://gitlab.cern.ch/atlas-phys/exot/ueh/EXOT-2017-19/DHNLAlgorithm) 
that was developed for the displaced heavy neutral lepton (DHNL) analysis. For installation instructions, see the "Getting Started" section below.

## Getting Started

To clone the project: 
```
setupATLAS
lsetup git
git clone --recursive  https://@gitlab.cern.ch:8443/atlas-phys/exot/ueh/EXOT-2017-19/DHNLNtupleAnalysis.git
```

There are two ways of setting up the I/O, depending on your preference for the I/O backend 

### Using NTupleAnalysisUtils (faster)  

The analysis code uses the package `NTupleAnalysisUtils` to load root files. For more details see [NTupleAnalysisUtils documentation](https://gitlab.cern.ch/Atlas-Inner-Tracking/NtupleAnalysisUtils_tutorial). To setup to run, use a recent release 22 analysis release, for example: 

To set up the environment (from scratch): 
```
mkdir build; cd build; asetup AnalysisBase,22.0.105 ; cd - 
```

To build the C++ backend (only needed the first time or when changing NtupleAnalysisUtils): 
```
cd build; 
cmake ../source; 
make ; 
source x*/setup.sh
```
This will make NTupleAnalysisUtils available, while the python scripts in this package will only need the python and ROOT versions that come with the release.

NB: The python release that comes with AnalysisBase doesn't support uproot3. Do not try to configure both uproot and AnalysisBase/NTAU at the same time. The releases will clash and nothing will run.

### Using uproot (historical default)

The analysis code historically uses the package `uproot` to load root files. This is the default. For more details see [uproot documentation](https://pypi.org/project/uproot/). To setup python to include the `uproot` package on lxplus do the following: 

```
source /cvmfs/sft.cern.ch/lcg/views/LCG_96c_LS/x86_64-centos7-gcc8-opt/setup.sh
```

If you are running on a local cluster the complier version might change (i.e. gcc9 -> gcc8). Alternatively you can setup python and install the uproot package locally yourself. 

### Running makeHistograms.py

To make histograms from DHNL ntuples for a specific channel (eg. uuu) do the following: 
```
cd python 
python makeHistograms.py -i path_to_dHNLntuple --config ../data/mc/config_mc_uuu.json [[--useNTAU ]]
```
One root file will be saved in DHNLNtupleAnalysis/output/ folder. 

Make sure your input file name contains the channel, mass, lifetime and mc campaign. Like this:
`MadGraph_uuu_10GeV_10mm_mc20e.root`.

For this example, the code will save the file DHNLNtupleAnalysis/output/histograms_mc16x_uuu.root. To see a full list of channels see this [section](#configs). (The x will be either a,d or e depending on what nutple was given as an input).

For a full list of configurable options for `makeHistograms.py` see this list of [options](#list-of-makehistogramspy-options).

## Quick Start Guide

This ntuple analysis code is designed to output histograms for analysis selection variables in the dHNL analysis. 

The `makeHistogram.py` file is the steering code for making histograms. To run the code, a configuration file is needed that provides the analysis code with a list of selections associated to a given channel name. For a list of the supported channel names see [List of configuration files](#configs).

To add a new channel, edit the corresponding config file and make a new channel that includes the cuts you wish to apply. 

For example if you edit `../data/mc/config_mc_uuu.json` and add "my_new_channel", when you run `makeHistograms.py -i inputFile --config ../data/mc/config_mc_uuu.json`, a new histograms file will be created called `histograms_my_new_channel.root` in the output directory. This output file will include histograms with your new cuts. 

N.B If you wish to apply a new set of cuts that are not implemented you may wish to make a new analysis class. See Toy Analysis Class in `analysis.py` as an example.

For a list of the currently implemented cuts see the Analysis class definition in `analysis.py`.

To add new histograms, add a new observable to the list in `observables.py`. Then fill the histogram in the corresponding function in `analysis.py`.

To add a new selection, make a new class in `selections.py`.


## Weights

**1. MC event weight**

MC event weights are computed in the `MC_event_weight` class in `helpers.py`. One weight is computed for each model in the interpretation.

In order properly calculate the cross section for the signal model, the mass and lifetime of the sample you are running is required to be in the name of the input string in the same format as the DAOD_RPVLL container names. If you see a warning such as:
```
"Can't determine the mass and lifetime of signal sample. MC mass-lifetime weight will be set to 1!!"
```
this is becuase your file is not appropriately named. Either rename your ntuple file following the convention from the DAOD_RPVLL conatiner name you used to make the ntuple or make due without MC event weighting. See the list of DAOD_RPVLL samples [here](https://twiki.cern.ch/twiki/pub/AtlasProtected/ExoticLongLivedHeavyNeutralLeptonRel21/MC16a_MC16d_MC16e_dHNL_WmuHNL_DAOD_RPVLLonly_updatedRECO.txt) for naming conventions.

A full lists of models and the descripton of what weights and trees should be used can be found [here](#HNL-models).

**2. Truth weight for spin correlations and leptons ordering bug fix**

Truth reweighting was introducted to add spin correlations with Pythia did not take into acount when simulating the truth distributions. This reweighting also fixed a lepton ordering bug that was also found in our official Pythia MC samples. For more details about these effects please see [this talk](https://indico.cern.ch/event/944478/contributions/3968769/attachments/2102276/3534533/MC_Reweighting_Sept15_20.pdf) and [this talk](https://indico.cern.ch/event/944479/contributions/3968779/attachments/2106444/3542592/MC_Reweighting_Sept22_20.pdf). 

In `selections.py` a class called `MCEventType` is used to determine if the current event is an LNC or LNV decay and accordingly calcute a weight based on the truth distributions from the particular event.

An overall "model event weight" is calculated in `analysis.py` the function called `calculate_event_weight`. The model weight is the product of both the MC event weight and the LNC/V weight.


**2. Scale factors**

Event-level scale factors have been implemented in the framework (also see [Systematics](## Systematics) selction below). The nominal scale factor is the product of three weights: 

`scale_factor = self.lepton_reco_sf['nominal'] * self.lepton_trig_sf['nominal'] * self.tree['weight_pileup']`

This scale factor is saved in mini-tree as the variable `SF_nominal`. 

**Note:** Some scale factors default to 1 and are only updated once the corresponding selection has been made. For example, the lepton trigger scale factor is only defined after the trigger matched lepton is selected in the `trigger_matched_medium_lepton_cut`. This means that in order to save the correct scale factor in the mini tree you must save the mini-trees after this cut. 


## Configs

Various configs are included in the `data/` folder. Here you will find various configs that are set up to run event selection for: 
- HNL signal samples (mc)
- inverted prompt lepton validation region (CR)
- inverted mhnl control region (inverted_mhnl)
- inverted tri-lepton mass region (inverted_mlll)
- kshort analysis (for tracking studies)
- same-sign background selection (SSbkg)
- unblinded signal region (unblind_SR)

Within the different folders there are configs that are setup to run the event selection for the six signal channels: 

- uuu: Prompt muon, 2-displaced muons
- uue: Prompt muon, 1-displaced muon & 1-displaced electron 
- uee: Prompt muon, 2-displaced electrons
- euu: Prompt electron, 2-displaced muons
- eeu: Prompt electron, 1-displaced muon & 1-displaced electron 
- eee: Prompt electron, 2-displaced electrons

## Systematics 

Systematics types are broadly separated into scale factor and tree systematics.

- SF systematics can be handled by providing a different weight for each scale factor.
- Tree systematics have different quantities for physics variables (muon pt, etc.) so the selections can be different. This requires a different tree for each.

The systematics can be run using --doSystematics as a command-line option.

A directory structure is used to separate different systematics. The nominal directory contains the tree without systematic variations. The other variations are kept in directories such as nominalMUON_ID__1down (see pictures in comment below). Only the nominal tree contains the scale factor variation ntuples (e.g. nominalMUON_ID__1down does not contain the ntuple MUON_EFF_RECO_SYS__1down).

The event weight must be multiplied by a scale factor. The nominal scale factor is SF_nominal but for systematic variations these would be MUON_EFF_RECO_SYS__1down etc. So each variable must be weighted by the event weight and the scale factor. e.g. HNLm * DV_weight_LNC_only * SF_nominal for the nominal variation and HNLm * DV_weight_LNC_only * MUON_EFF_RECO_SYS__1down for the muon reconstruction systematic variation.

An example file of 1000 events can be found here: [https://cernbox.cern.ch/index.php/s/GN0CtC4HrRnzPFc](https://cernbox.cern.ch/index.php/s/GN0CtC4HrRnzPFc).


## Guide To Using The Mini-Ntuples
Mini-ntuples are saved to the output file when you run `makeHistograms.py`. Mini-ntuples are designed to store the full tree information (not just the binned histogram) for each variable. Having access to a micro-ntuple means you can quickly re-bin or plot correlations after some selections are applied. By default the mini-ntuples are saved after the mHNL cut (the fianl SR selection). 

Inside of the histogram output file you will find the mini-ntuple trees in: `nominal/VSI_Leptons_Mod/ntuples_DVtype_VSI_Leptons_Mod`


If you want to change the cut to a different point in the event selection so that you can study the impact of the selection in data and mc (e.g. save ntuples after the `DVtype` cut), you can run: 

```
python makeHistograms.py -i path_to_dHNLntuple --config ../data/mc/config_mc_uuu.json --saveNtuples DVtype
```

Or if you want to save mini-trees after every cut: 

```
python makeHistograms.py -i path_to_dHNLntuple --config ../data/mc/config_mc_uuu.json --saveNtuples allcuts
```

N.B This "allcuts" option makes the histogram output files much larger, especially when running on data. But the feature if available if you would like to use it. 

Here is a list of cuts that you can update the code using the `-s` or `--saveNtuples` option: 

- all 
- DVtype
- cosmic
- DVlep_pt
- matveto (only for e-ee and u-ee channels)
- trkqual
- trig_match
- mvis
- mDV
- Zmass_veto
- mHNL

## HNL Models
In the interpretation of the analysis there are six different models. These models vary as to whether both LNC and LNV decays are allowed and what mixing angles are allowed and how strong the mixing to the various flavours are. 

**Mixing Models**
- **Single-flavour mixing**: HNL only mixes with muon or electron neutrinos 
  - Depending on mixing (electron or muon) only certain channels are non-zero. 
  - Muon-only mixing weights are saved for channels with prompt muons 
  - Electron-only mixing weights are saved for channels with prompt electrons 
- **Inverted heirarchy (IH) mixing**: Equal mixing with all three flavours of neutrinos
- **Normal heirarchy (NH) mixing**: Roughly equal mixing with muon and tau neutrinos. Mixing with electron neutrinos is supressed.

**HNL Models**
- **One Dirac HNL**: one HNL with only LNC decays.
- **One Majorana HNL**: one HNL with 50% LNC and 50% LNV decays.
- **Quasi-Dirac pair**: two HNLs depending on the mass splitting between the two HNLs, LNC and LNV decays have different contributions.
  - "Majorana limit": limit where HNL decays 50% LNC and 50% LNV.
  - "Dirac limit"   : limit where HNL decays 100% LNC.


To correctly compute the predicted number of events in the signal region for the different models you must select _both_ the correct model weight and the correct mini-tree. A summary is as follows: 

1. **One Dirac HNL with single-flavour mixing**: LNC mini-tree + `model_weight_one_dirac_hnl_LNC_single_flavour_mixing`
2. **Quasi-Dirac pair _"Dirac limit"_ with IH mixing**: LNC mini-tree + `model_weight_quasi_dirac_pair_LNC_ih_mixing`
3. **Quasi-Dirac pair _"Dirac limit"_ with NH mixing**: LNC mini-tree + `model_weight_quasi_dirac_pair_LNC_nh_mixing`
4. **One Majorana HNL with single-flavour mixing**: LNC_plus_LNV mini-tree + `model_weight_quasi_dirac_pair_LNCplusLNV_ih_mixing`
5. **Quasi-Dirac pair _"Majorana limit"_ with IH mixing**: LNC_plus_LNV mini-tree + `model_weight_quasi_dirac_pair_LNCplusLNV_ih_mixing`
6. **Quasi-Dirac pair _"Majorana limit"_ with NH mixing**: LNC_plus_LNV mini-tree + `model_weight_quasi_dirac_pair_LNCplusLNV_nh_mixing`

**Combinations** 

Applying the weights to the mini-trees listed above in each channel will give you the singal contribution in each respective channels. To get an overall combined signal strength for the different models, you will need to combine the various channels together. Here is a summary of the channels that should be combined for the various models: 

1. **One Dirac HNL with single-flavour mixing**: 
  - muon-only mixing: u-uu + u-ue + u-ee
  - electron-only mixing: e-ee + e-eu + e-uu
2. **Quasi-Dirac pair _"Dirac limit"_ with IH mixing**: u-uu + u-ue + u-ee + e-ee + e-eu + e-uu
3. **Quasi-Dirac pair _"Dirac limit"_ with NH mixing**: u-uu + u-ue + u-ee + e-ee + e-eu + e-uu
4. **One Majorana HNL with muon-only mixing**: 
  - muon-only mixing: u-uu + u-ue + u-ee
  - electron-only mixing: e-ee + e-eu + e-uu
5. **Quasi-Dirac pair _"Majorana limit"_ with IH mixing**: u-uu + u-ue + u-ee + e-ee + e-eu + e-uu
6. **Quasi-Dirac pair _"Majorana limit"_ with NH mixing**: u-uu + u-ue + u-ee + e-ee + e-eu + e-uu

These combinations assume that the accpetance to tau channels is zero.

## Study the displaced vertex efficiency using different VSI configurations

By default, the DHNLAlgorithm will only run on the `VSI_Leptons_Mod` container, which uses the custom VSI configuration to run the displaced vertex reconstruction. 
However, if your DHNL ntuples has been produced with different vertex containers you may be able to run the event selection to study the efficiency of the differen vertex configurations.
To run with a different vertex configuration, add the name of the new configuration to the `vtx_containers` list in the `DHNLNtupleAnalysis` configs located in 'data/'
For example, to run on VSI, VSI_Leptons and VSI_Leptons_Mod modify the config like this:
```
{
"OS_ee":{
        "vtx_containers" :  ["VSI", "VSI_Leptons", "VSI_LeptonsMod"],
        "selections":["CR", "nDV", "fidvol","2track","OS", "ee","2-veryveryloose","cosmicveto","DVmass" ]
        }
}
```

### Output File Structure 

When you run makeHistograms.py on a dHNL nutple you will get an output file located in the `DHNLNtupleAnalysis/ouput/` folder. This ouput file will contain variaous folders full of different histograms that were filled when you ran the event selection cuts in `analysis.py`. When you run on HNL signal samples you should get the following folder structure: 

```
├── nominal
│    ├── VSI_Leptons_Mod_ntuples_LNC_plus_LNV_mHNL <-- mini-ntuples from VSI Lepton Mod DVs that are saved after the mHNL cut with both LNC and LNV decays
│    ├── VSI_Leptons_Mod_ntuples_LNC_mHNL <-- mini-ntuples from VSI Lepton Mod DVs that are saved after the mHNL cut with LNC decays only
│    ├── VSI_Leptons_Mod_ntuples_LNV_mHNL <-- mini-ntuples from VSI Lepton Mod DVs that are saved after the mHNL cut with LNV decays only
│    ├── VSI_Leptons_Mod <-- folder that contains all histograms filled by VSI DVs
│       ├──truth <-- folder containing all truth histograms
│           ├──all <-- all truth histograms
│               ├──LNC <--LNC truth histograms 
│               └──LNV <-- LNV truth histograms
|           └── presel <-- truth histograms after preselection 
│       ├── all <-- histograms from all events and all DVs
│       ├── DVtype <-- histograms from events and all DVs after DVtype selection is applied
│           ├──LNC <-- LNC truth histograms 
│           └──LNV <-- LNV truth histograms
│       ... (more folders the same structure for all the differnt cuts)
│       └── CutFlow <-- VSI Leptons Mod cutflows
```


## List of makeHistograms.py Options

| **Option** | **Action** |
| ---------- | ---------- |
| `--config` | input config file for makeHisotgrams.py (required) |
| `-i` or `--input` | path to input DHNL ntuple that was produced with DHNL Algorithm (required) |
| `-o` or `--output` | output directory to store histograms. |
| `-f` or `--force` | overwrite previous histograms output file if it exists. (default: False) |
| `-a` or `--analysis` | name of the analysis you want to run (default: Full run 2 analysis) |
| `-s` or `--saveNtuples` | name of cut after which you want to save the micro-ntuples. (default: DVtype) |
| `-d` or `--debug`| Debug level. Options included are CRITICAL, ERROR, WARNING, INFO, DEBUG (default: INFO) |
| `--output_file` | Overrides the output filename and output directory |
| `--nevents` | Max. number of events to be processed |
| `--skipEvents` | Start processing the DHNL ntuple file a this event index. Event index starts a 0 and runs to n_max_events -1. |
| `--weight` | Use this flag to override the dHNL signal weight calculation for this sample. |
| `--notHNLmc` | Use this flag when running on mc that is not an HNL signal file. This will turn off HNL specific selections (e.g. HNL truth infomation.). |
| `--doSystematics` | Enable the filling of the systematic trees. Warning: due to the large number of systematics enabling doSystematics will make the hisotgram run run about 16 times more slowly.|



## Making Plots

Please see `DHNLPlotting` repository for making plots. [DHNLPlotting](https://gitlab.cern.ch/dtrischu/dhnlplotting/-/tree/master).

## Note About Colour Logs

This analysis code uses a python package called `coloredlogs` to output colourful logs to easily flag warning and error messages. To check if your python environment has the `coloredlogs` package installed run: 

```
python
>>> import coloredlogs
```

If you see the error message `ImportError: No module named coloredlogs`, then you can try to install the package via pip (or however you usually install python pacakges): 

```
pip install coloredlogs
```

## Running DHNL Algorithm on HTCondor 
If you are interested in running `DHNLNtupleAnalysis` on the lxplus batch system, Christian has put together a few scripts for submitting jobs on the HT condor system. You can check out his scripts [here](https://gitlab.cern.ch/cappelt/hnlntupleanalysis_htcondor).


# Shuffled Background Estimate
The DHNL analysis uses an object shuffling method to estimate the background from two leptons that randomly cross. The data-driven background model is obtained from a sample of "shuffled events."
This sample is created by combining each opposite-sign (OS) displaced vertex (DV) in the validation region (VR) with each prompt lepton found in a non-VR event that contains a same-sign (SS) DV satisfying loose requirements: DV mass > 1 GeV with no lepton-identification criteria imposed on its displaced leptons.
For each channel, the shuffled sample has at least 2,000 times the number of events in the "unshuffled" data sample, in which the DV and prompt lepton are from the same event.
As with unshuffled events, each shuffled event has a tri-lepton mass and HNL mass that is computed using the three objects (prompt lepton + two displaced leptons).

The shuffled background estimate can be run using the `shuffled_background_estimate.py` script. This script requires as input mini-ntuples for both SS events from the SR and OS/SS events from the VR. These mini ntuples are created by the DHNLNtupleAnalysis script (see instructions above). When running the script you must select a channel (e.g. uuu). This will then select the corresponding flavour of prompt lepton and DVs from the data events to create the sample of shuffled events. 

To run the main shuffled background estimate:
```
python shuffled_background_estimate.py --channel {channel} --run --VR_charge OS --select_plep_from_loose_regionB --select_DV_from_VR
```

To run the DV systematic: 
```
python shuffled_background_estimate.py --channel {channel} --run --VR_charge SS --select_plep_from_loose_regionB --select_DV_from_VR
```

To run the prompt lepton systematic: 
```
python shuffled_background_estimate.py --channel {channel} --run --VR_charge OS --select_plep_from_loose_regionB --select_DV_from_VR --select_plep_from_diff_channel

```

The output of the shuffled background estimate in a mini-ntuple in the same format as the main DHNLNtupleAnalysis mini-ntuple containing the relevent information about the shuffled events (i.e. not all mini-ntuple tress will be avaliable, but the necessary ones for the events selection/ statistical analysis code will be avaliable).


## List of shuffled_background_estimate.py Options

| **Option** | **Action** |
| ---------- | ---------- |
| `--channel` | The channel you want to produce the shuffled sample for (e.g. uuu, uue, eeu, eee, uee or euu). Required. |
| `--run ` | make the shuffled event sample. Default: False |
| `-p` or `--plot` | Plot the output of the shuffled background estimate. This will make validation plots as well.  |
| `--VR_charge` | The charge selection for the DV in the Validation region. Default: "OS" (Alternatively choose "SS") |
| `--select_plep_from_loose_regionB` | Select prompt leptons from the looser prompt lepton SR. This option will give you more prompt leptons to shuffle with and increase the background estimate statistics. Default: True |
| `--use_loose_DVs` | This flag will use a loose DV selection. This option is included for testing purposes only. Default: False |
| `--select_DV_from_VR`| Select DVs from the Validation Region. By default the code will select SR SS DVs. Default:False |
| `--select_plep_from_diff_channel` | This option will select prompt leptons from a different SV type event (e.g. u-ue SS events will be used to select prompt muons for u-uu channel. )This option is used to run prompt lepton systematic. Default: False |
| `--shuffle_PV` | This flag will also mix PV location, which changes the HNL flight direction used to compute the mhnl. This option is included for testing purposes only. Default: False.  |