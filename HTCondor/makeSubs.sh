#!/bin/bash
#
# make the submission scripts and submit them to the HTCondor batch system
# creates executable.sh and submission.sh for each data file
# submission directory has to be AFS, data can be read from EOS
# questions: christian.appelt@cern.ch

# define input data directory
dataDirectory=/nfs/dust/atlas/user/appelt/R22/data/2018/user.cappelt.data18_13TeV.periodAllYear.physics_Main.PhysCont.DAOD_LLP1.grp18_v01_p5332_tree.root

# define directory for submission files and output ntuples, logs (will be created)
DirForSubs=/nfs/dust/atlas/user/appelt/R22/data_SS_2018_uuu_VSI_Lep_Mod_2

# directory of NTupleAnalysis, change makeHistorgrams.py configuration in line 40 (output in line 42 has to be renamed as well)
realPWD=/nfs/dust/atlas/user/appelt/R22


submitJobs="y"
i=1

mkdir $DirForSubs
mkdir $DirForSubs/ntuples
mkdir $DirForSubs/logs

for eachfile in "$dataDirectory"/* #loop over files in data directory
 do
   echo "processing file $eachfile"
   file=$DirForSubs/executable-$i.sh
   file1=$DirForSubs/submission-$i.sub

   echo "#!/bin/bash" >> "$file"

   echo "inputFile=$eachfile" >> "$file"

   echo "tempDir=\$(mktemp -d -p /tmp/)" >> "$file"
   echo "cd \$tempDir" >> "$file"
   echo "cp -r $realPWD/DHNLNtupleAnalysis/* ." >> "$file"
   echo "mkdir build" >> "$file"
   echo "cd build" >> "$file"
   echo "export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase" >> "$file"
   echo "source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh" >> "$file"
   echo "asetup AnalysisBase,master,latest" >> "$file"
   echo "cmake .." >> "$file"
   echo "make" >> "$file"
   echo "source x*/setup.sh" >> "$file"

   echo "cd ../python" >> "$file"
   echo "python makeHistograms.py -i \$inputFile --config ../data/SSbkg/config_SSbkg_uuu.json --output_file ../output/output.root --useNTAU -f"  >> "$file"

   echo "cp -r ../output/output.root $DirForSubs/ntuples/output_$i.root"  >> "$file"
   echo "cd $realPWD/"  >> "$file"
   echo "rm -rf \$tempDir"  >> "$file"
   chmod +x "$file"

   echo "executable = $file" >> $file1
   echo "should_transfer_files = Yes" >> $file1
   echo "when_to_transfer_output = ON_EXIT" >> $file1
   echo "output     = $DirForSubs/logs/$i.out" >> $file1
   echo "error      = $DirForSubs/logs/$i.err" >> $file1
   echo "log        = $DirForSubs/logs/$i.log" >> $file1
   echo "universe   = vanilla" >> $file1
   echo "+JobFlavour= \"microcentury\"" >> $file1
   echo "RequestCpus = 1" >> $file1
   echo "queue" >> $file1
   chmod +x $file1

   # submit
   if [ $submitJobs = "y" ]; then
       condor_submit $file1
   fi
   (( i = i + 1 ))
 done

