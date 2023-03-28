DirForSubs=/nfs/dust/atlas/user/appelt/R22/data_SS_2018_eee_VSI_Lep_Mod_2
outputName=data_SS_2018.root


if test -f "$DirForSubs/ntuples/$outputName"; then
    echo "combined file already exists"
    return 1
fi

echo "adding files"
hadd $DirForSubs/ntuples/$outputName $DirForSubs/ntuples/output_*
sleep 2s
echo "cleaning up"
rm $DirForSubs/ntuples/output_*
rm -rf $DirForSubs/logs
rm $DirForSubs/executable-*
rm $DirForSubs/submission-*
echo "done, have fun with $DirForSubs/ntuples/$outputName"