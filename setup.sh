export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
source /cvmfs/sft.cern.ch/lcg/views/LCG_96c_LS/x86_64-centos7-gcc8-opt/setup.sh
lsetup git

# Installing python directory
# get directory of setup.sh script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export WorkDir=$DIR
# add python to PYTHONPATH and PATH variable if not already there
[[ ":$PYTHONPATH:" != *":$DIR/python:"* ]] && PYTHONPATH="$DIR/python:${PYTHONPATH}"
[[ ":$PATH:" != *":$DIR/python:"* ]] && PATH="$DIR/python:${PATH}"
