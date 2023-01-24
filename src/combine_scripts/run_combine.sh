#!/bin/sh
#ulimit -s unlimited
#set -e
cd /afs/cern.ch/work/r/rymuelle/public/nanoAODzPrime/higgscombine/CMSSW_10_2_13/src
export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/r/rymuelle/public/nanoAODzPrime/CMSSW_12_1_0/src/bff_plotting_code_v3/combine_data/ 

combine -M AsymptoticLimits "$1" --run blind 

