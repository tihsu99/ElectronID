#!/bin/bash
cd /user/tomc/eleIdTuning/CMSSW_9_2_10/src/ElectronID
source $VO_CMS_SW_DIR/cmsset_default.sh
eval `scram runtime -sh`
root -b -q 'fourPointOptimization.C(false)'  &> endcap.log
