export WORKINGDIR=/home/llr/cms/ndaci/SKWork/macro/skEfficiency/analyzers/tagAndProbe/Jobs/Study2012/Run2012D_PRV1_up/
cd /home/llr/cms/ndaci/WorkArea/HTauTau/Analysis/CMSSW_5_3_4/src/
export SCRAM_ARCH=slc5_amd64_gcc462
source /opt/exp_soft/cms/cmsset_default.sh
eval `scram runtime -sh`
cd $WORKINGDIR
touch start/submit_54.C.start
root -b macros/submit_54.C
touch done/submit_54.C.done
