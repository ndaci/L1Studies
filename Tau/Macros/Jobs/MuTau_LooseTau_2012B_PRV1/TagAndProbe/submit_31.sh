export WORKINGDIR=/home/llr/cms/ndaci/SKWork/macro/HTauTau/analyzers/tagAndProbe/Jobs/MuTau/LooseTau/Run2012B_PRV1/TagAndProbe/
cd /home/llr/cms/ndaci/WorkArea/MinCode/CMSSW_5_2_3_patch3/src/
export SCRAM_ARCH=slc5_amd64_gcc462
source /opt/exp_soft/cms/cmsset_default.sh
eval `scram runtime -sh`
cd $WORKINGDIR
touch start/submit_31.C.start
root -b macros/submit_31.C
touch done/submit_31.C.done
