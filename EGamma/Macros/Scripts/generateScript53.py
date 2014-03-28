#!/usr/bin/env python

import string
import sys

folder = str(sys.argv[1])
macro = str(sys.argv[2])

print 'export WORKINGDIR=' + folder
print 'cd /home/llr/cms/ndaci/WorkArea/HTauTau/Analysis/CMSSW_5_3_4/src/'
print 'export SCRAM_ARCH=slc5_amd64_gcc462'
print 'source /opt/exp_soft/cms/cmsset_default.sh'
print 'eval `scram runtime -sh`'
print 'cd $WORKINGDIR'
print 'touch start/' + macro + '.start'
print 'root -b macros/' + macro
print 'touch done/' + macro + '.done'
