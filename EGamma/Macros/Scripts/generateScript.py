#!/usr/bin/env python

import string
import sys

folder = str(sys.argv[1])
macro = str(sys.argv[2])

print 'export WORKINGDIR=' + folder
print 'cd /home/llr/cms/ndaci/WorkArea/MinCode/CMSSW_4_2_8_patch7/src'
print 'source /opt/exp_soft/cms/cmsset_default.sh'
print 'eval `scram runtime -sh`'
print 'cd $WORKINGDIR'
print 'touch start/' + macro + '.start'
print 'root -b macros/' + macro
print 'touch done/' + macro + '.done'
