#!/usr/bin/env python

import string
import sys

i         = str(sys.argv[1])
n         = str(sys.argv[2])
dirOut    = str(sys.argv[3])
dirIn     = str(sys.argv[4])
dirTo     = str(sys.argv[5])
tagData   = str(sys.argv[6])
nameChain = str(sys.argv[7])
json      = str(sys.argv[8])
data2011  = str(sys.argv[9])
debug     = str(sys.argv[10])
mapcheck  = str(sys.argv[11])
macro     = str(sys.argv[12])

print '{'
print 'gSystem->Load("/home/llr/cms/ndaci/SKWork/macro/skEfficiency/analyzers/tagAndProbe/MakePairs/'+macro+'");'
print 'makePairs(-1,-1,' + i + ',' + n + ',"' + dirOut + '","' + dirIn + '","' + dirTo + '","' + tagData + '", "' + nameChain + '",' + json + ',"' + data2011 + '",' + debug + ',' + mapcheck + ');'
print 'gSystem->Exit(0);'
print '}'
