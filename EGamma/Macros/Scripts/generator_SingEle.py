#!/usr/bin/env python

import string
import sys

i = str(sys.argv[1])
n = str(sys.argv[2])
dirOut = str(sys.argv[3])
dirIn = str(sys.argv[4])
tagData = str(sys.argv[5])
nameChain = str(sys.argv[6])
json = str(sys.argv[7])
data2011 = str(sys.argv[8])
debug = str(sys.argv[9])
macro = str(sys.argv[10])

print '{'
#print 'gSystem->Load("/home/llr/cms/ndaci/SKWork/macro/skEfficiency/analyzers/tagAndProbe/NormalOnly/makePairs_C.so");'
#print 'gSystem->Load("/home/llr/cms/ndaci/SKWork/macro/skEfficiency/analyzers/tagAndProbe/NormalOnly/'+macro+'");'
print 'gSystem->Load("/home/llr/cms/ndaci/SKWork/macro/skEfficiency/analyzers/tagAndProbe/MakePairs/'+macro+'");'
print 'makePairs(-1,' + i + ',' + n + ',"' + dirOut + '","' + dirIn + '","' + tagData + '", "' + nameChain + '",' + json + ',"' + data2011 + '",' + debug + ');'
print 'gSystem->Exit(0);'
print '}'
