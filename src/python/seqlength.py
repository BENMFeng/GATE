#!/usr/bin/env python
 
import os,sys
 
f = open(sys.argv[1],'rU')
header = f.readline()
header = header.rstrip(os.linesep)
sequence=''
for line in f:
  line = line.rstrip('\n')
  if(line[0] == '>'):
    header = header[1:]
    header = line
    print header, len(sequence)
    sequence = ''
  else:
    sequence += line
 
print header, len(sequence)