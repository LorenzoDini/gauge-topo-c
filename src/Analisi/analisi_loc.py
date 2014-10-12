#!/usr/bin/env python
import os, glob, sys

out_form="The output is latox, latoy, latoz, latot and then\n"
out_form+="<x>, V(<x^2>-<x>^2), 1-<x^4>/3<x^2>^2\n"
out_form+="where x = mean plaquette and polyakov loop\n"

try:
   sys.argv[6]
except: 
   print("Use: "+ sys.argv[0] + " -i input_file -o output_file -b block"); print(out_form); sys.exit(1)

for i in range(1,6,2):
   if(sys.argv[i]==r'-i'):
      file=sys.argv[i+1]
   elif(sys.argv[i]==r'-o'):
      output=sys.argv[i+1]
   elif(sys.argv[i]==r'-b'):
      blocco=int(sys.argv[i+1])



#data acquisition
in_file=open(file, 'r')
cont_file=in_file.readlines()
campione=len(cont_file)-1   ## sample dimension definition
in_file.close()



## print input file
out_file=open('input_analysis','w')
out_file.write("""
datafile %s  #input file
outfile  %s  #output file
block    %d  #block for jacknife
sample   %d  #number of lines
numobs   2  #number of observables
""" % (file, output, blocco, campione))
out_file.close()

## compile and execute
os.system('g++ -Wall --pedantic -O3 ANALISI_LOC.cc -o ANALISI')
os.system('./ANALISI input_analysis') 

## clean
os.system('rm ANALISI')
os.system('rm input_analysis') 
