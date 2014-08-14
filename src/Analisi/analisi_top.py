#!/usr/bin/env python
import os, glob, sys

from numpy import *
from scipy.optimize import fmin

out_form="The output is latox, latoy, latoz, latot and then\n"
out_form+="i, <TP>, <TP^2>/vol, b2, b4, <PL>, <PL^2>/vol, b2, b4\n"
out_form+="where TP=top. charge, PL=mean plaquette and measures are taken after i*cooling\n"

try:
   sys.argv[6]
except: 
   print "Use: "+ sys.argv[0] + " -i input_file -o output_file -b block"; print out_form; sys.exit(1)

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
num_top_charge=(int(len(cont_file[1].split()))-5)/2  ## number of topological charges measured on each lattice
in_file.close()



## compute the scale factor
if not os.path.isfile("./scale_factor.dat"):

  scale=empty(num_top_charge);

  for col1 in range(0, num_top_charge, 1):
     col=5+2*col1
     #function to minimize
     def residual(s):
       aux=0.0
       for i in range(1, len(cont_file)):
         helper=float(cont_file[i].split()[col])
         aux+=(helper*s-round(helper*s, 0))*(helper*s-round(helper*s, 0))  
       return aux

     scal=1.2;
     ris=fmin(residual, scal, ftol=1.0e-4)
     print ris
     scal=ris[0]
     scale[col1]=ris[0]

  ## print the rescaling factors
  out_file=open('scale_factor.dat','w')
  for i in range(0, num_top_charge, 1):
     out_file.write("%f\n" %(scale[i]) )
  out_file.close()

## print input file
out_file=open('input_analysis','w')
out_file.write("""
datafile %s  #input file
outfile  %s  #output file
block    %d  #block for jacknife
sample   %d  #number of lines
numobs   %d  #number of observables
""" % (file, output, blocco, campione, 3+2*num_top_charge))
out_file.close()

## compile and execute
os.system('g++ -Wall --pedantic -O3 ANALISI_TOP.cc -o ANALISI')
os.system('./ANALISI input_analysis') 

## clean
os.system('rm ANALISI')
os.system('rm input_analysis') 
