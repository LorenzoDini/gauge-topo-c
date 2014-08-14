#!/usr/bin/env python

########################################################################
# COPYRIGHT 2013 Claudio Bonati
# e-mail: claudio.bonati82@gmail.com
########################################################################
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
########################################################################

import os, glob, sys

out_form="The output is latox, latoy, latoz, latot and then\n";
out_form+="<O>, (<O^2>-<O>^2)/vol, b2, b4\n"+"for average plaquette, average Polyakov loop, top.charge without cooling\n"

try:
   sys.argv[6]
except: 
   print "Use: "+ sys.argv[0] + " -i input_files -o output_file -b block"; print out_form; sys.exit(1)

for i in range(1,6,2):
   if(sys.argv[i]==r'-i'):
      files=sys.argv[i+1]
   elif(sys.argv[i]==r'-o'):
      output=sys.argv[i+1]
   elif(sys.argv[i]==r'-b'):
      blocco=int(sys.argv[i+1])


# list of input files
filelist=glob.glob(files)

for file in filelist:
   
   #data acquisition
   in_file=open(file, 'r')
   cont_file=in_file.readlines()
   campione=len(cont_file)-1   ## sample dimension definition
   num_top_charge=(int(len(cont_file[1].split()))-5)/2  ## number of topological charges measured on each lattice
   in_file.close()

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
   os.system('g++ -Wall --pedantic -O3 ANALISI_LOC.cc -o ANALISI')
   os.system('./ANALISI input_analysis') 

   ## clean
   os.system('rm ANALISI')
   os.system('rm input_analysis') 
