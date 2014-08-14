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

import os, glob, sys, math

try:
   indir = sys.argv[1]
except: 
   print "Uso: "+ sys.argv[0] + " `dir_name`"; sys.exit(1)

# aggiunge il carattere finale "/" se non e' gia' presente
aux=indir[len(indir)-1:len(indir)]
if not aux == r'/':
   indir = indir+r'/'

# verifica se 'indir' e' una directory
if not os.path.isdir(indir):
   print "Uso: "+ sys.argv[0] + " `dir_name`"; sys.exit(1)

# lista dei file *.dat della cartella indir
filelist=glob.glob(indir+'*.dat')

input_file=open("input", 'w')

input_file.write("outfile  datimh.dat\n")
input_file.write("betamin  XXX\n")
input_file.write("betamax  XXX\n")
input_file.write("betastep XXX\n")
input_file.write("volume   XXX\n")
input_file.write("numfiles "+str(len(filelist))+"\n")

## unisce i files
for file in filelist:
  out_file=open(file.lstrip(indir), 'w')

  in_file=open(file, 'r')
  i=-1;

  for line in in_file:
     if(i==-1):
       beta=line.split()[4]
     if(i>=0):
       aux=line.split()
       if(len(aux)==7):
         plaq=(float(aux[0])+float(aux[1]))/2.0
         poly=math.sqrt(float(aux[2])**2+float(aux[3])**2)
         out_file.write(str(plaq)+" "+str(poly)+"\n")
       else:
         print "problem with the line "+str(i)+" of file "+in_file
     i+=1
    
  out_file.close() 

  input_file.write("./MH/"+file.lstrip(indir)+"  "+beta+"  "+str(i)+"  XXX \n")


os.system('g++ -Wall --pedantic -O3 ANALISIMH.cc -o ANALISIMH')
os.system('g++ -Wall --pedantic -O3 ANALISIMH_MAXPOLY.cc -o ANALISIMH_MAXPOLY')
