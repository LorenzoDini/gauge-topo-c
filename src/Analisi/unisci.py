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

README="Use: "+ sys.argv[0] + " dir_name + therm\n"+"with `therm` integer number. Put together all the data files in the directiory `dir_name' which are of the form\n"+"everything_1.dat everything_200.dat\n"+"throwing away the first `therm` lines"

# Use: ./unisci.py `dir_name` therm
# with therm integer number. Put together all the data files which are of the form 
# everything_1.dat everything_200.dat
# by throwing away the first `term` lines

try:
   indir = sys.argv[1]
except: 
  print README; sys.exit(1)


try:
   therm = int(sys.argv[2])
except: 
   print README; sys.exit(1)

# add a final "/" if is not present
aux=indir[len(indir)-1:len(indir)]
if not aux == r'/':
   indir = indir+r'/'

# verify that 'indir' is a directory
if not os.path.isdir(indir):
   print README; sys.exit(1)

# list of the files *.dat of the directory
filelist=glob.glob(indir+'*.dat')

## list files which differ only for the replica number
flist=[]
for file in filelist:
   for i in range(len(file)):
      if file[i]=='_':
         i_undscore=i
   aux=file[0:i_undscore]

   if not aux in flist:
      flist.append(aux)

## merge the files
for index in flist:
  flist2=glob.glob(index+'*.dat')
  aux=index+'.dat'
  out_file=open(aux.lstrip(indir), 'w')
  i=0
  for file in flist2:
     in_file=open(file, 'r')
     cont_file=in_file.readlines()

     if len(cont_file)>0:

        #for the first file write also the first line
        if i==0: 
           out_file.write(cont_file[0])
  
        for j in range(therm+1, len(cont_file), 1):
              out_file.write(cont_file[j])
        in_file.close()

        i+=1
  out_file.close() 

