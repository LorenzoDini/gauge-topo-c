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


import os, glob, sys, operator

try:
   file_dis = sys.argv[1]
   ord_colon = int(sys.argv[2])
except: 
   print "Uso: "+ sys.argv[0] + " `file_name` + column   (column>=0)"
   print "Sort the line of the file in such a way that the chosen column is growing"
   sys.exit(1)

file_ord=file_dis+'_ord'      #  sorted file
out_file=open(file_ord,'w')

in_file=open(file_dis,'r')
cont_file=in_file.readlines()

aux=[]
for i in range(0,len(cont_file),1):
  aux.append([i, (cont_file[i].split())[ord_colon]])
 
aux_sort=sorted(aux, key=operator.itemgetter(1))

for i in range(0,len(cont_file),1):
  out_file.write(cont_file[aux_sort[i][0]])
  
in_file.close()
out_file.close()
   
