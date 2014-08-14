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

import os, sys

try:
   file = sys.argv[1]
except: 
   print "Use: "+ sys.argv[0] + " `file_name`\n"; sys.exit(1)


# verifica se 'file' e' un file
if not os.path.isfile(file):
   print "Use: "+ sys.argv[0] + " `file_name`\n"; sys.exit(1)

in_file=open(file, 'r')
cont_file=in_file.readlines()
campione=len(cont_file)-1
numobs=len(cont_file[1].split())
in_file.close() 

out_file=open('input_autocorr_self','w')
out_file.write("""
datafile  %s   #data file
sample    %d   #number of lines
columns   %d   #number of columns
""" % (file, campione, numobs))
out_file.close()

os.system('g++ -Wall --pedantic -O3 AUTOCORR_SELF.cc -o AUTOCORR_SELF')
os.system('./AUTOCORR_SELF input_autocorr_self') 
os.system('rm AUTOCORR_SELF')
os.system('rm input_autocorr_self')
