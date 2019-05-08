#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
from optparse import OptionParser

parser = OptionParser(epilog='Generate the buttressing flux response number (theta_B) defined in Reese et al. (2018). This will create a file named as buttressingChangeRatio.txt in the current folder.')
parser.add_option("-i", "--index", dest="index_file", help="perturbcell_index file")
parser.add_option("-a", "--area", dest="area_file", help="perturbcell_area file")
parser.add_option("-m", "--method", dest="method", help="two options: mpas or albany")
parser.add_option("-f", "--folder", dest="folder", help="folder directory (name). For example, if the ctrl run folder name is ../test_ctrl, and the perturbation folder name is ../test_i, then set -f ../test")

for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]" 
options, args = parser.parse_args()

#-------------- constants --------------------
rho_i = 910.0
#---------------------------------------------


if options.method == "mpas":
    data_ctrl = Dataset(options.folder+'_ctrl'+'/globalStats.nc','r')
    glFlux_ctrl = data_ctrl.variables['groundingLineFlux'][1]

elif options.method == "albany":
    data_ctrl = open(options.folder+'_ctrl'+'/log.albany.0000.out','r')
    last_line = data_ctrl.readlines()[1151]
    glFlux_ctrl = float(last_line.split()[-1])
else:
    print ("wrong method")

perturbIndex = np.loadtxt(options.index_file)
area = np.loadtxt(options.area_file)

perturbN = len(area[:,0])

print ("{0:5d} cells perturbed".format(perturbN))

buttressingChangeRatio = np.zeros(perturbN)

if options.method == "mpas":
    for i in range(perturbN):
        dirN = options.folder+'_'+str(i+1)
        data_i = Dataset(dirN+'/globalStats.nc','r')
        glFlux_i = data_i.variables['groundingLineFlux'][1]
        area_i = area[i,1]
        glFlux_change = glFlux_i - glFlux_ctrl
        perturbMass_change = area_i*1*rho_i

        buttressingChangeRatio[i] = glFlux_change/perturbMass_change*100
        data_i.close()

elif options.method == "albany":
    for i in range(perturbN):
        dirN = options.folder+'_'+str(i+1)
        data_i = open(dirN+'/log.albany.0000.out','r')
        last_line = data_i.readlines()[1151]
        glFlux_i = float(last_line.split()[-1])
        glFlux_change = glFlux_i - glFlux_ctrl
        area_i = area[i,1]
        perturbMass_change = area_i*1*rho_i

        buttressingChangeRatio[i] = glFlux_change/perturbMass_change*100
        data_i.close()
else:
    print ("wrong method")

data_ctrl.close()

np.savetxt('buttressingChangeRatio.txt',buttressingChangeRatio)

