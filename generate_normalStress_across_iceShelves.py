#!/usr/bin/env python

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import subprocess
from optparse import OptionParser

parser = OptionParser(epilog='Get normal stress along specified direction. This script will generate a file named as new.nc that includes the normal stress variable called stress_direction (direction is the -d option we set). Note that only the calculation across ice shelves are correct as we use the plan-view (2HD) approximation for stress and strain rate calculation!')
parser.add_option("-i", "--input", dest="input_file", help="A nc file that includes strain rate, 3D velocity, temperature and thickness fields")
parser.add_option("-d", "--direction", dest="direction", help="direction for specified normal stress. Options are f (flow direction), af (across-flow firection), p1 (first principle stress direction), p2 (second principle stress), x (x direction) and y (y direction)")
parser.add_option("-v", "--visc", dest="viscosity", help="two options: temp (use temperature to recover) or const (set A_const in the constants section below).")

for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

#----- constants -------------
eps = 1e-30
n = 3.0
R = 8.31 # J mol-1 K-1
A01 = 3.985e-13
A02 = 1.916e3
# Pa-3 s-1
Q1 = 60e3
Q2 = 139e3
# J mol-1
rho_ocean = 1028.0
rho_ice = 910.0
gravity_acc = 9.8
A_const = 6.338e-25
#-----------------------------

stress_type = "stress_"+options.direction
subprocess.call(["cp",options.input_file,"new.nc"])
subprocess.call(["ncks","-v","thickness","new.nc", "stress.nc"])
subprocess.call(["ncks","-v","thickness","new.nc", "ux_mean.nc"])
subprocess.call(["ncks","-v","thickness","new.nc", "uy_mean.nc"])
subprocess.call(["ncks","-v","thickness","new.nc", "nx.nc"])
subprocess.call(["ncks","-v","thickness","new.nc", "ny.nc"])
subprocess.call(["ncrename","-v","thickness,"+stress_type,"stress.nc"])
subprocess.call(["ncrename","-v","thickness,ux_mean","ux_mean.nc"])
subprocess.call(["ncrename","-v","thickness,uy_mean","uy_mean.nc"])
subprocess.call(["ncrename","-v","thickness,nx","nx.nc"])
subprocess.call(["ncrename","-v","thickness,ny","ny.nc"])
subprocess.call(["ncks","-A","-v",stress_type,"stress.nc","new.nc"])
subprocess.call(["ncks","-A","-v","nx","nx.nc","new.nc"])
subprocess.call(["ncks","-A","-v","ny","ny.nc","new.nc"])
data = Dataset("new.nc",'r+')

exx = data.variables['exx'][0,:]
eyy = data.variables['eyy'][0,:]
exy = data.variables['exy'][0,:]
eyx = data.variables['eyx'][0,:]

#gh_x = data.variables['sfcElevGrad_x'][0,:]
#gh_y = data.variables['sfcElevGrad_y'][0,:]
# if we want to use the thickness gradient direction, see below for details

thickness = data.variables['thickness'][0,:]
#cellMask = data.variables['cellMask'][time,:]
#shelfIndex = np.where((cellMask&4)!=0)[0]

temperature = data.variables['temperature'][0,:,:]

ux = data.variables['uReconstructX'][0,:,:]
uy = data.variables['uReconstructY'][0,:,:]

mean_ux = np.mean(ux,axis=1)
mean_uy = np.mean(uy,axis=1)

data_ux['ux_mean'][0,:] = mean_ux
data_uy['uy_mean'][0,:] = mean_uy

meanT = np.mean(temperature,axis=1)
meanT[meanT==0] = 263.15 # for cells that contains invalid T data

A = A01*np.exp(-Q1/R/meanT)
A[meanT<=263.15] = A01*np.exp(-Q1/R/meanT[meanT<=263.15])
A[meanT>263.15] = A02*np.exp(-Q2/R/meanT[meanT>263.15])

if (options.viscosity == 'const'):
    A = A*0 + A_const
    print ("use a constant A value {0:5f}".format(A_const))
else:
    print ("calculate A using given ice temperature")


strainRate_e = np.sqrt(0.5*(exx**2+eyy**2) + 0.25*(exy+eyx)**2)
visc = 0.5*A**(-1.0/n)*(strainRate_e+eps)**((1.0-n)/n)

sigma_xx_d = 2*visc*exx
sigma_yy_d = 2*visc*eyy
sigma_xy = 2*visc*(exy+eyx)/2.0
sigma_yx = 2*visc*(eyx+exy)/2.0
# after tests, the real e_xy and e_yx equal (exy+eyx)/2 for the MALI outputs

sigma_xx = 2*sigma_xx_d+sigma_yy_d
sigma_yy = 2*sigma_yy_d+sigma_xx_d

N0 = 0.5*rho_ice*(1-rho_ice/rho_ocean)*gravity_acc*thickness

# flow direction
if options.direction == 'f':
    nx = mean_ux/np.sqrt(mean_ux**2+mean_uy**2+eps)
    ny = mean_uy/np.sqrt(mean_ux**2+mean_uy**2+eps)
    sigma_nn = ((sigma_xx*nx**2+(sigma_xy+sigma_yx)*nx*ny+sigma_yy*ny**2))
    data.variables[stress_type][0,:] = sigma_nn

# perpendicular-flow direction
elif options.direction == 'af':
    nx = mean_ux/np.sqrt(mean_ux**2+mean_uy**2+eps)
    ny = mean_uy/np.sqrt(mean_ux**2+mean_uy**2+eps)
    sigma_nn = ((sigma_xx*ny**2-(sigma_xy+sigma_yx)*nx*ny+sigma_yy*nx**2))
    data.variables[stress_type][0,:] = sigma_nn

# second principle stress direction
elif options.direction == 'p2':
    #theta1 = 0.5*np.arctan((sigma_xy+sigma_yx)/(sigma_xx-sigma_yy+eps))
    #theta2 = np.pi/2 + theta1
    #nx = np.cos(theta2)
    #ny = np.sin(theta2)
    #sigma_nn = ((sigma_xx*nx**2+(sigma_xy+sigma_yx)*nx*ny+sigma_yy*ny**2))
    sigma_nn = (sigma_xx+sigma_yy)/2.0-np.sqrt(((sigma_xx-sigma_yy)/2.0)**2+((sigma_xy+sigma_yx)/2)**2)
    data.variables[stress_type][0,:] = sigma_nn

# first principle stress direction
elif options.direction == 'p1':
    #theta1 = 0.5*np.arctan((sigma_xy+sigma_yx)/(sigma_xx-sigma_yy+eps))
    #theta2 = np.pi/2 + theta1
    #theta2 = 180.0/180*np.pi 
    #nx = np.cos(theta2)
    #ny = np.sin(theta2)
    #sigma_nn = ((sigma_xx*ny**2-(sigma_xy+sigma_yx)*nx*ny+sigma_yy*nx**2))
    sigma_nn = (sigma_xx+sigma_yy)/2.0+np.sqrt(((sigma_xx-sigma_yy)/2.0)**2+((sigma_xy+sigma_yx)/2)**2)
    data.variables[stress_type][0,:] = sigma_nn

# x direction
elif options.direction == 'x':
    nx = 1 
    ny = 0
    sigma_nn = ((sigma_xx*nx**2+(sigma_xy+sigma_yx)*nx*ny+sigma_yy*ny**2))
    data.variables[stress_type][0,:] = sigma_nn

# y direction
elif options.direction == 'y':
    nx = 1 
    ny = 0
    sigma_nn = ((sigma_xx*ny**2-(sigma_xy+sigma_yx)*nx*ny+sigma_yy*nx**2))
    data.variables[stress_type][0,:] = sigma_nn

# thickness gradient direction. This requires an extra thickness gradient variable which should be prepared by ourselves. An easy way is using Paraview to generate the gradient data, save it as a csv file and use the script called convert_csv_to_nc.py to convert it to a MPAS variable. If we want to use this option, or some other similar options, simply modify the script as you like.
#elif options.direction == 'flow_gh':
#    nx = mean_ux/np.sqrt(mean_ux**2+mean_uy**2+eps)
#    ny = mean_uy/np.sqrt(mean_ux**2+mean_uy**2+eps)
#    sigma_nn = (sigma_xx*nx+(sigma_xy+sigma_yx)/2*ny)*gh_x+(sigma_yy*ny+(sigma_xy+sigma_yx)/2*nx)*gh_y
#    data.variables[stress_type][0,:] = sigma_nn

else:
    print("wrong direction!")


#subprocess.call(["ncks","-v",stress_type,"new.nc",stress_type+".nc"])

data.close()
