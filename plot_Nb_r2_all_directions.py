#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from scipy import stats
from matplotlib.ticker import PercentFormatter
from optparse import OptionParser

parser = OptionParser(epilog='plot Nb-Nr correlation and Nb for all directions')
parser.add_option("-i", "--input", dest="input_file", help="input nc file")
parser.add_option("-x", "--index", dest="index_file", help="perturbcell_index file")
parser.add_option("-c", "--coord", dest="coord_file", help="perturbcell_coord file")
parser.add_option("-b", "--butCR", dest="butCR_file", help="buttressingChangeRatio file")

for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]" 
options, args = parser.parse_args()

#--------- constants ---------------
eps = 1.0e-30
n = 3.0
rho_i = 910.0
rho_w = 1028.0
g = 9.8
A = 6.338e-25
#-----------------------------------

thkGrad_x = np.loadtxt('../../data/gh_x.txt')
thkGrad_y = np.loadtxt('../../data/gh_y.txt')
thkGrad = np.sqrt(thkGrad_x**2+thkGrad_y**2)

index_P = np.loadtxt(options.index_file)
perturbN = len(index_P[:,0])
print ("{0:5d} cells perturbed".format(perturbN))
index_P = np.int16(index_P[0:perturbN,1])
data = Dataset(options.input_file,'r')

coord = np.loadtxt(options.coord_file)
xcoord_P = coord[0:perturbN,0]
ycoord_P = coord[0:perturbN,1]

exx = data.variables['exx'][0,:]
eyy = data.variables['eyy'][0,:]
exy = data.variables['exy'][0,:]
eyx = data.variables['eyx'][0,:]

thickness = data.variables['thickness'][0,:]

ux = data.variables['uReconstructX'][0,:,:]
uy = data.variables['uReconstructY'][0,:,:]

mean_ux = np.mean(ux,axis=1)
mean_uy = np.mean(uy,axis=1)

nx_flow = mean_ux/np.sqrt(mean_ux**2+mean_uy**2+eps)
ny_flow = mean_uy/np.sqrt(mean_ux**2+mean_uy**2+eps)


strainRate_e = np.sqrt(0.5*(exx**2+eyy**2) + 0.25*(exy+eyx)**2)
visc = 0.5*A**(-1.0/n)*(strainRate_e+eps)**((1.0-n)/n)

sigma_xx_d = 2*visc*exx
sigma_yy_d = 2*visc*eyy
sigma_xy = 2*visc*(exy+eyx)/2.0
sigma_yx = 2*visc*(eyx+exy)/2.0

sigma_xx = 2*sigma_xx_d+sigma_yy_d
sigma_yy = 2*sigma_yy_d+sigma_xx_d
        
index_like_gh = np.intersect1d(np.argwhere(thkGrad<7e-3), np.argwhere(xcoord_P<480e3))

butCR = np.loadtxt(options.butCR_file)
metricButCR = butCR**(1/n) 

r2_all = np.zeros(180)
r2_all1 = np.zeros(180)
theta1_mean = np.zeros(180)

theta_flow = np.arccos(nx_flow)
theta_principle = 0.5*np.arctan((sigma_xy+sigma_yx)/(sigma_xx-sigma_yy+eps))
theta_diff = theta_flow-theta_principle
theta_diff_P = theta_diff[index_P]*180.0/np.pi

for i in range(180):

    theta = i/180.0*np.pi

    theta_flow = np.arccos(nx_flow)

    theta_principle = 0.5*np.arctan(2*sigma_xy/(sigma_xx-sigma_yy))

    nx = np.cos(theta_flow+theta)
    ny = np.sin(theta_flow+theta)

    nx1 = np.cos(theta_principle+theta)
    ny1 = np.sin(theta_principle+theta)

    sigma_nn = sigma_xx*nx**2+(sigma_xy+sigma_yx)*nx*ny+sigma_yy*ny**2
    sigma_nn1 = sigma_xx*nx1**2+(sigma_xy+sigma_yx)*nx1*ny1+sigma_yy*ny1**2

    sigma_nn_P = sigma_nn[index_P]
    np.savetxt('stress_flow.txt',sigma_nn_P)
    sigma_nn_P1 = sigma_nn1[index_P]

    N0_P = 0.5*rho_i*(1-rho_i/rho_w)*g*thickness[index_P]

    sigma_nn_P_like = sigma_nn_P[index_like_gh]
    sigma_nn_P1_like = sigma_nn_P1[index_like_gh]
    N0_P_like = N0_P[index_like_gh]

    x=1-sigma_nn_P_like/N0_P_like
    x1=1-sigma_nn_P1_like/N0_P_like
    y=metricButCR[index_like_gh]
    mask = ~np.isnan(x) & ~np.isnan(y)
    mask1 = ~np.isnan(x1) & ~np.isnan(y)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x[mask],y[mask])
    slope, intercept, r_value1, p_value, std_err = stats.linregress(x1[mask1],y[mask1])

    r2_all[i] = r_value*r_value
    r2_all1[i] = r_value1*r_value1


plt.figure (1,figsize=(12,5))
plt.subplot(121)

plt.plot(np.arange(180), r2_all)
plt.plot(np.arange(180), r2_all1, 'r')
plt.xlabel('$\Delta\phi$ ($^\circ$)')
plt.ylabel('$r^2$')
plt.xlim(0,180)
#plt.ylim(0.3,1)
plt.legend(['flow','1st principle stress'])

plt.subplot(122)
#plt.plot(np.arange(180), r2_all1)
#plt.xlabel('$\Delta\phi$ ($^\circ$)')
#plt.xlim(0,180)
#plt.ylim(0.3,1)
plt.hist(theta_diff_P, bins=10, weights=np.ones(len(theta_diff_P))/len(theta_diff_P))
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.xlabel('$\Delta\phi$ ($^\circ$)')
plt.ylabel('percentange (%)')

plt.savefig('fig2.pdf')

plt.show()
