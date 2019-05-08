import numpy as np
from netCDF4 import Dataset
import subprocess
from optparse import OptionParser

parser = OptionParser(epilog='create thickness perturbed geometries. The thickness is perturbed cell by cell for a specified region. The script will generate many perturbed files (landice_grid_xxx.nc) and three geometry text files (index.txt, area.txt and coord.txt) in the perturbed_files directory')
parser.add_option("-i", "--input", dest="input_file", help="input unperturbed file")
parser.add_option("-p", "--perturbthk", dest="pThk", help="perturb ice thickness (1 m)")

for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]" 
options, args = parser.parse_args()

#---------- constants ----------------
rho_i = 910.0
rho_w = 1028.0
#-------------------------------------

perturb_thickness = float(options.pThk)

in_file = options.input_file

subprocess.call(["cp",in_file,"old.nc"]) 
subprocess.call(["ncks","-v","thickness","old.nc","tmp.nc"]) 
subprocess.call(["ncrename","-v","thickness,perturbMask","tmp.nc","perturbMask.nc"]) 
subprocess.call(["ncks","-A","-v","perturbMask","perturbMask.nc","old.nc"]) 
subprocess.call(["mkdir","perturbed_files"]) 

dataset = Dataset("old.nc",'r+')

xcoord = dataset.variables['xCell'][:]
ycoord = dataset.variables['yCell'][:]

thickness = dataset.variables['thickness'][0,:]
bedTopography = dataset.variables['bedTopography'][0,:]
perturbMask = dataset.variables['perturbMask'][0,:]
areaCell = dataset.variables['areaCell'][:]
minAreaCell = np.min(areaCell)

cellMask = bedTopography + thickness*rho_i/rho_w

perturbMask[:] = 0.0

nCells = len(dataset.dimensions['nCells'])

index1 = 0

index_array = np.zeros((10000,2))
area_array = np.zeros((10000,2))
coord_array = np.zeros((10000,2))
# we can also use append to avoid prescribe a large array before filling it out

for s in range(nCells):

    thicknessNew = np.copy(thickness)

    #if (areaCell[s]<minAreaCell*1.3) and ((cellMask[s]&4)==4) and ((cellMask[s]&2)==2):
    if (ycoord[s]>40e3) and (cellMask[s]<0) and (xcoord[s]<530e3):
        thicknessNew[s] = thickness[s] - perturb_thickness
        perturbMask[s] = 1
        index1 = index1 + 1
        index_array[index1-1,0] = index1
        index_array[index1-1,1] = s
        area_array[index1-1,0] = index1
        area_array[index1-1,1] = areaCell[s]
        coord_array[index1-1,0] = xcoord[s]
        coord_array[index1-1,1] = ycoord[s]
        out_file = "perturbed_files/landice_grid_"+str(index1)+".nc"
        subprocess.call(["cp","old.nc",out_file]) 
        dataset_tmp = Dataset(out_file,'r+')
        dataset_tmp['thickness'][0,:] = thicknessNew
        dataset_tmp.close()



dataset.variables['perturbMask'][0,:] = perturbMask

dataset.close()

np.savetxt('perturbed_files/index.txt', index_array[0:index1,:])
np.savetxt('perturbed_files/area.txt', area_array[0:index1,:])
np.savetxt('perturbed_files/coord.txt', coord_array[0:index1,:])

subprocess.call(["rm","old.nc"]) 

