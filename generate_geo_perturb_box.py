import numpy as np
from netCDF4 import Dataset
import subprocess
from optparse import OptionParser

def calculate_tri_area(x1,y1,x2,y2,x3,y3):
    return 0.5*abs((x1-x3)*(y2-y1)-(x1-x2)*(y3-y1))
# calculate the area of an triangle

def checkIfInRectangle(xi,yi,x1,y1,x2,y2):
    area1 = calculate_tri_area(xi,yi,x1,y1,x2,y1)
    area2 = calculate_tri_area(xi,yi,x1,y1,x1,y2)
    area3 = calculate_tri_area(xi,yi,x1,y2,x2,y2)
    area4 = calculate_tri_area(xi,yi,x2,y2,x2,y1)
    
    areaSum = area1 + area2 + area3 + area4

    area = abs((x2-x1)*(y2-y1))

    #print areaSum,area

    if abs(areaSum-area)<1e-3:
        return True
    else:
        return False
# check if a point (xi, yi) is in the box [(x1,y1),(x1,y2),(x2,y1),(x2,y2)], i.e, if the total area of four triangles equals to the area of the box

parser = OptionParser(epilog='create thickness perturbed geometries. The thickness is perturbed box by box (which contains a number of cells) for a specified region. The script will generate many perturbed files (landice_grid_xxx.nc) and three geometry text files (index.txt, area.txt and coord.txt) in the perturbed_files directory')
parser.add_option("-i", "--input", dest="input_file", help="input unperturbed file")
parser.add_option("-p", "--perturbthk", dest="pThk", help="perturb ice thickness (1 m)")

for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]" 
options, args = parser.parse_args()

#---------- constants ----------------
rho_i = 910.0
rho_o = 1028.0
dx = 20.0e3
# box size in x
dy = 20.0e3
# box size in y
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

xmin = min(xcoord)
xmax = max(xcoord)
ymin = min(ycoord)
ymax = max(ycoord)

Nx = int(np.ceil((xmax-xmin)/dx))
Ny = int(np.ceil((ymax-ymin)/dy))

mask = bedTopography + rho_i/rho_o*thickness

perturbMask[:] = 0.0

nCells = len(dataset.dimensions['nCells'])

index1 = 0

index_array = np.zeros((10000,2))
area_array = np.zeros((10000,2))
coord_array = np.zeros((10000,2))
# we can also use append to avoid prescribe a large array before filling it out

for i in range(Nx-1):
    for j in range(Ny-1):

        thicknessNew = np.copy(thickness)

        if_changed = 0
        area = 0

        x1 = xmin + i*dx
        x2 = xmin + (i+1)*dx
        y1 = ymin + j*dy
        y2 = ymin + (j+1)*dy

        index = i*(Ny-1)+j 
        x_mean = (x1+x2)/2.0
        y_mean = (y1+y2)/2.0

        for s in range(nCells):

            x = xcoord[s]
            y = ycoord[s]

            if checkIfInRectangle(x,y,x1,y1,x2,y2) and (mask[s]<=0) and (thickness[s]>10):
                # 10 is the dynamic thickness threshold
                thicknessNew[s] = thickness[s] - 1
                perturbMask[s] = 1
                if_changed = 1
                area = area + areaCell[s]

        if if_changed == 1:

            index1 = index1 + 1
            index_array[index1-1,0] = index1
            index_array[index1-1,1] = index
            area_array[index1-1,0] = index1
            area_array[index1-1,1] = area
            coord_array[index1-1,0] = x_mean
            coord_array[index1-1,1] = y_mean
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
