import numpy as np
from netCDF4 import Dataset
import subprocess

backgroundBmbFile = "backgroundBmb.nc"
# Rignot background bmb (converved from positive to negative using ncap2)
# we can get it by doing ncks -v floatingBasalMassBal mpas_initial_mesh_file.nc backgroundBmb.nc
bmbPerturb = 8.0
# Additional basal melt rate (m/a)
larmipRegionMaskFile = "larmipRegionMasks.nc"
# LARMIP region mask data, converted from LARMIP nc file to MPAS nc file using script interpolate_to_mpasli_grid.py (see interpolate_mask_to_mpas.sh)
regionMaskKeyWord = "larmipRegionMasks"
# the mask data name in the larmipRegionMaskFile

backgroundBmbData = Dataset(backgroundBmbFile,'r')

larmipRegionMaskData = Dataset(larmipRegionMaskFile,'r')
larmipRegionMask = larmipRegionMaskData.variables[regionMaskKeyWord][0,:]

nCells = len(backgroundBmbData.dimensions['nCells'])
nCells1 = len(larmipRegionMaskData.dimensions['nCells'])
assert nCells == nCells1

larmipRegionDic = {1:"EAIS", 2:"Ross", 3:"Amundsen", 4:"Weddell", 5:"Peninsula"}

for i in range(5):

    regionName = larmipRegionDic[i+1]
    print "processing region"+regionName

    newFileName = regionName + "Bmb.nc"

    subprocess.call(["cp",backgroundBmbFile,newFileName])

    newBmbdata = Dataset(newFileName,'r+')
    newBmb = newBmbdata.variables['floatingBasalMassBal'][0,:]

    for j in range(nCells):

        maskVal = larmipRegionMask[j]

        if abs((maskVal - (i+1))) < 1e-10:
            newBmb[j] = newBmb[j] - bmbPerturb*910/(365.0*24.0*3600.0)
            newBmbdata.variables['floatingBasalMassBal'][0,j] = newBmb[j]

    newBmbdata.close()

backgroundBmbData.close()
larmipRegionMaskData.close()
