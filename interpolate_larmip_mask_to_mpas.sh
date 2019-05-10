#!/bin/bash

inFile="./ais_2_to_20_km.nc"
#inFile="./ais20km.20180126.sfcMassBalFix.floatingBmbAdd.nc"
larmipRegionFile="LARMIP_regions_initMIPgrid_04.nc"

cp $inFile inFileCopy.nc

#ncrename -v x,x1 $larmipRegionFile
#ncrename -v y,y1 $larmipRegionFile
#ncrename -v regions,thk $larmipRegionFile

~/Apps/MPAS-Tools/grid_gen/landice_grid_tools/interpolate_to_mpasli_grid.py -d inFileCopy.nc -s $larmipRegionFile -m d

ncrename -v thickness,larmipRegionMasks inFileCopy.nc
ncks -O -v larmipRegionMasks inFileCopy.nc larmipRegionMasks.nc
