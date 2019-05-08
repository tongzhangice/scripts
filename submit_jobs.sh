#!/bin/bash

for i in {1..738..1}
do

    #cd /net/scratch1/tzhang/projects/buttressing_test/highres/larsenC_originalCells/test_$i
    #cd /net/scratch4/tzhang/projects/buttressing_test/mismip+/test_$i
    cd /scratch2/scratchdirs/tongz/projects/buttressing_test/mismip+_2km_newExe/test_$i
    sbatch batch-script 

done

