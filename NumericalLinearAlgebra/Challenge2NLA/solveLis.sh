#!/bin/bash

LIS_DATA_DIR=../lis-2.1.6/test/Data/Challenge2/
LIS_BIN_DIR=../lis-2.1.6/test/

# TASK3
echo "***************************************"
echo "* TASK3"
echo "***************************************"
# Copy ./data/ATA.mtx into the lis data directory 
cp ./data/ATA.mtx $LIS_DATA_DIR

# Compute the largest eigenvalue of A^TA (tol: 1.0e-8) => Power method
mpirun -n 4 $LIS_BIN_DIR/eigen1 $LIS_DATA_DIR/ATA.mtx $LIS_DATA_DIR/ATAeigvec.txt $LIS_DATA_DIR/ATAhist.txt -e pi -etol 1.e-8

echo "***************************************"
echo "* TASK4"
echo "***************************************"

# TASK4
# Find a μ yielding an acceleration of the previous eigensolver (power method)
mpirun -n 4 $LIS_BIN_DIR/eigen1 $LIS_DATA_DIR/ATA.mtx $LIS_DATA_DIR/ATAeigvec.txt $LIS_DATA_DIR/ATAhist.txt -e pi -etol 1.e-8 -shift 700
