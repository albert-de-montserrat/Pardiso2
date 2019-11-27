#!/bin/bash
#export MKL_INC=$MKLROOT/include
#export MKL_LIB=$MKLROOT/lib

icc -I$HDF5_INC -L$HDF5_LIB -lhdf5  -I$MKL_INC -L$MKL_LIB -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -mcmodel=medium -qopenmp -DH5_NO_DEPRECATED_SYMBOLS -o in3mg  in3mg.c
icc -I$HDF5_INC -L$HDF5_LIB -lhdf5  -I$MKL_INC -L$MKL_LIB -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -mcmodel=medium -qopenmp -DH5_NO_DEPRECATED_SYMBOLS -o i3mg   i3mg.c
icc -I$HDF5_INC -L$HDF5_LIB -lhdf5  -I$MKL_INC -L$MKL_LIB -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -mcmodel=medium -qopenmp -DH5_NO_DEPRECATED_SYMBOLS -o viz    viz.c
