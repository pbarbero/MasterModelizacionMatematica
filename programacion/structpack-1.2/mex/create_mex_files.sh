#! /bin/sh
# create_mex_files.sh
# Create the MEX files suitable for using the StructPack programs in Matlab
# or Octave
#
# Date: 28 March 2011
# Author: Pablo Martinez Naredo <pmnaredo@gmail.com>

prefix="/usr/local"
exec_prefix="${prefix}"

STRUCTPACK_INCLUDE_DIR="${prefix}/include"
STRUCTPACK_LIB_DIR="${exec_prefix}/lib"
STRUCTPACK_LIB_FLAGS="${STRUCTPACK_LIB_DIR}/libstructpack.a"

DFFTPACK_LIB_FLAGS=""

BLAS_LIB_FLAGS="-lblas"

MKL_LIB_DIR_FLAG=""
MKL_LIB_FLAGS=""

#OPENMP_FLAGS="-Wl,-fopenmp -lgomp"

# Flags for Octave
#
MEX="mkoctfile"
MEXFLAGS="--mex -v"

# Flags for Matlab
#
#MEX=mex
#MEXFLAGS=-v #-f $(MATLABROOT)/matopts.sh

export CC="gcc"
export DL_LD="gfortran"

case $DL_LD in
ifort)
  OPENMP_LIB_FLAGS="-liomp5"
  ;;
*)
  OPENMP_LIB_FLAGS="-lgomp"
  ;;
esac

EXTRA_LIB_FLAGS=""

# dts
#
$MEX $MEXFLAGS -I$STRUCTPACK_INCLUDE_DIR mex_dts.c \
    $STRUCTPACK_LIB_FLAGS \
    $DFFTPACK_LIB_FLAGS \
    $BLAS_LIB_FLAGS \
    $MKL_LIB_DIR_FLAG $MKL_LIB_FLAGS \
    $OPENMP_LIB_FLAGS \
    $EXTRA_LIB_FLAGS


# dpis
#
$MEX $MEXFLAGS -I$STRUCTPACK_INCLUDE_DIR mex_dpis.c \
    $STRUCTPACK_LIB_FLAGS \
    $DFFTPACK_LIB_FLAGS \
    $BLAS_LIB_FLAGS \
    $MKL_LIB_DIR_FLAG $MKL_LIB_FLAGS \
    $OPENMP_LIB_FLAGS \
    $EXTRA_LIB_FLAGS

# dt
#
$MEX $MEXFLAGS -I$STRUCTPACK_INCLUDE_DIR mex_dt.c \
    $STRUCTPACK_LIB_FLAGS \
    $DFFTPACK_LIB_FLAGS \
    $BLAS_LIB_FLAGS \
    $MKL_LIB_DIR_FLAG $MKL_LIB_FLAGS \
    $OPENMP_LIB_FLAGS \
    $EXTRA_LIB_FLAGS

# dtspg
#
$MEX $MEXFLAGS -I$STRUCTPACK_INCLUDE_DIR mex_dtspg.c \
    $STRUCTPACK_LIB_FLAGS \
    $DFFTPACK_LIB_FLAGS \
    $BLAS_LIB_FLAGS \
    $MKL_LIB_DIR_FLAG $MKL_LIB_FLAGS \
    $OPENMP_LIB_FLAGS \
    $EXTRA_LIB_FLAGS
