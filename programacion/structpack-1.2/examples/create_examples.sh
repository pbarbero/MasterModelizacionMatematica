#! /bin/sh
# create_examples.sh
# Build the example programs provided.
#
# Date: 7 November 2011
# Author: Pablo Martinez Naredo <pmnaredo@gmail.com>

prefix="/usr/local"
exec_prefix="${prefix}"

STRUCTPACK_INCLUDE_DIR="${prefix}/include"
STRUCTPACK_LIB_DIR="${exec_prefix}/lib"

DFFTPACK_LIB_FLAGS=""

BLAS_LIB_FLAGS="-lblas"
LAPACK_LIB_FLAGS="-llapack"

MKL_LIB_DIR_FLAG=""
MKL_CC_LIB_FLAGS=""
MKL_FC_LIB_FLAGS=""

EXTRA_CC_LIB_FLAGS=""
EXTRA_FC_LIB_FLAGS=""

CC="gcc"
CFLAGS="-O -Wall -std=c89 -DschedYes"

FC="gfortran"
FCFLAGS="-O -Wall -std=f2003 -fall-intrinsics -ffree-line-length-none -DschedYes"

OPENMP_CC_FLAGS="-fopenmp"
OPENMP_FC_FLAGS="-fopenmp"


# dts_example_c
#
$CC -I$STRUCTPACK_INCLUDE_DIR $CFLAGS -c -o dts_example_c.o \
    dts_example.c

$FC  -o dts_example_c dts_example_c.o \
    $STRUCTPACK_LIB_DIR/libstructpack.a \
    $DFFTPACK_LIB_FLAGS \
    $BLAS_LIB_FLAGS \
    $MKL_LIB_DIR_FLAG $MKL_FC_LIB_FLAGS \
    $EXTRA_FC_LIB_FLAGS \
    $OPENMP_FC_FLAGS


# dts_example_fortran
#
$FC -I$STRUCTPACK_INCLUDE_DIR $FCFLAGS -c -o dts_example_fortran.o \
    dts_example.F90

$FC -o dts_example_fortran dts_example_fortran.o \
    $STRUCTPACK_LIB_DIR/libstructpack.a \
    $DFFTPACK_LIB_FLAGS \
    $BLAS_LIB_FLAGS \
    $MKL_LIB_DIR_FLAG $MKL_FC_LIB_FLAGS \
    $EXTRA_FC_LIB_FLAGS \
    $OPENMP_FC_FLAGS


# dpis_example_c
#
$CC -I$STRUCTPACK_INCLUDE_DIR $CFLAGS -c -o dpis_example_c.o \
    dpis_example.c

$FC  -o dpis_example_c dpis_example_c.o \
    $STRUCTPACK_LIB_DIR/libstructpack.a \
    $DFFTPACK_LIB_FLAGS \
    $BLAS_LIB_FLAGS \
    $MKL_LIB_DIR_FLAG $MKL_FC_LIB_FLAGS \
    $EXTRA_FC_LIB_FLAGS \
    $OPENMP_FC_FLAGS


# dpis_example_fortran
#
$FC -I$STRUCTPACK_INCLUDE_DIR $FCFLAGS -c -o dpis_example_fortran.o \
    dpis_example.F90

$FC -o dpis_example_fortran dpis_example_fortran.o \
    $STRUCTPACK_LIB_DIR/libstructpack.a \
    $DFFTPACK_LIB_FLAGS \
    $BLAS_LIB_FLAGS \
    $MKL_LIB_DIR_FLAG $MKL_FC_LIB_FLAGS \
    $EXTRA_FC_LIB_FLAGS \
    $OPENMP_FC_FLAGS
    
    
    
    

# dt_example_c
#
$CC -I$STRUCTPACK_INCLUDE_DIR $CFLAGS -c -o dt_example_c.o \
    dt_example.c

$FC  -o dt_example_c dt_example_c.o \
    $STRUCTPACK_LIB_DIR/libstructpack.a \
    $DFFTPACK_LIB_FLAGS \
    $BLAS_LIB_FLAGS \
    $MKL_LIB_DIR_FLAG $MKL_FC_LIB_FLAGS \
    $EXTRA_FC_LIB_FLAGS \
    $OPENMP_FC_FLAGS


# dt_example_fortran
#
$FC -I$STRUCTPACK_INCLUDE_DIR $FCFLAGS -c -o dt_example_fortran.o \
    dt_example.F90

$FC -o dt_example_fortran dt_example_fortran.o \
    $STRUCTPACK_LIB_DIR/libstructpack.a \
    $DFFTPACK_LIB_FLAGS \
    $BLAS_LIB_FLAGS \
    $MKL_LIB_DIR_FLAG $MKL_FC_LIB_FLAGS \
    $EXTRA_FC_LIB_FLAGS \
    $OPENMP_FC_FLAGS


# dtspg_example_c
#
$CC -I$STRUCTPACK_INCLUDE_DIR $CFLAGS -c -o dtspg_example_c.o \
    dtspg_example.c

$FC  -o dtspg_example_c dtspg_example_c.o \
    $STRUCTPACK_LIB_DIR/libstructpack.a \
    $DFFTPACK_LIB_FLAGS \
    $BLAS_LIB_FLAGS $LAPACK_LIB_FLAGS \
    $MKL_LIB_DIR_FLAG $MKL_FC_LIB_FLAGS \
    $EXTRA_FC_LIB_FLAGS \
    $OPENMP_FC_FLAGS


# dtspg_example_fortran
#
$FC -I$STRUCTPACK_INCLUDE_DIR $FCFLAGS -c -o dtspg_example_fortran.o \
    dtspg_example.F90

$FC -o dtspg_example_fortran dtspg_example_fortran.o \
    $STRUCTPACK_LIB_DIR/libstructpack.a \
    $DFFTPACK_LIB_FLAGS \
    $BLAS_LIB_FLAGS $LAPACK_LIB_FLAGS \
    $MKL_LIB_DIR_FLAG $MKL_FC_LIB_FLAGS \
    $EXTRA_FC_LIB_FLAGS \
    $OPENMP_FC_FLAGS
    
