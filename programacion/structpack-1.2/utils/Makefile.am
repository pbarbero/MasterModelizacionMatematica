## Process this file with automake to produce Makefile.in

## Subdirectories
##
if BUILD_DFFTPACK
DFFTPACK_DIR = dfftpack
DFFTPACK_LIB = $(DFFTPACK_DIR)/libdfftpack.la
endif
SUBDIRS = $(DFFTPACK_DIR)

## Flags used in every rule in the Makefile
##
AM_CFLAGS = @openmp_cc_flags@
AM_FCFLAGS = @openmp_fc_flags@

## Libraries
##
noinst_LTLIBRARIES = libutils.la

if USE_MKL
dist_libutils_la_SOURCES = sp_utils.F90 sp_ctimer.F90 cutils.c ctimer.c sp_toegen.F90 csp_toegen.F90 
nodist_libutils_la_SOURCES = mkl_dfti.f90
else
libutils_la_SOURCES = sp_utils.F90 sp_ctimer.F90 cutils.c ctimer.c sp_toegen.F90 csp_toegen.F90 
endif
libutils_la_LIBADD = $(DFFTPACK_LIB)

## Header files
##
nodist_include_HEADERS = sp_utils.mod sp_ctimer.mod sp_toegen.mod csp_toegen.mod 
dist_noinst_HEADERS = blas.h cutils.h csp_toegen.h 

MOSTLYCLEANFILES = *.mod
CLEANFILES = mkl_dfti.f90

mkl_dfti.f90:
	$(LN_S) @mkl_include_dir@/mkl_dfti.f90
