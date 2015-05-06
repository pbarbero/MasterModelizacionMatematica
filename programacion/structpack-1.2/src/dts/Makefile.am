## Process this file with automake to produce Makefile.in

## Flags used in every rule in the Makefile
##
AM_CPPFLAGS = -I$(top_srcdir)/utils -I$(top_builddir)/utils \
	-I$(top_builddir)/src/dst -I$(top_builddir)/src/diagc @use_mkl_flag@
AM_FCFLAGS = @openmp_fc_flags@
AM_CFLAGS = @openmp_cc_flags@

## Libraries
##
noinst_LTLIBRARIES = libdts.la

libdts_la_CPPFLAGS = $(AM_CPPFLAGS) @f95_cpp_input_flag@
libdts_la_SOURCES = sp_dts.F90 csp_dts.F90 

## Programs
##
bin_PROGRAMS = dts

dts_SOURCES = affinity2.c affinity.c cmdline_dts.c dts.c cmdline_dts.h csp_dts.h
## Trick for using Fortran linker
##
nodist_EXTRA_dts_SOURCES = dummy.F90
dts_FCFLAGS =
## Static linking with the libtool libraries
##
dts_LDFLAGS = -static @nofor_main_flag@ @openmp_fc_flags@
dts_LDADD = \
	libdts.la \
	$(top_builddir)/src/diagc/libdiagc.la \
	$(top_builddir)/src/dst/libdst.la \
	$(top_builddir)/utils/libutils.la \
	@dfftpack_lib_flags@ \
	@blas_lib_flags@ \
	@mkl_fc_lib_flags@ \
	@intel_extra_fc_lib_flags@

## Header files that must be installed
##
nodist_include_HEADERS = sp_dts.mod csp_dts.mod
dist_include_HEADERS = csp_dts.h

## Other files that must be installed
##
mexdir = $(exec_prefix)/mex
dist_mex_DATA = mex_dts.c dts.m

## Other files that must be distributed
##
dist_noinst_DATA = README cmdline_dts.ggo

## Testing
##
check_PROGRAMS =
dist_check_SCRIPTS = test_dts.sh
TESTS = test_dts.sh

MOSTLYCLEANFILES = *.mod
