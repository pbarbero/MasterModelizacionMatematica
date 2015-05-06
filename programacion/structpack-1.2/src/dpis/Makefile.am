## Process this file with automake to produce Makefile.in

## Flags used in every rule in the Makefile
##
AM_CPPFLAGS = -I$(top_srcdir)/utils -I$(top_builddir)/utils \
	-I$(top_builddir)/src/dst @use_mkl_flag@
AM_FCFLAGS = @openmp_fc_flags@
AM_CFLAGS = @openmp_cc_flags@

## Libraries
##
noinst_LTLIBRARIES = libdpis.la

libdpis_la_CPPFLAGS = $(AM_CPPFLAGS) @f95_cpp_input_flag@
libdpis_la_SOURCES = sp_dpis.F90 csp_dpis.F90

## Programs
##
bin_PROGRAMS = dpis

dpis_SOURCES = cmdline_dpis.c dpis.c cmdline_dpis.h \
	csp_dpis.h
## Trick for using Fortran linker
##
nodist_EXTRA_dpis_SOURCES = dummy.F90
dpis_FCFLAGS =
## Static linking with the libtool libraries
##
dpis_LDFLAGS = -static @nofor_main_flag@ @openmp_fc_flags@
dpis_LDADD = \
	libdpis.la \
        $(top_builddir)/src/dst/libdst.la \
	$(top_builddir)/utils/libutils.la \
	@dfftpack_lib_flags@ \
	@blas_lib_flags@ \
	@mkl_fc_lib_flags@ \
	@intel_extra_fc_lib_flags@

## Header files that must be installed
##
nodist_include_HEADERS = sp_dpis.mod csp_dpis.mod
dist_include_HEADERS = csp_dpis.h

## Other files that must be installed
##
mexdir = $(exec_prefix)/mex
dist_mex_DATA = mex_dpis.c dpis.m

## Other files that must be distributed
##
dist_noinst_DATA = README cmdline_dpis.ggo

## Testing
##
check_PROGRAMS =
dist_check_SCRIPTS = test_dpis.sh
TESTS = test_dpis.sh

MOSTLYCLEANFILES = *.mod
