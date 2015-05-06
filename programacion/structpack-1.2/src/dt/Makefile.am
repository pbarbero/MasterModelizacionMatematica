## Process this file with automake to produce Makefile.in

## Flags used in every rule in the Makefile
##
AM_CPPFLAGS = -I$(top_srcdir)/utils -I$(top_builddir)/utils \
	-I$(top_builddir)/src/dst -I$(top_builddir)/src/diagc \
	-I$(top_builddir)/src/dctii -I$(top_builddir)/src/conv \
	-I$(top_builddir)/src/rft \
	@use_mkl_flag@
AM_FCFLAGS = @openmp_fc_flags@
AM_CFLAGS = @openmp_cc_flags@

## Libraries
##
noinst_LTLIBRARIES = libdt.la

libdt_la_CPPFLAGS = $(AM_CPPFLAGS) @f95_cpp_input_flag@
libdt_la_SOURCES = sp_dt.F90 csp_dt.F90

## Programs
##
bin_PROGRAMS = dt

dt_SOURCES = cmdline_dt.c dt.c cmdline_dt.h csp_dt.h
## Trick for using Fortran linker
##
nodist_EXTRA_dt_SOURCES = dummy.F90 
dt_FCFLAGS =
#dt_LDFLAGS = @nofor_main_flag@ @openmp_fc_flags@
# Static linking with the libtool libraries
dt_LDFLAGS = -static @nofor_main_flag@ @openmp_fc_flags@
dt_LDADD =  @lapack_lib_flags@ \
	libdt.la \
	$(top_builddir)/src/diagc/libdiagc.la \
	$(top_builddir)/src/dst/libdst.la \
	$(top_builddir)/src/dctii/libdctii.la \
	$(top_builddir)/utils/libutils.la \
	$(top_builddir)/src/conv/libconv.la \
	$(top_builddir)/src/rft/librft.la \
	@dfftpack_lib_flags@ \
	@blas_lib_flags@ \
	@mkl_fc_lib_flags@ \
	@intel_extra_fc_lib_flags@

## Header files that must be installed
##
nodist_include_HEADERS = sp_dt.mod csp_dt.mod
dist_include_HEADERS = csp_dt.h

## Other files that must be installed
##
mexdir = $(exec_prefix)/mex
dist_mex_DATA = mex_dt.c dt.m

## Other files that must be distributed
##
dist_noinst_DATA = README cmdline_dt.ggo

## Testing
##
check_PROGRAMS =
dist_check_SCRIPTS = test_dt.sh
TESTS = test_dt.sh

MOSTLYCLEANFILES = *.mod
