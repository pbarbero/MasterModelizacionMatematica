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
noinst_LTLIBRARIES = libdtspg.la

libdtspg_la_CPPFLAGS = $(AM_CPPFLAGS) @f95_cpp_input_flag@
libdtspg_la_SOURCES = sp_dtspg_.F90 csp_dtspg.F90

## Programs
##
bin_PROGRAMS = dtspg

dtspg_SOURCES = cmdline_dtspg.c dtspg.c cmdline_dtspg.h csp_dtspg.h
## Trick for using Fortran linker
##
nodist_EXTRA_dtspg_SOURCES = dummy.F90 
dtspg_FCFLAGS =
#dtspg_LDFLAGS = @nofor_main_flag@ @openmp_fc_flags@
# Static linking with the libtool libraries
dtspg_LDFLAGS = -static @nofor_main_flag@ @openmp_fc_flags@
dtspg_LDADD = @lapack_lib_flags@ \
	libdtspg.la \
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
nodist_include_HEADERS = sp_dtspg.mod csp_dtspg.mod
dist_include_HEADERS = csp_dtspg.h

## Other files that must be installed
##
mexdir = $(exec_prefix)/mex
dist_mex_DATA = mex_dtspg.c dtspg.m

## Other files that must be distributed
##
dist_noinst_DATA = README cmdline_dtspg.ggo

## Testing
##
check_PROGRAMS =
dist_check_SCRIPTS = test_dtspg.sh
TESTS = test_dtspg.sh

MOSTLYCLEANFILES = *.mod
