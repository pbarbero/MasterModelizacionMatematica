## Process this file with automake to produce Makefile.in

## Flags used in every rule in the Makefile
##
AM_CPPFLAGS = -I$(top_builddir)/utils -I$(top_builddir)/src/rft @f95_cpp_input_flag@ @use_mkl_flag@
AM_FCFLAGS = @openmp_fc_flags@

noinst_LTLIBRARIES = libconv.la

libconv_la_SOURCES = sp_conv.F90 csp_conv.F90

nodist_include_HEADERS = sp_conv.mod csp_conv.mod
dist_include_HEADERS = csp_conv.h

MOSTLYCLEANFILES = *.mod
