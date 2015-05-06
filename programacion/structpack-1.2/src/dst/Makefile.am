## Process this file with automake to produce Makefile.in

## Flags used in every rule in the Makefile
##
AM_CPPFLAGS = -I$(top_builddir)/utils @f95_cpp_input_flag@ @use_mkl_flag@
AM_FCFLAGS = @openmp_fc_flags@

noinst_LTLIBRARIES = libdst.la

libdst_la_SOURCES = sp_dst.F90

nodist_include_HEADERS = sp_dst.mod

MOSTLYCLEANFILES = *.mod
