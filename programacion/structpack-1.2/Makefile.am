## Process this file with automake to produce Makefile.in

## Subdirectories
##
SUBDIRS = utils src mex examples
DIST_SUBDIRS = $(SUBDIRS) doc

## Library versioning
## See the libtool manual for an explanation of the numbers
##
STRUCTPACK_CURRENT = 0
STRUCTPACK_REVISION = 0
STRUCTPACK_AGE = 0
STRUCTPACK_LT_VERSION = $(STRUCTPACK_CURRENT):$(STRUCTPACK_REVISION):$(STRUCTPACK_AGE)

## Flags used in every rule in the Makefile
##
AM_CPPFLAGS = -I$(top_srcdir)/utils -I$(top_builddir)/utils \
	-I$(top_builddir)/src/dctii \
	-I$(top_builddir)/src/dst \
	-I$(top_builddir)/src/rft \
	-I$(top_builddir)/src/conv \
	-I$(top_builddir)/src/diagc \
	-I$(top_builddir)/src/dts \
	-I$(top_builddir)/src/dpis \
	-I$(top_builddir)/src/dt \
	-I$(top_builddir)/src/dtspg \
	@use_mkl_flag@
AM_FCFLAGS = @openmp_fc_flags@
AM_CFLAGS = @openmp_cc_flags@

## Libraries
##
SUBLIBS = utils/libutils.la \
	src/dctii/libdctii.la \
        src/dst/libdst.la \
        src/rft/librft.la \
        src/conv/libconv.la \
        src/diagc/libdiagc.la \
	src/dts/libdts.la \
	src/dpis/libdpis.la \
	src/dt/libdt.la \
	src/dtspg/libdtspg.la

lib_LTLIBRARIES = libstructpack.la

libstructpack_la_SOURCES = structpack.F90
libstructpack_la_LDFLAGS = -version-info $(STRUCTPACK_LT_VERSION)
libstructpack_la_LIBADD = $(SUBLIBS)

include_HEADERS = structpack.h
nodist_include_HEADERS = structpack.mod

EXTRA_DIST =

MOSTLYCLEANFILES = *.mod
