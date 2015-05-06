#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#

AC_PREREQ(2.59)
# Ranilla 26-03-2013
AC_INIT([StructPack], [1.2], [structpack@gmail.com])
# Fin Ranilla 26-03-2013

AC_CONFIG_SRCDIR([src/dst/sp_dst.F90])

AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

# Default values
#
default_gcc_cflags="-O -Wall -std=c89"
default_gfortran_fcflags="-O -Wall -std=f2003 -fall-intrinsics -ffree-line-length-none"
default_icc_cflags="-O -std=c89"
default_ifort_fcflags="-O -auto -cpp -std03"

# Command line options
#
# --with-ifort
#
AC_ARG_WITH([ifort],
  [AS_HELP_STRING([--with-ifort],
    [use ifort instead of gfortran @<:@default=no@:>@])],
  [case $withval in
  yes | no)
    ;;
  *)
    AC_MSG_ERROR(['$withval' is not a valid argument for --with-ifort option])
    ;;
  esac],
  [with_ifort=no])

# --with-icc
#
AC_ARG_WITH([icc],
  [AS_HELP_STRING([--with-icc],
    [use icc instead of gcc @<:@default=no@:>@])],
  [case $withval in
  yes | no)
    ;;
  *)
    AC_MSG_ERROR(['$withval' is not a valid argument for --with-icc option])
    ;;
  esac],
  [with_icc=no])

# Ranilla 26-03-2013
# --with-sched
AC_ARG_WITH([sched],
  [AS_HELP_STRING([--with-sched],
    [use sched instead default settings @<:@default=yes@:>@])],
  [case $withval in
  yes | no)
    ;;
  *)
    AC_MSG_ERROR(['$withval' is not a valid argument for --with-sched option])
    ;;
  esac],
  [with_sched=yes])
# Fin Ranilla 26-03-2013

# --with-blas
#
AC_ARG_WITH([blas],
  [AS_HELP_STRING([--with-blas@<:@=FLAGS@:>@],
    [use generic BLAS and use FLAGS for linking with the BLAS library if given
      @<:@default=yes@:>@])],
  [with_blas_given=yes],
  [with_blas=yes])

# --with-lapack
#
AC_ARG_WITH([lapack],
  [AS_HELP_STRING([--with-lapack@<:@=FLAGS@:>@],
    [use generic LAPACK (requires BLAS) and use FLAGS for linking with the
      LAPACK library if given @<:@default=yes@:>@])],
  [with_lapack_given=yes],
  [with_lapack=yes])

# --with-dfftpack
#
AC_ARG_WITH([dfftpack],
  [AS_HELP_STRING([--with-dfftpack@<:@=FLAGS@:>@],
    [use dfftpack and use FLAGS for linking with the dfftpack library if given
      @<:@default=yes@:>@])],
  [],
  [with_dfftpack=yes])

# --with-mkl
#
AC_ARG_WITH([mkl],
  [AS_HELP_STRING([--with-mkl],
    [use MKL @<:@default=no@:>@])],
  [case $withval in
  yes | no)
    ;;
  *)
    AC_MSG_ERROR(['$withval' is not a valid argument for --with-mkl option])
    ;;
  esac],
  [with_mkl=no])

# Precious variables for MKL
#
AC_ARG_VAR([MKLINCLUDE], [MKL header files directory])
AC_ARG_VAR([MKLPATH], [MKL library files directory])


# Check inconsistencies between options
#
# Check whether MKL is used
#
if test "x$with_mkl" != xno; then
  # MKL BLAS is used
  #
  test "x$with_blas" != xno && test "x$with_blas_given" = xyes &&
    AC_MSG_WARN([--with-blas option will be ignored and MKL BLAS will be used])
  with_blas=no
  # MKL LAPACK is used
  #
  test "x$with_lapack" != xno && test "x$with_lapack_given" = xyes &&
    AC_MSG_WARN([--with-lapack option will be ignored and MKL LAPACK will be
      used])
  with_lapack=no
fi

# Check for which system
#
AC_CANONICAL_HOST

# Checks for programs.
#
# grep with extended regular expressions
#
AC_PROG_EGREP

# sed
#
m4_ifdef([AC_PROG_SED],
  [AC_PROG_SED],
  [AC_CHECK_PROG([SED], [sed], [sed], [:])])

# Fortran compiler
#
if test "x$with_ifort" = xyes; then
  fc="ifort"
else
  fc="gfortran"
fi
: ${FCFLAGS=""}
AC_PROG_FC([$fc])

openmp_fc_flags=
nofor_main_flag=
f95_cpp_input_flag=
case $FC in
gfortran*)
  test "x$FCFLAGS" = x && FCFLAGS="$default_gfortran_fcflags" 
  test "x$openmp_fc_flags" = x && openmp_fc_flags="-fopenmp"
  test "x$f95_cpp_input_flag" = x && f95_cpp_input_flag="-x f95-cpp-input"
  ;;
ifort)
  test "x$FCFLAGS" = x && FCFLAGS="$default_ifort_fcflags"
  test "x$openmp_fc_flags" = x && openmp_fc_flags="-openmp"
  test "x$nofor_main_flag" = x && nofor_main_flag="-nofor-main"
  ;;
*)
  ;;
esac
AC_SUBST([openmp_fc_flags])
AC_SUBST([nofor_main_flag])
AC_SUBST([f95_cpp_input_flag])

# Check for gfortran version
#
gfortran_version=
gfortran_version_major=
gfortran_version_minor=
case $FC in
gfortran*)
  AC_MSG_CHECKING([for $FC version])
  gfortran_version=$($FC --version |
    $EGREP -o ['[[:digit:]]+[.][[:digit:]]+([.][[:digit:]]+)?([[:space:]]|$)'])
  # Only first version number found is considered
  #
  for gfortran_version in $gfortran_version; do
    break
  done
  AC_MSG_RESULT([$gfortran_version])

  if test "x$gfortran_version" != x; then
    # Get major and minor components of the gfortran version number
    # Follows AX_SPLIT_VERSION macro from AC-Archive
    # Tom Howard <tomhoward@users.sf.net>
    #
    gfortran_version_major=$(echo "$gfortran_version" |
      $SED 's/\([[^.]][[^.]]*\).*/\1/')
    gfortran_version_minor=$(echo "$gfortran_version" |
      $SED 's/[[^.]][[^.]]*.\([[^.]][[^.]]*\).*/\1/')

    # Check major and minor components
    #
    test "$gfortran_version_major" -lt 4 || {
      test "$gfortran_version_major" -eq 4 &&
      test "$gfortran_version_minor" -lt 4; } &&
      AC_MSG_ERROR([gfortran version 4.4.0 or higher is required])
  else
    AC_MSG_WARN([cannot check $FC version])
  fi
  ;;
esac

# Fortran77 compiler
#
test "x$F77" = x && F77="$FC"
test "x$FFLAGS" = x && FFLAGS="$FCFLAGS"
AC_PROG_F77

# Determine the linker flags for the Fortran intrinsic and runtime libraries
#
AC_FC_LIBRARY_LDFLAGS

# C compiler
#
if test "x$with_icc" = xyes; then
  cc="icc"
else
  cc="gcc"
fi
: ${CFLAGS=""}
AC_PROG_CC([$cc])

openmp_cc_flags=
case $CC in
gcc*)
  test "x$CFLAGS" = x && CFLAGS="$default_gcc_cflags"
  test "x$openmp_cc_flags" = x && openmp_cc_flags="-fopenmp"
  ;;
icc)
  test "x$CFLAGS" = x && CFLAGS="$default_icc_cflags"
  test "x$openmp_cc_flags" = x && openmp_cc_flags="-openmp"
  ;;
*)
  ;;
esac
AC_SUBST([openmp_cc_flags])

# Ranilla y Pedro 28-03-2013
if test "x$with_sched" = xyes; then
  # Check for the presence of file sched.h
  # AC_CHECK_FILE([$KERNEL_INC/sched.h], [have_sched="yes"], [have_sched="no"])
  AC_CHECK_HEADERS([sched.h], [sched_flags="-DschedYes"], [], [])
  FCFLAGS="$FCFLAGS $sched_flags"
  FFLAGS="$FFLAGS $sched_flags"
  CFLAGS="$CFLAGS $sched_flags"
  #if test "x$have_sched" = no; then
  #  AC_MSG_ERROR([cannot find 'sched.h'])
  #else
  #  sched_flags="-DschedYes"
  #  FCFLAGS="$FCFLAGS $sched_flags"
  #  FFLAGS="$FFLAGS $sched_flags"
  #  CFLAGS="$CFLAGS $sched_flags"
  #fi
fi
# Fin Ranilla y Pedro 28-03-2013

# Libtool
#
AC_PROG_LIBTOOL

# Other program checks
#
AC_PROG_LN_S

# Libraries
#
# dfftpack
#
build_dfftpack=
dfftpack_lib_flags=
AC_FC_FUNC([dsint])
case $with_dfftpack in
no)
  ;;
yes)
  build_dfftpack="yes"
  ;;
*)
  dfftpack_lib_flags="$with_dfftpack"

  # Get Fortran linker name of function to check for
  #
  AC_FC_FUNC([dsint])

  save_LIBS="$LIBS"
  LIBS="$dfftpack_lib_flags $LIBS $FCLIBS"

  # Try linking using the specified flags
  #
  AC_MSG_CHECKING([for $dsint in $dfftpack_lib_flags])
  AC_LINK_IFELSE([AC_LANG_CALL([], [$dsint])],
    [lib_dfftpack_dsint="yes"],
    [lib_dfftpack_dsint="no"])
  AC_MSG_RESULT([$lib_dfftpack_dsint])

  LIBS="$save_LIBS"

  test "x$lib_dfftpack_dsint" = xno &&
    AC_MSG_ERROR([cannot link against '$dfftpack_lib_flags'])
  ;;
esac
AC_SUBST([dfftpack_lib_flags])
AM_CONDITIONAL([BUILD_DFFTPACK], [test "x$build_dfftpack" = xyes])


# BLAS
#
blas_lib_flags=
case $with_blas in
no)
  ;;
*)
  if test "x$with_blas" = xyes; then
    # Default flag
    #
    test "x$blas_lib_flags" = x && blas_lib_flags="-lblas"
  else
    test "x$blas_lib_flags" = x && blas_lib_flags="$with_blas"
  fi

  # Get Fortran linker name of function to check for
  #
  AC_FC_FUNC([daxpy])

  save_LIBS="$LIBS"
  LIBS="$blas_lib_flags $LIBS $FCLIBS"

  # Try linking with BLAS libraries
  #
  AC_MSG_CHECKING([for $daxpy in $blas_lib_flags])
  AC_LINK_IFELSE([AC_LANG_CALL([], [$daxpy])],
    [lib_blas_daxpy="yes"],
    [lib_blas_daxpy="no"])
  AC_MSG_RESULT([$lib_blas_daxpy])

  LIBS="$save_LIBS"

  test "x$lib_blas_daxpy" = xno &&
    AC_MSG_ERROR([cannot link against '$blas_lib_flags'])
  ;;
esac
AC_SUBST([blas_lib_flags])


# LAPACK + BLAS
#
lapack_lib_flags=
case $with_lapack in
no)
  ;;
*)
  if test "x$with_lapack" = xyes; then
    # Default flag
    #
    test "x$lapack_lib_flags" = x && lapack_lib_flags="-llapack"
  else
    test "x$lapack_lib_flags" = x && lapack_lib_flags="$with_lapack"
  fi

  # Get Fortran linker name of function to check for
  #
  AC_FC_FUNC([dgesv])

  save_LIBS="$LIBS"
  LIBS="$lapack_lib_flags $blas_lib_flags $LIBS $FCLIBS"

  # Try linking with LAPACK and BLAS libraries
  #
  AC_MSG_CHECKING([for $dgesv in $lapack_lib_flags $blas_lib_flags])
  AC_LINK_IFELSE([AC_LANG_CALL([], [$dgesv])],
    [lib_lapack_dgesv="yes"],
    [lib_lapack_dgesv="no"])
  AC_MSG_RESULT([$lib_lapack_dgesv])

  LIBS="$save_LIBS"

  test "x$lib_lapack_dgesv" = xno &&
    AC_MSG_ERROR([cannot link against '$lapack_lib_flags $blas_lib_flags'])
  ;;
esac
AC_SUBST([lapack_lib_flags])


# MKL
#
use_mkl_flag=
mkl_include_dir=
mkl_lib_dir=
mkl_include_dir_flag=
mkl_lib_dir_flag=
mkl_cc_lib_flags=
mkl_fc_lib_flags=

if test "x$with_mkl" = xyes; then
  # Check if environment variable MKLROOT has been defined
  #
  test "x$MKLROOT" = x &&
    AC_MSG_WARN([environment variable MKLROOT is not set])

  if test "x$MKLINCLUDE" = x; then
    mkl_include_dir="$MKLROOT/include"
  else
    mkl_include_dir="$MKLINCLUDE"
  fi

  # Check for the presence of file mkl_dfti.f90
  #
  AC_CHECK_FILE([$mkl_include_dir/mkl_dfti.f90],
    [have_mkl_dfti_f90="yes"],
    [have_mkl_dfti_f90="no"])
  if test "x$have_mkl_dfti_f90" = no; then
    AC_MSG_ERROR([cannot find 'mkl_dfti.f90', please specify the MKL include \
directory with MKLINCLUDE=DIR])
  fi

  # Set MKL link lines
  # Link lines are set according to the Intel MKL Link Line Advisor at
  # http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/
  #

  # Fortran link line
  #
  # Interface
  #
  mkl_iface_suffix=
  case $host_cpu in
  x86_64 | ia64)
    mkl_iface_suffix="_lp64"
    ;;
  *)
    ;;
  esac

  case $FC in
  gfortran*)
    mkl_fc_iface="-lmkl_gf${mkl_iface_suffix}"
    ;;
  *)
    mkl_fc_iface="-lmkl_intel${mkl_iface_suffix}"
    ;;
  esac
 
  # Threading
  #
  mkl_fc_thread="-lmkl_intel_thread"

  # Basic libraries
  #
  mkl_fc_base="$mkl_fc_iface $mkl_fc_thread -lmkl_core"

  # Other libraries (common)
  #
  mkl_extra="-liomp5 -lpthread -lm"

  # Complete Fortran line
  #
  mkl_fc_lib_flags="-Wl,--start-group $mkl_fc_base -Wl,--end-group $mkl_extra"


  # C link line
  #
  # Interface
  #
  mkl_cc_iface="-lmkl_intel${mkl_iface_suffix}"

  # Threading
  #
  mkl_cc_thread="-lmkl_intel_thread"

  # Basic libraries
  #
  mkl_cc_base="$mkl_cc_iface $mkl_cc_thread -lmkl_core"

  # Complete C line
  #
  mkl_cc_lib_flags="-Wl,--start-group $mkl_cc_base -Wl,--end-group $mkl_extra"


  # Try to link against MKL libraries
  #
  # Save CPPFLAGS
  #
  save_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="-I$mkl_include_dir $CPPFLAGS"

  # Save LDFLAGS
  #
  save_LDFLAGS="$LDFLAGS"
  test "x$MKLPATH" != x &&
    mkl_lib_dir=$MKLPATH &&
    LDFLAGS="-L$mkl_lib_dir $LDFLAGS"

  save_LIBS="$LIBS"
  LIBS="$mkl_cc_lib_flags $LIBS"
  AC_MSG_CHECKING([whether the MKL libraries are available])
  AC_LANG([C])
  AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[
        #include <mkl.h>
      ]],[[
        MKLVersion ver;
        mkl_get_version(&ver);
      ]])],
    [mkl_link_success="yes"],
    [mkl_link_success="no"])

  if test "x$mkl_link_success" = xno; then
   AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[
        #include <mkl.h>
      ]],[[
        MKLVersion ver;
        MKLGetVersion(&ver);
      ]])],
    [mkl_link_success="yes"],
    [mkl_link_success="no"])
  fi
  AC_MSG_RESULT([$mkl_link_success])
  LIBS="$save_LIBS"

  test "x$mkl_link_success" = xno &&
    AC_MSG_ERROR([cannot link against the MKL libraries, please specify the \
MKL lib directory with MKLPATH=DIR])

  AC_DEFINE([USE_MKL], [1], [Define to 1 if MKL is required])
  use_mkl_flag="-DUSE_MKL"
  mkl_include_dir_flag="-I$mkl_include_dir"
  test "x$mkl_lib_dir" != x &&
    mkl_lib_dir_flag="-L$mkl_lib_dir"

  # Restore CPPFLAGS and LDFLAGS
  #
  CPPFLAGS="$save_CPPFLAGS"
  LDFLAGS="$save_LDFLAGS"
fi
AC_SUBST([use_mkl_flag])
AC_SUBST([mkl_include_dir])
AC_SUBST([mkl_lib_dir])
AC_SUBST([mkl_include_dir_flag])
AC_SUBST([mkl_lib_dir_flag])
AC_SUBST([mkl_cc_lib_flags])
AC_SUBST([mkl_fc_lib_flags])
AM_CONDITIONAL([USE_MKL], [test "x$use_mkl_flag" != x])

# Other libraries
#
intel_extra_fc_lib_flags=
case $FC in
gfortran*)
  test "x$CC" = xicc && intel_extra_fc_lib_flags="-lirc -liomp5"
  ;;
ifort)
  intel_extra_cc_lib_flags="-liomp5"
  ;;
esac
AC_SUBST([intel_extra_cc_lib_flags])
AC_SUBST([intel_extra_fc_lib_flags])
  
# Checks for header files.
#
AC_CHECK_HEADERS([sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
#
#AC_C_CONST
#AC_HEADER_TIME
#AC_STRUCT_TM

# Checks for library functions.
#
AC_CHECK_FUNCS([gettimeofday])

AC_MSG_NOTICE([Build using MKLINCLUDE=$MKLINCLUDE])
AC_MSG_NOTICE([Build using MKLPATH=$MKLPATH])
AC_MSG_NOTICE([Build using FC=$FC])
AC_MSG_NOTICE([Build using FCFLAGS=$FCFLAGS])
AC_MSG_NOTICE([Build using LDFLAGS=$LDFLAGS])
AC_MSG_NOTICE([Build using F77=$F77])
AC_MSG_NOTICE([Build using FFLAGS=$FFLAGS])
AC_MSG_NOTICE([Build using FCLIBS=$FCLIBS])
AC_MSG_NOTICE([Build using CC=$CC])
AC_MSG_NOTICE([Build using CFLAGS=$CFLAGS])
AC_MSG_NOTICE([Build using CPPFLAGS=$CPPFLAGS])
AC_MSG_NOTICE([Build using CPP=$CPP])
AC_MSG_NOTICE([Build using openmp_fc_flags=$openmp_fc_flags])
AC_MSG_NOTICE([Build using sched_flags=$sched_flags])
AC_MSG_NOTICE([Build using nofor_main_flag=$nofor_main_flag])
AC_MSG_NOTICE([Build using f95_cpp_input_flag=$f95_cpp_input_flag])
AC_MSG_NOTICE([Build using openmp_cc_flags=$openmp_cc_flags])
AC_MSG_NOTICE([Build using dfftpack_lib_flags=$dfftpack_lib_flags])
AC_MSG_NOTICE([Build using blas_lib_flags=$blas_lib_flags])
AC_MSG_NOTICE([Build using lapack_lib_flags=$lapack_lib_flags])
AC_MSG_NOTICE([Build using use_mkl_flag=$use_mkl_flag])
AC_MSG_NOTICE([Build using mkl_include_dir=$mkl_include_dir])
AC_MSG_NOTICE([Build using mkl_lib_dir=$mkl_lib_dir])
AC_MSG_NOTICE([Build using mkl_include_dir_flag=$mkl_include_dir_flag])
AC_MSG_NOTICE([Build using mkl_lib_dir_flag=$mkl_lib_dir_flag])
AC_MSG_NOTICE([Build using mkl_cc_lib_flags=$mkl_cc_lib_flags])
AC_MSG_NOTICE([Build using mkl_fc_lib_flags=$mkl_fc_lib_flags])
AC_MSG_NOTICE([Build using intel_extra_cc_lib_flags=$intel_extra_cc_lib_flags])
AC_MSG_NOTICE([Build using intel_extra_fc_lib_flags=$intel_extra_fc_lib_flags])

AC_CONFIG_FILES([Makefile
	utils/Makefile
	utils/dfftpack/Makefile
	src/Makefile
	src/dctii/Makefile
	src/dst/Makefile
	src/rft/Makefile
	src/conv/Makefile
	src/diagc/Makefile
	src/dts/Makefile
	src/dpis/Makefile
	src/dt/Makefile
	src/dtspg/Makefile
	mex/Makefile
	mex/create_mex_files.sh
	doc/Makefile
	examples/Makefile
	examples/create_examples.sh])

AC_OUTPUT
