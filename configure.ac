#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
# Package name, version, email address for bug report
AC_INIT([DQGModel],[1.0.0],[ykawai@gfd-dennou.org])

# Set macro directory
AC_CONFIG_MACRO_DIR([m4])
# Set build directory
AC_CONFIG_AUX_DIR([build-aux])
# Check cross-compile
AC_CANONICAL_HOST
# Call automake
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
# Turn off shared libraries by default
#AC_DISABLE_SHARED([])
# Checks for programs.
AC_PROG_MAKE_SET
# Check fortran compiler
AC_MSG_NOTICE([finding fortran compiler])
case "$host" in
  *solaris*)
        F77S="xlf fort xlf95 ifort ifc efc pgf77 lf95 f95 f90 g95 gfortran f77 frt pgf95 fort77 fl32 af77 f90 xlf90 pgf90 epcf90 g77"
        F90S="xlf90 fort xlf95 ifort ifc efc pgf95 lf95 f95 f90 g95 gfortran xlf90 pgf90 epcf90"
        ;;
  *)
        F77S="xlf fort xlf95 ifort ifc efc pgf77 lf95 g95 gfortran f77 frt pgf95 f95 f90 fort77 fl32 af77 f90 xlf90 pgf90 epcf90 g77"
        F90S="xlf90 fort xlf95 ifort ifc efc pgf95 lf95 g95 gfortran f95 f90 f90 xlf90 pgf90 epcf90"
        ;;
esac
AC_PROG_FC($F90S)
if test "${FC+set}" = set && test "x$FC" = x; then
   AC_MSG_ERROR([Can't find F90 compiler.])
fi
AC_MSG_NOTICE([Your Fortran Compiler is $FC])
# Check the fortran >= 2003 compiler
AC_MSG_CHECKING([Check $FC can handle fortran >= 90 program.])
cat <<EOF >conftest.f90
program conftest
real, parameter :: a = 3.141592654
end program
EOF
ac_try='$FC -c ${FCFLAGS} conftest.f90'
if AC_TRY_EVAL(ac_try); then
  f90_check=OK
else
  f90_check=NG
fi
rm -f a.out conftest.*
AC_MSG_RESULT([$f90_check])
if test "x$f90_check" = xNG; then
  AC_MSG_ERROR([$FC is NOT a fortran >= 90 compiler])
fi
AC_MSG_NOTICE([Check $FC's module type])
AC_LANG_PUSH(Fortran)
AC_FC_SRCEXT(f90)
AC_LANG_POP(Fortran)
# @see m4/dc_fortran_check_modtype.m4
DC_FORTRAN_MODULE_TYPE
if test x"$ax_cv_fortran_modtype" = x"not found" ; then
  AC_MSG_ERROR([unable to find $FC's module type])
fi
AC_MSG_NOTICE([Fortran module type is $ax_cv_fortran_modtype])
AC_SUBST(F90MODTYPE, ${ax_cv_fortran_modtype})
case "${ax_cv_fortran_modtype:-undef}" in
intel.d)
  ax_cv_modext=.d
  ;;
HP.mod|std.mod|fqs.mod)
  ax_cv_modext=.mod
  ;;
hitachi.f90)
  ax_cv_modext=.f90
  ;;
esac
AC_SUBST(MODEXT, ${ax_cv_modext})
AC_MSG_NOTICE([Check $FC's module include flag])
AC_LANG_PUSH(Fortran)
AC_FC_SRCEXT(f90)
AC_LANG_POP(Fortran)
# @see m4/dc_fortran_check_modflag.m4
DC_FORTRAN_MODULE_FLAG
if test x"$ax_cv_fortran_modflag" = x"not found" ; then
  AC_MSG_ERROR([unable to find $FC module include flag])
fi
AC_SUBST(MODPATHFLAG, ${ax_cv_fortran_modflag})
# Check C/C++ compiler
AC_PROG_CC
AC_PROG_CXX

AC_SUBST(CPPHEADEREXT, ".hpp")
# Checks for libraries.
#
# * netcdf
#
DC_ARG_WITH([netcdf],
            [NetCDF library filename],
            [ac_cv_lib_netcdf],
  [AC_MSG_ERROR([specify netcdf library filename with --with-netcdf=])])
#
# * netcdff
#
DC_ARG_WITH([netcdff],
            [NetCDF Fortran library filename],
            [ac_cv_lib_netcdff],
            [ac_cv_lib_netcdff=""])
#
# * netcdf include
#
DC_ARG_WITH([netcdf-include],
            [NetCDF Fortran header filename],
            [ac_cv_netcdf_include],
            [ac_cv_netcdf_include=""])
#
# * Set NetCDF library, header...
#
DC_SET_LIBDIR_LIBNAME([$ac_cv_lib_netcdf], [netcdf_libdir], [netcdf_libname])
if test -n "$ac_cv_lib_netcdff"; then
  DC_SET_LIBDIR_LIBNAME([$ac_cv_lib_netcdff],
                        [netcdff_libdir],
                        [netcdff_libname])
fi
if test -n "$ac_cv_netcdf_include" && test -f "$ac_cv_netcdf_include"; then
  netcdf_incdir=`dirname $ac_cv_netcdf_include`
else
  try2=`dirname $netcdf_libdir`/include
  try2="$try2 `echo $netcdf_libdir | sed -e 's|/lib/|/include/|' -e 's|/lib$|/include|'`"
  for try3 in $try2; do
    if test -f "$try3/netcdf.inc"; then
      netcdf_incdir=$try3
      break
    fi
  done
fi
test -d "$netcdf_incdir" || AC_MSG_ERROR([netcdf include file not found. specify netcdf include filename with --with-netcdf-include=])
#
# * gtool5
#
DC_ARG_WITH([gtool5],
            [gtool5 library filename],
            [ac_cv_lib_gtool5],
  [AC_MSG_ERROR([specify gtool5 library filename with --with-gtool5= option.])])
DC_SET_LIBDIR_LIBNAME($ac_cv_lib_gtool5, gtool5_libdir, gtool5_libname)
DC_SET_MODDIR($ac_cv_lib_gtool5, gtool5_moddir)


#
# * --enable-mpi
#
DC_ARG_ENABLE([mpi], 
				[compile with MPI library [disable]], 
				[ac_cv_mpi], [ac_cv_mpi=default])

mpi_autosetting=no
case $ac_cv_mpi in
yes)
        SYSFFLAGS="${SYSFFLAGS} -DLIB_MPI"
        SYSCXXFLAGS="${SYSCXXFLAGS} -DLIB_MPI"
        ;;
no)
        ;;
*)
        build_mpicode=yes
        if test -n "${MPIFC}" ; then
            case "$FC" in
            *mpi*)
                SYSFFLAGS="${SYSFFLAGS} -DLIB_MPI"
                SYSCXXFLAGS="${SYSCXXFLAGS} -DLIB_MPI"
                mpi_autosetting=yes
                ;;
            *)
                CPPFLAGS=
            esac
        fi
        ;;
esac

#
# * set flags
#
if test -n "$ac_cv_lib_netcdff"; then
   SYSLDFLAGS="-L$netcdff_libdir $SYSLDFLAGS"
   SYSLDLIBS="-l$netcdff_libname $SYSLDLIBS"
fi
SYSLDFLAGS="-L$gtool5_libdir $SYSLDFLAGS -L$netcdf_libdir"
SYSLDLIBS="-l$gtool5_libname $SYSLDLIBS -l$netcdf_libname"
SYSFFLAGS="$SYSFFLAGS $ax_cv_fortran_modflag$netcdf_incdir $ax_cv_fortran_modflag$gtool5_moddir "

AC_SUBST(SYSLDLIBS)
AC_SUBST(SYSLDFLAGS)
AC_SUBST(SYSFFLAGS)

# Checks for header files.
# Checks for typedefs, structures, and compiler characteristics.
# Checks for library functions.
# * Calling libtool

LT_INIT
AC_CONFIG_FILES([
                 Makefile
                 src/Makefile
		 src/math/Makefile
                 src/data/Makefile
		 src/mesh/Makefile
		 src/fvm/Makefile
		 src/mesh2/Makefile
		 src/io/Makefile
		 test/Makefile
		 test/unitTests/Makefile
		 tool/Makefile
		 model/globalSWModel/Makefile
		 model/ogcm/Makefile
		 model/ogcm/setup/Makefile
		 model/ogcm/util/Makefile
		 model/ogcm/dynamics/Makefile
 ])

AC_OUTPUT

