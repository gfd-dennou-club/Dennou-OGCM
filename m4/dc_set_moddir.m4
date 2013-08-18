dnl== Set MODDIR from libfile
dnl
dnl Check existence of "libfile". And set MODDIR,
dnl
dnl usage: DC_SET_MODDIR(libfile, MODDIR)
dnl
AC_DEFUN([DC_SET_MODDIR],[
        if test ! -f $1 ; then
                AC_MSG_ERROR(specified library file \"$1\" is not exist)
        fi
        DC_SET_MODDIR_libdir=`dirname $1`
        DC_SET_MODDIR_libname=`basename $1 .a | sed 's/^lib//'`
        DC_SET_MODDIR_try1=""; DC_SET_MODDIR_try2=""
        if test -d ${DC_SET_MODDIR_try1:=$DC_SET_MODDIR_libdir/module} ; then
                $2=$DC_SET_MODDIR_try1
        elif test -d ${DC_SET_MODDIR_try2:=`dirname $DC_SET_MODDIR_libdir`/include} ; then
                $2=$DC_SET_MODDIR_try2
        else
                AC_MSG_ERROR($DC_SET_MODDIR_libname module directory not found)
        fi
])
