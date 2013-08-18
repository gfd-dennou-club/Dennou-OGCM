dnl == Check libfile and set LIBDIR and  LIBNAME
dnl
dnl  Check existence of "libfile" file, and set LIBDIR, LIBNAME
dnl  from "libfile". "libfile" file must have suffixes ".a" or ".so"
dnl
dnl  usage: DC_SET_LIBDIR_LIBNAME(libfile, LIBDIR, LIBNAME)
dnl
AC_DEFUN([DC_SET_LIBDIR_LIBNAME],[
        if test ! -f $1 ; then
                AC_MSG_ERROR(specified library file \"$1\" is not exist)
        fi
        $2=`dirname $1`
        case "$1" in
            *.a)
                $3=`basename $1 .a | sed 's/^lib//'`
                ;;
            *.so)
                $3=`basename $1 .so | sed 's/^lib//'`
                ;;
            *)
                AC_MSG_ERROR(specified library file \"$1\" have invalid suffix. Valid suffixes are \".a\" or \".so\")
                ;;
        esac
])

