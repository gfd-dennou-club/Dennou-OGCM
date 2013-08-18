# _DC_FORTRAN_MODULE_FLAG
#--------------------------------------------------------
# Authors:: Eizi TOYODA, Yasuhiro Morikawa, Youhei SASAKI
# Version:: $Id:$
# Copyright:: SPMODEL Development Group, All rights, reserved.
# License:: MIT-Like, See COPYRIGHT
#--------------------------------------------------------
AC_DEFUN([DC_FORTRAN_MODULE_FLAG],
[AC_CACHE_CHECK([Check Fortran's modules include flag],
                [ax_cv_fortran_modflag],
[AC_LANG_PUSH(Fortran)
i=0
while test \( -f tmpdir_$i \) -o \( -d tmpdir_$i \) ; do
  i=`expr $i + 1`
done
mkdir tmpdir_$i
cd tmpdir_$i
AC_COMPILE_IFELSE([module conftest_module
   contains
   subroutine conftest_routine
   character(len=4) :: a = "test"
   print '(a)', a
   end subroutine conftest_routine
   end module conftest_module
],[],[])
cd ..
ax_cv_fortran_modflag="not found"
for ax_flag in "-I" "-M" "-p"; do
  if test "$ax_cv_fortran_modflag" = "not found" ; then
    ax_save_FCFLAGS="$FCFLAGS"
    FCFLAGS="$ax_save_FCFLAGS ${ax_flag}tmpdir_$i"
    AC_COMPILE_IFELSE([program conftest_program
      use conftest_module
      call conftest_routine
      end program conftest_program
    ],[ax_cv_fortran_modflag="$ax_flag"],[])
    FCFLAGS="$ax_save_FCFLAGS"
  fi
done
rm -fr tmpdir_$i
AC_LANG_POP(Fortran)
])])
