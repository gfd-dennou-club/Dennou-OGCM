# _DC_FORTRAN_MODULE_TYPE
#---------------------------------------------------------------
# Authors:: Eizi TOYODA, Yasuhiro Morikawa, Youhei SASAKI
# Version:: $Id:$
# Copyright:: SPMODEL Development Group, All rights, reserved.
# License:: MIT-Like, See COPYRIGHT
#---------------------------------------------------------------
AC_DEFUN([DC_FORTRAN_MODULE_TYPE],
[AC_CACHE_CHECK([Check $FC's module type],
                [ax_cv_fortran_modtype],
[AC_LANG_PUSH(Fortran)
i=0
while test \( -f tmpdir_$i \) -o \( -d tmpdir_$i \) ; do
  i=`expr $i + 1`
done
mkdir tmpdir_$i
cd tmpdir_$i
cat <<EOF > conftes1.f90
module conftesa
logical :: b = .false.
end module conftesa
EOF
$FC -c conftes1.f90 1> /dev/null 2>&1
ax_cv_fortran_modtype="not found"
if test -f conftes1.d ; then
  ax_cv_fortran_modtype=intel.d
elif test -f CONFTESA.mod ; then
  ax_cv_fortran_modtype=HP.mod
elif test -f conftesa.mod ; then
  ax_cv_fortran_modtype=std.mod
else
  cat <<EOF > conftes2.f90
program conftes2
use conftesa, only: b
b = .true.
end program conftes2
EOF
  ln conftes1.f90 conftesa.f90
  if $FC -c conftes2.f90 1>/dev/null 2>&1 && test -f contes2.o  ; then
     ax_cv_fortran_modtype=hitachi.f90
  elif $FC -c -Am conftes1.f90 && $FC -c -Am conftes2.f90 1>/dev/null 2>&1 ;then
     ax_cv_fortran_modtype=fqs.mod
  fi
fi
cd ..
rm -fr tmpdir_$i
])])
