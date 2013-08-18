# _DC_ARG_ENABLE
#---------------------------------------------------------------
# Authors:: Eizi TOYODA, Yasuhiro Morikawa, Youhei SASAKI
# Version:: $Id:$
# Copyright:: SPMODEL Development Group, All rights, reserved.
# License:: MIT-Like, See COPYRIGHT
#---------------------------------------------------------------
AC_DEFUN([DC_ARG_ENABLE],
[AC_ARG_ENABLE($1, [ --enable-$1: $2],
[$3=$enableval],
[AC_CACHE_CHECK([$2], $3, [$4])]
)])
