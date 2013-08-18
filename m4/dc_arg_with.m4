# _DC_ARG_WITH
#---------------------------------------------------------------
# Authors:: Eizi TOYODA, Yasuhiro Morikawa, Youhei SASAKI
# Version:: $Id:$
# Copyright:: SPMODEL Development Group, All rights, reserved.
# License:: MIT-Like, See COPYRIGHT
#---------------------------------------------------------------
AC_DEFUN([DC_ARG_WITH],
[AC_ARG_WITH([$1],
[AC_HELP_STRING([--with-$1], [$2])],
[$3=$withval],
[AC_CACHE_CHECK([$2], $3, [$4])
])])
