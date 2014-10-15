SUBDIRS = src \
	test \
	tool \
	model/globalSWModel_FVM \
	model/globalSWModel_DG \
	model/ogcm \
	model/sogcm


ACLOCAL_AMFLAGS = -I m4

FC=@FC@
F90MODTYPE=@F90MODTYPE@
MODEXT=@MODEXT@
MODPATHFLAG=@MODPATHFLAG@
SYSFCFLAGS=@SYSFFLAGS@ 
SYSLDFLAGS=@SYSLDFLAGS@
SYSLDLIBS=@SYSLDLIBS@

SYSCXXFLAGS=@SYSCXXFLAGS@

#dist-hook:
#	sh $(top_srcdir)/dist_tool/generate-ChangeLog.sh       


test: check
	@echo ""
	@echo ""
	@echo " *** All unit tests are succeeded !! *** "
	@echo ""
	@echo ""


