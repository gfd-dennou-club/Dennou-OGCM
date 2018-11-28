ACLOCAL_AMFLAGS = -I m4

# These files get added to the distribution.
EXTRA_DIST =  COPYRIGHT ChangeLog VERSION
EXTRA_DIST += README INSTALL
EXTRA_DIST += doc 


FC=@FC@
F90MODTYPE=@F90MODTYPE@
MODEXT=@MODEXT@
MODPATHFLAG=@MODPATHFLAG@
SYSFFLAGS=@SYSFFLAGS@
SYSLDFLAGS=@SYSLDFLAGS@
SYSLDLIBS=@SYSLDLIBS@

SYSCXXFLAGS=@SYSCXXFLAGS@

#------------------------------------------------

SUBDIRS = \
	libsrc \
	test \
	tool \
	model/globalSWModel_FVM \
	model/globalSWModel_DG \
	model/ogcm \
	model/sogcm \
	model/dogcm

#	model/globalSWModel_DG2 

#------------------------------------------------

#dist-hook:
#	sh $(top_srcdir)/dist_tool/generate-ChangeLog.sh       


test: check
	@echo ""
	@echo ""
	@echo " *** All unit tests are succeeded !! *** "
	@echo ""
	@echo ""


install-data-hook:
	$(MKDIR_P) $(DESTDIR)$(includedir)
	$(install_sh_DATA) $(srcdir)/include/*@MODEXT@ $(DESTDIR)$(includedir)
