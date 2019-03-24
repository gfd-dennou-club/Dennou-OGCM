## This is a automake file, part of DCPAM.
# Copyright 2018, @see COPYRIGHT

# This directory stores libtool macros, put there by aclocal.
ACLOCAL_AMFLAGS = -I m4

# These files get added to the distribution.
EXTRA_DIST =  COPYRIGHT ChangeLog VERSION
EXTRA_DIST += README INSTALL
EXTRA_DIST += doc

FC=@FC@
F90MODTYPE=@F90MODTYPE@
MODEXT=@MODEXT@
MODPATHFLAG=@MODPATHFLAG@
W3M=@W3M@
RD2=@RD2@
SYSFFLAGS=@SYSFFLAGS@
SYSLDFLAGS=@SYSLDFLAGS@
SYSLDLIBS=@SYSLDLIBS@
SYSCXXFLAGS=@SYSCXXFLAGS@

#------------------------------------------------

SUBDIRS = \
	libsrc                  \
	test                    \
	tool                    \
	model/globalSWModel_FVM \
	model/globalSWModel_DG  \
	model/dogcm

#	model/globalSWModel_DG2 

#------------------------------------------------

#dist-hook:
#	sh $(top_srcdir)/dist_tool/generate-ChangeLog.sh


update-doc:
	( cd $(srcdir)/doc && make update )

	
test: check
	@echo ""
	@echo ""
	@echo " *** All unit tests for libsrc are succeeded !! *** "
	@echo ""
	@echo ""


install-data-hook:
	$(MKDIR_P) $(DESTDIR)$(includedir)
	$(install_sh_DATA) $(srcdir)/include/*@MODEXT@ $(DESTDIR)$(includedir)
	
clean-local:
	-@rm -f $(srcdir)/include/*@MODEXT@

maintainer-clean-local:
	rm -fr INSTALL INSTALL.ja
#	rm -fr doc/*.pdf doc/*.html doc/*.html.ja

