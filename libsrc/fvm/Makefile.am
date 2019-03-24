SUBDIRS = 

AM_FCFLAGS = $(SYSFFLAGS) \
	$(MODPATHFLAG)$(srcdir)/ \
	$(MODPATHFLAG)$(srcdir)/../math \
	$(MODPATHFLAG)$(srcdir)/../mesh \
	$(MODPATHFLAG)$(srcdir)/../data


CLEANFILES = *$(MODEXT)

lib_LTLIBRARIES = libDQGModelFVM_HexCGrid.la
libDQGModelFVM_HexCGrid_la_SOURCES = \
	Hex_C-Grid/fvMeshInfo_mod.f90 \
	Hex_C-Grid/fvCalculus_mod.f90



