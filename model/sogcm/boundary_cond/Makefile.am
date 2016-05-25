SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../libsrc/ \
        $(MODPATHFLAG)$(srcdir)/../util \
        $(MODPATHFLAG)$(srcdir)/../setup \
        $(MODPATHFLAG)$(srcdir)/../seaice \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

noinst_LTLIBRARIES= libDSOGCMBC.la 
dist_libDSOGCMBC_la_SOURCES= \
	BoundaryCondO_mod.f90

BoundaryCondO_mod.f90: \
	../seaice/VarSetSeaice_mod.f90 \
	../seaice/SeaIceBoundaryCond_mod.f90
