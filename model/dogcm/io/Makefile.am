SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../include/ \
        $(MODPATHFLAG)$(srcdir)/../util           \
        $(MODPATHFLAG)$(srcdir)/../admin          \
        $(MODPATHFLAG)$(srcdir)/../seaice         \
	$(MODPATHFLAG)$(srcdir)/../boundary       \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

noinst_LTLIBRARIES= libDOGCMIO.la 
dist_libDOGCMIO_la_SOURCES= \
	DOGCM_IO_History_mod.f90 \
	DOGCM_IO_Restart_mod.f90

DOGCM_IO_History_mod.f90 : \
	../admin/DOGCM_Admin_GridDef_mod.o        \
	../admin/DOGCM_Admin_GaussSpmGrid_mod.o   \
	../admin/DOGCM_Admin_Grid_mod.o           \
	../admin/DOGCM_Admin_TInteg_mod.o
