SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../include/ \
        $(MODPATHFLAG)$(srcdir)/../seaice         \
        $(MODPATHFLAG)$(srcdir)/../util           \
        $(MODPATHFLAG)$(srcdir)/../io             \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

noinst_LTLIBRARIES= libDOGCMAdmin.la 
dist_libDOGCMAdmin_la_SOURCES= \
	DOGCM_Admin_Constants_mod.f90    \
	DOGCM_Admin_GridDef_mod.f90      \
	DOGCM_Admin_GaussSpmGrid_mod.f90 \
	DOGCM_Admin_Grid_mod.f90         \
	DOGCM_Admin_TInteg_mod.f90       \
	DOGCM_Admin_BC_mod.f90           \
	DOGCM_Admin_Variable_mod.f90     \
	DOGCM_Admin_GovernEq_mod.f90

DOGCM_Admin_Grid_mod.o:                 \
	DOGCM_Admin_GaussSpmGrid_mod.o

DOGCM_Admin_Variable_mod.o:             \
	DOGCM_Admin_TInteg_mod.o        \
	../io/DOGCM_IO_History_mod.o 


