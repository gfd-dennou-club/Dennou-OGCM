SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../libsrc/ \
        $(MODPATHFLAG)$(srcdir)/../util          \
        $(MODPATHFLAG)$(srcdir)/../admin         \
        $(MODPATHFLAG)$(srcdir)/../io            \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

noinst_LTLIBRARIES= libDOGCMBoundary.la 
dist_libDOGCMBoundary_la_SOURCES= \
	DOGCM_Boundary_vars_mod.f90   \
	DOGCM_Boundary_common_mod.f90 \
	DOGCM_Boundary_spm_mod.f90    \
	DOGCM_Boundary_hspm_vfvm_mod.f90    \
	DOGCM_Boundary_driver_mod.f90

#BoundaryCondO_mod.f90: \
	../seaice/VarSetSeaice_mod.f90 \
	../seaice/SeaIceBoundaryCond_mod.f90
