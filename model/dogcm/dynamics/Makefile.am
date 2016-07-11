SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../include/ \
        $(MODPATHFLAG)$(srcdir)/../util           \
        $(MODPATHFLAG)$(srcdir)/../admin          \
	$(MODPATHFLAG)$(srcdir)/../boundary       \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

noinst_LTLIBRARIES= libDOGCMDynamics.la 
dist_libDOGCMDynamics_la_SOURCES= \
	HBEDiagnose_spm_mod.f90 \
	HBEBaroc_spm_mod.f90    \
	HBEBarot_spm_mod.f90    \
	DOGCM_Dyn_spm_mod.f90   \
	DOGCM_Dyn_driver_mod.f90


