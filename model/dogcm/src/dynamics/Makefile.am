SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../../include/ \
        $(MODPATHFLAG)$(srcdir)/../util              \
        $(MODPATHFLAG)$(srcdir)/../admin             \
	$(MODPATHFLAG)$(srcdir)/../boundary          \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

noinst_LTLIBRARIES= libDOGCMDynamics.la 
dist_libDOGCMDynamics_la_SOURCES= \
	hspm_vspm/HBEDiagnose_spm_mod.f90 \
	hspm_vspm/HBEBaroc_spm_mod.f90    \
	hspm_vspm/HBEBarot_spm_mod.f90    \
	hspm_vspm/DOGCM_Dyn_spm_mod.f90   \
	hspm_vfvm/HBEDiagnose_hspm_vfvm_mod.f90 \
	hspm_vfvm/HBEBaroc_hspm_vfvm_mod.f90    \
	hspm_vfvm/HBEBarot_hspm_vfvm_mod.f90    \
	hspm_vfvm/DOGCM_Dyn_hspm_vfvm_mod.f90   \
	DOGCM_Dyn_driver_mod.f90


