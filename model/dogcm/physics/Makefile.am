SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../src/ \
        $(MODPATHFLAG)$(srcdir)/../util     \
	$(MODPATHFLAG)$(srcdir)/../admin    \
	$(MODPATHFLAG)$(srcdir)/../boundary \
        $(MODPATHFLAG)$(srcdir)/../dynamics \
        $(MODPATHFLAG)$(srcdir)/../io \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

noinst_LTLIBRARIES= \
	libDOGCMLPhys.la   \
	libDOGCMVPhys.la   \
	libDOGCMPhysics.la

libDOGCMPhysics_la_SOURCES=\
	DOGCM_Phys_spm_mod.f90       \
	DOGCM_Phys_hspm_vfvm_mod.f90 \
	DOGCM_Phys_driver_mod.f90

libDOGCMPhysics_la_LIBADD = \
	libDOGCMLPhys.la \
	libDOGCMVPhys.la

dist_libDOGCMLPhys_la_SOURCES= \
	lphys/LPhys_DIFF_spm_mod.f90           \
	lphys/DOGCM_LPhys_DIFF_driver_mod.f90  \
	lphys/LPhys_RediGMHelper_mod.f90       \
	lphys/LPhys_RediGM_hspm_vfvm_mod.f90   \
	lphys/LPhys_RediGM_spm_mod.f90       

dist_libDOGCMVPhys_la_SOURCES= \
	vphys/DOGCM_VPhys_ConvAdjust_mod.f90  \
	vphys/DOGCM_VPhys_driver_mod.f90

#	vmixing/SGSSlowConvAdjust_mod.f90 
