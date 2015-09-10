SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../src/ \
        $(MODPATHFLAG)$(srcdir)/../../../src/math \
        $(MODPATHFLAG)$(srcdir)/../../../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../../../src/data \
        $(MODPATHFLAG)$(srcdir)/../../../src/fvm \
        $(MODPATHFLAG)$(srcdir)/../../../src/mesh2 \
        $(MODPATHFLAG)$(srcdir)/../../../src/io \
        $(MODPATHFLAG)$(srcdir)/../setup \
        $(MODPATHFLAG)$(srcdir)/../util \
        $(MODPATHFLAG)$(srcdir)/../dynamics \
        $(MODPATHFLAG)$(srcdir)/../experiment \
        $(MODPATHFLAG)$(srcdir)/../physics_driver \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

noinst_LTLIBRARIES= libDSOGCMSeaIce.la 
nodist_libDSOGCMSeaIce_la_SOURCES= \
	SeaIceConstants_mod.f90 \
	SeaIceBoundaryCond_mod.f90 \
	SeaIceThermDyn_Winton2000_mod.f90 \
	VarSetSeaice_mod.f90


#seaice_thermdyn_mom_SOURCES=\
	SeaIceThermDyn_mom_main.f90 \
	ice_thm.F90

#seaice_thermdyn_mom_LDADD=\
	$(LDADD) \
	libDSOGCMSEAICE.la
