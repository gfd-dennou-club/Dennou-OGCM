SUBDIRS =

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../src/ \
        $(MODPATHFLAG)$(srcdir)/../../../src/math \
        $(MODPATHFLAG)$(srcdir)/../../../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../../../src/data \
        $(MODPATHFLAG)$(srcdir)/../../../src/fvm \
        $(MODPATHFLAG)$(srcdir)/../../../src/mesh2 \
        $(MODPATHFLAG)$(srcdir)/../../../src/io \
        $(MODPATHFLAG)$(srcdir)/../../ogcm/util \
        $(MODPATHFLAG)$(srcdir)/../setup \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

noinst_LTLIBRARIES= libDSOGCMUtil.la 
nodist_libDSOGCMUtil_la_SOURCES= \
	SpmlUtil_mod.f90 \
	DiagnoseUtil_mod.f90 \
	EqState_Linear_mod.f90 \
	EOSDriver_mod.f90

DiagnoseUtil_mod.f90: \
	../setup/Constants_mod.o \
	../setup/GridSet_mod.o

EOSDriver_mod.f90: \
	../setup/Constants_mod.o \
	../setup/GridSet_mod.o \
	../setup/GovernEqSet_mod.o 


