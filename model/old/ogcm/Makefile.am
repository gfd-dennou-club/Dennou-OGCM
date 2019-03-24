SUBDIRS = setup util dynamics

all: copy_header

copy_header:
		@echo "Copy generated header files into $(srcdir)/include"
		@$(MKDIR_P) $(srcdir)/include
		@find . -name "*@MODEXT@" -exec $(install_sh_DATA) {} $(srcdir)/include \;

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../include/ \
        $(MODPATHFLAG)$(srcdir)/setup \
        $(MODPATHFLAG)$(srcdir)/util \
        $(MODPATHFLAG)$(srcdir)/dynamics \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/dynamics/libDOGCMDynamics.la \
	$(srcdir)/util/libDOGCMUtil.la \
	$(srcdir)/setup/libDOGCMSetup.la \
	$(srcdir)/../../libsrc/io/libDQGModelIO.la \
	$(srcdir)/../../libsrc/fvm/libDQGModelFVM_HexCGrid.la \
	$(srcdir)/../../libsrc/libDQGModelCore.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

bin_PROGRAMS = 

bin_PROGRAMS += dogcm 
dogcm_SOURCES =\
	GovernEqSolverDriver_mod.f90 \
	ogcm_main.f90

