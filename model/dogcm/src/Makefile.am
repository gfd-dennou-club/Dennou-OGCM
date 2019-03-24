SUBDIRS = \
	util admin io boundary seaice dynamics physics tint

################

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../../include/ \
        $(MODPATHFLAG)$(srcdir)/admin             \
        $(MODPATHFLAG)$(srcdir)/util              \
        $(MODPATHFLAG)$(srcdir)/io                \
	$(MODPATHFLAG)$(srcdir)/boundary          \
	$(MODPATHFLAG)$(srcdir)/dynamics          \
        $(MODPATHFLAG)$(srcdir)/physics           \
        $(MODPATHFLAG)$(srcdir)/seaice            \
        $(MODPATHFLAG)$(srcdir)/tint              \
        $(MODPATHFLAG)$(srcdir)/exp               \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/../../../libsrc/fvm/libDQGModelFVM_HexCGrid.la \
	$(srcdir)/../../../libsrc/libDQGModelCore.la             \
	$(SYSLDFLAGS) $(SYSLDLIBS)

#####################################################################

CLEANFILES = *$(MODEXT)
all: all-am copy_header

copy_header:
		@echo "Copy generated header files into $(srcdir)/include"
		@$(MKDIR_P) $(srcdir)/../include
		@find . -name "*@MODEXT@" -exec $(install_sh_DATA) {} $(srcdir)/../include \;


#######################################################################

noinst_LTLIBRARIES = libDOGCMCore.la

#-----------------------------------------------------------------

libDOGCMCore_la_SOURCES =

libDOGCMCore_la_LIBADD = \
	$(srcdir)/tint/libDOGCMTInt.la            \
	$(srcdir)/physics/libDOGCMPhysics.la      \
	$(srcdir)/dynamics/libDOGCMDynamics.la    \
	$(srcdir)/io/libDOGCMIO.la                \
	$(srcdir)/boundary/libDOGCMBoundary.la    \
	$(srcdir)/admin/libDOGCMAdmin.la          \
	$(srcdir)/util/libDOGCMUtil.la            \
	$(srcdir)/util/libSpmlUtil.la         	  \
	$(LDADD)


