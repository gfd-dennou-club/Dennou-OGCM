SUBDIRS = math mesh data fvm mesh2 io

AM_FCFLAGS = \
	$(MODPATHFLAG)$(srcdir)/ \
	$(MODPATHFLAG)$(srcdir)/../math \
	$(MODPATHFLAG)$(srcdir)/../mesh \
	$(MODPATHFLAG)$(srcdir)/../data \
	$(SYSFFLAGS)


#AM_CXXFLAGS = $(SYSCXXFLAGS) \
#	-I $(srcdir)/data/ \
#	-I $(srcdir)/util/ \
#	-I $(srcdir)/io/ \
#	-I $(srcdir)/dynamics/ \
#	-I $(srcdir)/experiment/ \
#	-I $(srcdir)/

#LDADD = \
#	$(SYSLDFLAGS) $(SYSLDLIBS) \
#	$(srcdir)/libDIGMOBaseLib.la 

CLEANFILES = *$(MODEXT)

#bin_PROGRAMS = IGModel-NHRM
#
#IGModel_NHRM_SOURCES = \
#	main.f90

lib_LTLIBRARIES = \
	libDQGModelCore.la

libDQGModelCore_la_SOURCES = 
libDQGModelCore_la_LIBADD = \
	$(srcdir)/math/libDQGModelMath.la \
	$(srcdir)/mesh/libDQGModelMesh.la \
	$(srcdir)/data/libDQGModelData.la \
	$(srcdir)/fvm/libDQGModelFVM.la \
	$(srcdir)/mesh2/libDQGModelMesh2.la
