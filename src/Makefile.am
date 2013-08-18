SUBDIRS = math mesh data

AM_FCFLAGS = $(SYSFFLAGS) \
	$(MODPATHFLAG)$(srcdir)/ 

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

lib_LTLIBRARIES = libDQGModelCore.la
libDQGModelCore_la_SOURCES = 
libDQGModelCore_la_LIBADD = \
	$(srcdir)/math/libDQGModelMath.la \
	$(srcdir)/mesh/libDQGModelMesh.la \
	$(srcdir)/data/libDQGModelData.la

