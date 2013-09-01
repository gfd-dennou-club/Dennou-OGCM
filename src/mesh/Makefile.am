SUBDIRS = 

AM_FCFLAGS = $(SYSFFLAGS) \
	$(MODPATHFLAG)$(srcdir)/ \
	$(MODPATHFLAG)$(srcdir)/../math 

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

#IGModel_NHRM_SOURCES = \
	main.f90

noinst_LTLIBRARIES = libDQGModelMesh.la
nodist_libDQGModelMesh_la_SOURCES = \
	PolyMesh_mod.f90 \
	SVoronoiGen_mod.f90 \
	SVoronoiGen2_mod.f90 \
	SCVoronoiGen_mod.f90