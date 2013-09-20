SUBDIRS = 

AM_FCFLAGS = $(SYSFFLAGS) \
	$(MODPATHFLAG)$(srcdir)/ \
	$(MODPATHFLAG)$(srcdir)/../math \
	$(MODPATHFLAG)$(srcdir)/../mesh \
	$(MODPATHFLAG)$(srcdir)/../data \
	$(MODPATHFLAG)$(srcdir)/../fvm 
 


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

noinst_LTLIBRARIES = libDQGModelMesh2.la
nodist_libDQGModelMesh2_la_SOURCES = \
	HexTriIcMesh_mod.f90

