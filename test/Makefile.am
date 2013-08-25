#test: check

SUBDIRS =unitTests apps

#AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../src/ \
        $(MODPATHFLAG)$(srcdir)/../src/math \
        $(MODPATHFLAG)$(srcdir)/../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../src/data \
	$(FCFLAGS) $(SYSFFLAGS)

#LDADD = \
	$(SYSLDFLAGS) $(SYSLDLIBS) \
	$(srcdir)/../src/libDQGModelCore.la 

#CLEANFILES = *$(MODEXT)

#TESTS = $(check_PROGRAMS)

#check_PROGRAMS = 

#check_PROGRAMS += advTest
#advTest_SOURCES =\
	advectTestOnSphere/SVoronoiGen_mod.f90 \
	advectTestOnSphere/HexTriIcMesh_mod.f90 \
	advectTestOnSphere/advTest.f90 