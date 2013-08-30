subdir =
test: check

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../../src/ \
        $(MODPATHFLAG)$(srcdir)/../../src/math \
        $(MODPATHFLAG)$(srcdir)/../../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../../src/data \
        $(MODPATHFLAG)$(srcdir)/../../src/fvm \
        $(MODPATHFLAG)$(srcdir)/../../src/io \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/../../src/io/libDQGModelIO.la \
	$(srcdir)/../../src/libDQGModelCore.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

TESTS = $(check_PROGRAMS)

check_PROGRAMS = 

check_PROGRAMS += advTest
advTest_SOURCES =\
	advectTestOnSphere/SVoronoiGen_mod.f90 \
	advectTestOnSphere/SVoronoiGen2_mod.f90 \
	advectTestOnSphere/SCVoronoiGen_mod.f90 \
	advectTestOnSphere/HexTriIcMesh_mod.f90 \
	advectTestOnSphere/advTest.f90 