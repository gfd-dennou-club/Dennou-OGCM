subdir =
test: check

AM_FCFLAGS = \
        $(MODPATHFLAG)$(srcdir)/../src/ \
        $(MODPATHFLAG)$(srcdir)/../src/math \
        $(MODPATHFLAG)$(srcdir)/../src/mesh \
        $(MODPATHFLAG)$(srcdir)/../src/data \
        $(MODPATHFLAG)$(srcdir)/../src/fvm \
        $(MODPATHFLAG)$(srcdir)/../src/mesh2 \
        $(MODPATHFLAG)$(srcdir)/../src/io \
	$(FCFLAGS) $(SYSFFLAGS)

LDADD = \
	$(srcdir)/../src/io/libDQGModelIO.la \
	$(srcdir)/../src/libDQGModelCore.la \
	$(SYSLDFLAGS) $(SYSLDLIBS)

CLEANFILES = *$(MODEXT)

TESTS = $(check_PROGRAMS)

bin_PROGRAMS = 

bin_PROGRAMS += scvgridGen
scvgridGen_SOURCES =\
	gridGen/gridGen.f90

bin_PROGRAMS += advTest 
advTest_SOURCES =\
	advectTestOnSphere/advTest.f90 