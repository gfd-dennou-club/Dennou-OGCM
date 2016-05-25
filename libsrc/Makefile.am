SUBDIRS = math mesh data fvm mesh2 io

all: copy_header

copy_header:
		@echo "Copy generated header files into $(srcdir)/include"
		@$(MKDIR_P) $(srcdir)/../include
		@find . -name "*@MODEXT@" -exec $(install_sh_DATA) {} $(srcdir)/../include \;

#AM_FCFLAGS = \
	$(MODPATHFLAG)$(srcdir)/ \
	$(MODPATHFLAG)$(srcdir)/../math \
	$(MODPATHFLAG)$(srcdir)/../mesh \
	$(MODPATHFLAG)$(srcdir)/../data \
	$(SYSFFLAGS)


CLEANFILES = *$(MODEXT)

lib_LTLIBRARIES = \
	libDQGModelCore.la

libDQGModelCore_la_SOURCES = 
libDQGModelCore_la_LIBADD = \
	$(srcdir)/math/libDQGModelMath.la \
	$(srcdir)/mesh/libDQGModelMesh.la \
	$(srcdir)/data/libDQGModelData.la \
	$(srcdir)/mesh2/libDQGModelMesh2.la
