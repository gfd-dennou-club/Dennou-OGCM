W3M=@W3M@
RD2=@RD2@
RD2OPTS=-r rd/rd2html-ext-lib --ref-extension --with-part=HTML:html
RD2OPTS+=--native-inline --out-code=utf8
RD2OPTS_JA=$(RD2OPS) --with-part=JA:rd
W3MOPTSJA =  -T text/html -I u -O u -cols 72 -dump
RD2OPTS_EN=$(RD2OPS) --with-part=EN:rd
W3MOPTSEN =  -T text/html -I u -O N -cols 72 -dump

#reference: $(srcdir)/../html/reference_manual/index.html
#$(srcdir)/../html/reference_manual/index.html:
#	( cd $(srcdir)/../ && $(RDOC) $(RDOCOPTS) --main $(MAINPAGE) $(F90FILES) )

update: gen-doc-local gen-description gen-usersguide

gen-doc-local: \
	$(srcdir)/../INSTALL    \
	$(srcdir)/../INSTALL.ja 


$(srcdir)/../INSTALL: INSTALL.rd
	$(RD2) $(RD2OPTS_EN) $< | $(W3M) $(W3MOPTSEN) > $@

$(srcdir)/../INSTALL.ja: INSTALL.rd
	$(RD2) $(RD2OPTS_JA) $< | $(W3M) $(W3MOPTSJA) > $@

gen-description:
	( cd $(srcdir)/description && make )

gen-usersguide:
	( cd $(srcdir)/users-guide && make )
		