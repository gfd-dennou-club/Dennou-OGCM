###########################################################
#
# Makefile for compiling tex document
#
###########################################################

# target file name
TARGET = dcpom_users-guide
VERFNAME = version.tex
VERTEMPL = version_template.tex

# commands to compile documents
LATEX    = platex
BIBTEX   = bibtex
DVIPDF   = dvipdfmx -p a4
CURL     = curl -f -s -S --create-dirs
EXTRACTBB  = extractbb
VERSION  = $(shell cat ../../VERSION)
IMAGE_URL = http://www.gfd-dennou.org/arch/dcpom/doc_figs
IMAGE_PATH = $(IMAGE_URL)/DCPOM_v$(VERSION)/users-guide

IMAGE_BASENAME = $(shell grep includegraphics *tex | ruby -n -e '/^\s*\\includegraphics.*\{(.+)\}/ =~ $$_.split(":")[1] && print(/\A([^.]+)(\..+)?\Z/=~File.basename($$1) && $$1+($$2 || ".eps"),"\n")')
IMAGES = $(patsubst %,figs/%,$(IMAGE_OGCM_BASENAME)) 


SRCS = ${TARGET}.tex \
	introduction.tex   \
	install.tex        \
	$(IMAGES)
	

all: ${TARGET}.pdf

clean:
	rm -f ${TARGET}.aux ${TARGET}.bbl ${TARGET}.blg ${TARGET}.dvi ${TARGET}.log ${TARGET}.out ${TARGET}.toc

allclean: clean
	rm -rf *.pdf figure

					
${TARGET}.pdf: ${TARGET}.dvi
	$(DVIPDF) ${TARGET}

${TARGET}.dvi: $(VERFNAME) ${SRCS} reference.bib
	$(LATEX)  ${TARGET}
	$(BIBTEX) ${TARGET}
	$(LATEX)  ${TARGET}
	$(LATEX)  ${TARGET}

$(VERFNAME): $(VERTEMPL)
	$(shell cat $< | sed s/"#VERSION#"/$(VERSION)/g > $@)
					

.SUFFIXES:
.SUFFIXES: .tex .bbl .dvi .pdf .bib .html .pdf .png .eps
	
figs/%.png:
	$(if $(IMAGE_URL),,$(error IMAGE_URL is empty))
	$(CURL) C$(IMAGE_PATH)/`basename $@` -o $@
	$(EXTRACTBB) $@

