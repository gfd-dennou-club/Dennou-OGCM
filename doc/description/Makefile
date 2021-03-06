###########################################################
#
# Makefile for compiling tex document
#
###########################################################

# target file name
TARGET = dcpom_description
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
IMAGE_PATH = $(IMAGE_URL)/DCPOM_v$(VERSION)/description

IMAGE_OGCM_BASENAME = $(shell grep includegraphics ogcm/*tex | ruby -n -e '/^\s*\\includegraphics.*\{(.+)\}/ =~ $$_.split(":")[1] && print(/\A([^.]+)(\..+)?\Z/=~File.basename($$1) && $$1+($$2 || ".eps"),"\n")')
IMAGE_SICE_BASENAME = $(shell grep includegraphics sice/*tex | ruby -n -e '/^\s*\\includegraphics.*\{(.+)\}/ =~ $$_.split(":")[1] && print(/\A([^.]+)(\..+)?\Z/=~File.basename($$1) && $$1+($$2 || ".eps"),"\n")')
IMAGES = $(patsubst %,ogcm/figs/%,$(IMAGE_OGCM_BASENAME)) \
				 $(patsubst %,sice/figs/%,$(IMAGE_SICE_BASENAME))


SRCS = ${TARGET}.tex \
	introduction.tex              \
	ogcm/OGCM_description.tex     \
	ogcm/model_basic.tex	        \
	ogcm/space_discretization.tex \
	ogcm/time_discretization.tex  \
	ogcm/OGCM_config_example.tex  \
	sice/SICE_description.tex     \
	sice/model_basic.tex	        \
	sice/space_discretization.tex \
	sice/time_discretization.tex  \
	sice/SICE_config_example.tex  \
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
	
ogcm/figs/%.png:
	$(if $(IMAGE_URL),,$(error IMAGE_URL is empty))
	$(CURL) C$(IMAGE_PATH)/ogcm/`basename $@` -o $@
	$(EXTRACTBB) $@

sice/figs/%.png:
	$(if $(IMAGE_URL),,$(error IMAGE_URL is empty))
	$(CURL) $(IMAGE_PATH)/sice/`basename $@` -o $@
	$(EXTRACTBB) $@
