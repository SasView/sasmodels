# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
PYTHON       ?= python
SPHINXOPTS    =
SPHINXBUILD   = "$(PYTHON)" -m sphinx
PAPER        ?= letter
BUILDDIR      = _build
ifdef ComSpec
    RMDIR   = rmdir /s/q
    CPDIR   = xcopy /q/e/i
    MKDIR   = mkdir
    PATHSEP = $(strip \)
    COPY    = copy
else
    RMDIR   = rm -rf
    CPDIR   = cp -rp
    MKDIR   = mkdir -p
    PATHSEP = $(strip /)
    COPY    = cp
endif

# Internal variables.
PAPEROPT_a4     = -D latex_elements.papersize=a4
PAPEROPT_letter = -D latex_elements.papersize=letter
ALLSPHINXOPTS   = -d _build/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

IMG_SOURCE = $(wildcard ../sasmodels/models/img/*)
IMG_DEST = $(patsubst ../sasmodels/models/img/%,model/img/%,$(IMG_SOURCE))
MODELS_PY = $(wildcard ../sasmodels/models/[a-z]*.py)
#MODELS_RST = $(patsubst ../sasmodels/models/%.py,model/%.rst,$(MODELS_PY))

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html      to make standalone HTML files"
	@echo "  dirhtml   to make HTML files named index.html in directories"
	@echo "  pickle    to make pickle files"
	@echo "  json      to make JSON files"
	@echo "  htmlhelp  to make HTML files and a HTML help project"
	@echo "  qthelp    to make HTML files and a qthelp project"
	@echo "  latex     to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  changes   to make an overview of all changed/added/deprecated items"
	@echo "  linkcheck to check all external links for integrity"
	@echo "  doctest   to run all doctests embedded in the documentation (if enabled)"

model/img/%: ../sasmodels/models/img/% model/img
ifdef ComSpec
	$(COPY) $(subst /,\,$<) $(subst /,\,$@)
else
	$(COPY) $< $@
endif

ref/models/index.rst: gentoc.py genmodel.py $(MODELS_PY)
	"$(PYTHON)" genmodel.py $(MODELS_PY)
	"$(PYTHON)" gentoc.py $(MODELS_PY)

.PHONY: help clean html dirhtml pickle json htmlhelp qthelp latex changes linkcheck doctest build

api: genapi.py
	-$(RMDIR) api
	"$(PYTHON)" genapi.py

model:
	$(MKDIR) model

model/img: model
ifdef ComSpec
	$(MKDIR) model\img
else
	$(MKDIR) model/img
endif

build: model $(IMG_DEST) api ref/models/index.rst
ifdef ComSpec
	rem cd ..\..
	rem python setup.py build
else
	#cd ../.. && python setup.py build
endif

upload:
	scp -r _build/html/* reflectometry.org:web/danse/docs/sasmodels
	scp _build/latex/SASModels.pdf reflectometry.org:web/danse/download

clean:
	-$(RMDIR) _build model api guide$(PATHSEP)models

html: build
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) _build/html
	@echo
	@echo "Build finished. The HTML pages are in _build/html."

dirhtml:
	$(SPHINXBUILD) -b dirhtml $(ALLSPHINXOPTS) _build/dirhtml
	@echo
	@echo "Build finished. The HTML pages are in _build/dirhtml."

pickle:
	$(SPHINXBUILD) -b pickle $(ALLSPHINXOPTS) _build/pickle
	@echo
	@echo "Build finished; now you can process the pickle files."

json:
	$(SPHINXBUILD) -b json $(ALLSPHINXOPTS) _build/json
	@echo
	@echo "Build finished; now you can process the JSON files."

htmlhelp:
	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) _build/htmlhelp
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in _build/htmlhelp."

qthelp:
	$(SPHINXBUILD) -b qthelp $(ALLSPHINXOPTS) _build/qthelp
	@echo
	@echo "Build finished; now you can run "qcollectiongenerator" with the" \
	      ".qhcp project file in _build/qthelp, like this:"
	@echo "# qcollectiongenerator _build/qthelp/Periodic Table.qhcp"
	@echo "To view the help file:"
	@echo "# assistant -collectionFile _build/qthelp/Periodic Table.qhc"

latex: build
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) _build/latex
	@echo
	@echo "Build finished; the LaTeX files are in _build/latex."
	@echo "Run \`make all-pdf' or \`make all-ps' in that directory to" \
	      "run these through (pdf)latex."

pdf: latex
	$(MAKE) -C _build/latex all-pdf
	$(COPY) _build/latex/SASModels.pdf _build/html

changes:
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) _build/changes
	@echo
	@echo "The overview file is in _build/changes."

linkcheck:
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) _build/linkcheck
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in _build/linkcheck/output.txt."

doctest:
	$(SPHINXBUILD) -b doctest $(ALLSPHINXOPTS) _build/doctest
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in _build/doctest/output.txt."
