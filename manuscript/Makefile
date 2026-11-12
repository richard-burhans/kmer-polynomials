# Build the k-mer polynomial manuscript: Markdown source, PDF output.
#   make figures  figures/make_figures.py -> figures/*.svg and *.pdf, from committed results files
#   make pdf      paper.md + references.bib + figures -> paper.pdf (and paper.tex, kept for linting and inspection)
#   make lint     writing-repo style linter on paper.md, and the reference validator on references.bib
#   make readme   copy the abstract from paper.md into ../README.md (between ABSTRACT markers)
#   make check    fail if ../README.md's abstract differs from paper.md's, or if the committed paper.pdf is stale
# Pandoc runs in a pinned image, so the build does not depend on a local TeX install.
PANDOC_IMAGE ?= pandoc/latex:3.5
WRITING ?= /home/rcb112/project/claude/bench/repos/writing
DOCKER = docker run --rm -u $(shell id -u):$(shell id -g) -v $(CURDIR):/data -w /data $(PANDOC_IMAGE)
PANDOC_FLAGS = --citeproc --bibliography=references.bib --pdf-engine=xelatex -V geometry:margin=1in --number-sections

FIGURES = $(addprefix figures/,$(addsuffix .pdf,fig1_encoding fig2_collapse fig3_speed fig4_sketch))
# matplotlib is pinned and fetched by uv, so the figures do not depend on a local install. The output is deterministic.
FIGURE_PY = SOURCE_DATE_EPOCH=0 nice -n 19 ionice -c3 uv run -q --no-project --with matplotlib==3.9.2 python

.PHONY: pdf tex figures lint readme check clean
pdf: paper.pdf
tex: paper.tex

paper.tex: paper.md references.bib $(FIGURES)
	$(DOCKER) $(PANDOC_FLAGS) -s paper.md -o paper.tex

# paper.pdf is COMMITTED so the draft can be read on GitHub. paper.pdf.sources records the sha256 of the files it was
# built from, and `make check` fails when the committed PDF is behind paper.md, references.bib or a figure.
SOURCES = paper.md references.bib $(FIGURES)
paper.pdf: $(SOURCES)
	$(DOCKER) $(PANDOC_FLAGS) paper.md -o paper.pdf
	sha256sum $(SOURCES) > paper.pdf.sources

figures:
	$(FIGURE_PY) figures/make_figures.py

$(FIGURES): figures/make_figures.py paper_bench_table.py
	$(FIGURE_PY) figures/make_figures.py

lint:
	python3 $(WRITING)/scripts/lint_manuscript_style.py paper.md
	python3 $(WRITING)/scripts/validate_references.py references.bib

readme:
	python3 sync_abstract.py

check:
	python3 sync_abstract.py --check
	@sha256sum --quiet -c paper.pdf.sources || { echo "paper.pdf is out of date: run make pdf and commit it"; exit 1; }
	@echo "paper.pdf matches its sources"

clean:
	rm -f paper.pdf paper.tex
