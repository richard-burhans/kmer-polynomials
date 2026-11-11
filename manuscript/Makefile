# Build the k-mer polynomial manuscript: Markdown source, PDF output.
#   make pdf      paper.md + references.bib -> paper.pdf (and paper.tex, kept for linting and inspection)
#   make lint     writing-repo style linter on paper.md, and the reference validator on references.bib
#   make readme   copy the abstract from paper.md into ../README.md (between ABSTRACT markers)
#   make check    fail if ../README.md's abstract differs from paper.md's
# Pandoc runs in a pinned image, so the build does not depend on a local TeX install.
PANDOC_IMAGE ?= pandoc/latex:3.5
WRITING ?= /home/rcb112/project/claude/bench/repos/writing
DOCKER = docker run --rm -u $(shell id -u):$(shell id -g) -v $(CURDIR):/data -w /data $(PANDOC_IMAGE)
PANDOC_FLAGS = --citeproc --bibliography=references.bib --pdf-engine=xelatex -V geometry:margin=1in --number-sections

.PHONY: pdf tex lint readme check clean
pdf: paper.pdf
tex: paper.tex

paper.tex: paper.md references.bib
	$(DOCKER) $(PANDOC_FLAGS) -s paper.md -o paper.tex

paper.pdf: paper.md references.bib
	$(DOCKER) $(PANDOC_FLAGS) paper.md -o paper.pdf

lint:
	python3 $(WRITING)/scripts/lint_manuscript_style.py paper.md
	python3 $(WRITING)/scripts/validate_references.py references.bib

readme:
	python3 sync_abstract.py

check:
	python3 sync_abstract.py --check

clean:
	rm -f paper.pdf paper.tex
