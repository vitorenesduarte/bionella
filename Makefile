.PHONY: site pdf

site:
	build-site.sh

pdf:
	jupyter nbconvert --to latex index.ipynb
	pdflatex index.tex
	rm index.aux index.log index.out index.tex index.toc
	google-chrome index.pdf

