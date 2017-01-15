.PHONY: site pdf

start-nb:
	jupyter notebook

site:
	python build-site.py

pdf:
	jupyter nbconvert --to latex index.ipynb
	pdflatex index.tex
	rm index.aux index.log index.out index.tex index.toc
	google-chrome index.pdf

