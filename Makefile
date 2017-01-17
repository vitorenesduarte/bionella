.PHONY: site pdf

site:
	python build-site.py

pdf:
	jupyter nbconvert --to latex index.ipynb
	pdflatex index.tex
	rm index.aux index.log index.out index.tex index.toc
	google-chrome index.pdf

first:
	python first.py

second:
	python second.py

third:
	python third.py

start-nb:
	jupyter notebook

