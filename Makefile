site:
	jupyter nbconvert --to html --template full index.ipynb
	google-chrome index.html

pdf:
	jupyter nbconvert --to latex index.ipynb
	pdflatex index.tex
	rm index.aux index.log index.out index.tex index.toc
	google-chrome index.pdf

