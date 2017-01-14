site:
	jupyter nbconvert --to html --template full index.ipynb
	google-chrome index.html
