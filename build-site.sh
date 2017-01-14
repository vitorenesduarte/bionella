#!/bin/bash
jupyter nbconvert --to html --template full index.ipynb

cd site/
for PAGE in $(ls *.ipynb); do
  jupyter nbconvert --to html --template full $PAGE
done
cd ../

google-chrome index.html
