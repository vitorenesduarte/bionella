#!/bin/bash
cd site/

for PAGE in $(ls *.ipynb); do
  jupyter nbconvert --to html --template full $PAGE
done

google-chrome index.html

