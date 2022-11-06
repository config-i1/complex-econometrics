#!/bin/sh

set -ev

# Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"
Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::pdf_book')"
# Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::epub_book')"

sleep 3

git commit -am "New version of the book compiled"

git push
#
# git add * || true
#
# git commit -am "New version of the book compiled" || true
#
# git push
