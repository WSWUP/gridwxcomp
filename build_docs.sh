#!/usr/bin/env bash

# after making and pushing changes
git checkout gh-pages
git rm -rf .
git checkout master docs gridwxcomp CHANGE_LOG.rst # get ALL docs source files
cd docs
make clean && make html # make docs from source
cd ..
mv ./docs/build/html/* ./
rm -rf docs/ gridwxcomp/ CHANGE_LOG.rst # remove docs source files
touch .nojekyll
git add -A
git commit -m "publishing updated docs"

# switch back
git checkout master


