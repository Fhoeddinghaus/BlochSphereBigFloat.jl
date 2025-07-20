#!/bin/sh

cd docs
git commit -am "Save uncommited changes (WIP)"
julia make.jl
git branch --delete --force gh-pages
git checkout --orphan gh-pages
git reset

rm -rf ../src
rm ../Manifest.toml
rm ../Project.toml
rm ../README.md
cp -r build/ ../

cd ..
rm -rf docs

git add --all
git commit -m "Rebuild GitHub pages"
git push -f origin gh-pages && git checkout main