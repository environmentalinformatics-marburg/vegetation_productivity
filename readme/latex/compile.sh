#!/bin/bash

pdflatex readme.tex
bibtex readme.aux
bibtex readme.aux
pdflatex readme.tex
pdflatex readme.tex

