#!/bin/bash

latex lvdb_table.tex
bibtex lvdb_table
latex lvdb_table.tex
latex lvdb_table.tex
pdflatex lvdb_table.tex