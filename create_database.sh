#!/bin/zsh

source ~/.zshrc

python code/combine_table_general.py

python code/create_latex_table.py

read -s -k '?Press any key to continue.'

cd table/

latex lvdb_table.tex
bibtex lvdb_table
latex lvdb_table.tex
latex lvdb_table.tex
pdflatex lvdb_table.tex