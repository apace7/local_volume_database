#!/bin/zsh

source ~/.zshrc

python scripts/combine_table_general.py

print 'creating latex tables'

python scripts/create_latex_table.py

read -s -k '?Press any key to create summary plots.'

python scripts/create_summary_plots.py

read -s -k '?Press any key to compile latex and create pdf.'

cd table/

latex lvdb_table.tex
bibtex lvdb_table
latex lvdb_table.tex
latex lvdb_table.tex
pdflatex lvdb_table.tex