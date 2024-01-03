# Contributing Guide

The package is built as a set of yaml files (1 per object) that are combined into a series of tables. 
There are 6 tables created directly from the yaml file:
- dwarf_mw
- dwarf_m31
- dwarf_local_field
- gc_ufsc
- gc_disk
- gc_harris

A 7th table is created that is the combined table of `dwarf_mw` `dwarf_m31` and `dwarf_local_field`
The `table` entry in the yaml file sets which table the object goes into. 

Contributions of updates to yaml file, new yaml file (new objects), and code updates are welcome.
