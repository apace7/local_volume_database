Tutorials \& Examples
=====================

`ipython notebook example <https://github.com/apace7/local_volume_database/blob/main/example_notebooks/example_plots.ipynb>`_ 


Example: How to add a new system the database
---------------------------------------------

Make a new yaml file for the system with the name ``key.yaml`` where key matches the ``key`` key in the yaml file in the ``data_input`` folder.  
(in general I make a copy of the code/example_yaml.yaml).  There are 4 required entries for each system: ``key``, ``table``, ``ra``, ``dec``.  The latter two keys are in the ``location`` collection. The ``table`` key is used to combine systems into the various tables. 

The new system can be added to the database tables by running the ``code/combine_table_general.py`` python script.

The new system can be added to the summary pdf tables by running the ``code/create_latex_table.py`` python script and running the latex scripts.

There is a bash/zsh shell script that will recreate the database and latex/pdf summary tables (there is a good chance that this script only works on my computer).

These are the current tables: 

* dwarf_mw
* dwarf_local_field
* dwarf_m31
* gc_harris
* gc_ufsc
* gc_disk
* gc_dwarf_hosted
* candidate
* local_field_distant
* misc

