How to Contribute to the Local Volume Database
==============================================

Community contributions to the local_volume_database are welcome. Some potential contributions include: adding new measurements, updating current systems, adding new systems (both newly discovered and more distant), and/or new properties of current systems. 

Example: Add or Update a Measurement
------------------------------------

This example shows how to add or update new measurements to the database. 
Here we will look at a recent paper on new Magellan/IMACS spectroscopic data of two ultra-faint dwarf galaxies, Aquarius II and Bootes II (Bruce et al 2023). 
`Link to paper on ADS <https://ui.adsabs.harvard.edu/abs/2023ApJ...950..167B/abstract>`_  

Bruce et al (2023) made new and improved measurements of the systemic velocity, velocity dispersion, mean metallicity, and metallicity dispersion of both galaxies and these updated quantities are included in LVDB. The updated results for Bootes II are shown in the code block below. 

.. code-block:: yaml

    velocity:
        ref_vlos: Bruce2023ApJ...950..167B
        vlos_sigma: 2.9
        vlos_sigma_em: 1.2
        vlos_sigma_ep: 1.6
        vlos_systemic: -130.4
        vlos_systemic_em: 1.1
        vlos_systemic_ep: 1.4
    metallicity_spectroscopic:
        metallicity_spectroscopic: -2.71
        metallicity_spectroscopic_em: 0.1
        metallicity_spectroscopic_ep: 0.11
        metallicity_spectroscopic_sigma_ul: 0.37
        ref_metallicity_spectroscopic: Bruce2023ApJ...950..167B

For any measurement in the database a reference should be included. 
The format for the database is author last name + ADS bibcode ( ``Bruce2023ApJ...950..167B`` for this example). The author last name should be stripped of special characters and spaces but capitalization is not changed. 
If the reference is not in the ``table/lvdb.bib`` file you should add it. 
Create the reference with "Export Citation" on ADS and update the entry to match the reference format in the database. 
The content of the database can include both published and unpublished papers on the arXiv.  Papers only on arXiv still have an ADS bibcode created.  The reference in the database can be updated to the published version later. 

Example: How to Add a New System the Database
---------------------------------------------

This example walks through the progress to add a new system to the database. 
Make a new yaml file in the ``data_input`` folder with the name ``new_system.yaml`` where new_system is the key of the new system.  This same key needs to be included as a yaml key named ``key`` (see example_yaml.yaml files).
The easiest way to verify that the correct YAML key names are  used is to copy  the example yaml file  ``code/example_yaml.yaml`` when creating the new YAML file.  Note that  ``example_yaml.yaml`` contains all possible properties/collections and keys for the database.  For the new system many of the YAML parameters will not be measured and should be deleted. 
There are 4 required entries for each system: ``key``, ``table``, ``ra``, ``dec``.  The latter two keys are in the ``location`` collection. The ``table`` key is used to combine systems into tables in the ``data/`` folder.  The ``key`` value needs to be unique and match the file name.

List of current options for the tables: 

* dwarf_mw
* dwarf_local_field
* dwarf_m31
* dwarf_local_field_distant
* gc_ufsc = gc_halo = gc_ambiguous (all these table values go into the same catalog)
* gc_harris
* gc_disk = gc_mw_new (these options both end up in the sane table/catalog)
* gc_dwarf_hosted
* candidate
* misc 
* gc_other

The new system will be added to the catalogs by running the ``create_database.sh`` script. 

Table/Catalog Construction 
---------------------------------------------

The bash/zsh shell script, ``create_database.sh``, creates the tables/catalogs, overview figures, and the summary pdf. 
The script runs three python scripts, ``scripts/combine_table_general.py``, ``scripts/create_latex_table.py``, and ``scripts/create_summary_plots.py``, and runs the LaTeX compilation of the pdf after. 
The first script, combines together all the systems, adds value-added columns and saves the catalogs (csv and fits files).
The second script, uses the output and the YAML files to create the input LaTeX data and citations files for the ``lvdb_table.pdf``. 
The third script, creates summary figures. This includes recreating all the figures in the overview paper with the latest input files and some additional figures in a single combined pdf file, ``paper_examples/overview_plots.pdf``. The plots in the combined pdf file are for checking the content of the database (looking for outliers, etc). 




Some Ideas for Contributions 
---------------------------------------

As stated earlier, community contributions are welcome and encouraged.  
Here is a short list of items that generally focus on expanding the scope of the database.  
The github issues are another list of potential contributions.
Some of these items have YAML keys that exist but are generally empty.

* Include gas kinematic properties. For example, peak rotation velocity and gas velocity dispersion.
* Statistics on RRL or other variable/rare stars in dwarf galaxies.
* Star formation history information.  For example, a quenching timescale could be included.
* Other star formation history tracers, FUV (GALEX, SWIFT etc), Halpha.
* kinematic information for globular clusters (average velocity dispersion and central velocity dispersion).
* open clusters.
* LMC/SMC/M31 clusters.
* Complete dwarf galaxy entries for systems beyond 3 Mpc. 
* Complete dwarf galaxy globular cluster systems and properties. 

