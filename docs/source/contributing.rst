How to Contribute to the Local Volume Database
==============================================

Community contributions to the local_volume_database are welcome. Examples contributions include: adding new measurements, updating current systems, adding new systems (both newly discovered and more distant), and/or new properties of current systems. 

Example: Add or update a measurement
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

Example: How to add a new system the database
---------------------------------------------

This example walks through the progress to add a new system to the database. 
Make a new yaml file in the ``data_input`` folder with the name ``new_system.yaml`` where new_system is the key of the new system.  This same key needs to be included as a yaml key named ``key`` (see example_yaml.yaml files).
The easiest why to find the parameter input options and names to simply make a copy of the example yaml file  ``code/example_yaml.yaml`` with the new_system name.  Note that  ``example_yaml.yaml`` contains all possible properties/collections and keys for the database.  For the new system many of these will not be measured and should be deleted or commented out. 
There are 4 required entries for each system: ``key``, ``table``, ``ra``, ``dec``.  The latter two keys are in the ``location`` collection. The ``table`` key is used to combine systems into tables in the ``data/`` folder. 

The new system will be added to the database tables by running the ``code/combine_table_general.py`` python script.

The new system will be added to the summary pdf tables by running the ``code/create_latex_table.py`` python script and running the latex scripts.

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
* misc (note that this key is used to collect objects that are not included in other tables)

.. How the database is constructed
.. ---------------------------------------------



Potential Contributions/Expansion
---------------------------------------

As stated earlier, community contributions are welcome and encouraged.  
Here is a short list of items that generally focus on expanding the scope of the database.  

* Include gas kinematic properties. For example, peak rotation velocity and gas velocity dispersion.
* Statistics on RRL or other variable/rare stars in dwarf galaxies
* Star formation history information.  For example, a quenching timescale could be included
* Other star formation history tracers, FUV (GALEX, SWIFT etc), Halpha
* kinematic information for globular clusters (average velocity dispersion and central velocity disperion)
* open clusters
* LMC/SMC/M31 clusters
* Complete dwarf galaxy entries for objects beyond 3 Mpc
* Complete dwarf galaxy globular cluster objects and properties

