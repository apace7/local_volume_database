Decription of yaml files 
===================================

There is an `example yaml file <https://github.com/apace7/local_volume_database/tree/main/data/>`_ in the /code/ folder. 
It includes all collections and keys in the database with a short descrition and units.  Not all keys are placed into the csv tables.
The yaml keys are **Bolded** below.

* **key** —- unique internal identifier. This should be the same as the name of the file (without .yaml) (required yaml key)
* **table** -- table to place system into (required yaml key)
* **location** -- center of the system (yaml collection)

  #. **ra** -- right ascension ICRS [degree]  (required yaml key)

  #. **dec** -- declination ICRS [degree] (required yaml key)

* **name_discovery**

  #. **name** -- name of system

  #. **other_name** -- list of additional names of the system

  #. **discovery_year** -- year of discovery. This may follow the arxiv year instead of the journal publication year.

  #. **host** -- host of system.