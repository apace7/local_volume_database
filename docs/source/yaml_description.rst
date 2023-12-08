Decription of yaml files 
===================================

There is an `example yaml file <https://github.com/apace7/local_volume_database/tree/main/data/>`_ in the /code/ folder. 
It includes all collections and keys in the database with a short descrition and units.  Not all keys are placed into the csv tables.
The yaml keys are **Bolded** below and the bullet points follow the yaml collection structure.  Errors columns are not included. 
The collections are split such that a single reference can describe the contents.

* **key** —- unique internal identifier. This should be the same as the name of the file (without .yaml) (required yaml key)
* **table** -- table to place system into (required yaml key)
* **location** -- center of the system (yaml collection)

  * **ra** -- right ascension ICRS [degree]  (required yaml key)

  * **dec** -- declination ICRS [degree] (required yaml key)

* **name_discovery**

  * **name** -- name of system

  * **other_name** -- list of additional names of the system

  * **ref_discovery** --- List of discovery references. There can be multiple discovery references due to independent discoveries made on similar    timescales. Follow-up confirmation studies (i.e. HST imaging for distant candidate dwarfs around local volume hosts). Re-discoveries of systems (i.e. globular clusters hosted by dwarf galaxies).

  * **discovery_year** -- year of discovery. This may follow the arxiv year instead of the journal publication year.

  * **host** -- host of system.

  * **confirmed_dwarf** -- 0/1 1 = confirmed dwarf galaxy.  

  * **confirmed_star_cluster** -- 0 or 1. 1 = confirmed star cluster.  

  * **confirmed_real** --

  * **false_positive** -- 

* **structure** -- yaml collection
  
  * **rhalf** -- [arcmin] 

  * **ellipticity**

  * **position_angle**

  * **ref_structure**

* **distance** -- yaml collection

  * **distance_modulus**

  * **ref_distance**

* **m_v** -- yaml collection

  * **apparent_magnitude_v** -- corrected for extinction

  * **ref_m_v**

* **velocity** -- stellar velocity/kinematics

  * **vlos_systemic** -- systemic heliocentric velocity of the system. Stellar velocities are preferred but some distant objects are from HI observations. [km/s]
  
  * **vlos_sigma** -- stellar velocity dispersion. [km/s]
  
  * **ref_vlos**

* **proper_motion**
  
  * **pmra** [mas/yr]

  * **pmdec** [mas/yr]

  * **ref_proper_motion** 