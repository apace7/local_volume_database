Decription of yaml files 
===================================

There is an `example yaml file <https://github.com/apace7/local_volume_database/blob/main/code/example_yaml.yaml>`_ in the /code/ folder. 
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
  
  * **rhalf** -- [arcmin]. Default units is arcmin if arcsec the **spatial_units** key needs to be set. 

  * **spatial_units** -- options = [arcmin, arcsec] sets the units for the input rhalf

  * **ellipticity**

  * **position_angle**

  * **ref_structure**

* **distance** -- yaml collection

  * **distance_modulus**

  * **distance_fixed_host** -- True/False. This option fixes the distance of the object to the distance of its host.  Commonly used for globular clusters hosted by dwarf galaxy and new (unconfirmed) satellites of more distant hosts (>3 Mpc)

  * **ref_distance**

* **m_v** -- yaml collection

  * **apparent_magnitude_v** -- corrected for extinction

  * **mean_ebv** -- Mean E(B-V) for reference.  This is not currently used in calculations. 

  * **ref_m_v**

* **velocity** -- stellar velocity/kinematics

  * **vlos_systemic** -- systemic heliocentric velocity of the system. Stellar velocities are preferred but some distant objects are from HI observations. [km/s]
  
  * **vlos_sigma** -- stellar velocity dispersion. [km/s]
  
  * **vlos_sigma_central** -- central stellar velocity dispersion. [km/s]

  * **ref_vlos**

* **proper_motion**
  
  * **pmra** [mas/yr]

  * **pmdec** [mas/yr]

  * **ref_proper_motion** 

* **spectroscopic_metallicity**

  * **metallicity_spectroscopic**

  * **metallicity_spectroscopic_sigma** -- metallicity dispersion

  * **ref_metallicity_spectroscopic**

* **metallicity_photometric**

  * **metallicity_photometric**

  * **metallicity_photometric_sigma**

  * **ref_metallicity_photometric**

* **structure_king**

  * **rcore** -- King core radius. [arcmin]. Default units is arcmin if arcsec the **spatial_units** key needs to be set. 

  * **rking** -- King limiting radius, commonly called the tidal radius. [arcmin]. Default units is arcmin if arcsec the **spatial_units** key needs to be set. 

  * **spatial_units** -- options = [arcmin, arcsec] sets the units for the input radial parameters.

  * **ellipticity** -- from King fit.

  * **position_angle** -- from King fit.

  * **ref_structure_king**

* **structure_sersic**

  * **n_sersic** -- Sersic powerlaw value.

  * **rad_sersic** -- Sersic radius [arcmin]. Default units is arcmin if arcsec the **spatial_units** key needs to be set. 

  * **spatial_units** -- options = [arcmin, arcsec] sets the units for the input radial parameter.

  * **ellipticity** -- from Sersic fit.

  * **position_angle** -- from Sersic fit.

  * **ref_structure_sersic**

* **structure_eff**

  * **gamma_eff** -- Powerlaw value from EFF profile (Elson, Fall & Freeman 1987).

  * **rad_eff** -- EFF scale radius [arcmin]. Default units is arcmin if arcsec the **spatial_units** key needs to be set. 

  * **spatial_units** -- options = [arcmin, arcsec] sets the units for the input radial parameter.

  * **ellipticity** -- from Sersic fit.

  * **position_angle** -- from Sersic fit.

  * **ref_structure_sersic**

* **flux_HI**

  * **flux_HI** -- HI flux in [km/s]

  **vlos_systemic_HI** -- Hi systemic velocity [km/s]

  **sigma_HI** -- velocity dispersion of HI gas [km/s]

  **vrot_HI** -- rotation velocity of HI gas [km/s]

* **age**

  * **age** -- [Gyr]

  * **ref_age**