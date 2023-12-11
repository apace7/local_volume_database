Usage & Decription of tables 
============================

.. _installation:

Installation (note that this doesn't work yet)
------------

To use local_volume_database, first install it using pip:

.. code-block:: console

   (.venv) $ pip install local_volume_database

Database content
----------------

The database is structured as individual yaml files for each system and combined tables as csv and fits files. 
The yaml files are located in data_input/ and the combined tables in data/. 

The tables can be directly loaded into jupyter notebooks without having to download the repository:

.. code-block:: python

   import astropy.table as table
   dsph_mw = table.Table.read('https://raw.githubusercontent.com/apace7/local_volume_database/main/data/dwarf_mw.csv')


Decription of tables 
--------------------

table descriptions (tab/sheet). These are available as csv and fits files. 

* dwarf_mw : Milky Way dwarf galaxies
* dwarf_m31: M31 dwarf galaxies
* dwarf_local_field: dwarf galaxies outside of MW/M31 to ~ 3 Mpc, mostly follows McConnachie 2012
* gc_ufsc: post Harris catalog star clusters in the MW halo (abs(b) > ~5-10)
* gc_disk: post-Harris catalog globular clusters at low Galactic latitude (abs(b) <10), some of these objects might be open clusters
* gc_harris: globular clusters in Harris catalog


pm_overview: key, reference, proper motion measurement, method (this includes most proper motion measurements of dwarf galaxies)

j_factor.csv: key, reference, angle, j-factor measurement [units are log10 GeV^2 cm^-5], notes (this includes some literature j-factor measurements, mostly from A. B. Pace)

Decription of table contents
----------------------------

columns:

* key: unique identifier for each system.  The yaml input files have the same name.
* host: host of system [MW, LMC, M31, etc]
* confirmed_real: system has been confirmed with either deeper photometry, follow-up spectroscopy, proper motion, or other methods
* confirmed_dwarf: (or confirmed_star_cluster) system has been confirmed to be dwarf galaxy (or star cluster) based on spectroscopy, and/or deeper photometry.
* ra: right ascension ICRS [degree]
* dec: declination ICRS [degree]
* rhalf: half-light radius (or plummer radius) in [arcmin]
* ellipticity: 1 - minor/major axis (or 1 - axis ratio)
* position_angle: N->E [degree] 
* distance_modulus [mag]
* distance: computed from distance_modulus [kpc] 
* rhalf_physical: rhalf * distance  [parsec] (computed from other columns)
* rhalf_sph_physical: rhalf * distance * sqrt(1-ellipticity) in [parsec] (computed from other columns)
* apparent_magnitude_v: apparent magnitude in V-band. Corrected for extinction. 
* M_V: absolute V-band magnitude, computed from distance_modulus and apparent_magnitude_V
* surface_brightness_rhalf: average surface brightness within spherically averaged half-light radius [mag arcsec^-2]
* vlos_systemic: heliocentric velocity of system [km/s]
* vlos_sigma: velocity dispersion in line-of-sight [km/s]
* metallicity: metallicity, spectroscopic preferred [dex]
* metallicity_type: lists whether `metallicity` column is photometric or spectroscopic
* metallicity_spectroscopic: spectroscopic metallicity [dex]
* metallicity_spectroscopic_sigma: spectroscopic metallicity dispersion [dex]
* pmra: proper motion in right ascension, includes cos(dec) term following Gaia [mas/yr]
* pmdec: proper motion in declination direction [mas/yr]
* rcore, rking: profile fits with king profile in arcmin
* rad_sersic, n_sersic: sersic profile parameters. rad_sersic in arcmin
* age: Gyr units
* metallicity_photometric: metallicity isochrone fitting (or non-spectroscopic metallicity)
* flux_HI: ( Jy km s^−1 ), flux in HI only included for dwarf galaxies
* ref: reference columns (ref_structure, ref_distance, ref_m_v, ref_vlos, ref_proper_motion) of author last name + ADS bibcode

error columns: 
* _em = error minus = minus 1 sigma (or 16% confidence interval) 
* _ep = error plus = plus 1 sigma (84% confidence interval)
* _ul = upper limit at 95% confidence interval (some are at 5sigma, 90% or 84%, but the goal is to make it consistent)

ref: reference of author last name + ADS bibcode

Decription of yaml files 
------------------------

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

* **spectroscopic_metallicity**

  * **metallicity_spectroscopic**

  * **metallicity_spectroscopic_sigma** -- metallicity dispersion

  * **ref_metallicity_spectroscopic**

* **metallicity_photometric**

  * **metallicity_photometric**

  * **metallicity_photometric_sigma**

  * **ref_metallicity_photometric**

* **structure_king**
* **structure_sersic**
* **structure_eff**
* **flux_HI**
* **age**