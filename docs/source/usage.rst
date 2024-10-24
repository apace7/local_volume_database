Usage & Decription of tables 
============================

For LVDB users that only want to  use  the combined catalogs/tables, installing the LVDB package is not required and the tables can be downloaded from the release page.

.. _installation:

Installation (note that this doesn't work yet)
------------

To use local_volume_database, first install it using pip:

.. code-block:: console

   (.venv) $ pip install local_volume_database

The `LVDBDIR` envirnment variable is used to point to the location of the input YAML files (/data_input/). 


Database content
----------------

The database is structured as individual yaml files for each system and combined tables as csv and fits files (descriptions of both below). 
The yaml files are located in `data_input/ <https://github.com/apace7/local_volume_database/tree/main/data_input>`_ and the combined csv tables are located in `data/ <https://github.com/apace7/local_volume_database/tree/main/data>`_. Fits files are located in the release pages.



The tables can be directly loaded into jupyter notebooks without having to download the repository from either the release page (recommend) or from the github:

.. code-block:: python

  import astropy.table as table
  ## release page, the version will need to be updated to the latest release
  dsph_all = table.Table.read('https://github.com/apace7/local_volume_database/releases/download/v0.0.2/dwarf_all.csv')
  ## latest github
  dsph_all = table.Table.read('https://raw.githubusercontent.com/apace7/local_volume_database/main/data/dwarf_all.csv')

There is also a pdf document in release page summarizing the contents and properties of each combined table (`pdf summary file <https://github.com/apace7/local_volume_database/releases/download/v0.0.2/lvdb_table.pdf>`_ Note to verify that you are using the latest release here). 


Decription of Catalogs/Tables 
--------------------

The following are the available tables (in csv and fits file formats). The fits file format is limited to github releases while the csv is include both in the release and the main catalog.

* **dwarf_mw** : Milky Way dwarf galaxies (the most distant dwarf galaxy is Eridanus II at ~ 350 kpc)
* **dwarf_m31**: M31 dwarf galaxies
* **dwarf_local_field**: dwarf galaxies outside of MW/M31 within the Local group to a distance of ~ 3 Mpc. This is an extension of galaxies from McConnachie 2012
* **dwarf_all** : combination of dwarf_mw, dwarf_m31, dwarf_local_field. Complete for known systems to ~ 3 Mpc.
* **dwarf_local_field_distant**: dwarf galaxies with distance > 3 Mpc. The limiting distance is set to ~10-40 Mpc (the approximate limits of HST/JWST). This table is not complete. 

* **gc_halo**: New star cluster-like systems in the Galactic halo or ambiguous compact stellar systems. These are faint star-cluster like systems (generally rhalf < 20 pc and M_V > -3 and at high Galactic latitudes abs(b) > ~5-10). A number of these systems are likely tidally stripped star clusters, tidally stripped dwarf galaxies, or the faintest dwarf galaxies. Many have are ambiguous classifications (as in the classification for dwarf galaxy and star cluster is ambiguous) and are difficult to classify. This does include several new Globular Cluster (Laevens 1/Crater I and Sagittarius II).
* **gc_disk**: post-Harris catalog globular clusters at low Galactic latitude (abs(b) <10), some of these systems might be open clusters, and some systems have not been confirmed
* **gc_harris**: globular clusters in Harris catalog (this excludes Koposov 1 and 2 which are in the gc_ufsc table)
* **gc_dwarf_hosted**: Globular clusters hosted by dwarf galaxies. This does not include the Sagittarius GCs which are in gc_harris. Incomplete catalog.

 (Note that older versions had  small differences between dwarf galaxy and star cluster catalogs)

There are two extra tables (data/pm_overview.csv and data/j_factor.csv). Both are collections of measurements (the other tables have one measurement per system). 

pm_overview.csv: key, reference, proper motion measurement, method (this includes most proper motion measurements of dwarf galaxies and the goal is to be complete for literature measurements).

j_factor.csv: key, reference, angle, j-factor measurement [units are log10 GeV^2 cm^-5], notes (this includes some literature j-factor measurements, mostly from A. B. Pace.  This is not complete.).

.. Decription of table contents
.. ----------------------------

Columns:

* key: unique identifier for each system.  The yaml input files have the same name.
* host: host of system [MW, LMC, M31, etc]
* confirmed_real: system has been confirmed with either deeper photometry, follow-up spectroscopy, proper motion, or other methods
* confirmed_dwarf: (or confirmed_star_cluster) system has been confirmed to be dwarf galaxy (or star cluster) based on spectroscopy, and/or deeper photometry.
* ra: right ascension ICRS [degree]
* dec: declination ICRS [degree]
* rhalf: major axis half-light radius (or plummer radius) in [arcmin]. Note that input yaml files can have arcsec or arcmin input units. 
* ellipticity: 1 - minor/major axis (or 1 - axis ratio)
* position_angle: N->E [degree] 
* distance_modulus [mag]
* apparent_magnitude_v: apparent magnitude in V-band. Corrected for extinction. Value added.
* vlos_systemic: heliocentric velocity of system [km/s]
* vlos_sigma: velocity dispersion in line-of-sight [km/s]
* metallicity_spectroscopic: spectroscopic metallicity [dex]
* metallicity_spectroscopic_sigma: spectroscopic metallicity dispersion [dex]
* pmra: systemic proper motion in right ascension, includes cos(dec) term following Gaia [mas/yr]
* pmdec: systemic proper motion in declination direction [mas/yr]
* rcore, rking: profile fits with king profile in arcmin
* rad_sersic, n_sersic: sersic profile parameters. rad_sersic in arcmin
* age: age of system [Gyr] 
* metallicity_photometric: metallicity from isochrone fitting (or non-spectroscopic metallicity such as metallicity sensitive narrowband imaging)
* flux_HI: flux in HI [Jy km s^−1]
* ref_ + x : reference columns such as ref_structure, ref_distance, ref_m_v, ref_vlos, ref_proper_motion.  All reference columns have the same format: author last name + ADS bibcode. 

Value-Added Columns:

* M_V: absolute V-band magnitude, computed from distance_modulus and apparent_magnitude_V
* mass_stellar: log10 stellar mass assuming M/L=2 and computed from M_V [log10 Msun]
* distance: heliocentric distance, computed from distance_modulus [kpc]
* ll: Galactic longitude
* bb: Galactic latitude
* sg_xx: Supergalactic coordinates, x [kpc]
* sg_yy: Supergalactic coordinates, y [kpc]
* sg_zz: Supergalactic coordinates, z [kpc] 
* distance_gc: 3D distance to Galactic center [kpc]
* distance_m31: 3D distance to M31 center [kpc]
* distance_lg: 3D distance to Local Group center [kpc] 
* distance_host: 3D distance to host galaxy [kpc]
* mass_HI: log10 HI mass computed from flux_HI and distance [log10 Msun] 
* metallicity: union of spectroscopic and photometric metallicity, spectroscopic preferred over photometric metallicity [dex]
* metallicity_type: lists whether `metallicity` column is photometric or spectroscopic. 
* velocity_gsr: Velocity in Galactic standard of rest frame [km/s]
* velocity_lg: Velocity of system relative to the Local Group centroid [km/s]
* mass_dynamical_wolf: Dynamical mass within 3D half-light radius using the dynamical mass estimator in `Wolf et al. 2010 <https://ui.adsabs.harvard.edu/abs/2010MNRAS.406.1220W/abstract>`_ [log10 Msun]. This column has errors and upper limit columns (em, ep, ul) using the errors from the half-light radius (rhalf), ellipticity, distance, and velocity dispersion (monte carlo errors). 
* rhalf_physical: half-light radius in physical units --  rhalf * distance  [parsec]. Includes monte carlo errors (distance and rhalf errors).
* rhalf_sph_physical: spherically averaged half-light radius (geometric mean); rhalf * distance * sqrt(1-ellipticity) in [parsec]. Includes monte carlo errors (distance, ellipticity, and rhalf errors).
* surface_brightness_rhalf: average surface brightness within spherically averaged half-light radius [mag arcsec^-2]
* ref_ + x : reference columns such as ref_structure, ref_distance, ref_m_v, ref_vlos, ref_proper_motion.  All reference columns have the same format: author last name + ADS bibcode. 

Many columns also have associated error columns. These follow the format of name + _em, + _ep + _ul (e.g., rhalf_em).

Error Columns: 

* _em = error minus = minus 1 sigma (or 16% confidence/credible interval) 
* _ep = error plus = plus 1 sigma (84% confidence/credible interval)
* _ll = lower limt at  5% confidence/credible interval 
* _ul = upper limit at 95% confidence/credible interval (some are at 5sigma, 90% or 84%, but the goal is to make it consistent)

The format for the reference columns is author last name + ADS bibcode. The author's last name has special characters removed but the capitalization is unchanged. 
There is an associated bibtex file (latex/lvdb.bib) that includes all references in the database. 

Decription of yaml files 
------------------------

There is an `example yaml file <https://github.com/apace7/local_volume_database/blob/main/code/example_yaml.yaml>`_ in the /code/ folder. 
It includes all collections and keys in the database.  Not all keys are included in the combined csv tables.
The yaml keys are **Bolded** below and the bullet points follow the yaml collection structure.  Errors columns are not included in the list below and some columns include upper limits in the combined table. 
The collections are split such that a single reference can describe the contents.

* **key** —- unique internal identifier. This should be the same as the name of the file (without .yaml) (required yaml key). All keys are lowercase in LVDB. Globular clusters and some dwarf galaxies are grouped by their host (for example, all LMC globular cluster keys have the prefix lmc_gc_ and many Centuarus A dwarf galaxy keys have the prefix cena_ ). 
* **table** -- table to place system into (required yaml key) list of possible tables [gc_harris, gc_dwarf_hosted, gc_disk, gc_halo=gc_ufsc, dwarf_mw , dwarf_local_field , dwarf_m31 , dwarf_local_field_distant, candidate, misc]. Systems in the candidate and misc tables are not combined into files. The candidate systems are included in the lvdb pdf summary. The misc systems are primarily bright host galaxies (MW, M31, Cen A) and only included for distance measurements (**distance_fixed_host**) and to link systems together. 
* **location** -- center of the system (yaml collection)

  * **ra** -- right ascension ICRS [degree]  (required yaml key)

  * **dec** -- declination ICRS [degree] (required yaml key)

* **name_discovery**

  * **name** -- name of system

  * **other_name** -- list of additional names of the system

  * **ref_discovery** --- List of discovery references. There can be multiple discovery references due to independent discoveries made on similar    timescales. Follow-up confirmation studies (i.e. HST imaging for distant candidate dwarfs around local volume hosts). Re-discoveries of systems (i.e. globular clusters hosted by dwarf galaxies).

  * **discovery_year** -- year of discovery. This may follow the arxiv year instead of the journal publication year.

  * **host** -- host of system.

  * **confirmed_dwarf** -- 0/1 -- 1 = confirmed dwarf galaxy.  

  * **confirmed_star_cluster** -- 0 or 1 -- 1 = confirmed star cluster.  

  * **confirmed_real** -- 1 = system is confirmed to be physical system.  This includes deeper imaging (i.e. HST), spectroscopic confirmation, and/or proper motion confirmation.

  * **false_positive** -- 1 = system is confirmed to be a false positive.  2 = system is confirmed to be background galaxy at much larger distances

  * **ref_false_positive** -- list of references that shows an system is a false positive. This could include new dwarf galaxy searches that do not recover the system. This includes dwarf galaxies candidates that are later shown to be background galaxies. 

  * **abbreviation** -- Common abbreviation for system (currently only for MW dwarf galaxies). 
  
  * **type** -- dSph, dIrr, NSC=Nuclear star cluster, GC=Globular Cluster

* **structure** -- yaml collection
  
  * **rhalf** -- elliptical half-light radius (or plummer radius) [arcmin]. This corresponds to the major axis. Default units is arcmin if arcsec the **spatial_units** key needs to be set. 

  * **spatial_units** -- options = [arcmin, arcsec] sets the units for the input radial parameter.

  * **spatial_model** -- options = [plummer, exponential, sersic, king, eff] model assumption for the primary model assumed to compute rhalf.  Included for reference.

  * **ellipticity** -- Ellipticity of the system, defined as 1 - b/a = 1- minor axis/major axis. 

  * **position_angle** -- position angle defined north to east [degree]

  * **ref_structure** -- reference

* **distance** -- yaml collection

  * **distance_modulus** --  distance modulus of the system. [mag] This quantity is used to compute the distance in kpc for each system.

  * **distance_fixed_host** -- True/False. This option fixes the distance of the system to the distance of its host.  Commonly used for globular clusters hosted by dwarf galaxy and new (unconfirmed) satellites of more distant hosts (>3 Mpc)

  * **ref_distance**

* **m_v** -- yaml collection

  * **apparent_magnitude_v** -- Apparent V-band magnitude of the system. This quantity is corrected for extinction. This quantity is combined with **distance_modulus** to compute the absolute V-band magnitude in the combined tables. 

  * **mean_ebv** -- Mean E(B-V) for reference.  This is not currently used in calculations. 

  * **ref_m_v** -- reference

* **velocity** -- stellar velocity/kinematics

  * **vlos_systemic** -- systemic heliocentric velocity of the system. Stellar velocities are preferred but some distant systems only have HI velocities. [km/s]
  
  * **vlos_sigma** -- stellar velocity dispersion. [km/s]. Sometimes called the global velocity dispersion.

  * **vlos_sigma_central** -- central stellar velocity dispersion. [km/s]. Primarily for globular clusters.
  
  * **ref_vlos** -- reference

* **proper_motion**
  
  * **pmra** -- systemic proper motion in the direction of right ascension (includes cosdec term) [mas/yr]

  * **pmdec** -- systemic proper motion in the direction of declination [mas/yr]

  * **ref_proper_motion** -- reference

* **spectroscopic_metallicity**

  * **metallicity_spectroscopic** -- mean metallicity

  * **metallicity_spectroscopic_sigma** -- metallicity dispersion

  * **ref_metallicity_spectroscopic** -- reference

* **metallicity_photometric**

  * **metallicity_photometric** -- photometric metallicity. This can include isochrone fitting or narrow band photometry.

  * **metallicity_photometric_sigma** -- metallicity dispersion from photometric measurements. Many for narrow band photometry. 

  * **ref_metallicity_photometric** -- reference

* **structure_king**

  * **rcore** -- King core radius [arcmin]. Default units is arcmin if arcsec the **spatial_units** key needs to be set. 

  * **rking** -- King limiting radius, sometimes referred to as the tidal radius [arcmin]. Default units is arcmin if arcsec the **spatial_units** key needs to be set. 

  * **spatial_units** -- options = [arcmin, arcsec] sets the units for the input radial parameter.
  
  * **ellipticity** and **position_angle** -- these are specfic to the King profile fit 

  * **ref_structure_king** -- reference

* **structure_sersic**

  * **n_sersic** -- Sersic powerlaw value.

  * **rad_sersic** -- Sersic radius [arcmin]. Default units is arcmin if arcsec the **spatial_units** key needs to be set. 

  * **spatial_units** -- options = [arcmin, arcsec] sets the units for the input radial parameter.

  * **ellipticity** -- from Sersic fit.

  * **position_angle** -- from Sersic fit.

  * **central_surface_brightness** -- central surface brightness of Sersic fit [mag/arcsec^2]

  * **ref_structure_sersic**

* **structure_eff**

  * **gamma_eff** -- Powerlaw value from EFF profile (Elson, Fall & Freeman 1987).

  * **rad_eff** -- EFF scale radius [arcmin]. Default units is arcmin if arcsec the **spatial_units** key needs to be set. 

  * **spatial_units** -- options = [arcmin, arcsec] sets the units for the input radial parameter.

  * **ellipticity** -- from EFF fit.

  * **position_angle** -- from EFF fit.

  * **ref_structure_sersic**

* **structure_plummer**

  * **rplummer** -- Plummer scale radius [arcmin]. Default units is arcmin if arcsec the **spatial_units** key needs to be set. 

  * **spatial_units** -- options = [arcmin, arcsec] sets the units for the input radial parameter.

  * **ellipticity** -- from Plummer fit.

  * **position_angle** -- from Plummer fit.

  * **ref_structure_plummer**

* **structure_exponential**

  * **rexponential** -- Exponential scale radius [arcmin]. Default units is arcmin if arcsec the **spatial_units** key needs to be set. 

  * **spatial_units** -- options = [arcmin, arcsec] sets the units for the input radial parameter.

  * **ellipticity** -- from EExponentialFF fit.

  * **position_angle** -- from Exponential fit.

  * **ref_structure_exponential**

* **flux_HI**

  * **flux_HI** -- [Jy km/s]

  * **vlos_systemic_HI** -- HI systemic velocity [km/s]

  * **sigma_HI** -- velocity dispersion of HI gas [km/s]

  * **vrot_HI** -- rotation velocity of HI gas [km/s]

  * **ref_flux_HI**

* **age**
  
  * **age** -- mean age of the systemic in [Gyr]. Mainly for star clusters. 

  * **ref_age** -- reference

* **star_formation_history**
  
  * **tau_50** -- time for 50 per cent of stellar mass to form [Gyr ago]

  * **tau_80** -- time where 80 per cent of stellar mass has formed, quenching time [Gyr ago]

  * **tau_90** -- time where 90 per cent of stellar mass has formed, quenching time [Gyr ago]

  * **ref_star_formation_history**

Citations to database and citations to the LVDB input
-----------------------------

The LVDB is set up to enable citations to the analysis and papers that serves as input to the LVDB. All reference columns (**ref_**) follow the same format of author last name (removed of special characters) + `NASA ADS bibcode <https://ui.adsabs.harvard.edu/>`_. There is a BibTeX file (`table/lvdb.bib <https://github.com/apace7/local_volume_database/blob/main/table/lvdb.bib>`_) with BibTeX entries from ADS with the key matching the LVDB reference column. There is an `ADS public library <https://ui.adsabs.harvard.edu/public-libraries/fVKkEJbdRyCmscCOwzsz6w>`_ that contains many of the input papers to the LVDB (with the goal to eventually contain all papers in the LVDB).  Papers replaced in the future will not be removed. 
The example notebook `example_notebooks/example_latex_citations.ipynb/ <https://github.com/apace7/local_volume_database/blob/main/example_notebooks/example_latex_citations.ipynb>`_ contains an example of creating a latex table with citations using the LVDB. 

As ADS bibcode are a fixed length of 19 characters, the ADS bibcode can be retrieved from the LVDB reference columns.  Other public tools such as  `adstex <https://github.com/yymao/adstex>`_ can be used to create bibtex files. 

Users of the LVDB are encouraged to cite the LVDB input of the systems studied in their analysis to give proper acknowledgment to the community.  

If you use this in your research please include a link to the github repository (https://github.com/apace7/local_volume_database) and cite the database paper (once it is written). 
An example in latex is: This work has made use of the Local Volume Database\footnote{\url{https://github.com/apace7/local_volume_database }}.
