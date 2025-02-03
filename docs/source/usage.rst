Usage & Decription of tables 
============================

For LVDB users that only want to  use  the combined catalogs/tables, installing the LVDB package is not required and the primary catalogs are located as attachments in the `GitHub release pages <https://github.com/apace7/local_volume_database/releases>`_.
The catalogs are included as csv and fits files and there is  a pdf summary file of the LVDB content.
For reproducibility, it is recommended to use one of the tagged release versions for scientific analysis.



.. _installation:

Installation 
------------

The LVDB package is installable locally:

.. code-block:: console

  git clone https://github.com/apace7/local_volume_database.git
  cd local_volume_database
  python -m build
  pip install .

The ``LVDBDIR`` envirnment variable is used to point to the location of the input YAML files (local_volume_database/data_input/). 
Pip package coming. As the package interacts with the YAML files, it is recommended to install locally.

Database Content
----------------

The database is structured as a collection of YAML files, where the properties of each system are located in an individual YAML file.
The YAML files are combined into catalogs that are available as csv and fits files on the release page. 
The YAML files are located in `data_input/ <https://github.com/apace7/local_volume_database/tree/main/data_input>`_ and the combined csv tables are located in `data/ <https://github.com/apace7/local_volume_database/tree/main/data>`_ or on the `release page <https://github.com/apace7/local_volume_database/releases>`_. Fits files are only located in the release pages.



The tables can be directly loaded into jupyter notebooks without having to download the repository.
An example of loading the tables remotely from either the release page (recommend) or  github is as follows:

.. code-block:: python

  import astropy.table as table
  ## release page, the version will need to be updated to the latest release
  dsph_all = table.Table.read('https://github.com/apace7/local_volume_database/releases/download/v1.0.0/dwarf_all.csv')
  ## latest github
  dsph_all = table.Table.read('https://raw.githubusercontent.com/apace7/local_volume_database/main/data/dwarf_all.csv')

Note that the version number will need to be changed to access the lasted release page while the github link will go to the version on the main branch.
There is also a pdf document (named lvdb_table.pdf) in release page summarizing the contents and properties of each combined table. 


Decription of Catalogs/Tables 
--------------------

The following are the available tables (in csv and fits file formats). The fits file is limited to the release pages while the csv file is included in both the release and main github.

* **dwarf_mw** : Milky Way dwarf galaxies (the most distant dwarf galaxy is Eridanus II at ~ 350 kpc).
* **dwarf_m31**: M31 dwarf galaxies.
* **dwarf_local_field**: dwarf galaxies outside of MW/M31 within the Local Field to a distance of ~ 3 Mpc. This is an extension of galaxies from McConnachie 2012 compilation.
* **dwarf_local_field_distant**: dwarf galaxies in the Local Volume with distance > 3 Mpc. The limiting distance is set to ~10-40 Mpc (the approximate limits of HST/JWST). This table is not complete to known systems (it is complete for known systems to a distance < 3.5 Mpc). 
* **dwarf_all** : combination of dwarf_mw, dwarf_m31, dwarf_local_field, dwarf_local_field_distant. Complete for known systems to ~ 3.5 Mpc. Note that earlier versions did not include dwarf_local_field_distant. 
* **gc_ambiguous**: systems with an ambiguous classification (referred to as ambiguous or hyper-faint compact stellar systems in the LVDB). These are all MW halo systems. 
* **gc_mw_new**: newly discovered globular clusters or candidate globular clusters (i.e. post-Harris catalog).  Many systems are at low Galactic latitude (abs(b) <10-20 deg) and candidate systems may be open clusters.
* **gc_harris**: globular clusters in Harris catalog (this excludes Koposov 1 and 2 which are in the gc_abmiguous table).
* **gc_dwarf_hosted**: globular clusters hosted by dwarf galaxies. This does not include the Sagittarius globular clusters which are in gc_harris. This catalog is incomplete for known systems.
* **gc_other**: for other globular clusters. (mostly for future work)
* **candidate**: known false-positive candidates, background galaxies, or low confidence candidates. **Only included in the release page.**
* **misc_host**: brighter galaxies that are hosts to the dwarf galaxies.  The catalog exists for completeness and for host information for dwarf galaxies.  The main properties compiled for these systems are phase-space information (ra,dec,distance, velocity) and overall stellar and gas mass. **Only included in the release page.**
* **comb_all**: the union of all tables.  Includes all systems in the LVDB.  This table has an additional column `table` that specifies the table origin of the system. **Only included in the release page.**





.. Decription of table contents
.. ----------------------------

Columns:

* key: unique identifier for each system.  The yaml input files have the same name.
* host: host [LVDB key] of system [MW, LMC, M31, etc]
* confirmed_real: system has been confirmed with either deeper photometry, follow-up spectroscopy, proper motion, or other methods (not a chance alignment of stars).
* confirmed_dwarf: (or confirmed_star_cluster) system has been confirmed to be dwarf galaxy (or star cluster) based on spectroscopy, and/or deeper photometry.
* ra: right ascension ICRS J2000.0 [degree]
* dec: declination ICRS J2000.0 [degree]
* rhalf: major axis of the half-light radius (or plummer radius) in [arcmin]. Note that input yaml files can have arcsec or arcmin input units but the combined catalogs are in arcmin. 
* ellipticity: 1 - minor/major axis (or 1 - axis ratio).
* position_angle: N->E [degree] 
* distance_modulus [mag]
* apparent_magnitude_v: apparent magnitude in V-band. Corrected for extinction. 
* vlos_systemic: heliocentric velocity of system [km/s]
* vlos_sigma: velocity dispersion in line-of-sight [km/s]
* metallicity_spectroscopic: spectroscopic metallicity [dex]
* metallicity_spectroscopic_sigma: spectroscopic metallicity dispersion [dex]
* pmra: systemic proper motion in right ascension, includes cos(dec) term following Gaia [mas/yr]
* pmdec: systemic proper motion in declination direction [mas/yr]
* rcore, rking: profile fits with king profile in arcmin
* rad_sersic, n_sersic: sersic profile parameters. rad_sersic in arcmin
* age: age of system [Gyr] 
* metallicity_isochrone: metallicity from isochrone or cmd fitting 
* flux_HI: flux in HI [Jy km s^−1]
* ref_ + x : reference columns such as ref_structure, ref_distance, ref_m_v, ref_vlos, ref_proper_motion.  All reference columns have the same format: author last name + ADS bibcode. 

Value-Added Columns:

* M_V: absolute V-band magnitude, computed from distance_modulus and apparent_magnitude_V
* mass_stellar: log10 stellar mass assuming M/L=2 and computed from M_V [log10 Msun]
* distance: heliocentric distance, computed from the distance_modulus column [kpc]
* ll: Galactic longitude [degree]
* bb: Galactic latitude [degree]
* sg_xx: Supergalactic coordinates, x [kpc]
* sg_yy: Supergalactic coordinates, y [kpc]
* sg_zz: Supergalactic coordinates, z [kpc] 
* distance_gc: 3D distance to Galactic center [kpc]
* distance_m31: 3D distance to M31 center [kpc]
* distance_lg: 3D distance to Local Group center [kpc] 
* distance_host: 3D distance to host galaxy [kpc]
* mass_HI: log10 HI mass computed from flux_HI and distance [log10 Msun] 
* metallicity: union of spectroscopic, photometric, and isochrone, spectroscopic preferred over photometric metallicity, and photometric over isochrone [dex]
* metallicity_type: lists whether `metallicity` column is photometric, isochrone or spectroscopic. 
* velocity_gsr: Velocity in Galactic standard of rest frame [km/s]
* velocity_lg: Velocity of system relative to the Local Group centroid [km/s]
* mass_dynamical_wolf: Dynamical mass within 3D half-light radius using the dynamical mass estimator in `Wolf et al. 2010 <https://ui.adsabs.harvard.edu/abs/2010MNRAS.406.1220W/abstract>`_ [log10 Msun]. This column has errors and upper limit columns (em, ep, ul) using the errors from the half-light radius (rhalf), ellipticity, distance, and velocity dispersion (monte carlo errors). 
* rhalf_physical: half-light radius in physical units --  rhalf * distance  [parsec]. Includes monte carlo errors (distance and rhalf errors).
* rhalf_sph_physical: azimuthally-averaged half-light radius (geometric mean); rhalf * distance * sqrt(1-ellipticity) in [parsec]. Includes monte carlo errors (distance, ellipticity, and rhalf errors).
* surface_brightness_rhalf: average surface brightness within azimuthally-averaged half-light radius [mag arcsec^-2]
* ref_ + x : reference columns such as ref_structure, ref_distance, ref_m_v, ref_vlos, ref_proper_motion.  All reference columns have the same format: author last name + ADS bibcode. 

Many columns also have associated error columns. These follow the format of name + _em, + _ep + _ul (e.g., rhalf_em). 

Error Columns: 

* _em = error minus = minus 1 sigma (or 16% confidence/credible interval) 
* _ep = error plus = plus 1 sigma (84% confidence/credible interval)
* _ll = lower limt at  5% confidence/credible interval 
* _ul = upper limit at 95% confidence/credible interval (some are at 5sigma, 90% or 84%, but the goal is to make it consistent)

The format for the reference columns is author last name + ADS bibcode. The author's last name has special characters removed but the capitalization is unchanged. 
There is an associated bibtex file (latex/lvdb.bib) that includes all references in the database. 

There are two extra tables: data/pm_overview.csv and data/j_factor.csv. The former is a compilation of systemic proper motion measurements for dwarf galaxies in the Local Group and the latter is a collection of J-factor measurements. Both are collections of measurements (the other tables have one measurement per system). 
The pm_overview table includes most proper motion literature measurements of dwarf galaxies and HFCSS. 
The j-factor table includes some literature j-factor measurements, mostly from A. B. Pace.  This is not complete for literature measurements.

pm_overview.csv: LVDB key, LVDB reference, ADS bibcode, proper motion measurements (full columns = pmra, pmra_em, pmra_ep, pmdec, pmdec_em, pmdec_ep, correlation) [the units are mas/yr expect for the unitless corrleation column], method [current options include=GAIA_EDR3, GAIA_DR2, Ground, HST_Ground, HST, GAIA_DR2_HST, HSC, SRG, GAIA_EDR3_HST, maser, GaiaHub], text citation, comments


j_factor.csv: LVDB key, LVDB reference, ADS bibcode, text citation, seleciton, angle [degree], j-factor measurement [units are log10 GeV^2 cm^-5] (full column names = logj, logj_em,	logj_ep,	logj_em05,	logj_ep95,	logj_ul95), use, comments

Decription of YAML Files 
------------------------

There is an `example yaml file <https://github.com/apace7/local_volume_database/blob/main/code/example_yaml.yaml>`_ in the /code/ folder. 
The example yaml file includes all collections and keys in the database.  Not all keys are included in the combined csv tables.
The yaml collections and keys are **Bolded**  and the bullet points follow the yaml collection structure.  Errors columns are not included in the list  and some columns include upper limits in the combined table. 
The collections are split such that a single reference can describe the contents.
Most keys are single entries and several keys are lists (specially other_name, ref_discovery, ref_false_positive).  

* **key** —- unique internal LVDB identifier (required yaml key). This should be the same as the name of the file (without .yaml).  All keys are lowercase in LVDB. Globular clusters and some dwarf galaxies are grouped by their host. For example, all LMC globular cluster keys have the prefix lmc_gc_ and many Centuarus A dwarf galaxy keys have the prefix cena_. Most new satellite systems will have a host prefix.
* **table** -- the table to place system into (required yaml key). The list of possible tables is: gc_harris, gc_dwarf_hosted, gc_disk=gc_mw_new, gc_halo=gc_ufsc=gc_abmiguous, dwarf_mw , dwarf_local_field , dwarf_m31 , dwarf_local_field_distant, candidate, misc, gc_other (there are several options that will place systems into the same table). The candidate and misc catalogs are only included in the release pages. The candidate systems are included in the lvdb pdf summary while the hosts/misc are not. The misc systems are primarily bright host galaxies (MW, M31, Cen A) and are partly included for distance measurements (**distance_fixed_host**) and to link systems together. 
* **location** -- yaml collection. center of the system 

  * **ra** -- right ascension ICRS [degree]  (required yaml key)

  * **dec** -- declination ICRS [degree] (required yaml key)

  * **ref_location** -- reference for center/location. Errors are supported for the center of the system. 

* **name_discovery** -- yaml collection

  * **name** -- name of system

  * **other_name** -- list of additional names of the system

  * **ref_discovery** --- List of discovery references. There can be multiple discovery references due to independent discoveries made on similar    timescales. Follow-up confirmation studies (i.e. HST imaging for distant candidate dwarfs around local volume hosts). Re-discoveries of systems (i.e. globular clusters hosted by dwarf galaxies).

  * **discovery_year** -- year of discovery. The year may be before the journal publication year due to an earlier arxiv submission.

  * **host** -- host galaxy of the system.

  * **confirmed_dwarf** -- Integer that denotes whether the system is confirmed to be a dwarf galaxy (options = 0,1). 1 = confirmed dwarf galaxy.  

  * **confirmed_star_cluster** -- Integer that denotes whether the system is confirmed to be a star cluster (options = 0,1).  1 = confirmed star cluster.  

  * **confirmed_real** -- Integer that denotes whether the system is confirmed  to be physical system (options = 0,1). 1 = confirmed system.  To confirm a system, deeper imaging (i.e. HST), spectroscopy, and/or proper motion/astrometry may be required. 

  * **false_positive** -- Integer that denotes whether the system is confirmed to a false positive or backkground galaxy (options = 0,1,2). 1 = system is confirmed to be a false positive.  2 = system is confirmed to be background galaxy at much larger distance (outside the Local Volume).

  * **ref_false_positive** -- list of references that shows an system is a false positive. This could include new dwarf galaxy searches that do not recover the system. This includes dwarf galaxies candidates that are later shown to be background galaxies. 

  * **abbreviation** -- Common abbreviation for system (currently only for MW dwarf galaxies). 
  
  * **type** -- Morphological type. This includes: dSph, dIrr, NSC=Nuclear star cluster, GC=Globular Cluster (this is not the full set of options). This key is generally incomplete.

  * **nme_lvg** -- exact name in the Catalog and Atlas of Local Volume galaxies (`LVG <https://www.sao.ru/lv/lvgdb/>`_). To enable a join on the LVG identifiers.

* **notes** -- List of notes in LaTeX. The notes are added to the summary pdf. 

* **structure** -- yaml collection
  
  * **rhalf** -- elliptical half-light radius [arcmin]. This corresponds to the major axis. The default units are arcmin if the **spatial_units** key is not included. 

  * **spatial_units** -- this key sets the units of the spatial parameter (rhalf here). The options are [arcmin, arcsec].

  * **spatial_model** -- options = [plummer, exponential, sersic, king, eff] model assumption for the primary model assumed to compute rhalf.  Included for reference.

  * **ellipticity** -- Ellipticity of the system, defined as 1 - b/a = 1- minor axis/major axis. 

  * **position_angle** -- position angle defined north to east [degree]

  * **diameter_holmberg** -- Holmberg isophote: projected major axis of galaxy at the isophotal level 25 mas/arcsec^2 in the B-band. Mainly included for systems without a half-light measurements (larger or brighter galaxies).

  * **ref_structure** -- reference

* **distance** -- yaml collection

  * **distance_modulus** --  distance modulus of the system. [mag] This quantity is used to compute the distance in kpc for each system.

  * **distance_fixed_host** -- True/False. This option fixes the distance of the system to the distance of its host.  Commonly used for globular clusters hosted by dwarf galaxy, systems without an independent distance measurement, and/or new candidate satellites in more distant systems (>3 Mpc).

  **distance_measurement_method** -- Refers to the method used for the distance measurement ['host', 'trgb', 'cmd', 'hb', 'rrl', 'sbf', 'nam']. 'hb' = horizontal branch, 'host' = distance fixed to the host (overlaps with **distance_fixed_host**), 'trgb' = tip of the red giant branch distance, 'sbf' = surface brightness fluctuation, 'rrl' = RR Lyrae, 'cmd' = color-magnitude diagram fitting, 'nam' = numerical action method based distance, 'btf' = baryonic Tully-Fisher distance, 'tf' = Tully-Fisher distance, 'sn' = supernova based distance

  * **ref_distance**

* **m_v** -- yaml collection

  * **apparent_magnitude_v** -- Apparent V-band magnitude of the system (Johnson-Kron-Cousins UBVRI photometric system). This quantity is corrected for extinction. This quantity is combined with **distance_modulus** to compute the absolute V-band magnitude in the combined tables. 

  * **apparent_magnitude_i** -- Apparent I-band magnitude of the system (Johnson-Kron-Cousins UBVRI photometric system). This quantity is corrected for extinction.

  * **apparent_magnitude_b** -- Apparent B-band magnitude of the system (Johnson-Kron-Cousins UBVRI photometric system). This quantity is corrected for extinction.

  * **mean_ebv** -- Mean E(B-V) for reference.  This is included for reference and is not used in calculations. 

  * **ref_m_v** -- Reference.

* **velocity** -- yaml collection. stellar velocity/kinematics

  * **vlos_systemic** -- systemic heliocentric velocity of the system. Stellar velocities are preferred but some distant systems only have HI velocities. [km/s]
  
  * **vlos_sigma** -- stellar velocity dispersion. [km/s]. Sometimes called the global velocity dispersion.

  * **vlos_sigma_central** -- central stellar velocity dispersion. [km/s]. Primarily for globular clusters.
  
  * **ref_vlos** -- reference

* **proper_motion** -- yaml collection
  
  * **pmra** -- systemic proper motion in the direction of right ascension (includes cosdec term) [mas/yr]

  * **pmdec** -- systemic proper motion in the direction of declination [mas/yr]

  * **pmra_pmdec_corr** -- correlation between pmra, pmdec, unitless [-1, 1]. 

  * **ref_proper_motion** -- reference

* **spectroscopic_metallicity** -- yaml collection

  * **metallicity_spectroscopic** -- mean metallicity

  * **metallicity_spectroscopic_sigma** -- metallicity dispersion

  * **ref_metallicity_spectroscopic** -- reference

* **metallicity_photometric** -- yaml collection

  * **metallicity_photometric** -- photometric metallicity. This generally is from metallicity sensistive photometry (Ca H&K, u-band). 

  * **metallicity_photometric_sigma** -- metallicity dispersion from photometric measurements. 

  * **ref_metallicity_photometric** -- reference

* **metallicity_isochrone** -- yaml collection

  * **metallicity_isochrone** -- metallicity from isochrone or color-magnitude diagram fitting. 

  * **metallicity_isochrone_sigma** -- metallicity dispersion from isochrone or color-magnitude diagram fitting. 

  * **ref_metallicity_isochrone** -- reference

* **structure_king** -- yaml collection

  * **rcore** -- King core radius [arcmin]. The default units are arcmin if the **spatial_units** key is not included. 

  * **rking** -- King limiting radius, sometimes referred to as the tidal radius [arcmin]. Default units is arcmin if arcsec the **spatial_units** key needs to be set. 

  * **spatial_units** -- this key sets the units of the spatial parameter. The options are [arcmin, arcsec].
  
  * **ellipticity** and **position_angle** -- these are specfic to the King profile fit 

  * **ref_structure_king** -- reference

* **structure_sersic** -- yaml collection

  * **n_sersic** -- Sersic powerlaw value.

  * **rad_sersic** -- Sersic radius [arcmin]. The default units are arcmin if the **spatial_units** key is not included. 

  * **spatial_units** -- this key sets the units of the spatial parameter. The options are [arcmin, arcsec].

  * **ellipticity** -- from Sersic fit.

  * **position_angle** -- from Sersic fit.

  * **central_surface_brightness** -- central surface brightness of Sersic fit [mag/arcsec^2]

  * **ref_structure_sersic**

* **structure_eff** -- yaml collection. EFF profile (Elson, Fall & Freeman 1987). Commonly used for globular clusters.

  * **gamma_eff** -- Powerlaw value from EFF profile (Elson, Fall & Freeman 1987).

  * **rad_eff** -- EFF scale radius [arcmin]. The default units are arcmin if the **spatial_units** key is not included. 

  * **spatial_units** -- this key sets the units of the spatial parameter. The options are [arcmin, arcsec].

  * **ellipticity** -- from EFF fit.

  * **position_angle** -- from EFF fit.

  * **ref_structure_sersic**

* **structure_plummer** -- yaml collection. 

  * **rplummer** -- Plummer scale radius [arcmin]. The default units are arcmin if the **spatial_units** key is not included. 

  * **spatial_units** -- this key sets the units of the spatial parameter. The options are [arcmin, arcsec].

  * **ellipticity** -- from Plummer fit.

  * **position_angle** -- from Plummer fit.

  * **ref_structure_plummer**

* **structure_exponential** -- yaml collection.

  * **rexponential** -- Exponential scale radius [arcmin]. The default units are arcmin if the **spatial_units** key is not included. 

  * **spatial_units** -- this key sets the units of the spatial parameter. The options are [arcmin, arcsec].

  * **ellipticity** -- from Exponential fit.

  * **position_angle** -- from Exponential fit.

  * **ref_structure_exponential**

* **flux_HI** -- yaml collection.

  * **flux_HI** -- [Jy km/s]

  * **vlos_systemic_HI** -- HI systemic velocity [km/s]

  * **sigma_HI** -- velocity dispersion of HI gas [km/s]

  * **vrot_HI** -- rotation velocity of HI gas [km/s]

  * **ref_flux_HI**

* **age** -- yaml collection.
  
  * **age** -- mean age of the systemic in [Gyr]. Mainly for star clusters. 

  * **ref_age** -- reference

* **star_formation_history** -- yaml collection. Mainly for dwarf galaxies.
  
  * **tau_50** -- time for 50 per cent of stellar mass to form [Gyr ago]

  * **tau_80** -- time for 80 per cent of stellar mass has formed, quenching time [Gyr ago]

  * **tau_90** -- time for 90 per cent of stellar mass has formed, quenching time [Gyr ago]

  * **ref_star_formation_history**

Citations to the LVDB and Citations to the LVDB Input
-----------------------------

The LVDB is set up to enable citations to the literature input of the LVDB. All reference columns (**ref_**) follow the same format of author last name (removed of special characters) + `NASA ADS bibcode <https://ui.adsabs.harvard.edu/>`_. There is a BibTeX file (`table/lvdb.bib <https://github.com/apace7/local_volume_database/blob/main/table/lvdb.bib>`_) with BibTeX entries from ADS with the key matching the LVDB reference column. There is an ADS public library (`Link <https://ui.adsabs.harvard.edu/public-libraries/fVKkEJbdRyCmscCOwzsz6w>`_) that contains the majority of the literature LVDB input.
The example notebook  contains an example of creating a latex table with citations using the LVDB (`example_notebooks/example_latex_citations.ipynb <https://github.com/apace7/local_volume_database/blob/main/example_notebooks/example_latex_citations.ipynb>`_). 
The LVDB package also contains a function that will output references (see `example_lvdb_package.ipynb <https://github.com/apace7/local_volume_database/blob/main/example_notebooks/example_lvdb_package.ipynb>`_).



As ADS bibcode are a fixed length of 19 characters, the ADS bibcode can be retrieved from the LVDB reference columns if users wish to use the ADS bibcodes instead.  Other public tools such as  `adstex <https://github.com/yymao/adstex>`_ package can be used to create bibtex files. 

Users of the LVDB are encouraged to cite the LVDB input (of the systems studied in their analysis) to give proper acknowledgment to the community.  The references could be included in a table or appendix. See Appendix A of this paper (`Cerny et al. 2024 <https://ui.adsabs.harvard.edu/abs/2024arXiv241000981C/abstract>`_) for an example of including internal LVDB references to the text of a paper.

If you use the LVDB in your research please include a link to the github repository (https://github.com/apace7/local_volume_database) and cite the LVDB overview paper (`Pace 2024 <https://ui.adsabs.harvard.edu/abs/2024arXiv241107424P/abstract>`_). 
An example in LaTeX that can be added to the acknowledgments section is: This work has made use of the Local Volume Database\footnote{\url{https://github.com/apace7/local_volume_database }} \citep{Pace2024arXiv241107424P}.

The bibtex of the LVDB paper is below:

.. code-block:: bibtex

  @ARTICLE{Pace2024arXiv241107424P,
    author = {{Pace}, Andrew B.},
        title = "{The Local Volume Database: a library of the observed properties of nearby dwarf galaxies and star clusters}",
    journal = {arXiv e-prints},
    keywords = {Astrophysics - Astrophysics of Galaxies},
        year = 2024,
        month = nov,
        eid = {arXiv:2411.07424},
        pages = {arXiv:2411.07424},
        doi = {10.48550/arXiv.2411.07424},
  archivePrefix = {arXiv},
    eprint = {2411.07424},
  primaryClass = {astro-ph.GA},
    adsurl = {https://ui.adsabs.harvard.edu/abs/2024arXiv241107424P},
    adsnote = {Provided by the SAO/NASA Astrophysics Data System}
  }

The LVDB releases are also indexed on `zenodo <https://doi.org/10.5281/zenodo.14076714>`_.

Link to the LVDB overview paper  on `arXiv <https://arxiv.org/abs/2411.07424>`_. and `ADS <https://ui.adsabs.harvard.edu/abs/2024arXiv241107424P/abstract>`_. 

.. The bibtex of the LVDB paper is below:

Extra Catalogs
-----------------------------

There are two additional catalogs included in the LVDB, pm_overview.csv and j_factor.csv. In constrast to other catalogs,  both these catalogs are compliations of measurements. pm_overview.csv compiles systemic proper motion measurements and j_factor.csv compiles J-factor  measurements (see appendix B of the LVDB overview paper for more details). 
The columns of the catalogs are described below.

pm_overview.csv  column description:

* key: LVDB key
* ref: ADS bibcode
* ref_cite: LVDB bibcode
* pmra: systemic proper motion [mas/yr] + (pmra_em and pmra_ep)
* pmdec: systemic proper motion [mas/yr] + (pmdec_em and pmdec_ep)
* pmra_pmdec_corr: correlation between errors [-1 to 1]
* method: options [GAIA_EDR3, GAIA_DR2, HST, GaiaHub, Ground, SRG, GAIA_EDR3_HST, HST_Ground, GAIA_DR2_HST, HSC, maser, Euclid + Gaia]
* citation: in text citation
* comments: notes

j_factor.csv column description: 

* key: LVDB key
* ref: ADS bibcode
* ref_cite: LVDB bibcode
* citation: in text citation
* selection: details on methodology
* angle: maximum angle [deg]
* logj: log_10 J-factor 
*	logj_em: 16% credible interval
*	logj_ep: 84% credible interval
* logj_em05: 5% credible interval
* logj_ep95: 95% credible interval
* logj_ul95: 95% upper limit
* use: value to use when there are multiple measurements in the same paper
* comments: noes