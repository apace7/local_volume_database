Decription of tables 
===================================

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
* rhalf: elliptical half-light radius (or plummer radius) in [arcmin]
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
