
Below is the example_yaml file for reference

- key: example_yaml
- table: misc (dwarf_mw, dwarf_m31, dwarf_local_field, gc_harris, gc_disk, gc_ufsc)
- location: ra, dec 
-- ra: 0
-- dec: 0

distance:
  distance_modulus: 
  
  ref_distance: ref
m_v:
  apparent_magnitude_v: 0
  apparent_magnitude_v_em: 0
  apparent_magnitude_v_ep: 0
  ref_m_v: ref
name_discovery:
  confirmed_dwarf: 0
  confirmed_star_cluster: 0
  confirmed_real: 0
  discovery_year: 0
  host: misc
  name: name
  other_name:
  - list
  ref_discovery:
  - list
notes:
  - list
proper_motion:
  pmra: 0
  pmra_em: 0
  pmra_ep: 0
  pmdec: 0
  pmdec_em: 0
  pmdec_ep: 0
  ref_proper_motion: ref
metallicity_spectroscopic:
  metallicity_spectroscopic: 0
  metallicity_spectroscopic_em: 0
  metallicity_spectroscopic_ep: 0
  metallicity_spectroscopic_sigma: 0
  metallicity_spectroscopic_sigma_em: 0
  metallicity_spectroscopic_sigma_ep: 0
  metallicity_spectroscopic_sigma_ul: 0
  ref_metallicity_spectroscopic: ref
metallicity_photometric:
  metallicity_photometric: 0
  metallicity_photometric_em: 0
  metallicity_photometric_ep: 0
  ref_metallicity_photometric: ref

structure:
## the ellipticity, position_angle here are the version in the standard table
  rhalf: 0
  rhalf_em: 0
  rhalf_ep: 0
  ellipticity: 0
  ellipticity_em: 0
  ellipticity_ep: 0
  ellipticity_ul: 0
  position_angle: 0
  position_angle_em: 0
  position_angle_ep: 0
  ref_structure: ref

structure_king:
## the ellipticity, position_angle here are not propagated to the standard table
  rcore: 0
  rcore_em: 0
  rcore_ep: 0
  rking: 0
  rking_em: 0
  rking_ep: 0
  ellipticity: 0
  ellipticity_em: 0
  ellipticity_ep: 0
  ellipticity_ul: 0
  position_angle: 0
  position_angle_em: 0
  position_angle_ep: 0
  ref_structure_king: ref

structure_sersic:
## the ellipticity, position_angle here are not propagated to the standard table
  n_sersic: 0
  n_sersic_em: 0
  n_sersic_ep: 0
  rad_sersic: 0
  rad_sersic_em: 0
  rad_sersic_ep: 0
  ellipticity: 0
  ellipticity_em: 0
  ellipticity_ep: 0
  ellipticity_ul: 0
  position_angle: 0
  position_angle_em: 0
  position_angle_ep: 0
  ref_structure_sersic: ref

velocity:
  vlos_systemic: 0
  vlos_systemic_em: 0
  vlos_systemic_ep: 0
  vlos_sigma: 0
  vlos_sigma_em: 0
  vlos_sigma_ep: 0
  ref_vlos: ref 

flux_HI:
  flux_HI:  0
  flux_HI_em: 0
  flux_HI_ep: 0
  flux_HI_ul: 0
  flux_HI_rotation_velocity: # not in standard table
  flux_HI_rotation_velocity_em: 
  flux_HI_rotation_velocity_ep: 
  flux_HI_dispersion: # not in standard table
  flux_HI_dispersion_em: 
  flux_HI_dispersion_ep: 
  ref_flux_HI: ref

age:
  age: 0
  age_em: 0
  age_ep: 0
  ref_age: ref
