import numpy as np
import matplotlib.pyplot as plt

import astropy.table as table

from astropy import units as u

from collections import Counter

import numpy.ma as ma
import yaml

dsph_mw = table.Table.read('data/dwarf_mw.csv')
items_to_join = dsph_mw['key']

## this is not the right way to do this
## this needs to be changed such that the table is made with masked columns
# np.ma.masked_all(len(items_to_join), dtype=type_for_column)
comb_table = table.Table(np.zeros((len(items_to_join), 3)),names=('key', 'ra', 'dec' ), dtype=('U100','f8', 'f8', ))
col_name = ['name','host','confirmed_real', 'confirmed_dwarf',  'rhalf', 'rhalf_em', 'rhalf_ep', 'position_angle', 'position_angle_em', 'position_angle_ep', 'ellipticity', 'ellipticity_em', 'ellipticity_ep', 'ellipticity_ul', 'ref_structure', 'distance_modulus', 'distance_modulus_em', 'distance_modulus_ep', 'ref_distance', 'apparent_magnitude_v', 'apparent_magnitude_v_em', 'apparent_magnitude_v_ep', 'ref_m_v', 'vlos_systemic', 'vlos_systemic_em', 'vlos_systemic_ep', 'vlos_sigma', 'vlos_sigma_em', 'vlos_sigma_ep', 'vlos_sigma_ul', 'ref_vlos', 'pmra', 'pmra_em', 'pmra_ep', 'pmdec', 'pmdec_em', 'pmdec_ep', 'ref_proper_motion', 'metallicity_spec', 'metallicity_spec_em', 'metallicity_spec_ep', 'metallicity_spec_sigma', 'metallicity_spec_sigma_em', 'metallicity_spec_sigma_ep', 'metallicity_spec_sigma_ul', 'ref_metallicity_spec', 'rcore', 'rcore_em', 'rcore_ep', 'rking', 'rking_em', 'rking_ep', 'ref_king', 'rad_sersic', 'rad_sersic_em', 'rad_sersic_ep', 'n_sersic', 'n_sersic_em', 'n_sersic_ep', 'ref_sersic', 'age', 'age_em', 'age_ep', 'ref_age', 'metallicity_photometric', 'metallicity_photometric_em', 'metallicity_photometric_ep', 'ref_metallicity_photometric', 'flux_HI', 'flux_HI_em', 'flux_HI_ep', 'flux_HI_ul', 'ref_flux_HI']
col_type = [ 'U100','U100','i1', 'i1','f8', 'f8', 'f8','f8', 'f8', 'f8','f8', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8' ,'f8', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'f8','U100']

for i,j in zip(col_name, col_type):
    comb_table[i] = np.ma.masked_all(len(items_to_join), dtype=j)

missing_key = []
for i in range(len(items_to_join)):
    with open("data_input/" + items_to_join[i] + '.yaml', 'r') as stream:
        try:
            stream_yaml = yaml.load(stream, Loader=yaml.Loader)
            comb_table['key'][i] = items_to_join[i]
            comb_table['ra'][i] = stream_yaml['location']['ra']
            comb_table['dec'][i] = stream_yaml['location']['dec']
            
            for y in ['structure', 'distance','m_v', 'velocity', 'proper_motion', 'spec_metallicity', 'structure_king', 'structure_sersic', 'age', 'metallicity_photometric', 'flux_HI', 'name_discovery']:
                if y in stream_yaml.keys():
                    x = list(stream_yaml[y].keys())
                    for list_key in range(len(x)):
                        name = x[list_key]

                        if name not in comb_table.dtype.names:
                            missing_key.append(name)

                            continue
                        else:
                            comb_table[name][i] = stream_yaml[y][x[list_key]]

        except yaml.YAMLError as exc:
            print(exc)

print(Counter(missing_key))

def dist_mod(mu, mu_em=0, mu_ep=0):
    def dm(x):
        return pow(10., x/5.+1.)/1000.
    return dm(mu), dm(mu)-dm(mu-mu_em), dm(mu+mu_ep)-dm(mu)

## value added columns
comb_table['M_V'] = comb_table['apparent_magnitude_v']-comb_table['distance_modulus']
comb_table['M_V_em'] = comb_table['apparent_magnitude_v_em']
comb_table['M_V_ep'] = comb_table['apparent_magnitude_v_ep']

d, dem, dep = dist_mod( comb_table['distance_modulus'], comb_table['distance_modulus_em'], comb_table['distance_modulus_ep'])
comb_table['distance'] = d
comb_table['distance_em'] = dem
comb_table['distance_ep'] = dep

comb_table['rhalf_physical'] = comb_table['distance']*1000.*comb_table['rhalf']/60./180.*np.pi
comb_table['rhalf_physical_sph'] = comb_table['rhalf_physical']*np.sqrt(1.-comb_table['ellipticity'])

comb_table['mass_HI'] = np.log10(235600 * comb_table['flux_HI']*(comb_table['distance']/1000.)**2 )

comb_table['metallicity'] = np.ma.masked_all(len(comb_table), dtype=float)
comb_table['metallicity_em'] = np.ma.masked_all(len(comb_table), dtype=float)
comb_table['metallicity_ep'] = np.ma.masked_all(len(comb_table), dtype=float)
comb_table['metallicity_type'] = np.ma.masked_all(len(comb_table), dtype='U100')
for i in range(len(comb_table)):
    if ma.is_masked(comb_table['spec_metallicity'][i])==False:
        comb_table['metallicity'][i] = comb_table['spec_metallicity'][i]
        comb_table['metallicity_em'][i] = comb_table['spec_metallicity_em'][i]
        comb_table['metallicity_ep'][i] = comb_table['spec_metallicity_ep'][i]
        comb_table['metallicity_type'][i] = 'spectroscopic'
    elif  ma.is_masked(comb_table['metallicity_photometric'][i])==False:
        comb_table['metallicity'][i] = comb_table['metallicity_photometric'][i]
        comb_table['metallicity_em'][i] = comb_table['metallicity_photometric_em'][i]
        comb_table['metallicity_ep'][i] = comb_table['metallicity_photometric_ep'][i]
        comb_table['metallicity_type'][i] = 'photometric'
    

comb_table.write('data_output/dwarf_mw.csv', overwrite=True)
