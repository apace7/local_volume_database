import numpy as np
import matplotlib.pyplot as plt

import astropy.table as table

from astropy import units as u

from collections import Counter

import numpy.ma as ma
import yaml

dsph_mw = table.Table.read('data/dwarf_mw.csv')
name_discovery = table.Table.read('data/name_discovery.csv')
notes = table.Table.read('data/notes.csv')

dsph_mw_name = table.join(dsph_mw, name_discovery, keys='key')
print(len(dsph_mw_name))

for i in range(len(dsph_mw_name)):
    d_t = {'key': str(dsph_mw_name['key'][i]) }
    
    location = {'ra':float(dsph_mw_name['ra'][i]), 'dec':float(dsph_mw_name['dec'][i])}
    d_t['location'] = location
    
    discovery = {'name':str(dsph_mw_name['name'][i]), 'discovery_year': int(dsph_mw_name['discovery_year'][i])}
    for prop in ['ref_discovery_1', 'ref_discovery_2', 'ref_discovery_3', 'other_name_1', 'other_name_2', 'other_name_3','other_name_4', 'other_name_5', 'other_name_6', 'table', 'host']:
        if prop in dsph_mw_name.dtype.names and ma.is_masked(dsph_mw_name[prop][i])==False:
            discovery[prop] = str(dsph_mw_name[prop][i])
    for prop in ['confirmed_real', 'confirmed_dwarf', 'confirmed_star_cluster']:
        if prop in dsph_mw_name.dtype.names and ma.is_masked(dsph_mw_name[prop][i])==False:
            discovery[prop] = int(dsph_mw_name[prop][i])
            
    d_t['name_discovery'] = discovery
    
    if ma.is_masked(dsph_mw_name['ref_structure'][i])==False:
        structure = {'ref': str(dsph_mw_name['ref_structure'][i]) }
        for prop in ['rhalf', 'rhalf_em', 'rhalf_ep', 'position_angle', 'position_angle_em', 'position_angle_ep', 'ellipticity', 'ellipticity_em', 'ellipticity_ep', 'ellipticity_ul',]:
            if prop in dsph_mw_name.dtype.names and ma.is_masked(dsph_mw_name[prop][i])==False:
                structure[prop] = float(dsph_mw_name[prop][i])
        d_t['structure'] = structure
    
    if ma.is_masked(dsph_mw_name['ref_distance'][i])==False:
        distance = {'ref': str(dsph_mw_name['ref_distance'][i]) }
        for prop in ['distance', 'distance_em', 'distance_ep', 'distance_modulus', 'distance_modulus_em', 'distance_modulus_ep',]:
            if prop in dsph_mw_name.dtype.names and ma.is_masked(dsph_mw_name[prop][i])==False:
                distance[prop] = float(dsph_mw_name[prop][i])
        d_t['distance'] = distance
    
    if ma.is_masked(dsph_mw_name['ref_mv'][i])==False:
        MV = {'ref':str(dsph_mw_name['ref_mv'][i]), }
        for prop_out,prop in zip( ['apparent_magnitude_v', 'apparent_magnitude_v_em', 'apparent_magnitude_v_ep'],['apparent_magnitude_V', 'apparent_magnitude_V_em', 'apparent_magnitude_V_ep']):
            if prop in dsph_mw_name.dtype.names and ma.is_masked(dsph_mw_name[prop][i])==False:
                MV[prop_out] = float(dsph_mw_name[prop][i])
        d_t['m_v'] = MV
    
    if ma.is_masked(dsph_mw_name['ref_vlos'][i])==False:
        vlos = {'ref': str(dsph_mw_name['ref_vlos'][i]) }
        for prop in ['vlos_systemic', 'vlos_systemic_em', 'vlos_systemic_ep', 'vlos_sigma', 'vlos_sigma_em', 'vlos_sigma_ep', 'vlos_sigma_ul',]:
            if prop in dsph_mw_name.dtype.names and ma.is_masked(dsph_mw_name[prop][i])==False:
                vlos[prop] = float(dsph_mw_name[prop][i])
        d_t['velocity'] = vlos
    
    if ma.is_masked(dsph_mw_name['ref_proper_motion'][i])==False:
        pm = {'ref': str(dsph_mw_name['ref_proper_motion'][i]) }
        for prop in ['pmra', 'pmra_em', 'pmra_ep', 'pmdec', 'pmdec_em', 'pmdec_ep', ]:
            if prop in dsph_mw_name.dtype.names and ma.is_masked(dsph_mw_name[prop][i])==False:
                pm[prop] = float(dsph_mw_name[prop][i])
        d_t['proper_motion'] = pm
    
    if ma.is_masked(dsph_mw_name['ref_metallicity'][i])==False:
        spec_metallicity = {'ref': str(dsph_mw_name['ref_metallicity'][i]) }
        for prop in ['metallicity', 'metallicity_em', 'metallicity_ep', 'metallicity_sigma', 'metallicity_sigma_em', 'metallicity_sigma_ep', ]:
            if prop in dsph_mw_name.dtype.names and ma.is_masked(dsph_mw_name[prop][i])==False:
                spec_metallicity['spec_'+prop] = float(dsph_mw_name[prop][i])
        d_t['spec_metallicity'] = spec_metallicity
    
    if ma.is_masked(dsph_mw_name['ref_king'][i])==False:
        king = {'ref': str(dsph_mw_name['ref_king'][i])}
        for prop in ['rcore', 'rcore_em', 'rcore_ep', 'rking', 'rking_em', 'rking_ep', 'king_concentration', 'king_concentration_em', 'king_concentration_ep', ]:
            if prop in dsph_mw_name.dtype.names and ma.is_masked(dsph_mw_name[prop][i])==False:
                king[prop] = float(dsph_mw_name[prop][i])
        d_t['structure_king'] = king
    
    if ma.is_masked(dsph_mw_name['ref_sersic'][i])==False:
        sersic = {'ref': str(dsph_mw_name['ref_sersic'][i])}
        for prop in ['rad_sersic', 'rad_sersic_em', 'rad_sersic_ep', 'n_sersic', 'n_sersic_em', 'n_sersic_ep',  ]:
            if prop in dsph_mw_name.dtype.names and ma.is_masked(dsph_mw_name[prop][i])==False:
                sersic[prop] = float(dsph_mw_name[prop][i])
        d_t['structure_sersic'] = sersic
    
    if ma.is_masked(dsph_mw_name['ref_age'][i])==False:
        age = {'ref': str(dsph_mw_name['ref_age'][i])}
        for prop in ['age', 'age_em', 'age_ep',   ]:
            if prop in dsph_mw_name.dtype.names and ma.is_masked(dsph_mw_name[prop][i])==False:
                age[prop] = float(dsph_mw_name[prop][i])
        d_t['age'] = age
    
    if ma.is_masked(dsph_mw_name['ref_metallicity_photometric'][i])==False:
        metallicity_photometric = {'ref': str(dsph_mw_name['ref_metallicity_photometric'][i])}
        for prop in ['metallicity_photometric', 'metallicity_photometric_em', 'metallicity_photometric_ep',   ]:
            if prop in dsph_mw_name.dtype.names and ma.is_masked(dsph_mw_name[prop][i])==False:
                metallicity_photometric[prop] = float(dsph_mw_name[prop][i])
        d_t['metallicity_photometric'] = metallicity_photometric
    
    if ma.is_masked(dsph_mw_name['ref_flux_HI'][i])==False:
        flux_HI = {'ref': str(dsph_mw_name['ref_flux_HI'][i])}
        for prop in ['flux_HI', 'flux_HI_em', 'flux_HI_ep', 'flux_HI_ul'  ]:
            if prop in dsph_mw_name.dtype.names and ma.is_masked(dsph_mw_name[prop][i])==False:
                flux_HI[prop] = float(dsph_mw_name[prop][i]) 
        d_t['flux_HI'] = flux_HI
    notes_temp = notes[notes['key']==dsph_mw_name['key'][i]]
    if len(notes_temp)>0:
        notes_list = []
        for kk in notes_temp['notes']:
#             print(kk)
            notes_list.append(str(kk))
        d_t['notes'] = notes_list
#     print(d_t)
    yaml_name = 'data_input/'+ dsph_mw_name['key'][i]+'.yaml'
    with open(yaml_name, 'w') as yaml_file:
        yaml.dump(d_t, yaml_file, default_flow_style=False)
