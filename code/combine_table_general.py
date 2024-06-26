## this takes yaml files as input and creates tables and saves them as csv (or other file types) 

import numpy as np
import matplotlib.pyplot as plt

import astropy.table as table

from astropy import units as u
import astropy.coordinates as coord

from collections import Counter

import numpy.ma as ma
import yaml

import os

path = "data_input/"
dir_list = os.listdir(path)
dir_list = [i for i in dir_list if i!='readme.md']
## this is to get a list of all "key"s which should correspond to file names (+.yaml)
table_list = []
for i in range(len(dir_list)):
    with open(path+ dir_list[i], 'r') as stream:
        try:
            stream_yaml = yaml.load(stream, Loader=yaml.Loader)
            if 'table' in stream_yaml.keys():
                table_list.append(stream_yaml['table'])
            else:
                print(stream_yaml['key'], "missing table")
        except yaml.YAMLError as exc:
            print(exc)
print("all tables")
print(Counter(table_list))

## we don't use every possible column for the "standard" table
col_name_dwarf = ['name','host','confirmed_real', 'confirmed_dwarf',  'rhalf', 'rhalf_em', 'rhalf_ep', 'position_angle', 'position_angle_em', 'position_angle_ep', 'ellipticity', 'ellipticity_em', 'ellipticity_ep', 'ellipticity_ul', 'ref_structure', 'distance_modulus', 'distance_modulus_em', 'distance_modulus_ep', 'ref_distance', 'apparent_magnitude_v', 'apparent_magnitude_v_em', 'apparent_magnitude_v_ep', 'ref_m_v', 'vlos_systemic', 'vlos_systemic_em', 'vlos_systemic_ep', 'vlos_sigma', 'vlos_sigma_em', 'vlos_sigma_ep', 'vlos_sigma_ul', 'ref_vlos', 'pmra', 'pmra_em', 'pmra_ep', 'pmdec', 'pmdec_em', 'pmdec_ep', 'ref_proper_motion', 'metallicity_spectroscopic', 'metallicity_spectroscopic_em', 'metallicity_spectroscopic_ep', 'metallicity_spectroscopic_sigma', 'metallicity_spectroscopic_sigma_em', 'metallicity_spectroscopic_sigma_ep', 'metallicity_spectroscopic_sigma_ul', 'ref_metallicity_spectroscopic', 'rcore', 'rcore_em', 'rcore_ep', 'rking', 'rking_em', 'rking_ep', 'ref_structure_king', 'rad_sersic', 'rad_sersic_em', 'rad_sersic_ep', 'n_sersic', 'n_sersic_em', 'n_sersic_ep', 'ref_structure_sersic', 'age', 'age_em', 'age_ep', 'ref_age', 'metallicity_photometric', 'metallicity_photometric_em', 'metallicity_photometric_ep', 'ref_metallicity_photometric', 'flux_HI', 'flux_HI_em', 'flux_HI_ep', 'flux_HI_ul', 'ref_flux_HI']
col_type_dwarf = [ 'U100','U100','i1', 'i1','f8', 'f8', 'f8','f8', 'f8', 'f8','f8', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8' ,'f8', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'f8','U100']
print(len(col_name_dwarf), len(col_type_dwarf))

col_name_gc = ['name','host','confirmed_real', 'confirmed_star_cluster',  'rhalf', 'rhalf_em', 'rhalf_ep', 'position_angle', 'position_angle_em', 'position_angle_ep', 'ellipticity', 'ellipticity_em', 'ellipticity_ep', 'ellipticity_ul', 'ref_structure', 'distance_modulus', 'distance_modulus_em', 'distance_modulus_ep', 'ref_distance', 'apparent_magnitude_v', 'apparent_magnitude_v_em', 'apparent_magnitude_v_ep', 'ref_m_v', 'vlos_systemic', 'vlos_systemic_em', 'vlos_systemic_ep', 'vlos_sigma', 'vlos_sigma_em', 'vlos_sigma_ep', 'vlos_sigma_ul', 'ref_vlos', 'pmra', 'pmra_em', 'pmra_ep', 'pmdec', 'pmdec_em', 'pmdec_ep', 'ref_proper_motion', 'metallicity_spectroscopic', 'metallicity_spectroscopic_em', 'metallicity_spectroscopic_ep', 'metallicity_spectroscopic_sigma', 'metallicity_spectroscopic_sigma_em', 'metallicity_spectroscopic_sigma_ep', 'metallicity_spectroscopic_sigma_ul', 'ref_metallicity_spectroscopic', 'rcore', 'rcore_em', 'rcore_ep', 'rking', 'rking_em', 'rking_ep', 'ref_structure_king', 'rad_sersic', 'rad_sersic_em', 'rad_sersic_ep', 'n_sersic', 'n_sersic_em', 'n_sersic_ep', 'ref_structure_sersic', 'age', 'age_em', 'age_ep', 'ref_age', 'metallicity_photometric', 'metallicity_photometric_em', 'metallicity_photometric_ep', 'ref_metallicity_photometric']
col_type_gc = [ 'U100','U100','i1', 'i1','f8', 'f8', 'f8','f8', 'f8', 'f8','f8', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8' ,'f8', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'U100', 'f8', 'f8', 'f8', 'U100',]
print(len(col_name_gc), len(col_type_gc))

## create gc tables
comb_gc_ufsc = table.Table(np.zeros((Counter(table_list)['gc_ufsc'], 3)),names=('key', 'ra', 'dec' ), dtype=('U100','f8', 'f8', ))
comb_gc_harris = table.Table(np.zeros((Counter(table_list)['gc_harris'], 3)),names=('key', 'ra', 'dec' ), dtype=('U100','f8', 'f8', ))
comb_gc_disk = table.Table(np.zeros((Counter(table_list)['gc_disk'], 3)),names=('key', 'ra', 'dec' ), dtype=('U100','f8', 'f8', ))
comb_gc_dwarf = table.Table(np.zeros((Counter(table_list)['gc_dwarf_hosted'], 3)),names=('key', 'ra', 'dec' ), dtype=('U100','f8', 'f8', ))

## add all the columns.  The columns are masked for missing entries 
for i,j in zip(col_name_gc, col_type_gc):
## the only difference between GC and dwarf is whether HI gas columns are included.  Have empty columns is fine.
# for i,j in zip(col_name_dwarf, col_type_dwarf):
    comb_gc_ufsc[i] = np.ma.masked_all(len(comb_gc_ufsc), dtype=j)
    comb_gc_harris[i] = np.ma.masked_all(len(comb_gc_harris), dtype=j)
    comb_gc_disk[i] = np.ma.masked_all(len(comb_gc_disk), dtype=j)
    comb_gc_dwarf[i] = np.ma.masked_all(len(comb_gc_dwarf), dtype=j)

## create dwarf tables
comb_dwarf_mw = table.Table(np.zeros((Counter(table_list)['dwarf_mw'], 3)),names=('key', 'ra', 'dec' ), dtype=('U100','f8', 'f8', ))
comb_dwarf_m31 = table.Table(np.zeros((Counter(table_list)['dwarf_m31'], 3)),names=('key', 'ra', 'dec' ), dtype=('U100','f8', 'f8', ))
comb_dwarf_lf = table.Table(np.zeros((Counter(table_list)['dwarf_local_field'], 3)),names=('key', 'ra', 'dec' ), dtype=('U100','f8', 'f8', ))
comb_dwarf_lf_distant = table.Table(np.zeros((Counter(table_list)['dwarf_local_field_distant'], 3)),names=('key', 'ra', 'dec' ), dtype=('U100','f8', 'f8', ))

for i,j in zip(col_name_dwarf, col_type_dwarf):
    comb_dwarf_mw[i] = np.ma.masked_all(len(comb_dwarf_mw), dtype=j)
    comb_dwarf_m31[i] = np.ma.masked_all(len(comb_dwarf_m31), dtype=j)
    comb_dwarf_lf[i] = np.ma.masked_all(len(comb_dwarf_lf), dtype=j)
    comb_dwarf_lf_distant[i] = np.ma.masked_all(len(comb_dwarf_lf_distant), dtype=j)

    # comb_gc_ufsc[i] = np.ma.masked_all(len(comb_gc_ufsc), dtype=j)

## distance modulus
def dist_mod(mu, mu_em=0, mu_ep=0):
    def dm(x):
        return pow(10., x/5.+1.)/1000.
    return dm(mu), dm(mu)-dm(mu-mu_em), dm(mu+mu_ep)-dm(mu)

def add_to_table(yaml_input, table_output, place, ):
    ## this function add input to tables
    missing_key=[]
    table_output['key'][place] = yaml_input['key']
    table_output['ra'][place] = yaml_input['location']['ra']
    table_output['dec'][place] = yaml_input['location']['dec']
    for y in ['structure', 'distance','m_v', 'velocity', 'proper_motion', 'metallicity_spectroscopic', 'structure_king', 'structure_sersic', 'age', 'metallicity_photometric', 'flux_HI', 'name_discovery']:
        if y in yaml_input.keys():
            x = list(yaml_input[y].keys())
            ## special case of distance fixed to the value in another file (ie host)
            if y == 'distance' and 'distance_fixed_host' in x and stream_yaml['name_discovery']['host'] + '.yaml' in dir_list:
                try:
                # if True:
                    # print("host missing for", stream_yaml['key'], path + stream_yaml['name_discovery']['host'] + '.yaml')
                    with open(path + stream_yaml['name_discovery']['host'] + '.yaml', 'r') as yaml_to_load:
                        stream_yaml_host = yaml.load(yaml_to_load, Loader=yaml.Loader)
                        # print("stream_yaml_host", stream_yaml_host)
                        x = list(stream_yaml_host['distance'].keys())
                        # print(x)
                        for list_key in range(len(x)):
                            name = x[list_key]
                            # print(name)
                            if name not in table_output.dtype.names:
                                missing_key.append(name)
                                continue
                            else:
                                table_output[name][place] = stream_yaml_host['distance'][x[list_key]]
                except:
                    print("host missing for", stream_yaml['key'])
                        # print(stream_yaml['name_discovery']['host'] + '.yaml' in dir_list)
                # continue
            ## dealing with units 
            elif y in ['structure', 'structure_king', 'structure_sersic'] and 'spatial_units' in x:
                unit_conversion = 1.
                if yaml_input[y]['spatial_units'].replace(' ', '') =='arcsec':
                    unit_conversion = 1./60. ## convert from arcsec to arcmin
                    # print(yaml_input['key'], unit_conversion)
                for list_key in range(len(x)):
                    name = x[list_key]
                    
                    if name not in table_output.dtype.names:
                        missing_key.append(name)
                        continue
                    else:
                        if name in ['rcore', 'rcore_em', 'rcore_ep', 'rking', 'rking_em', 'rking_ep',  'rad_sersic', 'rad_sersic_em', 'rad_sersic_ep','rhalf', 'rhalf_em', 'rhalf_ep']:
                            table_output[name][place] = yaml_input[y][x[list_key]] * unit_conversion
                            # print(name, yaml_input['key'])
                        else:
                            table_output[name][place] = yaml_input[y][x[list_key]]
            ## everything else
            else:
                x = list(yaml_input[y].keys())
                for list_key in range(len(x)):
                    name = x[list_key]
                    if name not in table_output.dtype.names:
                        missing_key.append(name)
                        continue
                    else:
                        table_output[name][place] = yaml_input[y][x[list_key]]
    return missing_key

def value_add(input_table, table_type='dwarf', **kwargs):
    spatial_units = kwargs.get("spatial_units", "arcmin")
    spatial_units_conversion = 60.
    if spatial_units == 'arcsec':
        spatial_units_conversion = 3600.
    ## value added columns
    input_table['M_V'] = input_table['apparent_magnitude_v']-input_table['distance_modulus']
    input_table['M_V_em'] = input_table['apparent_magnitude_v_em']
    input_table['M_V_ep'] = input_table['apparent_magnitude_v_ep']

    def lum(m_x, m_x_sun=4.83):
        return pow(10., -.4*(m_x - m_x_sun) )
    input_table['mass_stellar'] = np.log10(lum(input_table['M_V']) * 2.)

    d, dem, dep = dist_mod( input_table['distance_modulus'], input_table['distance_modulus_em'], input_table['distance_modulus_ep'])
    input_table['distance'] = d
    input_table['distance_em'] = dem
    input_table['distance_ep'] = dep

    input_table['rhalf_physical'] = input_table['distance']*1000.*input_table['rhalf']/spatial_units_conversion/180.*np.pi
    input_table['rhalf_sph_physical'] = np.ma.masked_all(len(input_table), dtype=float)
    for i in range(len(input_table)):
        if ma.is_masked(input_table['ellipticity'][i])==False:
            input_table['rhalf_sph_physical'][i] = input_table['rhalf_physical'][i]*np.sqrt(1.-input_table['ellipticity'][i])
        else:
            input_table['rhalf_sph_physical'][i] = input_table['rhalf_physical'][i]

    input_table['surface_brightness_rhalf'] = input_table['M_V'] + 19.78 + input_table['distance_modulus'] +  2.5 * np.log10(np.degrees(np.arctan(input_table['rhalf_sph_physical']/1000./input_table['distance']))**2)

    if table_type=='dwarf':
        input_table['mass_HI'] = np.log10(235600 * input_table['flux_HI']*(input_table['distance']/1000.)**2 )
        input_table['mass_HI_ul'] = np.log10(235600 * input_table['flux_HI_ul']*(input_table['distance']/1000.)**2 )
    
    input_table['metallicity'] = np.ma.masked_all(len(input_table), dtype=float)
    input_table['metallicity_em'] = np.ma.masked_all(len(input_table), dtype=float)
    input_table['metallicity_ep'] = np.ma.masked_all(len(input_table), dtype=float)
    input_table['metallicity_type'] = np.ma.masked_all(len(input_table), dtype='U100')
    for i in range(len(input_table)):
        if ma.is_masked(input_table['metallicity_spectroscopic'][i])==False:
            input_table['metallicity'][i] = input_table['metallicity_spectroscopic'][i]
            input_table['metallicity_em'][i] = input_table['metallicity_spectroscopic_em'][i]
            input_table['metallicity_ep'][i] = input_table['metallicity_spectroscopic_ep'][i]
            input_table['metallicity_type'][i] = 'spectroscopic'
        elif  ma.is_masked(input_table['metallicity_photometric'][i])==False:
            input_table['metallicity'][i] = input_table['metallicity_photometric'][i]
            input_table['metallicity_em'][i] = input_table['metallicity_photometric_em'][i]
            input_table['metallicity_ep'][i] = input_table['metallicity_photometric_ep'][i]
            input_table['metallicity_type'][i] = 'photometric'
    
    c_table_input = coord.SkyCoord(ra=input_table['ra']*u.deg, dec=input_table['dec']*u.deg,  frame='icrs',)
    input_table['ll'] = c_table_input.galactic.l.value
    input_table['bb'] = c_table_input.galactic.b.value
    
    return input_table


missing_key = []
missing_table = []
missing_table_key = []
place_dwarf_mw = 0
place_dwarf_m31 = 0
place_dwarf_local_field = 0
place_dwarf_local_field_distant = 0
place_gc_harris = 0
place_gc_disk = 0
place_gc_ufsc =0 
place_gc_dwarf =0 
example_keys= ['discovery_year', 'other_name', 'ref_discovery', 'type', 'spatial_units', 'central_surface_brightness', 'central_surface_brightness_em', 'central_surface_brightness_ep', 'false_positive', 'metallicity_photometric_sigma', 'mean_ebv', 'king_concentration', 'king_concentration_em', 'king_concentration_ep', 'abbreviation', 'vlos_sigma_central', 'vlos_sigma_central_em', 'vlos_sigma_central_ep', 'confirmed_star_cluster', 'vlos_systemic_HI', 'vlos_systemic_HI_em', 'vlos_systemic_HI_ep', 'sigma_HI', 'sigma_HI_em', 'sigma_HI_ep', 'vrot_HI', 'vrot_HI_em', 'vrot_HI_ep', 'ref_HI_kinematics',  'metallicity_photometric_sigma_em', 'metallicity_photometric_sigma_ep', 'apparent_magnitude_v_ul', 'age_ll']

## this add each galaxy/star cluster to the tables. 
for i in range(len(dir_list)):
    with open(path+ dir_list[i], 'r') as stream:
        stream_yaml = yaml.load(stream, Loader=yaml.Loader)
        if stream_yaml['table'] == 'dwarf_mw':
            miss = add_to_table(stream_yaml, comb_dwarf_mw, place_dwarf_mw)
            place_dwarf_mw+=1
        elif stream_yaml['table'] == 'dwarf_m31':
            miss = add_to_table(stream_yaml, comb_dwarf_m31, place_dwarf_m31)
            place_dwarf_m31+=1
        elif stream_yaml['table'] == 'dwarf_local_field':
            miss = add_to_table(stream_yaml, comb_dwarf_lf, place_dwarf_local_field)
            place_dwarf_local_field+=1 
        elif stream_yaml['table'] == 'dwarf_local_field_distant':
            miss = add_to_table(stream_yaml, comb_dwarf_lf_distant, place_dwarf_local_field_distant)
            place_dwarf_local_field_distant+=1 
        elif stream_yaml['table'] == 'gc_harris':
            miss = add_to_table(stream_yaml, comb_gc_harris, place_gc_harris)
            place_gc_harris+=1
        elif stream_yaml['table'] == 'gc_disk':
            miss = add_to_table(stream_yaml, comb_gc_disk, place_gc_disk)
            place_gc_disk+=1 
        elif stream_yaml['table'] == 'gc_ufsc':
            miss = add_to_table(stream_yaml, comb_gc_ufsc, place_gc_ufsc)
            place_gc_ufsc+=1 
        elif stream_yaml['table'] == 'gc_dwarf_hosted':
            miss = add_to_table(stream_yaml, comb_gc_dwarf, place_gc_dwarf)
            place_gc_dwarf+=1 
        else:
            missing_table.append(stream_yaml['table'])
            missing_table_key.append(stream_yaml['key'])
        if len(miss)>0:
            for iter in range(len(miss)):
                    if miss[iter] not in example_keys:
                        missing_key.append(miss[iter])
                        # print('extra key',stream_yaml['key'], miss)

# print("missing yaml entry", Counter(np.concatenate(missing_key).flat))
# print(missing_key)
print(Counter(missing_key).keys())
print()
missing_table = np.array(missing_table)
missing_table_key = np.array(missing_table_key)
print("missing table", Counter(missing_table))
for missing in Counter(missing_table).keys():
    temp = missing_table_key[missing_table==missing]
    print("objects missing (key)", missing, Counter(temp))

## save output
comb_gc_ufsc = value_add(comb_gc_ufsc, table_type='gc')
# comb_gc_ufsc = value_add(comb_gc_ufsc, table_type='dwarf')

comb_gc_harris = value_add(comb_gc_harris, table_type='gc')
comb_gc_disk = value_add(comb_gc_disk, table_type='gc')
comb_gc_dwarf = value_add(comb_gc_dwarf, table_type='gc')

comb_dwarf_mw = value_add(comb_dwarf_mw, table_type='dwarf')
comb_dwarf_m31 = value_add(comb_dwarf_m31, table_type='dwarf')
comb_dwarf_lf = value_add(comb_dwarf_lf, table_type='dwarf')
comb_dwarf_lf_distant = value_add(comb_dwarf_lf_distant, table_type='dwarf') # , spatial_units='arcsec'

comb_dwarf_mw.sort('key')
comb_dwarf_m31.sort('key')
comb_dwarf_lf.sort('key')
comb_dwarf_lf_distant.sort('key')
comb_dwarf_mw.write('data/dwarf_mw.csv', format='csv', overwrite=True)
comb_dwarf_m31.write('data/dwarf_m31.csv', format='csv',overwrite=True)
comb_dwarf_lf.write('data/dwarf_local_field.csv', format='csv',overwrite=True)
comb_dwarf_lf_distant.write('data/dwarf_local_field_distant.csv', format='csv',overwrite=True)

comb_gc_disk.sort('key')
comb_gc_harris.sort('key')
comb_gc_ufsc.sort('key')
comb_gc_dwarf.sort('key')
comb_gc_disk.write('data/gc_disk.csv', format='csv',overwrite=True)
comb_gc_harris.write('data/gc_harris.csv', format='csv',overwrite=True)
comb_gc_ufsc.write('data/gc_ufsc.csv', format='csv',overwrite=True)
comb_gc_dwarf.write('data/gc_dwarf_hosted.csv', format='csv',overwrite=True)

# comb_gc_ufsc.write('data/gc_ufsc.dat', format='ascii', overwrite=True)
# comb_gc_ufsc.write('data/gc_ufsc.mrt', format='ascii.mrt', overwrite=True)

comb_dwarf_mw.write('data/dwarf_mw.fits', format='fits', overwrite=True)
comb_dwarf_m31.write('data/dwarf_m31.fits', format='fits', overwrite=True)
comb_dwarf_lf.write('data/dwarf_local_field.fits', format='fits', overwrite=True)

comb_gc_disk.write('data/gc_disk.fits', format='fits', overwrite=True)
comb_gc_harris.write('data/gc_harris.fits', format='fits', overwrite=True)
comb_gc_ufsc.write('data/gc_ufsc.fits', format='fits', overwrite=True)
comb_gc_dwarf.write('data/gc_dwarf_hosted.fits', format='fits', overwrite=True)

comb_dwarf = table.vstack([comb_dwarf_mw, comb_dwarf_m31, comb_dwarf_lf])
comb_dwarf.write('data/dwarf_all.csv', format='csv', overwrite=True)
comb_dwarf.write('data/dwarf_all.fits', format='fits', overwrite=True)
