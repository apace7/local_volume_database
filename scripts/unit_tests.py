import numpy as np
import astropy.table as table
from astropy.io import ascii
import os

import bibtexparser
import local_volume_database as lvdb
from collections import Counter

lvdb_path = os.environ.get('LVDBDIR')
print("lvdb_path", lvdb_path)

comb_all = table.Table(ascii.read(lvdb_path + 'data/comb_all.csv'))

lvdb.add_column(comb_all,'structure_exponential','ref_structure_exponential', col_type='U50')
lvdb.add_column(comb_all,'structure_plummer','ref_structure_plummer', col_type='U50')
lvdb.add_column(comb_all,'structure_eff','ref_structure_eff', col_type='U50')
lvdb.add_column(comb_all,'star_formation_history','ref_star_formation_history', col_type='U50')

def unit_test():
    print("unit tests started")
    ra_dec_to_large = comb_all[np.logical_or.reduce((comb_all['ra']<0, comb_all['ra']>360., comb_all['dec']<-90, comb_all['dec']>90.))]
    if len(ra_dec_to_large)==0:
        print("ra dec checked")
    else:
        print("ra_dec_to_large", len(ra_dec_to_large))
        for i in range(len(ra_dec_to_large)):
            print(i, ra_dec_to_large['key'][i], ra_dec_to_large['ra'][i], ra_dec_to_large['dec'][i])
        print()

    print('checking references')
    ref_keys = ['ref_structure', 'ref_age', 'ref_distance', 'ref_m_v', 'ref_vlos', 'ref_proper_motion', 'ref_metallicity_spectroscopic', 'ref_structure_king', 'ref_structure_sersic', 'ref_metallicity_isochrone', 'ref_flux_HI', 'ref_metallicity_photometric','ref_structure_exponential','ref_structure_plummer','ref_structure_eff','ref_star_formation_history']

    all_reference = np.array([])
    for i in ref_keys:
        all_reference = np.ma.concatenate([comb_all[i], all_reference ])

    keep = np.zeros(len(all_reference), dtype=bool)
    for i in range(len(all_reference)):
        if np.ma.is_masked(all_reference[i])==False:
            keep[i] = True
        else:
            keep[i] = False
    all_reference_clean = all_reference[keep]
    print("number of references",len(all_reference_clean), len(all_reference))

    library = bibtexparser.parse_file(lvdb_path+'table/lvdb.bib')

    bib_keys = []
    for entry in library.entries:
        bib_keys.append(entry.key)
    print("bib_keys length", len(bib_keys))
    bad_reference = []
    for i in Counter(all_reference_clean).keys():
        if i not in bib_keys:
            bad_reference.append(i)
            # print(i)
    
    for bad in bad_reference:
        print(bad)
        for ref in ref_keys:
            test_array = comb_all[comb_all[ref]==bad]
            if len(test_array)>0:
                print(test_array['key'])
    print('references checked')




unit_test()
print("unit tests completed")