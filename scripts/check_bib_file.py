import astropy.table as table
import bibtexparser
from collections import Counter

dsph_mw = table.Table.read('data/dwarf_mw.csv')
dsph_m31 = table.Table.read('data/dwarf_m31.csv')
dsph_lf = table.Table.read('data/dwarf_local_field.csv')
ufsc = table.Table.read('data/gc_ufsc.csv')
gc_disk = table.Table.read('data/gc_disk.csv')
gc_harris = table.Table.read('data/gc_harris.csv')

x = "data/lvdb.bib"
bib_database = bibtexparser.parse_file(x)

ads_bibcode = []
for k in bib_database.entries:
    ads_bibcode.append(k.key)


missing = []
for test_table in [dsph_mw, dsph_m31, dsph_lf, ufsc, gc_disk, gc_harris]:
    for kk in ['ref_structure',  'ref_metallicity_photometric', 'ref_age', 'ref_structure_sersic', 'ref_structure_king', 'ref_metallicity_spectroscopic', 'ref_proper_motion', 'ref_vlos',  'ref_m_v', 'ref_distance']:
        for individual in range(len(test_table[kk])):
            if test_table[kk][individual] !='--':
                if test_table[kk][individual] not in ads_bibcode and test_table[kk][individual] not in [ 'internal', 'LVGDB', 'EDD']:
                    print( test_table['key'][individual], kk, test_table[kk][individual])
                    missing.append(test_table[kk][individual])

missing_gas = []
for test_table in [dsph_mw, dsph_m31, dsph_lf, ]:
    for kk in ['ref_flux_HI', ]:
        for individual in range(len(test_table[kk])):
            if test_table[kk][individual] !='--':
                if test_table[kk][individual] not in ads_bibcode and test_table[kk][individual] not in [ 'internal', 'LVGDB', 'EDD']:
                    print( test_table['key'][individual], kk, test_table[kk][individual])
                    missing_gas.append(test_table[kk][individual])

print(Counter(missing))
print(Counter(missing_gas))