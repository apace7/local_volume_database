import numpy as np
import astropy.table as table
from astropy.io import ascii
import os



lvdb_path = os.environ.get('LVDBDIR')
print("lvdb_path", lvdb_path)

comb_all = table.Table(ascii.read(lvdb_path + 'data/comb_all.csv'))

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

# import bibtexparser
# x = lvdb_path + "table/lvdb.bib"
# bib_database = bibtexparser.parse_file(x)


unit_test()
print("unit tests completed")