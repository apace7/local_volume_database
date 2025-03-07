import yaml
import astropy.coordinates as coord
from astropy import units as u
import astropy.table as table
import os.path
import numpy.ma as ma
import numpy as np
import corner

import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages

from collections import Counter
from adjustText import adjust_text

import warnings; warnings.filterwarnings('ignore')

import local_volume_database as lvdb
## load environment variable
lvdb_path = os.environ.get('LVDBDIR')
print("lvdb_path", lvdb_path)

plt.style.use(lvdb_path + 'code/std.mplstyle')
import matplotlib as mp
mp.rcParams['text.usetex'] = True

coord.galactocentric_frame_defaults.set('v4.0')
gc_frame = coord.Galactocentric()

from galpy.potential import KeplerPotential, MWPotential2014, vesc

dwarf_all = table.Table.read(lvdb_path + 'data/dwarf_all.csv')
# dsph_m31 = table.Table.read('../data/dwarf_m31.csv')
# dsph_lf = table.Table.read('../data/dwarf_local_field.csv')
# dsph_lf_distant = table.Table.read('../data/dwarf_local_field_distant.csv')
gc_ambiguous = table.Table.read(lvdb_path + 'data/gc_ambiguous.csv')
gc_mw_new = table.Table.read(lvdb_path + 'data/gc_mw_new.csv')
gc_harris = table.Table.read(lvdb_path + 'data/gc_harris.csv')
gc_dwarf = table.Table.read(lvdb_path + 'data/gc_dwarf_hosted.csv')


color_dsph_mw = 'tab:blue'
color_dsph_m31 = 'tab:orange'
color_dsph_lf = 'tab:green'
color_dsph_lf_distant = 'tab:red'

color_gc_disk = '#9EDAE5' #'tab:brown'  #'tab:purple'
color_gc_harris = 'tab:brown'  #'#9EDAE5' #'tab:cyan'
color_gc_ufcss = 'tab:olive'

color_gc_dwarf = 'tab:pink'
color_gc_lmc_smc = 'darkgreen'

label_dsph_mw = r'${\rm Dwarf~MW}$'
label_dsph_m31 = r'${\rm Dwarf~M31}$'
label_dsph_lf = r'${\rm Dwarf~LF}$'
label_dsph_lf_distant = r'${\rm Dwarf~LV}$'
label_gc_ufcss = r'${\rm Ambiguous/HFCSS}$'
label_gc_harris = r'${\rm GC~Harris}$'
label_gc_disk =r'${\rm GC~New~Disk/Bulge/Halo}$'
label_gc_dwarf =r'${\rm GC~Dwarf}$'
label_gc_lmc_smc =r'${\rm GC~LMC/SMC}$'

dsph_mw = dwarf_all[np.logical_or(dwarf_all['host']=='mw', dwarf_all['host']=='lmc')]
dsph_mw_extra = dwarf_all[dwarf_all['distance']<600]
dsph_m31 = dwarf_all[np.logical_or(dwarf_all['host']=='m_031',dwarf_all['host']=='m_033') ]

dsph_lf = dwarf_all[dwarf_all['distance']<=3e3]
dsph_lf = dsph_lf[dsph_lf['host']!='m_031']
dsph_lf = dsph_lf[dsph_lf['host']!='m_033']
dsph_lf = dsph_lf[dsph_lf['host']!='mw']
dsph_lf = dsph_lf[dsph_lf['host']!='lmc']
dsph_lf_distant = dwarf_all[dwarf_all['distance']>3e3]
print("number of MW / M31 / LF / LV dwarf galaxies", len(dsph_mw), len(dsph_m31), len(dsph_lf), len(dsph_lf_distant))

lvdb.add_column(dsph_lf,'name_discovery','discovery_year', col_type=int)
lvdb.add_column(dsph_mw,'name_discovery','discovery_year', col_type=int)
lvdb.add_column(dsph_mw,'name_discovery','abbreviation', col_type='U100')

lvdb.add_column(dsph_m31,'name_discovery','discovery_year', col_type=int)
lvdb.add_column(gc_ambiguous,'name_discovery','discovery_year', col_type=int)

def add_columns():
    dsph_mw_extra['xx'] = np.zeros(len(dsph_mw_extra), dtype='float')
    dsph_mw_extra['yy'] = np.zeros(len(dsph_mw_extra), dtype='float')
    dsph_mw_extra['zz'] = np.zeros(len(dsph_mw_extra), dtype='float')

    dsph_mw_extra['Vx'] = np.zeros(len(dsph_mw_extra), dtype='float')
    dsph_mw_extra['Vy'] = np.zeros(len(dsph_mw_extra), dtype='float')
    dsph_mw_extra['Vz'] = np.zeros(len(dsph_mw_extra), dtype='float')

    dsph_mw_extra['Lx'] = np.zeros(len(dsph_mw_extra), dtype='float')
    dsph_mw_extra['Ly'] = np.zeros(len(dsph_mw_extra), dtype='float')
    dsph_mw_extra['Lz'] = np.zeros(len(dsph_mw_extra), dtype='float')

    dsph_mw_extra['vrad'] = np.zeros(len(dsph_mw_extra), dtype=float)
    dsph_mw_extra['vtan'] = np.zeros(len(dsph_mw_extra), dtype=float)
    dsph_mw_extra['vrad_error'] = np.zeros(len(dsph_mw_extra), dtype=float)
    dsph_mw_extra['vtan_error'] = np.zeros(len(dsph_mw_extra), dtype=float)
    dsph_mw_extra['v3d'] = np.zeros(len(dsph_mw_extra), dtype=float)
    dsph_mw_extra['v3d_error'] = np.zeros(len(dsph_mw_extra), dtype=float)

    dsph_mw_extra['vgsr'] = np.zeros(len(dsph_mw_extra), dtype=float)
    dsph_mw_extra['vgsr_em'] = np.zeros(len(dsph_mw_extra), dtype=float)
    dsph_mw_extra['vgsr_ep'] = np.zeros(len(dsph_mw_extra), dtype=float)

    dsph_mw_extra['distance_gc'] = np.zeros(len(dsph_mw_extra), dtype=float)
    dsph_mw_extra['distance_gc_em'] = np.zeros(len(dsph_mw_extra), dtype=float)
    dsph_mw_extra['distance_gc_ep'] = np.zeros(len(dsph_mw_extra), dtype=float)

    for i in range(len(dsph_mw_extra)):
        
        icrs2 = coord.SkyCoord(ra=dsph_mw_extra['ra'][i]*u.deg,
                            dec=dsph_mw_extra['dec'][i]*u.deg,
                            distance=dsph_mw_extra['distance'][i]*u.kpc,
                            pm_ra_cosdec=dsph_mw_extra['pmra'][i]*u.mas/u.yr,
                            pm_dec=dsph_mw_extra['pmdec'][i]*u.mas/u.yr,
                            radial_velocity=dsph_mw_extra['vlos_systemic'][i]*u.km/u.s)

        icrs_err2 = coord.SkyCoord(ra=0*u.deg, dec=0*u.deg, distance=(dsph_mw_extra['distance_em'][i]+dsph_mw_extra['distance_ep'][i])/2.*u.kpc,
                                pm_ra_cosdec=(dsph_mw_extra['pmra_em'][i]+dsph_mw_extra['pmra_ep'][i])/2.*u.mas/u.yr,
                                pm_dec=(dsph_mw_extra['pmdec_em'][i]+dsph_mw_extra['pmdec_ep'][i])/2.*u.mas/u.yr,
                                radial_velocity=(dsph_mw_extra['vlos_systemic_em'][i]+dsph_mw_extra['vlos_systemic_ep'][i])/2.*u.km/u.s)
        
        n_samples = 1000
        dist = np.random.normal(icrs2.distance.value, icrs_err2.distance.value,
                                n_samples) * icrs2.distance.unit

        pm_ra_cosdec = np.random.normal(icrs2.pm_ra_cosdec.value,
                                        icrs_err2.pm_ra_cosdec.value,
                                        n_samples) * icrs2.pm_ra_cosdec.unit

        pm_dec = np.random.normal(icrs2.pm_dec.value,
                                icrs_err2.pm_dec.value,
                                n_samples) * icrs2.pm_dec.unit

        rv = np.random.normal(icrs2.radial_velocity.value, icrs_err2.radial_velocity.value,
                            n_samples) * icrs2.radial_velocity.unit

        ra = np.full(n_samples, icrs2.ra.degree) * u.degree
        dec = np.full(n_samples, icrs2.dec.degree) * u.degree
        
        icrs_samples = coord.SkyCoord(ra=ra, dec=dec, distance=dist,
                                pm_ra_cosdec=pm_ra_cosdec,
                                pm_dec=pm_dec, radial_velocity=rv)
        galcen_samples = icrs_samples.transform_to(gc_frame)
        
        
        r = np.sqrt(galcen_samples.x**2 + galcen_samples.y**2 + galcen_samples.z**2)
        vrad = galcen_samples.x*galcen_samples.v_x/r + galcen_samples.y*galcen_samples.v_y/r + galcen_samples.z*galcen_samples.v_z/r
        vtan = np.sqrt((r*galcen_samples.v_z - galcen_samples.z *vrad)**2+(galcen_samples.x*galcen_samples.v_y - galcen_samples.y*galcen_samples.v_x)**2) /np.sqrt(galcen_samples.x**2 + galcen_samples.y**2)
    #     print(dsph_mw_extra['key'][i], vrad.value, vtan.value)
        pers_quant = corner.quantile(vrad.value, [.5, .1587, .8413], )
        dsph_mw_extra['vrad'][i] = pers_quant[0]
        dsph_mw_extra['vrad_error'][i] = (pers_quant[2]-pers_quant[1])/2.
        
        pers_quant = corner.quantile(vtan.value, [.5, .1587, .8413], )
        dsph_mw_extra['vtan'][i] = pers_quant[0]
        dsph_mw_extra['vtan_error'][i] = (pers_quant[2]-pers_quant[1])/2.
    #     print("rad",dsph_mw_extra['key'][i], dsph_mw_extra['vrad'][i], dsph_mw_extra['vrad_error'][i])
    #     print(dsph_mw_extra['key'][i], dsph_mw_extra['vtan'][i], dsph_mw_extra['vtan_error'][i])
        pers_quant = corner.quantile(vtan.value, [.5, .1587, .8413], )
        dsph_mw_extra['vtan'][i] = pers_quant[0]
        dsph_mw_extra['vtan_error'][i] = (pers_quant[2]-pers_quant[1])/2.
    #     print(vrad.value**2 + vtan.value**2, galcen_samples.v_x.value**2+galcen_samples.v_y.value**2+galcen_samples.v_z.value**2)
        pers_quant = corner.quantile(np.sqrt(vrad.value**2 + vtan.value**2), [.5, .1587, .8413], )
        dsph_mw_extra['v3d'][i] = pers_quant[0]
        dsph_mw_extra['v3d_error'][i] = (pers_quant[2]-pers_quant[1])/2.
    #     print(dsph_mw_extra['v3d'][i], dsph_mw_extra['v3d_error'][i])
    #     print(dsph_mw_extra['key'][i], icrs3.x.value, icrs3.y.value, icrs3.z.value)
    #     print('\t', icrs3.v_x.value, icrs3.v_y.value, icrs3.v_z.value)
        
        pers_quant = corner.quantile(galcen_samples.radial_velocity.value, [.5, .1587, .8413], )
        dsph_mw_extra['vgsr'][i] = pers_quant[0]
        dsph_mw_extra['vgsr_ep'][i] = (pers_quant[2]-pers_quant[0])
        dsph_mw_extra['vgsr_em'][i] = (pers_quant[0]-pers_quant[1])
        
        pers_quant = corner.quantile(np.sqrt(galcen_samples.x.value**2 + galcen_samples.y.value**2 + galcen_samples.z.value**2), [.5, .1587, .8413], )
        dsph_mw_extra['distance_gc'][i] = pers_quant[0]
        dsph_mw_extra['distance_gc_ep'][i] = (pers_quant[2]-pers_quant[0])
        dsph_mw_extra['distance_gc_em'][i] = (pers_quant[0]-pers_quant[1])

add_columns()


## compuates average density within the half-light radius mass and errors
## this assumes Gaussian errors and averages the errors 
def compute_density_error(rhalf, rhalf_em, rhalf_ep, ellipticity, ellipticity_em, ellipticity_ep, distance, distance_em, distance_ep, sigma, sigma_em,sigma_ep, n=10000):
    if (ma.is_masked(sigma_em)==True or ma.is_masked(sigma_ep)==True) and ma.is_masked(sigma)==False:
        rh = distance*rhalf/180./60.*1000.*np.pi
        comb_mass = 930. * rh * sigma**2
        
        if ma.is_masked(ellipticity)==False:
            rh = rh*np.sqrt(1.-ellipticity)
            return [comb_mass*np.sqrt(1.-ellipticity)/(4.*np.pi/3.*(4./3.*rh)**3),0,0]
        else:
            return[comb_mass/(4.*np.pi/3.*(4./3.*rh)**3),0,0]
    else: 
        if ma.is_masked(rhalf_em)==False and ma.is_masked(rhalf)==False:
            x = np.random.normal(rhalf, (rhalf_em+rhalf_ep)/2., n)
        elif ma.is_masked(rhalf)==False:
            x=np.empty(n)
            x.fill(rhalf)
        else:
            x = np.zeros(n)
            return [np.ma.masked,np.ma.masked,np.ma.masked]
        
        if ma.is_masked(ellipticity_em)==False and ma.is_masked(ellipticity)==False:
            y = np.random.normal(ellipticity, (ellipticity_em+ellipticity_ep)/2., n)
        elif ma.is_masked(ellipticity)==False:
            y=np.empty(len(x))
            y.fill(ellipticity)
        else:
            y = np.zeros(len(x))
        
        if ma.is_masked(distance_em)==False and ma.is_masked(distance)==False:
            z = np.random.normal(distance, (distance_em+distance_ep)/2., n)
        elif ma.is_masked(distance)==False:
            z=np.empty(len(x))
            z.fill(distance)
        else:
            z = np.full(10000, distance)
        
        if ma.is_masked(sigma_em)==False and ma.is_masked(sigma_ep)==False:
            sig = np.random.normal(sigma, (sigma_em+sigma_ep)/2., n)
        elif ma.is_masked(sigma)==False:
            sig=np.empty(len(x))
            sig.fill(sigma)
        else:
            sig = np.zeros(len(x))
        
        x2 = x[np.logical_and(y>=0, y<1)]
        y2 = y[np.logical_and(y>=0, y<1)]
        z2 = z[np.logical_and(y>=0, y<1)]
        sig2 = sig[np.logical_and(y>=0, y<1)]
        
        rh_array = x2 *np.pi/180./60.*1000.*np.sqrt(1. - y2)* z2
        comb_mass = 930. * x2 *np.pi/180./60.*1000.*np.sqrt(1. - y2)* z2 * sig2**2
        comb_mass2 = comb_mass[~np.isnan(comb_mass)]
        rh_array2 = rh_array[~np.isnan(comb_mass)]

        comb_mass2_density = comb_mass2/(4./3.*np.pi * (4./3.*rh_array2)**3)
        mean_mass = corner.quantile(comb_mass2_density, [.5, .1587, .8413, 0.0227501, 0.97725])
        return [mean_mass[0], mean_mass[0]-mean_mass[1], mean_mass[2]-mean_mass[0]]
    
def add_density(input_array):
    input_array['density_dyn_mcmc'] = np.zeros(len(input_array), dtype=float)
    input_array['density_dyn_mcmc_em'] = np.zeros(len(input_array), dtype=float)
    input_array['density_dyn_mcmc_ep'] = np.zeros(len(input_array), dtype=float)
    input_array['density_dyn_mcmc_ul'] = np.zeros(len(input_array), dtype=float)

    for i in range(len(input_array)):
        y= compute_density_error(input_array['rhalf'][i], input_array['rhalf_em'][i], input_array['rhalf_ep'][i], input_array['ellipticity'][i], input_array['ellipticity_em'], input_array['ellipticity_ep'][i], input_array['distance'][i], input_array['distance_em'][i], input_array['distance_ep'][i],input_array['vlos_sigma'][i], input_array['vlos_sigma_em'][i], input_array['vlos_sigma_ep'][i])

        input_array['density_dyn_mcmc'][i] = y[0]
        input_array['density_dyn_mcmc_em'][i] = y[1]
        input_array['density_dyn_mcmc_ep'][i] = y[2]

        z= compute_density_error(input_array['rhalf'][i], input_array['rhalf_em'][i], input_array['rhalf_ep'][i], input_array['ellipticity'][i], input_array['ellipticity_em'], input_array['ellipticity_ep'][i], input_array['distance'][i], input_array['distance_em'][i], input_array['distance_ep'][i],input_array['vlos_sigma_ul'][i], np.ma.masked, np.ma.masked)
        input_array['density_dyn_mcmc_ul'][i] = z[0]
    # return input_array

add_density(dsph_m31)
add_density(dsph_mw)
add_density(dsph_lf)
add_density(gc_ambiguous)


## there are generally two plots created, one plot is added to the overview_plots.pdf and the second is saved in the same location as the paper_example/.ipynb notebooks
## axis limits are generally not included in overview_plots.pdf so that outliers can be identified.

with PdfPages(lvdb_path + 'paper_examples/overview_plots.pdf') as pdf:
    fig = plt.figure()

    ### Figure 1
    sorted_data = np.sort(dsph_mw['discovery_year']) 
    plt.step(sorted_data, np.arange(sorted_data.size), where='pre' , c=color_dsph_mw, label=label_dsph_mw)
    plt.plot( [np.max(sorted_data), np.max(sorted_data)+100], [len(sorted_data)-1, len(sorted_data)-1], c=color_dsph_mw)

    sorted_data = np.sort(dsph_m31['discovery_year']) 
    plt.step(sorted_data, np.arange(sorted_data.size), where='pre', c=color_dsph_m31, label=label_dsph_m31)
    plt.plot( [np.max(sorted_data), np.max(sorted_data)+100], [len(sorted_data)-1, len(sorted_data)-1], c=color_dsph_m31)

    sorted_data = np.sort(dsph_lf['discovery_year']) 
    plt.step(sorted_data, np.arange(sorted_data.size), where='pre' , c=color_dsph_lf, label=label_dsph_lf)
    plt.plot( [np.max(sorted_data), np.max(sorted_data)+100], [len(sorted_data)-1, len(sorted_data)-1], c=color_dsph_lf)

    sorted_data = np.sort(gc_ambiguous['discovery_year']) 
    plt.step(sorted_data, np.arange(sorted_data.size), where='pre' , c=color_gc_ufcss, label=label_gc_ufcss)
    plt.plot( [np.max(sorted_data), np.max(sorted_data)+100], [len(sorted_data)-1, len(sorted_data)-1], c=color_gc_ufcss)

    plt.xlim(1920, 2030)
    plt.ylim(0, 75)
    plt.xlabel(r'${\rm Discovery~Year}$')
    plt.ylabel(r'${\rm Cumulative~Number}$')
    ## add lines denoting eras (SDSS, current, LSST)
    plt.axvline(2005, c='k', ls=":", )
    plt.axvline(2015, c='k', ls=":", )
    plt.axvline(2025, c='k', ls=":", )

    plt.legend(loc=2)
    plt.tight_layout()
    plt.savefig(lvdb_path + 'paper_examples/milky_way_m31_lf_discovery_year.pdf')
    pdf.savefig()
    plt.close()
    print("figure 1 finished")

    ### Figure 2
    comb_ngc253 = dwarf_all[dwarf_all['host']=='ngc_0253']
    comb_cena = dwarf_all[dwarf_all['host']=='ngc_5128']

    for i,j,k in zip([dsph_mw, dsph_m31, dsph_lf, comb_ngc253], [color_dsph_mw, color_dsph_m31, color_dsph_lf, 'k',], [label_dsph_mw, label_dsph_m31, label_dsph_lf, r'{\rm NGC~253~Satellite}']):
        sorted_data = np.sort(i['M_V']) 
        plt.step(sorted_data, np.arange(sorted_data.size), label=k, c=j)

    plt.xlim(2.5, -20)
    plt.yscale('log')
    plt.xlabel(r'${\rm  M_V}$')
    plt.ylabel(r'${\rm N(<M_V)}$')
    plt.legend(loc=3)
    plt.tight_layout()
    plt.savefig(lvdb_path + 'paper_examples/cumulative_distribution.pdf')
    
    for i,j,k in zip([comb_cena], ['tab:red'], [r'{\rm Cen~A~Satellite}']):
        sorted_data = np.sort(i['M_V']) 
        plt.step(sorted_data, np.arange(sorted_data.size), label=k, c=j)
    plt.legend(loc=3)
    plt.tight_layout()
    pdf.savefig()
    plt.close()
    print("figure 2 finished")
    ### Figure 3 part 1
    c_dsph_mw = coord.SkyCoord(ra=dsph_mw['ra']*u.deg, dec=dsph_mw['dec']*u.deg,  frame='icrs', distance=dsph_mw['distance']*u.kpc, radial_velocity=dsph_mw['vlos_systemic']*u.km/u.s, pm_ra_cosdec=dsph_mw['pmra']*u.mas/u.yr,pm_dec=dsph_mw['pmdec']*u.mas/u.yr)

    fig = plt.figure(1,figsize=(16*1.2,6*1.2))
    ax = fig.add_subplot(111, projection='mollweide')
    cb = ax.scatter(c_dsph_mw.galactic.l.wrap_at(180*u.deg).rad, c_dsph_mw.galactic.b.rad  , c=c_dsph_mw.distance.value, label=r'${\rm Dwarf~MW}$', cmap='bwr', s=dsph_mw['mass_stellar']**2.5, zorder=10, ec='k')
    plt.colorbar(cb)
    ## colorbar is based on distance
    ## size of points is based on stellar mass (mass-to-light ratio =2)

    ## this labels the objects
    texts = [plt.text(c_dsph_mw.galactic.l.wrap_at(180*u.deg).rad[i], c_dsph_mw.galactic.b.rad[i], r'${\rm '+dsph_mw['abbreviation'][i]+'}$', ha='center', va='center') for i in range(len(c_dsph_mw))]

    adjust_text(texts, arrowprops=dict(arrowstyle='->', color='grey'))

    ## this adds ra, dec=0 grey line
    gal_long = np.arange(-180., 180., .1)
    c_gc = coord.SkyCoord(ra=gal_long*u.deg, dec=np.zeros(len(gal_long))*u.deg,  frame='icrs')
    temp = c_gc.transform_to('galactic')
    ax.plot(temp.galactic.l.wrap_at(180*u.deg).rad, temp.galactic.b.rad  , '.', c='k', ms=.1, rasterized=True,zorder=1)

    ax.grid(ls=':')

    plt.tight_layout()
    plt.savefig(lvdb_path + 'paper_examples/aitoff_MW_v3.pdf', transparent=True)
    pdf.savefig()
    plt.close()

    ### Figure 3 part 2
    c_dsph_m31 = coord.SkyCoord(ra=dsph_m31['ra']*u.deg, dec=dsph_m31['dec']*u.deg,  frame='icrs', distance=dsph_m31['distance']*u.kpc)
    c_dsph_lf = coord.SkyCoord(ra=dsph_lf['ra']*u.deg, dec=dsph_lf['dec']*u.deg,  frame='icrs', distance=dsph_lf['distance']*u.kpc)
    c_lf_distant = coord.SkyCoord(ra=dsph_lf_distant['ra']*u.deg, dec=dsph_lf_distant['dec']*u.deg,  frame='icrs', distance=dsph_lf_distant['distance']*u.kpc)
    fig = plt.figure(1,figsize=(16*.9,9*.9))
    ax = fig.add_subplot(111, projection='mollweide')

    ax.plot(c_dsph_m31.galactic.l.wrap_at(180*u.deg).rad, c_dsph_m31.galactic.b.rad  , 'o', label=r'${\rm Dwarf~M31}$', mec='k', c=color_dsph_m31)
    ax.plot(c_dsph_lf.galactic.l.wrap_at(180*u.deg).rad, c_dsph_lf.galactic.b.rad  , 'o', label=r'${\rm Dwarf~Local~Field}$', mec='k',c=color_dsph_lf)
    ax.plot(c_lf_distant.galactic.l.wrap_at(180*u.deg).rad, c_lf_distant.galactic.b.rad  , 'o', label=r'${\rm Dwarf~Local~Volume}$', mec='k',c=color_dsph_lf_distant)

    ax.grid(ls=':')
    ax.legend(loc=4)
    plt.tight_layout()
    plt.savefig(lvdb_path + 'paper_examples/aitoff_dwarf_LF.pdf',transparent=True)
    pdf.savefig()
    plt.close()

    ### Figure 3 part 2b
    fig = plt.figure(1,figsize=(16*1.2,6*1.2))
    ax = fig.add_subplot(111, projection='mollweide')

    ax.scatter(c_dsph_m31.galactic.l.wrap_at(180*u.deg).rad, c_dsph_m31.galactic.b.rad  , label=r'${\rm Dwarf~M31}$', c=dsph_m31['distance']/1000., cmap='bwr', vmin=0., vmax=11, s=dsph_m31['mass_stellar']**2)
    ax.scatter(c_dsph_lf.galactic.l.wrap_at(180*u.deg).rad, c_dsph_lf.galactic.b.rad , label=r'${\rm Dwarf~LF}$', c=dsph_lf['distance']/1000., cmap='bwr', vmin=0., vmax=11, s=dsph_lf['mass_stellar']**2)
    cb = ax.scatter(c_lf_distant.galactic.l.wrap_at(180*u.deg).rad, c_lf_distant.galactic.b.rad  ,  label=r'${\rm Dwarf~LF~Distant}$', c=dsph_lf_distant['distance']/1000., cmap='bwr', vmin=0., vmax=11, s=dsph_lf_distant['mass_stellar']**2)
    plt.colorbar(cb)

    ax.grid(ls=':')

    plt.tight_layout()
    pdf.savefig()
    plt.close()

    ### Figure 3 part 3
    gc_harris2= gc_harris[gc_harris['host']=='mw']
    c_gc_harris2 = coord.SkyCoord(ra=gc_harris2['ra']*u.deg, dec=gc_harris2['dec']*u.deg,  frame='icrs', distance=gc_harris2['distance']*u.kpc)
    
    # keep = np.zeros(len(gc_harris), dtype=bool)
    # for i in range(len(gc_harris)):
    #     if gc_harris['host'][i] in ['sagittarius_1']:
    #         keep[i]=True
    #     else:
    #         keep[i]=False
    # gc_dwarf_sgr = gc_harris[keep]
    gc_dwarf_sgr = gc_harris[gc_harris['host'] == 'sagittarius_1'] 
    c_gc_dwarf_sgr = coord.SkyCoord(ra=gc_dwarf_sgr['ra']*u.deg, dec=gc_dwarf_sgr['dec']*u.deg,  frame='icrs', distance=gc_dwarf_sgr['distance']*u.kpc)

    c_gc_ambiguous = coord.SkyCoord(ra=gc_ambiguous['ra']*u.deg, dec=gc_ambiguous['dec']*u.deg,  frame='icrs', distance=gc_ambiguous['distance']*u.kpc)

    c_gc_disk = coord.SkyCoord(ra=gc_mw_new['ra']*u.deg, dec=gc_mw_new['dec']*u.deg,  frame='icrs', distance=gc_mw_new['distance']*u.kpc)
    keep = np.zeros(len(gc_dwarf), dtype=bool)
    for i in range(len(gc_dwarf)):
        if gc_dwarf['host'][i] in ['lmc', 'smc']:
            keep[i]=True
        else:
            keep[i]=False
    lmc_smc = gc_dwarf[keep]
    c_lmc_smc = coord.SkyCoord(ra=lmc_smc['ra']*u.deg, dec=lmc_smc['dec']*u.deg,  frame='icrs', distance=lmc_smc['distance']*u.kpc)

    keep = np.zeros(len(gc_dwarf), dtype=bool)
    for i in range(len(gc_dwarf)):
        if gc_dwarf['host'][i] not in ['lmc', 'smc'] and gc_dwarf['distance'][i] <400:
            keep[i]=True
        else:
            keep[i]=False
    gc_dwarf_mw = gc_dwarf[keep]
    c_gc_dwarf_mw = coord.SkyCoord(ra=gc_dwarf_mw['ra']*u.deg, dec=gc_dwarf_mw['dec']*u.deg,  frame='icrs', distance=gc_dwarf_mw['distance']*u.kpc)
    print("figure 3, number of LMC/SMC cluster", len(gc_dwarf_mw), len(lmc_smc))

    fig = plt.figure(1,figsize=(16*.9,9*.9))
    ax = fig.add_subplot(111, projection='mollweide')

    ax.plot(c_gc_ambiguous.galactic.l.wrap_at(180*u.deg).rad, c_gc_ambiguous.galactic.b.rad  , 'o', label=r'${\rm Ambiguous/HFCSS}$', c=color_gc_ufcss, mec='k', ms=6)
    ax.plot(c_gc_disk.galactic.l.wrap_at(180*u.deg).rad, c_gc_disk.galactic.b.rad  , 's', label=r'${\rm GC~New~Bulge/Disk/Halo}$', c=color_gc_disk, mec='k', ms=6)
    # color_gc_disk
    ax.plot(c_gc_harris2.galactic.l.wrap_at(180*u.deg).rad, c_gc_harris2.galactic.b.rad  , 's', label=r'${\rm GC~Harris}$', c=color_gc_harris, mec='k', ms=6)
    #tab:cyan
    # color_gc_harris
    ax.plot(c_lmc_smc.galactic.l.wrap_at(180*u.deg).rad, c_lmc_smc.galactic.b.rad  , 's', label=r'${\rm GC~LMC/SMC}$', c=color_gc_lmc_smc, mec='k', ms=6)
    ax.plot(c_gc_dwarf_sgr.galactic.l.wrap_at(180*u.deg).rad, c_gc_dwarf_sgr.galactic.b.rad  , 's', label=r'${\rm GC~Fornax/Sgr/Dwarf~Hosted}$', c=color_gc_dwarf, mec='k', ms=6)
    ax.plot(c_gc_dwarf_mw.galactic.l.wrap_at(180*u.deg).rad, c_gc_dwarf_mw.galactic.b.rad  , 's', c=color_gc_dwarf, mec='k', ms=6)
    ax.grid(ls=':')
    ax.legend(loc=4)
    plt.tight_layout()
    plt.savefig(lvdb_path + 'paper_examples/aitoff_mw_gc.pdf')
    pdf.savefig()
    plt.close()

    print("figure 3 finished")

    ### Figure 4

    fig, ax = plt.subplots(1,2, figsize=(10,5))
    ax[0].plot(dwarf_all['sg_xx']/1000., dwarf_all['sg_yy']/1000., '.',)

    ax[1].plot(dwarf_all['sg_yy']/1000., dwarf_all['sg_zz']/1000., '.')
    # ax[0].legend(loc=2)
    ax[0].set_xlabel(r'${\rm SGX~(Mpc)}$')
    ax[0].set_ylabel(r'${\rm SGY~(Mpc)}$')
    ax[1].set_xlabel(r'${\rm SGY~(Mpc)}$')
    ax[1].set_ylabel(r'${\rm SGZ~(Mpc)}$')

    ## this makes the axis limits have the same min/max and can see structure
    xmax = np.max([abs(ax[0].get_xlim()[0]),abs(ax[0].get_xlim()[1]), abs(ax[0].get_ylim()[0]), abs(ax[0].get_ylim()[1])])
    ax[0].set_xlim(-xmax, xmax)
    ax[0].set_ylim(-xmax, xmax)

    xmax = np.max([abs(ax[1].get_xlim()[0]),abs(ax[1].get_xlim()[1])])
    ymax = np.max([abs(ax[1].get_ylim()[0]),abs(ax[1].get_ylim()[1])])
    ax[1].set_xlim(-xmax, xmax)
    ax[1].set_ylim(-ymax, ymax)

    plt.tight_layout()
    pdf.savefig()
    plt.close()


    keep = np.zeros(len(dwarf_all), dtype=bool)
    for i in range(len(dwarf_all)):
        if np.ma.is_masked(dwarf_all['host'][i])==False:
            keep[i] = True
        else:
            keep[i] = False
    close_has_host = dwarf_all[keep]
    isolated_close = dwarf_all[~keep]
    print("figure4: host/isolated numbers", len(close_has_host), len(dwarf_all), len(isolated_close))

    fig, ax = plt.subplots(1,2, figsize=(10,5))
    # ax[0].plot(close_has_host['sg_xx']/1000., close_has_host['sg_yy']/1000., '.',  c='tab:blue')
    ax[0].plot(isolated_close['sg_xx']/1000., isolated_close['sg_yy']/1000., '.', c='tab:orange')

    # ax[1].plot(close_has_host['sg_yy']/1000., close_has_host['sg_zz']/1000., '.', c='tab:blue')
    ax[1].plot(isolated_close['sg_yy']/1000., isolated_close['sg_zz']/1000., '.', c='tab:orange')
    # ax[0].legend(loc=2)
    ax[0].set_xlabel(r'${\rm SGX~(Mpc)}$')
    ax[0].set_ylabel(r'${\rm SGY~(Mpc)}$')
    ax[1].set_xlabel(r'${\rm SGY~(Mpc)}$')
    ax[1].set_ylabel(r'${\rm SGZ~(Mpc)}$')

    for host in Counter(close_has_host['host']).keys():
        satellite = close_has_host[close_has_host['host']==host]
        if len(satellite)<5:
            ax[0].plot(satellite['sg_xx']/1000., satellite['sg_yy']/1000., '.', c='tab:orange')
            ax[1].plot(satellite['sg_xx']/1000., satellite['sg_zz']/1000., '.', c='tab:orange')
        else:
            ax[0].plot(satellite['sg_xx']/1000., satellite['sg_yy']/1000.,)
            ax[1].plot(satellite['sg_xx']/1000., satellite['sg_zz']/1000.,)

    ## this makes the axis limits have the same min/max and can see structure
    xmax = np.max([abs(ax[0].get_xlim()[0]),abs(ax[0].get_xlim()[1]), abs(ax[0].get_ylim()[0]), abs(ax[0].get_ylim()[1])])
    ax[0].set_xlim(-xmax, xmax)
    ax[0].set_ylim(-xmax, xmax)

    xmax = np.max([abs(ax[1].get_xlim()[0]),abs(ax[1].get_xlim()[1])])
    ymax = np.max([abs(ax[1].get_ylim()[0]),abs(ax[1].get_ylim()[1])])
    ax[1].set_xlim(-xmax, xmax)
    ax[1].set_ylim(-ymax, ymax)

    plt.tight_layout()
    # plt.savefig('supergalactic_overview.pdf')
    pdf.savefig()
    plt.close()
    print("Figure 4 finished")

    ### Figure 5

    ## radial velocity of M31 in GSR
    rv_gsr_m31 = -113.66

    fig, ax = plt.subplots(1,3, figsize=(15,5))
    ax[0].errorbar( dsph_mw_extra['distance_gc'], dsph_mw_extra['v3d']*dsph_mw_extra['vrad']/abs(dsph_mw_extra['vrad']), fmt='o', yerr=dsph_mw_extra['v3d_error'], xerr=[dsph_mw_extra['distance_gc_em'], dsph_mw_extra['distance_gc_ep']],mec = 'k', mew=.75)

    ## for simplicity adds same escape velocity and virial radius line for MW and M31
    rad = np.arange(1,500,1)
    for i in [0,1, 2]:
        ax[i].plot(rad, vesc(MWPotential2014,rad/8.)*220., c='k')
        ax[i].plot(rad, -vesc(MWPotential2014,rad/8.)*220., c='k')
        ax[i].set_xlim(1,500)
        ax[i].set_ylim(-700,700)
        ax[i].axvline(300, c='k',ls=':')

    ax[0].set_xlabel(r'$d_{\rm {\rm GC}}~({\rm kpc})$')
    ax[0].set_ylabel(r'$v_{\rm {\rm 3D}}\times v_{\rm rad}/ | v_{\rm rad} |~({\rm km~s^{-1}})$')

    ax[1].set_xlabel(r'$d_{\rm {\rm GC}}~({\rm kpc})$')
    ax[1].set_ylabel(r'$\sqrt{3}~v_{\rm {\rm gsr}}~({\rm km~s^{-1}})$')

    ax[1].errorbar( dsph_mw_extra['distance_gc'], np.sqrt(3)*dsph_mw_extra['vgsr'], fmt='o', yerr=[np.sqrt(3)*dsph_mw_extra['vgsr_em'],np.sqrt(3)*dsph_mw_extra['vgsr_ep']], xerr= [dsph_mw_extra['distance_gc_em'], dsph_mw_extra['distance_gc_ep']],mec = 'k', mew=.75)

    ## velocity relative to M31
    ax[2].errorbar(dwarf_all['distance_m31'], np.sqrt(3)*(dwarf_all['velocity_gsr'] - rv_gsr_m31), fmt='o', yerr=(dwarf_all['vlos_systemic_em']+dwarf_all['vlos_systemic_ep'])/2., mec='k')

    ax[2].set_xlabel(r'${\rm d_{\rm M31}~(Mpc)}$')
    ax[2].set_ylabel(r'$ \sqrt{3}~ v_{\rm M31}~({\rm km~s^{-1}})$')
    ax[2].set_xlim(0,500)
    ax[2].set_ylim(-500, 500)

    plt.tight_layout()
    pdf.savefig()
    plt.savefig(lvdb_path + 'paper_examples/rad_velocity_3panel.pdf')
    plt.close()
    print("figure 5 finished ")

    ### Figure 6
    ## simple Hubble flow model
    ## numbers should match Penarrubia et al 2014  
    ## https://ui.adsabs.harvard.edu/abs/2014MNRAS.443.2204P/abstract 
    def velocity_relation(r, A = 76.5, B = 46.6 ):
        return  A*r - B/np.sqrt(r) 
    Gravity = 4.301e-6

    plt.plot(dwarf_all['distance_lg']/1e3, np.sqrt(3.)*dwarf_all['velocity_lg'] , 'o', mec='k')
    plt.xlabel(r'${\rm d_{\rm LG}~(Mpc)}$')
    plt.ylabel(r'$ \sqrt{3}~ v_{\rm LG}~({\rm km~s^{-1}})$')

    plt.axvline(1.06, c='k',ls=':')

    plt.text((1.06*1.1)/10.5, .9, r'$R_{LG}$', fontsize=20, horizontalalignment='left', transform=plt.gca().transAxes)

    rad2 = np.arange(.1,3,.1)
    plt.plot(rad2, +vesc(KeplerPotential(amp=Gravity*5e12, normalize=False),rad2*1000.), c='k')
    plt.plot(rad2, -vesc(KeplerPotential(amp=Gravity*5e12, normalize=False),rad2*1000.), c='k')
    rad3 = np.arange(1, 12.1, .1)
    plt.plot(rad3, velocity_relation(rad3)*np.sqrt(3), c='k')

    plt.tight_layout()
    pdf.savefig()

    plt.xlim(0.1, 10.)
    plt.savefig(lvdb_path + 'paper_examples/rad_velocity_LG2.pdf')
    plt.close()
    print("figure 6 finished")

    ### Figure 7 
    mp.rcParams['ytick.right'] = False

    def const_mu(muV, rhalf):
        return muV - 36.57 - 2.5 * np.log10(2.*np.pi*rhalf**2)
    x = np.arange(1, 1e4, 1)
    for mu in [24, 26, 28, 30, 32]:
        plt.plot(x,  [const_mu(mu, i/1000.) for i in x], c='k', lw=2, ls=':')

    ## plot data
    plt.errorbar(dsph_mw['rhalf_sph_physical'], dsph_mw['M_V'], fmt='o', label=r'${\rm Dwarf~MW}$', c=color_dsph_mw, mec='k', mew=0.75)

    plt.plot(gc_ambiguous['rhalf_sph_physical'], gc_ambiguous['M_V'], 'o', c=color_gc_ufcss, mec='k', mew=0.75,label=r'${\rm Ambiguous/HFCSS}$')
    plt.plot(gc_mw_new['rhalf_sph_physical'], gc_mw_new['M_V'], 's', c=color_gc_disk, mec='k', mew=0.75,label=r'${\rm GC~New~Bulge/Disk/Halo}$')
    plt.plot(gc_harris['rhalf_sph_physical'], gc_harris['M_V'], 's', c=color_gc_harris, mec='k', mew=0.75,label=r'${\rm GC~Harris}$', ms=6)


    plt.gca().set_xscale('log') 
    plt.gca().invert_yaxis()
    plt.gca().set_xlabel(r'$R_{1/2}~({\rm pc})$')
    plt.gca().set_ylabel(r'$M_V$')
    # plt.gca().invert_yaxis()

    plt.ylim(2.5, -20)
    # plt.xlim(1, 4e3)
    plt.legend(loc=2)

    ## functions for twin y-axis
    def lum(x):
        m_x_sun=4.83
        return -.4*(x - m_x_sun) + np.log10(2.)
    def lum_inverse(x):
        m_x_sun=4.83
        return m_x_sun - (x )/0.4 + np.log10(2.)

    secax = plt.gca().secondary_yaxis('right', functions=(lum, lum_inverse))
    secax.set_ylabel(r'$\log_{10}{(M_{\star}/M_{\odot})}$', fontsize=15)
    # secax.invert_yaxis()

    plt.tight_layout()
    pdf.savefig()
    plt.xlim(1, 4e3)
    plt.savefig(lvdb_path + 'paper_examples/mw_rhalf_m_v_m_star.pdf')
    plt.close()
    

    plt.errorbar(dsph_mw['distance'], dsph_mw['M_V'], fmt='o', label=r'${\rm Dwarf~MW}$', c=color_dsph_mw, mec='k', mew=0.75)

    plt.plot(gc_ambiguous['distance'], gc_ambiguous['M_V'], 'o', c=color_gc_ufcss, mec='k', mew=0.75,label=r'${\rm Ambiguous/HFCSS}$')
    plt.plot(gc_mw_new['distance'], gc_mw_new['M_V'], 's', c=color_gc_disk, mec='k', mew=0.75,label=r'${\rm GC~New~Bule/Disk/Halo}$')
    plt.plot(gc_harris['distance'], gc_harris['M_V'], 's', c=color_gc_harris, mec='k', mew=0.75,label=r'${\rm GC~Harris}$')
    plt.gca().set_xscale('log') 
    plt.gca().invert_yaxis()

    plt.gca().set_xlabel(r'$d_{\odot}~({\rm kpc})$')
    plt.gca().set_ylabel(r'$M_V$')

    # plt.gca().invert_yaxis()
    plt.ylim(2.5, -20)
    # plt.legend(loc=2)

    def lum(x):
        m_x_sun=4.83
        return -.4*(x - m_x_sun) + np.log10(2.)


    def lum_inverse(x):
        m_x_sun=4.83
        return m_x_sun - (x )/0.4 + np.log10(2.)

    secax = plt.gca().secondary_yaxis('right', functions=(lum, lum_inverse))

    secax.set_ylabel(r'$\log_{10}{(M_{\star}/M_{\odot})}$', fontsize=15)
    plt.tight_layout()

    pdf.savefig()
    plt.savefig(lvdb_path + 'paper_examples/mw_distance_m_v_m_star.pdf')
    plt.close()
    print ("Figure 7 finished")

    #### Figure 8 
    x = np.arange(1, 1e4, 1)
    for mu in [24, 26, 28, 30, 32]:
        plt.plot(x,  [const_mu(mu, i/1000.) for i in x], c='k', lw=2, ls=':')

    ## plot data
    plt.errorbar(dsph_mw['rhalf_sph_physical'], dsph_mw['M_V'], fmt='o', label=r'${\rm Dwarf~MW}$', c=color_dsph_mw, mec='k', mew=0.75)
    plt.plot(dsph_m31['rhalf_sph_physical'], dsph_m31['M_V'], 'o', label=r'${\rm Dwarf~M31}$', c=color_dsph_m31, mec='k', mew=0.75,)
    plt.plot(dsph_lf['rhalf_sph_physical'],dsph_lf['M_V'], 'o', label=label_dsph_lf, c=color_dsph_lf, mec='k', mew=0.75,)
    plt.plot(dsph_lf_distant['rhalf_sph_physical'],dsph_lf_distant['M_V'], 'o', label=label_dsph_lf_distant, c=color_dsph_lf_distant, mec='k', mew=0.75,)

    plt.gca().set_xscale('log')
    plt.gca().invert_yaxis()
    plt.gca().set_xlabel(r'$R_{1/2}~({\rm pc})$')
    plt.gca().set_ylabel(r'$M_V$')

    plt.ylim(-0.1, -20)
    # plt.xlim(10, 5e3)
    # plt.gca().invert_yaxis()

    plt.legend(loc=2)

    ## twin y-axis
    def lum(x):
        m_x_sun=4.83
        return -.4*(x - m_x_sun) + np.log10(2.)
    def lum_inverse(x):
        m_x_sun=4.83
        return m_x_sun - (x )/0.4 + np.log10(2.)

    secax = plt.gca().secondary_yaxis('right', functions=(lum, lum_inverse))
    secax.set_ylabel(r'$\log_{10}{(M_{\star}/M_{\odot})}$', fontsize=15)
    plt.tight_layout()
    pdf.savefig()
    plt.xlim(10, 5e3)
    plt.savefig(lvdb_path + 'paper_examples/dwarf_rhalf_m_v_m_star.pdf')
    plt.close()

    x = np.arange(1, 1e4, 1)
    for mu in [24, 26, 28, 30, 32]:
        plt.plot(x,  [const_mu(mu, i/1000.) for i in x], c='k', lw=2, ls=':')

    ## plot data
    plt.errorbar(dsph_mw['rhalf_sph_physical'], dsph_mw['M_V'], fmt='o', label=r'${\rm Dwarf~MW}$', c=color_dsph_mw, mec='k', mew=0.75)
    plt.plot(dsph_m31['rhalf_sph_physical'], dsph_m31['M_V'], 'o', label=r'${\rm Dwarf~M31}$', c=color_dsph_m31, mec='k', mew=0.75,)
    plt.plot(dsph_lf['rhalf_sph_physical'],dsph_lf['M_V'], 'o', label=label_dsph_lf, c=color_dsph_lf, mec='k', mew=0.75,)
    # plt.plot(dsph_lf_distant['rhalf_sph_physical'],dsph_lf_distant['M_V'], 'o', label=label_dsph_lf_distant, c=color_dsph_lf_distant, mec='k', mew=0.75,)

    plt.gca().set_xscale('log')
    plt.gca().invert_yaxis()
    plt.gca().set_xlabel(r'$R_{1/2}~({\rm pc})$')
    plt.gca().set_ylabel(r'$M_V$')

    # plt.ylim(-0.1, -20)
    # 
    plt.legend(loc=2)

    ## twin y-axis
    def lum(x):
        m_x_sun=4.83
        return -.4*(x - m_x_sun) + np.log10(2.)
    def lum_inverse(x):
        m_x_sun=4.83
        return m_x_sun - (x )/0.4 + np.log10(2.)

    secax = plt.gca().secondary_yaxis('right', functions=(lum, lum_inverse))
    secax.set_ylabel(r'$\log_{10}{(M_{\star}/M_{\odot})}$', fontsize=15)
    plt.tight_layout()
    pdf.savefig()
    # plt.xlim(10, 5e3)
    # plt.savefig('dwarf_rhalf_m_v_m_star_LF.pdf')
    plt.close()


    x = np.arange(1, 1e4, 1)
    for mu in [24, 26, 28, 30, 32]:
        plt.plot(x,  [const_mu(mu, i/1000.) for i in x], c='k', lw=2, ls=':')

    ## plot data
    # plt.errorbar(dsph_mw['rhalf_sph_physical'], dsph_mw['M_V'], fmt='o', label=r'${\rm Dwarf~MW}$', c=color_dsph_mw, mec='k', mew=0.75)
    # plt.plot(dsph_m31['rhalf_sph_physical'], dsph_m31['M_V'], 'o', label=r'${\rm Dwarf~M31}$', c=color_dsph_m31, mec='k', mew=0.75,)
    # plt.plot(dsph_lf['rhalf_sph_physical'],dsph_lf['M_V'], 'o', label=label_dsph_lf, c=color_dsph_lf, mec='k', mew=0.75,)
    plt.plot(dsph_lf_distant['rhalf_sph_physical'],dsph_lf_distant['M_V'], 'o', label=label_dsph_lf_distant, c=color_dsph_lf_distant, mec='k', mew=0.75,)

    plt.gca().set_xscale('log')
    plt.gca().invert_yaxis()
    plt.gca().set_xlabel(r'$R_{1/2}~({\rm pc})$')
    plt.gca().set_ylabel(r'$M_V$')

    # plt.ylim(-0.1, -20)
    # plt.xlim(10, 5e3)
    plt.legend(loc=2)

    ## twin y-axis
    def lum(x):
        m_x_sun=4.83
        return -.4*(x - m_x_sun) + np.log10(2.)
    def lum_inverse(x):
        m_x_sun=4.83
        return m_x_sun - (x )/0.4 + np.log10(2.)

    secax = plt.gca().secondary_yaxis('right', functions=(lum, lum_inverse))
    secax.set_ylabel(r'$\log_{10}{(M_{\star}/M_{\odot})}$', fontsize=15)
    plt.tight_layout()
    pdf.savefig()
    plt.close()

    mp.rcParams['ytick.right'] = True
    print ("Figure 8 finished")

    ### Figure 9

    fig, ax = plt.subplots(2,3, figsize=(18, 12))
    def apy_lum(m_x, m_x_sun=4.83):
        return pow(10., -0.4*(m_x - m_x_sun) )
    for obj, color, label in zip([dsph_mw, dsph_m31, dsph_lf, gc_ambiguous], [color_dsph_mw,color_dsph_m31, color_dsph_lf, color_gc_ufcss], [label_dsph_mw, label_dsph_m31, label_dsph_lf, label_gc_ufcss]):
    # for obj, color, label in zip([dsph_mw, dsph_m31, dsph_lf, ], [color_dsph_mw,color_dsph_m31, color_dsph_lf, ], [label_dsph_mw, label_dsph_m31, label_dsph_lf, ]):
        ax[0][1].errorbar(obj['M_V'], obj['vlos_sigma'], fmt='o', yerr=[obj['vlos_sigma_em'], obj['vlos_sigma_ep']], xerr=[obj['M_V_em'], obj['M_V_ep']], c=color, label=label,  mec='k', mew=0.75)
        ax[0][1].errorbar(obj['M_V'], obj['vlos_sigma_ul'], fmt='o', yerr=obj['vlos_sigma_ul']*.4, uplims=True, c=color, xerr=[obj['M_V_em'], obj['M_V_ep']],  mec='k', mew=0.75)
        
        ax[0][0].errorbar(obj['rhalf_sph_physical'], obj['vlos_sigma'], fmt='o', yerr=[obj['vlos_sigma_em'], obj['vlos_sigma_ep']], xerr=[obj['rhalf_sph_physical_em'], obj['rhalf_sph_physical_ep']], c=color,  mec='k', mew=0.75)
        ax[0][0].errorbar(obj['rhalf_sph_physical'], obj['vlos_sigma_ul'], fmt='o', yerr=obj['vlos_sigma_ul']*.4, uplims=True, c=color, xerr=[obj['rhalf_sph_physical_em'], obj['rhalf_sph_physical_ep']],  mec='k', mew=0.75 )
        
        ax[1][0].errorbar(obj['M_V'], obj['mass_dynamical_wolf'], fmt='o', yerr=[obj['mass_dynamical_wolf_em'], obj['mass_dynamical_wolf_ep']], xerr=[obj['M_V_em'], obj['M_V_ep']], c=color,  mec='k', mew=0.75)
        ax[1][0].errorbar(obj['M_V'], obj['mass_dynamical_wolf_ul'], fmt='o', yerr=obj['mass_dynamical_wolf_ul']*.4, uplims=True, c=color, xerr=[obj['M_V_em'], obj['M_V_ep']],  mec='k', mew=0.75 )
        
        ax[1][2].errorbar(obj['M_V'], obj['mass_dynamical_wolf']/apy_lum(obj['M_V']), fmt='o', yerr=[obj['mass_dynamical_wolf_em']/apy_lum(obj['M_V']), obj['mass_dynamical_wolf_ep']/apy_lum(obj['M_V'])], xerr=[obj['M_V_em'], obj['M_V_ep']], c=color,  mec='k', mew=0.75)
        ax[1][2].errorbar(obj['M_V'], obj['mass_dynamical_wolf_ul']/apy_lum(obj['M_V']), fmt='o', yerr=obj['mass_dynamical_wolf_ul']/apy_lum(obj['M_V'])*.4, uplims=True, c=color, xerr=[obj['M_V_em'], obj['M_V_ep']],  mec='k', mew=0.75 )
        
        ax[0][2].errorbar(obj['rhalf_sph_physical'], obj['density_dyn_mcmc'], fmt='o', yerr=[obj['density_dyn_mcmc_em'], obj['density_dyn_mcmc_ep']], xerr=[obj['rhalf_sph_physical_em'], obj['rhalf_sph_physical_ep']], c=color,  mec='k', mew=0.75)
        ax[0][2].errorbar(obj['rhalf_sph_physical'], obj['density_dyn_mcmc_ul'], fmt='o', yerr=obj['density_dyn_mcmc_ul']*.4, uplims=True, c=color, xerr=[obj['rhalf_sph_physical_em'], obj['rhalf_sph_physical_ep']],  mec='k', mew=0.75 )
        
        ax[1][1].errorbar(obj['rhalf_sph_physical'], obj['mass_dynamical_wolf'], fmt='o', yerr=[obj['mass_dynamical_wolf_em'], obj['mass_dynamical_wolf_ep']], xerr=[obj['rhalf_sph_physical_em'], obj['rhalf_sph_physical_ep']], c=color,  mec='k', mew=0.75)
        ax[1][1].errorbar(obj['rhalf_sph_physical'], obj['mass_dynamical_wolf_ul'], fmt='o', yerr=obj['mass_dynamical_wolf_ul']*.4, uplims=True, c=color, xerr=[obj['rhalf_sph_physical_em'], obj['rhalf_sph_physical_ep']],  mec='k', mew=0.75 )
        

    ax[0][0].set_xscale('log')
    ax[0][0].set_ylim(0.5, 45)
    ax[0][1].set_yscale('log')
    ax[0][0].set_yscale('log')
    ax[0][1].set_ylim(0.5, 45)
    ax[0][1].set_xlim(3, -20)
    ax[1][0].set_xlim(3, -20)
    ax[0][0].set_xlim(1, 4e3)

    ax[0][0].set_xlabel(r'$R_{1/2}~({\rm pc})$')
    ax[0][1].set_xlabel(r'$M_V$')
    ax[0][1].set_ylabel(r'$\sigma_{\rm los}~({\rm km~s^{-1}})$')
    ax[0][0].set_ylabel(r'$\sigma_{\rm los}~({\rm km~s^{-1}})$')
    ax[0][1].legend(loc=4)

    ax[1][0].set_yscale('log')
    ax[1][0].set_xlabel(r'$M_V$')
    ax[1][0].set_ylabel(r'$M_{\rm dyn}~(r=r_{1/2})~({\rm M_{\odot}})$')

    ax[1][2].set_yscale('log')
    ax[1][2].set_xlabel(r'$M_V$')
    ax[1][2].set_ylabel(r'$M_{\rm dyn}/L~(r=r_{1/2})~({\rm M_{\odot}~L_{\odot}^{-1}})$')
    ax[1][2].invert_xaxis()

    ax[0][2].set_xscale('log')
    ax[0][2].set_yscale('log')
    ax[0][2].set_xlabel(r'$R_{1/2}~({\rm pc})$')
    ax[0][2].set_ylabel(r'$\rho_{\rm dyn}~(r=r_{1/2})~({\rm M_{\odot}~pc^{-3}})$')

    ax[1][1].set_xscale('log')
    ax[1][1].set_yscale('log')
    ax[1][1].set_xlabel(r'$R_{1/2}~({\rm pc})$')
    ax[1][1].set_ylabel(r'$M_{\rm dyn}~(r=r_{1/2})~({\rm M_{\odot}})$')


    plt.tight_layout()
    plt.savefig(lvdb_path + 'paper_examples/stellar_kinematics_6panal.pdf')
    pdf.savefig()
    plt.close()
    print("figure 9 finished")

    ### Figure 10
    def lum(m_x, m_x_sun=4.83):
        return pow(10., -0.4*(m_x - m_x_sun) )
    fig, ax = plt.subplots(1,2, figsize=(12+3, 6))
    ## stellar mass/luminosity - stellar metallicity 

    ## for the legend
    ax[0].plot([], [], 'o', mfc='None', mec='k', lw=3, zorder=1000, label=r'${\rm Spectroscopic~[Fe/H]}$')

    for obj, color, label, marker in zip([dsph_mw,  gc_ambiguous, gc_mw_new, gc_harris], [color_dsph_mw, color_gc_ufcss,color_gc_disk,color_gc_harris ], [label_dsph_mw, label_gc_ufcss, label_gc_disk, label_gc_harris], ['o', 'o', 's', 's']):
        ax[0].plot(obj['M_V'], obj['metallicity'], marker, c=color, lw=3, zorder=500,label=label)
        obj2= obj[obj['metallicity_type']=='spectroscopic']
        ax[0].plot(obj2['M_V'], obj2['metallicity'], marker, mfc='None', mec='k', lw=2, zorder=1000)
        
        ax[1].plot(obj['M_V'], obj['metallicity_spectroscopic_sigma'], marker, c=color, lw=2, zorder=500)
        obj2= obj[obj['metallicity_type']=='spectroscopic']
        ax[1].plot(obj2['M_V'], obj2['metallicity_spectroscopic_sigma'], marker, mfc='None', mec='k', lw=2, zorder=1000)
        ax[1].errorbar(obj2['M_V'], obj2['metallicity_spectroscopic_sigma_ul'], fmt='_',yerr=.05, uplims=True, c=color, )
    ax[0].invert_xaxis()
    ax[0].set_ylabel(r'${\rm [Fe/H]}$')

    ## luminoisty-metallicity relation from Simon 2019
    ## https://ui.adsabs.harvard.edu/abs/2019ARA%26A..57..375S/abstract
    x = np.arange( -20,3, .1)
    ax[0].plot(x,-1.68 + 0.29 * np.log10(lum(x)/1e6), c='k', lw=2, label=r'${\rm Simon~2019}$', zorder=1)

    ax[0].set_xlabel(r'$M_V$')
    ax[1].set_xlabel(r'$M_V$')
    
    ax[1].set_ylabel(r'$\sigma_{\rm [Fe/H]}$')
    ax[0].legend(loc=(.65, 0),)

    
    pdf.savefig()
    ax[0].set_xlim(1, -15)
    ax[1].set_xlim(1, -15)
    plt.savefig(lvdb_path + 'paper_examples/metallicity_overview_mw.pdf')
    plt.close()
    print("figure 10 finished")
    ### Figure 11

    fig, ax = plt.subplots(1,2, figsize=(12+3, 6))
    ## mass/luminosity- metallicity 

    ax[0].plot([], [], 'o', mfc='None', mec='k', lw=3, zorder=1000, label=r'${\rm Spectroscopic~[Fe/H]}$')
    for obj, color, label  in zip([dsph_mw,  dsph_m31, dsph_lf, dsph_lf_distant], [color_dsph_mw, color_dsph_m31,color_dsph_lf,color_dsph_lf_distant ], [label_dsph_mw, label_dsph_m31, label_dsph_lf, label_dsph_lf_distant]):
        ax[0].plot(obj['M_V'], obj['metallicity'], 'o', c=color, lw=2, zorder=500, label=label)
        obj2= obj[obj['metallicity_type']=='spectroscopic']
        ax[0].plot(obj2['M_V'], obj2['metallicity'], 'o', mfc='None', mec='k', lw=2, zorder=1000)
        
        ax[1].plot(obj['M_V'], obj['metallicity_spectroscopic_sigma'], 'o', c=color, lw=3, zorder=500)
        obj2= obj[obj['metallicity_type']=='spectroscopic']
        ax[1].plot(obj2['M_V'], obj2['metallicity_spectroscopic_sigma'], 'o', mfc='None', mec='k', lw=2, zorder=1000)
        ax[1].errorbar(obj2['M_V'], obj2['metallicity_spectroscopic_sigma_ul'], fmt='_',yerr=.05, uplims=True, c=color, )
    ax[0].invert_xaxis()
    ax[0].set_ylabel(r'${\rm [Fe/H]}$')

    ax[0].legend(loc=2)
    ## luminoisty-metallicity relation from Simon 2019
    ## https://ui.adsabs.harvard.edu/abs/2019ARA%26A..57..375S/abstract
    x = np.arange( -20,3, .1)
    ax[0].plot(x,-1.68 + 0.29 * np.log10(lum(x)/1e6), c='k', lw=2, label=r'${\rm Simon~2019}$', zorder=1)
    ax[0].set_xlabel(r'$M_V$')
    ax[1].set_xlabel(r'$M_V$')
    # ax[0].legend(loc=(1,0))
    
    ax[1].set_ylabel(r'$\sigma_{\rm [Fe/H]}$')

    plt.tight_layout()
    pdf.savefig()
    ax[0].set_xlim(1, -17)
    ax[1].set_xlim(1, -17)
    plt.savefig(lvdb_path + 'paper_examples/metallicity_overview_dwarf.pdf')
    plt.close()
    print("figure 11 finished")
    ### Figure 12

    keep = np.zeros(len(gc_dwarf), dtype=bool)
    for i in range(len(gc_dwarf)):
        if gc_dwarf['host'][i] in ['lmc', 'smc']:
            keep[i]=True
        else:
            keep[i]=False
    lmc_smc = gc_dwarf[keep]
    gc_dwarf2 = gc_dwarf[~keep]

    keep = np.zeros(len(gc_dwarf), dtype=bool)
    for i in range(len(gc_dwarf)):
        if gc_dwarf['host'][i] not in ['lmc', 'smc'] :
            keep[i]=True
        else:
            keep[i]=False
    gc_dwarf_mw = gc_dwarf[keep]
    print("figure 12, number in LMC/SMC and othw dwarf galaxies", len(gc_dwarf_mw), len(lmc_smc), len(gc_dwarf2), len(gc_dwarf))

    fig, ax = plt.subplots(2,2, figsize=(16*.8, 16*.8))

    for obj, color, label in zip([gc_harris2, gc_mw_new, gc_ambiguous, gc_dwarf2, gc_dwarf_sgr, lmc_smc], [color_gc_harris,color_gc_disk,  color_gc_ufcss, color_gc_dwarf,color_gc_dwarf, color_gc_lmc_smc], [label_gc_harris, label_gc_disk, label_gc_ufcss, label_gc_dwarf, '', label_gc_lmc_smc]):
        ax[0][0].errorbar( obj['rhalf_sph_physical'], obj['M_V'],fmt='s',  c=color, label=label,  mec='k', mew=0.75)
        
        ax[0][1].errorbar( obj['M_V'], obj['metallicity'],fmt='s',  c=color, label=label,  mec='k', mew=0.75)
        
        ax[1][0].errorbar( obj['rhalf_sph_physical'], obj['metallicity'],fmt='s',  c=color, label=label,  mec='k', mew=0.75)

        ax[1][1].errorbar( obj['age'], obj['metallicity'],fmt='s',  c=color, label=label,  mec='k', mew=0.75)
    
    ax[0][0].set_xlabel(r'$R_{1/2}~({\rm pc})$')
    ax[0][0].set_ylabel(r'$M_V$')

    ax[0][1].set_ylabel(r'${\rm [Fe/H]}$')
    ax[0][1].set_xlabel(r'$M_V$')
   
    ax[1][0].set_ylabel(r'${\rm [Fe/H]}$')
    ax[1][0].set_xlabel(r'$R_{1/2}~({\rm pc})$')

    ax[1][1].set_ylabel(r'${\rm [Fe/H]}$')
    ax[1][1].set_xlabel(r'${\rm Age}~({\rm Gyr})$')

    plt.tight_layout()
    ax[0][0].legend(loc=4)
    pdf.savefig()
    ax[0][0].set_xlim(0,20)
    ax[0][0].set_ylim(3,-12)
    ax[0][1].set_xlim(3,-12)
    ax[0][1].set_ylim(-3, .5)
    ax[1][0].set_xlim(0,20)
    ax[1][0].set_ylim(-3, .5)
    ax[1][1].set_ylim(-3, .5)
    plt.savefig(lvdb_path + 'paper_examples/gc_summary.pdf')
    plt.close()
    print("figure 12 finished")

    ### Figure 13

    dwarf_all['distance_GC_M31'] = np.zeros(len(dwarf_all), dtype=float)
    for i in range(len(dwarf_all)):
        dwarf_all['distance_GC_M31'][i] = np.min([dwarf_all['distance_gc'][i], dwarf_all['distance_m31'][i]])
    dwarf_all['ratio_HI_mstar_ul'] = 235600*dwarf_all['flux_HI_ul']*(dwarf_all['distance']/1000.)**2/10**dwarf_all['mass_stellar']
    comb_lg = dwarf_all[dwarf_all['distance']<3000]

    plt.errorbar(comb_lg['distance_GC_M31'], 10**comb_lg['mass_HI']/10**comb_lg['mass_stellar'], fmt='o', label=r'${\rm Satellite}$', c='tab:blue', zorder=1000)
    plt.errorbar(comb_lg['distance_GC_M31'],comb_lg['ratio_HI_mstar_ul'], fmt='.', uplims=True, yerr=comb_lg['ratio_HI_mstar_ul']/2., c='tab:orange')

    plt.gca().set_yscale('log')

    plt.gca().set_ylabel(r'$M_{H I}/M_{\star}$')
    plt.gca().set_xlabel(r'$D_{\rm min(MW,~M31)}~({\rm kpc})$')

    plt.tight_layout()
    plt.savefig(lvdb_path + 'paper_examples/distance_ratio_gas_stellar.pdf')
    pdf.savefig()
    plt.close()
    print("figure 13 finished")

    print("distance vs M_V")
    plt.figure(figsize = (6.4*2, 4.8*2))
    plt.plot(dwarf_all['distance']/1e3, dwarf_all['M_V'], 'o')
    max_x = plt.gca().get_xlim()[1]
    min_x=plt.gca().get_xlim()[0]
    plt.gca().set_xlim(min_x, max_x)
    x = np.arange(0.1, max_x, .1)
    plt.plot(x, 20. - lvdb.dist_mod_kpc(x*1e3), c='k',ls=':')
    plt.plot(x, 19. - lvdb.dist_mod_kpc(x*1e3), c='k',ls=':')
    plt.plot(x, 21. - lvdb.dist_mod_kpc(x*1e3), c='k', ls=':')
    
    plt.gca().invert_yaxis()
    plt.gca().set_xlabel(r'$d~({\rm Mpc})$')
    plt.gca().set_ylabel(r'$M_V$')
    plt.tight_layout()
    pdf.savefig()
    plt.close()

    plt.figure(figsize = (6.4*2, 4.8*2))
    plt.gca().set_xlim(min_x, max_x)
    plt.plot(x, 20. - lvdb.dist_mod_kpc(x*1e3), c='k',ls=':')
    plt.plot(x, 19. - lvdb.dist_mod_kpc(x*1e3), c='k',ls=':')
    plt.plot(x, 21. - lvdb.dist_mod_kpc(x*1e3), c='k', ls=':', label=r'$V=19,20,21~{\rm mag}$')
    plt.gca().invert_yaxis()
    plt.gca().set_xlabel(r'$d~({\rm Mpc})$')
    plt.gca().set_ylabel(r'$M_V$')
    plt.plot(isolated_close['distance']/1e3, isolated_close['M_V'], 'o', c='gray', label=r'${\rm isolated}$')
    for host in Counter(close_has_host['host']).keys():
        satellite = close_has_host[close_has_host['host']==host]
        if len(satellite)<5:
            plt.plot(satellite['distance']/1e3, satellite['M_V'], 'o', c='k')
        else:
            plt.plot(satellite['distance']/1000., satellite['M_V'],'o')
    plt.plot([],[], 'o', c='k', label=r'${\rm satellite,~N<5}$')
    plt.plot([],[], 'o', c='tab:blue', label=r'${\rm satellite,~N>5~(colored~points)}$')
    plt.legend(loc=4)
    plt.tight_layout()
    pdf.savefig()
    plt.close()

    plt.figure(figsize = (6.4*2, 4.8*2))
    plt.gca().set_xlim(min_x, max_x)
    plt.plot(x, 20. - lvdb.dist_mod_kpc(x*1e3), c='k',ls=':')
    plt.plot(x, 19. - lvdb.dist_mod_kpc(x*1e3), c='k',ls=':')
    plt.plot(x, 21. - lvdb.dist_mod_kpc(x*1e3), c='k', ls=':')
    plt.gca().invert_yaxis()
    plt.gca().set_xlabel(r'$d~({\rm Mpc})$')
    plt.gca().set_ylabel(r'$M_V$')
    plt.plot(isolated_close['distance']/1e3, isolated_close['M_V'], '.', c='gray')
    for host in Counter(close_has_host['host']).keys():
        satellite = close_has_host[close_has_host['host']==host]
        if len(satellite)<5:
            plt.plot(satellite['distance']/1e3, satellite['M_V'], '.', c='k')
        else:
            plt.plot(satellite['distance']/1000., satellite['M_V'],)
    plt.tight_layout()
    pdf.savefig()
    plt.close()
    print("distance vs M_V finished")

    print("surface brightness plots")

    plt.figure(figsize = (6.4*2, 4.8*2))

    plt.gca().set_ylabel(r'$\mu~({\rm mag~arcsec^{-2}})$')
    plt.gca().set_xlabel(r'$M_V$')
    plt.plot(dwarf_all['M_V'], dwarf_all['surface_brightness_rhalf'], 'o')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    pdf.savefig()
    plt.close()

    ## second plots
    plt.figure(figsize = (6.4*2, 4.8*2))
    plt.gca().set_ylabel(r'$\mu~({\rm mag~arcsec^{-2}})$')
    plt.gca().set_xlabel(r'$d~({\rm Mpc})$')
    plt.plot(dwarf_all['distance']/1e3, dwarf_all['surface_brightness_rhalf'], 'o')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    pdf.savefig()
    plt.close() 

    print("surface brightness plots finished")