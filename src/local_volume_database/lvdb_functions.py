import yaml
import astropy.coordinates as coord
from astropy import units as u
import astropy.table as table
import os.path
import numpy.ma as ma
import numpy as np

import matplotlib.pyplot as plt
import matplotlib._color_data as mcd
from collections import Counter
import corner



## load environment variable
lvdb_path = os.environ.get('LVDBDIR')
# print(lvdb_path)

plt.style.use(lvdb_path+'/code/std.mplstyle')
import matplotlib as mp
mp.rcParams['text.usetex'] = True

def get_notes(key, **kwargs):
    ## input is lvdb key
    ## outputs notes in LVDB YAML file
    print_output = kwargs.get("print_output", True)
    path = kwargs.get('path', lvdb_path + '/data_input/')
    with open(path+ key + '.yaml', 'r') as stream:
        try:
            stream_yaml = yaml.load(stream, Loader=yaml.Loader)
            if 'notes' in stream_yaml.keys():
                if print_output:
                    print("notes:\t", key)
    #                 print(stream_yaml['notes'])
                    for i in range(len(stream_yaml['notes'])):
                        print(stream_yaml['notes'][i])
                return stream_yaml['notes']
            else:
                if print_output:
                    print("no notes", key)
        except:
            if print_output:
                print("no key: ", key)

coord.galactocentric_frame_defaults.set('v4.0')
gc_frame = coord.Galactocentric()

def add_coord(table_input):
    ## add Galactic coordinates to table
    c_table_input = coord.SkyCoord(ra=table_input['ra']*u.deg, dec=table_input['dec']*u.deg,  frame='icrs',)
    table_input['ll'] = c_table_input.galactic.l.value
    table_input['bb'] = c_table_input.galactic.b.value

def latex_name(name):
    ## adds some special latex characters to some systems 
    name = name.replace('_', '\\_')
    if len(name)<6:
        return name
    if name[:6] == 'Bootes':
        name = name.replace('oo', 'o\\"{o}')
    if len(name)<5:
        return name
    if name[:5] == 'Munoz':
        name = name.replace('n', '{\\~n}')
    return name

def add_column(table, yaml_key, yaml_name, **kwargs):
    """
    adds columns to the table.  
    Parameters
    ----------
    table : input/output table
    yaml_key : overall key  in YAML file 
    yaml_name : name in the nested key strucutre in the YAML file
    col_type (optional) : python column type, default is float
    table_yaml_name (optional) : if the table output needs to have the name of the column changed from `yaml_name`
    """
    table_yaml_name = kwargs.get('table_yaml_name', yaml_name)
    col_type = kwargs.get('col_type', float)
    table[table_yaml_name] = np.ma.masked_all(len(table), dtype=col_type)
    #np.zeros(len(table), dtype=col_type)
    path = kwargs.get('path', lvdb_path + '/data_input/')
    for i in range(len(table)):
        k = table['key'][i]
        with open(path+ k +'.yaml', 'r') as stream:
            try:
                stream_yaml = yaml.load(stream, Loader=yaml.Loader)

                if yaml_key in stream_yaml.keys() and  yaml_name in stream_yaml[yaml_key].keys():
                    table[table_yaml_name][i] = stream_yaml[yaml_key][yaml_name]
                # else:
                #     print(stream_yaml['key'], "discovery_year")
            except yaml.YAMLError as exc:
                print(exc)
    # return table[table_yaml_name]
                
def make_latex_value(value, em, ep, **kwargs):
    ## for rounding when making latex table
    ul = kwargs.get('ul', False)
    n = kwargs.get('n', 2)
    str_out = ' '
    if ul!=False:
        str_out = '$<' + "{:0.{prec}f}".format(ul, prec=n) +'$'
    elif ma.is_masked(value)==False:
#         if ma.is_masked(value)==False:
#             return str_out
        str_out = '$'+"{:0.{prec}f}".format(value, prec=n)
        if ma.is_masked(em)==False and ma.is_masked(ep)== False :
            if "{:0.{prec}f}".format(em, prec=n) == "{:0.{prec}f}".format(ep, prec=n):
                str_out+='\\pm'+"{:0.{prec}f}".format(em, prec=n)
            else:
                str_out+='_{-'+"{:0.{prec}f}".format(em,prec=n) + '}^{+' +"{:0.{prec}f}".format(ep,prec=n) +'}'
        str_out+='$'
    return str_out

def add_year(table, **kwargs):
    ## initial version of add_column()
    table['year'] = np.zeros(len(table), dtype=int)
    path = kwargs.get('path', lvdb_path +'data_input/')
    for i in range(len(table)):
        k = table['key'][i]
        with open(path+ k +'.yaml', 'r') as stream:
            try:
                stream_yaml = yaml.load(stream, Loader=yaml.Loader)

                if 'discovery_year' in stream_yaml['name_discovery'].keys():
                    table['year'][i] = stream_yaml['name_discovery']['discovery_year']
                # else:
                #     print(stream_yaml['key'], "discovery_year")
            except yaml.YAMLError as exc:
                print(exc)
    # return table['year']

reference_column = ['ref_structure', 'ref_distance', 'ref_vlos', 'ref_proper_motion', 'ref_metallicity_spectroscopic', 'ref_structure_king', 'ref_age', 'ref_structure_sersic', 'ref_metallicity_isochrone', 'ref_flux_HI']
def get_citations(systems, style='individual', reference_column=reference_column):
    ## output references in citep{}
    if style=='individual':
        for i in range(len(systems)):
            cite_temp = []
            for j in reference_column:
                
                if ma.is_masked(systems[j][i])==False:
                    cite_temp.append(systems[j][i])
            ## unique entries
            cite_temp2 = np.unique(cite_temp)
            out =''
            for temp in cite_temp2:
                out+=temp+', '
#         x=x.replace('&', '\string&')
            print(systems['name'][i], "\\citep{"+out[:-2]+"}" +',')
    if style=='all':
        cite_temp = []
        for i in range(len(systems)):
            
            for j in reference_column:
                
                if ma.is_masked(systems[j][i])==False:
                    cite_temp.append(systems[j][i])
            ## unique entries
        cite_temp2 = np.unique(cite_temp)
        out =''
        for temp in cite_temp2:
            out+=temp+', '
#         x=x.replace('&', '\string&')
        print("\\citep{"+out[:-2]+"}")

def lum(m_x, m_x_sun=4.83):
    ## M_V -> L_V, default is V-band (absolute magnitude to luminosity conversion)
	return pow(10., -.4*(m_x - m_x_sun) )
def lum_inverse(lum, m_x_sun=4.83 ):
    ## Luminosity to absoltue magnitude conversion, default is V-band
    return np.log10(lum)/(-.4)+m_x_sun
def dm(x):
    ## distance modulus equation
    return pow(10., x/5.+1.)/1000.
def dist_mod(mu, mu_e=0, mu_em=0, mu_ep=0):
    ## distance modulus to distance, has option to mcmc over errors
    ## output is in kpc
	# def dm(x):
	# 	return pow(10., x/5.+1.)/1000.

	if mu_e > 0:
		sample = []
		for i in range(100000):
			x = np.random.normal(mu, mu_e)
			sample.append(dm(x))
		a = corner.quantile(np.array(sample), [.5, .1587, .8413])

		return (a[0], a[0]-a[1], a[2]-a[0], (a[2]-a[1])/2.)
	else:
		return dm(mu)
def dist_mod_kpc(dist):
    return (np.log10(dist*1000.) -1.)*5.

pm_data = table.Table.read(lvdb_path+'/data/pm_overview.csv')

def plot_proper_motion_galaxy_3panel(key, pm_overview = pm_data,  add=[], **kwargs):
    #key = LVDB key, system to examine
    ## pm_overview = input proper motion array.  Default is to load the local LVDB version but this option allows for other input
    ## add = by hand point to include (new measurement to compare to, )
    # add = [pmra, pmra_error, pmdec, pmdec_error] (for plotting only)
    error_option = kwargs.get('error_option', False)
    save_file_name = kwargs.get('save_file_name', '')
    non_gaia = kwargs.get('non_gaia',False)## include 3rd plot with non-Gaia based proper motion measurements
    
    if non_gaia:
        fig, ax = plt.subplots(1,3,figsize=(18,5))
    else:
        fig, ax = plt.subplots(1,2,figsize=(12,5))

    pm_overview2 = pm_overview[pm_overview['key']==key]
    print("number of measurements and mehtods:",len(pm_overview2), Counter(pm_overview2['method']))
    print_pm = kwargs.get('print_pm',True)
    if print_pm:
        for kk in range(len(pm_overview2)):
            print(pm_overview2['ref_cite'][kk], pm_overview2['method'][kk], pm_overview2['pmra'][kk], pm_overview2['pmra_em'][kk], pm_overview2['pmdec'][kk], pm_overview2['pmdec_em'][kk])
    keep = np.zeros(len(pm_overview2), dtype=bool)
    exclude = kwargs.get('exclude', [])
    for i in range(len(pm_overview2)):
        if pm_overview2['ref_cite'][i] in exclude:
            keep[i]=False
        else:
            keep[i]=True
    print("excluded measurements:",len(pm_overview2)-np.sum(keep), len(exclude))
    pm_overview2 = pm_overview2[keep]
    
    temp_color = list(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    for i in ['black', 'lightgreen', 'magenta', 'navy', 'gold', 'maroon']:
        temp_color.append(i)
    tot = len(list(mcd.XKCD_COLORS))
    x = np.random.choice(tot, tot)
    for x2 in x:
        temp_color.append(list(mcd.XKCD_COLORS)[x2])
    xkcd = 0

    if len(add)>0:
        for jj in range(len(add)):
            for i in [0,1]:
                ax[i].errorbar(add[jj][0], add[jj][3], fmt='*',xerr=[[add[jj][1]],[add[jj][2]]], yerr=[[add[jj][4]],[add[jj][5]]],label=add[jj][6], c=temp_color[xkcd], ms=10, zorder=1000)
            if non_gaia:
                ax[2].errorbar(add[jj][0], add[jj][3], fmt='*',  xerr=[[add[jj][1]],[add[jj][2]]], yerr=[[add[jj][4]],[add[jj][5]]], label=add[jj][6], c=temp_color[xkcd], ms=10, zorder=1000)
            xkcd+=1
        
    for kk in range(len(pm_overview2)):
        if non_gaia:
            ax[2].errorbar(pm_overview2['pmra'][kk], pm_overview2['pmdec'][kk], fmt='o', xerr=[[pm_overview2['pmra_em'][kk]],[pm_overview2['pmra_ep'][kk]]],yerr=[[pm_overview2['pmdec_em'][kk]], [pm_overview2['pmdec_ep'][kk]]], label=pm_overview2['citation'][kk],c=temp_color[xkcd])
            
        if pm_overview2['method'][kk] in ['GAIA_EDR3', 'GAIA_DR2']:
            ax[1].errorbar(pm_overview2['pmra'][kk], pm_overview2['pmdec'][kk], fmt='o', xerr=[[pm_overview2['pmra_em'][kk]], [pm_overview2['pmra_ep'][kk]]], yerr=[[pm_overview2['pmdec_em'][kk]], [pm_overview2['pmdec_ep'][kk]]], label=pm_overview2['citation'][kk],c=temp_color[xkcd])
            if pm_overview2['method'][kk]=='GAIA_EDR3':
                ax[0].errorbar(pm_overview2['pmra'][kk], pm_overview2['pmdec'][kk], fmt='o', xerr=[[pm_overview2['pmra_em'][kk]], [pm_overview2['pmra_ep'][kk]]],  yerr=[[pm_overview2['pmdec_em'][kk]], [pm_overview2['pmdec_ep'][kk]]], label=pm_overview2['citation'][kk],c=temp_color[xkcd])
        xkcd+=1
    
#     ax[3].axis("off")
    for i in [0,1]:
        ax[i].set_xlabel(r'$\mu_{\alpha *}~({\rm mas~yr^{-1}})$', fontsize=15)
        ax[i].set_ylabel(r'$\mu_{\delta }~({\rm mas~yr^{-1}})$', fontsize=15)
        xvals,yvals = ax[i].axes.get_xlim(),ax[i].axes.get_ylim()

        xrange = xvals[1]-xvals[0]
        yrange = yvals[1]-yvals[0]
        ax[i].set_aspect((xrange/yrange)*.9, adjustable='box')
#         ax[i].set_box_aspect(1)
    ax[0].set_title(r'${\rm {\it Gaia}~EDR3}$')
    ax[1].set_title(r'${\rm All~{\it Gaia}}$')
    if non_gaia:
        ax[2].set_xlabel(r'$\mu_{\alpha *}~({\rm mas~yr^{-1}})$', fontsize=15)
        ax[2].set_ylabel(r'$\mu_{\delta }~({\rm mas~yr^{-1}})$', fontsize=15)
        ax[2].set_title(r'${\rm All~Measurements}$')
        xvals,yvals = ax[2].axes.get_xlim(),ax[2].axes.get_ylim()

        xrange = xvals[1]-xvals[0]
        yrange = yvals[1]-yvals[0]
        ax[2].set_aspect((xrange/yrange)*.9, adjustable='box')
#     ax.set_box_aspect(1)
    
    if non_gaia:
        lgd = ax[2].legend(fontsize=13, loc=(1,0))
    else:
        lgd = ax[1].legend(fontsize=13, loc=(1,0))
        
    plt.tight_layout()
    if len(save_file_name)>0:
        plt.savefig(save_file_name)
        
    plt.show()