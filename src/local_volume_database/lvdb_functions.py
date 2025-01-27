import yaml
import astropy.coordinates as coord
from astropy import units as u
import astropy.table as table
import os.path
import numpy.ma as ma
import numpy as np

## load environment variable
lvdb_path = os.environ.get('LVDBDIR')
# print(lvdb_path)

def get_notes(key, **kwargs):
    ## input is lvdb key
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
	return pow(10., -.4*(m_x - m_x_sun) )
def lum_inverse(lum, m_x_sun=4.83 ):
    return np.log10(lum)/(-.4)+m_x_sun
def dm(x):
		return pow(10., x/5.+1.)/1000.
def dist_mod(mu, mu_e=0, mu_em=0, mu_ep=0):
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