import corner
import numpy as np
import matplotlib.pyplot as plt
# %matplotlib inline

import astropy.table as table

from collections import Counter

import numpy.ma as ma
import yaml

from astropy import units as u
import astropy.coordinates as coord

## this makes a->z plus aa -> zz for labeling citations
import string
string.ascii_lowercase
long_list = list(string.ascii_lowercase)
for i in list(string.ascii_lowercase):
    for j in list(string.ascii_lowercase):
        long_list.append(i+j)

import os.path

path = 'data_input/'

def make_latex_value(value, em, ep, **kwargs):
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
def add_year(table):
    table['year'] = np.zeros(len(table), dtype=int)
    path = '/Users/apace/Documents/local_volume_database/data_input/'
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

def add_column(table, yaml_key, yaml_name, **kwargs):
    table_yaml_name = kwargs.get('table_yaml_name', yaml_name)
    col_type = kwargs.get('col_type', 'float')
    table[table_yaml_name] = np.ma.masked_all(len(table), dtype=col_type)
    #np.zeros(len(table), dtype=col_type)
    path = '/Users/apace/Documents/local_volume_database/data_input/'
    for i in range(len(table)):
        k = table['key'][i]
        with open(path+ k +'.yaml', 'r') as stream:
            try:
                stream_yaml = yaml.load(stream, Loader=yaml.Loader)

                if yaml_key in stream_yaml.keys() and  table_yaml_name in stream_yaml[yaml_key].keys():
                    table[table_yaml_name][i] = stream_yaml[yaml_key][table_yaml_name]
                # else:
                #     print(stream_yaml['key'], "discovery_year")
            except yaml.YAMLError as exc:
                print(exc)
    # return table[table_yaml_name]

def latex_name(name):
    if len(name)<6:
        return name
    if name[:6] == 'Bootes':
        name = name.replace('oo', 'o\\"{o}')
    if len(name)<5:
        return name
    if name[:5] == 'Munoz':
        name = name.replace('n', '{\\~n}')
    return name

coord.galactocentric_frame_defaults.set('v4.0')
gc_frame = coord.Galactocentric()

def add_coord(table_input):
    c_table_input = coord.SkyCoord(ra=table_input['ra']*u.deg, dec=table_input['dec']*u.deg,  frame='icrs',)
    table_input['ll'] = c_table_input.galactic.l.value
    table_input['bb'] = c_table_input.galactic.b.value

def create_latex_table_name_discovery(output, input_table, **kwargs):
    classification_column = kwargs.get("classification_column", 'confirmed_dwarf')
    classification_output = kwargs.get("classification_output", 'Dwarf Galaxy')
    missing_host = []
    with open(output, 'w+') as f:
        for i in range(len(input_table)):
            end_line = '\\\\'
            if i == len(input_table)-1:
                end_line=''
            k = input_table['key'][i]
            name = latex_name(input_table['name'][i])

            out_str = '' + name + ' & '
            other_name = []
            ref = []
            host = ''
            
            with open(path+ k +'.yaml', 'r') as stream:
                try:
                    stream_yaml = yaml.load(stream, Loader=yaml.Loader)
        #             out_str = ''
                    if 'other_name' in stream_yaml['name_discovery'].keys():
                        for other_name_list in range(len(stream_yaml['name_discovery']['other_name'])):
                            other_name.append(stream_yaml['name_discovery']['other_name'][other_name_list])
        #             else:
        #                 print(stream_yaml['key'], "missing table")
                    if 'ref_discovery' in stream_yaml['name_discovery'].keys():
                        for other_name_list in range(len(stream_yaml['name_discovery']['ref_discovery'])):
                            x = stream_yaml['name_discovery']['ref_discovery'][other_name_list]
                            x=x.replace('&', '\string&')
                            ref.append(x)
                    if 'host' in stream_yaml['name_discovery'].keys():
                        host_key = stream_yaml['name_discovery']['host']
                        if host_key in ['MW', 'LF']:
                            host = host_key
                        else:
                            if os.path.isfile(path + host_key +'.yaml'):  
                                with open(path + host_key +'.yaml', 'r') as stream:
                                    stream_yaml_key = yaml.load(stream, Loader=yaml.Loader)
                                    host = stream_yaml_key['name_discovery']['name']
                            else:
                                host = host_key
                                host=host.replace('_', ' ')
                                missing_host.append(stream_yaml['key'])
                                print('no host key', host_key)
                except yaml.YAMLError as exc:
                    print(exc)
            # 
            c = coord.SkyCoord(ra=input_table['ra'][i]*u.degree, dec=input_table['dec'][i]*u.degree)
            x = c.to_string('hmsdms')
            x1 = c.ra.to_string(unit=u.hourangle, sep=":", precision=1, alwayssign=False, pad=True)
            x2 = c.dec.to_string(sep=":", precision=1, alwayssign=True, pad=True)
            place =0
            if len(other_name)>0:
                out_str += other_name[place] + ' & '
            else:
                out_str +=  ' & '
            out_str += x1 + ' & ' + x2  + ' & '
            out_str += host + ' & '
            if len(ref)>0:
                out_str += "\\citet{" +ref[place] +'}' + ' & '
            else:
                out_str +=  ' & '
            place +=1
        #     print(k,  x1, x2)
    #         print( out_str+ ' \\\\')
            if input_table['confirmed_real'][i]==1:
                out_str +=  '  & '
            else:
                out_str +=  ' Cand. & '

            if input_table[classification_column][i]==1:
                out_str +=  classification_output
            elif input_table[classification_column][i]==0:
                out_str +=  '  '        

            # if :
            out_str += end_line+'\n'
            # else:
            #     out_str += '\n'
            f.write(out_str )
            while len(other_name)>place or len(ref) > place:
                out_str2 = ' & '
                if len(other_name)>place:
                    name = latex_name(other_name[place])
                    out_str2 +=  name
                out_str2 += ' &&&& '
                if len(ref)>place:
                    out_str2 += "\\citet{" + ref[place]+'}'
                out_str2 += end_line
                place+=1
                f.write( out_str2+'\n')
    print(Counter(missing_host))

def create_latex_table_structure(output, output_citations, input_table, **kwargs):
    citations = []
    letter = []
    spatial_units = kwargs.get("spatial_units", "arcmin")
    spatial_units_conversion = 60.
    spatial_convert_factor = kwargs.get("spatial_convert_factor", 1.) ## to be used to convert csv values from arcmin to arsec
    round_rhalf = 1
    round_distance = 1
    if spatial_units == 'arcsec':
        spatial_units_conversion = kwargs.get("spatial_units_conversion", 3600.) 
        round_rhalf = 0
        round_distance = 0
    with open(output, 'w+') as f:
        for i in range(len(input_table)):
            letter_to_list = []

            ## this adds combines all the citations per object that this table is using
            cite_temp = []
            if ma.is_masked(input_table['ref_structure'][i])==False:
                cite_temp.append(input_table['ref_structure'][i])
            if ma.is_masked(input_table['ref_distance'][i])==False:
                cite_temp.append(input_table['ref_distance'][i])
            if ma.is_masked(input_table['ref_m_v'][i])==False:
                cite_temp.append(input_table['ref_m_v'][i])

            ## unique entries
            cite_temp2 = np.unique(cite_temp)
#             print(input_table['key'][i], cite_temp2)
            ## this checks if a citation has already been used and pulls it, otherwise it finds the next letter to assign to a citation 
            for tt in cite_temp2:
                if  isinstance(tt,str)==False:
                    continue
#                 if not tt:
#                     continue
#                 if len(tt)<5:
#                     continue
                if tt in  citations:
                    letter_to_list.append(letter[citations.index(tt)])
                else:
                    citations.append(tt)
                    letter_to_list.append(long_list[len(letter)])
                    letter.append(long_list[len(letter)])

            letter_to_list_string = ""
            if len(letter_to_list)>0:
                for kk in letter_to_list:
                    letter_to_list_string+=kk  +','
                letter_to_list_string = letter_to_list_string[:-1]

            x = np.random.normal(input_table['rhalf'][i], (input_table['rhalf_em'][i]+input_table['rhalf_ep'][i])/2., 1000)
            if ma.is_masked(input_table['ellipticity_em'][i])==False and ma.is_masked(input_table['ellipticity'][i])==False:
                y = np.random.normal(input_table['ellipticity'][i], (input_table['ellipticity_em'][i]+input_table['ellipticity_ep'][i])/2., 1000)
            else:
                y = np.zeros(len(x))
            if ma.is_masked(input_table['distance_em'][i])==False and ma.is_masked(input_table['distance_ep'][i])==False:
                z = np.random.normal(input_table['distance'][i], (input_table['distance_em'][i]+input_table['distance_ep'][i])/2., 1000)
            else:
                z = np.full(1000, input_table['distance'][i])
            x2 = x[np.logical_and(y>=0, y<1)]
            y2 = y[np.logical_and(y>=0, y<1)]
            z2 = z[np.logical_and(y>=0, y<1)]
            comb = x2 *np.pi/180./spatial_units_conversion*1000.*np.sqrt(1. - y2)* z2
            comb2 = comb[~np.isnan(comb)]
            if len(comb2)==0:
                str_rhalf=''
                if ma.is_masked(input_table['rhalf_em'][i])==True or ma.is_masked(input_table['rhalf_ep'][i])==True and ma.is_masked(input_table['rhalf'][i])==False:
                    rh = input_table['distance'][i]*input_table['rhalf'][i]/180./spatial_units_conversion*1000.*np.pi
                    if ma.is_masked(input_table['ellipticity'][i])==False:
                        rh = rh*np.sqrt(1.-input_table['ellipticity'][i])
                    str_rhalf=make_latex_value(rh, input_table['rhalf_em'][i], input_table['rhalf_em'][i], n=round_rhalf)  
            else:
                mean = corner.quantile(comb2, [.5, .1587, .8413, 0.0227501, 0.97725])
                str_rhalf = make_latex_value(mean[0], mean[0]-mean[1], mean[2]-mean[0], n=round_rhalf)
        #     rhalf_pc = mean[0]
        #     rhalf_pc_error = (mean[2]-mean[1])/2.
            end_line = '\\\\'
            if i == len(input_table)-1:
                end_line=''
            ## output each row of our table, plus the citations at the end of the line
            name = latex_name(input_table['name'][i])
            f.write(name + '&'+"{:0.4f}".format(input_table['ra'][i])+'&'+"{:0.4f}".format(input_table['dec'][i])+'&'+
                  make_latex_value(spatial_convert_factor*input_table['rhalf'][i], spatial_convert_factor*input_table['rhalf_em'][i], spatial_convert_factor*input_table['rhalf_ep'][i], n=2)+'&'+
                 make_latex_value(input_table['ellipticity'][i], input_table['ellipticity_em'][i], input_table['ellipticity_ep'][i], ul=input_table['ellipticity_ul'][i], n=2)+ '& '+
                  make_latex_value(input_table['position_angle'][i], input_table['position_angle_em'][i], input_table['position_angle_ep'][i], n=1)+ '& '+
                  str_rhalf + '& '+
                  make_latex_value(input_table['distance_modulus'][i], input_table['distance_modulus_em'][i], input_table['distance_modulus_ep'][i], n=2)+' & '+
                 make_latex_value(input_table['distance'][i], input_table['distance_em'][i], input_table['distance_ep'][i], n=round_distance)+ ' & '+ make_latex_value(input_table['apparent_magnitude_v'][i], np.ma.masked, np.ma.masked, n=1)+ ' & '+
                   make_latex_value(input_table['M_V'][i], input_table['M_V_em'][i], input_table['M_V_ep'][i], n=1)+ ' & '+ letter_to_list_string+
                  end_line+'\n')
    with open(output_citations, 'w+') as f:
        for i,j in zip(letter, citations):
        #     print(i, j)
            j = j.replace('&', '\string&')
            f.write( "("+i+") \citet{"+j+"}\n",)

def create_latex_table_kinematics(output, output_citations, input_table, **kwargs):
    citations = []
    letter = []
    add_age = kwargs.get("add_age", False)
    with open(output, 'w+') as f:
        for i in range(len(input_table)):
    #         print(input_table['key'][i])
            letter_to_list = []

            ## this adds combines all the citations per object that this table is using
            cite_temp = []
            if ma.is_masked(input_table['ref_vlos'][i])==False:
                cite_temp.append(input_table['ref_vlos'][i])
            if ma.is_masked(input_table['ref_proper_motion'][i])==False:
                cite_temp.append(input_table['ref_proper_motion'][i])
            if ma.is_masked(input_table['ref_metallicity_spectroscopic'][i])==False:
                cite_temp.append(input_table['ref_metallicity_spectroscopic'][i])
            if add_age == True and ma.is_masked(input_table['ref_age'][i])==False:
                cite_temp.append(input_table['ref_age'][i])
            ## unique entries
    #         if not cite_temp:
    # #             print(input_table['key'][i])
    #             f.write(input_table['name'][i]+
    #               '\\\\'+'\n')
    #             continue

            cite_temp2 = np.unique(cite_temp)
    #         print(cite_temp, cite_temp2)
            ## this checks if a citation has already been used and pulls it, otherwise it finds the next letter to assign to a citation 
    #         if not cite_temp:
            for tt in cite_temp2:
                if  isinstance(tt,str)==False:
                    continue
                if tt in  citations:
                    letter_to_list.append(letter[citations.index(tt)])
                else:
                    citations.append(tt)
                    letter_to_list.append(long_list[len(letter)])
                    letter.append(long_list[len(letter)])

            letter_to_list_string = ""
            if len(letter_to_list)>0:
                for kk in letter_to_list:
                    letter_to_list_string+=kk  +','
                letter_to_list_string = letter_to_list_string[:-1]
            end_line = '\\\\'
            if i == len(input_table)-1:
                end_line=''
            name = latex_name(input_table['name'][i])
            age_str = ''
            if add_age:
                age_str = make_latex_value(input_table['age'][i], input_table['age_em'][i], input_table['age_ep'][i], n=1)+'&'
            f.write(name + '&'+"{:0.4f}".format(input_table['ll'][i])+'&'+"{:0.4f}".format(input_table['bb'][i])+'&'+ make_latex_value(input_table['vlos_systemic'][i], input_table['vlos_systemic_em'][i], input_table['vlos_systemic_ep'][i], n=1)+ '& '+
               make_latex_value(input_table['vlos_sigma'][i], input_table['vlos_sigma_em'][i], input_table['vlos_sigma_ep'][i], ul=input_table['vlos_sigma_ul'][i], n=2)+ '& '+
                    make_latex_value(input_table['metallicity_spectroscopic'][i], input_table['metallicity_spectroscopic_em'][i], input_table['metallicity_spectroscopic_ep'][i], n=2)+'&'+
                 make_latex_value(input_table['metallicity_spectroscopic_sigma'][i], input_table['metallicity_spectroscopic_sigma_em'][i], input_table['metallicity_spectroscopic_sigma_ep'][i], ul=input_table['metallicity_spectroscopic_sigma_ul'][i],  n=2)+ '& '+ age_str + 
                  make_latex_value(input_table['pmra'][i], input_table['pmra_em'][i], input_table['pmra_ep'][i], n=3)+'&'+
                 make_latex_value(input_table['pmdec'][i], input_table['pmdec_em'][i], input_table['pmdec_ep'][i],  n=3)+ '& '+
                  letter_to_list_string+
                  end_line+'\n')
    with open(output_citations, 'w+') as f:
        for i,j in zip(letter, citations):
        #     print(i, j)
            j = j.replace('&', '\string&')
            f.write( "("+i+") \citet{"+j+"}\n",)

def create_latex_table_mass(output, output_citations, input_table):
    citations = []
    letter = []

    with open(output, 'w+') as f:
        for i in range(len(input_table)):
            letter_to_list = []
            cite_temp = []
            if ma.is_masked(input_table['ref_m_v'][i])==False:
                cite_temp.append(input_table['ref_m_v'][i])
            if ma.is_masked(input_table['ref_vlos'][i])==False:
                cite_temp.append(input_table['ref_vlos'][i])
            if ma.is_masked(input_table['ref_structure'][i])==False:
                cite_temp.append(input_table['ref_structure'][i])
            if ma.is_masked(input_table['ref_flux_HI'][i])==False:
                cite_temp.append(input_table['ref_flux_HI'][i])


            cite_temp2 = np.unique(cite_temp)
            for tt in cite_temp2:
                if  isinstance(tt,str)==False:
                    continue
                if tt in  citations:
                    letter_to_list.append(letter[citations.index(tt)])
                else:
                    citations.append(tt)
                    letter_to_list.append(long_list[len(letter)])
                    letter.append(long_list[len(letter)])

            letter_to_list_string = ""
            if len(letter_to_list)>0:
                for kk in letter_to_list:
                    letter_to_list_string+=kk  +','
                letter_to_list_string = letter_to_list_string[:-1]

            # input_table
            def lum(m_x, m_x_sun=4.83):
                return pow(10., -.4*(m_x - m_x_sun) )
            if ma.is_masked(input_table['M_V'][i])==False:
                mstar = lum(input_table['M_V'][i]) * 2.
                mstar_str = '$'+"{:0.1e}".format(mstar)[:3] + '\\times 10^{'+str(int("{:0.1e}".format(mstar).split('e')[1]))+'}'+'$'
            else:
                mstar_str=''
        #.replace('+','')    
            mass_to_light_str = ''
            m_dyn_str = ''
            if ma.is_masked(input_table['vlos_sigma'][i])==False and ma.is_masked(input_table['rhalf_sph_physical'][i])==False:
                mdyn = 930 * input_table['rhalf_sph_physical'][i] * input_table['vlos_sigma'][i]**2
                m_dyn_str = '$'+"{:0.1e}".format(mdyn)[:3] + '\\times 10^{'+str(int("{:0.1e}".format(mdyn).split('e')[1]))+'}'+'$'
                mass_to_light = mdyn/(mstar/2.)
                mass_to_light_str = '$'+"{:0.1f}".format(mass_to_light)+'$'
            elif ma.is_masked(input_table['vlos_sigma_ul'][i])==False and ma.is_masked(input_table['rhalf_sph_physical'][i])==False:
                mdyn_ul = 930 * input_table['rhalf_sph_physical'][i] * input_table['vlos_sigma_ul'][i]**2
                m_dyn_str = '$<'+"{:0.1e}".format(mdyn_ul)[:3] + '\\times 10^{'+str(int("{:0.1e}".format(mdyn_ul).split('e')[1]))+'}'+'$'
                mass_to_light = mdyn/(mstar/2.)
                mass_to_light_str = '$<'+"{:0.1f}".format(mass_to_light)+'$'
        #         print(input_table['key'][i], m_dyn_str)
            else:
                m_dyn_str = ''
                mass_to_light_str=''
            m_HI_str=''
            m_HI_m_star_str=''
            if ma.is_masked(input_table['mass_HI'][i])==False:
                mdyn = 10.**input_table['mass_HI'][i]
                m_HI_str = '$'+"{:0.1e}".format(mdyn)[:3] + '\\times 10^{'+str(int("{:0.1e}".format(mdyn).split('e')[1]))+'}'+'$'
                m_HI_m_star = mdyn/mstar
                m_HI_m_star_str = '$'+"{:0.1f}".format(m_HI_m_star)+'$'
            elif ma.is_masked(input_table['mass_HI_ul'][i])==False:
                temp = 10.**input_table['mass_HI_ul'][i]
                m_HI_str = '$<'+"{:0.1e}".format(temp)[:3] + '\\times 10^{'+str(int("{:0.1e}".format(temp).split('e')[1]))+'}'+'$'
                m_HI_m_star = temp/mstar
                ul_str = abs(int("{:0.1e}".format(m_HI_m_star).split('e')[1]))
                fmt = "{:0."+str(ul_str)+"f}"
                m_HI_m_star_str = '$<'+ fmt.format(m_HI_m_star)+'$'
            else:
                m_HI_str = ''
                m_HI_m_star_str=''

            end_line = '\\\\'
            if i == len(input_table)-1:
                end_line=''
        #     mstar_s = "{:0.1e}".format(mstar)
        #     mstar_split = mstar_s.split('e')
        #     mstar_split[1]
            name = latex_name(input_table['name'][i])

            f.write(name +' & ' +  mstar_str +' & ' + m_dyn_str +' & ' + mass_to_light_str +' & ' + m_HI_str+' & ' + m_HI_m_star_str+ '& '+ letter_to_list_string + end_line + '\n')

    with open(output_citations, 'w+') as f:
        for i,j in zip(letter, citations):
        #     print(i, j)
            j = j.replace('&', '\string&')
            f.write( "("+i+") \citet{"+j+"}\n",)

dsph_mw = table.Table.read('data/dwarf_mw.csv')
dsph_m31 = table.Table.read('data/dwarf_m31.csv')
dsph_lf = table.Table.read('data/dwarf_local_field.csv')
dsph_lf_distant = table.Table.read('data/dwarf_local_field_distant.csv')
gc_ufsc = table.Table.read('data/gc_ufsc.csv')
gc_disk = table.Table.read('data/gc_disk.csv')
gc_harris = table.Table.read('data/gc_harris.csv')
gc_dwarf = table.Table.read('data/gc_dwarf_hosted.csv')

# dsph_mw['year'] = add_year(dsph_mw)
# dsph_m31['year'] = add_year(dsph_m31)
# dsph_lf['year'] = add_year(dsph_lf)
# dsph_lf_distant['year'] = add_year(dsph_lf_distant)

for tab in [dsph_mw, dsph_m31,  dsph_lf, dsph_lf_distant]:
    # tab['year'] = add_year(tab)
    add_year(tab)
    add_coord(tab)

# gc_ufsc['year'] = add_year(gc_ufsc)
# gc_disk['year'] = add_year(gc_disk)
# gc_harris['year'] = add_year(gc_harris)
# gc_dwarf['year'] = add_year(gc_dwarf)
for tab in [gc_ufsc, gc_disk,  gc_harris, gc_dwarf]:
    # tab['year'] = add_year(tab)
    add_year(tab)
    # tab['age'] = add_column(tab, 'age', 'age')
    # tab['age_em'] = add_column(tab, 'age', 'age_em')
    # tab['age_ep'] = add_column(tab, 'age', 'age_ep')
    # tab['ref_age'] = add_column(tab, 'age', 'ref_age', col_type='U100')
    add_column(tab, 'age', 'age')
    add_column(tab, 'age', 'age_em')
    add_column(tab, 'age', 'age_ep')
    add_column(tab, 'age', 'ref_age', col_type='U100')
    add_coord(tab)

# add_coord(dsph_m31)
# add_coord(dsph_mw)
# add_coord(dsph_lf)
# add_coord(dsph_lf_distant)

# add_coord(gc_ufsc)
# add_coord(gc_harris)
# add_coord(gc_disk)

# add_coord(gc_dwarf)

dsph_mw.sort(['year', 'key'])
create_latex_table_name_discovery('table/table_data/dwarf_mw_name_discovery_data.tex', dsph_mw)
dsph_mw2 = dsph_mw[dsph_mw['key']!='LMC']
dsph_mw2 = dsph_mw2[dsph_mw2['key']!='SMC']
dsph_mw2.sort('key')
create_latex_table_structure('table/table_data/dwarf_mw_structure_data.tex', 'table/table_data/dwarf_mw_structure_citations.tex', dsph_mw2)
create_latex_table_kinematics('table/table_data/dwarf_mw_kinematics_data.tex', 'table/table_data/dwarf_mw_kinematics_citations.tex', dsph_mw2)
create_latex_table_mass('table/table_data/dwarf_mw_mass_data.tex', 'table/table_data/dwarf_mw_mass_citations.tex', dsph_mw2)


dsph_m31.sort(['year', 'key'])
create_latex_table_name_discovery('table/table_data/dwarf_m31_name_discovery_data.tex', dsph_m31)
dsph_m31.sort('key')
create_latex_table_structure('table/table_data/dwarf_m31_structure_data.tex', 'table/table_data/dwarf_m31_structure_citations.tex', dsph_m31)
create_latex_table_kinematics('table/table_data/dwarf_m31_kinematics_data.tex', 'table/table_data/dwarf_m31_kinematics_citations.tex', dsph_m31)
create_latex_table_mass('table/table_data/dwarf_m31_mass_data.tex', 'table/table_data/dwarf_m31_mass_citations.tex', dsph_m31)

dsph_lf.sort(['year', 'key'])
create_latex_table_name_discovery('table/table_data/dwarf_lf_name_discovery_data.tex', dsph_lf)
dsph_lf.sort('key')
create_latex_table_structure('table/table_data/dwarf_lf_structure_data.tex', 'table/table_data/dwarf_lf_structure_citations.tex', dsph_lf)
create_latex_table_kinematics('table/table_data/dwarf_lf_kinematics_data.tex', 'table/table_data/dwarf_lf_kinematics_citations.tex', dsph_lf)
create_latex_table_mass('table/table_data/dwarf_lf_mass_data.tex', 'table/table_data/dwarf_lf_mass_citations.tex', dsph_lf)

dsph_lf_distant.sort(['year', 'key'])
create_latex_table_name_discovery('table/table_data/dwarf_lf_distant_name_discovery_data.tex', dsph_lf_distant)
dsph_lf_distant.sort(['host', 'key'])
create_latex_table_structure('table/table_data/dwarf_lf_distant_structure_data.tex', 'table/table_data/dwarf_lf_distant_structure_citations.tex', dsph_lf_distant, spatial_convert_factor=60, spatial_units='arcsec', spatial_units_conversion=60.)
create_latex_table_kinematics('table/table_data/dwarf_lf_distant_kinematics_data.tex', 'table/table_data/dwarf_lf_distant_kinematics_citations.tex', dsph_lf_distant)
create_latex_table_mass('table/table_data/dwarf_lf_distant_mass_data.tex', 'table/table_data/dwarf_lf_distant_mass_citations.tex', dsph_lf_distant)

gc_ufsc.sort(['year', 'key'])
create_latex_table_name_discovery('table/table_data/gc_ufsc_name_discovery_data.tex', gc_ufsc, classification_column='confirmed_star_cluster', classification_output='Star Cluster')
gc_ufsc.sort('key')
create_latex_table_structure('table/table_data/gc_ufsc_structure_data.tex', 'table/table_data/gc_ufsc_structure_citations.tex', gc_ufsc, )
create_latex_table_kinematics('table/table_data/gc_ufsc_kinematics_data.tex', 'table/table_data/gc_ufsc_kinematics_citations.tex', gc_ufsc, add_age=True)
# create_latex_table_mass('table/table_data/gc_ufsc_mass_data.tex', 'table/table_data/gc_ufsc_mass_citations.tex', gc_ufsc)

gc_disk.sort(['year', 'key'])
create_latex_table_name_discovery('table/table_data/gc_disk_name_discovery_data.tex', gc_disk, classification_column='confirmed_star_cluster', classification_output='Star Cluster')
gc_disk.sort('key')
create_latex_table_structure('table/table_data/gc_disk_structure_data.tex', 'table/table_data/gc_disk_structure_citations.tex', gc_disk)
create_latex_table_kinematics('table/table_data/gc_disk_kinematics_data.tex', 'table/table_data/gc_disk_kinematics_citations.tex', gc_disk, add_age=True)

gc_harris.sort(['year', 'key'])
create_latex_table_name_discovery('table/table_data/gc_harris_name_discovery_data.tex', gc_harris, classification_column='confirmed_star_cluster', classification_output='Star Cluster')
gc_harris.sort('key')
create_latex_table_structure('table/table_data/gc_harris_structure_data.tex', 'table/table_data/gc_harris_structure_citations.tex', gc_harris)
create_latex_table_kinematics('table/table_data/gc_harris_kinematics_data.tex', 'table/table_data/gc_harris_kinematics_citations.tex', gc_harris, add_age=True)

gc_dwarf.sort(['year', 'key'])
create_latex_table_name_discovery('table/table_data/gc_dwarf_name_discovery_data.tex', gc_dwarf, classification_column='confirmed_star_cluster', classification_output='Star Cluster')
gc_dwarf.sort(['host','key'])
create_latex_table_structure('table/table_data/gc_dwarf_structure_data.tex', 'table/table_data/gc_dwarf_structure_citations.tex', gc_dwarf, spatial_convert_factor=60, spatial_units='arcsec', spatial_units_conversion=60.)
create_latex_table_kinematics('table/table_data/gc_dwarf_kinematics_data.tex', 'table/table_data/gc_dwarf_kinematics_citations.tex', gc_dwarf, add_age=True)


dir_list = os.listdir(path)
dir_list = [i for i in dir_list if i!='readme.md']

table_list = []
for i in range(len(dir_list)):
    with open(path+ dir_list[i], 'r') as stream:
        try:
            stream_yaml = yaml.load(stream, Loader=yaml.Loader)
            if 'table' in stream_yaml.keys():
                if stream_yaml['table'] == 'candidate':

                    if 'discovery_year' in stream_yaml['name_discovery'].keys():
                        y = stream_yaml['name_discovery']['discovery_year']
                    else:
                        y=0
                    if 'name' in stream_yaml['name_discovery'].keys():
                        n = stream_yaml['name_discovery']['name']
                    else:
                        n=0
                    table_list.append((stream_yaml['key'], y, n))
        except yaml.YAMLError as exc:
            print(exc)
print("candidate list", len(table_list))

table_list_sort = sorted(table_list, key=lambda student: student[1])

output = 'table/table_data//candidate_name_discovery_data.tex'
with open(output, 'w+') as f:
    for candidate in range(len(table_list_sort)):
        yaml_name = table_list_sort[candidate][0]
#         print(yaml_name)
        end_line = '\\\\'
        if candidate == len(table_list_sort)-1:
            end_line=''
        k = yaml_name
        
        other_name = []
        ref = []
        ref_fp = []
        host = ''
        with open(path+ yaml_name +'.yaml', 'r') as stream:
            try:
                stream_yaml = yaml.load(stream, Loader=yaml.Loader)
                c = coord.SkyCoord(ra=stream_yaml['location']['ra']*u.degree, dec=stream_yaml['location']['dec']*u.degree)
                x = c.to_string('hmsdms')
                x1 = c.ra.to_string(unit=u.hourangle, sep=":", precision=1, alwayssign=False, pad=True)
                x2 = c.dec.to_string(sep=":", precision=1, alwayssign=True, pad=True)
                name = ''
                if 'name' in stream_yaml['name_discovery'].keys():
                    name= latex_name(stream_yaml['name_discovery']['name'])
    #             out_str = ''
                if 'other_name' in stream_yaml['name_discovery'].keys():
                    for other_name_list in range(len(stream_yaml['name_discovery']['other_name'])):
                        other_name.append(stream_yaml['name_discovery']['other_name'][other_name_list])
    #             else:
    #                 print(stream_yaml['key'], "missing table")
                if 'ref_discovery' in stream_yaml['name_discovery'].keys():
                    for other_name_list in range(len(stream_yaml['name_discovery']['ref_discovery'])):
                        x = str(stream_yaml['name_discovery']['ref_discovery'][other_name_list])
                        x=x.replace('&', '\string&')
                        ref.append(x)
                if 'ref_false_positive' in stream_yaml['name_discovery'].keys():
                    for other_name_list in range(len(stream_yaml['name_discovery']['ref_false_positive'])):
                        x = str(stream_yaml['name_discovery']['ref_false_positive'][other_name_list])
                        x=x.replace('&', '\string&')
                        ref_fp.append(x)
                if 'host' in stream_yaml['name_discovery'].keys():
                    host_key = stream_yaml['name_discovery']['host']
                    if host_key in ['MW', 'LF']:
                        host = host_key
                    else:
                        if os.path.isfile(path + host_key +'.yaml'):  
                            with open(path + host_key +'.yaml', 'r') as stream:
                                stream_yaml_key = yaml.load(stream, Loader=yaml.Loader)
                                host = stream_yaml_key['name_discovery']['name']
                        else:
                            host = host_key
                            host=host.replace('_', ' ')
                            # print('no host key', host_key)
                out_str = '' + name + ' & '
                place =0
                if len(other_name)>0:
                    out_str += other_name[place] + ' & '
                else:
                    out_str +=  ' & '
                out_str += x1 + ' & ' + x2  + ' & '
                out_str += host + ' & '
                if len(ref)>0:
                    out_str += "\\citet{" +ref[place] +'}' + ' & '
                else:
                    out_str +=  ' & '
                
            #     print(k,  x1, x2)
        #         print( out_str+ ' \\\\')
                fp =0
                if 'false_positive' in stream_yaml['name_discovery'].keys():
                    fp = stream_yaml['name_discovery']['false_positive']
#                     print(stream_yaml['name_discovery']['false_positive'])
                if fp==1:
                    out_str +=  ' FP & '
                else:
                    out_str +=  '  & '

#                 if stream_yaml_key[classification_column][i]==1:
#                     out_str +=  classification_output
#                 elif input_table[classification_column][i]==0:
#                     out_str +=  '  '        
                if len(ref_fp)>0:
                    out_str += "\\citet{" + ref_fp[place] +'}'

                place +=1
                out_str += end_line+'\n'
                # else:
                #     out_str += '\n'
                f.write(out_str )
                while len(other_name)>place or len(ref) > place or len(ref_fp)>place:
                    out_str2 = ' & '
                    if len(other_name)>place:
                        name = latex_name(other_name[place])
                        out_str2 +=  name
                    out_str2 += ' &&&& '
                    if len(ref)>place:
                        out_str2 += "\\citet{" + ref[place]+'}'
                    out_str2 += ' && '
                    if len(ref_fp)>place:
                        out_str2 += "\\citet{" + ref_fp[place] +'}'
                    out_str2 += end_line
                    place+=1
                    f.write( out_str2+'\n')
            except yaml.YAMLError as exc:
                    print(exc)

citations = []
letter = []
# 
output = 'table/table_data/candidate_structure_data.tex'
output_citations = 'table/table_data/candidate_structure_citations.tex'

with open(output, 'w+') as f:
    for candidate in range(len(table_list_sort)):
        letter_to_list = []
        yaml_name = table_list_sort[candidate][0]
        # print(yaml_name, path+ yaml_name +'.yaml')
        
        with open(path+ yaml_name +'.yaml', 'r') as stream:
            try:
                # yaml_name = table_list_sort[candidate][0]
                stream_yaml = yaml.load(stream, Loader=yaml.Loader)
                ## this adds combines all the citations per object that this table is using
                cite_temp = []
                if 'structure' in stream_yaml.keys() and 'ref_structure' in stream_yaml['structure'].keys():
                    if ma.is_masked(stream_yaml['structure']['ref_structure'])==False:
                        cite_temp.append(stream_yaml['structure']['ref_structure'])
                if 'm_v' in stream_yaml.keys() and 'ref_m_v' in stream_yaml['m_v'].keys():
                    if ma.is_masked(stream_yaml['m_v']['ref_m_v'])==False:
                        cite_temp.append(stream_yaml['m_v']['ref_m_v'])
                if 'distance' in stream_yaml.keys() and 'ref_distance' in stream_yaml['distance'].keys():
                    if ma.is_masked(stream_yaml['distance']['ref_distance'])==False:
                        cite_temp.append(stream_yaml['distance']['ref_distance'])

                ## unique entries
                cite_temp2 = np.unique(cite_temp)
        #             print(input_table['key'][i], cite_temp2)
                ## this checks if a citation has already been used and pulls it, otherwise it finds the next letter to assign to a citation 
                for tt in cite_temp2:
                    if  isinstance(tt,str)==False:
                        continue
        #                 if not tt:
        #                     continue
        #                 if len(tt)<5:
        #                     continue
                    if tt in  citations:
                        letter_to_list.append(letter[citations.index(tt)])
                    else:
                        citations.append(tt)
                        letter_to_list.append(long_list[len(letter)])
                        letter.append(long_list[len(letter)])

                letter_to_list_string = ""
                if len(letter_to_list)>0:
                    for kk in letter_to_list:
                        letter_to_list_string+=kk  +','
                    letter_to_list_string = letter_to_list_string[:-1]
                rh_str = ''
                str_rhalf = ''
                if 'structure' in stream_yaml.keys() and 'rhalf' in stream_yaml['structure'].keys():
                    rh_str = make_latex_value(stream_yaml['structure']['rhalf'], np.ma.masked,np.ma.masked, n=2)
                def dm(x):
                    return pow(10., x/5.+1.)/1000.
                dist_str = ''
                dm_str = ''
                if 'distance' in stream_yaml.keys() and 'distance' in stream_yaml.keys()  and 'distance_modulus' in stream_yaml['distance'].keys():
                    d = dm(stream_yaml['distance']['distance_modulus'])
                    dm_str = make_latex_value(stream_yaml['distance']['distance_modulus'], np.ma.masked,np.ma.masked, n=2)
                    dist_str = make_latex_value(d, np.ma.masked,np.ma.masked, n=1)
                if 'distance' in stream_yaml.keys() and 'structure' in stream_yaml.keys() and 'rhalf' in stream_yaml['structure'].keys() and 'distance_modulus' in stream_yaml['distance'].keys():
                    d = dm(stream_yaml['distance']['distance_modulus'])
                    rh = d*stream_yaml['structure']['rhalf']/180./60.*1000.*np.pi
                    str_rhalf = make_latex_value(rh, np.ma.masked,np.ma.masked, n=1)
                v_str = ''
                mv_str = ''
                if 'distance' in stream_yaml.keys() and 'distance_modulus' in stream_yaml['distance'].keys() and 'apparent_magnitude_v' in stream_yaml['m_v'].keys():
                    mv_str = make_latex_value(stream_yaml['m_v']['apparent_magnitude_v']-stream_yaml['distance']['distance_modulus'], np.ma.masked,np.ma.masked, n=1)
                if 'm_v' in stream_yaml.keys() and 'apparent_magnitude_v' in stream_yaml['m_v'].keys():
                    v_str = make_latex_value(stream_yaml['m_v']['apparent_magnitude_v'], np.ma.masked,np.ma.masked, n=1) 
                end_line = '\\\\'
                if candidate == len(table_list_sort)-1:
                    end_line=''
                ## output each row of our table, plus the citations at the end of the line
                name = latex_name(stream_yaml['name_discovery']['name'])
                f.write(name + '&'+"{:0.4f}".format(stream_yaml['location']['ra'])+'&'+"{:0.4f}".format(stream_yaml['location']['dec'])+'&'+
                      rh_str+'&'+
                      str_rhalf + '& '+
                      dm_str +' & '+
                     dist_str + ' & '+ v_str+ ' & '+
                       mv_str + ' & '+ letter_to_list_string+
                      end_line+'\n')
            except yaml.YAMLError as exc:
                print(exc)
with open(output_citations, 'w+') as f:
    for i,j in zip(letter, citations):
    #     print(i, j)
        j = j.replace('&', '\string&')
        f.write( "("+i+") \citet{"+j+"}\n",)
        