Tutorials \& Examples
=====================

Example: Accessing the database and creating a figure 
---------------------------------------------

This example creates a figure comparing half-light radius and absolute magnitude of dwarf galaxies and star clusters. 

Load required packages.

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    import astropy.table as table

load the data from github

.. code-block:: python

    ## can be loaded from the release page (recommended/perfered method). Note that the release tag needs to be selected.
    dsph_all = table.Table.read('https://github.com/apace7/local_volume_database/releases/download/v1.0.0/dwarf_all.csv')

    ## load from main github branch
    dsph_mw = table.Table.read('https://raw.githubusercontent.com/apace7/local_volume_database/main/data/dwarf_mw.csv')
    dsph_m31 = table.Table.read('https://raw.githubusercontent.com/apace7/local_volume_database/main/data/dwarf_m31.csv')
    dsph_lf = table.Table.read('https://raw.githubusercontent.com/apace7/local_volume_database/main/data/dwarf_local_field.csv')
    gc_halo = table.Table.read('https://raw.githubusercontent.com/apace7/local_volume_database/main/data/gc_ambiguous.csv')
    gc_disk = table.Table.read('https://raw.githubusercontent.com/apace7/local_volume_database/main/data/gc_mw_new.csv')
    gc_harris = table.Table.read('https://raw.githubusercontent.com/apace7/local_volume_database/main/data/gc_harris.csv')

make a figure

.. code-block:: python

    def const_mu(muV, rhalf):
        return muV - 36.57 - 2.5 * np.log10(2.*np.pi*rhalf**2)
    x = np.arange(1, 1e4, 1)
    for mu in [24, 26, 28, 30, 32]:
        plt.plot(x,  [const_mu(mu, i/1000.) for i in x], c='k', lw=2, ls=':')

    plt.errorbar(dsph_mw['rhalf_sph_physical'], dsph_mw['M_V'], fmt='o', label=r'${\rm Dwarf~MW}$', )
    plt.plot(dsph_m31['rhalf_sph_physical'], dsph_m31['M_V'], 'o', label=r'${\rm Dwarf~M31}$')
    plt.plot(dsph_lf['rhalf_sph_physical'],dsph_lf['M_V'], 'o', label=r'${\rm Dwarf~LF}$')
    plt.plot(gc_halo['rhalf_sph_physical'], gc_halo['M_V'], 'o',label=r'${\rm UFSC}$')
    plt.plot(gc_disk['rhalf_sph_physical'], gc_disk['M_V'], 'o',label=r'${\rm GC~New}$')
    plt.plot(gc_harris['rhalf_sph_physical'], gc_harris['M_V'], 'o',label=r'${\rm GC~Harris}$')
    plt.gca().set_xscale('log')
    plt.gca().invert_yaxis()
    plt.gca().set_xlabel(r'$r_{1/2}~({\rm pc})$')
    plt.gca().set_ylabel(r'$M_V$')
    plt.legend(loc=(2))
    plt.ylim(3, -20)
    plt.xlim(1, 4e3)
    plt.show()

Figure from the sample (this does include custom matplotlib stylefile and latex for the captions).
I note that the name for one type of system in the catalog (labeled UFSC = ultra-faint star clusters) has been changed several times. In the majority of the documentation these are referred to as, new star cluster-like systems in Galactic halo and/or ambiguous systems (as in the nature of dwarf galaxy and star cluster is ambiguous).

.. figure:: /figures/example_rhalf_MV.png
   :width: 60%
   :align: center
   :alt: Example

   Example figure

Example Jupyter Notebooks
---------------------------------------------

In addition to the example above, there are two folders with example ipython notebooks:

* `general example notebooks <https://github.com/apace7/local_volume_database/blob/main/example_notebooks/>`_ 
* `LVDB overview plots  <https://github.com/apace7/local_volume_database/tree/main/paper_examples>`_


Some Recommendations 
---------------------------------------------

For scientific analysis, I would recommendation fixing the LVDB version to a tagged GitHub release.
Alternatively, a specific commit can be useed instead of the latest github version.  

For example, this loads an older version of the data/dwarf_all.csv table for either a tagged release or a specific commit.

.. code-block:: python

    ## tagged GitHub release
    dsph_mw = table.Table.read('https://github.com/apace7/local_volume_database/releases/download/v1.0.0/dwarf_all.csv')
    ## specific commit
    dsph_mw = table.Table.read('https://raw.githubusercontent.com/apace7/local_volume_database/3a473c7f839f228a5702fa0293cebfea5fe3bcb6/data/dwarf_all.csv')

Interactive Website
---------------------------------------------

An interactive version of the LVDB is available here: `link <https://lvd-interactive.streamlit.app/>`_ (made by Katya Gozman).
