Tutorials \& Examples
=====================

`ipython notebook example <https://github.com/apace7/local_volume_database/blob/main/example_notebooks/example_plots.ipynb>`_ 

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

    dsph_mw = table.Table.read('https://raw.githubusercontent.com/apace7/local_volume_database/main/data/dwarf_mw.csv')
    dsph_m31 = table.Table.read('https://raw.githubusercontent.com/apace7/local_volume_database/main/data/dwarf_m31.csv')
    dsph_lf = table.Table.read('https://raw.githubusercontent.com/apace7/local_volume_database/main/data/dwarf_local_field.csv')
    ufsc = table.Table.read('https://raw.githubusercontent.com/apace7/local_volume_database/main/data/gc_ufsc.csv')
    gc_disk = table.Table.read('https://raw.githubusercontent.com/apace7/local_volume_database/main/data/gc_disk.csv')
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
    plt.plot(ufsc['rhalf_sph_physical'], ufsc['M_V'], 'o',label=r'${\rm UFSC}$')
    plt.plot(gc_disk['rhalf_sph_physical'], gc_disk['M_V'], 'o',label=r'${\rm GC~New}$')
    plt.plot(gc_harris['rhalf_sph_physical'], gc_harris['M_V'], 'o',label=r'${\rm GC~Harris}$')
    plt.gca().set_xscale('log')
    plt.gca().invert_yaxis()
    plt.gca().set_xlabel(r'$r_{1/2}~({\rm pc})$')
    plt.gca().set_ylabel(r'$M_V$')
    plt.legend(loc=(2))
    plt.ylim(3, -20)
    plt.xlim(1, 4e3)
    ## plt.tight_layout()
    ## plt.savefig('/example_rhalf_MV.png')
    plt.show()

Figure from the sample (this does include custom matplotlib stylefile and latex for the captions).

.. figure:: /figures/example_rhalf_MV.png
   :width: 60%
   :align: center
   :alt: Example

   Example figure

Some Recommendations 
---------------------------------------------

For detailed analysis, I would recommendation fixing the version of the tables instead of the current version.  
Each table can be loaded from a specific commit. For example, this loads an older version of the data/dwarf_mw.csv table.

.. code-block:: python

    dsph_mw = table.Table.read('https://raw.githubusercontent.com/apace7/local_volume_database/3a473c7f839f228a5702fa0293cebfea5fe3bcb6/data/dwarf_mw.csv')