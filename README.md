# **local_volume_database** 

### DESCRIPTION:

Database of Local Volume dwarf galaxies and  star clusters. The database is complete for dwarf galaxies within ~3 Mpc. 
The planned limiting distance is ~10+ Mpc (i.e. resolved stars with HST or JWST).  
The star cluster collections are currently limited to old star clusters/globular clusters in the Milky Way halo and globular clusters/star clusters hosted by low mas dwarf galaxies.

The main data tables are located in the release pages in `csv` and `fits` file formats. The data tables are also available in the [`data/`](data/) folder in `csv` format. 
A summary pdf document of the database content is available in the release page. 

Documentation: [readthedocs](https://local-volume-database.readthedocs.io/en/latest/index.html).

Overview paper: [arXiv](https://arxiv.org/abs/2411.07424) and  [ADS](https://ui.adsabs.harvard.edu/abs/2024arXiv241107424P/abstract).

The tables can be directly loaded into Jupyter notebooks without having to download the repository:

    import astropy.table as table
    ## from the release page (recommended)
    dwarf_mw = table.Table.read('https://github.com/apace7/local_volume_database/releases/download/v1.0.0/dwarf_all.csv')
    # or from a  github branch
    dwarf_mw = table.Table.read('https://raw.githubusercontent.com/apace7/local_volume_database/main/data/dwarf_mw.csv') 


### ACKNOWLEDGEMENT:

If you use this in your research please include a link to the github repository (https://github.com/apace7/local_volume_database) and cite the database paper. It is highly encouraged to cite the input references of the database if they are used in your work. An example in latex is: This work has made use of the Local Volume Database\footnote{\url{https://github.com/apace7/local_volume_database }} \citep{Pace2024arXiv241107424P}. The BibTex of the citation is available here [ADS]{https://ui.adsabs.harvard.edu/abs/2024arXiv241107424P/exportcitation}.

The LVDB releases are also indexed on [zenodo](https://doi.org/10.5281/zenodo.14076714).

### INSTALLATION and USE:

Users that only want to use the catalogs do not need to install the package. 
The catalogs are available as csv and fits files in the release pages.
The package is focused on helper functions and yaml file access. 


The package is installable locally via pip:

```
git clone https://github.com/apace7/local_volume_database.git
cd local_volume_database
python -m build
pip install .
```

The `LVDBDIR` environment variable is used to point to the location of the input YAML files (/data_input/). 
