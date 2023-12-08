Usage
=====

.. _installation:

Installation (note that this doesn't work yet)
------------

To use local_volume_database, first install it using pip:

.. code-block:: console

   (.venv) $ pip install local_volume_database

Database content
----------------

The tables can be directly loaded into jupyter notebooks without having to download the repository:

>>> import astropy.table as table
>>> dwarf_mw = table.Table.read('https://raw.githubusercontent.com/apace7/local_volume_database/main/data/dwarf_mw.csv')



