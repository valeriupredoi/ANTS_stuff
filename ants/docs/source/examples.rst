
Example usages
==============

Regrid
------

Without filling missing points (notice ``--target-grid`` argument):

.. code-block::

    ancil_general_regrid.py --ants-config path/to/config/file --target-grid <grid filenames> --output <output filename> <source filenames>
..

With filling missing points for consistency with the land sea mask (notice ``--target-lsm`` argument):


.. code-block::

    ancil_general_regrid.py --ants-config path/to/config/file --target-lsm <lsm filename> --output <output filename> <source filenames>
..

Please note that neither horizontal nor vertical regridding schemes are set by default -
these must be configured via the configuration file.
For examples, see the `Regridding schemes`_ section below.

More details can be found with:


.. code-block::

    ancil_general_regrid.py --help
..

ancil_2anc
----------

.. code-block::

  ancil_2anc.py source_file -o output_file --grid-staggering=6
..

This will result in an ancillary file at `output_file`, and a NetCDF file at
`output_file.nc`.  The ancillary file will have the fixed length header
grid staggering set as a value of 6.  All other pieces of metadata will be
determined from the content of the `source_file`.

Configuration
=============

More details on the available configuration options can be found :class:`here <ants.config.GlobalConfiguration>`.

.. code-block::

    --ants-config <configuration filename>
..

to the command.

Regridding schemes
------------------

 For horizontal regridding, the valid scheme names are:

    **Linear, TwoStage, ConservativeESMF, AreaWeighted, Nearest**

 For vertical regridding, the valid scheme names are:

    **Linear, Conservative, Nearest**

 For vertical regridding it is also possible to specify the extrapolation mode.  By default, there is no extrapolation.  To enable extrapolation, use one of the valid mode names:

    **Nearest, Linear**

 As an example, to specify regridding using ESMF's first order conservative regridding for the horizontal, linear for the vertical regridding, and nearest neighbour vertical extrapolation, use:


.. code-block::

    [ants_regridding_horizontal]
    scheme=ConservativeESMF

    [ants_regridding_vertical]
    scheme=Linear
    extrapolation_mode=Nearest
..

Decomposition
-------------

By default (as of ANTS 1.0), decomposition is disabled. To enable it, the size of the decomposed chunks must be specified. For example, for an 8x6 split, use:

.. code-block::

    [ants_decomposition]
    x_split=8
    y_split=6
..

It is possible, though discouraged, to have ANTS decompose to 800Mb chunks. To do this, set both splits to 'automatic'. This option may be removed in a future version of ANTS.

Decomposition can also be disabled by setting both splits to 0:

.. code-block::

    [ants_decomposition]
    x_split=0
    y_split=0
..

