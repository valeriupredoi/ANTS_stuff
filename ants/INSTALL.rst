===============
Installing ANTS
===============

If you do not have Conda installed, you will need to follow these
`installation instructions <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_
to get started.  Here is the `getting started <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html>`_ guide on using conda.

Stable installation
===================

A stable installation of ANTS uses the same package versions for all
installations of a particular ANTS version.  To create a stable ANTS
environment, including ANTS, using the `environment.lock` file from an ANTS
working copy::

    $ conda create -p <path/to/install/to/ants_x.y> --file environment.lock
    $ conda activate <ants_x.y>
    $ python setup.py install

Some dependencies are not conda/PyPI installable.  In the activated
environment, follow the respective package installation instructions. These
packages are mule, um_spiral_search, and optionally um_utils (See
`Dependencies`_).  Please consult your local UM support team as they may have
suitable versions of these dependencies available that can be used.

A quick verification of the completed installation by activating the environment (if it's not
already active) and running::

    $ python -c 'import ants; print(ants.__file__)'

The returned path should be inside the newly created conda environment.

A more thorough verification can be done by running the test workflow.  To run
the tests in a stable install, run the quick verification command above and cd
to the ants installation directory (i.e. the path returned by the previous
command without the ``__init__.py``).  Then run::

    $ cp -r <working/copy>/rose-stem/sources ./ants/tests/.
    $ cp -r <working/copy>/rose-stem/KGO ./ants/tests/.
    $ python -m pytest --continue-on-collection-errors .

If the two copy commands are omitted, the tests will still run, but around 20%
will error fail due to missing data files.

Developer installation
======================

A developer installation of ANTS only uses the same package versions for key
dependencies, and lets other versions differ.  The key dependency versions are
restricted to ensure that the trunk tests pass, and may not be the same as the
previous stable release.  A developer install needs a conda environment with
all of the ANTS dependencies installed, and a working copy of ANTS for
development work.

For the conda environment, create the ANTS environment without ANTS and using
the `environment.yml` file.  This may use different versions of dependencies
than the previous release::

    $ conda env create -n <ants_x.y> -f environment.yml
    $ conda activate <ants_x.y>

Some dependencies are not conda/PyPI installable.  In the activated
environment, follow the respective package installation instructions. These
packages are mule, um_spiral_search, and optionally um_utils (See
`Dependencies`_).  Please consult your local UM support team as they may have
suitable versions of these dependencies available that can be used.

To checkout a working copy of ANTS use::

   $ fcm co https://code.metoffice.gov.uk/svn/ancil/ants/<branch>
   $ cd <branch>

Typical branches would include `trunk`, `tags/vx.x.x` (tag releases) or
`branches/dev/<userid>/<branch_name>` (development branches).
See the `fcm user guide <http://metomi.github.io/fcm/doc/user_guide/>`_ for more details.

Ensure this checkout is on your PYTHONPATH::

    $ export PYTHONPATH=<PATH_TO_ANTS_CHECKOUT>/lib:${PYTHONPATH}

Verify the completed installation by activating the environment (if it's not
already active) and :ref:`running the tests <running-tests>`.

The returned ANTS location should be inside the working copy, while the
returned iris location should be inside the newly created conda environment.

Most ANTS applications are in contrib.  To check these out, use::

   $ fcm co https://code.metoffice.gov.uk/svn/ancil/contrib/<branch>
   $ cd <branch>

The branch information for contrib includes `trunk`, `tags/vx.x.x` (tag releases)
or `branches/dev/<userid>/<branch_name>` (development branches), the same as
the core ANTS branches.

Dependencies
============

Build and runtime dependencies
------------------------------
Python 3
    ANTS supports Python 3.
    The minimum version of Python required is 3.7.0.

Iris v2.3
    Iris is a Python package for analysing and visualising meteorological and
    oceanographic data sets.  See the `iris documentation
    <https://scitools-iris.readthedocs.io/en/latest/getting_started.html>`_
    for more details.  Note that ANTS is using v2.3 of iris, so there may be
    discrepancies between iris behaviour in ANTS and the documentation for
    more recent versions of iris.

mule 2020.01.1 or later
    A package for working with UM FieldsFiles and variants thereof.
    See `<https://code.metoffice.gov.uk/doc/um/mule/latest>`_.

numba
    JIT compiler that translates a subset of Python and NumPy code into fast
    machine code.

PyKDTree
    Kd-tree implementation for nearest neighbour search in Python.  See
    `<https://github.com/storpipfugl/pykdtree>`_.


Optional dependencies
---------------------

Additional ANTS functionality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
f90nml
    Provides namelist reading support, commonly utilised for defining target
    grid information, see `<https://f90nml.readthedocs.io/en/latest>`_.
    Loading UM grid target grid definitions will require this.

gdal
    Provides raster loading support, commonly utilised format for storing
    large geolocated data such as satellite data.  Uses within ANTS include
    the pre-processing of ITE and IGBP datasets in the derivation of
    vegetation fraction ancillaries.
    See `<https://gdal.org/api/index.html#python-api>`_.

python-stratify
    Vectorized interpolators that are especially useful for Nd vertical
    interpolation/stratification of atmospheric and oceanographic datasets.
    See `<https://github.com/SciTools/python-stratify>`_.
    Ancillaries defined on a set of vertical levels will require stratify,
    such as aerosols and ozone ancillaries.

esmpy
    ESMPy is a Python interface to the Earth System Modeling Framework (ESMF)
    regridding utility.
    See `<https://earthsystemmodeling.org/esmpy>`_.

um_spiral_search
    This "um_spiral_search" module provides a Python extension from the
    SHUMlib spiral search library.
    `<https://code.metoffice.gov.uk/trac/um/browser/mule/trunk/um_spiral_search/README>`_.


ANTS testing
^^^^^^^^^^^^
filelock
    Python package for software testing.  A platform independent file lock.

black
    Python package for software testing.  Python style guide checker.

isort
    Python package for software testing.  Python source code import ordering checker.

flake8
    Python package for software testing.  Python source code style and error checker.

um_utils
    UM file utilities based on the Mule API.  Python package for software
    testing.
    See `<https://code.metoffice.gov.uk/trac/um/browser/mule/trunk/um_utils/README>`_.

nccmp
    Tool for comparing the contents of netCDF files.

Rose
    A toolkit for writing, editing and running application configurations. Used for configuring test workflow applications. See `<http://metomi.github.io/rose/doc/html/index.html>`_.

cylc
    A workflow engine for cycling systems. Used for running the test workflows. See `<https://cylc.github.io/>`_.
