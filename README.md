### ANTS Stuff

Stuff related to the MO package ANTS

### Resources

- [ANTS User Guide](https://code.metoffice.gov.uk/doc/ancil/ants/latest/index.html)
- [ANTS on MOSRS](https://code.metoffice.gov.uk/trac/ancil)
- [ANTS codebase at CodeMO](https://code.metoffice.gov.uk/svn/ancil/ants/)
- [ANTS contrib](https://code.metoffice.gov.uk/trac/ancil/browser/contrib/)

### 16 October 2023

Trying to install and get ANTS to work with latest Python and ``iris``.

Work: local, own machine, using

```bash
mamba 1.4.9
conda 23.5.2
```

#### Installation Steps

- use ``environment.yml`` but with no hard pins
- run the env creation ``mamba env create -n new_ants -f environment.yml``
- env solves well and now we have:

```
python 3.11.6 hab00c5b_0_cpython conda-forge
iris   3.7.0  pyha770c72_0       conda-forge
```
- use ``pip`` to install:

```
pip install -e .
...
Successfully built ANTS
Installing collected packages: ANTS
Successfully installed ANTS-1.2.0.dev0
```
- run tests with ``pytest``: 185 errors, all ``E   ModuleNotFoundError: No module named 'mule'``;
- install ``mule`` from Conda with ``mamba install -c coecms mule``
- run ``pytest`` again - one single error:
```
lib/ants/fileformats/_ugrid.py:335: in _construct_mesh
    node_latitudes = self.cubes.extract(
E   TypeError: CubeList.extract() got an unexpected keyword argument 'strict'
```
- remove all instances of `strict` from `iris.Constraint`
- adjust all instances of resulted `CubeList` instead of `Cube` (to point to a single cube, that is the only element of the cubes list)
- handle `0000-00-00 00:00:00` times that cftime doesn't like: in `lib/ants/fileformats/pp/__init__.py` modufy the time converter:

```python
    def int_time(time):
        # Converts a time of "1-02-03 04:05:06" to 10203040506
        if time.month == 0 and time.day == 0:
            time = cftime.datetime(time.year, 1, 1,
                                   time.hour, time.minute,
                                   time.second)
        return int(time.strftime("%Y%m%d%H%M%S").strip())
```
- `cubes.extract_strict()` -> `cubes.extract_cube()`
- iris exceptions of ConstraintMismatchError can be circumvented, and just check on the result,
  if null, apply the `except` clause from previous try/except implement;
- please don't assume `nccmp` is a given in the system!
- `numpy.float` etc deprecated (I haven't yet fixed that)
- in total:
```
M       environment.yml
M       lib/ants/coord_systems.py
M       lib/ants/fileformats/_ugrid.py
M       lib/ants/fileformats/pp/__init__.py
M       lib/ants/tests/fileformats/ancil/test_integration.py
M       lib/ants/tests/fileformats/ugrid/test_integration.py
```
- after this first pass API changes: 122 failed, 1034 passed, 165 skipped, 10 xfailed, 576 warnings in 55.57s 
