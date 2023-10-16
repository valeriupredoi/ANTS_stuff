### ANTS Stuff

Stuff related to the MO package ANTS

### Resources

- [ANTS User Guide](https://code.metoffice.gov.uk/doc/ancil/ants/latest/index.html)
- [ANTS on MOSRS](https://code.metoffice.gov.uk/trac/ancil)
- [ANTS codebase at CodeMO](https://code.metoffice.gov.uk/svn/ancil/ants/)

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

