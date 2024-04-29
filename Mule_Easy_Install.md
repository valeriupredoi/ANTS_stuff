## Install Mule Supereasy

- install mamba: in a conda environment: `conda install -c conda-forge mamba`
- grab latest `trunk` for `mule`: `svn checkout https://code.metoffice.gov.uk/svn/um/mule/trunk/mule`
- go inside the mule: `cd mule`
- create conda environment file:

```yaml
---
name: mule
channels:
  - conda-forge
  - coecms
  - main

dependencies:
  - blas
  - intel-openmp
  - libgfortran-ng
  - libgfortran5
  - libmo_unpack
  - mkl
  - mkl-service
  - mkl_fft
  - mkl_random
  - mo_pack
  - numpy
  - numpy-base
  - pytest
  - python >=3.9
  - shumlib
  - six
  - tbb
```
- create `mule` conda environment: `mamba env create -n mule -f environment.yml`
- activate upon creation: `conda activate mule`
- pip-install `mule`: `pip install -e .`
- test installation:

```bash
cd lib/mule
pytest
```

- done!

