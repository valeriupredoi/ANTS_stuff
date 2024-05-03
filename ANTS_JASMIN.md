## Module on JASMIN

- Module name: `ants` (2.0.0)

```
[valeriu@sci2 ~]$ module load ants
[valeriu@sci2 ~]$ module show ants
-------------------------------------------------------------------------------------------------------------------------------------------------
   /apps/jasmin/modulefiles/ants/2.0.0:
-------------------------------------------------------------------------------------------------------------------------------------------------
prepend_path("PATH","/apps/jasmin/community/ANTS/miniconda3/envs/ants2/bin:/apps/jasmin/community/ANTS/miniconda3/envs/ants2/lib/python3.10/site-pack
ages:/apps/jasmin/community/ANTS/miniconda3/bin")
setenv("ESMFMKFILE","/apps/jasmin/community/ANTS/miniconda3/envs/ants2/lib/esmf.mk")
```
- Module location: `/apps/jasmin/modulefiles/ants/2.0.0`
- Conda environment location: `/apps/jasmin/community/ANTS/miniconda3/envs/ants2`
- Conda/mamba specs:
```
[valeriu@sci2 ~]$ mamba --version
mamba 1.5.8
conda 24.4.0
```
- `ants` version:
```
[valeriu@sci2 ~]$ ants-version 
[INFO] ANTS version loaded was:
ants: /apps/jasmin/community/ANTS/ants_2.0.0/lib/ants/__init__.py (version 2.0.0)
[INFO] Iris version loaded was:
iris: /apps/jasmin/community/ANTS/miniconda3/envs/ants2/lib/python3.10/site-packages/iris/__init__.py (version 3.7.1)
```
- Python interpreter:
```
[valeriu@sci2 ~]$ python -V && which python
Python 3.10.13
/apps/jasmin/community/ANTS/miniconda3/envs/ants2/bin/python
```
- check of imports:
```
[valeriu@sci2 ~]$ python -c "import ants"
[valeriu@sci2 ~]$ python -c "from ants import *"
```

## Admin area (``rwx`` rights on source dir)

- loading environment:

```
[valeriu@sci2 ~]$ cd /apps/jasmin/community/ANTS
[valeriu@sci2 ANTS]$ source conda_base.sh 
(base) [valeriu@sci2 ANTS]$ conda activate ants2
```
- important dependencies:

```
ants                      2.0.0                    pypi_0    pypi
iris                      3.7.1              pyha770c72_0    conda-forge
mule                      2023.8.1                  dev_0    <develop>
shumlib                   2023.06.1            h3218e01_0    coecms
esmf                      8.4.2           mpi_mpich_h2a0de38_103    conda-forge
um-spiral-search          2023.8.1                  dev_0    <develop>
esmpy                     8.4.2              pyhc1e730c_4    conda-forge
numpy                     1.26.0          py310hb13e2d6_0    conda-forge
pytest                    8.2.0              pyhd8ed1ab_0    conda-forge
um-utils                  2023.8.1                 pypi_0    pypi
```
- tests:
  - `ants`: 1352 passed, 10 xfailed in 33.15s
  - `mule`: 264 passed, 100 warnings in 0.78s (lots of Numpy 2.0 deprecation warnings)
  - `um_spiral_search`: 6 passed in 0.24s
  - `um_utils`: 62 passed in 0.60s

- notes on dependencies: `esmf==8.4.2` is required for `ants` tests to pass

## Install steps

```
cd /apps/jasmin/community/ANTS
svn checkout https://code.metoffice.gov.uk/svn/ancil/ants/tags/2.0.0/
svn checkout https://code.metoffice.gov.uk/svn/um/mule/trunk/mule
svn checkout https://code.metoffice.gov.uk/svn/um/mule/trunk/um_spiral_search/
Checked out revision 123808
svn checkout https://code.metoffice.gov.uk/svn/um/mule/trunk/um_utils/
Checked out revision 123808
cp -r 2.0.0 /apps/jasmin/community/ANTS/ants_2.0.0
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
/apps/jasmin/community/ANTS/miniconda3
source conda_base.sh
conda -V (conda 24.3.0)
conda update -n base conda
conda install -c conda-forge mamba
(base) [valeriu@sci2 ANTS]$ mamba -V                                                                                                                 
mamba 1.5.8                                                                                                                                          
conda 24.4.0
cd ants_2.0.0
mamba env create -n ants2 -f environment.yml
conda activate ants2
cd ../mule
mamba update -n ants2 -f environment.yml
  - ca-certificates   2024.2.2  hbcca054_0  conda-forge     Cached
  + ca-certificates  2024.3.11  h06a4308_0  pkgs/main       Cached
pip install -e .
Successfully installed mule-2023.8.1
cd lib/mule
pytest
264 passed, 116 warnings in 2.82s
cd ../../../ants_2.0.0
pip install -e .
Successfully built ANTS
Installing collected packages: ANTS
Successfully installed ANTS-2.0.0
cd lib/ants
pytest
162 errorrs: UserWarning:  No module named 'um_spiral_search'
cd ../../../um_spiral_search
add header file https://github.com/metomi/shumlib/blob/master/shum_spiral_search/src/c_shum_spiral_search.h
mamba install -c coecms shumlib
pip install -e .
cd um_utils
python setup.py build
python setup.py install
cd ../../ants_2.0.0
mamba install -c conda-forge esmf=8.4.2
cd lib/ants
pytest
========================================================= 1352 passed, 10 xfailed in 44.47s =========================================================
[WARNING] yaksa: 10 leaked handle pool objects
(that's mpich leaky, no bother)
```
