

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
