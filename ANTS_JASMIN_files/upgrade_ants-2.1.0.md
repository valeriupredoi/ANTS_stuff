## Date: 7 March 2025

## Install type
Overwrite, ie full upgrade to `ants==2.1.0` with replacement of `ants==2.0.0`

## Dependencies changed

- `gdal==3.9.1`
- `geovista`
- `iris-esmf-regrid>=0.9`
- `iris-sample-data`
- `nccmp`
- `pykdtree`
- `pytest-cov`
- `pytest-xdist`
- `python=3.10.14`
- `ruff`
- `sphinx=7.3.7`

## Install triggers
```
cwd=/apps/jasmin/community/ANTS
cd ants_2.1.0/
pip install -e .
cd lib/ants
pytest
====== 1142 passed, 11 xfailed in 28.17s =========
```

## Mod module
```
cp /apps/jasmin/modulefiles/ants/2.0.0 /apps/jasmin/modulefiles/ants/2.1.0
```

## Final tests
```
[valeriu@sci-vm-01 ~]$ ants-version
[INFO] ANTS version loaded was:
ants: /apps/jasmin/community/ANTS/ants_2.1.0/lib/ants/__init__.py (version 2.1.0)
[INFO] Iris version loaded was:
iris: /apps/jasmin/community/ANTS/miniconda3/envs/ants2/lib/python3.10/site-packages/iris/__init__.py (version 3.7.1)
```
