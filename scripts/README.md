### ANTS

Environment created from stock yml file.

ANTS installed from ``trunk/`` at 25 October 2023 ``ants=ANTS-1.2.0.dev0``.

```
pip install -e .
```

Mule installed from conda `mamba install -c coecms mule`.

Tests run: 

```
(ants) valeriu@valeriu-PORTEGE-Z30-C:~/ANTS_trunk$ ants-version 
/home/valeriu/ANTS_trunk/lib/ants/analysis/_merge.py:35: UserWarning:  No module named 'um_spiral_search'
Unable to import "spiral", proceeding without the capabilities it provides.  See install.rst
  warnings.warn(msg.format(str(_SPIRAL_IMPORT_ERROR)))
[INFO] ANTS version loaded was:
ants: /home/valeriu/ANTS_trunk/lib/ants/__init__.py (version 1.2.0dev)
[INFO] Iris version loaded was:
iris: /home/valeriu/miniconda3/envs/ants/lib/python3.7/site-packages/iris/__init__.py (version 2.3.0)
```
and
```
pytest >
7 failed, 1291 passed, 22 skipped, 11 xfailed, 2 warnings in 73.81s (0:01:13)
```

(all fails due to missing ``nccmp``)

### 2to3 NitrogenDeposition

#### Conversion

```
2to3-3.7 NitrogenDeposition_py2/* -o NitrogenDeposition -n -w
```

### Package NitrogenDeposition

``NitrogenDeposition`` is now a barebones package installable with ``pip install -e .``
