## Package: NitrogenDeposition

- module: ``nitrogendeposition``
- install: ``pip install -e .``
- Python3-ed (``python=3.7``)

## Tests passing

black==black, 21.12b0 (compiled: no)

```
pip install -e .
black --check .
pytest
nitrogendeposition/cmip6_ndep_jasmin_vn21.py --help
nitrogendeposition/cmip6_ndep_jasmin_vn23_ssp585.py --help
nitrogendeposition/cmip6_ndep_jasmin_vn25_past_and_future.py --help
nitrogendeposition/cmip6_ndep_jasmin_vn26_past_and_future.py  --help
```
