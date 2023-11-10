## Package: GreenhouseGases

- module: ``greenhousegases``
- install: ``pip install -e .``
- Python3-ed (``python=3.7``)

## Tests passing

black==black, 21.12b0 (compiled: no)

```
pip install -e .
black --check .
pytest
greenhousegases/GHG_UKCA.py --help
greenhousegases/GHG_radiation.py --help
```
