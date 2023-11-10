## Package: GreenhouseGases

- module: ``greenhousegases``
- install: ``pip install -e .``
- Python3-ed (``python=3.7``)

## Tests passing

```
pip install -e .
flake8
pytest
greenhousegases/GHG_UKCA.py --help
greenhousegases/GHG_radiation.py --help
```
