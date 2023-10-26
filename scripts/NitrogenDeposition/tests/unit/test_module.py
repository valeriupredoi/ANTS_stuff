import pytest

from pathlib import Path


def test_module():
    """Run a basic module health test."""
    print("This is NitrogenDeposition with Python 3.7")
    test_file = str(Path(__file__).resolve().parent)
    print(f"Unit tests path: {test_file}")
    from NitrogenDeposition import cmip6_ndep_jasmin
    from NitrogenDeposition import cmip6_ndep_jasmin_vn16
    from NitrogenDeposition import cmip6_ndep_jasmin_vn17
    from NitrogenDeposition import cmip6_ndep_jasmin_vn21
    from NitrogenDeposition import cmip6_ndep_jasmin_vn23_ssp585
    from NitrogenDeposition import cmip6_ndep_jasmin_vn25_past_and_future
    from NitrogenDeposition import cmip6_ndep_jasmin_vn26_past_and_future
    from NitrogenDeposition import cmip6_ndep
