from pathlib import Path


def test_module():
    """Run a basic module health test."""
    print("This is nitrogendeposition with Python 3.7")
    test_file = str(Path(__file__).resolve().parent)
    print(f"Unit tests path: {test_file}")
    from nitrogendeposition import cmip6_ndep_jasmin
    from nitrogendeposition import cmip6_ndep_jasmin_vn16
    from nitrogendeposition import cmip6_ndep_jasmin_vn17
    from nitrogendeposition import cmip6_ndep_jasmin_vn21
    from nitrogendeposition import cmip6_ndep_jasmin_vn23_ssp585
    from nitrogendeposition import cmip6_ndep_jasmin_vn25_past_and_future
    from nitrogendeposition import cmip6_ndep_jasmin_vn26_past_and_future
    from nitrogendeposition import cmip6_ndep
