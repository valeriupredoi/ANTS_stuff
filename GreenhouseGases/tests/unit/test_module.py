import pytest

from pathlib import Path


def test_module():
    """Run a basic module health test."""
    print("This is greenhousegases with Python 3.7")
    test_file = str(Path(__file__).resolve().parent)
    print(f"Unit tests path: {test_file}")
    from greenhousegases import GHG_radiation
    from greenhousegases import ghg_rcp_ncplot
    from greenhousegases import ghg_rcp_plot
    from greenhousegases import ghg_rcp
    from greenhousegases import GHG_reform
    from greenhousegases import ghg_ukca_ncplot
    from greenhousegases import GHG_UKCA
