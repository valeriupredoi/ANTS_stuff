import pytest

from pathlib import Path


def test_module():
    """Run a basic module health test."""
    print("This is GreenhouseGases with Python 3.7")
    test_file = str(Path(__file__).resolve().parent)
    print(f"Unit tests path: {test_file}")
    from GreenhouseGases import GHG_radiation
    from GreenhouseGases import ghg_rcp_ncplot
    from GreenhouseGases import ghg_rcp_plot
    from GreenhouseGases import ghg_rcp
    from GreenhouseGases import GHG_reform
    from GreenhouseGases import ghg_ukca_ncplot
    from GreenhouseGases import GHG_UKCA
