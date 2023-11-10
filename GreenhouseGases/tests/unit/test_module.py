from pathlib import Path


def test_module():
    """Run a basic module health test."""
    print("This is greenhousegases with Python 3.7")
    test_file = str(Path(__file__).resolve().parent)
    print(f"Unit tests path: {test_file}")
    from greenhousegases import GHG_radiation  # noqa
    from greenhousegases import ghg_rcp_ncplot  # noqa
    from greenhousegases import ghg_rcp_plot  # noqa
    from greenhousegases import ghg_rcp  # noqa
    from greenhousegases import GHG_reform  # noqa
    from greenhousegases import ghg_ukca_ncplot  # noqa
    from greenhousegases import GHG_UKCA  # noqa
