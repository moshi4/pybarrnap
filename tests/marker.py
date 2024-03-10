import pytest

from pybarrnap.utils import is_cmscan_installed

skipif_cmscan_not_installed = pytest.mark.skipif(
    not is_cmscan_installed(),
    reason="cmscan(infernal) is not installed.",
)
