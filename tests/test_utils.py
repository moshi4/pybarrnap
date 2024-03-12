from pybarrnap.utils import array_split


def test_array_split():
    """Test `array_split()`"""
    assert array_split([1], 1) == [[1]]
    assert array_split([1, 2, 3, 4, 5, 6], 3) == [[1, 2], [3, 4], [5, 6]]
    assert array_split([1, 2, 3, 4, 5], 2) == [[1, 2, 3], [4, 5]]
    assert array_split([1, 2, 3], 5) == [[1], [2], [3]]
