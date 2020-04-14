import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_almost_equal
from pyconmech.frame_analysis import numpy_stiffness

@pytest.fixture
def tol():
    return 1e-8

@pytest.mark.local_stiff
def test_local_stiffness_matrix(tol):
    pass