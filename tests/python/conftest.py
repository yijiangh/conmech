import os
import pytest

@pytest.fixture
def test_data_dir():
    root_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(root_dir, '..', 'test_data')
    return json_path

@pytest.fixture
def n_attempts():
    return 10