import os
import pytest

def pytest_addoption(parser):
    parser.addoption('--viewer', action='store_true', help='enable viewer for visualization.')

@pytest.fixture
def test_data_dir():
    root_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(root_dir, '..', 'test_data')
    return json_path

@pytest.fixture
def n_attempts():
    return 10

@pytest.fixture
def viewer(request):
    return request.config.getoption("--viewer")