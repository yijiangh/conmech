import os
import pytest

def pytest_addoption(parser):
    parser.addoption('--viewer', action='store_true', help='enable viewer for visualization.')
    parser.addoption(
        "--engine",
        choices={'cpp', 'numpy'},
        action="append",
        default=[],
        help="list of checker engines",
    )
    parser.addoption('--db', action='store_true', help='debug mode for tests')

@pytest.fixture
def test_data_dir():
    root_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(root_dir, 'test_data')
    return json_path

@pytest.fixture
def n_attempts():
    return 10

@pytest.fixture
def viewer(request):
    return request.config.getoption("--viewer")

def pytest_generate_tests(metafunc):
    if "engine" in metafunc.fixturenames:
        engine_opts = metafunc.config.getoption("engine")
        if len(engine_opts) > 0:
            metafunc.parametrize("engine", engine_opts)
        else:
            metafunc.parametrize("engine", ['cpp', 'numpy'])

@pytest.fixture
def debug(request):
    return request.config.getoption("--debug")