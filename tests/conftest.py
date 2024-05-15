"""
Add the '--runslow' command line option to pytest else skip slow markers
"""
import pytest
import base64
import logging
import os

import ee


@pytest.fixture(scope="session", autouse=True)
def test_init():
    logging.basicConfig(level=logging.DEBUG, format='%(message)s')
    logging.getLogger('googleapiclient').setLevel(logging.ERROR)
    logging.debug('Test Setup')

    # For GitHub Actions authenticate using private key environment variable
    if 'EARTHENGINE_TOKEN' in os.environ:
        print('Writing privatekey.json from environmental variable ...')
        content = base64.b64decode(os.environ['EARTHENGINE_TOKEN']).decode('ascii')
        EE_KEY_FILE = 'privatekey.json'
        with open(EE_KEY_FILE, 'w') as f:
            f.write(content)
        ee.Initialize(ee.ServiceAccountCredentials('', key_file=EE_KEY_FILE))
    else:
        ee.Initialize()



def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
