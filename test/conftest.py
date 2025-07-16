import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--run-optional", action="store_true", default=False, help="run optional tests"
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run-optional"):
        # --run-optional given in cli: do not skip optional tests
        return
    skip_optional = pytest.mark.skip(reason="need --run-optional option to run")
    for item in items:
        if "optional" in item.keywords:
            item.add_marker(skip_optional)
