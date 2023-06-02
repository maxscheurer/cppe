def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "needs_polarizationsolver: marks tests that need the polarizationsolver Python library (deselect with '-m \"not needs_polarizationsolver\"')",
    )
