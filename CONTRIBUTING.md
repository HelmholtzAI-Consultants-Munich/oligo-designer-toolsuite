# Contribution Guidelines for Oligo Designer Toolsuite

Thank you for your interest in contributing to Oligo Designer Toolsuite! This guide explains how to set up your environment, follow our style, add features, write tests, and submit pull requests.

## General Guidelines

1. **Fork and clone**: Fork the repository and clone your fork locally. Always work in a new feature branch.
2. **Keep in sync**: Regularly sync your branch with `upstream` (the main repository) to avoid large conflicts.
3. **Commit messages**: Write clear, imperative commit messages (e.g., "Add X", "Fix Y").
4. **Pull requests**: Open a PR from your branch to this repository. Provide a concise description, screenshots/logs when relevant, and link any related issues.

## Environment Setup

1. **Supported Python and platforms**: Tested with Python 3.10–3.12 on Linux (x64) and macOS (x64, arm64).
2. **Recommended**: Use a conda environment for a stable installation.
   - If your institution does not support Anaconda, use miniforge: `https://github.com/conda-forge/miniforge`.
3. **Create and activate an environment** (choose one):
   - Using conda:
     ```bash
     conda create -n odt python=3.12 -y
     conda activate odt
     ```
   - Using `venv`:
     ```bash
     python3 -m venv .venv
     source .venv/bin/activate
     ```
4. **Install the project (editable) and dev tools**:
   ```bash
   pip install -e .[dev]
   ```
5. **Install pre-commit hooks**:
   ```bash
   pip install pre-commit
   pre-commit install
   ```

## Code Style and Tooling

We use the following tools (configured in `pyproject.toml` and `.pre-commit-config.yaml`):

- **Black**: code formatting (line length 110)
- **isort**: import sorting (Black profile)
- **autoflake**: remove unused imports/variables
- **pre-commit**: runs the above hooks automatically

Before committing, run:
```bash
pre-commit run --all-files
```

## Implementing Modules

1. **Inherit from the appropriate base class**: When adding a new class, inherit from the relevant base class for the submodule. This enforces a consistent structure. For example, a property filter should inherit from `BasePropertyFilter` and must implement the `apply` method.
2. **Design for clarity**: Prefer descriptive names and small, focused functions/classes.
3. **Unit tests**: Add tests in `tests/test_<module_name>.py` that cover success and error cases.

## Implementing Pipelines

1. **Reuse existing modules**: Compose pipelines using existing modules where possible. If new steps are required, add them to the library so they can be reused.
2. **Configuration**: Ensure new pipelines accept configuration via YAML when appropriate and document the fields.
3. **CLI entry points**: If adding a new CLI, expose it via `project.scripts` in `pyproject.toml` and include usage examples.

## Testing

- Run the test suite:
  ```bash
  pytest
  ```
- Add focused unit tests for new behavior and update fixtures if needed.
- Ensure tests pass locally before opening a PR.

## Documentation

Please document public APIs and user-facing behavior using reStructuredText (reST). This standard integrates with Sphinx for consistent docs.

Example reST-style docstring:
```
def function(arg1, arg2):
    """
    A brief description of the function.

    :param arg1: Description of arg1.
    :type arg1: int
    :param arg2: Description of arg2.
    :type arg2: str
    :returns: Description of return value.
    :rtype: bool
    """
    pass
```

To build the documentation locally:
```bash
pip install -r docs/requirements.txt
cd docs
make html
```
The HTML will be generated in `docs/_build/html`.

## Pull Request Checklist

- [ ] Code formatted with Black and imports sorted with isort
- [ ] `pre-commit run --all-files` passes with no changes
- [ ] Tests added/updated and `pytest` passes locally
- [ ] Documentation added/updated (docstrings and, if relevant, Sphinx pages)
- [ ] Clear PR title and description; linked to related issues

## Thank you!
By following these guidelines, you help maintain the quality and reliability of Oligo Designer Toolsuite. We appreciate your contributions!
