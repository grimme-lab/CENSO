# CENSO – Agent Guidelines

CENSO is a Python CLI tool for automated refinement and energetic sorting of
Conformer–Rotamer Ensembles (CREs) at the DFT level.

---

## Repository Layout

```
src/censo/
├── cli/          - CLI parsing and user interface
├── config/       - Pydantic-based configuration models
│   └── parts/    - Per-part configs (prescreening, screening, etc.)
├── ensembleopt/  - Workflow steps (prescreening, screening, optimization, refinement)
├── processing/   - Processors for external programs (ORCA, xtb, Turbomole)
├── properties/   - Property calculators (NMR, optical rotation, UV/Vis)
├── parallel/     - Dask cluster orchestration
├── assets/       - Static data files (DFA, solvents JSON)
└── scripts/      - Auxiliary scripts (c2anmr, nmrplot, uvvisplot)
test/
├── unit/         - Unit tests (mirror of src structure)
└── integration/  - Integration tests (require external programs)
```

---

## Development Setup

```bash
pip install -e .[dev]       # install in editable mode with all dev dependencies
pip install -e .[scripts]   # additionally install numpy/matplotlib/pandas for scripts
```

Requires **Python >= 3.12**.

---

## Build / Lint / Test Commands

### Running Tests

```bash
pytest                          # run all tests (skips optional/integration tests)
pytest -svv                     # verbose output with captured output shown
pytest test/unit/               # unit tests only
pytest test/unit/test_ensemble.py                        # single test file
pytest test/unit/test_ensemble.py::TestEnsembleDataReadOutput  # single class
pytest test/unit/test_ensemble.py::TestEnsembleDataReadOutput::test_read_output_success  # single test
pytest -k "test_read_output"    # pattern-match test name
pytest --run-optional           # include optional/integration tests
pytest --keep-log               # keep censo.log after test run
tox                             # test against Python 3.12 and 3.13 (CI matrix)
tox -e py312                    # test against specific Python version
```

Custom pytest markers (defined in `pyproject.toml`):
- `optional` – skipped by default; run with `--run-optional`
- `requires_xtb` / `requires_orca` / `requires_turbomole` / `requires_cosmotherm` –
  skipped when the corresponding external program is unavailable

### Formatting

```bash
black src/ test/                # auto-format (line length: black default 88)
```

### Linting

```bash
ruff check src/ test/           # lint and check imports
ruff check --fix src/ test/     # auto-fix safe issues
```

### Type Checking

```bash
mypy src/                       # type-check the source tree
# mypy is run with --explicit-package-bases --check-untyped-defs
```

### Pre-commit (run before every commit)

```bash
pre-commit run --all-files      # runs trailing-whitespace, EOF, black, ruff, pyupgrade, mypy
pre-commit run --files <file>   # run only on specific files
```

The full pre-commit pipeline (`.pre-commit-config.yaml`):
1. `trailing-whitespace`, `end-of-file-fixer`, `check-toml`, `check-yaml`,
   `check-added-large-files` (max 10 MB), `debug-statements`
2. `ruff` (with `--fix`)
3. `black`
4. `pyupgrade`
5. `mypy` (`--explicit-package-bases --check-untyped-defs`)

CI runs pre-commit and tox (py312 + py313) on every push/PR.

---

## Code Style Guidelines

### General

- Formatter: **black** (default settings, line length 88).
- Linter/import sorter: **ruff**.
- Type checker: **mypy** with `--check-untyped-defs`.
- Target language version: **Python 3.12+**; use modern syntax freely (e.g.,
  `list[str]`, `X | Y`, `type Alias = ...`, generic classes `class Foo[T]:`).

### Naming Conventions

| Symbol | Convention |
|---|---|
| Functions, methods, variables | `lower_snake_case` |
| Classes | `UpperCamelCase` |
| Constants / module-level enums | `UPPER_SNAKE_CASE` |
| Private instance attributes | `__double_underscore` (name-mangled) |
| Protected / internal helpers | `_single_underscore` |
| Logger instance | `logger = setup_logger(__name__)` at module level |

### Imports

- Order enforced by ruff: stdlib → third-party → local (relative).
- Use **relative imports** within the `censo` package (e.g., `from ..molecules import ...`).
- Prefer `from pathlib import Path` over `os.path` for filesystem operations.
- Use `from collections.abc import Callable` (not `typing.Callable`).

### Type Annotations

- Annotate all public function signatures and class attributes.
- Use built-in generics directly (`list[str]`, `dict[str, int]`, `tuple[int, ...]`).
- Union types: `X | Y` syntax (not `Union[X, Y]`).
- `Any` is allowed where truly dynamic; prefer narrower types when possible.
- Pydantic models use `model_config = ConfigDict(...)` class variable (v2 API).

### Docstrings

- Required on all public classes and public functions/methods.
- Module-level docstring at the top of each file describing its purpose.
- Parameter/return documentation uses Sphinx-style (`:param name:`, `:type name:`,
  `:return:`, `:rtype:`).

```python
def prescreening(ensemble: EnsembleData, config: PartsConfig) -> dict[str, Any]:
    """
    Short description of what this does.

    :param ensemble: The conformer ensemble to process.
    :type ensemble: EnsembleData
    :param config: Configuration for this part.
    :type config: PartsConfig
    :return: JSON-serializable results dictionary.
    :rtype: dict
    """
```

### Configuration / Pydantic Models

- All configuration classes extend `GenericConfig` (which extends `pydantic.BaseModel`).
- Use `ConfigDict(str_to_lower=True, str_strip_whitespace=True, validate_default=True,
  use_attribute_docstrings=True)`.
- Enum fields should extend both `str` and `Enum` for JSON-safe serialization.
- Centralize defaults in `src/censo/config/parts/`.

### Error Handling

- Raise concrete exceptions (`ValueError`, `RuntimeError`, `TypeError`) with descriptive
  messages rather than bare `Exception`.
- Use Pydantic `ValidationError` for configuration input errors; surface them via
  `print_validation_errors()` from `censo.utilities`.
- Log errors/warnings with the module-level `logger` (never `print()` for diagnostics).
- User-facing output uses `printf()` from `censo.utilities` (wraps `print` with
  consistent formatting).

### Logging

- Obtain a logger per module: `logger = setup_logger(__name__)` from `censo.logging`.
- Do not call `logging.basicConfig()` or configure handlers outside `censo.logging`.

### Parallel / Dask

- Parallel work is dispatched via `censo.parallel.execute()`; do not call
  `client.submit()` / `client.map()` directly in workflow code.
- Keep objects that cross the Dask boundary small and picklable (e.g., pass
  `GeometryData` rather than full `MoleculeData`).

### Tests

- Tests live in `test/unit/` mirroring the `src/censo/` structure.
- Test classes use `UpperCamelCase` prefixed with `Test`; test functions are
  `lower_snake_case` prefixed with `test_`.
- Use `pytest.approx` for floating-point comparisons.
- Fixtures (shared data) go in `test/unit/fixtures/`; the `fixtures_path` fixture
  is defined in `test/unit/conftest.py`.
- The `tmp_wd` autouse fixture changes the working directory to a temp path for every
  test—do not rely on the project root as CWD inside tests.
- Mark slow or external-tool-dependent tests with `@pytest.mark.optional` or the
  appropriate `requires_*` marker.
- Do not modify test files unless explicitly approved.
