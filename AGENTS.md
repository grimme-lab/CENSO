# CENSO – Agent Guidelines

CENSO is a Python CLI tool for automated refinement and energetic sorting of
Conformer–Rotamer Ensembles (CREs) at the DFT level.

---

## Build / Lint / Test Commands

### Running Tests

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

CI runs pre-commit and tox (py312 + py313) on every push/PR.

---

## Code Style Guidelines

### Naming Conventions

| Symbol | Convention |
|---|---|
| Functions, methods, variables | `lower_snake_case` |
| Classes | `UpperCamelCase` |
| Constants / module-level enums | `UPPER_SNAKE_CASE` |
| Private instance attributes | `__double_underscore` (name-mangled) |
| Protected / internal helpers | `_single_underscore` |
| Logger instance | `logger = setup_logger(__name__)` at module level |

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

### Error Handling

- Raise concrete exceptions (`ValueError`, `RuntimeError`, `TypeError`) with descriptive
  messages rather than bare `Exception`.
- Use Pydantic `ValidationError` for configuration input errors; surface them via
  `print_validation_errors()` from `censo.utilities`.
- Log errors/warnings with the module-level `logger` (never `print()` for diagnostics).
- User-facing output uses `printf()` from `censo.utilities` (wraps `print` with
  consistent formatting).
