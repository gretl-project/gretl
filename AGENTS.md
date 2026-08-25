# AGENTS.md

Gretl (GNU Regression, Econometrics and Time-series Library) - C/C++ econometrics software with CLI, GUI, and plugin architecture.

## Build

```bash
./configure
make
make check
sudo make install
```

Quick CI build: `./configure -q --enable-quiet-build --disable-nls`

Build specific targets: `make clean lib cli gui plugin po share tests doc xdg addons`

### Build dependencies

Parallel make respects these dependencies (see `Makefile.in:129-136`):
- `lib` must build first (depends on `buildstamp`)
- `cli`, `share`, `gui`, `editor`, `plugin`, `tests` depend on `lib`
- `addons` depends on `cli`

### LAPACK configuration

Auto-detects `libopenblas` or `liblapack+libblas`. For non-standard locations:
```bash
LAPACK_LIBS='-L/opt/openblas/lib -lopenblas' ./configure
```

## Tests

```bash
make check                           # NIST tests + unit tests

# Individual test suites (from unittests/)
./run_tests.sh --practice            # Integration tests (no assertions)
./run_tests.sh --fundamentals        # Math fundamentals
./run_tests.sh --commands            # Gretl commands
./run_tests.sh --functions           # Function unit tests
./run_tests.sh --plots               # Plot tests
./run_tests.sh --all                 # All above

# Single test script
gretlcli -b -e -q path/to/script.inp
```

Test files use `.inp` extension (hansl scripting language). Tests use `assert()` or the assertion package.

## Key directories

| Directory | Purpose |
|-----------|---------|
| `lib/src/` | Core libgretl library (~220 C files) |
| `cli/` | Command-line interface (gretlcli) |
| `gui/` | GTK-based GUI |
| `plugin/` | Dynamically loaded plugins |
| `addons/` | User function packages (hansl) |
| `tests/` | NIST regression validation |
| `unittests/` | Hansl script-based tests |
| `cephes/`, `minpack/`, `rng/` | Math libraries |

## Conventions

- Version defined in `lib/src/version.h`
- Plugins are `.so`/`.dll` loaded at runtime
- Test scripts named `run_*.inp`
- `gretlcli` flags: `-b` batch, `-q` quiet, `-e` exit on error
