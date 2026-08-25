# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Gretl (GNU Regression, Econometrics and Time-series Library) is a C/C++ econometrics
package: a core library (`libgretl`), a CLI client, a GTK GUI, dynamically loaded
plugins, and a collection of function packages ("addons") written in gretl's own
scripting language, hansl.

An `AGENTS.md` at the repo root already covers build/test commands, key directories,
and conventions in detail — read it first. This file adds architecture that spans
multiple files/directories and isn't obvious from a directory listing alone.

## Architecture

### Layering

```
addons/   hansl (.inp) function packages, compiled to .gfn via gretlcli -t pkg.inp
gui/      GTK client (gretl.c entry point) -> calls into lib/
cli/      command-line client (gretlcli) -> calls into lib/
plugin/   .so/.dll modules, dlopen'd at runtime by lib/ for optional/heavy features
lib/src/  libgretl: the engine everything else is built on
```

`lib/src/` is the core (~110 C files, headers in the same dir, all pulled in via
`lib/src/libgretl.h`). Both `cli` and `gui` are thin front ends over it, so a fix or
feature typically lives in `lib/src/` first, with `cli/`/`gui/` only adding
presentation. `plugin/` holds functionality that's optional, GPL-incompatible-license-
adjacent, or heavyweight (import filters for Stata/SPSS/SAS/Excel/ODS, ARMA/GARCH/VAR
estimators, kalman filter, SVM, lpsolve, mailer, zip handling, etc.) — these are
discovered and loaded at runtime rather than linked directly, see `gretl_plugin.c`/
`plugins.h`.

### Core data model (`lib/src/libgretl.h`)

- **`DATASET`**: the dataset in memory — `Z` is the `double**` data matrix (series
  indexed by variable number), `varname`/`varinfo` hold names/metadata, `t1`/`t2`
  bound the current (sub)sample, `submask`/`restriction` record subsampling. Nearly
  every estimation and data-handling function takes a `DATASET *` argument.
- **`MODEL`** (`gretl_model.h`): result of an estimation command — coefficients,
  standard errors, test results (`ModelTest`), diagnostics.
- **`gretl_matrix`** (`gretl_matrix.c`/`.h`): the numeric matrix type used throughout,
  backed by LAPACK/BLAS (openblas or lapack+blas, auto-detected by `configure`).
- **`gretl_bundle`** / **`gretl_array`**: hansl's associative container and array
  types, used both internally and as the standard way user functions and plugins
  return structured results.

### The hansl language (gretl's scripting language)

Scripts (`.inp` files) and function packages are written in hansl. Its implementation
is split across:

- `interact.c` / `gretl_commands.c`: parses and dispatches gretl *commands* (`ols`,
  `smpl`, `open`, ...) — the imperative statement layer.
- `genmain.c` / `genlex.c` / `gensyntax.c` / `geneval.c` (the largest file in the
  tree): the `genr`/expression engine — lexer, parser, and evaluator for hansl
  expressions (matrix algebra, function calls, series arithmetic).
  `genfuncs.c` implements the built-in function library called from there.
  `gretl_bfgs.c`, `nls.c`-style numerical routines etc. are invoked from user scripts
  via this layer.
- `gretl_func.c`: user-defined *functions* — packaging hansl code into callable units
  with typed parameters/returns, as used by function packages (addons).
- `monte_carlo.c`: `loop` blocks (the `for`/`while`/index-loop construct).

`cli/` and `gui/console.c` both feed user input through this same command/expression
pipeline, so behavior should be consistent between the two front ends by construction
— don't special-case one without checking the other.

### Tests

Two independent test trees, both hansl-based (see `AGENTS.md` for directory layout
and the Given/When/Then test-writing convention used in `unittests/`):
- `tests/`: NIST regression-validation datasets plus longer-running "practice"
  scripts (integration smoke tests with no assertions).
- `unittests/`: `assert()`-based unit tests, organized by fundamentals/commands/
  functions/plots, run via `unittests/run_tests.sh`.

Addon packages also carry their own `tests/` subdirectories, run via their Makefile's
`check` target.
