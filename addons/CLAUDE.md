# CLAUDE.md — addons/

Each subdirectory (e.g. `regls`, `extra`, `SVAR`, `dbnomics`, `gig`, `kmeans`, ...) is
a self-contained hansl function package: `<name>.inp` (source), `<name>.spec`
(package metadata), `pkg.inp.in` (build driver), and a `tests/` subdir. `make` in an
addon directory invokes `gretlcli -t pkg.inp` to compile the `.inp` sources into a
`.gfn` package file, which is what gretl actually loads. When editing an addon, treat
the `.inp` files as the source of truth; the `.gfn` is a build artifact.
