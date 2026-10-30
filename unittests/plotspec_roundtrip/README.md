# Round-trip test for the plotspec reader and writer

Checks `plotspec_read_file()` (`lib/src/plotspec_read.c`) against
`plotspec_print()` (`lib/src/plotspec.c`), the two halves of what the GUI
plot editor needs.

    make check

This generates a set of plots with `gretlcli` (`gen_plots.inp`), then
`roundtrip.c` takes each through read -> write -> read -> write and checks:

1. the two intermediate specs have the same structure and the two written
   scripts are identical (reader and writer agree);
2. if gnuplot is installed, the original script and the first written script
   draw the same thing on the text-only "dumb" terminal (nothing was lost on
   the first pass, which check 1 cannot see).

Plot types the GUI editor does not support (`cant_edit()`) are skipped.
Known failures can be listed with `-x name` (see `roundtrip.c`); there are none at present.

    make check-asan

additionally compiles the reader (and `plotspec.c`) into the test binary with
AddressSanitizer, UBSan and leak detection, repeats the round trip, and then
feeds the reader tens of thousands of damaged variants of the plots
(truncated, lines deleted, blanked or cut short, random byte changes). The
reader may reject these, but must not misuse memory or leak.

Use `./roundtrip -v -k plots/foo.plt` to see successes and keep the
intermediate files for inspection. Needs `build/` in the top-level source
directory, with libgretl and gretlcli built.
