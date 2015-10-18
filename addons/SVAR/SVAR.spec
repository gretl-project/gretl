author = Riccardo "Jack" Lucchetti and Sven Schreiber
email = r.lucchetti@univpm.it
version = 1.2
date = 2015-10-18
description = Structural VARs
public = SVAR_setup SVAR_restrict SVAR_ident SVAR_estimate \
    SVAR_cumulate SVAR_boot SVAR_hd SVAR_coint \
    GetShock IRFplot IRFsave FEVD \
    GUI_SVAR GUI_plot FEVDplot FEVDsave HDplot HDsave \
    SVAR_bundle_print

gui-main = GUI_SVAR

bundle-plot = GUI_plot
bundle-print = SVAR_bundle_print
label = Structural VARs

help = SVAR.pdf
sample-script = examples/simple_C.inp
min-version = 1.10.0
data-requirement = needs-time-series-data
data-files = examples
