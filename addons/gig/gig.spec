author = Riccardo "Jack" Lucchetti and Stefano Balietti
email = r.lucchetti@univpm.it
tags = C22
version = @VERSION@
date = 2025-03-19
description = An assortment of univariate GARCH models
public = GUI_gig \
    gig_setup gig_set_dist gig_set_pq gig_set_vQR gig_set_vcvtype \
    gig_print gig_estimate \
    gig_var_fcast \
    gig_plot gig_dplot gig_vfgraph \
    gig_bundle_print GUI_gig_plot

bundle-print = gig_bundle_print
bundle-plot = GUI_gig_plot
plot-precheck = gig_plot_precheck
gui-main = GUI_gig
label = gig
tags = C22
help = gig.pdf
sample-script = examples/example1.inp
min-version = @VERSION@
data-requirement = needs-time-series-data
data-files = examples
