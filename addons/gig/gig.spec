author = Riccardo "Jack" Lucchetti and Stefano Balietti
email = r.lucchetti@univpm.it
version = 2.16
date = 2014-12-12
description = An assortment of univariate GARCH models
public = GUI_gig \
    gig_setup gig_set_dist gig_set_pq gig_set_vQR gig_set_vcvtype \
    gig_print gig_estimate \
    gig_plot gig_dplot \
    gig_bundle_print GUI_gig_plot

bundle-print = gig_bundle_print
bundle-plot = GUI_gig_plot
gui-main = GUI_gig
label = gig
help = gig.pdf
sample-script = examples/example1.inp
min-version = 1.9.92
data-requirement = needs-time-series-data
data-files = examples
