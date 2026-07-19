author = Artur Tarassow and Allin Cottrell
email = atecon@posteo.de
version = @VERSION@
date = 2026-07-12
description = K-means clustering algorithm
tags = C13 C52
min-version = @VERSION@
public = kmeans_fit kmeans_predict kmeans_plot kmeans_summary \
    kmeansGUI kmeans_scree kmeans_screeplot kmeans_crosstab \
    kmeans_sil_samples kmeans_sil_score kmeans_sil_plot
help = kmeans.pdf
sample-script = kmeans_sample.inp
depends = PairPlot

gui-main = kmeansGUI
bundle-plot = GUI_kmeans_plot
plot-precheck = kmeans_plot_precheck

label = K-Means
menu-attachment = MAINWIN/View
data-files = examples
