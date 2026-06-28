author = Artur Tarassow and Allin Cottrell
email = atecon@posteo.de
version = @VERSION@
date = 2026-06-24
description = K-means clustering algorithm
tags = C13 C52
min-version = @VERSION@
public = kmeans_fit kmeans_predict kmeans_plot kmeans_summary \
    kmeansGUI kmeans_scree kmeans_screeplot \
    kmeans_sil_samples kmeans_sil_score kmeans_sil_plot
help = kmeans_help.md
sample-script = kmeans_sample.inp
depends = PairPlot
gui-main = kmeansGUI
#menu-only = kmeansGUI
label = K-Means
menu-attachment = MAINWIN/View
data-files = examples
