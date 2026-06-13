# kmeans Package

This package includes functionalities to identify unknown clusters of multidimensional data using the well known (at least in the machine-learning field) k-means algorithm.

The k-means algorithm divides a set of `N` samples `X` into `k` disjoint clusters `C`, each described by the mean of the samples in the cluster. The means are called the cluster centroids.

The objective is to minimize some loss. For instance, the objective is to minimize "inertia", or within-cluster sum-of-squares criterion in case of the Euclidean distance function.

For more information see:

https://scikit-learn.org/stable/modules/clustering.html#kmeans

Additionally, the package provides functions to compute silhouette samples and scores, which are measures of how similar an observation is to its own cluster compared to other clusters. The silhouette score is a measure of how well clusters are separated; higher values indicate better-defined clusters. The best value is 1 and the worst value is -1. Values near 0 indicate overlapping clusters. Negative values generally indicate that a sample has been assigned to the wrong cluster, as a different cluster is more similar. The reference for silhouette analysis is:

Rousseeuw, P. J. (1987). Silhouettes: A graphical aid to the interpretation and validation of cluster analysis. *Journal of Computational and Applied Mathematics*, 20, 53-65, URL: https://www.sciencedirect.com/science/article/pii/0377042787901257

Please ask questions and report bugs on the Gretl mailing list if possible. Alternatively, create an issue ticket on the github repo (see below).
Source code and test script(s) can be found here: https://github.com/atecon/kmeans


## GUI access

Currently, only k-means estimation but not silhouette analysis can be performed via the dialog box, which can be opened via `View -> k-Means`.

Two choices are available:

1. If you activate the checkbox "Show scree plot" the scree plot will be shown, and the parameter `n_clusters` is the maximum number of clusters to evaluate and plot. The returned bundle will contain only a matrix with two columns: the first column is the number of clusters, and the second column is the within-cluster sum of squares (*inertia*). Each row corresponds to a different number of clusters.

2. If you activate the checkbox "Show pair-plot" the pair-plot will be shown, and the parameter `n_clusters` is the number of clusters to estimate. The returned bundle will contain detailed information about the estimation.


# Public Functions

## kmeans_fit

```
kmeans_fit (list xlist, int n_clusters[2::2], bundle opts[null])
```

Execute the kmeans algorithm and estimate the clusters. It is required that `xlist` does not contain missing values. So, please make sure to clean the data before calling this function.

**Arguments:**

- `xlist`: list, Features (regressors) to train the model.
- `n_clusters`: int, Number of assumed clusters (default: 2)
- `opts`: bundle, Optional parameters affecting the kmeans algorithm. **You can pass the following parameters:**

    * `algorithm`: string, kmeans algorithm to use. Currently, only `full` is supported (classical EM-style algorithm).
    * `distance_type`: string, Name of the distance metric applied (default: `euclidean`). For more distance metrics, see gretl's built-in function `distance()`.
    * `initializer`: string, Method for initialization. Either `random`: Choose `n_clusters` observations (rows) at random from data for the initial centroids. Or `pca`: Try to pick data points that are as far apart as possible by means of PCA.
    * `max_iter`: int, Maximum number of iterations of the kmeans algorithm to run.
    * `n_draws`: int, Number of time the kmeans algorithm will be run with different centroid seeds. The final results will be the best output of `n_draws` consecutive runs in terms of inertia.
    * `tolerance`:  scalar, Minimum improvement of the `within_variation_total` (Sum of the squared distances across all clusters) required before early stopping the algorithm (default: 1e-4)
    * `verbose`: int, Level of verbosity: `0`: don't print anything, `1`: print some details, `2`: print more details (default: `0`)


**Return:** Bundle holding various items.

- `between_variation`: scalar, Between cluster sum of squares = `total_ssq - within_variation_total`
- `centroids`: matrix, Estimated mean values (centroids) for each feature (columns) and for each cluster (rows).
- `cluster_id`: matrix, Estimated cluster ID for each observation for the best draw minimizing `inertia`.
- `distances`: matrix, Estimated distance for the best draw minimizing `inertia`.
- `error`: int, Error code. In case of no error `FALSE`, otherwise positive integer.
- `nobs`: int, Number of non-missing observations used for training.
- `pointsize`: scalar, Size of points being plotted when calling the `kmeans_plot()` function.
- `total_ssq`: scalar, Sum of the squared distances of the features from its mean values
- `use_circles`: bool, Plot circles instead of point when calling the `kmeans_plot()` function.
- `within_variation_total`: scalar, Sum of the squared distances across all clusters (*inertia*).
- `within_variation_avg`: scalar, Sum of the average squared distances across all clusters.


## kmeans_predict

```
kmeans_predict (list xlist, bundle Model)
```

Predict cluster belonging based on the estimated model.

**Arguments:**

- `xlist`: list, Features (regressors) used for predicting cluster belonging.
- `Model`: bundle, Model object returned by the `kmeans_fit()` function.

**Return:** Series holding the predicted cluster ID for each observation.


## kmeans_summary

```
kmeans_summary (bundle Model)
```

Print summarizing information on estimation step after having applied the `kmeans_fit()` function.

**Arguments:**

- Model: bundle, Bundle returned by the `kmeans_fit()` function.

**Return:** Nothing.


## kmeans_plot

```
kmeans_plot (list xlist, bundle self[null])
```

Factorized scatter plot estimated clusters for each 2-dimensional combination of features. This function calls the user-defined package "PairPlot" which must be installed.

**Arguments:**

- `xlist`: list, Features (regressors) used for plotting.
- `self`: bundle, Bundle for manipulating the plot. **Note** Here you can also pass options accepted by the "PairPlot" package which is used in the background.

**Return:** Nothing.


## kmeans_screeplot

```
kmeans_screeplot (list xlist, int max_clusters, string filename[null], bundle self[null])
```

This function plots the scree plot for the kmeans algorithm. The scree plot shows the within-cluster sum of squares (*inertia*) for different numbers of clusters (from 1 to `max_clusters`). The optimal number of clusters is the one where the within-cluster sum of squares starts to decrease more slowly.

**Arguments:**

- `xlist`: list, Features (regressors) used for plotting.
- `max_clusters`: int, Maximum number of clusters to plot (default: 3).
- `filename`: string, Name of the file to save the plot (optional). If not provided, the plot will be shown on the screen.
- `self`: bundle, Bundle for manipulating the plot (optional). **Note, accepted options are:**

    * The same as for the `kmeans_fit()` function for the `opts` bundle. Relevant for controlling the k-means algorithm.
    * `verbose`: int, If `2` show summary for each model, else show nothing.
    * `fontsize`: scalar, Font size for the plot.
    * `linewidth`: scalar, Line width for the plot.

**Return:** Matrix with two columns: the first column is the number of clusters, and the second column is the within-cluster sum of squares (*inertia*). Each row corresponds to a different number of clusters.


## kmeans_sil_samples

```
kmeans_sil_samples (bundle Model, list xlist)
```

Compute the silhouette sample value for each observation in the fitted kmeans model. The silhouette value measures how similar an observation is to its own cluster compared to other clusters.

**Reference:** https://en.wikipedia.org/wiki/Silhouette_(clustering)

**Arguments:**

- `Model`: bundle, Model object returned by the `kmeans_fit()` function.
- `xlist`: list, Features (regressors) used for clustering.

**Return:** Series holding the silhouette value for each observation.


## kmeans_sil_score

```
kmeans_sil_score (bundle Model, list xlist, bool global[TRUE])
```

Compute the mean silhouette score for all samples in the fitted kmeans model if the optional parameter `global` is set to `TRUE` (default). The silhouette score is a measure of how well clusters are separated; higher values indicate better-defined clusters. The best value is 1 and the worst value is -1. Values near 0 indicate overlapping clusters. Negative values generally indicate that a sample has been assigned to the wrong cluster, as a different cluster is more similar. If `global` is set to `FALSE`, the mean silhouette score is computed for each cluster separately.

**Reference:** https://scikit-learn.org/stable/modules/generated/sklearn.metrics.silhouette_score.html

**Arguments:**

- `Model`: bundle, Model object returned by the `kmeans_fit()` function.
- `xlist`: list, Features (regressors) used for clustering.
- `global`: bool, If TRUE (default) compute the mean silhouette score for all samples, otherwise compute it for each cluster separately.

**Return:** Scalar value representing the mean silhouette score if `global` is set to TRUE, otherwise a matrix is returned  (columns refer to the cluster and the score).


## kmeans_sil_plot

```
kmeans_sil_plot (bundle Model, list xlist, bundle opts[null])
```

Create a silhouette plot for the fitted kmeans model. Each cluster is plotted separately, showing the silhouette value for each observation. Plot options can be passed via the `opts` bundle. The dashed vertical line shows the mean Silhouette score.

**Reference:** https://en.wikipedia.org/wiki/Silhouette_(clustering)

**Arguments:**

- `Model`: bundle, Model object returned by the `kmeans_fit()` function.
- `xlist`: list, Features (regressors) used for clustering.
- `opts`: bundle, Optional parameters for controlling the plot. The following options are supported:
    * `filename`: string, Name of the file to save the plot. If not provided, the plot will be shown on the screen.
    * `width`: scalar, Width of the plot in pixels (default: 600).
    * `height`: scalar, Height of the plot in pixels (default: 900).
    * `columns`: scalar, Number of columns in the grid layout for clusters (default: 1).
    * `base_margin`: scalar, Base margin for top and bottom plot margins (default: 2.5).
    * `pointsize`: scalar, Size of points (default: 0.75)

**Return:** Nothing.


# Changelog

* **v0.6 (August 2025)**
    * Add new silhouette functions: `kmeans_sil_samples()`, `kmeans_sil_score()`, and `kmeans_sil_plot()`
    * Bugifx: Correctly handle the case for the "mahalanobis" distance metric in the `compute_distance()` function. (Kudos to Allin Cottrell for reporting and providing a fix)
    * GUI dialog: Add option to plot silhouette analysis and general improvements
    * Raise minimum version of gretl to 2023c

* **v0.5 (February 2025)**
    * Bugfix: in case of a single cluster, switch to "random" initializer as "pca" would fail; throw error if pca initializer is called for a single cluster
    * Improvement: Improve robustness if empty clusters are estimated (e.g., due to bad random initialization of centroids)

* **v0.4 (January 2025)**
    * Add the new function `kmeans_screeplot()` for plotting the scree plot (method to determine the optimal number of clusters).
    * Change default values for parameters `pointsize` and `fontsize` in the `kmeans_plot()` function. Now set by gretl's default values.
    * Bugfix: The distance measure is now taken into account when calling the function via the GUI. Before, the default distance measure was always used.

* **v0.3 (February 2024)**
    * Add GUI dialog
    * Move to markdown-based help file
    * Internal improvements

* **v0.2 (July 2022)**
    * Fix bug that arises if the sample range is restricted, and you're trying to coerce a column vector that's not the full length of the dataset into a series on adding it to a bundle.
    * Returned objects `cluster_id` and `distances` when calling the `kmeans_fit()` function are of type matrix instead of series, now.

* **v0.1 (February 2022)**
    * initial release
