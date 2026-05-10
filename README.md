# TreeAngle

`TreeAngle` is an R package for computing the **TreeAngle** statistic, a geometric summary of association for paired generation-structured data.

The package was motivated by cell lineage gene expression analysis, but the method can also be used for other paired data organized by generation.

## Installation

### From GitHub

```r
# install.packages("remotes")
remotes::install_github("LiuSX-69/TreeAngle")
```

### From a local source tarball

```r
install.packages("TreeAngle_<version>.tar.gz", repos = NULL, type = "source")
```

## Main arguments

| argument | description |
|---|---|
| `data` | paired data, supplied as a generation-wise list, a matrix/data frame with `x`, `y`, and optionally `generation`, or a direct two-column object |
| `alpha` | coverage proportion used to define the angle |
| `method` | rule used to choose the final angle |
| `normalize` | whether to apply generation-wise normalization |
| `target_sd` | target marginal standard deviation after normalization |
| `tau` | generation shift parameter used in normalization |

### `alpha`

In this package, `alpha` denotes the coverage proportion.

```r
alpha = 0.95
```

means that the angle is estimated to enclose approximately 95% of the points.
In some manuscript notation, `alpha` may instead denote a tail probability. Under that convention, package-level `alpha = 0.95` corresponds to a tail probability of about `0.05`.

### `method`

The `method` argument specifies how the final reported angle is chosen from these valid candidates.
Available options are:
- `"mean"`: average of the valid candidate angles; package default;
- `"min"`: smallest valid candidate angle, giving the tightest admissible angle;
- `"median"`: median of the valid candidate angles;
- `"equal"`: candidate chosen by the equal-allocation rule in the original implementation;
- `"neigh"`: candidate chosen by the neighborhood-adjusted rule in the original implementation.
- `"all"`: return all five original outputs:
  `anglemin`, `anglemedian`, `angleequal`, `anglemean`, and `angleneigh`.

### Output
For `method = "mean"`, `"min"`, `"median"`, `"equal"`, or `"neigh"`,
`treeangle()` returns a plain `5 x 2` numeric matrix.
The rows are:
| row | meaning |
|---|---|
| `angle` | angle summary |
| `pointup` | point on the upper boundary ray |
| `pointdown` | point on the lower boundary ray |
| `coefup` | slope and intercept of the upper boundary line |
| `coefdown` | slope and intercept of the lower boundary line |
The columns are:
| column | meaning |
|---|---|
| `deg_x_slope` | angle in degrees for row `angle`; x-coordinate for point rows; slope for coefficient rows |
| `tan_y_intercept` | tangent of the angle for row `angle`; y-coordinate for point rows; intercept for coefficient rows |
For `method = "all"`, `treeangle()` returns a list with five such matrices:
`anglemin`, `anglemedian`, `angleequal`, `anglemean`, and `angleneigh`.


### `normalize`

`normalize` controls whether generation-wise normalization is applied.

- `normalize = NULL`: automatic behavior;
- `normalize = TRUE`: force normalization;
- `normalize = FALSE`: skip normalization.

For generation-structured data, normalization scales the paired coordinates generation by generation before computing the angle.

The default target marginal standard deviation is:

```r
target_sd = 1
```

When normalization is used:

- each generation must contain at least two observations;
- both coordinates must have non-zero empirical standard deviation within each generation;
- generation labels for matrix/data frame input must be consecutive positive integers starting from 1.

Valid generation labels:

```r
1, 2, 3, 4
```

Invalid generation labels:

```r
2, 3, 4
1, 2, 4
```

If the data have already been normalized or are already on a common scale, use:

```r
normalize = FALSE
```
## Data format and example data


In the package data-frame interface, the recommended format is a data frame
with three columns:

- `generation`: generation index;
- `x`: first paired coordinate;
- `y`: second paired coordinate.

Each row should represent one paired observation or one paired increment.

The easiest way to try the package is to first generate simulated data with
`simdata()` and then pass the result to `treeangle()`.

```r
dat <- simdata(
  generations = 5,
  branching = 2,
  observations_per_node = 1,
  mu = 7,
  sigma2 = 0.02,
  rho = 0.5,
  decay = "power",
  seed = 1
)

```

The generated data frame has the columns:

```r
names(dat)
#> [1] "generation" "x"          "y"
```

The generation sizes follow the specified branching structure:

```r
table(dat$generation)
#>  1  2  3  4  5 
#>  2  4  8 16 32
```

Now compute the TreeAngle statistic:

```r
res <- treeangle(
  dat,
  alpha = 0.95,
  method = "mean",
  normalize = TRUE
)

res
```

The result is a plain numeric matrix with row and column names.

```r
rownames(res)
#> [1] "angle"     "pointup"   "pointdown" "coefup"    "coefdown"

colnames(res)
#> [1] "deg_x_slope"     "tan_y_intercept"
```

The rows and columns have the following meanings:

| row | first column: `deg_x_slope` | second column: `tan_y_intercept` |
|---|---|---|
| `angle` | angle in degrees | tangent of the angle |
| `pointup` | x-coordinate of the upper boundary point | y-coordinate of the upper boundary point |
| `pointdown` | x-coordinate of the lower boundary point | y-coordinate of the lower boundary point |
| `coefup` | slope of the upper boundary line | intercept of the upper boundary line |
| `coefdown` | slope of the lower boundary line | intercept of the lower boundary line |

To return all five original summaries, use:

```r
all_res <- treeangle(
  dat,
  alpha = 0.95,
  method = "all",
  normalize = TRUE
)

names(all_res)
#> [1] "anglemin"    "anglemedian" "angleequal"  "anglemean"   "angleneigh"
```

`TreeAngle` does not reconstruct paired observations from raw lineage or tree files. If starting from raw lineage or tree data, first construct the paired observations or paired increments, then pass the resulting table to `treeangle()`.


### Decay options

- `decay = "power"` uses `rho^i`;
- `decay = "linear"` uses the manuscript-style linear rule  
  `rho * (1 - (i - 1) / generations)`;


## Accepted input objects

`treeangle()` accepts several input representations, but the recommended interface is a data frame with columns `generation`, `x`, and `y`.
Other supported formats, such as generation-wise lists or two-column numeric matrices, are mainly provided for compatibility with the original implementation.

## Citation

If you use `TreeAngle`, please cite the associated manuscript:

> Mao, S., Fan, X., Hu, J., and Liu, S. (2025).  
> *A Geometric Statistic for Quantifying the Spatial-Temporal Association in Cell Lineage Gene Expression Data.*