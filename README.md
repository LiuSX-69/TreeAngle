# TreeAngle

**TreeAngle** is an R package for computing the geometric **TreeAngle**
statistic from paired generation-structured data, especially data arising from
cell lineage studies.

The package is written for clean-environment use:
- no absolute file paths,
- no `source()` calls,
- no hidden workspace dependence.

## What the package does

- computes the TreeAngle statistic from paired two-dimensional data;
- optionally applies the manuscript normalization step when generation
  information is available;
- supports five summary rules:
  - `"mean"` (recommended),
  - `"min"`,
  - `"median"`,
  - `"equal"`,
  - `"neigh"`;
- generates synthetic manuscript-style data with `simdata()`.

---

## Important input note

`treeangle()` expects **paired 2D observations**, optionally grouped by
generation.

In the manuscript, these inputs are typically **generation-wise paired
increments after preprocessing**, rather than a raw nested tree object.

So the package accepts:

1. a list of generation-specific two-column matrices/data frames;
2. a matrix/data frame with columns `x`, `y`, and `generation`;
3. a direct two-column matrix/data frame.

---

## Installation

### From GitHub

```r
# install.packages("remotes")
remotes::install_github("LiuSX-69/TreeAngle")
```

### From a local source package

```r
install.packages("TreeAngle_1.0.1.tar.gz", repos = NULL, type = "source")
```

---

## Main functions

```r
treeangle(
  data,
  alpha = 0.95,
  method = "mean",
  normalize = NULL,
  target_sd = 0.01,
  tau = 0.1
)

simdata(...)
```

---

## `treeangle()` arguments

### `alpha`
Coverage proportion used to construct the angle.

- `alpha = 0.95` means the estimated angle encloses about 95% of points;
- about `1 - alpha = 0.05` of points are allowed outside.

### `method`
Summary rule used to choose the final angle.

Available values are:

- `"mean"`: average of all valid candidate angles; recommended;
- `"min"`: smallest candidate angle;
- `"median"`: median candidate angle;
- `"equal"`: central candidate pair;
- `"neigh"`: neighborhood-adjusted candidate around the minimum.

### `normalize`
Controls whether normalization is applied.

- `normalize = NULL` (default):
  - generation-structured input -> normalize automatically;
  - direct two-column input -> do not normalize.
- `normalize = TRUE`: force normalization; generation information is required.
- `normalize = FALSE`: skip normalization.

### `target_sd`
Target marginal standard deviation used in normalization.

Default:

```r
target_sd = 0.01
```

### `tau`
Generation shift parameter used in normalization.

Default:

```r
tau = 0.1
```

These defaults come from the original code/manuscript setting.

---

## Input formats

### 1. Generation-wise list input

```r
library(TreeAngle)

sim <- simdata(
  generations = 6,
  rho = 0.7,
  decay = "linear",
  format = "list",
  seed = 123
)

res <- treeangle(
  sim,
  alpha = 0.95,
  method = "mean"
)

print(res)
res$angle_degrees
```

### 2. Matrix/data frame with generation column

```r
sim_mat <- simdata(
  generations = 6,
  rho = 0.7,
  decay = "power",
  format = "matrix",
  seed = 123
)

res <- treeangle(
  sim_mat,
  alpha = 0.95,
  method = "mean",
  normalize = TRUE,
  target_sd = 0.01,
  tau = 0.1
)

print(res)
```

### 3. Direct two-column input

```r
xy <- do.call(rbind, sim$data)

res <- treeangle(
  xy,
  alpha = 0.95,
  method = "median",
  normalize = FALSE
)

print(res)
```

---

## Output

`treeangle()` returns a `treeangle_result` object.

Important components include:

- `angle_degrees`
- `tan_angle`
- `upper_point`
- `lower_point`
- `upper_line`
- `lower_line`
- `method_table`

The legacy 5 x 2 matrix form can be recovered by:

```r
as.matrix(res)
```

You can inspect all five method summaries from:

```r
res$method_table
```

---

## Simulating data

`simdata()` generates generation-wise paired Gaussian data that can be passed
directly to `treeangle()`.

```r
sim <- simdata(
  generations = 6,
  branching = 2,
  observations_per_node = 1,
  mu = 7,
  sigma2 = 0.02,
  rho = 0.5,
  decay = "power",
  format = "list",
  seed = 1
)

print(sim)

treeangle(
  sim,
  alpha = 0.95,
  method = "mean"
)
```

### Decay options

- `decay = "power"` uses `rho^i`, which matches the manuscript simulation setup.
- `decay = "linear"` reproduces the old root-level script logic more closely.

---

## Reproducing old script behavior

Your old root-level scripts used a linear generation decay and path-based
simulation outside the package. The package version replaces that workflow with:

```r
sim <- simdata(
  generations = 6,
  rho = 0.7,
  decay = "linear",
  format = "matrix",
  seed = 1
)

treeangle(
  sim,
  alpha = 0.95,
  method = "mean",
  normalize = TRUE,
  target_sd = 0.01,
  tau = 0.1
)
```

---

## Clean-environment check

Before resubmission, run in a fresh R session:

```r
devtools::document()
devtools::check()
```

This will regenerate `NAMESPACE` and `man/` and check for missing dependencies,
namespace problems, and undocumented arguments.

---

## Reference

If you use **TreeAngle**, please cite the associated manuscript:

> Mao, S., Fan, X., Hu, J., and Liu, S. (2025).  
> A Geometric Statistic for Quantifying the Spatial-Temporal Association in Cell Lineage Gene Expression Data.
```