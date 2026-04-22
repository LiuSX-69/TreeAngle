# TreeAngle

`TreeAngle` is an R package for computing the **TreeAngle** statistic, a geometric summary of association for paired generation-structured data. The package was motivated by cell lineage gene expression analysis, but it can also be used for other paired data organized by generation.

The package is intentionally lightweight and uses minimal dependencies, which makes it easier to install and test in a clean R environment.

---

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

---

## Quick start

```r
library(TreeAngle)

sim <- simdata(seed = 1)
res <- treeangle(sim, alpha = 0.95, method = "mean")

print(res)
as.matrix(res)
```

---

## Important note on `alpha`

In the manuscript, the symbol `alpha` may denote a **tail probability**, so the corresponding coverage is `1 - alpha`.

In this package, for user convenience and consistency with the original simulation scripts, the argument `alpha` denotes the **coverage proportion**.

Therefore:

- `alpha = 0.95` means that the estimated angle encloses approximately **95%** of points;
- this corresponds to a manuscript-style tail probability of about **0.05**.

---

## Main function

```r
treeangle(
  data,
  alpha = 0.95,
  method = "mean",
  normalize = NULL,
  target_sd = 0.01,
  tau = 0.1
)
```

### Arguments

#### `data`
Input data. Supported formats are:

1. a **list** of generation-specific two-column matrices or data frames;
2. a **matrix/data frame** with columns `x`, `y`, and `generation`;
3. a direct **two-column** matrix/data frame.

Objects returned by `simdata()` can be passed directly to `treeangle()`.

#### `alpha`
Coverage proportion used to define the angle.

- `alpha = 0.95` means that about 95% of points are enclosed;
- about 5% of points are allowed outside.

#### `method`
Summary rule used to choose the final angle.

Available values are:

- `"mean"`: average of valid candidate angles; recommended;
- `"min"`: smallest candidate angle;
- `"median"`: median candidate angle;
- `"equal"`: central candidate rule;
- `"neigh"`: neighborhood-adjusted candidate rule.

#### `normalize`
Controls whether normalization is applied.

- `NULL` (default): normalize automatically for generation-structured input and skip normalization for direct two-column input;
- `TRUE`: force normalization;
- `FALSE`: skip normalization.

#### `target_sd`
Target marginal standard deviation used in normalization. Default is `0.01`.

#### `tau`
Generation shift parameter used in normalization. Default is `0.1`.

---

## Input formats

### 1. Generation-wise list input

```r
sim <- simdata(
  generations = 6,
  rho = 0.7,
  decay = "linear",
  format = "list",
  seed = 123
)

res <- treeangle(sim, alpha = 0.95, method = "mean")
print(res)
```

### 2. Matrix/data frame with generation column

```r
sim <- simdata(
  generations = 6,
  rho = 0.7,
  decay = "power",
  format = "matrix",
  seed = 123
)

res <- treeangle(
  sim,
  alpha = 0.95,
  method = "median",
  normalize = TRUE
)

print(res)
```

### 3. Direct two-column input

```r
sim <- simdata(format = "xy", seed = 123)

res <- treeangle(
  sim,
  alpha = 0.95,
  method = "min",
  normalize = FALSE
)

print(res)
```

---

## Output

`treeangle()` returns an object of class `treeangle_result`.

Important components include:

- `angle_degrees`
- `tan_angle`
- `upper_point`
- `lower_point`
- `upper_line`
- `lower_line`
- `method_table`
- `result_matrix`

The legacy matrix representation can be recovered by:

```r
as.matrix(res)
```

All five method summaries can be inspected with:

```r
res$method_table
```

---

## Simulating data

`simdata()` generates paired Gaussian data generation by generation and is designed to work directly with `treeangle()`.

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

treeangle(sim, alpha = 0.95, method = "mean")
```

### Decay options

- `decay = "power"` uses `rho^i`
- `decay = "linear"` uses a generation-wise linear decay

---

## Reproducibility and package checks

Before submission or resubmission, run the package in a **clean R session**.

```r
devtools::document()
devtools::test()
devtools::check()
```

A minimal smoke test is:

```r
library(TreeAngle)

sim <- simdata(seed = 1)
res <- treeangle(sim, alpha = 0.95, method = "mean")

stopifnot(
  inherits(res, "treeangle_result"),
  is.finite(res$angle_degrees)
)
```

---

## Citation

If you use `TreeAngle`, please cite the associated manuscript:

> Mao, S., Fan, X., Hu, J., and Liu, S. (2025).  
> *A Geometric Statistic for Quantifying the Spatial-Temporal Association in Cell Lineage Gene Expression Data.*
