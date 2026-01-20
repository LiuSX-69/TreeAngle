# TreeAngle: A Geometric Statistic for Quantifying Spatial-Temporal Association in Cell Lineage Trees

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**TreeAngle** is an R package that implements a geometric statistic—the **Minimal Enclosing Angle ($\Delta \theta$)**—designed to quantify the spatial-temporal association in tree-shaped datasets, specifically for cell lineage gene expression data.

Based on the methodology proposed by **Mao, Fan, and Hu (2025)**, this approach provides a robust, topology-aware alternative to traditional vector-based correlation measures (e.g., Pearson correlation) by accounting for the unique dependencies inherent in cell division and differentiation processes.

## Theoretical Background

In developmental biology, gene expression data measured on cell lineage trees (e.g., *C. elegans* embryogenesis) possess inherent structural properties:
*   **Temporal Dependence:** Expression values accumulate along developmental paths.
*   **Spatial Dependence:** Hierarchical relationships exist between parent and daughter cells.
*   **Correlation Damping:** Association strength may vary across generations.

Traditional methods that "flatten" the tree into vectors fail to capture these dynamics. **TreeAngle** addresses this by:
1.  Mapping paired gene expression **increments** to a 2D plane.
2.  Constructing a **Quantile Ellipse** (based on a significance level $\alpha$).
3.  Calculating the minimal angle $\Delta \theta$ that encloses the ellipse.

**Interpretation:**
*   **Small Angle ($\Delta \theta$)**: Indicates **strong association** (high synchrony in expression trends).
*   **Large Angle ($\Delta \theta$)**: Indicates **weak association** or independence.

## Installation

You can install the development version of TreeAngle from GitHub using the `devtools` package:

```r
# install.packages("devtools")
devtools::install_github("LiuSX-69/TreeAngle")
```


## Usage Example

The core function `calangle` calculates the geometric statistic from 2D data points.

### 1. Basic Calculation
Here is how to calculate the angle for a simulated dataset representing paired gene expression increments.

```r
library(TreeAngle)

# --- Simulation ---
# Generate synthetic data representing paired increments (normalized)
set.seed(123)

# Case A: Strong Association (points cluster tightly along the diagonal)
# Represents gene pairs with high co-expression
data_strong <- matrix(rnorm(200, mean=0, sd=1), ncol=2)
data_strong[, 2] <- 0.9 * data_strong[, 1] + 0.1 * rnorm(100) 

# Case B: Weak Association (points are dispersed)
# Represents unrelated genes
data_weak <- matrix(rnorm(200, mean=0, sd=1), ncol=2)

# --- Calculation ---
# Calculate the geometric statistic using default parameters (alpha=0.95)
result_strong <- calangle(data_strong)
result_weak   <- calangle(data_weak)

# --- Interpretation ---
# The first element [1,1] of the result matrix is the Angle statistic.
print(paste("Strong Association Angle:", round(result_strong[1, 1], 2)))
print(paste("Weak Association Angle:  ", round(result_weak[1, 1], 2)))
```

### 2. Advanced Options
You can customize the quantile threshold (`alpha`) and the estimation method based on your specific analysis requirements.

```r
# Use a custom alpha (e.g., enclosing 90% of data instead of 95%)
# Use 'min' method to find the absolute minimal candidate angle
res_custom <- calangle(data_strong, alpha = 0.90, method = "min")

# The output is a matrix containing:
# Row 1: The Angle and Tangent values
# Rows 2-5: Coordinates of the defining points and lines (useful for plotting)
print(res_custom)
```

## Parameters

| Argument | Default | Description | Correspondence in Paper |
| :--- | :--- | :--- | :--- |
| `data` | - | A matrix or data frame of 2 columns ($X, Y$). | Represents the paired **increments** of gene expression (after normalization). |
| `alpha` | `0.95` | Numeric value between 0 and 1. | Corresponds to the **$(1-\alpha)$ quantile ellipse**. `0.95` filters out the top 5% outliers to ensure robust estimation. |
| `method` | `"mean"` | Character string. Options: `"min"`, `"median"`, `"equal"`, `"mean"`, `"neigh"`. | The strategy to aggregate candidate angles from the dataset (see **Supplementary Fig S.5**). The paper recommends `"mean"` for the estimator $\Delta \hat{\theta}$. |
| `start` | `c(0,0)` | Numeric vector of length 2. | The **root point** ($x_0, y_0$) from which the rays emanate. Typically set to the origin for normalized increments. |

## Reference

If you use **TreeAngle** in your research, please cite the original paper:

> **Mao, S., Fan, X., & Hu, J. (2025).** "A Geometric Statistic for Quantifying the Spatial-Temporal Association in Cell Lineage Gene Expression Tree."