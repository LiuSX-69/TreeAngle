# TreeAngle: A Geometric Statistic for Quantifying Spatial-Temporal Association in Cell Lineage Trees

**TreeAngle** is an R package that implements a geometric statistic—the **Minimal Enclosing Angle ($\Delta \theta$)**—designed to quantify the spatial-temporal association in tree-shaped datasets, specifically for cell lineage gene expression data.

Based on the methodology proposed by **Mao, Fan, Hu, and Liu (2025)**, this approach provides a robust, topology-aware alternative to traditional vector-based correlation measures (e.g., Pearson correlation) by accounting for the unique dependencies inherent in cell division and differentiation processes.

## Theoretical Background

In developmental biology, gene expression data measured on cell lineage trees (e.g., *C. elegans* embryogenesis) possess inherent structural properties:
*   **Temporal Dependence:** Expression values accumulate along developmental paths.
*   **Spatial Dependence:** Hierarchical relationships exist between parent and daughter cells.
*   **Correlation Damping:** Association strength between genes often weakens as cells differentiate across generations.

Traditional methods that "flatten" the tree into vectors fail to capture these dynamics. **TreeAngle** addresses this by employing the **TD-DeltaTheta Algorithm** (Algorithm 1 in the paper), which consists of:

1.  **Generation-Dependent Normalization:**
    *   **Scaling:** Data at each generation $i$ is scaled using the standard deviation of the **last generation** ($T$) as a global reference to preserve relative variance trends.
    *   **Shifting:** Data is shifted away from the origin dynamically. The shift distance increases with generation depth $i$ according to a harmonic series sum ($\sum 1/k$), separating later generations further in the geometric space.
2.  **Geometric Statistic:** Calculating the minimal angle $\Delta \theta$ that encloses the normalized data's quantile ellipse.

## Installation

You can install the development version of TreeAngle from GitHub using the `devtools` package:

```r
# install.packages("devtools")
devtools::install_github("LiuSX-69/TreeAngle")
```
## Usage Example

The core function `treeangle` calculates the geometric statistic. To correctly apply the normalization algorithm described in the paper, the input data must be structured to distinguish between generations.

### 1. Data Preparation (List Format)
The most convenient way to pass data is a **List**, where each element represents the data matrix ($N \times 2$) for one generation.

Here, we simulate data using `MASS::mvrnorm` to mimic gene co-expression that **damps (weakens)** over generations, matching the simulation logic in the original study.

```r
library(TreeAngle)
library(MASS) 

# --- Simulation: Generating Correlated Tree Data ---
set.seed(123)
tree_data_list <- list()
n_samples <- 100

# Generation 1: Strong Association (rho = 0.9)
sigma1 <- matrix(c(1, 0.9, 0.9, 1), 2, 2)
tree_data_list[[1]] <- mvrnorm(n_samples, mu = c(0,0), Sigma = sigma1)

# Generation 2: Damped Association (rho = 0.7)
sigma2 <- matrix(c(1, 0.7, 0.7, 1), 2, 2)
tree_data_list[[2]] <- mvrnorm(n_samples, mu = c(2,2), Sigma = sigma2)

# Generation 3 : Weak Association (rho = 0.5)
sigma3 <- matrix(c(1, 0.5, 0.5, 1), 2, 2) 
tree_data_list[[3]] <- mvrnorm(n_samples, mu = c(4,4), Sigma = sigma3)

# Structure: List of 3 matrices
str(tree_data_list) 
```

### 2. Calculation
Calculate the angle using the strict Algorithm 1 normalization.

```r

result <- treeangle(tree_data_list, normalize = TRUE, tau = 0.1)
print(paste("Tree Association Angle:", round(result[1, 1], 4)))
print(result)
```

### 3. Alternative Input (Matrix with Column)
Alternatively, you can provide a single matrix with a 3rd column indicating the **Generation**.

```r

mat_data <- rbind(
  cbind(tree_data_list[[1]], 1),
  cbind(tree_data_list[[2]], 2),
  cbind(tree_data_list[[3]], 3)
)

result_mat <- treeangle(mat_data, normalize = TRUE)
print(result_mat[1,1])
```

## Parameters

| Argument | Default | Description | Correspondence in Paper |
| :--- | :--- | :--- | :--- |
| `data` | - | **List** of matrices (where `data[[i]]` is Gen $i$) <br>OR **Matrix** with 3 columns (X, Y, Gen). | Represents the paired **increments** of gene expression across different tree depths. Using a List is recommended to naturally represent generational data structure. |
| `alpha` | `0.95` | Numeric (0-1). | Corresponds to the **$(1-\alpha)$ quantile ellipse**. `0.95` filters out the top 5% outliers to ensure robust estimation. |
| `method` | `"mean"` | String. | The strategy to aggregate candidate angles. The paper recommends `"mean"` for the estimator $\Delta \hat{\theta}$. Options: `"min"`, `"median"`, `"equal"`, `"mean"`, `"neigh"`. |
| `normalize`| `TRUE` | Logical. | Whether to apply the **Algorithm 1** normalization (scaling by last gen SD, shifting by harmonic series). |
| `target_sd`| `0.01` | Numeric. | The target standard deviation after scaling. This matches the term $\sqrt{10^{-4}}$ used in the original simulation code. It standardizes the spread of data clouds. |
| `tau` | `0.1` | Numeric. | A hyperparameter $\tau$ that controls the magnitude of the **Dynamic Mean Shift** $\mu_i^\ast$. <br>Formula: $\mu_i^\ast = \sqrt{-2\ln(1-\alpha)\sigma_{\text{target}}^2 + \tau \sum_{k=1}^i \frac{1}{k}}$ <br>Larger `tau` moves later generations further from the origin, reducing overlap. |
| `start` | `c(0,0)` | Numeric vector. | The **root point** ($x_0, y_0$) from which the rays emanate. Typically set to the origin for normalized increments. |

## Reference

If you use **TreeAngle** in your research, please cite the original paper:

> **Mao, S., Fan, X., Hu, J., & Liu, S. (2025).** "A Geometric Statistic for Quantifying the Spatial-Temporal Association in Cell Lineage Gene Expression Tree."
