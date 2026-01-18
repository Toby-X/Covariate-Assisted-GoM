# Covariate-assisted Grade of Membership Models via Shared Latent Geometry

This repository contains the R code for the paper "Covariate-assisted Grade of Membership Models via Shared Latent Geometry".

## Installation

Install the required packages using the following command in your R console:
```R
install.packages(readLines("requirements.txt"))
```

## Usage

To run the covariate-assisted GoM, use the following command in your R console:
```R
source("cogom.r")
```

The code for the three simulations in the paper is in the `simulations` directory. The code and processed data for the real data analysis is in the `real_data` directory.

### Quick Start

**Input Format:**

$\mathbf{R}$: $N \times J$ binary response matrix.

$\mathbf{X}$: $N \times W$ covariate matrix (continuous).

$K$: Number of latent profiles.

```R
# 1. Load your data (Example with random generation not following the model assumptions)
N <- 200; J <- 50; W <- 20; K <- 3
R_mat <- matrix(rbinom(N * J, 1, 0.5), nrow = N) # Your response matrix
X_mat <- matrix(rnorm(N * W), nrow = N)          # Your covariate matrix

# 2. Select tuning parameter alpha (optional cross-validation step)
# alpha_seq <- seq(0, 1, length=20)
# alpha_opt <- find_best_alpha(alpha_seq, R_mat, X_mat, K)
alpha_opt <- 0.5

# 3. Run the estimator
res <- gom.cov(R_mat, X_mat, K, alpha_opt)

# 4. Access estimates
# Pi: N x K Mixed Membership Matrix
# Theta: J x K Response Probability Matrix
print(head(res$Pi))    
print(head(res$Theta)) 
```