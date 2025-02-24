# Load necessary package for GitHub installation
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# Function to install GitHub packages if not installed
install_github_if_missing <- function(repo, pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    devtools::install_github(repo)
  }
}

## Install not installed Cran packages
cran_packages <- c(
  "tidyverse",
  "here",
  "janitor",
  "rsample",
  "MASS",
  "imputeR",
  "palmerpenguins",
  "LassoSIR",
  "TULIP",
  "rstudioapi",
  "miscset",
  "ggh4x"
)
not_installed <- cran_packages[!(cran_packages %in% rownames(installed.packages()))]

if (length(not_installed) > 0) {
  install.packages(not_installed)
}

# List of GitHub packages (repo = "username/repository", package_name = "package")
github_packages <- list(
  list(repo = "D3r1kBoonstra/sdr", pkg = "sdr"),
  list(repo = "ywwry66/QDA-by-Projection-R-Package", pkg = "QDA-by-Projection-R-Package")
)


# Install missing GitHub packages
lapply(github_packages, function(x) install_github_if_missing(x$repo, x$pkg))
