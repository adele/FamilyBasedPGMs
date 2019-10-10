## FamilyBasedPGMs: An R Package for Learning Genetic and Environmental Graphical Models from Family Data

### Overview

This package provides methods for learning, from observational Gaussian family data (i.e., Gaussian data clusterized in families), Gaussian undirected and directed acyclic PGMs describing linear relationships among multiple phenotypes and a decomposition of the learned PGM into unconfounded genetic and environmental PGMs. 

The structure learning is based on zero partial correlation tests, derived in the work by Ribeiro and Soler, entitled 
"Learning Genetic and Environmental Graphical Models from Family Data" (submitted for publication). These tests are based on univariate polygenic linear mixed, with two components of variance: the polygenic or family-specific random effect, which models the phenotypic variability across the families, and the environmental or subject-specific error, which models phenotypic variability after removing the familial aggregation effect.

Particularly, for causal structure learning (learning of the structure of directed acyclic PGMs), these partial correlation tests are used as d-separation oracles in the IC/PC algorithm.

**Keywords:** Structure learning; Causal inference; Covariance matrix decomposition; Polygenic mixed model; Zero partial correlation test; Confounded residuals.

### Reference Manual

Documentation of the methods provided by FamilyBasedPGMs can be found at: https://github.com/adele/FamilyBasedPGMs/blob/master/manual.pdf

### Example

An example is provided at https://github.com/adele/FamilyBasedPGMs/blob/master/vignettes/familybasedpgms-example.pdf.


### Installation

First, install the following R packages:

```r
source("https://bioconductor.org/biocLite.R")

# The following packages must be installed together with RGBL
biocLite("graph")
biocLite("RBGL")
biocLite("Rgraphviz")

install.packages(c("stringi", "curl", "pcalg"), dependencies=TRUE)
```
Note: if you are asked to update packages, then press "a" for all.

Then, proceed with the installation of the FamilyBasedPGMs R package. 

You can download the latest tar.gz file with the source code of the FamilyBasedPGMs R package, available at https://github.com/adele/FamilyBasedPGMs/releases/latest, and install it with the following command, where `path_to_file` represents the full path and file name of the tar.gz file:

```r
install.packages(path_to_file, repos=NULL, type="source", dependencies=TRUE)
```

All releases are available at https://github.com/adele/FamilyBasedPGMs/releases. If you want a specific version of the FamilyBasedPGMs R package, for example, v1.0, you can install it directly from the URL:

```r
install.packages("https://github.com/adele/FamilyBasedPGMs/releases/download/v1.0/FamilyBasedPGMs_1.0.tar.gz", repos=NULL, method="libcurl", dependencies=TRUE)
```

Alternatively, you can install the development version directly from GitHub. Make sure you have the devtools R package installed. If not, install it with `install.packages("devtools", dependencies=TRUE)`.

```r
devtools::install_github("adele/FamilyBasedPGMs", dependencies="Import")
```

### Installation Trobleshooting on Debian and Ubuntu

If any of the following packages is missing, simply install the missing package(s) from the repository of your Linux distribution:

```console
sudo apt-get install libgmp10 libgmp-dev libatlas3-base libv8-dev libcurl4-openssl-dev libmagick++-dev
```

If you are getting an error when loading the shared library libgfortran.so.4, then you may need to install GCC 7 including the Fortran part, that includes libgfortran 4. See [here](https://stackoverflow.com/questions/46516394/how-to-install-libgfortran-so-4-on-ubuntu-16-06) for instructions.

```console
sudo add-apt-repository ppa:jonathonf/gcc-7.1
sudo apt-get update

sudo apt-get install gcc-7 g++-7 gfortran-7
```
