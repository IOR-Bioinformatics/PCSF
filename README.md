

PCSF: an R-Package for Network-Based Interpretation of High-throughput Data
===============================

The PCSF package performs fast and user-friendly network analysis of high-throughput data. Using interaction networks as a template, it determines high-confidence subnetworks relevant to the data, which potentially leads to predictions of functional units. It also interactively visualizes the resulting subnetwork with functional enrichment analysis.

Contact: Murodzhon Akhmedov [murodzhon.akhmedov@irb.usi.ch]



Reference:
--------------------
[A divide and conquer matheuristic algorithm for the Prize-collecting Steiner Tree Problem.](http://www.sciencedirect.com/science/article/pii/S0305054815003019)
Akhmedov M, Kwee I, and Montemanni R (2016). *Computers and Operations Research*, 70, 18-25.

A fast Prize-collecting Steiner Forest algorithm for Functional Analyses in Biological Networks.
Akhmedov M, LeNail A, Bertoni F, Kwee I, Fraenkel E and Montemanni R (2017). *Lecture Notes in Computer Science*, to appear.


System Requirements:
--------------------
1. R (>= 3.1.0)

2. Boost C++ library: http://www.boost.org



Installation:
--------------------

1. The PCSF package depends on the following R-packages: 

 - `BH` and `igraph` - for efficient graph handling and calculations,
 - `httr`, `methods`, `org.Hs.eg.db`, and `topGO` - to perform enrichment analysis,
 - `Rcpp`  - to employ C++ source code within R,
 - `visNetwork` - for visualization.


2. In order to compile the source, Windows users should install the `Rtools` package by the following [link](https://cran.r-project.org/bin/windows/Rtools/) that installs GCC and CMake.


3. The PCSF package and its dependencies can be installed on Mac OS, Linux and Windows by running the following commands in the R console.

```
source("http://bioconductor.org/biocLite.R")
biocLite("topGO")
install.packages("devtools", dependencies=TRUE)
devtools::install_github("IOR-Bioinformatics/PCSF", repos=BiocInstaller::biocinstallRepos(),
                         dependencies=TRUE, type="source", force=TRUE)
```


Comments:
--------------------

#### Test environments

* Mac OS X (10.12.4), R 3.4.0
* Ubuntu (16.04), R 3.2.3
* Windows 7, R 3.4.1

#### R CMD check results

There were no ERRORs, WARNINGs or NOTEs. 
