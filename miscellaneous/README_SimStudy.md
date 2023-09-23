# Vignette for codes used for the simulation study Alam and Staicu, 20xx

This file illustrates how to execute the codes in files
*utilities_generate.R*, *utilities_main.R*, and *utilities_estimate.R*
to reproduce the simulation study.

- First, load all the functions in files *utilities_generate.R* and
  *utilities_main.R*

``` r
source("utilities_generate.R")
source("utilities_main.R")
```

- Note that, all the functions in files *utilities_generate.R* and
  *utilities_main.R* are available in *sLFDA* package. It can be
  installed from github via

``` r
devtools::install_github("https://github.com/msalam14/sLFDA")
```

- All the simulation results, except the quantile trajectory estimation
  via the proposed *sLFDA* approach, are obtained by excuting the codes
  written in *utilities_estimate.R* script.

``` r
source("utilities_estimate.R")
```

<div>

> **Note**
>
> We perform the simulation on HPC; particularly, parallel jobs were run
> on distributed memory. To do so, we use the *Rmpi* package in the
> script *utilities_estimate.R*. For a given sample size and skewness
> type, this script will save $4$ files; in total $36$ files will be
> stored for $9$ scenarios ($3$ sample sizes and $3$ skewness scenarios)

</div>

- The *R* script named as *utilities_results.R* contains all the
  necessary codes to summarize simulation results stored in the $36$
  files by the script *utilities_estimate.R*. This script also generate
  the Tables provided in the original manuscript and supplementary
  material.
