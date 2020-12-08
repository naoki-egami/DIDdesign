DIDdesign: Double Difference-in-Differences
===========================================

<!-- badges: start -->

[![R build
status](https://github.com/naoki-egami/DIDdesign/workflows/R-CMD-check/badge.svg)](https://github.com/naoki-egami/DIDdesign/actions)
<!-- badges: end -->

Authors:

-   [Naoki Egami](https://naokiegami.com/)
-   [Soichiro Yamauchi](https://soichiroy.github.io/)

Reference:

-   Egami and Yamauchi (2019) “How to Improve the
    Difference-in-Differences Design with Multiple Pre-treatment
    Periods.”

Installation Instructions:
--------------------------

-   Downloading the most recent version of `DIDdesign` from Github

    ``` r
    ## need to install `devtools` if necessary
    require(devtools)
    install_github("naoki-egami/DIDdesign", dependencies = TRUE)
    ```

The Standard Difference-in-Differences Design
---------------------------------------------

### Panel Data

``` r
## load package 
require(DIDdesign)
```

    ## Loading required package: DIDdesign

``` r
## load data 
data(anzia2012)


## estimate treatment effect
set.seed(1234)
fit_panel <- did(
  formula = lnavgsalary_cpi ~ oncycle | teachers_avg_yrs_exper +
                        ami_pc + asian_pc + black_pc + hisp_pc,
  data    = anzia2012,
  id_unit = "district",
  id_time = "year",
  option  = list(n_boot = 20, parallel = FALSE)
)

fit_panel$estimates
```

    ## # A tibble: 3 x 3
    ##   estimator  estimate std.error
    ##   <chr>         <dbl>     <dbl>
    ## 1 Double-DID -0.00409   0.00277
    ## 2 DID        -0.00620   0.00280
    ## 3 sDID       -0.00439   0.00408

``` r
## treatment effect on Y_{t+1}
fit_panel_t1 <- did(
  formula = lnavgsalary_cpi ~ oncycle | teachers_avg_yrs_exper +
                        ami_pc + asian_pc + black_pc + hisp_pc,
  data    = anzia2012,
  id_unit = "district",
  id_time = "year",
  option  = list(n_boot = 20, lead = 1, parallel = FALSE)
)

fit_panel_t1$estimates
```

    ## # A tibble: 3 x 3
    ##   estimator  estimate std.error
    ##   <chr>         <dbl>     <dbl>
    ## 1 Double-DID -0.00512   0.00232
    ## 2 DID        -0.00620   0.00255
    ## 3 sDID       -0.00439   0.00281

### Repeated Cross-sectional Data

Staggered Adoption Design
-------------------------
