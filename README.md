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

## load data
data(anzia2012)
```

``` r
## estimate treatment effect
set.seed(1234)
fit_panel <- did(
  formula  = lnavgsalary_cpi ~ oncycle | teachers_avg_yrs_exper + 
                                          ami_pc + asian_pc + black_pc + hisp_pc,
  data     = anzia2012,
  id_unit  = "district",
  id_time  = "year",
  design   = "did",
  is_panel = TRUE,
  option   = list(n_boot = 100, parallel = TRUE, lead = 0:2)
)
```

`did()` function takes the following arguments:

-   `formula`: A formula specifying variables. It should follow the form
    of `outcome ~ treatment | covariates`.
    -   `treatment` should be time-varying, that is, `treatment` takes
        zero for everyone before the treatment assignment, and takes 1
        for units who are treated. See the example for how the treatment
        variable should be coded.
    -   `covariates` can be omitted as `outcome ~ treatment`.
-   `data`: A data frame. This can be either `data.frame` or `tibble`/
-   `id_unit`: A variable name in the data that uniquely identifies
    units (e.g., individuals or states).
-   `id_time`: A variable name in the data that uniquely identifies time
    (e.g., year).
-   `design`: Design option. It should be `"did"` when the standard DID
    design is used.
-   `is_panel`: A boolean argument to indicate the type of the data.
    When the dataset is panel (i.e., same observations are measured
    repeately overtime), it should take `TRUE`. See the next section for
    how to analyze the repeated cross-section data.
-   `option`: A list of options.
    -   `n_boot`: Number of bootstrap iterations to estimate weighting
        matrix.
    -   `parallel`: A boolean argument. If `TRUE`, bootstrap is
        conducted in parallel using `future` package.
    -   `lead`: A vector of non-negative lead parameter. For example,
        when `lead = c(0, 1)`, treatment effect when the treatment is
        assigned (`lead = 0`) as well as one-time ahead effect
        (`lead = 1`) will be estimated. Default is `lead = 0`.

``` r
## view the estimates
summary(fit_panel)
```

    ## # A tibble: 9 x 7
    ##   estimator   lead estimate std.error statistic p_value ddid_weight
    ##   <chr>      <int>    <dbl>     <dbl>     <dbl>   <dbl>       <dbl>
    ## 1 Double-DID     0 -0.00388   0.00264    -1.47  0.142        NA    
    ## 2 DID            0 -0.00620   0.00275    -2.26  0.0241       -0.286
    ## 3 sDID           0 -0.00439   0.00435    -1.01  0.312         1.29 
    ## 4 Double-DID     1 -0.00735   0.00346    -2.13  0.0336       NA    
    ## 5 DID            1 -0.0115    0.00392    -2.94  0.00332       0.504
    ## 6 sDID           1 -0.00313   0.00391    -0.800 0.424         0.496
    ## 7 Double-DID     2 -0.00648   0.00422    -1.53  0.125        NA    
    ## 8 DID            2 -0.0115    0.00505    -2.27  0.0230        0.613
    ## 9 sDID           2  0.00145   0.00457     0.318 0.750         0.387

`summary()` function can be used to view estimates.

### Repeated Cross-sectional Data

``` r
data(malesky2014)

set.seed(1234)
ff_rcs <- did(
 formula = transport ~ treatment + post_treat | factor(city),
 data    = malesky2014,
 id_time = 'year',
 is_panel= FALSE,
 option  = list(n_boot = 10, id_cluster = "tinh", parallel = TRUE)
)

summary(ff_rcs)
```

Staggered Adoption Design
-------------------------
