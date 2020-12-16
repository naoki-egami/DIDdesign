# DIDdesign: Double Difference-in-Differences

<!-- badges: start -->

[![R build
status](https://github.com/naoki-egami/DIDdesign/workflows/R-CMD-check/badge.svg)](https://github.com/naoki-egami/DIDdesign/actions)
<!-- badges: end -->

Authors:

  - [Naoki Egami](https://naokiegami.com/)
  - [Soichiro Yamauchi](https://soichiroy.github.io/)

Reference:

  - Egami and Yamauchi (2019) “Using Multiple Pre-treatment Periods to
    Improve Difference-in-Differences and Staggered Adoption Designs.”

## Installation Instructions:

  - Downloading the most recent version of `DIDdesign` from Github
    
    ``` r
    ## need to install `devtools` if necessary
    require(devtools)
    install_github("naoki-egami/DIDdesign", dependencies = TRUE)
    ```

## Table of Contents

1.  Basic DID design with Panel Data
2.  Basic DID design with Repeated Cross-Section Data
3.  Staggered adoption design

## The Basic Difference-in-Differences Design with Panel Data

``` r
## load package
require(DIDdesign)

## load data
data(anzia2012)
```

### Step 1: Assess the parallel trends assumption

``` r
## check parallel trends
set.seed(1234)
check_panel <- did_check(
  formula = lnavgsalary_cpi ~ oncycle | teachers_avg_yrs_exper +
                       ami_pc + asian_pc + black_pc + hisp_pc,
  data    = anzia2012,
  id_unit = "district",
  id_time = "year",
  option  = list(n_boot = 200, parallel = TRUE, lag = 1:3)
)

check_panel$estimate
#> # A tibble: 3 x 5
#>    estimate   lag std.error EqCI95_LB EqCI95_UB
#>       <dbl> <int>     <dbl>     <dbl>     <dbl>
#> 1 -0.00361      1   0.00265  -0.00796   0.00796
#> 2  0.00326      2   0.00231  -0.00705   0.00705
#> 3 -0.000434     3   0.00271  -0.00489   0.00489
```

``` r
plot(check_panel)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" style="display: block; margin: auto;" />

If a user wish to create custom figures, data used to generate the above
plot are available via `check_panel$plot[[1]]$dat_plot` and
`check_panel$plot[[2]]$dat_plot`.

### Step 2: Estimate the treatment effect with the double DID estimator

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
  option   = list(n_boot = 200, parallel = TRUE, lead = 0:2)
)
```

`did()` function takes the following arguments:

  - `formula`: A formula specifying variables. It should follow the form
    of `outcome ~ treatment | covariates`.
      - `treatment` should be time-varying, that is, `treatment` takes
        zero for everyone before the treatment assignment, and takes 1
        for units who are treated. See the example for how the treatment
        variable should be coded.
      - `covariates` can be omitted as `outcome ~ treatment`.
  - `data`: A data frame. This can be either `data.frame` or `tibble`.
  - `id_unit`: A variable name in the data that uniquely identifies
    units (e.g., individuals or states).
  - `id_time`: A variable name in the data that uniquely identifies time
    (e.g., year).
  - `design`: Design option. It should be `"did"` when the standard DID
    design is used.
  - `is_panel`: A boolean argument to indicate the type of the data.
    When the dataset is panel (i.e., same observations are measured
    repeately overtime), it should take `TRUE`. See the next section for
    how to analyze the repeated cross-section data.
  - `option`: A list of options.
      - `n_boot`: Number of bootstrap iterations to estimate weighting
        matrix.
      - `parallel`: A boolean argument. If `TRUE`, bootstrap is
        conducted in parallel using `future` package.
      - `lead`: A vector of non-negative lead parameter. For example,
        when `lead = c(0, 1)`, treatment effect when the treatment is
        assigned (`lead = 0`) as well as one-time ahead effect (`lead
        = 1`) will be estimated. Default is `lead = 0`.

<!-- end list -->

``` r
## view the estimates
summary(fit_panel)
#> # A tibble: 9 x 6
#>   estimator   lead estimate std.error statistic p_value
#>   <chr>      <int>    <dbl>     <dbl>     <dbl>   <dbl>
#> 1 Double-DID     0 -0.00412   0.00261    -1.58  0.115  
#> 2 DID            0 -0.00620   0.00265    -2.34  0.0194 
#> 3 sDID           0 -0.00439   0.00440    -0.999 0.318  
#> 4 Double-DID     1 -0.00674   0.00323    -2.09  0.0370 
#> 5 DID            1 -0.0115    0.00365    -3.16  0.00160
#> 6 sDID           1 -0.00313   0.00393    -0.796 0.426  
#> 7 Double-DID     2 -0.00508   0.00430    -1.18  0.238  
#> 8 DID            2 -0.0115    0.00492    -2.33  0.0196 
#> 9 sDID           2  0.00145   0.00490     0.297 0.766
```

`summary()` function can be used to view estimates.

``` r
require(patchwork)
plot(fit_panel) + plot(fit_panel, check_panel)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" style="display: block; margin: auto;" />

## The Standard DID Design with Repeated Cross-sectional Data

Sometimes, each period consists of different units, instead of repeated
observations of the same units. `did()` can handle such “repeated
cross-sectional” data by setting `is_panel = FALSE`. As an example, we
analyze `malesky2014` dataset (see `?malesky2014` for more details on
this dataset).

``` r
## load data
data(malesky2014)

## check parallel trends
set.seed(1234)
check_rcs <- did_check(
  formula = transport ~ treatment + post_treat | factor(city),
  data    = malesky2014,
  id_time = 'year',
  is_panel= FALSE,
  option  = list(n_boot = 200, parallel = TRUE, id_cluster = "tinh", lag = 1)
)

## summary
check_rcs$estimate
#> # A tibble: 1 x 5
#>   estimate   lag std.error EqCI95_LB EqCI95_UB
#>      <dbl> <dbl>     <dbl>     <dbl>     <dbl>
#> 1  -0.0609     1    0.0485    -0.141     0.141
```

``` r
plot(check_rcs)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" style="display: block; margin: auto;" />

``` r
## estimate ATT
ff_rcs <- did(
  formula = transport ~ treatment + post_treat | factor(city),
  data    = malesky2014,
  id_time = 'year',
  is_panel= FALSE,
  option  = list(n_boot = 200, parallel = TRUE, id_cluster = "tinh", lead = 0)
)
```

  - `formula` now includes `treatment` variable as well as the
    post-treatment time indicator.
  - `id_cluster`: A variable used to cluster the standard errors.

<!-- end list -->

``` r
summary(ff_rcs)
#> # A tibble: 3 x 6
#>   estimator   lead estimate std.error statistic p_value
#>   <chr>      <dbl>    <dbl>     <dbl>     <dbl>   <dbl>
#> 1 Double-DID     0    0.239    0.0909     2.63  0.00842
#> 2 DID            0    0.101    0.105      0.957 0.338  
#> 3 sDID           0    0.169    0.138      1.22  0.221
```

## Staggered Adoption Design

`DIDdesign` supports the staggered adoption design where units receive
the treatment at different periods of time. As an example, we analyze
`paglayan2019` dataset in the package (see `?paglayan2019` for more
details about this dataset).

``` r
## data
require(dplyr)
require(tibble)

## format dataset
paglayan2019 <- paglayan2019 %>%
  filter(!(state %in% c("WI", "DC"))) %>%
  mutate(id_time = year,
         id_subject = as.numeric(as.factor(state)),
         log_expenditure = log(pupil_expenditure + 1),
         log_salary      = log(teacher_salary + 1))
```

<img src="man/figures/README-sa_plot-1.png" width="100%" />

As we can see in the above plot, states receive the treatment at
different years ranging from 1965 at earliest to 1987 at latest (and
some of the states never receive the treatment).

## Step 1: Assess the pre-treatment parallel trends

``` r
set.seed(1234)
check_sa <- did_check(
  formula = log_expenditure ~ treatment,
  data    = paglayan2019,
  id_unit = "id_subject",
  id_time = "id_time",
  design  = "sa",
  option  = list(n_boot = 200, parallel = TRUE, thres = 1, lag = 1:5)
)

## view estimates
check_sa$estimate
#> # A tibble: 5 x 5
#>   estimate   lag std.error EqCI95_LB EqCI95_UB
#>      <dbl> <int>     <dbl>     <dbl>     <dbl>
#> 1 -0.00267     1   0.00864   -0.0169    0.0169
#> 2 -0.0124      2   0.00886   -0.0270    0.0270
#> 3  0.00227     3   0.0110    -0.0204    0.0204
#> 4 -0.00758     4   0.0119    -0.0271    0.0271
#> 5 -0.0107      5   0.00894   -0.0254    0.0254
```

``` r
plot(check_sa)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" style="display: block; margin: auto;" />

## Step 2: Estimate staggered-adoption average treatment effect

`did()` function can handle the staggered adoption design by setting the
`design` argument to `design = "sa"`.

``` r
## estimate time-weighted SA-ATE
set.seed(1234)
fit_sa <- did(
  formula = log_expenditure ~ treatment,
  data    = paglayan2019,
  id_unit = "id_subject",
  id_time = "id_time",
  design  = "sa",
  option  = list(n_boot = 200, lead = 0:9, thres = 2, parallel = TRUE)
)
```

In addition to options described in the previous section, there is one
additional argument specific to the staggered adoption design.

  - `thres` parameter in the option control the minimum number of
    treated units for a particular time to be included in the treatment
    effect estimation. For example if `thres = 2`, the effect for
    Tennessee will be removed from the time-average effect because it’s
    the only unit who received the treatment in 1972 (i.e., the number
    of treated units in 1972 is less than the threshold).

<!-- end list -->

``` r
head(summary(fit_sa))
#> # A tibble: 6 x 6
#>   estimator      lead  estimate std.error statistic p_value
#>   <chr>         <int>     <dbl>     <dbl>     <dbl>   <dbl>
#> 1 SA-Double-DID     0  0.00869     0.0206    0.422    0.673
#> 2 SA-DID            0  0.0149      0.0162    0.921    0.357
#> 3 SA-sDID           0  0.00903     0.0195    0.463    0.643
#> 4 SA-Double-DID     1 -0.0152      0.0163   -0.936    0.349
#> 5 SA-DID            1  0.000916    0.0151    0.0606   0.952
#> 6 SA-sDID           1 -0.0199      0.0206   -0.970    0.332
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="80%" />
