DIDdesign: Double Difference-in-Differences
===========================================

Authors:

-   [Naoki Egami](https://scholar.princeton.edu/negami/home)
-   Soichiro Yamauchi

Reference

Installation
------------

-   Downloading `DIDdesign` from Github

    ``` r
    ## need to install `devtools` if necessary
    require(devtools)
    install_github("naoki-egami/DIDdesign", dependencies = TRUE)
    ```

Example
-------

### Load Sample Data

``` r
# load package 
require(DIDdesign)
```

    ## Loading required package: DIDdesign

``` r
# load sample data 
data("anzia2012")
```

### Visualizing Data

### Nonparametric Estimator

``` r
# estimate ATT 
fit_np <- did(lnavgsalary_cpi ~ oncycle, data = anzia2012,
    id_subject = "district", id_time = "year",
    post_treatment = c(2007, 2008, 2009),
    method = "nonparametric",
    se_boot = TRUE, n_boot  = 100, boot_min = FALSE,
    select = "HQIC")
```

    ## 
    ## ... estimating treatment effect for  2007  ...
    ## ... bootstraping to compute standard errors ...
    ## ... computing the standard DiD estimate ...
    ## 
    ## ... estimating treatment effect for  2008  ...
    ## ... bootstraping to compute standard errors ...
    ## ... computing the standard DiD estimate ...
    ## 
    ## ... estimating treatment effect for  2009  ...
    ## ... bootstraping to compute standard errors ...
    ## ... computing the standard DiD estimate ...

``` r
# summary function 
summary(fit_np)
```

    ##                                2007             2008             2009
    ## 1   D-DiD      ATT           -0.007           -0.012           -0.012
    ## 2           95% CI [-0.013, -0.002]  [-0.02, -0.003]  [-0.021, 0.002]
    ## 3              BIC          -18.018          -18.018          -18.018
    ## 4             HQIC           -9.147           -9.147           -9.147
    ## 5         Selected               M1               M1               M1
    ## 6 Std-DiD      ATT           -0.007           -0.011           -0.011
    ## 7           95% CI [-0.013, -0.001] [-0.019, -0.003] [-0.021, -0.003]

``` r
# plot results 
plot(fit_np, full = TRUE)
```

![](README_files/figure-markdown_github/unnamed-chunk-2-1.png)

### Parametric Estimator (without Covariates)

``` r
fit_p1 <- did(lnavgsalary_cpi ~ oncycle, data = anzia2012,
                            id_subject = "district", id_time = "year",
                            post_treatment = c(2007, 2008, 2009),
                            method = "parametric",
                            se_boot = FALSE, select = "HQIC")
```

    ## 
    ## ... estimating treatment effect for  2007  ...
    ## 
    ## ... estimating treatment effect for  2008  ...
    ## 
    ## ... estimating treatment effect for  2009  ...

``` r
# summary function 
summary(fit_p1)
```

    ##                                2007             2008             2009
    ## 1   D-DiD      ATT           -0.007           -0.012           -0.012
    ## 2           95% CI [-0.008, -0.006] [-0.013, -0.011] [-0.013, -0.011]
    ## 3              BIC          -18.018          -18.018          -18.018
    ## 4             HQIC           -9.147           -9.147           -9.147
    ## 5         Selected               M1               M1               M1
    ## 6 Std-DiD      ATT           -0.007           -0.011           -0.011
    ## 7           95% CI [-0.013, -0.001] [-0.019, -0.003] [-0.022, -0.001]

``` r
# plot function 
plot(fit_p1, full = TRUE, ylim = c(-0.025, 0.02))
```

![](README_files/figure-markdown_github/unnamed-chunk-4-1.png)
