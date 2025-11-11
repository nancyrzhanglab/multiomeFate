# Purpose

`{r, include = FALSE} knitr::opts_chunk$set( collapse = TRUE, comment = "#>" )`

This repository contains all the functions to perform CYFER and the
downstream analysis, for the paper “Temporal and Clonal Resolution of
Cellular Evolution Under Stress”. See the companion GitHub package
<https://github.com/nancyrzhanglab/multiomeFate_analysis> for all the
analyses performed in the paper.

This code was developed and tested primarily on R 4.3.2. on a Macbook
(macOS Sonoma 14.2.1) equipped with an i7 processor.

# Installation

This package can be installed through `devtools` in R.

``` r
library("devtools")
devtools::install_github("nancyrzhanglab/multiomeFate")
```

The package itself depends on several packages. These include `Seurat`,
`cowplot`, `ggplot2`, `ggrepel`, `ggtern`, `scCustomize`, and `plyr`.
See the last section of this README to see where (i.e., CRAN,
Bioconductor, or GitHub) to download all such packages.

After installation of all the dependencies, the installation of the
`multiomeFate` package itself is fast (less than 2 minutes).

``` r
> devtools::session_info()
─ Session info ───────────────────
 setting  value
 version  R version 4.3.2 (2023-10-31)
 os       macOS Sonoma 14.2.1
 system   aarch64, darwin20
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       America/Los_Angeles
 date     2025-06-11
 rstudio  2024.12.0+467 Kousa Dogwood (desktop)
 pandoc   NA
 quarto   1.5.57 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto

─ Packages ───────────────────────
 !  package          * version    date (UTC) lib source
    abind              1.4-8      2024-09-12 [1] CRAN (R 4.3.3)
    bayesm             3.1-6      2023-09-23 [1] CRAN (R 4.3.1)
    beeswarm           0.4.0      2021-06-01 [1] CRAN (R 4.3.0)
    brio               1.1.5      2024-04-24 [1] CRAN (R 4.3.1)
    cachem             1.1.0      2024-05-16 [1] CRAN (R 4.3.3)
    circlize           0.4.16     2024-02-20 [1] CRAN (R 4.3.3)
    cli                3.6.5      2025-04-23 [1] CRAN (R 4.3.3)
    cluster            2.1.8.1    2025-03-12 [1] CRAN (R 4.3.3)
    codetools          0.2-20     2024-03-31 [1] CRAN (R 4.3.1)
    colorspace         2.1-1      2024-07-26 [1] CRAN (R 4.3.3)
    compositions       2.0-8      2024-01-31 [1] CRAN (R 4.3.1)
    cowplot            1.1.3      2024-01-22 [1] CRAN (R 4.3.1)
    data.table         1.17.4     2025-05-26 [1] CRAN (R 4.3.3)
    deldir             2.0-4      2024-02-28 [1] CRAN (R 4.3.1)
    DEoptimR           1.1-3-1    2024-11-23 [1] CRAN (R 4.3.3)
    desc               1.4.3      2023-12-10 [1] CRAN (R 4.3.1)
    devtools         * 2.4.5      2022-10-11 [1] CRAN (R 4.3.0)
    dichromat          2.0-0.1    2022-05-02 [1] CRAN (R 4.3.3)
    digest             0.6.37     2024-08-19 [1] CRAN (R 4.3.3)
    dotCall64          1.2        2024-10-04 [1] CRAN (R 4.3.3)
    dplyr              1.1.4      2023-11-17 [1] CRAN (R 4.3.1)
    ellipsis           0.3.2      2021-04-29 [1] CRAN (R 4.3.0)
    farver             2.1.2      2024-05-13 [1] CRAN (R 4.3.3)
    fastDummies        1.7.5      2025-01-20 [1] CRAN (R 4.3.3)
    fastmap            1.2.0      2024-05-15 [1] CRAN (R 4.3.3)
    fitdistrplus       1.2-2      2025-01-07 [1] CRAN (R 4.3.3)
    forcats            1.0.0      2023-01-29 [1] CRAN (R 4.3.0)
    fs                 1.6.6      2025-04-12 [1] CRAN (R 4.3.3)
    future             1.40.0     2025-04-10 [1] CRAN (R 4.3.3)
    future.apply       1.11.3     2024-10-27 [1] CRAN (R 4.3.3)
    generics           0.1.4      2025-05-09 [1] CRAN (R 4.3.3)
    ggbeeswarm         0.7.2      2023-04-29 [1] CRAN (R 4.3.0)
    ggplot2            3.5.2      2025-04-09 [1] CRAN (R 4.3.3)
    ggprism            1.0.6      2025-05-17 [1] CRAN (R 4.3.3)
    ggrastr            1.0.2      2023-06-01 [1] CRAN (R 4.3.0)
    ggrepel            0.9.6      2024-09-07 [1] CRAN (R 4.3.3)
    ggridges           0.5.6      2024-01-23 [1] CRAN (R 4.3.1)
    ggtern             3.5.0      2024-03-24 [1] CRAN (R 4.3.1)
    GlobalOptions      0.1.2      2020-06-10 [1] CRAN (R 4.3.0)
    globals            0.17.0     2025-04-16 [1] CRAN (R 4.3.3)
    glue               1.8.0      2024-09-30 [1] CRAN (R 4.3.3)
    goftest            1.2-3      2021-10-07 [1] CRAN (R 4.3.0)
    gridExtra          2.3        2017-09-09 [1] CRAN (R 4.3.0)
    gtable             0.3.6      2024-10-25 [1] CRAN (R 4.3.3)
    hexbin             1.28.5     2024-11-13 [1] CRAN (R 4.3.3)
    htmltools          0.5.8.1    2024-04-04 [1] CRAN (R 4.3.2)
    htmlwidgets        1.6.4      2023-12-06 [1] CRAN (R 4.3.1)
    httpuv             1.6.16     2025-04-16 [1] CRAN (R 4.3.3)
    httr               1.4.7      2023-08-15 [1] CRAN (R 4.3.0)
    ica                1.0-3      2022-07-08 [1] CRAN (R 4.3.0)
    igraph             2.1.4      2025-01-23 [1] CRAN (R 4.3.3)
    irlba              2.3.5.1    2022-10-03 [1] CRAN (R 4.3.2)
    janitor            2.2.1      2024-12-22 [1] CRAN (R 4.3.3)
    jsonlite           2.0.0      2025-03-27 [1] CRAN (R 4.3.3)
    KernSmooth         2.23-26    2025-01-01 [1] CRAN (R 4.3.3)
    later              1.4.2      2025-04-08 [1] CRAN (R 4.3.3)
    latex2exp          0.9.6      2022-11-28 [1] CRAN (R 4.3.0)
    lattice            0.22-7     2025-04-02 [1] CRAN (R 4.3.3)
    lazyeval           0.2.2      2019-03-15 [1] CRAN (R 4.3.0)
    lifecycle          1.0.4      2023-11-07 [1] CRAN (R 4.3.1)
    listenv            0.9.1      2024-01-29 [1] CRAN (R 4.3.1)
    lmtest             0.9-40     2022-03-21 [1] CRAN (R 4.3.0)
    lubridate          1.9.4      2024-12-08 [1] CRAN (R 4.3.3)
    magrittr           2.0.3      2022-03-30 [1] CRAN (R 4.3.0)
    MASS               7.3-60.0.1 2024-01-13 [1] CRAN (R 4.3.1)
    Matrix             1.6-5      2024-01-11 [1] CRAN (R 4.3.2)
    matrixStats        1.5.0      2025-01-07 [1] CRAN (R 4.3.3)
    memoise            2.0.1      2021-11-26 [1] CRAN (R 4.3.0)
    mime               0.13       2025-03-17 [1] CRAN (R 4.3.3)
    miniUI             0.1.2      2025-04-17 [1] CRAN (R 4.3.3)
 VP multiomeFate     * 1.0.0      2025-06-04 [?] load_all() (on disk 0.0.0.100)
    nlme               3.1-168    2025-03-31 [1] CRAN (R 4.3.3)
    paletteer          1.6.0      2024-01-21 [1] CRAN (R 4.3.1)
    parallelly         1.45.0     2025-06-02 [1] CRAN (R 4.3.3)
    patchwork          1.3.0      2024-09-16 [1] CRAN (R 4.3.3)
    pbapply            1.7-2      2023-06-27 [1] CRAN (R 4.3.0)
    pillar             1.10.2     2025-04-05 [1] CRAN (R 4.3.3)
    pkgbuild           1.4.8      2025-05-26 [1] CRAN (R 4.3.3)
    pkgconfig          2.0.3      2019-09-22 [1] CRAN (R 4.3.0)
    pkgload            1.4.0      2024-06-28 [1] CRAN (R 4.3.3)
    plotly             4.10.4     2024-01-13 [1] CRAN (R 4.3.1)
    plyr               1.8.9      2023-10-02 [1] CRAN (R 4.3.1)
    png                0.1-8      2022-11-29 [1] CRAN (R 4.3.0)
    polyclip           1.10-7     2024-07-23 [1] CRAN (R 4.3.3)
    profvis            0.4.0      2024-09-20 [1] CRAN (R 4.3.3)
    progressr          0.15.1     2024-11-22 [1] CRAN (R 4.3.3)
    promises           1.3.3      2025-05-29 [1] CRAN (R 4.3.3)
    proto              1.0.0      2016-10-29 [1] CRAN (R 4.3.0)
    purrr              1.0.4      2025-02-05 [1] CRAN (R 4.3.3)
    R6                 2.6.1      2025-02-15 [1] CRAN (R 4.3.3)
    RANN               2.6.2      2024-08-25 [1] CRAN (R 4.3.3)
    RColorBrewer       1.1-3      2022-04-03 [1] CRAN (R 4.3.0)
    Rcpp               1.0.14     2025-01-12 [1] CRAN (R 4.3.3)
    RcppAnnoy          0.0.22     2024-01-23 [1] CRAN (R 4.3.1)
    RcppHNSW           0.6.0      2024-02-04 [1] CRAN (R 4.3.1)
    rematch2           2.1.2      2020-05-01 [1] CRAN (R 4.3.0)
    remotes            2.5.0      2024-03-17 [1] CRAN (R 4.3.3)
    reshape2           1.4.4      2020-04-09 [1] CRAN (R 4.3.0)
    reticulate         1.42.0     2025-03-25 [1] CRAN (R 4.3.3)
    rlang              1.1.6      2025-04-11 [1] CRAN (R 4.3.3)
    robustbase         0.99-4-1   2024-09-27 [1] CRAN (R 4.3.3)
    ROCR               1.0-11     2020-05-02 [1] CRAN (R 4.3.0)
    rprojroot          2.0.4      2023-11-05 [1] CRAN (R 4.3.1)
    RSpectra           0.16-2     2024-07-18 [1] CRAN (R 4.3.3)
    rstudioapi         0.17.1     2024-10-22 [1] CRAN (R 4.3.3)
    Rtsne              0.17       2023-12-07 [1] CRAN (R 4.3.1)
    scales             1.4.0      2025-04-24 [1] CRAN (R 4.3.3)
    scattermore        1.2        2023-06-12 [1] CRAN (R 4.3.0)
    scCustomize        3.0.1      2024-12-18 [1] CRAN (R 4.3.3)
    sctransform        0.4.2      2025-04-30 [1] CRAN (R 4.3.3)
    sessioninfo        1.2.3      2025-02-05 [1] CRAN (R 4.3.3)
    Seurat             5.2.1      2025-01-24 [1] CRAN (R 4.3.3)
    SeuratObject       5.1.0      2025-04-22 [1] CRAN (R 4.3.3)
    shape              1.4.6.1    2024-02-23 [1] CRAN (R 4.3.1)
    shiny              1.10.0     2024-12-14 [1] CRAN (R 4.3.3)
    snakecase          0.11.1     2023-08-27 [1] CRAN (R 4.3.0)
    sp                 2.2-0      2025-02-01 [1] CRAN (R 4.3.3)
    spam               2.11-1     2025-01-20 [1] CRAN (R 4.3.3)
    spatstat.data      3.1-6      2025-03-17 [1] CRAN (R 4.3.3)
    spatstat.explore   3.4-3      2025-05-21 [1] CRAN (R 4.3.3)
    spatstat.geom      3.4-1      2025-05-20 [1] CRAN (R 4.3.3)
    spatstat.random    3.4-1      2025-05-20 [1] CRAN (R 4.3.3)
    spatstat.sparse    3.1-0      2024-06-21 [1] CRAN (R 4.3.3)
    spatstat.univar    3.1-3      2025-05-08 [1] CRAN (R 4.3.3)
    spatstat.utils     3.1-4      2025-05-15 [1] CRAN (R 4.3.3)
    stringi            1.8.7      2025-03-27 [1] CRAN (R 4.3.3)
    stringr            1.5.1      2023-11-14 [1] CRAN (R 4.3.1)
    survival           3.8-3      2024-12-17 [1] CRAN (R 4.3.3)
    tensor             1.5        2012-05-05 [1] CRAN (R 4.3.0)
    tensorA            0.36.2.1   2023-12-13 [1] CRAN (R 4.3.1)
    testthat         * 3.2.3      2025-01-13 [1] CRAN (R 4.3.3)
    tibble             3.3.0      2025-06-08 [1] CRAN (R 4.3.3)
    tidyr              1.3.1      2024-01-24 [1] CRAN (R 4.3.1)
    tidyselect         1.2.1      2024-03-11 [1] CRAN (R 4.3.1)
    timechange         0.3.0      2024-01-18 [1] CRAN (R 4.3.1)
    urlchecker         1.0.1      2021-11-30 [1] CRAN (R 4.3.0)
    usethis          * 3.1.0      2024-11-26 [1] CRAN (R 4.3.3)
    uwot               0.2.3      2025-02-24 [1] CRAN (R 4.3.3)
    vctrs              0.6.5      2023-12-01 [1] CRAN (R 4.3.1)
    vipor              0.4.7      2023-12-18 [1] CRAN (R 4.3.1)
    viridisLite        0.4.2      2023-05-02 [1] CRAN (R 4.3.0)
    withr              3.0.2      2024-10-28 [1] CRAN (R 4.3.3)
    xtable             1.8-4      2019-04-21 [1] CRAN (R 4.3.0)
    zoo                1.8-14     2025-04-10 [1] CRAN (R 4.3.3)

 [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library

 * ── Packages attached to the search path.
 V ── Loaded and on-disk version mismatch.
 P ── Loaded and on-disk path mismatch.

──────────────────────────────────
```
