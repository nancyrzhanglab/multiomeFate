# Plot a simplex (ternary) diagram

Creates a ternary plot for three-way compositional data using ggplot2,
without requiring the ggtern package.

## Usage

``` r
plot_simplex(
  df,
  x_col,
  y_col,
  z_col,
  color_col = NULL,
  size_col = NULL,
  col_palette = NULL,
  xlab = "x",
  ylab = "y",
  zlab = "z",
  title = "title"
)
```

## Arguments

- df:

  A data frame where each row is a point to plot.

- x_col:

  Character name of the column in `df` for the first axis (bottom-left
  corner).

- y_col:

  Character name of the column in `df` for the second axis (bottom-right
  corner).

- z_col:

  Character name of the column in `df` for the third axis (top corner).

- color_col:

  Character name of the column in `df` for point colors. Default `NULL`.

- size_col:

  Character name of the column in `df` for point sizes. Default `NULL`.

- col_palette:

  Named vector of colors passed to `scale_color_manual`. Default `NULL`.

- xlab:

  Label for the first axis (bottom-left corner). Default `"x"`.

- ylab:

  Label for the second axis (bottom-right corner). Default `"y"`.

- zlab:

  Label for the third axis (top corner). Default `"z"`.

- title:

  Plot title. Default `"title"`.

## Value

A ggplot2 object.
