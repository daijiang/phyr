# Deprecated functions in phyr

These functions still work but will be removed in a future version of
phyr; use their replacements:

- `communityPGLMM()`:
  [`pglmm()`](https://daijiang.github.io/phyr/reference/pglmm.md)

- `communityPGLMM.show.re()` and `communityPGLMM.plot.re()`:
  [`pglmm_plot_ranef()`](https://daijiang.github.io/phyr/reference/pglmm-plot-re.md)
  (or
  [`pglmm_plot_re()`](https://daijiang.github.io/phyr/reference/pglmm-plot-re.md))

The functions `communityPGLMM.profile.LRT()`,
`communityPGLMM.matrix.structure()` and
`communityPGLMM.predicted.values()`, deprecated since phyr 1.0.3, were
removed in phyr 1.9.0: use
[`pglmm_profile_LRT()`](https://daijiang.github.io/phyr/reference/pglmm-profile-LRT.md),
[`pglmm_matrix_structure()`](https://daijiang.github.io/phyr/reference/pglmm-matrix-structure.md)
and
[`pglmm_predicted_values()`](https://daijiang.github.io/phyr/reference/pglmm-predicted-values.md).

## Usage

``` r
communityPGLMM(...)

communityPGLMM.show.re(...)

communityPGLMM.plot.re(...)
```

## Arguments

- ...:

  Arguments passed to the replacement function.

## Value

The value of the replacement function.
