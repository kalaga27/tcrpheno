# tcrpheno

`tcrpheno` is an R package that applies a logistic regression model to the amino acid sequences of T-cell receptor (TCR) chains to predict T cell phenotypes. The model generates scores that correlate with distinct T cell fates, offering insights into the functional potential of T cells based on their receptor sequences.
Citation

If you use `tcrpheno` in your research, please cite the publication:

> Lagattuta, K. et al. The T cell receptor sequence influences the likelihood of 
T cell memory formation. Cell Reports. 2025 Jan 28;44(1)


## Installation

You can install tcrpheno from Github with:

``` r
remotes::install_github("kalaga27/tcrpheno")
```

## The `score_tcrs()` Function

This is the core function of the package. It takes a data frame of TCR 
sequences and calculates the four phenotype scores.

```r
score_tcrs(data, chain)
```

### Arguments:
* `data`: A data frame containing the TCR information. The data frame must 
include a cell identifier in the first column. The order of the remaining 
columns does not matter, but specific column names are required. To see the 
expected format, load the example dataset with ```data("tcrpheno_clones")```.
  * For paired-chain analysis, the following columns are required: `TCRA_cdr3aa`, 
  `TCRA_vgene`, `TCRA_jgene`, `TCRB_cdr3aa`, `TCRB_vgene`, `TCRB_jgene`.
  * For single-chain analysis, only the columns corresponding to that chain 
  are needed.
* `chain`: A character string specifying the TCR chains available in the data.
  * "ab": For paired alpha and beta TCR data.
  * "a": For alpha chain data only.
  * "b": For beta chain data only.
  
## Output Scores

The `score_tcrs()` function returns a data frame with four columns, each 
corresponding to a predicted phenotype score:
* TCRinnate: Higher scores suggest a greater likelihood of the T cell adopting 
an innate-like, PLZF-high phenotype (e.g., MAIT and iNKT cells).
* TCR.CD8: Higher scores indicate a predisposition towards a CD8+ T cell fate 
over a CD4+ fate.
* TCRreg: Higher scores point to an increased probability of the T cell 
becoming a regulatory T cell (Treg).
* TCRmem: Higher scores suggest a T cell is more likely to differentiate into 
a memory cell rather than remaining naive.

# Reporting Issues

We welcome community feedback and contributions!

## Reporting Issues

If you encounter a bug, have a question, or want to suggest an enhancement,
please open an issue on the [GitHub repository](https://github.com/kalaga27/tcrpheno). 
When reporting a bug, please include a minimal reproducible example so we can quickly identify 
and resolve the problem.
