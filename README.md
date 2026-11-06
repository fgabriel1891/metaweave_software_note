# MetaWeave: Ecography Software Note

This repository contains the manuscript source and supporting analyses for
the MetaWeave Software Note. The article uses MetaWeave 1.0.0.

## Contents

- `article.qmd` and `references.bib`: manuscript and bibliography.
- `images/`: the workflow figure, its DiagrammeR source, and the short
  in-paper example figure.
- `supplementary-material/`: S1 synthetic benchmark, S2 palm–mammal case,
  S3 Ecuador assembly case, and the inputs and saved outputs for each.
- `supplementary-material/scripts/synthetic_benchmark.R`: full benchmark
  simulation, inference and scoring script, sourced in S1.
- `supplementary-material/scripts/fit_ecuador_sbm.R`: script used to fit the
  Ecuador regional model from the interaction data included with S3.
- `images/` also contains the manuscript's benchmark and case-study figures.

## Render

Install Quarto, R, `knitr` and `rmarkdown`. From the repository root:

```sh
quarto render article.qmd --to docx
```

The article uses existing figures and does not execute its R code during
rendering. To render the supplements, work in `supplementary-material/`:

```sh
quarto render S1-synthetic-benchmark.qmd --to docx
quarto render S2-palm-mammal.qmd --to docx
quarto render S3-ecuador-assembly.qmd --to docx
```

The supplements require MetaWeave 1.0.0 or later and `tidyverse`; S1 and S2
also need `terra`, and S2 needs `patchwork`. Fitting the displayed SBM
examples requires `cassandRa`; repeating S3 simulations requires `vegan`
and `bipartite`. Main simulations and fits are shown without executing when
rendering, and supplied results are loaded for figures and tables. See the
[supplement README](supplementary-material/README.md) for details.

To regenerate the workflow SVG from the project root, install `DiagrammeR`
and `DiagrammeRsvg`, then run `source("images/create-workflow.R")` in R.

## Sources 

The Ecuador data are archived at <https://zenodo.org/records/5537283>.
The palm interaction records are from <https://doi.org/10.5061/dryad.rd46vq3>.
The package and its license are available at
<https://github.com/fgabriel1891/metaweave>. 


