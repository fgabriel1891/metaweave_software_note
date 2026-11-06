---

editor_options: 
  markdown: 
    wrap: 72
---

# MetaWeave Ecography supplementary analyses

Render each QMD from this directory, for example:

``` sh
quarto render S1-synthetic-benchmark.qmd --to docx
quarto render S2-palm-mammal.qmd --to docx
quarto render S3-ecuador-assembly.qmd --to docx
```

The `data/` folder is part of the supplement. Keep it beside the QMDs. All three use MetaWeave 1.0.0 or later and `tidyverse`, including `ggplot2`. S1 and S2 also use `terra`; S2 uses `patchwork` for the map layout. S3 uses `vegan` and `bipartite` when running the displayed ensemble code. The displayed SBM fitting chunks require `cassandRa` only if run manually. Quarto, `knitr` and `rmarkdown` are needed to render the documents.

Both HTML and Word are configured in each QMD. Use `--to html` for an HTML preview. Code is visible in both formats. The fitting and benchmark generation chunks marked `eval: false` show how archived quantities were constructed; they are displayed without being executed. S3 loads saved simulation draws by default. These choices keep the precomputed steps visible while avoiding repeated fits during rendering.

| Document | Question | Included input |
|------------------------|------------------------|------------------------|
| S1 | How do block and species-pair rules recover known local probabilities? | Recorded simulation seed, model parameters, known-process probabilities; ranges are regenerated in the QMD |
| S2 | Where do regional records and an SBM make different spatial predictions? | Regional interaction CSV, matched palm and mammal range rasters, regional matrix, SBM parameters |
| S3 | What do regional interaction structure and species composition contribute to local network metrics? | Ecuador interaction records, pooled SBM fit, 49-draw results |

S1's `scripts/synthetic_benchmark.R` regenerates the ranges and the complete benchmark from the original seed. S1 displays a `source()` call for that script, while its default render reads the saved results. S2 starts from the matched range rasters used by the analysis; it does not rebuild them from their original range sources. S3 starts from a saved SBM fitted to the six-site pooled metaweb; its QMD contains the complete local reconstruction and simulation code. The fit was produced with `cassandRa` 0.2.0, and the source fit script is retained at `scripts/fit_ecuador_sbm.R`. Run it with Rscript to write a new fit to `supplementary-material/results/` without replacing the included fit. The original Ecuador data are at <https://zenodo.org/records/5537283>.

The main computational chunks are displayed with `eval: false`. Hidden chunks load saved results for rendering. The displayed code can be run interactively to repeat the analysis; the S3 simulations take much longer than rendering the saved results.
