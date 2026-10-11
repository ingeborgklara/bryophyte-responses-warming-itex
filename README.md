# Bryophyte community responses to long-term warming manipulations across contrasting tundra habitats

Data and analysis code for:

> Jónsdóttir, I. S., Klarenberg, I. & Russi Colmenares, A. J. Bryophyte community
> responses to long-term warming manipulations vary greatly between contrasting
> tundra habitats. Preprint: [10.32942/X2XM4Q](https://doi.org/10.32942/X2XM4Q)

The study compares bryophyte community structure, taxonomic composition,
diversity, bryophyte functional group (BFG) composition and selected functional
traits across five tundra habitats spanning three Arctic bioclimatic sub-zones,
and how each responded to decadal experimental warming with open top chambers
following [ITEX](https://www.gvsu.edu/itex/) protocols.

| Code | Habitat | Region | Established |
|------|---------|--------|-------------|
| BH | *Betula nana* heath, Audkúluheidi | Iceland (oro-Arctic) | 1996 |
| MH | *Racomitrium lanuginosum* moss heath, Thingvellir | Iceland (sub-Arctic) | 1995 |
| DH | *Dryas octopetala* heath, Endalen | Svalbard (high Arctic) | 2002 |
| CH | *Cassiope tetragona* heath, Endalen | Svalbard (high Arctic) | 2002 |
| SB | Snowbed, Endalen | Svalbard (high Arctic) | 2002 |

## Contents

```
├── bryophytes_itex_analysis.R      all analyses, run top to bottom
├── data/raw/
│   └── bryophytes_itex_data.xlsx   all data, one sheet per type
├── output/
│   ├── figures/                    written by the script
│   ├── tables/
│   └── models/                     fitted MCMCglmm objects (.rds)
├── CITATION.cff
├── LICENSE                         code, MIT
└── LICENSE-data                    data, CC BY 4.0
```

## Running the analyses

```r
source("bryophytes_itex_analysis.R")
```

Paths resolve with `here::here()` relative to the repository root, so the project
runs from anywhere without editing.

Fitting all MCMCglmm models takes several hours. Each fitted model is saved to
`output/models/` as `.rds`, so to regenerate figures and tables without
refitting, set the flag at the top of the script:

```r
refit_models <- FALSE
```

Each model is preceded by an explicit `set.seed()` drawn from the `seeds` vector
in section 1, so posterior summaries reproduce exactly rather than approximately.

### Requirements

R ≥ 4.3 and:

```r
install.packages(c("here", "readxl", "dplyr", "tidyr", "tibble", "MASS",
                   "MCMCglmm", "broom.mixed", "vegan", "ggplot2", "ggpubr",
                   "ggh4x", "patchwork", "cowplot", "openxlsx"))
```

Versions used for the submitted results are written to `sessionInfo.txt` when the
script finishes.

Two loading notes: `MASS::select()` masks `dplyr::select()` — `MASS` is needed for
`stepAIC()`, so call `dplyr::select()` explicitly. `plot_grid()` is always called
as `cowplot::plot_grid()`.

## Where each result comes from

The script is one file in twelve numbered sections.

| Output | Section |
|--------|---------|
| Figure 2 — cover and layer depth | 5 |
| Figure 3, Tables S3, S4 — vascular plant regressions | 6 |
| Figure 4, Table S5 — species PCA and PERMANOVA | 7 |
| Figure 5 — richness and Shannon diversity | 8 |
| Figure 6, Table S6 — BFG PCA and PERMANOVA | 9 |
| Figure 7, Figure S2 — bryophyte functional groups | 10 |
| Figures 8, 9, 10 — growth, biomass, N₂ fixation, WHC, colony density | 11 |
| Posterior summaries for all models | 12 |

Figure 1 is a map and is not produced by this code (credit: Arto Vitikka, Arctic
Centre, University of Lapland).

**Not included here:** Figure S1 (relative abundance of major growth forms) draws
on the wider plant community datasets of Jónsdóttir et al. (2005, 2023) rather
than on the bryophyte data archived here; figure credit Katrín Björnsdóttir.
Figure S3 (annual growth by the cranked wire method) was analysed separately with
a linear mixed model.

## Data

One workbook, `data/raw/bryophytes_itex_data.xlsx`, with a `README` sheet
describing every column.

| Sheet | Contents | Years |
|-------|----------|-------|
| `point_intercept` | One row per intercepted taxon per point; 100 points per plot, so one hit = 1% cover | 2014–2019 |
| `canopy_height` | Vascular plant canopy height, several measurements per plot | 2014–2019 |
| `layer_depth` | Bryophyte layer depth from 25 pin insertions per plot | 2014–2019 |
| `shoot_growth` | Mesh-bag transplant shoot growth and dry weight | 2015 |
| `n_fixation` | Acetylene reduction assay | 2013–2014 |
| `whc_colony_density` | Colony water holding capacity and density | 2020 |
| `species_matrix` | Plot × taxon matrix of bryophyte hits | 1995–2019 |
| `bfg_matrix` | Plot × functional group matrix, sensu Lett et al. (2022) | 1995–2019 |

Every sheet shares the same five identifiers: `SITE`, `HABITAT`, `TREATMENT`
(`CTL` or `OTC`), `PLOT` and `YEAR`.

The point-intercept sheet covers 2014 onwards only; composition for the earlier
surveys is held in the species and BFG matrices, which therefore cannot be
rebuilt from it.

The Betula heath also carries a grazed (`GRA`) treatment from a separate
experiment. It is excluded from every analysis in the paper and is not included
here.

## Statistical approach

Bayesian generalized linear mixed models were fitted with **MCMCglmm** (Hadfield
2010) using cell-means parameterisation (`~ GROUP - 1`), so each habitat ×
treatment combination is estimated directly rather than as a contrast against a
reference level. Continuous responses use a Gaussian family; counts (cover,
richness, BFG hits) use Poisson with a log link, back-transformed with `exp()`
applied to the posterior draws before summarising, which is what the figures
plot.

Warming effects reported in the paper are based on whether the 95% credible
intervals of the control and OTC posteriors overlap. The results workbook also
reports the posterior of the per-draw OTC − CTL difference and whether its
interval excludes zero; this is provided for transparency and is not an
additional test.

PERMANOVAs use `vegan::adonis2()` on Bray–Curtis dissimilarities with year and
treatment as marginal terms (`by = "margin"`) and 10,000 permutations, fitted
separately per habitat. Year enters as a number, giving a one-degree-of-freedom
linear trend through time. Multivariate dispersion is tested with `betadisper()`
followed by `permutest()`. The moss heath (MH) is excluded from the ordinations
and PERMANOVAs: its bryophyte layer is a single species, and so a single
functional group.

Vascular plant regressions are ordinary least squares, reduced by backward
selection on AIC with a ΔAIC > 2 retention rule (Burnham & Anderson 2002). Tables
report the full models; the retained terms are printed to the console as the
selection runs.

Biological soil crust (`CRUST`) is one of the columns of `species_matrix` and is
therefore counted in the richness, Shannon, PERMANOVA and PCA results, as in the
submitted analysis. It occurs in the three Svalbard habitats only. The functional
group ordinations use all columns of `bfg_matrix` except `Lichen`; `CRUST` is
retained there too.

## Licence

- Code: MIT (`LICENSE`)
- Data: CC BY 4.0 (`LICENSE-data`)

## Acknowledgements

Fieldwork and laboratory work were supported by Vísindasjóður (Rannís), UNIS, The
Icelandic Research Fund, Svalbardstiftelsen and the University of Iceland
Research Fund. Logistical support was provided by Landsvirkjun, Thingvellir
National Park and UNIS. See the paper for the full list of people who contributed
to the community analyses.
