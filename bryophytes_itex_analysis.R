## #####################################################################
##
##  Bryophyte community responses to long-term warming manipulations
##  vary greatly between contrasting tundra habitats
##
##  Jonsdottir, I. S., Klarenberg, I. & Russi Colmenares, A. J.
##
##  Analysis code for all figures and tables in the paper and its
##  supporting information.
##
##  Habitats: BH = Betula nana heath (Audkuluheidi, Iceland)
##            MH = Racomitrium moss heath (Thingvellir, Iceland)
##            DH = Dryas heath (Endalen, Svalbard)
##            CH = Cassiope heath (Endalen, Svalbard)
##            SB = snowbed (Endalen, Svalbard)
##
##  Treatments: CTL = control, OTC = open top chamber (warmed).
##  The Betula heath also carries a grazed (GRA) treatment from a
##  separate experiment. It is excluded from every analysis in the paper
##  and is not present in the published data.
##
##  CONTENTS
##    1  Setup: packages, paths, seeds, plotting defaults
##    2  Helper functions
##    3  Load data
##    4  Derive per-plot variables
##    5  Bryophyte cover and layer depth          -> Figure 2
##    6  Vascular plant regressions               -> Figure 3, Tables S3, S4
##    7  Species composition: PERMANOVA and PCA   -> Figure 4, Table S5
##    8  Richness and Shannon diversity           -> Figure 5
##    9  BFG composition: PERMANOVA and PCA       -> Figure 6, Table S6
##   10  BFG hit counts                           -> Figure 7, Figure S2
##   11  Functional traits                        -> Figures 8, 9, 10
##   12  Results workbook
##
##  Figure 1 is a map and is not produced here.
##  Figures S1 and S3 are analysed elsewhere; see README.md.
##
##  Fitting all models takes several hours. Set refit_models <- FALSE to
##  load saved posteriors from output/models/ and regenerate figures and
##  tables only.
##
## #####################################################################

refit_models <- TRUE


## #####################################################################
##  1  SETUP
## #####################################################################

library(here)          # project-relative paths
library(readxl)        # read raw .xlsx data
library(dplyr)         # data manipulation
library(tidyr)         # pivot_longer / pivot_wider
library(tibble)        # rownames_to_column, used by the PCA sections
library(MASS)          # stepAIC() for backward selection (Table S4)
library(MCMCglmm)      # Bayesian GLMMs (Hadfield 2010)
library(broom.mixed)   # tidy() for MCMCglmm objects
library(vegan)         # adonis2, betadisper, permutest, rda, diversity
library(ggplot2)       # figures
library(ggpubr)        # theme_pubr(), used throughout
library(ggh4x)         # facet_nested() / facet_nested_wrap()
library(patchwork)     # multi-panel assembly (Figures 2, 4, 6)
library(cowplot)       # plot_grid(); always called as cowplot::plot_grid()
library(openxlsx)      # results workbooks

## MASS::select() masks dplyr::select(). Call dplyr::select() explicitly.

## ---- Paths ----------------------------------------------------------
dir_raw     <- here("data", "raw")
dir_figures <- here("output", "figures")
dir_tables  <- here("output", "tables")
dir_models  <- here("output", "models")

for (d in c(dir_figures, dir_tables, dir_models)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

data_file <- file.path(dir_raw, "bryophytes_itex_data.xlsx")

if (!file.exists(data_file)) {
  stop("Data file not found: ", data_file)
}

## ---- Reproducibility ------------------------------------------------
## One seed per model, set immediately before each MCMCglmm() call.
seeds <- c(
  cover = 1001, depth = 1002, richness = 1003, shannon = 1004,
  growth = 1005, biomass_ice = 1006, biomass_sval = 1007,
  whc = 1008, colony_density = 1009,
  nfix_iceland = 1010, nfix_svalbard = 1011,
  bfg_BT = 1012, bfg_We = 1013, bfg_SU = 1014, bfg_TU = 1015
)

seed_permanova <- 42     # as used for the submitted PERMANOVA results
n_perm         <- 10000

## ---- Conventions ----------------------------------------------------
habitat_levels   <- c("BH", "MH", "DH", "CH", "SB")

## Habitats entering the ordinations and PERMANOVAs. MH is excluded: its
## bryophyte layer is a single species, and so a single functional group,
## which leaves nothing to ordinate and makes the PERMANOVA unreliable.
ordination_habitats <- c("BH", "DH", "CH", "SB")
treatment_levels <- c("CTL", "OTC")
trt_cols         <- c(CTL = "#2166ac", OTC = "#cb181d")

region_of <- function(h) ifelse(h %in% c("BH", "MH"), "Iceland", "Svalbard")

## Survey years, shared across all PCA panels so the colour scale is
## identical between habitats (used with drop = FALSE on the scales).
all_years <- c("1995","1997","1998","2000","2003","2007","2009","2014","2015","2019")

## Binned year palette: adjacent survey years share a colour so that the
## Icelandic and Svalbard panels are visually comparable.
year_cols <- c(
  "1995" = "#440154FF", "1997" = "#440154FF", "1998" = "#440154FF",
  "2000" = "#31688EFF", "2003" = "#31688EFF",
  "2007" = "#1F9E89FF", "2009" = "#1F9E89FF",
  "2014" = "#6DCD59FF", "2015" = "#6DCD59FF",
  "2019" = "#FDE725FF"
)

## Palette for the stacked BFG bars (Figure 7)
bfg_palette <- c("#DBA04C","#514E1D","#7E8E2C","#B8AD62","#457D60",
                 "#55588B","#664151","#8C7B6B","#843810","#BD5630FF")

theme_set(theme_pubr())

save_fig <- function(plot, filename, width, height, dpi = 400) {
  ggsave(file.path(dir_figures, filename), plot = plot,
         width = width, height = height, dpi = dpi)
}


## #####################################################################
##  2  HELPER FUNCTIONS
## #####################################################################

#' Split a wide community matrix into metadata and taxon columns by name
#'
#' Replaces positional indexing such as Moss_SPECIES[, 6:32], which fails
#' silently if the column order of the published file differs by one.
split_community <- function(df, meta_cols, drop_cols = character(0)) {
  absent <- setdiff(meta_cols, names(df))
  if (length(absent)) stop("Metadata columns not found: ", paste(absent, collapse = ", "))

  taxa_cols <- setdiff(names(df), c(meta_cols, drop_cols))
  spp <- as.data.frame(df[, taxa_cols, drop = FALSE])

  ## Hits read in as text where a column is sparse; coerce, then treat
  ## absences as the structural zeros they are.
  spp[] <- lapply(spp, function(x) { x <- as.numeric(x); x[is.na(x)] <- 0; x })

  list(meta = as.data.frame(df[, meta_cols, drop = FALSE]), spp = spp)
}

#' PERMANOVA with marginal terms plus dispersion tests for both factors
#'
#' by = "margin" tests year and treatment against each other. Dispersion
#' is tested per factor with betadisper() + permutest().
permanova_habitat <- function(spp, env, label) {
  empty <- rowSums(spp) == 0
  if (any(empty)) {
    message(label, ": dropped ", sum(empty), " plot(s) with zero total hits")
    spp <- spp[!empty, , drop = FALSE]
    env <- env[!empty, , drop = FALSE]
  }
  spp <- spp[, colSums(spp) > 0, drop = FALSE]
  d <- vegan::vegdist(spp, method = "bray")

  set.seed(seed_permanova)
  bd_trt  <- vegan::betadisper(d, group = as.factor(env$TREATMENT))
  bd_year <- vegan::betadisper(d, group = as.factor(env$YEAR))
  p_trt   <- vegan::permutest(bd_trt,  permutations = 999)$tab[1, "Pr(>F)"]
  p_year  <- vegan::permutest(bd_year, permutations = 999)$tab[1, "Pr(>F)"]

  ## YEAR enters as a NUMBER, not a factor: the published models fit a
  ## linear trend through time with one degree of freedom. Passing it as
  ## a factor costs one df per extra survey year and changes every term.
  env$YEAR <- as.numeric(as.character(env$YEAR))

  set.seed(seed_permanova)
  ad  <- vegan::adonis2(d ~ YEAR + TREATMENT, data = env,
                        permutations = n_perm, by = "margin")
  res <- as.data.frame(ad)

  homog <- rep(NA_real_, nrow(res))
  homog[match("YEAR",      rownames(res))] <- round(p_year, 3)
  homog[match("TREATMENT", rownames(res))] <- round(p_trt,  3)

  data.frame(
    Habitat       = c(label, rep(NA, nrow(res) - 1)),
    Term          = rownames(res),
    DF            = res$Df,
    SumOfSqs      = round(res$SumOfSqs, 4),
    R2            = round(res$R2, 4),
    F             = round(res$F, 4),
    P             = signif(res$`Pr(>F)`, 4),
    Homogeneity_P = homog,
    row.names     = NULL
  )
}

#' PCA of a community matrix, with site scores joined to the metadata
pca_scores <- function(spp, env) {
  pca <- vegan::rda(spp)
  u   <- as.data.frame(pca$CA$u)
  out <- dplyr::inner_join(tibble::rownames_to_column(env),
                           tibble::rownames_to_column(u), by = "rowname")
  out$YEAR <- factor(as.character(out$YEAR), levels = all_years)
  list(
    scores = out,
    var1   = round(pca$CA$eig[1] / sum(pca$CA$eig) * 100, 1),
    var2   = round(pca$CA$eig[2] / sum(pca$CA$eig) * 100, 1)
  )
}

#' One PCA panel: plots coloured by year, shaped by treatment, with
#' treatment centroids joined chronologically by arrows.
pca_panel <- function(p, title, annotation = NULL, xlim, ylim) {
  centroids <- p$scores %>%
    dplyr::group_by(TREATMENT, YEAR) %>%
    dplyr::summarise(PC1 = mean(PC1, na.rm = TRUE),
                     PC2 = mean(PC2, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(TREATMENT, as.numeric(as.character(YEAR)))

  g <- ggplot(p$scores, aes(x = PC1, y = PC2,
                            fill = YEAR, color = YEAR, shape = TREATMENT)) +
    geom_point(size = 2.5, alpha = 0.5) +
    scale_fill_manual(values = year_cols, name = "Year", drop = FALSE) +
    scale_color_manual(values = year_cols, name = "Year", drop = FALSE) +
    scale_shape_manual(values = c(CTL = 21, OTC = 24),
                       labels = c(CTL = "Control", OTC = "OTC"),
                       name = "Treatment") +
    geom_path(data = centroids, aes(group = TREATMENT),
              arrow = arrow(length = unit(0.15, "cm"), type = "closed", ends = "last"),
              alpha = 0.5, show.legend = FALSE, color = "black") +
    geom_point(data = centroids, size = 3, alpha = 1, color = "black") +
    theme_classic() +
    coord_equal(xlim = xlim, ylim = ylim) +
    guides(fill  = guide_legend(override.aes = list(shape = 21)),
           shape = guide_legend(override.aes = list(fill = "black"))) +
    labs(x = paste0("PC1 (", p$var1, "%)"),
         y = paste0("PC2 (", p$var2, "%)")) +
    ggtitle(title) +
    theme(aspect.ratio = 1)

  if (!is.null(annotation)) {
    g <- g + annotate("text", x = xlim[1], y = ylim[2], hjust = 0, vjust = 1,
                      size = 4, label = annotation)
  }
  g
}

#' Coefficient table for one lm, in the layout of Table S4
lm_table <- function(fit, response, habitat_year, note = NA_character_) {
  s  <- summary(fit)
  co <- as.data.frame(coef(s))
  ci <- as.data.frame(confint(fit))

  data.frame(
    Response  = c(response, rep(NA, nrow(co) - 1)),
    Group     = c(habitat_year, rep(NA, nrow(co) - 1)),
    Predictor = rownames(co),
    Estimate  = round(co[["Estimate"]], 3),
    SE        = round(co[["Std. Error"]], 3),
    CI_low    = round(ci[[1]], 3),
    CI_high   = round(ci[[2]], 3),
    Statistic = round(co[["t value"]], 2),
    P         = signif(co[["Pr(>|t|)"]], 3),
    N         = c(length(fit$residuals), rep(NA, nrow(co) - 1)),
    R2        = c(round(s$r.squared, 3), rep(NA, nrow(co) - 1)),
    R2_adj    = c(round(s$adj.r.squared, 3), rep(NA, nrow(co) - 1)),
    Note      = c(note, rep(NA, nrow(co) - 1)),
    row.names = NULL
  )
}

#' Fit a full model, run backward AIC selection, and report both
#'
#' The paper reports the full model (Table S4) with AIC selection as a
#' check on which predictors are retained, following the rule that a
#' predictor is dropped only where removal raises AIC by more than two.
fit_and_select <- function(formula, data, response, habitat_year, note = NA_character_) {
  fit  <- lm(formula, data = data)
  step <- MASS::stepAIC(fit, direction = "backward", trace = FALSE)
  message(habitat_year, " (", response, ") retained: ",
          paste(attr(terms(step), "term.labels"), collapse = " + "))
  lm_table(fit, response, habitat_year, note)
}

#' Tidy a cell-means coefficient name into a readable label
tidy_cell <- function(x) {
  x <- gsub("Loc_Year_TRT|Loc_TRT_plot|Loc_TRT|SITE_YEAR_TRT|SITE_TRT_Month|SITE_TRT",
            "", x)
  x <- gsub("^[_.]+|[_.]+$", "", x)

  repl <- c("DRY-L" = "DH", "CAS-L" = "CH", "BIS-L" = "SB",
            "DRY" = "DH", "CAS" = "CH", "BIS" = "SB",
            "AUD" = "BH", "THI" = "MH",
            "Audkuluheidi" = "BH", "Thingvellir" = "MH",
            "Endalen_Dryas" = "DH", "Endalen_Cassiope" = "CH",
            "Endalen_Snowbed" = "SB",
            "EndalenDryas" = "DH", "EndalenCassiope" = "CH",
            "EndalenSnowbed" = "SB",
            "DH" = "DH", "CH" = "CH", "SB" = "SB", "BH" = "BH", "MH" = "MH")

  parts <- unlist(strsplit(x, "[_.]"))
  parts <- parts[nzchar(parts)]
  parts <- vapply(parts, function(p) if (p %in% names(repl)) repl[[p]] else p, "")
  trimws(paste(parts, collapse = " "))
}

#' Per-treatment posterior intervals and the OTC - CTL contrast
#'
#' Inference in the paper uses the credible-interval overlap criterion.
#' The difference columns are provided for transparency, not as a second
#' test.
#'
#' Poisson models are back-transformed with exp() applied to the draws
#' before summarising, which is how the figures are plotted. Applying
#' exp() to an already-summarised value would not match them.
mcmc_diff <- function(model, family = c("gaussian", "poisson"),
                      ctl_label = "CTL", otc_label = "OTC",
                      model_name = "") {

  family <- match.arg(family)
  sol <- as.matrix(model$Sol)
  cn  <- colnames(sol)

  ## Remove the treatment token entirely. Leaving a placeholder in its
  ## position lets the token leak into the printed label.
  stem_of <- function(x) {
    s <- gsub(ctl_label, "", x, fixed = TRUE)
    s <- gsub(otc_label, "", s, fixed = TRUE)
    gsub("__+", "_", s)
  }

  ctl_cols <- cn[grepl(ctl_label, cn, fixed = TRUE)]
  otc_cols <- cn[grepl(otc_label, cn, fixed = TRUE)]
  if (length(ctl_cols) == 0 || length(otc_cols) == 0) {
    stop(model_name, ": no CTL/OTC columns. First coefficients: ",
         paste(utils::head(cn), collapse = " | "))
  }

  ctl_map  <- setNames(ctl_cols, vapply(ctl_cols, stem_of, ""))
  otc_map  <- setNames(otc_cols, vapply(otc_cols, stem_of, ""))
  paired   <- intersect(names(ctl_map), names(otc_map))
  unpaired <- setdiff(union(names(ctl_map), names(otc_map)), paired)

  ## Unpaired cells are expected where a treatment was not measured,
  ## e.g. N2 fixation in Iceland in July (control plots only).
  if (length(unpaired)) {
    message("  ", model_name, ": dropped (no CTL/OTC partner): ",
            paste(vapply(unpaired, tidy_cell, ""), collapse = "; "))
  }
  if (length(paired) == 0) stop(model_name, ": no CTL/OTC pairs.")

  back <- if (family == "poisson") exp else identity

  out <- lapply(paired, function(s) {
    ctl <- back(as.numeric(sol[, ctl_map[[s]]]))
    otc <- back(as.numeric(sol[, otc_map[[s]]]))
    cq  <- quantile(ctl, c(.025, .975))
    oq  <- quantile(otc, c(.025, .975))
    d   <- otc - ctl
    dq  <- quantile(d, c(.025, .975))

    data.frame(
      cell        = tidy_cell(s),
      CTL_mean    = mean(ctl), CTL_low = cq[1], CTL_high = cq[2],
      OTC_mean    = mean(otc), OTC_low = oq[1], OTC_high = oq[2],
      CIs_overlap = if (!(cq[1] > oq[2] || oq[1] > cq[2])) "yes" else "no",
      diff_mean   = mean(d), diff_low = dq[1], diff_high = dq[2],
      diff_excludes_zero = if (dq[1] > 0 || dq[2] < 0) "yes" else "no",
      diff_type   = "OTC - CTL (figure units)",
      row.names   = NULL
    )
  })

  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}

#' Fit a model, or load it from output/models/ when refit_models is FALSE
fit_or_load <- function(name, seed, expr) {
  path <- file.path(dir_models, paste0(name, ".rds"))
  if (!refit_models && file.exists(path)) {
    message("Loading ", name)
    return(readRDS(path))
  }
  message("Fitting ", name)
  set.seed(seed)
  model <- force(expr)
  saveRDS(model, path)
  model
}


## #####################################################################
##  3  LOAD DATA
##
##  All data live in one workbook, one sheet per data type. Identifiers
##  are already harmonised across sheets:
##
##    SITE       AUDKULUHEIDI | THINGVELLIR | ENDALEN
##    HABITAT    BH | MH | DH | CH | SB
##    TREATMENT  CTL | OTC
##    PLOT       plot number within habitat
##    YEAR       survey or sampling year
##
##  Grazed (GRA) plots are a separate treatment and are not in the data.
##  See the README sheet of the workbook for column definitions.
## #####################################################################

read_sheet <- function(sheet) {
  d <- as.data.frame(read_excel(data_file, sheet = sheet))
  d$HABITAT   <- factor(d$HABITAT, levels = habitat_levels)
  d$TREATMENT <- factor(d$TREATMENT, levels = treatment_levels)
  d$PLOT      <- as.character(d$PLOT)
  if ("YEAR" %in% names(d)) d$YEAR <- as.character(d$YEAR)
  if (anyNA(d$HABITAT) || anyNA(d$TREATMENT)) {
    stop(sheet, ": unexpected habitat or treatment codes")
  }
  d
}

## ---- Point-intercept data -------------------------------------------
## One row per intercepted taxon per point; 100 points per plot, so one
## hit = 1% cover. Source of bryophyte cover, litter and the shrub
## variables. Covers 2014-2019; the earlier surveys are summarised in
## the species and BFG matrices.
Point_Int <- read_sheet("point_intercept")
Point_Int$SITEPLOT <- interaction(Point_Int$SITE, Point_Int$HABITAT,
                                  Point_Int$TREATMENT, Point_Int$PLOT,
                                  Point_Int$YEAR, sep = "_")

## ---- Canopy height ---------------------------------------------------
Canopy_height <- read_sheet("canopy_height")
Canopy_height$SITEPLOT <- interaction(Canopy_height$SITE, Canopy_height$HABITAT,
                                      Canopy_height$TREATMENT, Canopy_height$PLOT,
                                      Canopy_height$YEAR, sep = "_")

## ---- Bryophyte layer depth -------------------------------------------
## 25 pin measurements per plot. Iceland 2014 and 2019, Svalbard 2015.
Moss_depth <- read_sheet("layer_depth")
Moss_depth$Loc_Year_TRT <- interaction(Moss_depth$HABITAT, Moss_depth$YEAR,
                                       Moss_depth$TREATMENT, sep = "_")
Moss_depth$PLOTEFFECT   <- interaction(Moss_depth$PLOT, Moss_depth$TREATMENT, sep = "_")
Moss_depth$SITEPLOT     <- interaction(Moss_depth$SITE, Moss_depth$HABITAT,
                                       Moss_depth$TREATMENT, Moss_depth$PLOT,
                                       Moss_depth$YEAR, sep = "_")

## ---- Annual shoot growth and shoot biomass ---------------------------
## Mesh-bag transplants: R. lanuginosum in Iceland (10 shoots per bag,
## trimmed to 30 mm), S. uncinata in Svalbard (15 shoots, 20 mm), three
## bags per plot. DH is absent: its bryophyte layer is too shallow for
## the method, which is expected, not missing data.
Moss_growth <- read_sheet("shoot_growth")
Moss_growth$Loc_TRT <- interaction(Moss_growth$HABITAT, Moss_growth$TREATMENT, sep = "_")
Moss_growth$PLOT_ID <- interaction(Moss_growth$HABITAT, Moss_growth$TREATMENT,
                                   Moss_growth$PLOT, sep = "_")
Moss_growth$BAG_ID  <- interaction(Moss_growth$HABITAT, Moss_growth$TREATMENT,
                                   Moss_growth$PLOT, Moss_growth$Position, sep = "_")

## ---- Water holding capacity and colony density -----------------------
## R. lanuginosum colonies of standard volume (98.17 cm3), two cores per
## plot, Icelandic habitats only, June 2020.
Moss_WHC <- read_sheet("whc_colony_density")
Moss_WHC$Loc_TRT      <- interaction(Moss_WHC$HABITAT, Moss_WHC$TREATMENT, sep = "_")
Moss_WHC$Loc_TRT_plot <- interaction(Moss_WHC$HABITAT, Moss_WHC$TREATMENT,
                                     Moss_WHC$PLOT, sep = "_")
Moss_WHC$COLEFFECT    <- interaction(Moss_WHC$TREATMENT, Moss_WHC$Colony, sep = "_")

## Note: WHC uses COLEFFECT (colony) as its random effect while colony
## density uses Loc_TRT_plot (plot). Both are kept exactly as fitted for
## the submitted results.

## ---- N2 fixation (acetylene reduction assay) -------------------------
## Iceland: R. lanuginosum, June/July/August 2014, July controls only.
## Svalbard: S. uncinata, early August 2013.
N_fixation <- read_sheet("n_fixation")
N_fixation$REGION <- region_of(N_fixation$HABITAT)

## Iceland was sampled in three months, Svalbard in one, so the two
## regions are modelled separately with different cell structures.
N_fixation$SITE_TRT_Month <- ifelse(
  N_fixation$REGION == "Iceland",
  as.character(interaction(N_fixation$HABITAT, N_fixation$TREATMENT,
                           N_fixation$MONTH, sep = "_")),
  as.character(interaction(N_fixation$HABITAT, N_fixation$TREATMENT, sep = "_"))
)
N_fixation$PLOTEFFECT <- interaction(N_fixation$TREATMENT, N_fixation$PLOT, sep = "_")

## ---- Bryophyte functional groups (wide matrix) -----------------------
## Species assigned to the twelve BFGs of Lett et al. (2022); eight are
## represented here. Lichen and biological soil crust are reported in
## Figure 7; crust is also carried in the ordination, matching the
## submitted analysis.
Moss_BFG <- read_sheet("bfg_matrix")

bfg_meta_cols <- c("SITE", "HABITAT", "TREATMENT", "PLOT", "YEAR")
bfg_all_cols  <- setdiff(names(Moss_BFG), bfg_meta_cols)

## The ordination excludes lichens only; crust is retained.
bfg_ord_cols <- setdiff(bfg_all_cols, "Lichen")

Moss_BFG_long <- Moss_BFG %>%
  tidyr::pivot_longer(cols = all_of(bfg_all_cols),
                      names_to = "BFG", values_to = "HITS") %>%
  as.data.frame()

Moss_BFG_long$HITS <- as.numeric(Moss_BFG_long$HITS)
Moss_BFG_long$BFG  <- factor(Moss_BFG_long$BFG)

## ---- Bryophyte species (wide matrix) ---------------------------------
## Note: CRUST is one of the columns, so biological soil crust is counted
## in richness, Shannon diversity, the species PERMANOVAs and the species
## PCA, as in the submitted analysis.
Moss_SPECIES <- read_sheet("species_matrix")

species_meta_cols <- c("SITE", "HABITAT", "TREATMENT", "PLOT", "YEAR")

message("Loaded: ", nrow(Point_Int), " point-intercept records | ",
        nrow(Moss_depth), " depth measurements | ",
        nrow(Moss_growth), " growth records | ",
        nrow(Moss_WHC), " colonies | ",
        nrow(N_fixation), " ARA samples")


## #####################################################################
##  4  DERIVE PER-PLOT VARIABLES
##
##  Point-intercept records are summed to plot level, canopy height is
##  averaged per plot, and the vascular plant table is assembled.
##  Structural zeros are filled by column NAME after each join: column
##  positions shift with every merge, so positional filling silently
##  writes zeros into the wrong variable.
## #####################################################################

#' Sum point-intercept hits per plot for a set of growth forms or taxa
sum_hits_per_plot <- function(df, col, keep, value_name) {
  df %>%
    dplyr::filter(.data[[col]] %in% keep) %>%
    dplyr::group_by(SITEPLOT, SITE, HABITAT, TREATMENT, PLOT, YEAR) %>%
    dplyr::summarise(!!value_name := sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    as.data.frame()
}

Point_Int_litter  <- sum_hits_per_plot(Point_Int, "SPPORIG", "LITTER", "Litter abundance")
Point_Int_allshrub<- sum_hits_per_plot(Point_Int, "GBROAD", c("ESHRUB","DSHRUB"), "Shrub abundance")
Point_Int_dshrub  <- sum_hits_per_plot(Point_Int, "GBROAD", "DSHRUB", "Deciduous shrub abundance")
Point_Int_eshrub  <- sum_hits_per_plot(Point_Int, "GBROAD", "ESHRUB", "Evergreen shrub abundance")
Point_Int_moss    <- sum_hits_per_plot(Point_Int, "GBROAD", "MOSS",   "Moss cover")
Point_Int_lichen  <- sum_hits_per_plot(Point_Int, "GBROAD", "LICHEN", "Lichen cover")

Canopy_height_plot <- Canopy_height %>%
  dplyr::group_by(SITEPLOT, SITE, HABITAT, TREATMENT, PLOT, YEAR) %>%
  dplyr::summarise(HEIGHT = mean(HEIGHT, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame()

Moss_depth_plot <- Moss_depth %>%
  dplyr::group_by(SITEPLOT, SITE, HABITAT, TREATMENT, PLOT, YEAR) %>%
  dplyr::summarise(`Moss depth` = mean(Depth_cm, na.rm = TRUE),
                   .groups = "drop") %>%
  as.data.frame()

## ---- Assemble the vascular plant table -------------------------------
join_keys <- c("SITEPLOT", "SITE", "HABITAT", "TREATMENT", "PLOT", "YEAR")

vasc_plant_comm <- Point_Int %>%
  dplyr::distinct(dplyr::across(all_of(join_keys))) %>%
  as.data.frame()

for (tbl in list(Point_Int_litter, Canopy_height_plot, Point_Int_allshrub,
                 Point_Int_dshrub, Point_Int_eshrub, Point_Int_moss,
                 Moss_depth_plot)) {
  vasc_plant_comm <- merge(vasc_plant_comm, tbl, by = join_keys, all.x = TRUE)
}

## Absence from the point-intercept survey is zero hits, not missing.
## Canopy height and moss depth are genuinely absent where not measured,
## so they keep their NAs.
abundance_vars <- c("Litter abundance", "Shrub abundance",
                    "Deciduous shrub abundance", "Evergreen shrub abundance",
                    "Moss cover", "Lichen cover")
present <- intersect(abundance_vars, names(vasc_plant_comm))
vasc_plant_comm[present] <- lapply(vasc_plant_comm[present],
                                   function(x) { x <- as.numeric(x); x[is.na(x)] <- 0; x })

vasc_plant_comm$SUBSITEYEAR <- factor(paste(vasc_plant_comm$HABITAT,
                                            vasc_plant_comm$YEAR, sep = "."))

## ---- Richness and Shannon diversity ----------------------------------
## Computed once here, by name. The original indexed species columns as
## [, 6:32]; if the published file's column order differs by one, that
## silently changes both metrics and Figure 5 stops reproducing.
species_parts <- split_community(Moss_SPECIES, species_meta_cols)

Moss_SPECIES$Richness <- vegan::specnumber(species_parts$spp)
Moss_SPECIES$Shannon  <- vegan::diversity(species_parts$spp, index = "shannon")
Moss_SPECIES$Island   <- region_of(Moss_SPECIES$HABITAT)
Moss_SPECIES$SITE_TRT <- interaction(Moss_SPECIES$HABITAT, Moss_SPECIES$TREATMENT,
                                     sep = "_")

stopifnot(!anyNA(Moss_SPECIES$Richness), !anyNA(Moss_BFG_long$HITS))


## #####################################################################
##  5  BRYOPHYTE COVER AND LAYER DEPTH            -> Figure 2
## #####################################################################

## ---- Cover: Poisson cell-means, habitat x year x treatment -----------
Point_Int_moss$Loc_Year_TRT <- interaction(Point_Int_moss$HABITAT,
                                           Point_Int_moss$TREATMENT,
                                           Point_Int_moss$YEAR, sep = "_")

prior_R  <- list(R = list(V = 1, nu = 0.002))
prior_G1 <- list(G = list(G1 = list(V = 1, nu = 0.002)),
                 R = list(V = 1, nu = 0.002))

iter <- 400
mosscovermcmcglmm <- fit_or_load("cover", seeds[["cover"]],
  MCMCglmm(`Moss cover` ~ Loc_Year_TRT - 1, rcov = ~ units,
           data   = Point_Int_moss,
           prior  = prior_R,
           nitt   = 13000 * iter, burnin = 3000 * iter, thin = 10 * iter,
           family = "poisson", verbose = FALSE))

mosscover_df <- broom.mixed::tidy(mosscovermcmcglmm, conf.int = TRUE, conf.level = 0.95)
mosscover_df <- subset(mosscover_df, effect == "fixed")

## Derive habitat, treatment and year from the coefficient names rather
## than from a hand-written vector. The original assigned these by
## position, which breaks if a cell is dropped or the level order shifts.
cover_parts <- do.call(rbind, strsplit(sub("Loc_Year_TRT", "", mosscover_df$term), "_"))
mosscover_df$location  <- factor(cover_parts[, 1], levels = habitat_levels)
mosscover_df$treatment <- cover_parts[, 2]
mosscover_df$year      <- cover_parts[, 3]
mosscover_df$Island    <- region_of(mosscover_df$location)

## ---- Layer depth: Gaussian, plot random ------------------------------
mossdepth <- fit_or_load("depth", seeds[["depth"]],
  MCMCglmm(Depth_cm ~ Loc_Year_TRT - 1,
           random = ~ PLOTEFFECT, rcov = ~ units,
           data   = Moss_depth,
           prior  = prior_G1,
           nitt   = 130000, burnin = 30000, thin = 100,
           family = "gaussian", verbose = FALSE))

mossdepth_df <- broom.mixed::tidy(mossdepth, conf.int = TRUE, conf.level = 0.95)
mossdepth_df <- subset(mossdepth_df, effect == "fixed")

depth_parts <- do.call(rbind, strsplit(sub("Loc_Year_TRT", "", mossdepth_df$term), "_"))
mossdepth_df$location  <- factor(depth_parts[, 1], levels = habitat_levels)
mossdepth_df$year      <- depth_parts[, 2]
mossdepth_df$treatment <- depth_parts[, 3]
mossdepth_df$Island    <- region_of(mossdepth_df$location)

## Depth cannot be negative; the lower credible bound is truncated at 0
## for plotting only, never for inference.
mossdepth_df$conf.low2 <- pmax(mossdepth_df$conf.low, 0)

## ---- Figure 2 --------------------------------------------------------
cover_panel <- function(region) {
  ggplot(subset(mosscover_df, Island == region),
         aes(x = location, y = exp(estimate),
             ymin = exp(conf.low), ymax = exp(conf.high), color = treatment)) +
    geom_pointrange(position = position_dodge(width = 0.3), size = 0.3) +
    xlab("") + ylab("Bryophyte cover (%)") +
    facet_wrap(~year, nrow = 1) +
    theme_pubr() +
    theme(strip.background = element_rect(color = "white", fill = "white", linetype = "blank"),
          strip.text = element_text(size = 12, colour = "black"),
          legend.position = "none") +
    scale_color_manual(values = trt_cols) +
    ylim(0, 110)
}

depth_panel <- function(region) {
  ggplot(subset(mossdepth_df, Island == region),
         aes(x = location, y = estimate,
             ymin = conf.low2, ymax = conf.high, color = treatment)) +
    geom_pointrange(position = position_dodge(width = 0.3), size = 0.3) +
    xlab("") + ylab("Bryophyte layer depth (cm)") +
    facet_wrap(~year, nrow = 1) +
    theme_pubr() +
    theme(strip.background = element_rect(color = "white", fill = "white", linetype = "blank"),
          strip.text = element_text(size = 12, colour = "black"),
          legend.position = "none") +
    scale_color_manual(values = trt_cols) +
    ylim(0, 18)
}

fig2 <- (cover_panel("Iceland") | cover_panel("Svalbard")) /
        (depth_panel("Iceland") | depth_panel("Svalbard")) +
  plot_annotation(tag_levels = "a", tag_suffix = ")")

save_fig(fig2, "Figure2_cover_and_depth.png", width = 10, height = 7)


## #####################################################################
##  6  VASCULAR PLANT REGRESSIONS      -> Figure 3, Tables S3, S4
##
##  OLS regressions of cover and layer depth on canopy height, litter,
##  and total / deciduous / evergreen shrub abundance, by habitat and
##  year. Backward selection on AIC, dropping a predictor only where
##  removal raises AIC by more than two (Burnham & Anderson 2002).
##  Table S4 reports the full models; the retained terms are printed to
##  the console as the selection runs.
## #####################################################################

vp <- function(group, extra = TRUE) {
  subset(vasc_plant_comm, SUBSITEYEAR == group)
}

tableS4 <- rbind(
  ## --- Bryophyte cover ---
  fit_and_select(`Moss cover` ~ `Litter abundance` + HEIGHT +
                   `Shrub abundance` + `Deciduous shrub abundance`,
                 vp("BH.2014"), "Bryophyte cover", "BH 2014"),

  fit_and_select(`Moss cover` ~ `Litter abundance` + HEIGHT +
                   `Shrub abundance` + `Deciduous shrub abundance`,
                 vp("BH.2019"), "Bryophyte cover", "BH 2019"),

  fit_and_select(`Moss cover` ~ `Litter abundance` + HEIGHT,
                 vp("MH.2014"), "Bryophyte cover", "MH 2014"),

  fit_and_select(`Moss cover` ~ HEIGHT + `Litter abundance`,
                 vp("MH.2019"), "Bryophyte cover", "MH 2019"),

  ## One MH 2019 plot has unusually high litter. The relationship is
  ## reported with and without it; both remain significant.
  fit_and_select(`Moss cover` ~ HEIGHT + `Litter abundance`,
                 subset(vp("MH.2019"), `Litter abundance` < 80),
                 "Bryophyte cover", "MH 2019",
                 note = "outlier removed (litter abundance 88)"),

  fit_and_select(`Moss cover` ~ HEIGHT + `Litter abundance` +
                   `Deciduous shrub abundance` + `Shrub abundance`,
                 vp("DH.2015"), "Bryophyte cover", "DH 2015", note = "model NS"),

  fit_and_select(`Moss cover` ~ HEIGHT + `Litter abundance` +
                   `Deciduous shrub abundance` + `Shrub abundance`,
                 vp("CH.2015"), "Bryophyte cover", "CH 2015", note = "model NS"),

  fit_and_select(`Moss cover` ~ HEIGHT + `Litter abundance` +
                   `Shrub abundance` + `Evergreen shrub abundance`,
                 vp("SB.2015"), "Bryophyte cover", "SB 2015", note = "model NS"),

  ## --- Bryophyte layer depth ---
  fit_and_select(`Moss depth` ~ `Litter abundance` + HEIGHT +
                   `Shrub abundance` + `Deciduous shrub abundance`,
                 vp("BH.2014"), "Bryophyte depth", "BH 2014"),

  fit_and_select(`Moss depth` ~ `Litter abundance` + HEIGHT +
                   `Shrub abundance` + `Deciduous shrub abundance`,
                 vp("BH.2019"), "Bryophyte depth", "BH 2019"),

  fit_and_select(`Moss depth` ~ `Litter abundance` + HEIGHT,
                 vp("MH.2014"), "Bryophyte depth", "MH 2014"),

  fit_and_select(`Moss depth` ~ `Litter abundance` + HEIGHT,
                 vp("MH.2019"), "Bryophyte depth", "MH 2019", note = "model NS"),

  fit_and_select(`Moss depth` ~ HEIGHT + `Deciduous shrub abundance`,
                 vp("DH.2015"), "Bryophyte depth", "DH 2015", note = "model NS"),

  fit_and_select(`Moss depth` ~ `Deciduous shrub abundance`,
                 vp("CH.2015"), "Bryophyte depth", "CH 2015"),

  fit_and_select(`Moss depth` ~ `Litter abundance` + HEIGHT +
                   `Shrub abundance` + `Deciduous shrub abundance`,
                 vp("CH.2015"), "Bryophyte depth", "CH 2015",
                 note = "full model, NS"),

  fit_and_select(`Moss depth` ~ `Litter abundance` + HEIGHT +
                   `Shrub abundance` + `Deciduous shrub abundance`,
                 vp("SB.2015"), "Bryophyte depth", "SB 2015", note = "model NS")
)

write.xlsx(tableS4, file.path(dir_tables, "TableS4_regressions.xlsx"),
           overwrite = TRUE)

## ---- Table S3: vascular plant variables, mean +- SD ------------------
## Last survey year per habitat: 2019 in Iceland, 2015 in Svalbard.
vasc_plant_lastyear <- subset(
  vasc_plant_comm,
  ((HABITAT %in% c("BH","MH") & YEAR == "2019") |
     (HABITAT %in% c("DH","CH","SB") & YEAR == "2015"))
)

var_levels <- c("Canopy height (cm)",
                "Deciduous shrub abundance (hits)",
                "Evergreen shrub abundance (hits)",
                "Shrub abundance (hits)",
                "Litter abundance (hits)")

vasc_long <- vasc_plant_lastyear %>%
  dplyr::select(HABITAT, YEAR, TREATMENT, HEIGHT,
                `Deciduous shrub abundance`, `Evergreen shrub abundance`,
                `Shrub abundance`, `Litter abundance`) %>%
  tidyr::pivot_longer(cols = -c(HABITAT, YEAR, TREATMENT),
                      names_to = "variable", values_to = "value") %>%
  dplyr::mutate(
    value    = as.numeric(value),
    variable = dplyr::recode(variable,
      "HEIGHT"                    = "Canopy height (cm)",
      "Deciduous shrub abundance" = "Deciduous shrub abundance (hits)",
      "Evergreen shrub abundance" = "Evergreen shrub abundance (hits)",
      "Shrub abundance"           = "Shrub abundance (hits)",
      "Litter abundance"          = "Litter abundance (hits)"),
    variable = factor(variable, levels = var_levels)
  )

vasc_summary <- vasc_long %>%
  dplyr::group_by(HABITAT, YEAR, TREATMENT, variable) %>%
  dplyr::summarise(mean = mean(value, na.rm = TRUE),
                   sd   = sd(value,  na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(cell = ifelse(is.nan(mean), "—",
                              sprintf("%.2f ± %.2f", mean, sd)))

plot_n <- vasc_plant_lastyear %>%
  dplyr::group_by(HABITAT, YEAR, TREATMENT) %>%
  dplyr::summarise(n = dplyr::n_distinct(PLOT), .groups = "drop")

tableS3 <- vasc_summary %>%
  dplyr::mutate(Island = region_of(HABITAT)) %>%
  dplyr::select(Island, HABITAT, YEAR, TREATMENT, variable, cell) %>%
  tidyr::pivot_wider(names_from = variable, values_from = cell) %>%
  dplyr::left_join(plot_n, by = c("HABITAT","YEAR","TREATMENT")) %>%
  dplyr::relocate(n, .after = TREATMENT) %>%
  dplyr::arrange(HABITAT, TREATMENT) %>%
  dplyr::rename(Habitat = HABITAT, Year = YEAR, Treatment = TREATMENT) %>%
  as.data.frame()

write.xlsx(tableS3, file.path(dir_tables, "TableS3_vascular_plant_variables.xlsx"),
           overwrite = TRUE)

## ---- Figure 3: the four significant relationships --------------------
reg_panel <- function(group, yvar, xvar, ylab, xlab, title, filter_expr = TRUE) {
  d <- subset(vasc_plant_comm, SUBSITEYEAR == group)
  d <- d[eval(substitute(filter_expr), d, parent.frame()), , drop = FALSE]
  ggplot(d, aes(y = .data[[yvar]], x = .data[[xvar]])) +
    geom_point() +
    stat_smooth(method = "lm") +
    theme_pubr() +
    ylab(ylab) + xlab(xlab) + ggtitle(title)
}

fig3a <- reg_panel("BH.2019", "Moss cover", "Litter abundance",
                   "Bryophyte cover (%)", "Litter abundance (hits)", "BH 2019")
fig3b <- reg_panel("BH.2019", "Moss depth", "Deciduous shrub abundance",
                   "Bryophyte depth (cm)", "Deciduous shrub abundance (hits)", "BH 2019")
fig3c <- reg_panel("MH.2019", "Moss cover", "Litter abundance",
                   "Bryophyte cover (%)", "Litter abundance (hits)", "MH 2019")
fig3d <- reg_panel("CH.2015", "Moss depth", "Deciduous shrub abundance",
                   "Bryophyte depth (cm)", "Deciduous shrub abundance (hits)", "CH 2015")

fig3 <- (fig3a | fig3b) / (fig3c | fig3d) +
  plot_annotation(tag_levels = "a", tag_suffix = ")")

save_fig(fig3, "Figure3_vascular_plant_relationships.png", width = 7, height = 6, dpi = 300)


## #####################################################################
##  7  SPECIES COMPOSITION: PERMANOVA AND PCA   -> Figure 4, Table S5
##
##  MH is excluded (see ordination_habitats in section 1).
## #####################################################################

species_subset <- function(h) {
  keep <- Moss_SPECIES$HABITAT == h
  split_community(Moss_SPECIES[keep, , drop = FALSE],
                  meta_cols = species_meta_cols,
                  drop_cols = c("Richness", "Shannon", "Island", "SITE_TRT"))
}

tableS5 <- do.call(rbind, lapply(ordination_habitats, function(h) {
  parts <- species_subset(h)
  permanova_habitat(parts$spp, parts$meta, h)
}))

write.xlsx(tableS5, file.path(dir_tables, "TableS5_permanova_species.xlsx"),
           overwrite = TRUE)

## ---- Figure 4: species PCA per habitat -------------------------------
## Axis limits are fixed per panel, as in the submitted figure.
species_pca_limits <- list(
  BH = list(x = c(-0.30, 0.15), y = c(-0.20, 0.40), ann = "Year *\nTreatment *"),
  DH = list(x = c(-0.50, 0.30), y = c(-0.50, 0.50), ann = "Year *"),
  CH = list(x = c(-0.40, 0.40), y = c(-0.30, 0.40), ann = "Year *"),
  SB = list(x = c(-0.30, 0.40), y = c(-0.40, 0.40), ann = "Treatment *")
)

species_panels <- lapply(ordination_habitats, function(h) {
  parts <- species_subset(h)
  p     <- pca_scores(parts$spp, parts$meta)
  lim   <- species_pca_limits[[h]]
  pca_panel(p, title = h, annotation = lim$ann, xlim = lim$x, ylim = lim$y)
})

fig4 <- Reduce(`|`, species_panels) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 12, face = "plain"))

save_fig(fig4, "Figure4_PCA_species.png", width = 13, height = 4.5, dpi = 300)


## #####################################################################
##  8  RICHNESS AND SHANNON DIVERSITY             -> Figure 5
##
##  Final year of analysis in each region: 2019 in Iceland (2014 also
##  retained for the Icelandic habitats) and 2015 in Svalbard.
## #####################################################################

diversity_data <- subset(Moss_SPECIES, YEAR %in% c("2014","2015","2019"))

richnessmcmcglmm <- fit_or_load("richness", seeds[["richness"]],
  MCMCglmm(Richness ~ SITE_TRT - 1, rcov = ~ units,
           data   = diversity_data,
           prior  = prior_R,
           nitt   = 500000, burnin = 100000, thin = 400,
           family = "poisson", verbose = FALSE))

iter <- 400
shannonmcmcglmm <- fit_or_load("shannon", seeds[["shannon"]],
  MCMCglmm(Shannon ~ SITE_TRT - 1, rcov = ~ units,
           data   = diversity_data,
           prior  = prior_R,
           nitt   = 13000 * iter, burnin = 3000 * iter, thin = 10 * iter,
           verbose = FALSE))

#' Tidy a SITE_TRT cell-means model into a plotting frame
tidy_site_trt <- function(model) {
  d <- broom.mixed::tidy(model, conf.int = TRUE, conf.level = 0.95)
  d <- subset(d, effect == "fixed")
  parts <- do.call(rbind, strsplit(sub("SITE_TRT", "", d$term), "_"))
  d$location  <- factor(parts[, 1], levels = habitat_levels)
  d$treatment <- parts[, 2]
  d$Island    <- region_of(d$location)
  d
}

mossrichness_df <- tidy_site_trt(richnessmcmcglmm)
mossshannon_df  <- tidy_site_trt(shannonmcmcglmm)

diversity_panel <- function(df, ylab, ylims, backtransform = FALSE) {
  tf <- if (backtransform) exp else identity
  ggplot(df, aes(x = location, y = tf(estimate),
                 ymin = tf(conf.low), ymax = tf(conf.high), color = treatment)) +
    geom_pointrange(position = position_dodge(width = 0.3), size = 0.3) +
    xlab("") + ylab(ylab) +
    facet_nested_wrap(~Island, nrow = 1, scales = "free_x") +
    theme_pubr() +
    theme(strip.background = element_rect(color = "white", fill = "white", linetype = "blank"),
          strip.text = element_text(size = 12, colour = "black"),
          legend.position = "none") +
    scale_color_manual(values = trt_cols) +
    ylim(ylims)
}

fig5 <- cowplot::plot_grid(
  diversity_panel(mossrichness_df, "Richness", c(0, 12), backtransform = TRUE),
  diversity_panel(mossshannon_df,  "Shannon diversity", c(-0.2, 2)),
  ncol = 1, align = "hv", axis = "tblr", labels = c("a)", "b)"))

save_fig(fig5, "Figure5_richness_diversity.png", width = 5, height = 6)


## #####################################################################
##  9  BFG COMPOSITION: PERMANOVA AND PCA     -> Figure 6, Table S6
##
##  MH is excluded (see ordination_habitats in section 1).
## #####################################################################

bfg_subset <- function(h) {
  keep <- Moss_BFG$HABITAT == h
  parts <- split_community(Moss_BFG[keep, , drop = FALSE],
                           meta_cols = bfg_meta_cols,
                           drop_cols = "Lichen")
  parts
}

tableS6 <- do.call(rbind, lapply(ordination_habitats, function(h) {
  parts <- bfg_subset(h)
  permanova_habitat(parts$spp, parts$meta, h)
}))

write.xlsx(tableS6, file.path(dir_tables, "TableS6_permanova_bfg.xlsx"),
           overwrite = TRUE)

## ---- Figure 6: BFG PCA per habitat (MH is single-group) --------------
bfg_pca_limits <- list(
  BH = list(x = c(-0.30, 0.20), y = c(-0.30, 0.25), ann = "Year *\nTreatment *"),
  DH = list(x = c(-0.30, 0.40), y = c(-0.40, 0.40), ann = "Year *\nTreatment *"),
  CH = list(x = c(-0.40, 0.40), y = c(-0.30, 0.40), ann = "Year *"),
  SB = list(x = c(-0.50, 0.40), y = c(-0.30, 0.50), ann = NULL)
)

bfg_panels <- lapply(ordination_habitats, function(h) {
  parts <- bfg_subset(h)
  p     <- pca_scores(parts$spp, parts$meta)
  lim   <- bfg_pca_limits[[h]]
  pca_panel(p, title = h, annotation = lim$ann, xlim = lim$x, ylim = lim$y)
})

fig6 <- Reduce(`|`, bfg_panels) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 12, face = "plain"))

save_fig(fig6, "Figure6_PCA_bfg.png", width = 13, height = 4.5, dpi = 300)


## #####################################################################
## 10  BFG HIT COUNTS                   -> Figure 7, Figure S2
##
##  Figure 7 is the stacked mean cover of every cryptogam group,
##  including lichen and crust. Figure S2 models the four most abundant
##  BFGs. Habitat sets differ by group because not every group occurs
##  everywhere: BT in all five, We in all but MH, SU in Svalbard only,
##  TU in the snowbed only.
## #####################################################################

## ---- Figure 7: stacked bars with upper-only error bars ---------------
bfg_stack <- subset(Moss_BFG_long, YEAR %in% c("2015","2019")) %>%
  dplyr::group_by(HABITAT, TREATMENT, BFG) %>%
  dplyr::summarise(mean = mean(HITS, na.rm = TRUE),
                   se   = sd(HITS, na.rm = TRUE) / sqrt(sum(!is.na(HITS))),
                   .groups = "drop") %>%
  dplyr::mutate(HABITAT = factor(HABITAT, levels = habitat_levels),
                BFG     = factor(as.character(BFG)))

## Cumulative top of each segment within a habitat x treatment bar.
## ggplot stacks top-to-bottom in fill order, hence rev(cumsum(rev(.))).
bfg_stack <- bfg_stack %>%
  dplyr::arrange(HABITAT, TREATMENT, BFG) %>%
  dplyr::group_by(HABITAT, TREATMENT) %>%
  dplyr::mutate(ytop = rev(cumsum(rev(mean))),
                ymax = ytop + se) %>%
  dplyr::ungroup()

## geom_errorbar() cannot suppress its lower cap, so the upper-only bars
## are drawn as two segments. Both need inherit.aes = FALSE, or y is
## inherited from the column aesthetic and the segments land wrong.
fig7 <- ggplot(bfg_stack, aes(x = TREATMENT, y = mean, fill = BFG)) +
  geom_col(position = "stack", linewidth = 0.7, colour = NA) +
  geom_segment(aes(x = TREATMENT, xend = TREATMENT, y = ytop, yend = ymax),
               inherit.aes = FALSE, colour = "grey20", linewidth = 0.5) +
  geom_segment(aes(x = as.numeric(factor(TREATMENT)) - 0.1,
                   xend = as.numeric(factor(TREATMENT)) + 0.1,
                   y = ymax, yend = ymax),
               inherit.aes = FALSE, colour = "grey20", linewidth = 0.5) +
  facet_nested_wrap(~HABITAT, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = rev(bfg_palette)) +
  xlab("Treatment") + ylab("Cover (hits)") +
  theme_pubr() +
  theme(strip.background = element_rect(color = "white", fill = "white", linetype = "blank"),
        strip.text = element_text(size = 12, colour = "black"),
        legend.title = element_blank(),
        legend.position = "right")

save_fig(fig7, "Figure7_BFG_stacked.png", width = 8, height = 6, dpi = 600)

## ---- Figure S2: Poisson cell-means per BFG ---------------------------
iter <- 400

fit_bfg_hits <- function(grp, habs, seed) {
  fit_or_load(paste0("bfg_", grp), seed, {
    d <- subset(Moss_BFG_long,
                BFG == grp & HABITAT %in% habs & YEAR %in% c("2015","2019"))
    d$SITE_TRT   <- interaction(d$HABITAT, d$TREATMENT, sep = "_")
    d$PLOTEFFECT <- interaction(d$HABITAT, d$TREATMENT, d$PLOT, sep = "_")
    MCMCglmm(HITS ~ SITE_TRT - 1, random = ~ PLOTEFFECT, rcov = ~ units,
             data   = droplevels(d),
             prior  = prior_G1,
             nitt   = 13000 * iter, burnin = 3000 * iter, thin = 10 * iter,
             family = "poisson", verbose = FALSE)
  })
}

## TU occurs only in the snowbed and is sparse there. Check before
## trusting the model; if mostly zeros, report raw mean +- SE instead.
print(with(subset(Moss_BFG_long, BFG == "TU" & HABITAT == "SB" &
                    YEAR %in% c("2015","2019")),
           table(TREATMENT, HITS > 0)))

mod_BT <- fit_bfg_hits("BT", habitat_levels,            seeds[["bfg_BT"]])
mod_We <- fit_bfg_hits("We", c("BH","DH","CH","SB"),    seeds[["bfg_We"]])
mod_SU <- fit_bfg_hits("SU", c("DH","CH","SB"),         seeds[["bfg_SU"]])
mod_TU <- fit_bfg_hits("TU", "SB",                      seeds[["bfg_TU"]])

bfg_models <- list(BT = mod_BT, We = mod_We, SU = mod_SU, TU = mod_TU)

bfg_hits_df <- do.call(rbind, Map(function(model, grp) {
  d <- tidy_site_trt(model)
  d$BFG <- grp
  d
}, bfg_models, names(bfg_models)))

bfg_hits_df$BFG <- factor(bfg_hits_df$BFG, levels = c("BT","We","SU","TU"))

figS2 <- ggplot(bfg_hits_df,
                aes(x = location, y = exp(estimate),
                    ymin = exp(conf.low), ymax = exp(conf.high),
                    color = treatment)) +
  geom_pointrange(position = position_dodge(width = 0.3), size = 0.3) +
  facet_grid(BFG ~ Island, scales = "free", space = "free_x", drop = TRUE) +
  scale_color_manual(values = trt_cols) +
  xlab("") + ylab("Cover (hits)") +
  theme_pubr() +
  theme(strip.background = element_rect(color = "white", fill = "white", linetype = "blank"),
        strip.text = element_text(size = 11, colour = "black"),
        legend.position = "none")

save_fig(figS2, "FigureS2_BFG_posteriors.png", width = 8, height = 8)


## #####################################################################
## 11  FUNCTIONAL TRAITS              -> Figures 8, 9, 10
##
##  Shoot growth is a single joint Iceland + Svalbard model with bag
##  nested within plot. DH is absent from the growth and biomass data
##  because its bryophyte layer is too shallow to measure.
## #####################################################################

## ---- Shoot growth ----------------------------------------------------
prior_G2 <- list(G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002)),
                 R = list(V = 1, nu = 0.002))

mossgrowth_all <- fit_or_load("growth", seeds[["growth"]],
  MCMCglmm(Growth_mm ~ Loc_TRT - 1,
           random = ~ PLOT_ID + BAG_ID, rcov = ~ units,
           data   = Moss_growth,
           prior  = prior_G2,
           nitt   = 300000, burnin = 50000, thin = 250,
           family = "gaussian", verbose = FALSE))

tidy_loc_trt <- function(model) {
  d <- broom.mixed::tidy(model, conf.int = TRUE, conf.level = 0.95)
  d <- subset(d, effect == "fixed")
  stem <- sub("Loc_TRT", "", d$term)
  d$treatment <- ifelse(grepl("_CTL$", stem), "CTL", "OTC")
  d$location  <- factor(sub("_(CTL|OTC)$", "", stem), levels = habitat_levels)
  d$Island    <- region_of(d$location)
  d
}

mossgrowth_df <- tidy_loc_trt(mossgrowth_all)

trait_panel <- function(df, ylab, ylims, facet = TRUE) {
  g <- ggplot(df, aes(x = location, y = estimate,
                      ymin = conf.low, ymax = conf.high, color = treatment)) +
    geom_pointrange(position = position_dodge(width = 0.3), size = 0.3) +
    xlab("") + ylab(ylab) +
    scale_color_manual(values = trt_cols) +
    theme_pubr() +
    ylim(ylims) +
    theme(legend.position = "none")
  if (facet) {
    g <- g + facet_nested_wrap(~Island, nrow = 1, scales = "free_x") +
      theme(strip.background = element_rect(color = "white", fill = "white", linetype = "blank"),
            strip.text = element_text(size = 12, colour = "black"))
  }
  g
}

mossgrowthplot <- trait_panel(mossgrowth_df, "Shoot growth (mm)", c(0, 11))

## ---- Shoot biomass ---------------------------------------------------
## Fitted separately per region, as in the submitted analysis.
iter <- 200
mossbiom10Ice <- fit_or_load("biomass_iceland", seeds[["biomass_ice"]],
  MCMCglmm(Biomass_10shoots_g ~ Loc_TRT - 1, rcov = ~ units,
           data  = droplevels(subset(Moss_growth, HABITAT %in% c("BH","MH"))),
           prior = prior_R,
           nitt  = 13000 * iter, burnin = 3000 * iter, thin = 10 * iter,
           family = "gaussian", verbose = FALSE))

iter <- 100
mossbiom10Sval <- fit_or_load("biomass_svalbard", seeds[["biomass_sval"]],
  MCMCglmm(Biomass_10shoots_g ~ Loc_TRT - 1, rcov = ~ units,
           data  = droplevels(subset(Moss_growth, HABITAT %in% c("DH","CH","SB"))),
           prior = prior_R,
           nitt  = 13000 * iter, burnin = 3000 * iter, thin = 10 * iter,
           family = "gaussian", verbose = FALSE))

mossbiom_df <- rbind(tidy_loc_trt(mossbiom10Ice), tidy_loc_trt(mossbiom10Sval))
mossbiomplot <- trait_panel(mossbiom_df, "Shoot biomass (g)", c(0, 0.2))

fig8 <- cowplot::plot_grid(mossgrowthplot, mossbiomplot,
                           axis = "tblr", align = "hv",
                           labels = c("a)", "b)"), hjust = -3, vjust = 1)

save_fig(fig8, "Figure8_growth_biomass.png", width = 10, height = 5)

## ---- N2 fixation -----------------------------------------------------
iter <- 1000
mossnfixice <- fit_or_load("nfix_iceland", seeds[["nfix_iceland"]],
  MCMCglmm(Nfix_umol_g_day ~ SITE_TRT_Month - 1, random = ~ PLOTEFFECT, rcov = ~ units,
           data  = subset(N_fixation, REGION == "Iceland"),
           prior = prior_G1,
           nitt  = 13000 * iter, burnin = 3000 * iter, thin = 10 * iter,
           verbose = FALSE))

mossnfixsval <- fit_or_load("nfix_svalbard", seeds[["nfix_svalbard"]],
  MCMCglmm(Nfix_umol_g_day ~ SITE_TRT_Month - 1, random = ~ PLOTEFFECT, rcov = ~ units,
           data  = subset(N_fixation, REGION == "Svalbard"),
           prior = prior_G1,
           nitt  = 13000 * iter, burnin = 3000 * iter, thin = 10 * iter,
           verbose = FALSE))

tidy_nfix <- function(model, region) {
  d <- broom.mixed::tidy(model, conf.int = TRUE, conf.level = 0.95)
  d <- subset(d, effect == "fixed")
  parts <- strsplit(sub("SITE_TRT_Month", "", d$term), "_")
  d$location  <- factor(vapply(parts, `[`, "", 1), levels = habitat_levels)
  d$treatment <- vapply(parts, `[`, "", 2)
  d$MONTH     <- if (region == "Iceland") vapply(parts, `[`, "", 3) else ""
  d$Island    <- region
  ## Fixation rates cannot be negative; truncate the lower bound for
  ## plotting only.
  d$conf.low2 <- pmax(d$conf.low, 0)
  d
}

nfix_ice  <- tidy_nfix(mossnfixice,  "Iceland")
nfix_sval <- tidy_nfix(mossnfixsval, "Svalbard")

nfix_ice$MONTH <- factor(nfix_ice$MONTH, levels = c("JUNE", "JULY", "AUGUST"))

nfix_panel <- function(df, ylab, nested) {
  g <- ggplot(df, aes(x = location, y = estimate,
                      ymin = conf.low2, ymax = conf.high, color = treatment)) +
    geom_pointrange(position = position_dodge(width = 0.3), size = 0.3) +
    xlab("") + ylab(ylab) +
    theme_pubr() +
    theme(strip.background = element_rect(color = "white", fill = "white", linetype = "blank"),
          strip.text = element_text(size = 12, colour = "black"),
          legend.position = "none") +
    scale_color_manual(values = c(CTL = "#08306b", OTC = "#cb181d"))
  if (nested) {
    g + facet_nested(~ Island + MONTH, scales = "free", space = "free")
  } else {
    g + facet_nested(~ Island, scales = "free", space = "free")
  }
}

fig9 <- cowplot::plot_grid(
  nfix_panel(nfix_ice,  "N fixation (ethylene g-1 day-1)", nested = TRUE),
  nfix_panel(nfix_sval, "", nested = FALSE),
  ncol = 2, rel_widths = c(3.5, 2), labels = c("a)", "b)"))

save_fig(fig9, "Figure9_nitrogen_fixation.png", width = 7.5, height = 5)

## ---- Water holding capacity and colony density -----------------------
WHC <- fit_or_load("whc", seeds[["whc"]],
  MCMCglmm(WHC_g_g ~ Loc_TRT - 1, random = ~ COLEFFECT, rcov = ~ units,
           data  = Moss_WHC,
           prior = prior_G1,
           nitt  = 200000, burnin = 40000, thin = 160,
           family = "gaussian", verbose = FALSE))

## Colony density is on the order of 0.05 g cm-3, small enough that the
## default inverse-Wishart prior with V = 1 inflates the variance. A
## parameter-expanded prior is used instead.
prior_PX <- list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)),
                 R = list(V = 1, nu = 0.002))

Biomcol <- fit_or_load("colony_density", seeds[["colony_density"]],
  MCMCglmm(Colony_density_g_cm3 ~ Loc_TRT - 1,
           random = ~ Loc_TRT_plot, rcov = ~ units,
           data   = Moss_WHC,
           prior  = prior_PX,
           nitt   = 250000, burnin = 40000, thin = 120,
           family = "gaussian", verbose = FALSE))

mosswhcplot  <- trait_panel(tidy_loc_trt(WHC),
                            "Water holding capacity (g/g)", c(0, 8), facet = FALSE)
Biomcolplot  <- trait_panel(tidy_loc_trt(Biomcol),
                            "Colony density (g cm-3)", c(0, 0.08), facet = FALSE)

fig10 <- cowplot::plot_grid(mosswhcplot, Biomcolplot,
                            align = "hv", axis = "tblr",
                            labels = c("a)", "b)"))

save_fig(fig10, "Figure10_WHC_colony_density.png", width = 10, height = 5)


## #####################################################################
## 12  RESULTS WORKBOOK
##
##  Per-treatment posterior intervals and the OTC - CTL contrast for
##  every model, one sheet each. Rows where the credible intervals do
##  not overlap are highlighted.
## #####################################################################

model_specs <- list(
  list(obj = mossdepth,         fam = "gaussian", sheet = "Layer depth"),
  list(obj = mossgrowth_all,    fam = "gaussian", sheet = "Shoot growth"),
  list(obj = mossbiom10Ice,     fam = "gaussian", sheet = "Biomass Iceland"),
  list(obj = mossbiom10Sval,    fam = "gaussian", sheet = "Biomass Svalbard"),
  list(obj = WHC,               fam = "gaussian", sheet = "WHC"),
  list(obj = Biomcol,           fam = "gaussian", sheet = "Colony density"),
  list(obj = shannonmcmcglmm,   fam = "gaussian", sheet = "Shannon"),
  list(obj = mossnfixice,       fam = "gaussian", sheet = "Nfix Iceland"),
  list(obj = mossnfixsval,      fam = "gaussian", sheet = "Nfix Svalbard"),
  list(obj = mosscovermcmcglmm, fam = "poisson",  sheet = "Cover"),
  list(obj = richnessmcmcglmm,  fam = "poisson",  sheet = "Richness"),
  list(obj = mod_BT,            fam = "poisson",  sheet = "BFG BT"),
  list(obj = mod_We,            fam = "poisson",  sheet = "BFG We"),
  list(obj = mod_SU,            fam = "poisson",  sheet = "BFG SU"),
  list(obj = mod_TU,            fam = "poisson",  sheet = "BFG TU")
)

wb  <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "bottom")
hln <- createStyle(fgFill = "#FCE4D6")
num <- createStyle(numFmt = "0.000")

for (spec in model_specs) {
  tbl <- tryCatch(mcmc_diff(spec$obj, family = spec$fam, model_name = spec$sheet),
                  error = function(e) {
                    message("FAILED ", spec$sheet, ": ", conditionMessage(e)); NULL
                  })
  if (is.null(tbl)) next

  addWorksheet(wb, spec$sheet)
  writeData(wb, spec$sheet, tbl, headerStyle = hdr)
  addStyle(wb, spec$sheet, num, rows = 2:(nrow(tbl) + 1), cols = c(2:7, 9:11),
           gridExpand = TRUE, stack = TRUE)

  nr <- which(tbl$CIs_overlap == "no") + 1
  if (length(nr)) {
    addStyle(wb, spec$sheet, hln, rows = nr, cols = 1:ncol(tbl),
             gridExpand = TRUE, stack = TRUE)
  }
  setColWidths(wb, spec$sheet, cols = 1:ncol(tbl), widths = "auto")
}

saveWorkbook(wb, file.path(dir_tables, "moss_treatment_effects.xlsx"), overwrite = TRUE)


## ---------------------------------------------------------------------
##  Session information
## ---------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), here("sessionInfo.txt"))
message("Done. Figures in output/figures, tables in output/tables.")
