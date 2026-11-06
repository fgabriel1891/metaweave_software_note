# Fit one bipartite SBM to the union of all six Ecuador networks.
# Run with Rscript before 01_run.R. Requires cassandRa.

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(arg) != 1L) stop("Run this file with Rscript.")
here <- dirname(normalizePath(sub("^--file=", "", arg), mustWork = TRUE))
if (!requireNamespace("cassandRa", quietly = TRUE)) {
  stop("Install cassandRa before fitting the pooled SBM.")
}
input <- file.path(dirname(here), "data", "ecuador", "dataSet1.RData")
if (!file.exists(input)) stop("Missing extracted Zenodo data: ", input)
env <- new.env(parent = emptyenv())
load(input, envir = env)
dat <- env$dataSet1$N
needed <- c("site", "plantCode", "animalCode", "frequency")
if (!is.data.frame(dat) || !all(needed %in% names(dat))) {
  stop("The archived interaction table has unexpected columns.")
}
dat <- dat[!is.na(dat$frequency) & dat$frequency > 0, ]
if (!nrow(dat) || anyNA(dat[c("site", "plantCode", "animalCode")])) {
  stop("Invalid interaction records.")
}
p_pool <- sort(unique(as.character(dat$plantCode)))
a_pool <- sort(unique(as.character(dat$animalCode)))
sites <- sort(unique(as.character(dat$site)))
if (length(sites) != 6L) stop("Expected six Ecuador sites.")

# A regional link is 1 when recorded at one or more sites, regardless of
# frequency or the number of sites where it was recorded.
metaweb <- matrix(0L, length(p_pool), length(a_pool),
                  dimnames = list(p_pool, a_pool))
metaweb[cbind(match(as.character(dat$plantCode), p_pool),
              match(as.character(dat$animalCode), a_pool))] <- 1L
stopifnot(all(metaweb %in% 0:1), sum(metaweb) > 0L)

seed <- 5537283L
set.seed(seed)
cass_network <- cassandRa::CreateListObject(metaweb)
fit <- cassandRa::FitSBM(list = cass_network, n_SBM = 10, G = NULL)
# The existing palm-mammal workflow uses SBM1 as the selected best fit.
if (is.null(fit$SBM1)) stop("cassandRa did not return SBM1; inspect its API.")
best <- fit$SBM1
theta <- best$Omega_rs
if (!is.matrix(theta) || !is.numeric(theta) || !length(theta) ||
    any(!is.finite(theta)) || any(theta < 0 | theta > 1)) {
  stop("The fitted block probabilities are invalid.")
}
plant_block <- as.integer(best$SB_H)
animal_block <- as.integer(best$SB_W)
if (length(plant_block) != nrow(metaweb) ||
    length(animal_block) != ncol(metaweb) ||
    anyNA(plant_block) || anyNA(animal_block) ||
    any(!plant_block %in% seq_len(nrow(theta))) ||
    any(!animal_block %in% seq_len(ncol(theta)))) {
  stop("Fitted block assignments do not cover the complete regional metaweb.")
}
rownames(theta) <- paste0("plant_block_", seq_len(nrow(theta)))
colnames(theta) <- paste0("animal_block_", seq_len(ncol(theta)))
row_lookup <- data.frame(species = rownames(metaweb),
                         guild = rownames(theta)[plant_block])
column_lookup <- data.frame(species = colnames(metaweb),
                            guild = colnames(theta)[animal_block])
stopifnot(!anyNA(row_lookup$guild), !anyNA(column_lookup$guild),
          identical(row_lookup$species, rownames(metaweb)),
          identical(column_lookup$species, colnames(metaweb)))

out <- file.path(dirname(here), "results")
dir.create(out, showWarnings = FALSE)
saveRDS(list(version = 1L, source_doi = "10.5281/zenodo.5537283",
             input_md5 = unname(tools::md5sum(input)),
             metaweb = metaweb, seed = seed, starts = 10L,
             block_selection = "G = NULL (cassandRa default)",
             selected_fit = "SBM1", fit_fields = names(best),
             cassandra_version = as.character(utils::packageVersion("cassandRa")),
             row_lookup = row_lookup, column_lookup = column_lookup,
             theta = theta),
        file.path(out, "pooled_sbm_fit.rds"))
message("Saved fitted pooled SBM: ", file.path(out, "pooled_sbm_fit.rds"))
message("Blocks: ", nrow(theta), " plant, ", ncol(theta), " animal")
