# Synthetic benchmark for MetaWeave
#
# This script exercises the package's current public implementation:
#   * aligned terra::SpatRaster distribution groups;
#   * block_model();
#   * probability_matrix_model();
#   * a context-dependent backend built with new_inference_model();
#   * run_spatial_inference(), summarize_networks(), and map_metrics().

set.seed(20260809)

if (!requireNamespace("terra", quietly = TRUE)) {
  stop("The 'terra' package is required.")
}

library(metaweave)


inv_logit <- function(x) 1 / (1 + exp(-x))
clamp_probability <- function(x, epsilon = 1e-8) {
  pmin(pmax(x, epsilon), 1 - epsilon)
}

output_dir <- file.path(getwd(), "synthetic-benchmark-results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Generate a landscape and two bipartite species groups
grid_size <- 15L
n_producers <- 24L
n_consumers <- 18L
n_producer_blocks <- 3L
n_consumer_blocks <- 3L

template <- terra::rast(
  nrows = grid_size,
  ncols = grid_size,
  xmin = 0,
  xmax = 1,
  ymin = 0,
  ymax = 1,
  crs = "EPSG:4326"
)

xy <- terra::xyFromCell(template, seq_len(terra::ncell(template)))
environment <- as.numeric(scale(xy[, 1], center = TRUE, scale = TRUE))

producer_names <- sprintf("producer_%02d", seq_len(n_producers))
consumer_names <- sprintf("consumer_%02d", seq_len(n_consumers))

producer_lookup <- data.frame(
  species = producer_names,
  guild = paste0(
    "P",
    rep(seq_len(n_producer_blocks), length.out = n_producers)
  )
)

consumer_lookup <- data.frame(
  species = consumer_names,
  guild = paste0(
    "C",
    rep(seq_len(n_consumer_blocks), length.out = n_consumers)
  )
)

producer_traits <- data.frame(
  species = producer_names,
  optimum = runif(n_producers, -1.5, 1.5),
  breadth = runif(n_producers, 0.65, 1.1),
  intercept = rnorm(n_producers, 1.1, 0.25)
)

consumer_traits <- data.frame(
  species = consumer_names,
  optimum = runif(n_consumers, -1.5, 1.5),
  breadth = runif(n_consumers, 0.65, 1.1),
  intercept = rnorm(n_consumers, 0.9, 0.25)
)

simulate_distribution_stack <- function(template, environment, traits) {
  layers <- lapply(seq_len(nrow(traits)), function(i) {
    probability <- inv_logit(
      traits$intercept[i] -
        (environment - traits$optimum[i])^2 /
          (2 * traits$breadth[i]^2)
    )
    presence <- stats::rbinom(length(probability), 1L, probability)
    raster <- terra::rast(template)
    terra::values(raster) <- ifelse(presence == 1L, 1, NA_real_)
    names(raster) <- traits$species[i]
    raster
  })
  terra::rast(layers)
}

producer_stack <- simulate_distribution_stack(
  template,
  environment,
  producer_traits
)

consumer_stack <- simulate_distribution_stack(
  template,
  environment,
  consumer_traits
)

distributions <- list(
  producers = new_distribution_collection(
    producer_stack,
    guild = "producers",
    metadata = producer_traits
  ),
  consumers = new_distribution_collection(
    consumer_stack,
    guild = "consumers",
    metadata = consumer_traits
  )
)

# Generate regional interaction structure
theta <- matrix(
  c(
    0.70, 0.22, 0.08,
    0.15, 0.68, 0.25,
    0.10, 0.28, 0.72
  ),
  nrow = n_producer_blocks,
  byrow = TRUE,
  dimnames = list(
    paste0("P", seq_len(n_producer_blocks)),
    paste0("C", seq_len(n_consumer_blocks))
  )
)

producer_block <- producer_lookup$guild[
  match(producer_names, producer_lookup$species)
]
consumer_block <- consumer_lookup$guild[
  match(consumer_names, consumer_lookup$species)
]

sbm_probability <- theta[producer_block, consumer_block, drop = FALSE]
dimnames(sbm_probability) <- list(producer_names, consumer_names)

# Add species-pair heterogeneity around the block-level expectation. This is
# the regional probability surface used by probability_matrix_model().
pair_effect <- matrix(
  stats::rnorm(n_producers * n_consumers, 0, 0.55),
  nrow = n_producers,
  ncol = n_consumers,
  dimnames = list(producer_names, consumer_names)
)

regional_probability <- inv_logit(
  stats::qlogis(clamp_probability(sbm_probability)) + pair_effect
)

dimnames(regional_probability) <- list(producer_names, consumer_names)

# Each pair has an environmental optimum. Local interaction probabilities
# depart from the regional matrix as environmental context changes.
pair_optimum <- matrix(
  stats::runif(n_producers * n_consumers, -1.5, 1.5),
  nrow = n_producers,
  ncol = n_consumers,
  dimnames = list(producer_names, consumer_names)
)

context_probability <- function(rows, columns, site_environment) {
  base <- regional_probability[rows, columns, drop = FALSE]
  optimum <- pair_optimum[rows, columns, drop = FALSE]
  context_match <- exp(-0.5 * ((site_environment - optimum) / 0.8)^2)
  inv_logit(
    stats::qlogis(clamp_probability(base)) - 0.9 + 1.8 * context_match
  )
}

# Construct models through MetaWeave's public API

sbm <- block_model(
  row_lookup = producer_lookup,
  column_lookup = consumer_lookup,
  theta = theta,
  row_group = "producers",
  column_group = "consumers"
)

regional_matrix <- probability_matrix_model(
  probability_matrix = regional_probability,
  row_group = "producers",
  column_group = "consumers",
  name = "regional species-pair model"
)

context_predict <- function(model, assemblage) {
  par <- model$parameters
  rows <- assemblage$species[[par$row_group]]
  columns <- assemblage$species[[par$column_group]]

  rows <- rows[rows %in% rownames(par$regional_probability)]
  columns <- columns[columns %in% colnames(par$regional_probability)]

coordinates <- as.numeric(assemblage$coordinates)

if (length(coordinates) < 1L || !is.finite(coordinates[[1L]])) {
  stop("The context-dependent model requires a finite x coordinate.")
}

site_environment <- as.numeric(scale(
  coordinates[[1L]],
  center = par$x_center,
  scale = par$x_scale
))

  base <- par$regional_probability[rows, columns, drop = FALSE]
  optimum <- par$pair_optimum[rows, columns, drop = FALSE]
  context_match <- exp(
    -0.5 * ((site_environment - optimum) / par$context_width)^2
  )
  probability <- inv_logit(
    stats::qlogis(clamp_probability(base)) +
      par$context_intercept +
      par$context_strength * context_match
  )

  new_ecological_network(
    probability,
    type = "probability",
    metadata = list(
      cell_id = assemblage$cell_id,
      model = model$name,
      environment = site_environment
    ),
    row_group = par$row_group,
    column_group = par$column_group
  )
}

context_model <- new_inference_model(
  predict = context_predict,
  name = "context-dependent species-pair model",
  parameters = list(
    regional_probability = regional_probability,
    pair_optimum = pair_optimum,
    row_group = "producers",
    column_group = "consumers",
    x_center = mean(xy[, 1]),
    x_scale = stats::sd(xy[, 1]),
    context_width = 0.8,
    context_intercept = -0.9,
    context_strength = 1.8
  ),
  subclass = "context_probability_model"
)

models <- list(
  SBM = sbm,
  regional_matrix = regional_matrix,
  context_model = context_model
)


# Run MetaWeave end to end


results <- lapply(models, function(model) {
  run_spatial_inference(
    distributions = distributions,
    model = model,
    min_species = 1L,
    threshold = NULL
  )
})

# All models must operate on the same localized assemblages.
reference_index <- results[[1L]]$result$index[, c("cell_id", "x", "y")]
for (i in seq_along(results)[-1L]) {
  stopifnot(identical(
    reference_index,
    results[[i]]$result$index[, c("cell_id", "x", "y")]
  ))
}

# Generate known local truth and evaluate probability recovery

evaluate_one_network <- function(network, site_environment, seed) {
  rows <- rownames(network$matrix)
  columns <- colnames(network$matrix)
  truth <- context_probability(rows, columns, site_environment)

  set.seed(seed)
  observed <- matrix(
    stats::rbinom(length(truth), 1L, as.numeric(truth)),
    nrow = nrow(truth),
    ncol = ncol(truth),
    dimnames = dimnames(truth)
  )

  predicted <- network$matrix
  predicted_safe <- clamp_probability(predicted)

  data.frame(
    probability_rmse = sqrt(mean((predicted - truth)^2)),
    probability_mae = mean(abs(predicted - truth)),
    brier = mean((predicted - observed)^2),
    log_loss = -mean(
      observed * log(predicted_safe) +
        (1 - observed) * log(1 - predicted_safe)
    ),
    true_expected_links = sum(truth),
    predicted_expected_links = sum(predicted),
    observed_links = sum(observed),
    true_connectance = mean(truth),
    predicted_connectance = mean(predicted)
  )
}

site_results <- do.call(
  rbind,
  lapply(names(results), function(model_name) {
    collection <- results[[model_name]]$result
    rows <- lapply(seq_along(collection$networks), function(i) {
      network <- collection$networks[[i]]
      site_environment <- as.numeric(scale(
        collection$index$x[i],
        center = mean(xy[, 1]),
        scale = stats::sd(xy[, 1])
      ))
      metrics <- evaluate_one_network(
        network,
        site_environment,
        seed = 100000L + collection$index$cell_id[i]
      )
      cbind(
        data.frame(
          model = model_name,
          cell_id = collection$index$cell_id[i],
          x = collection$index$x[i],
          y = collection$index$y[i],
          environment = site_environment
        ),
        metrics
      )
    })
    do.call(rbind, rows)
  })
)

rmse <- function(observed, predicted) {
  sqrt(mean((observed - predicted)^2, na.rm = TRUE))
}

model_summary <- do.call(
  rbind,
  lapply(split(site_results, site_results$model), function(x) {
    data.frame(
      model = unique(x$model),
      probability_rmse = mean(x$probability_rmse),
      probability_mae = mean(x$probability_mae),
      mean_brier = mean(x$brier),
      mean_log_loss = mean(x$log_loss),
      expected_links_rmse = rmse(
        x$true_expected_links,
        x$predicted_expected_links
      ),
      connectance_rmse = rmse(
        x$true_connectance,
        x$predicted_connectance
      )
    )
  })
)

rownames(model_summary) <- NULL
model_summary <- model_summary[
  order(model_summary$probability_rmse),
]

# The custom context backend is the data-generating model and should recover
# the known probabilities up to numerical precision.
context_error <- model_summary$probability_rmse[
  model_summary$model == "context_model"
]
stopifnot(length(context_error) == 1L, context_error < 1e-10)


# Save reproducible outputs
utils::write.csv(
  site_results,
  file.path(output_dir, "site_level_performance.csv"),
  row.names = FALSE
)

utils::write.csv(
  model_summary,
  file.path(output_dir, "model_summary.csv"),
  row.names = FALSE
)

terra::writeRaster(
  results$context_model$spatial$data,
  file.path(output_dir, "context_model_metrics.tif"),
  overwrite = TRUE
)


model_colours <- c(
  SBM = "#D55E00",
  regional_matrix = "#0072B2",
  context_model = "#009E73"
)

grDevices::png(
  file.path(output_dir, "benchmark_diagnostics.png"),
  width = 1800,
  height = 850,
  res = 150
)

graphics::par(mfrow = c(1, 2), mar = c(7, 4, 2, 1))
graphics::boxplot(
  probability_rmse ~ model,
  data = site_results,
  col = model_colours[sort(unique(site_results$model))],
  las = 2,
  xlab = "",
  ylab = "Site-level probability RMSE"
)

graphics::plot(
  NA,
  xlim = range(site_results$true_connectance),
  ylim = range(site_results$predicted_connectance),
  xlab = "True expected connectance",
  ylab = "Predicted expected connectance"
)
graphics::abline(0, 1, lty = 2, col = "grey50")

for (model_name in names(results)) {
  keep <- site_results$model == model_name
  graphics::points(
    site_results$true_connectance[keep],
    site_results$predicted_connectance[keep],
    pch = 16,
    col = grDevices::adjustcolor(
      model_colours[model_name],
      alpha.f = 0.45
    )
  )
}

graphics::legend(
  "topleft",
  legend = names(results),
  col = model_colours[names(results)],
  pch = 16,
  bty = "n"
)

grDevices::dev.off()

saveRDS(
  list(
    settings = list(
      seed = 20260809,
      grid_size = grid_size,
      n_producers = n_producers,
      n_consumers = n_consumers
    ),
    distributions = distributions,
    producer_lookup = producer_lookup,
    consumer_lookup = consumer_lookup,
    theta = theta,
    regional_probability = regional_probability,
    pair_optimum = pair_optimum,
    results = results,
    site_results = site_results,
    model_summary = model_summary
  ),
  file.path(output_dir, "benchmark_objects.rds")
)

print(model_summary)
message("Benchmark outputs written to: ", normalizePath(output_dir))


