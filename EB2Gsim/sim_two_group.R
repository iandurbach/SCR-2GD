# Two-group extension of Efford & Boulanger (2019)
# - Fix spacing at 2 * sigma1
# - Vary grid size (nx, ny) as in the original study
# - Vary sigma2 while holding a0_2 constant (lambda0_2 adjusts with sigma2)
# - Fit two-group models (D ~ g, lambda0 ~ g, sigma ~ g) and summarise RSE

library(secr)
library(secrdesign)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(ggplot2)
library(foreach)
library(doParallel)

lambda_from_a0 <- function(a0, sigma) {
  a0 / (2 * pi * sigma^2)
}

build_two_group_plan <- function(sigma1 = 20,
                                 sigma2 = c(2, 5, 10, 20, 40, 60, 100),
                                 grids = list(c(6, 6), c(8, 8), c(10, 10)),
                                 a0 = 400,
                                 D1 = 10,
                                 D2 = 10,
                                 nocc = 5,
                                 nrepeats = 1,
                                 detector = "proximity") {
  spacing <- 2 * sigma1
  lambda01 <- lambda_from_a0(a0, sigma1)

  trapset <- lapply(grids, function(g) {
    make.grid(nx = g[1], ny = g[2], spacing = spacing, detector = detector)
  })
  names(trapset) <- sprintf("grid_%sx%s", vapply(grids, `[`, numeric(1), 1), vapply(grids, `[`, numeric(1), 2))

  mask_buffer <- 4 * max(c(sigma1, sigma2))
  mask_spacing <- sigma1 / 2
  maskset <- lapply(trapset, make.mask, spacing = mask_spacing, buffer = mask_buffer, type = "trapbuffer")

  scen_list <- list()
  scenario_id <- 1

  for (s2 in sigma2) {
    lambda02 <- lambda_from_a0(a0, s2)
    for (g_idx in seq_along(trapset)) {
      base <- make.scenarios(
        D = D1,
        detectfn = "HHN",
        lambda0 = lambda01,
        sigma = sigma1,
        trapsindex = g_idx,
        noccasions = nocc,
        nrepeats = nrepeats,
        groups = c("g1", "g2")
      )
      is_g2 <- base$group == "g2"
      base$D[is_g2] <- D2
      base$lambda0[is_g2] <- lambda02
      base$sigma[is_g2] <- s2

      base$scenario <- scenario_id
      base$scenario_label <- sprintf("%s_sigma2_%s", names(trapset)[g_idx], s2)
      base$sigma1 <- sigma1
      base$sigma2 <- s2
      base$lambda01 <- lambda01
      base$lambda02 <- lambda02
      base$a0 <- a0
      base$spacing <- spacing
      base$grid_label <- names(trapset)[g_idx]

      scen_list[[length(scen_list) + 1]] <- base
      scenario_id <- scenario_id + 1
    }
  }

  scenarios <- bind_rows(scen_list)
  list(scenarios = scenarios, trapset = trapset, maskset = maskset)
}

tidy_pred <- function(pred_list, scenario_label, sigma1, sigma2, grid_label, repl_id) {
  if (is.null(pred_list)) return(data.frame())
  map_dfr(seq_along(pred_list), function(idx) {
    nm <- names(pred_list)[idx]
    grp <- str_match(nm, "g\\s*=\\s*(\\S+)")[, 2]
    pars <- c("D", "lambda0", "sigma")
    df <- pred_list[[idx]]
    df$parameter <- pars
    df$group <- grp
    df$scenario_label <- scenario_label
    df$grid_label <- grid_label
    df$sigma1 <- sigma1
    df$sigma2 <- sigma2
    df$repl <- repl_id
    df
  })
}

count_groups <- function(capthist) {
  if (isTRUE(getOption("EB2Gsim.debug", FALSE))) {
    message("count_groups debug: is.null=", is.null(capthist), " is.array=", is.array(capthist))
  }
  if (is.null(capthist) || !is.array(capthist)) {
    return(data.frame(group = NA_character_, n = NA_real_, r = NA_real_, r_spatial = NA_real_))
  }
  if (dim(capthist)[1] == 0) {
    return(data.frame(group = NA_character_, n = NA_real_, r = NA_real_, r_spatial = NA_real_))
  }
  detections <- rowSums(apply(capthist, c(1, 2), sum))
  det_by_detector <- apply(capthist, c(1, 3), sum)
  spatial_recaps_ind <- pmax(rowSums(det_by_detector > 0) - 1, 0)
  grp <- covariates(capthist)$group
  if (is.null(grp)) grp <- rep("all", length(detections))
  counts_total <- aggregate(detections, list(group = grp), function(x) c(n = sum(x > 0), r = sum(pmax(x - 1, 0))))
  counts_spatial <- aggregate(spatial_recaps_ind, list(group = grp), function(x) c(r_spatial = sum(x)))
  if (isTRUE(getOption("EB2Gsim.debug", FALSE))) {
    message("count_groups: head grp=", paste(head(grp), collapse = ","), " detections length=", length(detections))
    message("count_groups total matrix: ", paste(capture.output(print(counts_total)), collapse = " | "))
    message("count_groups spatial matrix: ", paste(capture.output(print(counts_spatial)), collapse = " | "))
  }
  data.frame(
    group = counts_total$group,
    n = counts_total$x[, "n"],
    r = counts_total$x[, "r"],
    r_spatial = counts_spatial$x[, "r_spatial"]
  )
}

extract_pred_and_counts <- function(fit) {
  if (isTRUE(getOption("EB2Gsim.debug", FALSE))) {
    message("extract_pred_and_counts called with class: ", paste(class(fit), collapse = ", "),
            " capthist NULL? ", is.null(fit$capthist))
  }
  if (inherits(fit, "try-error") || is.null(fit$capthist)) {
    return(list(pred = NULL, counts = data.frame(group = NA_character_, n = NA_real_, r = NA_real_, approx_rse = NA_real_)))
  }
  counts <- try(count_groups(fit$capthist), silent = TRUE)
  if (inherits(counts, "try-error")) {
    counts <- data.frame(group = NA_character_, n = NA_real_, r = NA_real_)
  }
  pred <- try(predict(fit), silent = TRUE)
  if (inherits(pred, "try-error")) pred <- NULL
  counts$approx_rse <- 1 / pmin(sqrt(counts$n), sqrt(counts$r))
  list(pred = pred, counts = counts)
}

run_two_group <- function(plan,
                          nrepl = 200,
                          seed = 123,
                          fit_model = list(detectfn = "HHN", model = list(D ~ g, lambda0 ~ g, sigma ~ g), groups = "group"),
                          ncores = NULL) {
  scenarios_split <- split(plan$scenarios, plan$scenarios$scenario)
  if (!is.null(ncores) && ncores > 1) {
    cl <- makeCluster(ncores)
    on.exit(stopCluster(cl), add = TRUE)
    registerDoParallel(cl)
  }

  sims_list <- foreach(
    scen = scenarios_split,
    .packages = c("secr", "secrdesign"),
    .export = c("extract_pred_and_counts", "count_groups", "tidy_pred")
  ) %dopar% {
    run.scenarios(
      nrepl = nrepl,
      scenarios = scen,
      trapset = plan$trapset,
      maskset = plan$maskset,
      fit = TRUE,
      fit.args = fit_model,
      extractfn = extract_pred_and_counts,
      seed = seed,
      multisession = FALSE,
      ncores = NULL
    )
  }

  # stitch back into one "sims" structure
  sims <- list(
    output = lapply(sims_list, function(x) x$output[[1]]),
    scenarios = plan$scenarios,
    trapset = plan$trapset,
    maskset = plan$maskset
  )

  scenario_meta <- plan$scenarios %>%
    distinct(scenario, scenario_label, sigma2, grid_label, spacing, a0, lambda01, lambda02, sigma1)
  scenario_group_meta <- plan$scenarios %>%
    select(scenario, scenario_label, grid_label, sigma1, sigma2, group, D_true = D, lambda0_true = lambda0, sigma_true = sigma) %>%
    distinct()

  per_repl_rows <- list()

  tidy_list <- map2(
    sims$output,
    split(scenario_meta, scenario_meta$scenario),
    function(scenario_outputs, meta_row) {
      preds <- list()
      counts <- list()
      meta_groups <- scenario_group_meta %>% filter(scenario_label == meta_row$scenario_label)
      for (repl_id in seq_along(scenario_outputs)) {
        out <- scenario_outputs[[repl_id]]
        pred_df <- tidy_pred(out$pred, meta_row$scenario_label, meta_row$sigma1, meta_row$sigma2, meta_row$grid_label, repl_id)
        cnt_df <- mutate(out$counts,
          scenario_label = meta_row$scenario_label,
          sigma2 = meta_row$sigma2,
          grid_label = meta_row$grid_label,
          repl = repl_id,
          spacing = meta_row$spacing,
          sigma1 = meta_row$sigma1
        )
        preds[[repl_id]] <- pred_df
        counts[[repl_id]] <- cnt_df

        if (nrow(cnt_df) > 0 && nrow(pred_df) > 0) {
          pred_wide <- pred_df %>%
            select(group, parameter, estimate, SE.estimate, repl) %>%
            pivot_wider(names_from = parameter, values_from = c(estimate, SE.estimate), names_sep = "_")
          combined <- cnt_df %>%
            left_join(meta_groups, by = c("scenario_label", "grid_label", "sigma1", "sigma2", "group")) %>%
            left_join(pred_wide, by = c("group", "repl"))
          combined <- combined %>%
            mutate(
              approx_rse = ifelse(!is.na(n) & !is.na(r) & n > 0 & r > 0, 1 / pmin(sqrt(n), sqrt(r)), NA_real_),
              rse_sim = SE.estimate_D / estimate_D
            )
          per_repl_rows[[length(per_repl_rows) + 1]] <<- combined
        }
      }
      list(pred = bind_rows(preds), counts = bind_rows(counts))
    }
  )

  per_repl_df <- if (length(per_repl_rows)) bind_rows(per_repl_rows) else tibble()
  if (nrow(per_repl_df) && !"r_spatial" %in% names(per_repl_df)) per_repl_df$r_spatial <- NA_real_
  tidy_pred_df <- bind_rows(map(tidy_list, "pred"))
  tidy_counts_df <- bind_rows(map(tidy_list, "counts"))

  empirical_rse <- if (nrow(per_repl_df) == 0) {
    tibble()
  } else {
    per_repl_df %>%
      group_by(scenario_label, grid_label, sigma1, sigma2, group) %>%
      summarise(
        mean_rse = mean(rse_sim, na.rm = TRUE),
        sd_rse = sd(rse_sim, na.rm = TRUE),
        n_rep = n(),
        .groups = "drop"
      )
  }

  approx_counts <- if (nrow(per_repl_df) == 0) {
    tibble()
  } else {
    per_repl_df %>%
      filter(!is.na(n)) %>%
      group_by(scenario_label, grid_label, sigma1, sigma2, group) %>%
      summarise(
        mean_n = mean(n, na.rm = TRUE),
        mean_r = mean(r, na.rm = TRUE),
        mean_r_spatial = mean(r_spatial, na.rm = TRUE),
        approx_rse = mean(approx_rse, na.rm = TRUE),
        .groups = "drop"
      )
  }

  list(
    sims = sims,
    predictions = tidy_pred_df,
    counts = tidy_counts_df,
    per_repl = per_repl_df,
    empirical_rse = empirical_rse,
    approx_counts = approx_counts
  )
}

plot_two_group_rse <- function(results, outfile = NULL) {
  combined <- results$empirical_rse %>%
    left_join(results$approx_counts, by = c("scenario_label", "grid_label", "sigma2", "group")) %>%
    mutate(sigma_ratio = sigma2 / sigma1)

  if (nrow(combined) == 0) stop("No data available to plot")

  plt <- ggplot(combined, aes(x = sigma2, color = group)) +
    geom_line(aes(y = mean_rse), size = 1) +
    geom_point(aes(y = mean_rse), size = 2) +
    geom_line(aes(y = approx_rse), linetype = "dashed") +
    geom_point(aes(y = approx_rse), shape = 1) +
    facet_wrap(~grid_label, scales = "free_y") +
    labs(
      x = expression(sigma[2]),
      y = "RSE(D-hat)",
      color = "Group",
      title = "Empirical vs approximate RSE(D-hat) for two-group simulations",
      subtitle = "Solid lines: empirical from simulations; dashed: 1/min(sqrt(n), sqrt(r))"
    ) +
    theme_minimal()

  if (!is.null(outfile)) {
    ggsave(outfile, plt, width = 10, height = 6, dpi = 300)
  }
  plt
}
