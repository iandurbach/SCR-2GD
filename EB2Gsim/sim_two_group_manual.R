# Two-group simulation using explicit sim.popn -> sim.capthist -> secr.fit steps
# - Constant D, no spatial covariates
# - Group-specific lambda0/sigma (via a0 = 2 * pi * lambda0 * sigma^2)
# - Trap spacing fixed at 2 * sigma1; grid size varies
# - Parallelised across scenarios with foreach/doParallel (secr itself runs serially)

library(secr)
library(secrdesign)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(foreach)
library(doParallel)

lambda_from_a0 <- function(a0, sigma, nocc = 1) {
  a0 / (2 * pi * sigma^2 * nocc)
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

build_two_group_plan_manual <- function(sigma1 = 20,
                                        sigma2 = c(2, 5, 10, 20, 40, 60, 100),
                                        grids = list(c(6, 6), c(8, 8), c(10, 10)),
                                        a0_1 = 400,
                                        a0_2 = 400,
                                        lambda0_1 = NULL,
                                        lambda0_2 = NULL,
                                        D1 = NULL,
                                        D2 = NULL,
                                        Dmult = 4000,
                                        nocc = 5,
                                        detector = "proximity",
                                        D_model = c("group", "common")) {
  D_model <- match.arg(D_model)
  spacing <- 2 * sigma1
  lambda01 <- if (!is.null(lambda0_1)) lambda0_1 else lambda_from_a0(a0_1, sigma1, nocc)

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
    lambda02 <- if (!is.null(lambda0_2)) lambda0_2 else lambda_from_a0(a0_2, s2, nocc)
    d1_val <- if (is.null(D1)) Dmult / (sigma1^2) else D1
    d2_val <- if (is.null(D2)) Dmult / (s2^2) else D2
    for (g_idx in seq_along(trapset)) {
      scen_list[[length(scen_list) + 1]] <- data.frame(
        scenario = scenario_id,
        scenario_label = sprintf("%s_sigma2_%s", names(trapset)[g_idx], s2),
        grid_label = names(trapset)[g_idx],
        trap_id = g_idx,
        mask_id = g_idx,
        spacing = spacing,
        sigma1 = sigma1,
        sigma2 = s2,
        lambda01 = lambda01,
        lambda02 = lambda02,
        a0_1 = a0_1,
        a0_2 = a0_2,
        D1 = d1_val,
        D2 = d2_val,
        nocc = nocc,
        detector = detector,
        D_model = D_model,
        stringsAsFactors = FALSE
      )
      scenario_id <- scenario_id + 1
    }
  }

  scenarios <- bind_rows(scen_list)
  list(scenarios = scenarios, trapset = trapset, maskset = maskset, a0_1 = a0_1, a0_2 = a0_2)
}

count_by_group <- function(capthist) {
  detections <- rowSums(apply(capthist, c(1, 2), sum))
  det_by_detector <- apply(capthist, c(1, 3), sum)
  spatial_recaps_ind <- pmax(rowSums(det_by_detector > 0) - 1, 0)
  grp <- covariates(capthist)$group
  if (is.null(grp)) grp <- rep("all", length(detections))

  groups <- unique(grp)
  do.call(
    rbind,
    lapply(groups, function(g) {
      idx <- grp == g
      data.frame(
        group = g,
        n = sum(detections[idx] > 0),
        r = sum(pmax(detections[idx] - 1, 0)),
        r_spatial = sum(spatial_recaps_ind[idx])
      )
    })
  )
}

sim_and_fit <- function(scen_row, trapset, maskset, nrepl, seed, trace_fit = FALSE) {
  trap <- trapset[[scen_row$trap_id]]
  mask <- maskset[[scen_row$mask_id]]

  out_rows <- list()
  for (i in seq_len(nrepl)) {
    set.seed(seed + scen_row$scenario * 10000 + i)

    pop1 <- sim.popn(D = scen_row$D1, core = mask, buffer = 0)
    pop2 <- sim.popn(D = scen_row$D2, core = mask, buffer = 0)
    covariates(pop1)$group <- "g1"
    covariates(pop2)$group <- "g2"

    ch1 <- sim.capthist(
      traps = trap,
      popn = pop1,
      noccasions = scen_row$nocc,
      detectpar = list(lambda0 = scen_row$lambda01, sigma = scen_row$sigma1),
      detectfn = "HHN"
    )
    ch2 <- sim.capthist(
      traps = trap,
      popn = pop2,
      noccasions = scen_row$nocc,
      detectpar = list(lambda0 = scen_row$lambda02, sigma = scen_row$sigma2),
      detectfn = "HHN"
    )
    covariates(ch1)$group <- "g1"
    covariates(ch2)$group <- "g2"

    ch <- rbind(ch1, ch2)
    counts <- count_by_group(ch)

    fit <- try(
      secr.fit(
        capthist = ch,
        mask = mask,
        detectfn = "HHN",
        model = if (scen_row$D_model == "group") list(D ~ g, lambda0 ~ g, sigma ~ g) else list(D ~ 1, lambda0 ~ g, sigma ~ g),
        groups = "group",
        trace = trace_fit
      ),
      silent = !trace_fit
    )

    if (!inherits(fit, "try-error")) {
      preds <- predict(fit)
      if (scen_row$D_model == "group") {
        pred_df <- tidy_pred(preds, scen_row$scenario_label, scen_row$sigma1, scen_row$sigma2, scen_row$grid_label, i)
        pred_wide <- pred_df %>%
          select(group, parameter, estimate, SE.estimate) %>%
          pivot_wider(names_from = parameter, values_from = c(estimate, SE.estimate), names_sep = "_")

        combined <- counts %>%
          left_join(pred_wide, by = "group") %>%
          mutate(
            scenario_label = scen_row$scenario_label,
            grid_label = scen_row$grid_label,
            spacing = scen_row$spacing,
            sigma1 = scen_row$sigma1,
            sigma2 = scen_row$sigma2,
            lambda01 = scen_row$lambda01,
            lambda02 = scen_row$lambda02,
            D1 = scen_row$D1,
            D2 = scen_row$D2,
            repl = i,
            D_true = ifelse(group == "g1", scen_row$D1, scen_row$D2),
            lambda0_true = ifelse(group == "g1", scen_row$lambda01, scen_row$lambda02),
            sigma_true = ifelse(group == "g1", scen_row$sigma1, scen_row$sigma2),
            approx_rse = ifelse(!is.na(n) & !is.na(r) & n > 0 & r > 0, 1 / pmin(sqrt(n), sqrt(r)), NA_real_),
            rse_sim = SE.estimate_D / estimate_D
          )
      } else {
        beta <- coef(fit)
        V <- try(vcov(fit), silent = TRUE)
        if (inherits(V, "try-error")) {
          V <- matrix(NA_real_, nrow = nrow(beta), ncol = nrow(beta), dimnames = list(rownames(beta), rownames(beta)))
        }
        D_est <- exp(beta["D", "beta"])
        se_logD <- if (!is.na(V["D", "D"])) sqrt(V["D", "D"]) else NA_real_
        SE_D <- if (!is.na(se_logD)) D_est * se_logD else NA_real_

        g_levels <- levels(covariates(ch)$group)
        pred_group_df <- tidy_pred(preds, scen_row$scenario_label, scen_row$sigma1, scen_row$sigma2, scen_row$grid_label, i)
        pred_wide <- pred_group_df %>%
          filter(parameter %in% c("lambda0", "sigma")) %>%
          select(group, parameter, estimate, SE.estimate) %>%
          pivot_wider(names_from = parameter, values_from = c(estimate, SE.estimate), names_sep = "_")

        total_n <- sum(counts$n, na.rm = TRUE)
        total_r <- sum(counts$r, na.rm = TRUE)
        approx_total <- ifelse(total_n > 0 & total_r > 0, 1 / pmin(sqrt(total_n), sqrt(total_r)), NA_real_)

        combined <- counts %>%
          left_join(pred_wide, by = "group") %>%
          mutate(
            estimate_D = D_est,
            SE.estimate_D = SE_D,
            scenario_label = scen_row$scenario_label,
            grid_label = scen_row$grid_label,
            spacing = scen_row$spacing,
            sigma1 = scen_row$sigma1,
            sigma2 = scen_row$sigma2,
            lambda01 = scen_row$lambda01,
            lambda02 = scen_row$lambda02,
            D1 = scen_row$D1,
            D2 = scen_row$D2,
            repl = i,
            D_true = ifelse(group == "g1", scen_row$D1, scen_row$D2),
            lambda0_true = ifelse(group == "g1", scen_row$lambda01, scen_row$lambda02),
            sigma_true = ifelse(group == "g1", scen_row$sigma1, scen_row$sigma2),
            approx_rse = approx_total,
            rse_sim = ifelse(!is.na(SE_D) & !is.na(D_est), SE_D / D_est, NA_real_)
          )
      }
      out_rows[[length(out_rows) + 1]] <- combined
    } else {
      counts$scenario_label <- scen_row$scenario_label
      counts$grid_label <- scen_row$grid_label
      counts$spacing <- scen_row$spacing
      counts$sigma1 <- scen_row$sigma1
      counts$sigma2 <- scen_row$sigma2
      counts$lambda01 <- scen_row$lambda01
      counts$lambda02 <- scen_row$lambda02
      counts$D1 <- scen_row$D1
      counts$D2 <- scen_row$D2
      counts$repl <- i
      counts$D_true <- ifelse(counts$group == "g1", scen_row$D1, scen_row$D2)
      counts$lambda0_true <- ifelse(counts$group == "g1", scen_row$lambda01, scen_row$lambda02)
      counts$sigma_true <- ifelse(counts$group == "g1", scen_row$sigma1, scen_row$sigma2)
      if (scen_row$D_model == "group") {
        counts$approx_rse <- ifelse(!is.na(counts$n) & !is.na(counts$r) & counts$n > 0 & counts$r > 0, 1 / pmin(sqrt(counts$n), sqrt(counts$r)), NA_real_)
      } else {
        total_n <- sum(counts$n, na.rm = TRUE)
        total_r <- sum(counts$r, na.rm = TRUE)
        approx_total <- ifelse(total_n > 0 & total_r > 0, 1 / pmin(sqrt(total_n), sqrt(total_r)), NA_real_)
        counts$approx_rse <- approx_total
      }
      counts$rse_sim <- NA_real_
      counts$estimate_D <- NA_real_
      counts$SE.estimate_D <- NA_real_
      counts$estimate_lambda0 <- NA_real_
      counts$SE.estimate_lambda0 <- NA_real_
      counts$estimate_sigma <- NA_real_
      counts$SE.estimate_sigma <- NA_real_
      out_rows[[length(out_rows) + 1]] <- counts
    }
  }

  bind_rows(out_rows)
}

run_two_group_manual <- function(plan,
                                 nrepl = 100,
                                 seed = 123,
                                 ncores = NULL,
                                 trace_fit = FALSE) {
  scenarios_split <- split(plan$scenarios, plan$scenarios$scenario)

  if (!is.null(ncores) && ncores > 1) {
    cl <- makeCluster(ncores)
    on.exit(stopCluster(cl), add = TRUE)
    registerDoParallel(cl)

    per_scen <- foreach(
      scen_df = scenarios_split,
      .packages = c("secr", "secrdesign", "dplyr", "tidyr", "stringr"),
      .export = c("count_by_group", "sim_and_fit", "tidy_pred")
    ) %dopar% {
      scen_row <- scen_df[1, ]
      sim_and_fit(scen_row, plan$trapset, plan$maskset, nrepl, seed, trace_fit = trace_fit)
    }
  } else {
    per_scen <- lapply(scenarios_split, function(scen_df) {
      scen_row <- scen_df[1, ]
      sim_and_fit(scen_row, plan$trapset, plan$maskset, nrepl, seed, trace_fit = trace_fit)
    })
  }

  per_repl <- bind_rows(per_scen)

  empirical_rse <- per_repl %>%
    filter(!is.na(rse_sim)) %>%
    group_by(scenario_label, grid_label, sigma1, sigma2, group) %>%
    summarise(
      mean_rse = mean(rse_sim, na.rm = TRUE),
      sd_rse = sd(rse_sim, na.rm = TRUE),
      n_rep = n(),
      .groups = "drop"
    )

  approx_counts <- per_repl %>%
    filter(!is.na(n)) %>%
    group_by(scenario_label, grid_label, sigma1, sigma2, group) %>%
    summarise(
      mean_n = mean(n, na.rm = TRUE),
      mean_r = mean(r, na.rm = TRUE),
      mean_r_spatial = mean(r_spatial, na.rm = TRUE),
      approx_rse = mean(approx_rse, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    per_repl = per_repl,
    empirical_rse = empirical_rse,
    approx_counts = approx_counts
  )
}
