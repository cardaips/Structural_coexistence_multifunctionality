# functions to run and predict more efficiently the models with multimembership random factors

# run and plot the models ####

multimembership_model <- function(formula, presence_matrix, data) {
  # Omega and differential only
  lmod <- lFormula(formula, data = data)
  # I modified the random factor by introducing "speciesmat" which is a presence/absence matrix of PaNDiv species in my communities of 5 plant species.
  lmod$reTrms$Zt <- lmod$reTrms$Ztlist[[1]] <- Matrix(t(presence_matrix))
  # preparing to run the optimisation of the lmer
  devfun <- do.call(mkLmerDevfun, lmod)
  # optimising per se
  opt <- optimizeLmer(devfun)
  # run the optimisation and save the model into a variable called m1, m2,... to m5 -> important !
  model <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
  return(model)
}

# plotting the multimembership so that model output look beautiful
plot_multi <- function(model) {
  plot <- modelplot(model,
    coef_omit = "SD|Intercept",
    coef_rename = c(
      "omega" = "ND",
      "structural.niche" = "ND",
      "non.logged.theta" = "FD",
      "structural.fitness" = "FD",
      "differential" = "ID",
      "indirect.interactions" = "ID",
      "non.logged.min.distance" = "min. distance to exclusion",
      "min.distance" = "min. distance to exclusion",
      "threshold.continuous" = "threshold",
      "nitrogen1" = "Nitrogen"
    )
  ) +
    aes(
      shape = ifelse(p.value > 0.05,
        "Not significant",
        "Significant"
      ),
      color = "black"
    ) +
    facet_grid(~model) +
    scale_shape_manual(
      name = "Significance",
      values = c(4, 19),
      labels = c("Not significant", "Significant"),
      drop = FALSE
    ) +
    scale_color_manual(values = rep("black", length(model))) +
    guides(colour = "none") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +
    theme_minimal()+
    theme(text= element_text(size = 15, family = "serif"))
  return(plot)
}

# predict with the model ####

# predict for minimum dstance to exclusion
predict_multifunctionality_dist <- function(model, new.data) {
  # useful prediction function
  pfun <- function(.) {
    predict(., newdata = new.data, re.form = ~0, type = "response")
  }

  # summarise output of bootstrapping
  sumBoot <- function(merBoot) {
    return(
      data.frame(
        fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs = .5, na.rm = TRUE))),
        lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs = .025, na.rm = TRUE))),
        upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs = .975, na.rm = TRUE)))
      )
    )
  }

  bootdata <- lme4::bootMer(model, pfun, nsim = 1000, type = "parametric")

  sum <- sumBoot(bootdata)


  new.data$fit <- sum$fit
  new.data$lwr <- sum$lwr
  new.data$upr <- sum$upr

  return(new.data)
}

# predict for each coex mechanism
predict_multifunctionality_coex <- function(model, new.data) {
  # useful prediction function
  pfun <- function(.) {
    predict(., newdata = new.data, re.form = ~0, type = "response")
  }

  # summarise output of bootstrapping
  sumBoot <- function(merBoot) {
    return(
      data.frame(
        fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs = .5, na.rm = TRUE))),
        lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs = .025, na.rm = TRUE))),
        upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs = .975, na.rm = TRUE)))
      )
    )
  }

  bootdata <- lme4::bootMer(model, pfun, nsim = 100, type = "parametric")
  sum <- sumBoot(bootdata)

  new.data$fit <- sum$fit

  return(new.data)
}
