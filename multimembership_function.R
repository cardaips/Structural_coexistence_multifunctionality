# functions to run more efficiently the models with multimembership random factors

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
  plot <- modelplot(model, coef_omit = "SD",
                    coef_rename = c("omega"="ND", 
                                    "non.logged.theta"="FD",
                                    "differential"="ID",
                                    "non.logged.min.distance"="min. distance to exclusion",
                                    "nitrogen1" = "Nitrogen")) +
    aes(
      shape = ifelse(p.value > 0.05,
        "Not significant",
        "Significant"
      ),
      color="black"
      ) +
    facet_grid(~model) +
    scale_shape_manual(
      name = "Significance",
      values = c(4, 19),
      labels = c("Not significant", "Significant"),
      drop = FALSE
    ) +
    scale_color_manual(values = rep("black",length(model)))+
    guides(colour="none")+
  geom_vline(xintercept = 0, linetype = "dashed") +
    theme_minimal()
  return(plot)
}
