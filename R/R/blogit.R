

## Binary logistic (blogit) regression with emphasis on visualization.
## Available models include glm(...) and bayesglm(...) from `arm`
##+package. For some reason that perfect sepperation occurs, the 
##+coefficients are inflated. In this case, switch to use `bayesglm`
##+is recommended.

blogit <-
  function(...,
           clip = c(0, 8),
           user_specified_variables,
           plotting_for_sig_only = FALSE,
           main = NULL, main.xycoord = c(0.3, 1.05),
           use.glm=FALSE)
  {
  if (use.glm) {
    fit = glm(...)
  } else {
    fit <- arm::bayesglm(...)
  }
  r <- arm::display(fit, digits=6, detail=TRUE)
  if (plotting_for_sig_only && all(r$p.value[2:length(r$p.value)] > 0.05))
    return(fit)
  m <- signif(exp(r$coef), 2)
  l <- signif(exp(r$coef - 2.96 * r$se), 2)
  u <- signif(exp(r$coef + 2.96 * r$se), 2)
  p <- signif(r$p.value, 2)
  if (any(is.infinite(u))) {
    warning("Inf found for coefficients, set use.glm = FALSE is recommended.")
  }
  if (missing(user_specified_variables)) {
    variable_labels <- names(m)
  } else {
    variable_labels <- c(NA, user_specified_variables)
  }
  tabletext <- list()
  tabletext[[1]] <- c("Var", "OR", "2.5% CI", "97.5% CI", "P-value")
  for (i in 2:length(m)) {
    tabletext[[i]] <- c(variable_labels[i], m[i], l[i], u[i], p[i])
  }
  tabletext <- do.call("rbind", tabletext)
  is.summary = c(1, rep(0, 20))
  m[1] <- NA
  l[1] <- NA
  u[1] <- NA
  forestplot(
    tabletext,
    m, l, u,
    clip = clip,
    zero = 1,
    is.summary = is.summary,
    boxsize = 0.5,
    align = rep('l',5),
    xlab = "OR (95% CI)",
    col = meta.colors(box = "black", line = "black", summary = "black", zero = "black")
  )
  if (!is.null(main)) {
    x <- main.xycoord[1]
    y <- main.xycoord[2]
    text(x, y, main, xpd = TRUE)
  }
  invisible(fit)
}

library(arm)
library(rmeta)

## An example of blogit in regressing SMG mutation status versus mutational exposures
setwd("/Users/lixiangchun/Public/WorkSpace/Project/DigestiveSystemCancer/Analysis/iCGA_PanelOfNormals/StomachPublishedDataOnlyReanalyses/mutated_genes_vs_mutational_signatures")
load('SMG_vs_MSig.RData')
#b=bayesglm(as.formula(f2), data=regular, family=binomial)
fit = blogit(
  ERBB4 ~ msig1 + msig2 + msig3 + msig4 + msig5 + msig6 + msig7,
  data = regular,
  family = binomial,
  clip = c(0, 16), main="ERBB4 vs. mutational exposures",
  use.glm = FALSE
)
