###############################################################################
# Regression, factorial ANOVA (RCBD) and Tukey for TLA/LAI, SPAD and ExG
# Data: dados/tabela_v2_2026-10-01.csv (table sent by Filipe on 2026-10-01)
# Design: RCBD, 2 (biochar 0, 12 Mg ha-1) x 3 (lime 0, 75, 100% of the
# recommended rate), 4 blocks. LAI = TLA x 10 plants m-2 / 10 000.
# Run from paper3-lai-drone/:  Rscript analise/anova/anova_tukey_regressao.R
###############################################################################
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(emmeans); library(multcomp)
  library(multcompView); library(car)
})
invisible(Sys.setlocale("LC_ALL", "C.UTF-8"))
SAIDA <- "resultados/anova_tukey"; dir.create(SAIDA, recursive = TRUE, showWarnings = FALSE)
ALFA <- 0.05

d <- read.csv("dados/tabela_v2_2026-10-01.csv")
d$LAI <- d$TLA * 10 / 1e4
d$Biochar <- factor(d$BC, levels = c("BC0", "BC12"), labels = c("0", "12"))
d$Lime    <- factor(d$L, levels = c("L0", "L75", "L100"), labels = c("0", "75", "100"))
d$Block   <- factor(d$BL)

vars <- list(
  LAI  = list(lab = expression("LAI (m"^2*" m"^-2*")"),            txt = "LAI (m2 m-2)",  dig = 2),
  TLA  = list(lab = expression("Leaf area (cm"^2*" plant"^-1*")"), txt = "TLA (cm2 plant-1)", dig = 0),
  SPAD = list(lab = "SPAD index",                                   txt = "SPAD index",   dig = 1),
  ExG  = list(lab = "ExG (DN)",                                     txt = "ExG (DN)",     dig = 1))
log <- c(sprintf("n = %d plots. LAI computed from TLA x 10 plants m-2.", nrow(d)),
         sprintf("Plots where reported LAI differs from TLA x 10/10000 by > 0.01: %s",
                 paste(d$Parcela[abs(d$LAI_informado - d$LAI) > 0.01], collapse = ", ")))

# ---- 1. Regressions -------------------------------------------------------
pares <- list(c("LAI", "ExG"), c("LAI", "SPAD"), c("SPAD", "ExG"), c("TLA", "ExG"))
reg <- do.call(rbind, lapply(pares, function(p) {
  f <- lm(reformulate(p[2], p[1]), d); s <- summary(f); ct <- cor.test(d[[p[2]]], d[[p[1]]])
  data.frame(Response = p[1], Predictor = p[2],
             Equation = sprintf("y = %.4g %+.4g x", coef(f)[1], coef(f)[2]),
             r = ct$estimate, r_CI95 = sprintf("%.2f to %.2f", ct$conf.int[1], ct$conf.int[2]),
             R2 = s$r.squared, R2_adj = s$adj.r.squared, F = s$fstatistic[1],
             p = coef(s)[2, 4], RMSE = sigma(f), n = nrow(d))
}))
fm <- lm(LAI ~ ExG + SPAD, d); sm <- summary(fm)
reg <- rbind(reg, data.frame(Response = "LAI", Predictor = "ExG + SPAD",
  Equation = sprintf("y = %.4g %+.4g ExG %+.4g SPAD", coef(fm)[1], coef(fm)[2], coef(fm)[3]),
  r = sqrt(sm$r.squared), r_CI95 = "(multiple R)", R2 = sm$r.squared, R2_adj = sm$adj.r.squared,
  F = sm$fstatistic[1],
  p = pf(sm$fstatistic[1], sm$fstatistic[2], sm$fstatistic[3], lower.tail = FALSE),
  RMSE = sigma(fm), n = nrow(d)))
reg_out <- reg; reg_out[, c("r", "R2", "R2_adj", "F", "RMSE")] <- round(reg_out[, c("r", "R2", "R2_adj", "F", "RMSE")], 3)
reg_out$p <- signif(reg_out$p, 3)
write.csv(reg_out, file.path(SAIDA, "Table_regression.csv"), row.names = FALSE)
log <- c(log, "", "REGRESSIONS", capture.output(print(reg_out, row.names = FALSE)),
         sprintf("Partial p in LAI ~ ExG + SPAD: ExG p = %.4f; SPAD p = %.4f; VIF = %.2f",
                 coef(sm)[2, 4], coef(sm)[3, 4], vif(fm)[1]))

tema <- theme_classic(base_size = 10, base_family = "sans") +
  theme(axis.text = element_text(colour = "black"), axis.ticks = element_line(colour = "black"),
        legend.position = "none", plot.tag = element_text(face = "bold"))
pt <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
painel_reg <- function(y, x) {
  f <- lm(reformulate(x, y), d); s <- summary(f)
  lab <- sprintf("y = %.3g %s %.3gx\nr = %.2f; R² = %.2f; %s", coef(f)[1],
                 ifelse(coef(f)[2] < 0, "−", "+"), abs(coef(f)[2]),
                 cor(d[[x]], d[[y]]), s$r.squared, pt(coef(s)[2, 4]))
  ggplot(d, aes(.data[[x]], .data[[y]])) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black", fill = "grey80", linewidth = 0.5) +
    geom_point(aes(shape = Biochar, fill = Lime), size = 2.2, stroke = 0.4) +
    scale_shape_manual(values = c(21, 24)) +
    scale_fill_manual(values = c("white", "grey60", "black")) +
    labs(x = vars[[x]]$lab, y = vars[[y]]$lab, subtitle = lab) +
    tema + theme(plot.subtitle = element_text(size = 7.5))
}
g_reg <- (painel_reg("LAI", "ExG") | painel_reg("LAI", "SPAD") | painel_reg("SPAD", "ExG")) +
  plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")
g_reg <- g_reg + plot_layout(guides = "collect") &
  guides(shape = guide_legend(title = expression("Biochar (Mg ha"^-1*")")),
         fill = guide_legend(title = "Lime (% of rate)", override.aes = list(shape = 21)))
ggsave(file.path(SAIDA, "Fig_regression.png"), g_reg, width = 180, height = 75, units = "mm", dpi = 600)
ggsave(file.path(SAIDA, "Fig_regression.pdf"), g_reg, width = 180, height = 75, units = "mm")

# ---- 2. ANOVA, assumptions, Tukey -----------------------------------------
letras <- function(em) {
  cl <- cld(em, Letters = letters, adjust = "tukey", alpha = ALFA, reversed = TRUE)  # highest mean = a
  cl$.group <- gsub(" ", "", cl$.group); as.data.frame(cl)
}
resumo <- function(v, g) {
  a <- aggregate(d[[v]], d[g], function(x) c(m = mean(x), sd = sd(x), n = length(x)))
  a <- cbind(a[g], as.data.frame(a$x))
  a$ci95 <- qt(0.975, a$n - 1) * a$sd / sqrt(a$n); a
}
tab_anova <- list(); tab_medias <- list(); figs <- list()
for (v in names(vars)) {
  f <- lm(reformulate(c("Block", "Biochar * Lime"), v), d)
  av <- anova(f)
  sw <- shapiro.test(residuals(f)); lv <- leveneTest(reformulate("Biochar:Lime", v), d)
  p_int <- av["Biochar:Lime", "Pr(>F)"]
  cv <- 100 * sqrt(av["Residuals", "Mean Sq"]) / mean(d[[v]])
  tab_anova[[v]] <- data.frame(Variable = vars[[v]]$txt, Source = rownames(av), df = av$Df,
    SS = signif(av$`Sum Sq`, 4), MS = signif(av$`Mean Sq`, 4), F = round(av$`F value`, 2),
    p = signif(av$`Pr(>F)`, 3), CV_pct = c(rep(NA, nrow(av) - 1), round(cv, 1)),
    Shapiro_p = c(rep(NA, nrow(av) - 1), round(sw$p.value, 3)),
    Levene_p = c(rep(NA, nrow(av) - 1), round(lv$`Pr(>F)`[1], 3)))
  dig <- vars[[v]]$dig
  if (p_int < ALFA) {
    # Interaction: uppercase = biochar within each lime rate; lowercase = lime within each biochar rate
    up <- letras(emmeans(f, ~ Biochar | Lime)); up$.group <- toupper(up$.group)
    lo <- letras(emmeans(f, ~ Lime | Biochar))
    s <- resumo(v, c("Biochar", "Lime"))
    s <- merge(s, up[, c("Biochar", "Lime", ".group")], by = c("Biochar", "Lime"))
    names(s)[names(s) == ".group"] <- "up"
    s <- merge(s, lo[, c("Biochar", "Lime", ".group")], by = c("Biochar", "Lime"))
    names(s)[names(s) == ".group"] <- "lo"
    s$lab <- paste0(s$up, s$lo)
    s$Lime <- factor(s$Lime, levels = levels(d$Lime)); s <- s[order(s$Biochar, s$Lime), ]
    tab_medias[[v]] <- data.frame(Variable = vars[[v]]$txt, Effect = "Biochar x Lime",
      Level = paste0("BC", s$Biochar, " L", s$Lime),
      Mean_SD = sprintf(paste0("%.", dig, "f ± %.", dig, "f"), s$m, s$sd),
      CI95 = sprintf(paste0("± %.", dig, "f"), s$ci95), Tukey = s$lab)
    s$ytxt <- s$m + s$sd + 0.04 * max(s$m + s$sd)
    figs[[v]] <- ggplot(s, aes(Lime, m, fill = Biochar)) +
      geom_col(position = position_dodge(0.8), width = 0.7, colour = "black", linewidth = 0.3) +
      geom_errorbar(aes(ymin = m - sd, ymax = m + sd), position = position_dodge(0.8), width = 0.2, linewidth = 0.3) +
      geom_text(aes(y = ytxt, label = lab, group = Biochar), position = position_dodge(0.8), size = 2.8, vjust = 0) +
      scale_fill_manual(values = c("white", "grey55"), name = expression("Biochar (Mg ha"^-1*")")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
      labs(x = "Lime (% of recommended rate)", y = vars[[v]]$lab,
           subtitle = sprintf("B: %s; L: %s; B × L: %s", pt(av["Biochar", "Pr(>F)"]),
                              pt(av["Lime", "Pr(>F)"]), pt(p_int))) +
      tema + theme(legend.position = "top", plot.subtitle = element_text(size = 7.5))
  } else {
    # No interaction: main effects in separate panels
    pan <- list(); linhas <- list()
    for (fa in c("Biochar", "Lime")) {
      em <- emmeans(f, as.formula(paste("~", fa)))
      cl <- letras(em); s <- resumo(v, fa)
      s <- merge(s, cl[, c(fa, ".group")], by = fa)
      s[[fa]] <- factor(s[[fa]], levels = levels(d[[fa]])); s <- s[order(s[[fa]]), ]
      p_fa <- av[fa, "Pr(>F)"]
      if (p_fa >= ALFA) s$.group <- ""    # no letters when the F test is not significant
      linhas[[fa]] <- data.frame(Variable = vars[[v]]$txt, Effect = fa,
        Level = paste0(ifelse(fa == "Biochar", "BC", "L"), s[[fa]]),
        Mean_SD = sprintf(paste0("%.", dig, "f ± %.", dig, "f"), s$m, s$sd),
        CI95 = sprintf(paste0("± %.", dig, "f"), s$ci95),
        Tukey = ifelse(s$.group == "", "ns", s$.group))
      s$ytxt <- s$m + s$sd + 0.04 * max(s$m + s$sd)
      pan[[fa]] <- ggplot(s, aes(.data[[fa]], m)) +
        geom_col(width = 0.6, fill = ifelse(fa == "Biochar", "grey80", "grey55"), colour = "black", linewidth = 0.3) +
        geom_errorbar(aes(ymin = m - sd, ymax = m + sd), width = 0.15, linewidth = 0.3) +
        geom_text(aes(y = ytxt, label = .group), size = 3, vjust = 0) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
        labs(x = if (fa == "Biochar") expression("Biochar (Mg ha"^-1*")") else "Lime (% of recommended rate)",
             y = vars[[v]]$lab, subtitle = sprintf("%s (F = %.2f)", pt(p_fa), av[fa, "F value"])) +
        tema + theme(plot.subtitle = element_text(size = 7.5))
    }
    tab_medias[[v]] <- do.call(rbind, linhas)
    figs[[v]] <- (pan$Biochar | pan$Lime) + plot_layout(widths = c(2, 3))
  }
  log <- c(log, "", sprintf("ANOVA %s (interaction p = %.4f -> %s)", v, p_int,
                            ifelse(p_int < ALFA, "interaction plot", "main effects")),
           capture.output(print(av)),
           sprintf("CV = %.1f%%; Shapiro-Wilk p = %.3f; Levene p = %.3f", cv, sw$p.value, lv$`Pr(>F)`[1]))
}
TA <- do.call(rbind, tab_anova); TM <- do.call(rbind, tab_medias)
write.csv(TA, file.path(SAIDA, "Table_ANOVA.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
write.csv(TM, file.path(SAIDA, "Table_means_Tukey.csv"), row.names = FALSE, fileEncoding = "UTF-8")
log <- c(log, "", "MEANS AND TUKEY", capture.output(print(TM, row.names = FALSE)))
writeLines(log, file.path(SAIDA, "log_analysis.txt")); cat(log, sep = "\n")

for (v in names(vars)) {
  ggsave(file.path(SAIDA, sprintf("Fig_bars_%s.png", v)), figs[[v]], width = 140, height = 75, units = "mm", dpi = 600)
}
g_all <- wrap_plots(figs[c("LAI", "SPAD", "ExG")], ncol = 1) + plot_annotation(tag_levels = "a")
ggsave(file.path(SAIDA, "Fig_bars_LAI_SPAD_ExG.png"), g_all, width = 140, height = 210, units = "mm", dpi = 600)
ggsave(file.path(SAIDA, "Fig_bars_LAI_SPAD_ExG.pdf"), g_all, width = 140, height = 210, units = "mm")
