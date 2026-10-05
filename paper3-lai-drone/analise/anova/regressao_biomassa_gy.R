###############################################################################
# Regressions of SPAD, leaf area (TLA/LAI) and ExG on aboveground biomass with
# ear (AGBce) and grain yield (GY). Data: Livro2.xlsx (n = 24 plots).
# Run from paper3-lai-drone/:  Rscript analise/anova/regressao_biomassa_gy.R
###############################################################################
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
invisible(Sys.setlocale("LC_ALL", "C.UTF-8"))
SAIDA <- "resultados/anova_tukey"; dir.create(SAIDA, recursive = TRUE, showWarnings = FALSE)

d <- merge(read.csv("dados/preliminar_exg_tla.csv"), read.csv("dados/livro2_biomassa_gy.csv"), by = "Parcela")
names(d)[names(d) == "TLA_cm2_planta"] <- "TLA"; names(d)[names(d) == "ExG_DN_QGIS"] <- "ExG"
d$LAI <- d$TLA * 10 / 1e4
d$Biochar <- factor(d$Biochar, levels = c("BC0", "BC12"))
d$Lime    <- factor(d$Calcario, levels = c("L0", "L75", "L100"))
stopifnot(nrow(d) == 24)

lab <- list(
  SPAD  = "SPAD index",
  TLA   = expression("Leaf area (cm"^2*" plant"^-1*")"),
  LAI   = expression("LAI (m"^2*" m"^-2*")"),
  ExG   = "ExG (DN)",
  AGBce = "Biomass with ear",          # aboveground, with ear [VERIFICAR unit]
  GY    = expression("Grain yield (kg ha"^-1*")"))
num <- function(x) sub("^-", "−", sub("\\.$", "", trimws(formatC(x, digits = 3, format = "fg", flag = "#"))))
sg  <- function(x) ifelse(x < 0, " − ", " + ")
pt  <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)

# ---- table: linear and quadratic fits ---------------------------------------
pares <- expand.grid(x = c("SPAD", "TLA", "LAI", "ExG"), y = c("AGBce", "GY"), stringsAsFactors = FALSE)
tab <- do.call(rbind, lapply(seq_len(nrow(pares)), function(i) {
  x <- pares$x[i]; y <- pares$y[i]
  f1 <- lm(reformulate(x, y), d); s1 <- summary(f1); ct <- cor.test(d[[x]], d[[y]])
  f2 <- lm(as.formula(sprintf("%s ~ %s + I(%s^2)", y, x, x)), d); s2 <- summary(f2)
  data.frame(Response = y, Predictor = x,
    Equation = paste0("y = ", num(coef(f1)[1]), sg(coef(f1)[2]), num(abs(coef(f1)[2])), "x"),
    r = round(ct$estimate, 3), r_CI95 = sprintf("%.2f to %.2f", ct$conf.int[1], ct$conf.int[2]),
    R2 = round(s1$r.squared, 3), R2_adj = round(s1$adj.r.squared, 3), F = round(s1$fstatistic[1], 2),
    p = signif(coef(s1)[2, 4], 3), RMSE = signif(sigma(f1), 4),
    R2_quadratic = round(s2$r.squared, 3), p_quadratic_term = signif(coef(s2)[3, 4], 3), n = nrow(d))
}))
# extra models: leaf area x SPAD (canopy chlorophyll proxy) and two-predictor models
d$LAIxSPAD <- d$LAI * d$SPAD
extra <- do.call(rbind, lapply(c("AGBce", "GY"), function(y) do.call(rbind, lapply(
  c("LAIxSPAD", "TLA + SPAD", "TLA + ExG", "TLA + SPAD + ExG"), function(x) {
    f <- lm(as.formula(paste(y, "~", x)), d); s <- summary(f)
    data.frame(Response = y, Model = x, R2 = round(s$r.squared, 3), R2_adj = round(s$adj.r.squared, 3),
               p_model = signif(pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE), 3),
               AIC = round(AIC(f), 1))
  }))))
write.csv(tab, file.path(SAIDA, "Table_regression_biomass_yield.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(extra, file.path(SAIDA, "Table_regression_biomass_yield_extra.csv"), row.names = FALSE, fileEncoding = "UTF-8")
print(tab, row.names = FALSE); print(extra, row.names = FALSE)

# ---- figure: 2 rows (AGBce, GY) x 3 columns (SPAD, TLA, ExG) -----------------
tema <- theme_classic(base_size = 10, base_family = "sans") +
  theme(axis.text = element_text(colour = "black"), axis.ticks = element_line(colour = "black"),
        plot.tag = element_text(face = "bold"), plot.subtitle = element_text(size = 7), plot.margin = margin(4, 8, 4, 4))
painel <- function(y, x) {
  f <- lm(reformulate(x, y), d); s <- summary(f)
  sub <- sprintf("y = %s%s%sx\nr = %.2f; R² = %.2f; %s", num(coef(f)[1]), sg(coef(f)[2]), num(abs(coef(f)[2])),
                 cor(d[[x]], d[[y]]), s$r.squared, pt(coef(s)[2, 4]))
  ggplot(d, aes(.data[[x]], .data[[y]])) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black", fill = "grey80", linewidth = 0.5) +
    geom_point(aes(shape = Biochar, fill = Lime), size = 2.2, stroke = 0.4) +
    scale_shape_manual(values = c(21, 24)) + scale_fill_manual(values = c("white", "grey60", "black")) +
    labs(x = lab[[x]], y = lab[[y]], subtitle = sub) + tema
}
g <- (painel("AGBce", "SPAD") | painel("AGBce", "TLA") | painel("AGBce", "ExG")) /
     (painel("GY", "SPAD")    | painel("GY", "TLA")    | painel("GY", "ExG")) +
  plot_annotation(tag_levels = "a") + plot_layout(guides = "collect") &
  theme(legend.position = "bottom") &
  guides(shape = guide_legend(title = "Biochar", order = 1),
         fill = guide_legend(title = "Lime", order = 2, override.aes = list(shape = 21)))
ggsave(file.path(SAIDA, "Fig_regression_biomass_yield.png"), g, width = 190, height = 150, units = "mm", dpi = 600)
ggsave(file.path(SAIDA, "Fig_regression_biomass_yield.pdf"), g, width = 190, height = 150, units = "mm")
