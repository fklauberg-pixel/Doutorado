###############################################################################
# LAI do milho x ExG extraído no QGIS (um índice por parcela)
# Regressão linear, quadrática e Random Forest, com validação cruzada
#
# Dados: dados/preliminar_exg_tla.csv (Livro2.xlsx, conferido com EXG.csv)
#   TLA_cm2_planta = média de todas as folhas de 5 plantas/parcela (CI-202)
#   ExG_DN_QGIS    = ExG por parcela calculado no QGIS (DN 0-255)
# Voo e medição no mesmo dia (12/12/2022). Densidade 10 plantas m-2.
#
# Rodar a partir da pasta paper3-lai-drone/:  Rscript analise/R/lai_exg_qgis.R
###############################################################################
suppressPackageStartupMessages({library(randomForest); library(ggplot2)})
set.seed(2026)

DENSIDADE  <- 10
N_ARVORES  <- 500
N_REP_7030 <- 500
N_REP_KF   <- 100
SAIDA <- "resultados/exg_qgis"
dir.create(SAIDA, recursive = TRUE, showWarnings = FALSE)

dados <- read.csv("dados/preliminar_exg_tla.csv")
dados$LAI   <- dados$TLA_cm2_planta * DENSIDADE / 10000
dados$ExG   <- dados$ExG_DN_QGIS
dados$bloco <- factor(dados$Bloco)
dados$trat  <- factor(paste0(dados$Biochar, dados$Calcario),
                      levels = c("BC0L0", "BC0L75", "BC0L100", "BC12L0", "BC12L75", "BC12L100"))
n <- nrow(dados); stopifnot(n == 24)

metricas <- function(obs, pred) {
  c(R2    = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2),
    RMSE  = sqrt(mean((obs - pred)^2)),
    rRMSE = 100 * sqrt(mean((obs - pred)^2)) / mean(obs),
    MAE   = mean(abs(obs - pred)),
    vies  = mean(pred - obs))
}

ajustar_prever <- function(tr, te, modelo) {
  switch(modelo,
    Nulo       = rep(mean(tr$LAI), nrow(te)),
    Linear     = predict(lm(LAI ~ ExG, tr), te),
    Quadratica = predict(lm(LAI ~ ExG + I(ExG^2), tr), te),
    Exponencial = exp(predict(lm(log(LAI) ~ ExG, tr), te)),
    RF         = predict(randomForest(LAI ~ ExG, data = tr, ntree = N_ARVORES,
                                      nodesize = 3), te))
}
modelos <- c("Nulo", "Linear", "Quadratica", "Exponencial", "RF")

# ---- 1. Descritiva e ajuste com todos os dados -----------------------------
ct <- cor.test(dados$ExG, dados$LAI)
fit_lin <- lm(LAI ~ ExG, dados)
fit_quad <- lm(LAI ~ ExG + I(ExG^2), dados)
desc <- c(
  sprintf("n = %d parcelas; densidade = %d plantas m-2", n, DENSIDADE),
  sprintf("LAI: média %.3f, DP %.3f, CV %.1f%%, mín %.3f, máx %.3f",
          mean(dados$LAI), sd(dados$LAI), 100 * sd(dados$LAI) / mean(dados$LAI),
          min(dados$LAI), max(dados$LAI)),
  sprintf("ExG: média %.2f, DP %.2f, mín %.2f, máx %.2f",
          mean(dados$ExG), sd(dados$ExG), min(dados$ExG), max(dados$ExG)),
  sprintf("r de Pearson = %.3f (IC95%% %.3f a %.3f), p = %.4f",
          ct$estimate, ct$conf.int[1], ct$conf.int[2], ct$p.value),
  sprintf("Linear (todos os dados): LAI = %.4f + %.5f ExG; R2 = %.3f; p = %.4f",
          coef(fit_lin)[1], coef(fit_lin)[2], summary(fit_lin)$r.squared,
          summary(fit_lin)$coefficients[2, 4]),
  sprintf("Quadrática (todos os dados): R2 = %.3f; p do termo ExG^2 = %.4f",
          summary(fit_quad)$r.squared, summary(fit_quad)$coefficients[3, 4]))
# ExG e LAI respondem aos tratamentos? (ANOVA em DBC)
av_lai <- anova(lm(LAI ~ bloco + trat, dados))
av_exg <- anova(lm(ExG ~ bloco + trat, dados))
desc <- c(desc,
  sprintf("ANOVA DBC, efeito de tratamento: LAI F = %.2f, p = %.4f; ExG F = %.2f, p = %.4f",
          av_lai["trat", "F value"], av_lai["trat", "Pr(>F)"],
          av_exg["trat", "F value"], av_exg["trat", "Pr(>F)"]))
# Correlação entre médias de tratamento (remove ruído de parcela)
med <- aggregate(cbind(LAI, ExG) ~ trat, dados, mean)
ctm <- cor.test(med$ExG, med$LAI)
desc <- c(desc, sprintf("Médias dos 6 tratamentos: r = %.3f, p = %.4f", ctm$estimate, ctm$p.value))
writeLines(desc, file.path(SAIDA, "descritiva.txt")); cat(desc, sep = "\n")

# ---- 2. Validação cruzada ---------------------------------------------------
cv_pooled <- function(folds) {           # folds: lista de índices de teste
  sapply(modelos, function(m) {
    pred <- rep(NA_real_, n)
    for (te in folds) pred[te] <- ajustar_prever(dados[-te, ], dados[te, ], m)
    pred
  })
}
# 2a. Deixa um bloco fora (principal)
f_bl <- split(seq_len(n), dados$bloco)
pred_bl <- cv_pooled(f_bl)
res_bl <- t(apply(pred_bl, 2, function(p) metricas(dados$LAI, p)))
# 2b. LOOCV
pred_lo <- cv_pooled(as.list(seq_len(n)))
res_lo <- t(apply(pred_lo, 2, function(p) metricas(dados$LAI, p)))
# 2c. 5-fold repetido
kf <- lapply(seq_len(N_REP_KF), function(i) {
  f <- split(sample(n), rep(1:5, length.out = n))
  p <- cv_pooled(f)
  t(apply(p, 2, function(pp) metricas(dados$LAI, pp)))
})
res_kf <- Reduce(`+`, kf) / N_REP_KF
# 2d. 500 partições 70/30 (Du et al., 2022), iguais para todos os modelos
part <- replicate(N_REP_7030, sample(n, round(0.7 * n)), simplify = FALSE)
r7030 <- lapply(modelos, function(m) t(sapply(part, function(tr)
  metricas(dados$LAI[-tr], ajustar_prever(dados[tr, ], dados[-tr, ], m)))))
names(r7030) <- modelos
res_7030 <- t(sapply(r7030, colMeans))
pct_melhor_nulo <- sapply(modelos, function(m)
  100 * mean(r7030[[m]][, "RMSE"] < r7030[["Nulo"]][, "RMSE"]))

tab <- rbind(
  data.frame(esquema = "deixa_um_bloco_fora", modelo = modelos, res_bl, pct_RMSE_menor_que_nulo = NA),
  data.frame(esquema = "LOOCV", modelo = modelos, res_lo, pct_RMSE_menor_que_nulo = NA),
  data.frame(esquema = "5fold_x100", modelo = modelos, res_kf, pct_RMSE_menor_que_nulo = NA),
  data.frame(esquema = "70_30_x500", modelo = modelos, res_7030,
             pct_RMSE_menor_que_nulo = pct_melhor_nulo))
tab[, 3:7] <- round(tab[, 3:7], 4)
write.csv(tab, file.path(SAIDA, "tab_desempenho.csv"), row.names = FALSE)
print(tab, row.names = FALSE)

# ---- 3. Rendimento (GY do Livro2) ------------------------------------------

if (file.exists("dados/preliminar_gy.csv")) {
  g <- read.csv("dados/preliminar_gy.csv")
  dados <- merge(dados, g, by = "Parcela")
  rg <- c(sprintf("r(LAI, GY) = %.3f, p = %.4f", cor(dados$LAI, dados$GY), cor.test(dados$LAI, dados$GY)$p.value),
          sprintf("r(ExG, GY) = %.3f, p = %.4f", cor(dados$ExG, dados$GY), cor.test(dados$ExG, dados$GY)$p.value))
  writeLines(rg, file.path(SAIDA, "rendimento.txt")); cat(rg, sep = "\n")
}

# ---- 4. Figuras -------------------------------------------------------------
tema <- theme_bw(base_size = 11)
g1 <- ggplot(dados, aes(ExG, LAI)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.6) +
  geom_point(aes(shape = trat), size = 2.4) +
  scale_shape_manual(values = c(1, 2, 0, 16, 17, 15)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 3.5,
           label = sprintf("r = %.2f (p = %.3f)\nR2 ajuste = %.2f",
                           ct$estimate, ct$p.value, summary(fit_lin)$r.squared)) +
  labs(x = "ExG (QGIS)", y = expression("LAI (m"^2*" m"^-2*")"), shape = "Tratamento") + tema
ggsave(file.path(SAIDA, "fig_ExG_LAI.png"), g1, width = 6, height = 4.2, dpi = 300)

df <- data.frame(LAI = dados$LAI, bloco = dados$bloco, pred_bl[, c("Linear", "RF")])
df <- rbind(data.frame(LAI = df$LAI, bloco = df$bloco, modelo = "Linear", pred = df$Linear),
            data.frame(LAI = df$LAI, bloco = df$bloco, modelo = "RF", pred = df$RF))
rot <- data.frame(modelo = c("Linear", "RF"),
                  lab = sprintf("R² = %.2f\nRMSE = %.3f", res_bl[c("Linear", "RF"), "R2"],
                                res_bl[c("Linear", "RF"), "RMSE"]))
lim <- range(c(df$LAI, df$pred))
g2 <- ggplot(df, aes(LAI, pred)) +
  geom_abline(linetype = 2, colour = "grey50") +
  geom_point(aes(shape = bloco), size = 2.2) +
  geom_text(data = rot, aes(x = -Inf, y = Inf, label = lab), hjust = -0.1, vjust = 1.2, size = 3.3) +
  facet_wrap(~ modelo) + coord_equal(xlim = lim, ylim = lim) +
  labs(x = expression("LAI medido, CI-202 (m"^2*" m"^-2*")"),
       y = expression("LAI estimado (m"^2*" m"^-2*")"), shape = "Bloco") + tema
ggsave(file.path(SAIDA, "fig_obs_pred_bloco.png"), g2, width = 7.5, height = 3.8, dpi = 300)

d7 <- do.call(rbind, lapply(modelos, function(m) data.frame(modelo = m, R2 = r7030[[m]][, "R2"])))
d7$modelo <- factor(d7$modelo, levels = modelos)
g3 <- ggplot(d7, aes(modelo, R2)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_boxplot(outlier.size = 0.5) + labs(x = NULL, y = "R² no teste (500 partições 70/30)") + tema
ggsave(file.path(SAIDA, "fig_boxplot_R2_7030.png"), g3, width = 5.5, height = 4, dpi = 300)
