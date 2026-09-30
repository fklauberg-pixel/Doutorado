###############################################################################
# Estimativa da área foliar do milho (TLA / LAI) com índices RGB de drone
# Regressão linear x Random Forest, com validação cruzada
#
# Experimento: DBC, 4 blocos x 6 tratamentos (biochar de açaí x calcário)
#   BC0L0, BC0L75, BC0L100, BC12L0, BC12L75, BC12L100  -> 24 parcelas
# Área foliar: todas as folhas de 5 plantas/parcela (CI-202), média por parcela
# Imagem: 1 voo RGB (VT-R1), ortomosaico do Pix4D
#
# Base metodológica: Du et al. (2022) Sci Rep 12:15937
#   - LAI = área foliar por planta x densidade (plantas m-2)
#   - bandas normalizadas r = R/(R+G+B) etc.
#   - índices RGB; seleção de VIs com |r| > 0,7; RF com 500 árvores
#   - ULR, MLR e RF comparados por R2 e RMSE em dados de teste
#
# Versão revisada do script do Filipe (30/09/2026). Mudanças marcadas com [REV]:
#   [REV1] R2 de validação = 1 - SQres/SQtot (o cor()^2 ignora viés e escala)
#   [REV2] modelo nulo (média do treino) como linha de base
#   [REV3] seleção de VIs poda índices redundantes (|r| > 0,95 entre si)
#   [REV4] mtry do RF escolhido pelo erro OOB dentro do treino; RF com todos
#          os VIs também avaliado; importância vem do mesmo modelo validado
#   [REV5] as 500 partições 70/30 são as MESMAS para todos os modelos
#          (comparação pareada: % de partições em que o RF vence)
#   [REV6] na seção 8, o melhor VI é escolhido dentro de cada treino
#   [REV7] extração ignora pixels transparentes (banda alfa do Pix4D)
###############################################################################

# ---- 0. Pacotes -------------------------------------------------------------
pkgs <- c("terra", "randomForest", "ggplot2", "dplyr", "tidyr")
novos <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(novos)) install.packages(novos)
invisible(lapply(pkgs, library, character.only = TRUE))
set.seed(2026)

# ---- 1. Parâmetros ----------------------------------------------------------
USAR_DADOS_DEMO <- TRUE     # TRUE = testa o script com dados simulados
                            # FALSE = usa seus arquivos (seção 2)

DENSIDADE <- 10             # plantas m-2 (2 plantas/cova, 0,80 x 0,25 m)
                            # troque por uma coluna 'estande' se tiver contagem
LIMIAR_R  <- 0.70           # seleção de VIs, como em Du et al. (2022)
LIMIAR_REDUND <- 0.95       # [REV3] |r| máximo entre VIs escolhidos
N_ARVORES <- 500            # idem
N_REP_7030 <- 500           # repetições da partição 70/30 (idem)

# ---- 2. Extração dos índices do ortomosaico (rode uma vez) -----------------
# Requer: ortomosaico RGB do Pix4D (GeoTIFF) e um vetor com os polígonos das
# parcelas (GeoPackage/shapefile) com a coluna 'parcela'.
# Estratégia: média de R, G, B na área útil da parcela -> bandas normalizadas
# -> índices. Calcular os índices sobre as médias evita a instabilidade de
# índices razão (VARI, RGRI) pixel a pixel, quando o denominador ~ 0.
extrair_bandas <- function(orto_path, parcelas_path, recuo_m = 0.5) {
  orto <- terra::rast(orto_path)
  rgb  <- orto[[1:3]]
  names(rgb) <- c("R", "G", "B")
  if (terra::nlyr(orto) >= 4) {                   # [REV7] alfa = 0 -> NA
    rgb <- terra::mask(rgb, orto[[4]], maskvalues = 0)
  }
  parc <- terra::vect(parcelas_path)
  parc <- terra::project(parc, terra::crs(rgb))
  util <- terra::buffer(parc, width = -recuo_m)   # exclui bordadura
  ext  <- terra::extract(rgb, util, fun = mean, na.rm = TRUE, ID = TRUE)
  data.frame(parcela = util$parcela, R = ext$R, G = ext$G, B = ext$B)
}
# bandas <- extrair_bandas("ortomosaico_VT_R1.tif", "parcelas.gpkg", 0.5)
# write.csv(bandas, "bandas_parcelas.csv", row.names = FALSE)

# ---- 3. Índices RGB ---------------------------------------------------------
# Atenção: como r + g + b = 1, ExG = 3g - 1 e ExR, ExGR, TGI e CIVE também são
# combinações lineares de r e g. Eles carregam no máximo duas informações
# independentes; por isso a poda de redundância [REV3].
calc_indices <- function(d) {
  S <- d$R + d$G + d$B
  r <- d$R / S; g <- d$G / S; b <- d$B / S
  ExG <- 2 * g - r - b
  ExR <- 1.4 * r - g
  data.frame(
    ExG   = ExG,                                   # Woebbecke et al. 1995
    ExR   = ExR,                                   # Meyer & Neto 2008
    ExGR  = ExG - ExR,                             # Meyer & Neto 2008
    GLI   = (2 * g - r - b) / (2 * g + r + b),     # Louhaichi et al. 2001 (= VDVI)
    VARI  = (g - r) / (g + r - b),                 # Gitelson et al. 2002
    NGRDI = (g - r) / (g + r),                     # Tucker 1979; Hunt et al. 2005
    MGRVI = (g^2 - r^2) / (g^2 + r^2),             # Bendig et al. 2015
    RGBVI = (g^2 - b * r) / (g^2 + b * r),         # Bendig et al. 2015
    RGRI  = r / g,                                 # razão vermelho/verde
    TGI   = g - 0.39 * r - 0.61 * b,               # Hunt et al. 2011 (forma simplificada)
    CIVE  = 0.441 * r - 0.811 * g + 0.385 * b + 18.78745  # Kataoka et al. 2003
  )
}
VIS <- c("ExG", "ExR", "ExGR", "GLI", "VARI", "NGRDI",
         "MGRVI", "RGBVI", "RGRI", "TGI", "CIVE")

# ---- 4. Dados ---------------------------------------------------------------
trats <- c("BC0L0", "BC0L75", "BC0L100", "BC12L0", "BC12L75", "BC12L100")

if (USAR_DADOS_DEMO) {
  # Dados SIMULADOS só para testar o fluxo. Não interprete estes resultados.
  dados <- expand.grid(trat = trats, bloco = paste0("B", 1:4),
                       stringsAsFactors = FALSE)
  dados$parcela <- seq_len(nrow(dados))
  efeito <- c(0, 300, 450, 500, 800, 900)[match(dados$trat, trats)]
  dados$TLA <- 5200 + efeito + rnorm(24, 0, 350)            # cm2 planta-1
  G <- 120 + 0.012 * dados$TLA + rnorm(24, 0, 4)             # DN simulados
  dados$R <- 95 - 0.006 * dados$TLA + rnorm(24, 0, 3)
  dados$G <- G
  dados$B <- 70 - 0.003 * dados$TLA + rnorm(24, 0, 3)
  dados$rend <- 3000 + 0.9 * dados$TLA + rnorm(24, 0, 500)   # kg ha-1
} else {
  # Planilha com UMA linha por parcela:
  #   parcela, bloco, trat, TLA (média das 5 plantas, cm2), rend (kg ha-1)
  # Confira a numeração: o Livro2.xlsx e o matriz_consolidada.csv numeram as
  # parcelas de formas diferentes. A coluna 'parcela' precisa ser a mesma dos
  # polígonos.
  campo  <- read.csv("dados_campo_milho.csv")
  bandas <- read.csv("bandas_parcelas.csv")        # saída da seção 2
  dados  <- merge(campo, bandas, by = "parcela")
}

dados$LAI <- dados$TLA * DENSIDADE / 10000          # m2 m-2
dados <- cbind(dados, calc_indices(dados))
dados$bloco <- factor(dados$bloco)
dados$trat  <- factor(dados$trat, levels = trats)
stopifnot(nrow(dados) == 24)

# Com densidade constante, LAI = TLA x 0,001: r, R2 e rRMSE são idênticos.
# Para comparar com Du et al., converta RMSE e MAE para LAI (x DENSIDADE/1e4).
RESP <- "TLA"

# ---- 5. Métricas ------------------------------------------------------------
# [REV1] R2 de validação (coeficiente de eficiência). Pode ser negativo: aí o
# modelo erra mais do que simplesmente prever a média. r2_pearson fica só como
# referência: ele dá 1,0 para uma predição com erro sistemático de +400 cm2.
metricas <- function(obs, pred) {
  c(R2         = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2),
    r2_pearson = if (sd(pred) > 0) cor(obs, pred)^2 else NA_real_,
    RMSE       = sqrt(mean((obs - pred)^2)),
    rRMSE      = 100 * sqrt(mean((obs - pred)^2)) / mean(obs),
    MAE        = mean(abs(obs - pred)),
    vies       = mean(pred - obs))
}

# ---- 6. Correlação de Pearson (descritiva, todos os dados) -----------------
cor_tab <- data.frame(
  VI = VIS,
  r  = sapply(VIS, function(v) cor(dados[[v]], dados[[RESP]])),
  p  = sapply(VIS, function(v) cor.test(dados[[v]], dados[[RESP]])$p.value)
)
cor_tab <- cor_tab[order(-abs(cor_tab$r)), ]
print(cor_tab, digits = 3, row.names = FALSE)

# Redundância entre índices (ExG x GLI x CIVE têm |r| ~ 1)
print(round(cor(dados[, VIS]), 2))

# ---- 7. Modelos -------------------------------------------------------------
# A seleção de VIs é refeita DENTRO de cada treino. Selecionar com todos os
# 24 pontos e depois validar vazaria informação do teste (otimismo).
selecionar_vis <- function(treino, resp, limiar = LIMIAR_R, max_k = 3) {
  r <- sapply(VIS, function(v) abs(cor(treino[[v]], treino[[resp]])))
  cand <- names(sort(r, decreasing = TRUE))
  cand <- c(cand[r[cand] >= limiar], cand[1])[1:max(1, sum(r >= limiar))]
  sel <- character(0)                               # [REV3]
  for (v in cand) {
    if (!length(sel) || all(abs(cor(treino[[v]], treino[sel])) <= LIMIAR_REDUND)) sel <- c(sel, v)
    if (length(sel) == max_k) break
  }
  sel
}

# [REV4] RF com mtry escolhido pelo erro out-of-bag, só com o treino
ajustar_rf <- function(treino, vars, resp) {
  f <- reformulate(vars, resp)
  grade <- unique(pmax(1, round(length(vars) * c(1/3, 1/2, 2/3, 1))))
  oob <- sapply(grade, function(m)
    tail(randomForest(f, data = treino, ntree = N_ARVORES, mtry = m)$mse, 1))
  randomForest(f, data = treino, ntree = N_ARVORES, mtry = grade[which.min(oob)],
               importance = TRUE)
}

ajustar_prever <- function(treino, teste, modelo, resp = RESP) {
  if (modelo == "Nulo") return(rep(mean(treino[[resp]]), nrow(teste)))  # [REV2]
  if (modelo == "RF_todos") return(predict(ajustar_rf(treino, VIS, resp), teste))
  sel <- selecionar_vis(treino, resp)
  if (modelo == "ULR") return(predict(lm(reformulate(sel[1], resp), treino), teste))
  if (modelo == "MLR") return(predict(lm(reformulate(sel, resp), treino), teste))
  if (modelo == "RF")  return(predict(ajustar_rf(treino, sel, resp), teste))
  stop("modelo desconhecido")
}

modelos <- c("Nulo", "ULR", "MLR", "RF", "RF_todos")

# 7a. Validação principal: deixa um BLOCO inteiro de fora (4 folds, 18/6).
#     Respeita o DBC: cada teste traz os 6 tratamentos de um bloco não visto.
cv_blocos <- function(modelo, resp = RESP) {
  pred <- rep(NA_real_, nrow(dados))
  for (bl in levels(dados$bloco)) {
    te <- dados$bloco == bl
    pred[te] <- ajustar_prever(dados[!te, ], dados[te, ], modelo, resp)
  }
  pred
}

pred_bl <- sapply(modelos, cv_blocos)
res_bl  <- t(sapply(modelos, function(m) metricas(dados[[RESP]], pred_bl[, m])))
cat("\n== Validação deixando um bloco de fora (TLA, cm2 planta-1) ==\n")
print(round(res_bl, 3))

# 7b. Comparável a Du et al. (2022): 500 partições aleatórias 70/30
#     [REV5] as mesmas partições para todos os modelos
n <- nrow(dados)
particoes <- replicate(N_REP_7030, sample(n, round(0.7 * n)), simplify = FALSE)
cv_7030 <- function(modelo, resp = RESP) {
  t(sapply(particoes, function(tr) {
    p <- ajustar_prever(dados[tr, ], dados[-tr, ], modelo, resp)
    metricas(dados[[resp]][-tr], p)
  }))
}

res_7030 <- lapply(modelos, cv_7030)
names(res_7030) <- modelos
resumo_7030 <- do.call(rbind, lapply(modelos, function(m) {
  x <- res_7030[[m]]
  data.frame(modelo = m,
             R2_media   = mean(x[, "R2"]),
             R2_p05     = quantile(x[, "R2"], 0.05),
             R2_p95     = quantile(x[, "R2"], 0.95),
             RMSE_media = mean(x[, "RMSE"]),
             RMSE_p05   = quantile(x[, "RMSE"], 0.05),
             RMSE_p95   = quantile(x[, "RMSE"], 0.95),
             pct_RF_vence = if (m == "RF") NA else
               100 * mean(res_7030[["RF"]][, "RMSE"] < x[, "RMSE"]))
}))
cat("\n== 500 partições 70/30 (teste com 7 parcelas) ==\n")
cat("pct_RF_vence = % das partições em que o RF teve RMSE menor que o modelo da linha\n")
print(resumo_7030, digits = 3, row.names = FALSE)

# 7c. Cada índice sozinho, com validação por bloco (tabela do paper)
tab_vi <- do.call(rbind, lapply(VIS, function(v) {
  pred <- rep(NA_real_, nrow(dados))
  for (bl in levels(dados$bloco)) {
    te <- dados$bloco == bl
    pred[te] <- predict(lm(reformulate(v, RESP), dados[!te, ]), dados[te, ])
  }
  data.frame(VI = v, t(metricas(dados[[RESP]], pred)))
}))
tab_vi <- tab_vi[order(tab_vi$RMSE), ]
cat("\n== Cada VI isolado, validação por bloco ==\n")
print(tab_vi, digits = 3, row.names = FALSE)

# 7d. Importância das variáveis no RF_todos ajustado com todos os dados
#     [REV4] mesma especificação (todos os VIs, mtry por OOB) validada em 7a/7b
rf_final <- ajustar_rf(dados, VIS, RESP)
imp <- importance(rf_final, type = 1)                 # %IncMSE
imp <- data.frame(VI = rownames(imp), IncMSE = imp[, 1])
print(imp[order(-imp$IncMSE), ], row.names = FALSE, digits = 3)

# ---- 8. Área foliar -> rendimento de grãos ---------------------------------
cv_lm_bloco <- function(formula) {
  resp <- all.vars(formula)[1]
  pred <- rep(NA_real_, nrow(dados))
  for (bl in levels(dados$bloco)) {
    te <- dados$bloco == bl
    pred[te] <- predict(lm(formula, dados[!te, ]), dados[te, ])
  }
  metricas(dados[[resp]], pred)
}
# [REV6] VI escolhido dentro de cada treino, sem olhar o bloco de teste
cv_vi_rend <- function() {
  pred <- rep(NA_real_, nrow(dados))
  for (bl in levels(dados$bloco)) {
    te <- dados$bloco == bl
    v <- selecionar_vis(dados[!te, ], "rend", max_k = 1)
    pred[te] <- predict(lm(reformulate(v, "rend"), dados[!te, ]), dados[te, ])
  }
  metricas(dados$rend, pred)
}
dados$TLA_drone <- pred_bl[, "RF"]    # TLA estimada pelo drone (fora da amostra)

rend_tab <- rbind(
  TLA_medida  = cv_lm_bloco(rend ~ TLA),
  TLA_drone   = cv_lm_bloco(rend ~ TLA_drone),
  VI_direto   = cv_vi_rend()
)
cat("\n== Rendimento de grãos (kg ha-1), validação por bloco ==\n")
print(round(rend_tab, 3))

# ---- 9. Figuras -------------------------------------------------------------
tema <- theme_bw(base_size = 12)
mod_fig <- c("ULR", "MLR", "RF", "RF_todos")

# 9a. Observado x predito (validação por bloco)
df_pred <- tidyr::pivot_longer(
  data.frame(TLA = dados$TLA, bloco = dados$bloco, pred_bl[, mod_fig]),
  cols = all_of(mod_fig), names_to = "modelo", values_to = "pred")
df_pred$modelo <- factor(df_pred$modelo, levels = mod_fig)
rot <- data.frame(modelo = factor(mod_fig, levels = mod_fig),
                  lab = sprintf("R² = %.2f\nRMSE = %.0f cm²",
                                res_bl[mod_fig, "R2"], res_bl[mod_fig, "RMSE"]))
lim <- range(c(df_pred$TLA, df_pred$pred))

g1 <- ggplot(df_pred, aes(TLA, pred)) +
  geom_abline(linetype = 2, colour = "grey50") +
  geom_point(aes(shape = bloco), size = 2.4) +
  geom_text(data = rot, aes(x = -Inf, y = Inf, label = lab),
            hjust = -0.1, vjust = 1.2, size = 3.5) +
  facet_wrap(~ modelo, nrow = 1) + coord_equal(xlim = lim, ylim = lim) +
  labs(x = expression("TLA medida (cm"^2*" planta"^-1*")"),
       y = expression("TLA estimada (cm"^2*" planta"^-1*")"),
       shape = "Bloco") + tema
ggsave("fig_obs_pred_TLA.png", g1, width = 11, height = 3.8, dpi = 300)

# 9b. Distribuição do R2 nas 500 partições (como a Fig. 6 de Du et al.)
df_box <- do.call(rbind, lapply(modelos, function(m)
  data.frame(modelo = m, R2 = res_7030[[m]][, "R2"],
             RMSE = res_7030[[m]][, "RMSE"])))
df_box$modelo <- factor(df_box$modelo, levels = modelos)
g2 <- ggplot(df_box, aes(modelo, R2)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_boxplot(outlier.size = 0.6) +
  labs(x = NULL, y = "R² no teste (500 partições 70/30)") + tema
ggsave("fig_boxplot_R2.png", g2, width = 5.5, height = 4, dpi = 300)

# 9c. TLA medida x rendimento
g3 <- ggplot(dados, aes(TLA, rend)) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.6) +
  geom_point(aes(colour = trat), size = 2.4) +
  labs(x = expression("TLA (cm"^2*" planta"^-1*")"),
       y = expression("Rendimento de grãos (kg ha"^-1*")"),
       colour = "Tratamento") + tema
ggsave("fig_TLA_rendimento.png", g3, width = 6, height = 4.2, dpi = 300)

# ---- 10. Exporta tabelas ----------------------------------------------------
write.csv(cor_tab,     "tab_correlacao_VI_TLA.csv", row.names = FALSE)
write.csv(tab_vi,      "tab_VI_isolado_cv_bloco.csv", row.names = FALSE)
write.csv(data.frame(modelo = rownames(res_bl), res_bl),
          "tab_modelos_cv_bloco.csv", row.names = FALSE)
write.csv(resumo_7030, "tab_modelos_cv_7030.csv", row.names = FALSE)
write.csv(data.frame(alvo = rownames(rend_tab), rend_tab),
          "tab_rendimento_cv_bloco.csv", row.names = FALSE)
write.csv(imp,         "tab_importancia_RF.csv", row.names = FALSE)
