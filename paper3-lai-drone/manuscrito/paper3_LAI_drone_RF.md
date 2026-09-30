# Estimativa do índice de área foliar do milho com imagens RGB de drone de baixo custo e Random Forest na Amazônia Central

*Título em inglês (proposta):* Estimating maize leaf area index from low-cost UAV RGB imagery: random forest versus conventional regression in Central Amazonia

**Autores:** Filipe Eduardo Danielli [PREENCHER coautores e afiliações]

**Periódico-alvo:** [PREENCHER — sugestões no README: *Drones*, *Precision Agriculture*, *Smart Agricultural Technology*]

> **Convenções deste rascunho** (seguem o CLAUDE.md do projeto):
> `[RESULTADO PENDENTE: ...]` marca número que só sai da análise; `[PREENCHER: ...]` marca dado de campo que só o autor tem;
> `[VERIFICAR: ...]` marca ponto a conferir na fonte; `[CITAÇÃO NECESSÁRIA: ...]` marca afirmação sem referência verificada.
> Nenhum resultado abaixo foi inventado. Sobre o grau de conferência das referências, ver a nota no início da lista.

---

## Resumo

O índice de área foliar (LAI) governa a interceptação de radiação, a transpiração e a produtividade do milho, mas medi-lo em campo é lento e trabalhoso. Veículos aéreos não tripulados (VANTs) com câmeras RGB de consumo oferecem uma alternativa barata, desde que índices de vegetação no visível consigam representar a variação do LAI. Este estudo avaliou se 13 índices RGB, extraídos de ortomosaicos com distância de amostragem no solo (GSD) de 0,84 cm obtidos com um DJI Mini 2, estimam o LAI do milho medido com o scanner a laser portátil CI-202 em 24 parcelas de um experimento fatorial de biochar de açaí e calcário na Amazônia Central. Comparamos regressão linear simples (ULR), regressão linear múltipla (MLR) e Random Forest (RF), seguindo o delineamento de Du et al. (2022), mas com seleção de variáveis e ajuste de hiperparâmetros feitos dentro de cada partição de validação cruzada, para evitar vazamento de informação em amostra pequena. [RESULTADO PENDENTE: melhor índice e seu r; R², RMSE e rRMSE de validação de ULR, MLR e RF; diferença RF − linear e sua significância pelo teste t corrigido]. [RESULTADO PENDENTE: frase-conclusão, que depende de o RF superar ou não a regressão linear com n = 24].

**Palavras-chave:** índice de área foliar; VANT; índices de vegetação no visível; Random Forest; validação cruzada; *Zea mays*; biochar.

---

## 1. Introdução

**O LAI é a variável-chave do dossel e medi-lo em campo é caro.** O LAI, definido como a área foliar verde de uma face por unidade de área de solo, controla a interceptação de luz, as trocas de água e carbono e, por consequência, o acúmulo de biomassa e a produtividade das culturas (Gower et al., 1999; Fang et al., 2019). Os métodos diretos, destrutivos ou com medidores de área foliar, são exatos mas lentos; os métodos ópticos indiretos aceleram a coleta, mas dependem de pressupostos sobre a distribuição das folhas e exigem calibração (Jonckheere et al., 2004; Yan et al., 2019). Em experimentos com muitas parcelas, o custo de medir o LAI limita a frequência de amostragem justamente quando os tratamentos mais divergem.

**VANTs com câmera RGB tornam o monitoramento do dossel acessível.** Sensores multiespectrais e hiperespectrais embarcados em VANTs estimam o LAI com boa acurácia, mas custam uma ordem de grandeza a mais que câmeras RGB de consumo [CITAÇÃO NECESSÁRIA: comparação de custo entre sensores RGB e multiespectrais em VANT]. Índices calculados só com as bandas do visível, como o ExG (Woebbecke et al., 1995), o NGRDI (Tucker, 1979; Hunt et al., 2005), o VARI (Gitelson et al., 2002), o GLI (Louhaichi et al., 2001) e os índices MGRVI e RGBVI (Bendig et al., 2015), respondem à fração de cobertura verde e à cor do dossel, duas propriedades ligadas ao LAI. Li et al. (2019) mostraram, em arroz, que índices de cor extraídos de imagens RGB de VANT estimam o LAI, e que combinar cor e textura melhora a estimativa.

**O aprendizado de máquina tem superado a regressão convencional, mas em bases de dados grandes.** O Random Forest (Breiman, 2001) combina muitas árvores de regressão treinadas em amostras *bootstrap*, lida bem com preditores colineares e capta relações não lineares, o que explica sua adoção ampla em sensoriamento remoto (Belgiu & Drăguţ, 2016). Para o LAI do milho, Du et al. (2022) compararam regressão linear simples e múltipla, rede neural e RF a partir de 14 índices RGB de imagens de VANT com GSD de 0,8 cm, usando 264 amostras de dois anos e quatro estádios; o RF obteve o melhor desempenho no teste (R² = 0,71; RMSE = 0,25 no período completo) [VERIFICAR: valores na Tabela de resultados de Du et al., 2022]. Magalhães & Rossi (2024) relataram R² = 0,98 com RF e índices RGB em milho para silagem, e observaram saturação dos índices RGB a partir de LAI ≈ 3 m² m⁻². Hussain et al. (2025) também encontraram o RF como melhor algoritmo para LAI de milho-doce a partir de índices multiespectrais de VANT.

**Duas lacunas motivam este trabalho.** A primeira é geográfica: os estudos acima vêm de sistemas temperados ou subtropicais, e não encontramos avaliação de índices RGB de VANT para o LAI do milho em solos altamente intemperizados da Amazônia Central [VERIFICAR: busca sistemática em Scopus/Web of Science antes da submissão]. A segunda é metodológica: com poucas dezenas de parcelas, típico de experimentos agronômicos, a vantagem do RF sobre a regressão linear não está garantida, e selecionar índices ou ajustar hiperparâmetros com os mesmos dados usados para validar infla a acurácia estimada (Varma & Simon, 2006). Du et al. (2022), por exemplo, escolheram os índices com r > 0,7 antes de separar treino e teste [VERIFICAR: se a seleção foi feita com o conjunto completo]. O experimento de biochar e calcário cria um gradiente de crescimento do milho que oferece variação de LAI em condições reais de manejo.

**Objetivo e hipóteses.** Avaliar a capacidade de índices RGB, derivados de imagens de um VANT de baixo custo com GSD submétrico (< 1 cm), de estimar o LAI do milho medido com o CI-202, comparando regressão convencional e Random Forest. Testamos duas hipóteses: (H1) ao menos um índice RGB se correlaciona significativamente com o LAI de campo; (H2) o RF estima o LAI com menor erro de validação cruzada que a melhor regressão linear.

---

## 2. Material e métodos

### 2.1 Área de estudo e delineamento experimental

O experimento foi conduzido na Estação Experimental de Fruticultura Tropical do Instituto Nacional de Pesquisas da Amazônia (INPA), km 41 da rodovia AM-010, Manaus, Amazonas (2°37′11,8″ S; 60°02′28,8″ W). [PREENCHER: tipo de solo (Latossolo Amarelo / Ferralsol), clima (Köppen), precipitação e temperatura médias do período].

O delineamento foi em blocos casualizados, em esquema fatorial 2 × 3, com quatro repetições: duas doses de biochar de resíduo de açaí (*Euterpe oleracea* Mart.) — 0 (BC0) e 12 t ha⁻¹ (BC12) — e três doses de calcário dolomítico — 0, 75 e 100% da dose recomendada (L0, L75, L100), totalizando 24 parcelas de 3 × 4 m (12 m²). O milho foi semeado no espaçamento de 0,80 m entre linhas e 0,25 m entre covas, com duas plantas por cova, o que resulta em 10 plantas m⁻². [PREENCHER: híbrido/cultivar, data de semeadura, adubação de base e cobertura]. O efeito dos tratamentos não é o objeto deste artigo: aqui eles servem para gerar variação de LAI entre parcelas.

### 2.2 LAI de referência com o CI-202

A área foliar foi medida de forma não destrutiva com o medidor portátil a laser CI-202 (CID Bio-Science, Camas, WA, EUA). Em cada parcela, todas as folhas de cinco plantas foram medidas, e a área foliar da parcela é a média dessas cinco plantas. A medição foi feita em 12 de dezembro de 2022, no mesmo dia do voo, o que elimina a defasagem entre a verdade de campo e a imagem. [VERIFICAR: ano da medição; o Filipe escreveu 12/12/2023, mas o voo e o processamento do Pix4D são de 12/12/2022]. [PREENCHER: dias após a semeadura; estádio fenológico; critério de escolha das plantas (ex.: aleatórias na área útil)]. O LAI de cada parcela foi calculado como

LAI = AF × D / 10 000,

em que AF é a área foliar média por planta (cm² planta⁻¹) e D é a densidade de plantas (plantas m⁻²), seguindo Du et al. (2022).



### 2.3 Aquisição e processamento das imagens

As imagens foram obtidas com um VANT DJI Mini 2 (DJI, Shenzhen, China), câmera FC7303 com sensor CMOS de 1/2,3″ e resolução de 4000 × 2250 pixels, em 12 de dezembro de 2022, o voo mais próximo da medição da área foliar. [VERIFICAR: o artigo 2 descreve outro voo, de 17/12/2022, com GSD de 1,7 cm; alinhar os dois textos para não parecer contradição entre os artigos]. [PREENCHER: altura de voo, horário, condição de céu, sobreposições frontal e lateral, modo de exposição/ISO, número de imagens]. As imagens foram processadas no Pix4D Mapper [PREENCHER: versão] para gerar o ortomosaico RGB, com GSD médio de 0,84 cm, segundo o relatório de processamento do Pix4D. [PREENCHER: uso ou não de pontos de controle; sistema de referência final]. [VERIFICAR: o relatório de processamento inicial disponível indica 11 de 25 imagens calibradas (44%); confirmar se houve reprocessamento com todas as imagens antes de extrair os índices.]

### 2.4 Extração dos índices de vegetação

Os polígonos das 24 parcelas foram desenhados sobre o ortomosaico e reduzidos por uma bordadura interna de 0,5 m para excluir efeito de borda e plantas das parcelas vizinhas. Em cada parcela, os valores digitais R, G e B de cada pixel foram convertidos em coordenadas cromáticas normalizadas, r = R/(R+G+B), g = G/(R+G+B) e b = B/(R+G+B), como em Du et al. (2022), e 13 índices foram calculados pixel a pixel (Tabela 1). Tomamos a média de cada índice (i) sobre todos os pixels da parcela, que incorpora a fração de solo exposto, e (ii) só sobre os pixels de vegetação, separados do solo pelo limiar de Otsu (1979) aplicado ao ExG. A fração de pixels de vegetação foi registrada como cobertura do dossel (CC). A abordagem (i), igual à de Du et al. (2022), é a análise principal; a (ii) é uma análise de sensibilidade.

**Tabela 1.** Índices de vegetação no visível usados neste estudo (r, g e b normalizados). [VERIFICAR: conferir cada fórmula contra a Tabela 2 de Du et al., 2022.]

| Índice | Fórmula | Referência |
|---|---|---|
| ExG | 2g − r − b | Woebbecke et al. (1995) |
| ExR | 1,4r − g | Meyer & Neto (2008) |
| ExGR | ExG − ExR | Meyer & Neto (2008) |
| NGRDI | (g − r)/(g + r) | Tucker (1979); Hunt et al. (2005) |
| RGRI | r/g | Verrelst et al. (2008) |
| MGRVI | (g² − r²)/(g² + r²) | Bendig et al. (2015) |
| RGBVI | (g² − b·r)/(g² + b·r) | Bendig et al. (2015) |
| VDVI (= GLI) | (2g − r − b)/(2g + r + b) | Louhaichi et al. (2001) [CITAÇÃO NECESSÁRIA: Wang et al. (2015) para o VDVI] |
| VARI | (g − r)/(g + r − b) | Gitelson et al. (2002) |
| VEG | g/(r^0,667 · b^0,333) | [CITAÇÃO NECESSÁRIA: Hague et al. (2006), Precision Agriculture] |
| CIVE | 0,441r − 0,811g + 0,385b + 18,78745 | [CITAÇÃO NECESSÁRIA: Kataoka et al. (2003)] |
| GR | g/r | [CITAÇÃO NECESSÁRIA: referência do GR usada por Du et al., 2022] |
| BI | √[(r² + g² + b²)/3] | [VERIFICAR: forma do BI em Du et al. (2022), que citam Camargo Neto (2004)] |
| CC | fração de pixels com ExG acima do limiar de Otsu | Otsu (1979) |

O VDVI e o GLI têm a mesma expressão algébrica; mantivemos apenas um deles para não duplicar preditores.

### 2.5 Modelos de estimativa

Comparamos cinco modelos, todos ajustados apenas com os dados de treino de cada partição:

1. **Nulo:** a média do LAI no treino, linha de base que qualquer modelo útil precisa superar.
2. **Regressão linear simples (ULR):** o índice com maior |r| no treino como único preditor, como em Du et al. (2022).
3. **Regressão linear múltipla (MLR):** índices com |r| > 0,7 no treino (critério de Du et al., 2022). Como os índices RGB são fortemente colineares, descartamos o índice com |r| > 0,95 com outro já selecionado de maior correlação com o LAI e limitamos o modelo a três preditores, cerca de um preditor para cada seis a sete observações de treino.
4. **RF com índices selecionados (RF_sel):** mesmos preditores do MLR.
5. **RF com todos os índices (RF_todos):** os 13 índices e a CC.

O RF foi implementado com o `RandomForestRegressor` do scikit-learn [PREENCHER: versão] com 500 árvores, como em Du et al. (2022). A fração de preditores sorteados em cada nó (*max_features*: 0,33; 0,5; 0,75; 1,0) e o número mínimo de observações por folha (*min_samples_leaf*: 1, 2, 3, 5) foram escolhidos pelo menor erro *out-of-bag* dentro do conjunto de treino, sem acesso aos dados de validação.

### 2.6 Validação e comparação dos modelos

Com 24 parcelas, um único teste 70/30 deixaria sete parcelas para validar, e a métrica dependeria fortemente de quais parcelas caíram no teste. Por isso usamos quatro esquemas: (i) validação cruzada *k*-fold com k = 5 repetida 100 vezes (esquema principal); (ii) *leave-one-out* (LOOCV); (iii) deixar um bloco inteiro de fora, que testa a transferência entre blocos com posições diferentes no campo; e (iv) 500 partições aleatórias 70/30, reproduzindo o esquema de Du et al. (2022) para permitir comparação direta.

Em cada repetição, calculamos o coeficiente de determinação de validação, R² = 1 − Σ(y − ŷ)²/Σ(y − ȳ)², a raiz do erro quadrático médio (RMSE), o RMSE relativo à média (rRMSE, %), o erro absoluto médio (MAE) e o viés médio (ŷ − y). Note que o R² de validação pode ser negativo quando o modelo erra mais que a média. Comparamos o RF com os modelos lineares pela diferença pareada de RMSE nas mesmas partições e pelo teste t corrigido para reamostragem de Nadeau & Bengio (2003), que leva em conta a sobreposição entre conjuntos de treino. A importância dos preditores no RF final foi estimada por permutação (50 repetições; queda no R²).

As análises foram feitas em Python [PREENCHER: versão] com scikit-learn, rasterio e geopandas; o código está em [PREENCHER: repositório/DOI Zenodo].

---

## 3. Resultados

> **Nota ao autor (remover antes da submissão):** as Seções 3.1–3.3 usam por enquanto só o ExG já extraído no QGIS (valores digitais 0–255, média por parcela), a pedido do autor. Saídas em `resultados/exg_qgis/`, geradas por `analise/R/lai_exg_qgis.R`. Quando os demais índices forem extraídos do ortomosaico, estas seções ganham as tabelas completas.

### 3.1 Variação do LAI e do ExG

O LAI medido com o CI-202 variou de 2,03 a 3,74 m² m⁻² (média 2,96; DP 0,34; CV 11,5%), e o ExG variou de 30,8 a 73,9 (média 53,1; DP 11,7). Os dois responderam aos tratamentos de forma diferente (ANOVA em blocos, esquema fatorial; `anova_fatorial.txt`). O LAI respondeu ao biochar (F = 15,48; p = 0,001), mas não ao calcário (F = 1,43; p = 0,27). O ExG respondeu principalmente ao calcário (F = 29,78; p < 0,001) e, em menor grau, ao biochar (F = 7,26; p = 0,017). O caso mais claro é o tratamento BC12L0: teve LAI médio de 3,15 m² m⁻², entre os maiores do experimento, mas ExG médio de 46,2, o segundo menor. Nenhuma interação biochar × calcário foi significativa para as duas variáveis (p ≥ 0,24).

### 3.2 Relação entre ExG e LAI

O ExG correlacionou-se positivamente com o LAI (r = 0,41; IC 95% 0,01–0,70; p = 0,044; Figura 1, `fig_ExG_LAI.png`). A regressão linear ajustada com as 24 parcelas foi LAI = 2,324 + 0,0121·ExG (R² = 0,17). O termo quadrático não foi significativo (p = 0,091). Com as médias dos seis tratamentos, a correlação foi r = 0,65 (p = 0,16; n = 6).

### 3.3 Desempenho dos modelos

Em validação cruzada, nenhum modelo baseado no ExG estimou o LAI melhor do que a média do conjunto de treino (Tabela 3, `tab_desempenho.csv`). Na validação principal, deixando um bloco inteiro de fora, o R² de validação foi −0,22 para a regressão linear, −0,15 para a quadrática, −0,25 para a exponencial e −0,30 para o Random Forest. O modelo nulo teve R² de −0,17, e os RMSE ficaram todos entre 0,36 e 0,38 m² m⁻² (rRMSE de 12,0% a 12,8%). No LOOCV e no 5-fold repetido 100 vezes, as regressões linear e quadrática ficaram marginalmente acima do nulo: RMSE de 0,334–0,344 contra 0,347–0,349 m² m⁻², com R² entre −0,07 e −0,01. O Random Forest foi o pior modelo em todos os esquemas: RMSE de 0,358 no LOOCV e 0,367 no 5-fold, e nas 500 partições 70/30 superou o modelo nulo em só 38% das partições, contra 58% da regressão linear (Figura 2, `fig_obs_pred_bloco.png`; Figura 3, `fig_boxplot_R2_7030.png`).

**Tabela 3.** Desempenho de validação na estimativa do LAI (m² m⁻²) a partir do ExG.

| Esquema | Modelo | R² | RMSE | rRMSE (%) |
|---|---|---|---|---|
| Deixa um bloco fora | Nulo | −0,17 | 0,360 | 12,1 |
| | Linear | −0,22 | 0,367 | 12,4 |
| | Quadrática | −0,15 | 0,357 | 12,0 |
| | Random Forest | −0,30 | 0,380 | 12,8 |
| LOOCV | Nulo | −0,09 | 0,347 | 11,7 |
| | Linear | −0,04 | 0,340 | 11,5 |
| | Quadrática | −0,01 | 0,334 | 11,3 |
| | Random Forest | −0,16 | 0,358 | 12,1 |
| 5-fold × 100 | Nulo | −0,10 | 0,349 | 11,8 |
| | Linear | −0,07 | 0,344 | 11,6 |
| | Quadrática | −0,05 | 0,341 | 11,5 |
| | Random Forest | −0,22 | 0,367 | 12,4 |

O ExG e o LAI correlacionaram-se com o rendimento de grãos com intensidade semelhante: r = 0,53 (p = 0,008) e r = 0,55 (p = 0,006), respectivamente (`rendimento.txt`).

### 3.4 Importância dos índices

[RESULTADO PENDENTE: Figura 4 (`fig_importancia_RF.png`) e frequência com que cada índice foi selecionado nas partições (`tab_frequencia_selecao_MLR.csv`)].

### 3.5 Análise de sensibilidade: só pixels de vegetação

[RESULTADO PENDENTE: repetir 3.2–3.3 com `--sufixo veg` e comparar com a análise principal].

---

## 4. Discussão

> Estrutura proposta; o conteúdo de cada parágrafo depende dos resultados.

### 4.1 Quais índices RGB carregam informação sobre o LAI

**Resultado com o ExG (a desenvolver):** o ExG mede sobretudo o verdor do dossel, e o verdor respondeu ao calcário, enquanto a área foliar respondeu ao biochar (Seção 3.1). Uma hipótese é que o calcário dolomítico, ao fornecer Ca e Mg e corrigir a acidez, aumentou o teor de clorofila e a cor verde das folhas sem aumentar na mesma proporção a área foliar [CITAÇÃO NECESSÁRIA: Mg, calagem e clorofila em milho; conferir com os teores foliares de Mg e o SPAD do próprio experimento]. Se essa hipótese se confirmar, a cor e a quantidade de folhas se desacoplam quando os tratamentos alteram a nutrição, e índices de verdor não servem como substitutos diretos do LAI nessas condições. Isso também explica por que o ExG prevê o rendimento tão bem quanto o LAI: o rendimento depende das duas coisas.


[RESULTADO PENDENTE]. Pontos a discutir: se os índices de razão verde-vermelho (NGRDI, MGRVI, VARI) se destacam, como no estádio de enchimento de grãos em Du et al. (2022) [VERIFICAR]; o papel da cobertura do dossel (CC) versus a cor das folhas; e se a média de todos os pixels (que mistura solo e planta) supera a média só da vegetação, o que indicaria que a fração de cobertura, e não a cor, conduz a relação.

### 4.2 Random Forest versus regressão convencional com amostra pequena

[RESULTADO PENDENTE]. Pontos a discutir: comparar com o ganho de R² de 0,64 para 0,71 do RF sobre a MLR em Du et al. (2022), obtido com 264 amostras [VERIFICAR]; se o RF não superar a regressão linear aqui, a explicação provável é o tamanho amostral e a baixa dimensionalidade efetiva dos índices RGB, que são quase redundantes entre si; o RF também não extrapola além da amplitude de treino, limitação relevante com LAI de amplitude estreita (Belgiu & Drăguţ, 2016). Discutir o quanto a seleção de variáveis fora da validação cruzada pode ter inflado acurácias publicadas (Varma & Simon, 2006).

### 4.3 Limitações

- **Amplitude do LAI.** A área foliar por planta variou pouco entre parcelas (CV de 11,5% nos dados preliminares). Com erro de medição do CI-202 e da amostragem de plantas da ordem dessa variação, a relação índice-LAI fica atenuada [CITAÇÃO NECESSÁRIA: erro de amostragem do LAI com poucas plantas por parcela].
- **Uma data.** Du et al. (2022) usaram quatro estádios em dois anos; aqui há um voo. Não é possível avaliar a transferência entre estádios.
- **Saturação.** Índices RGB tendem a saturar em LAI elevado (Magalhães & Rossi, 2024).
- **Radiometria.** Câmeras de consumo aplicam balanço de branco e compressão JPEG automáticos; sem painel de calibração, os valores digitais não são reflectância [CITAÇÃO NECESSÁRIA: efeito do processamento interno de câmeras RGB de consumo nos índices].
- **Escala.** A área foliar vem de poucas plantas por parcela; o índice, de toda a parcela.

### 4.4 Implicações práticas

[RESULTADO PENDENTE: o que o RMSE obtido significa para o manejo; custo do Mini 2 frente a sensores multiespectrais; potencial de acrescentar altura do dossel a partir do modelo digital de superfície do mesmo voo].

---

## 5. Conclusões

[RESULTADO PENDENTE: responder H1 e H2 em duas ou três frases, com os números principais].

---

## Referências

> Ano, periódico e DOI foram conferidos nesta sessão contra a página do artigo, a lista de referências de Du et al. (2022) ou uma busca bibliográfica. As entradas marcadas com [VERIFICAR] têm DOI, páginas ou autores ainda não conferidos na fonte. Antes da submissão, importe tudo por DOI no gerenciador de referências (Zotero/Mendeley) para fixar títulos e listas completas de autores.

- Belgiu, M., & Drăguţ, L. (2016). Random forest in remote sensing: A review of applications and future directions. *ISPRS Journal of Photogrammetry and Remote Sensing*, 114, 24–31. [VERIFICAR: DOI 10.1016/j.isprsjprs.2016.01.011]
- Bendig, J., Yu, K., Aasen, H., Bolten, A., Bennertz, S., Broscheit, J., Gnyp, M. L., & Bareth, G. (2015). Combining UAV-based plant height from crop surface models, visible, and near infrared vegetation indices for biomass monitoring in barley. *International Journal of Applied Earth Observation and Geoinformation*, 39, 79–87. [VERIFICAR: DOI e lista de autores]
- Breiman, L. (2001). Random forests. *Machine Learning*, 45, 5–32. https://doi.org/10.1023/A:1010933404324
- CID Bio-Science, Inc. (s.d.). *CI-202 portable laser leaf area meter: User's guide* (Rev. 6). Camas, WA.
- Du, L., Yang, H., Song, X., Wei, N., Yu, C., Wang, W., & Zhao, Y. (2022). Estimating leaf area index of maize using UAV-based digital imagery and machine learning methods. *Scientific Reports*, 12, 15937. https://doi.org/10.1038/s41598-022-20299-0
- Fang, H., Baret, F., Plummer, S., & Schaepman-Strub, G. (2019). An overview of global leaf area index (LAI): Methods, products, validation, and applications. *Reviews of Geophysics*, 57, 739–799. https://doi.org/10.1029/2018RG000608 [VERIFICAR: páginas]
- Gitelson, A. A., Kaufman, Y. J., Stark, R., & Rundquist, D. (2002). Novel algorithms for remote estimation of vegetation fraction. *Remote Sensing of Environment*, 80, 76–87. https://doi.org/10.1016/S0034-4257(01)00289-9
- Gower, S. T., Kucharik, C. J., & Norman, J. M. (1999). Direct and indirect estimation of leaf area index, fAPAR, and net primary production of terrestrial ecosystems. *Remote Sensing of Environment*, 70, 29–51. https://doi.org/10.1016/S0034-4257(99)00056-5
- Hunt, E. R., Cavigelli, M., Daughtry, C. S. T., McMurtrey, J. E., & Walthall, C. L. (2005). Evaluation of digital photography from model aircraft for remote sensing of crop biomass and nitrogen status. *Precision Agriculture*, 6, 359–378. https://doi.org/10.1007/s11119-005-2324-5
- Hussain, S., Teshome, F. T., Tulu, B. B., Awoke, G. W., Hailegnaw, N. S., & Bayabil, H. K. (2025). Leaf area index (LAI) prediction using machine learning and UAV based vegetation indices. *European Journal of Agronomy*, 168, 127557. https://doi.org/10.1016/j.eja.2025.127557
- Jonckheere, I., Fleck, S., Nackaerts, K., Muys, B., Coppin, P., Weiss, M., & Baret, F. (2004). Review of methods for in situ leaf area index determination: Part I. Theories, sensors and hemispherical photography. *Agricultural and Forest Meteorology*, 121, 19–35. [VERIFICAR: DOI 10.1016/j.agrformet.2003.08.027]
- Li, S., Yuan, F., Ata-Ul-Karim, S. T., Zheng, H., Cheng, T., Liu, X., Tian, Y., Zhu, Y., Cao, W., & Cao, Q. (2019). Combining color indices and textures of UAV-based digital imagery for rice LAI estimation. *Remote Sensing*, 11, 1763. https://doi.org/10.3390/rs11151763 [VERIFICAR: lista de autores]
- Louhaichi, M., Borman, M. M., & Johnson, D. E. (2001). Spatially located platform and aerial photography for documentation of grazing impacts on wheat. *Geocarto International*, 16(1), 65–70. https://doi.org/10.1080/10106040108542184
- Magalhães, L. P., & Rossi, F. (2024). Use of indices in RGB and random forest regression to measure the leaf area index in maize. *Agronomy*, 14, 750. https://doi.org/10.3390/agronomy14040750
- Meyer, G. E., & Neto, J. C. (2008). Verification of color vegetation indices for automated crop imaging applications. *Computers and Electronics in Agriculture*, 63, 282–293. https://doi.org/10.1016/j.compag.2008.03.009
- Nadeau, C., & Bengio, Y. (2003). Inference for the generalization error. *Machine Learning*, 52, 239–281. [VERIFICAR: DOI 10.1023/A:1024068626366]
- Otsu, N. (1979). A threshold selection method from gray-level histograms. *IEEE Transactions on Systems, Man, and Cybernetics*, 9(1), 62–66. https://doi.org/10.1109/TSMC.1979.4310076
- Tucker, C. J. (1979). Red and photographic infrared linear combinations for monitoring vegetation. *Remote Sensing of Environment*, 8, 127–150. https://doi.org/10.1016/0034-4257(79)90013-0
- Varma, S., & Simon, R. (2006). Bias in error estimation when using cross-validation for model selection. *BMC Bioinformatics*, 7, 91. https://doi.org/10.1186/1471-2105-7-91
- Verrelst, J., Schaepman, M. E., Koetz, B., & Kneubühler, M. (2008). Angular sensitivity analysis of vegetation indices derived from CHRIS/PROBA data. *Remote Sensing of Environment*, 112, 2341–2353. https://doi.org/10.1016/j.rse.2007.11.001
- Woebbecke, D. M., Meyer, G. E., Von Bargen, K., & Mortensen, D. A. (1995). Color indices for weed identification under various soil, residue, and lighting conditions. *Transactions of the ASAE*, 38(1), 259–269. https://doi.org/10.13031/2013.27838
- Yan, G., Hu, R., Luo, J., Weiss, M., Jiang, H., Mu, X., Xie, D., & Zhang, W. (2019). Review of indirect optical measurements of leaf area index: Recent advances, challenges, and perspectives. *Agricultural and Forest Meteorology*, 265, 390–411. [VERIFICAR: DOI 10.1016/j.agrformet.2018.11.033 e páginas]

**Leituras complementares citadas por Du et al. (2022), úteis para a discussão** (DOIs obtidos da lista de referências de Du et al.; conferir antes de citar): Houborg & McCabe (2018), *ISPRS J. Photogramm.*, 10.1016/j.isprsjprs.2017.10.004 (LAI com Cubist e RF); Han et al. (2019), *Plant Methods*, 10.1186/s13007-019-0394-z (biomassa de milho com aprendizado de máquina); Maimaitijiang et al. (2020), *Remote Sensing*, 10.3390/rs12091357; Osco et al. (2020), *Remote Sensing*, 10.3390/rs12193237 (milho, VANT, aprendizado de máquina, Brasil); Yue et al. (2018), *Remote Sensing*, 10.3390/rs10071138; Rodriguez-Galiano et al. (2012), *ISPRS J. Photogramm.*, 10.1016/j.isprsjprs.2011.11.002.
