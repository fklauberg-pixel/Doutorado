# Paper 3 — LAI do milho por drone RGB: Random Forest vs regressão convencional

Base metodológica: Du et al. (2022), *Estimating leaf area index of maize using UAV-based
digital imagery and machine learning methods*, Scientific Reports 12:15937
(https://doi.org/10.1038/s41598-022-20299-0).

Verdade de campo: área foliar com o scanner a laser CI-202. Imagens: DJI Mini 2, GSD 0,84 cm
(relatório Pix4D de 12/12/2022). 24 parcelas, fatorial biochar (BC0, BC12) × calcário (L0, L75, L100), 4 blocos.

## Estrutura

```
paper3-lai-drone/
├── config.yaml                      # caminhos e parâmetros (editar aqui)
├── analise/
│   ├── indices.py                   # fórmulas dos 13 índices RGB
│   ├── 00_preliminar_exg.py         # ExG já existente x área foliar (roda hoje)
│   ├── 01_extrair_indices.py        # ortomosaico + parcelas -> índices por parcela
│   ├── 02_modelagem.py              # ULR, MLR, Random Forest, validação cruzada, figuras
│   └── teste_sintetico.py           # testa o pipeline com dados falsos
├── dados/
│   ├── preliminar_exg_tla.csv       # extraído de Livro2.xlsx (conferido com EXG.csv)
│   └── lai_ci202_MODELO.csv         # formato esperado do CSV de campo
├── manuscrito/paper3_LAI_drone_RF.md
└── resultados/preliminar/           # saída do 00_preliminar_exg.py
```

## Arquivos que faltam para rodar a análise completa

Coloque em `dados/` (rasters e vetores ficam fora do git, ver `.gitignore`):

1. **`ortomosaico_rgb.tif`** — ortomosaico RGB georreferenciado exportado do Pix4D
   (`3_dsm_ortho/2_mosaic/*_transparent_mosaic_group1.tif`). No material compartilhado só
   existe a etapa 1 (`1_initial`) e o `orthomosaic_preview.png`, que não serve para extrair
   índices. O relatório mostra só 11 de 25 imagens calibradas: vale reprocessar antes.
2. **`parcelas.gpkg`** (ou `.shp`) — polígonos das 24 parcelas com a coluna `Parcela`
   (mesma numeração do Livro2.xlsx), em CRS métrico (ex.: SIRGAS 2000 / UTM 20S, EPSG:31980).
   O `B.gpkg` compartilhado tem um único polígono.
3. **`lai_ci202.csv`** — uma linha por parcela, como em `lai_ci202_MODELO.csv`. A coluna
   `TLA_cm2_planta` já existe no Livro2.xlsx; faltam o número de plantas medidas e a data.
4. **Densidade de plantas** (plantas m⁻²) em `config.yaml`, para converter área foliar
   por planta em LAI. Sem ela o pipeline roda com cm² planta⁻¹ (r e R² não mudam).
5. Para o texto de Métodos: híbrido, data de semeadura, espaçamento, data da medição com o
   CI-202, altura e horário do voo. Estão marcados como `[PREENCHER]` no manuscrito.

## Como rodar

```bash
pip install -r requirements.txt
python analise/00_preliminar_exg.py                 # já roda com os dados atuais
python analise/01_extrair_indices.py                # precisa de 1, 2
python analise/02_modelagem.py                      # precisa de 3; ~30 min com 4 núcleos
python analise/02_modelagem.py --sufixo veg         # sensibilidade: só pixels de vegetação
python analise/teste_sintetico.py /tmp/teste_lai    # verificação do código com dados falsos
```

Use `--rapido` no `02_modelagem.py` para um teste de ~1 minuto.

## Decisões de método (diferenças em relação a Du et al., 2022)

- **Seleção de índices e ajuste do RF dentro de cada partição de treino.** Du et al.
  selecionaram índices com r > 0,7 antes de dividir os dados; com n = 24 isso infla a
  acurácia (Varma & Simon, 2006).
- **Validação cruzada repetida em vez de só 70/30.** Com 24 parcelas, um teste 70/30 tem 7
  parcelas. O esquema 70/30 × 500 de Du et al. também roda, para comparação direta.
- **Deixar um bloco fora** como teste extra de transferência espacial.
- **Modelo nulo** (média) como linha de base: R² de validação negativo significa pior que a média.
- **Teste t corrigido** (Nadeau & Bengio, 2003) para dizer se a diferença RF − linear é real.
- **Cobertura do dossel (CC)** acrescentada aos 13 índices.

## Resultado preliminar (ExG já calculado no QGIS × área foliar do CI-202)

Ver `resultados/preliminar/preliminar_exg_tla.txt`. Resumo: r = 0,41 (IC 95% 0,01–0,70),
R² de ajuste = 0,17, mas R² de validação LOOCV = −0,04. O ExG sozinho, calculado sobre a
parcela inteira, não prevê a área foliar melhor que a média. A área foliar por planta
varia pouco entre parcelas (CV = 11,5%), o que limita qualquer modelo.

## Periódicos sugeridos

*Drones* (MDPI), *Precision Agriculture* (Springer), *Smart Agricultural Technology*
(Elsevier). Com n = 24 e uma data, *Drones* ou *Smart Agricultural Technology* são mais realistas;
*Precision Agriculture* pede mais datas ou mais parcelas.
