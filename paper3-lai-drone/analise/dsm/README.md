# DSM, altura e grade de parcelas (voos 12/12 e 17/12/2022)

Scripts usados no relatório `relatorio_DSM_altura_12x17` (pasta compartilhada PAPER3-LAI-DRONE/dsm_altura).

- `detect.py`: densidade de vegetação (ExG) reamostrada do ortomosaico.
- `gridfit.py`: ajusta a grade 4 × 6 das parcelas (centro, rotação, passo) por otimização sobre a densidade; salva `<tag>_grid.npy`.
- `overlay.py`: desenha a grade sobre o ortomosaico para conferência.
- `extract.py <tag> <ortomosaico> <dsm>`: por parcela, índices RGB, cobertura (Otsu do ExG) e altura (DSM − chão local, percentil 2 dos pixels de solo).
- `simulacao.py`: regressão linear × Random Forest com n = 24 (LOOCV).

Ortomosaico e DSM precisam ser os originais do Pix4D (mesma grade de pixels).
