"""Índices de vegetação no visível (RGB) usados por Du et al. (2022, Sci. Rep. 12:15937).

Todas as fórmulas operam sobre as coordenadas cromáticas normalizadas
r = R/(R+G+B), g = G/(R+G+B), b = B/(R+G+B), como em Du et al. (2022).
Para índices de razão (NGRDI, MGRVI, RGBVI, VARI, GLI...) a normalização não
altera o valor; para ExG, ExR, ExGR e CIVE ela altera a escala.

[VERIFICAR] Conferir cada fórmula contra a Tabela 2 de Du et al. (2022) antes da
submissão. A leitura automática da tabela deixou dúvida em RGRI e BI.
"""
import numpy as np

EPS = 1e-9


def chromatic(R, G, B):
    R, G, B = (np.asarray(x, dtype="float64") for x in (R, G, B))
    s = R + G + B
    s = np.where(s == 0, np.nan, s)
    return R / s, G / s, B / s


def _div(a, b):
    b = np.where(np.abs(b) < EPS, np.nan, b)
    return a / b


def compute_indices(R, G, B):
    """Recebe arrays de DN (0-255) e devolve dict {nome: array}."""
    r, g, b = chromatic(R, G, B)
    exg = 2 * g - r - b
    exr = 1.4 * r - g
    return {
        # Woebbecke et al. (1995)
        "ExG": exg,
        # Meyer & Neto (2008)
        "ExR": exr,
        "ExGR": exg - exr,
        # Tucker (1979); Hunt et al. (2005)
        "NGRDI": _div(g - r, g + r),
        # Verrelst et al. (2008). [VERIFICAR] Du et al. podem usar outra forma.
        "RGRI": _div(r, g),
        # Bendig et al. (2015)
        "MGRVI": _div(g**2 - r**2, g**2 + r**2),
        "RGBVI": _div(g**2 - b * r, g**2 + b * r),
        # Wang et al. (2015); algebricamente idêntico ao GLI (Louhaichi et al., 2001)
        "VDVI": _div(2 * g - r - b, 2 * g + r + b),
        # Gitelson et al. (2002)
        "VARI": _div(g - r, g + r - b),
        # Hague et al. (2006)
        "VEG": _div(g, np.power(r, 0.667) * np.power(b, 0.333)),
        # Kataoka et al. (2003)
        "CIVE": 0.441 * r - 0.811 * g + 0.385 * b + 18.78745,
        # Razão verde/vermelho
        "GR": _div(g, r),
        # Índice de brilho. [VERIFICAR] forma exata em Du et al. (2022)
        "BI": np.sqrt((r**2 + g**2 + b**2) / 3.0),
    }


# VDVI e GLI são a mesma expressão; mantemos só VDVI para não duplicar preditores.
INDEX_NAMES = list(compute_indices(np.ones(1), np.ones(1), np.ones(1)).keys())
