import numpy as np
import pandas as pd
import statsmodels.api as sm

def calc_risk_diff(exposed, n_exp, control, n_ctrl):
    p1 = exposed / n_exp
    p0 = control / n_ctrl
    rd = p1 - p0
    se = np.sqrt(p1 * (1 - p1) / n_exp + p0 * (1 - p0) / n_ctrl)
    return rd, rd - 1.96 * se, rd + 1.96 * se

COVS = [
    "age_years",
    "sex_group",
    "indication",
    "patient_source_group",
    "hp_model",
    "endoscopist_volume_group",
    "platform_group",
    "study_phase",
]

def design_matrix(work):
    parts = [pd.Series(1.0, index=work.index, name="const"), work[["g"]].astype(float)]
    for col in COVS:
        if str(work[col].dtype) in ("object", "category", "bool"):
            d = pd.get_dummies(work[col].astype("category"), prefix=col, drop_first=True)
            parts.append(d.astype(float))
        else:
            parts.append(work[[col]].astype(float))
    return pd.concat(parts, axis=1)

def poisson_cluster_fit(frame, outcome):
    cols = [outcome, "g", "platform_session_cluster"] + COVS
    work = frame[cols].copy().dropna()
    x = design_matrix(work)
    res = sm.GLM(work[outcome].astype(float), x, family=sm.families.Poisson()).fit(
        cov_type="cluster",
        cov_kwds={"groups": work["platform_session_cluster"], "use_correction": True, "df_correction": True},
    )
    return res, x

def adj_risk_diff(res, x):
    x1 = x.copy()
    x0 = x.copy()
    x1["g"] = 1.0
    x0["g"] = 0.0
    beta = np.asarray(res.params, dtype=float)
    cov = np.asarray(res.cov_params(), dtype=float)
    mu1 = np.exp(x1.to_numpy(dtype=float) @ beta)
    mu0 = np.exp(x0.to_numpy(dtype=float) @ beta)
    d = float(mu1.mean() - mu0.mean())
    grad = (x1.to_numpy(dtype=float) * mu1[:, None]).mean(0) - (x0.to_numpy(dtype=float) * mu0[:, None]).mean(0)
    se = float(np.sqrt(max(float(grad @ cov @ grad), 0.0)))
    return d, d - 1.96 * se, d + 1.96 * se


file_name = r""

df = pd.read_csv(file_name)
df["g"] = (df["cade"] == "A").astype(float)
a = df[df["cade"] == "A"]
b = df[df["cade"] == "M"]
print("A组人数", len(a))
print("B组人数", len(b))
print()

rate_list = [
    ("outcome_neoplasm", "胃肿瘤检出率"),
    ("outcome_atrophy_im", "肠化或萎缩性胃炎检出率"),
    ("outcome_polyp", "胃息肉检出率"),
    ("outcome_gastritis", "胃炎检出率"),
    ("biopsy_any", "活检率"),
    ("biopsy_rate_cardia_fundus", "贲门-胃底活检率"),
    ("biopsy_rate_body", "胃体活检率"),
    ("biopsy_rate_angulus", "胃角活检率"),
    ("biopsy_rate_antrum_pylorus", "胃窦-幽门活检率"),
    ("biopsy_rate_remnant", "残胃活检率"),
]

for col, label in rate_list:
    ac = int(a[col].sum())
    bc = int(b[col].sum())
    rd, lo, hi = calc_risk_diff(ac, len(a), bc, len(b))
    res, x = poisson_cluster_fit(df, col)
    adr, alo, ahi = adj_risk_diff(res, x)
    print(label)
    print("A组事件数", ac)
    print("B组事件数", bc)
    print("绝对风险差", rd)
    print("绝对风险差置信区间", lo, hi)
    print("聚类调整绝对风险差", adr)
    print("聚类调整绝对风险差置信区间", alo, ahi)
    print()
