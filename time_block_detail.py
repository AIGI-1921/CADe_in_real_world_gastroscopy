import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

def calc_rr(exposed, n_exp, control, n_ctrl):
    if exposed == 0 or control == 0:
        exposed += 0.5
        control += 0.5
        n_exp += 1
        n_ctrl += 1
    rr = (exposed / n_exp) / (control / n_ctrl)
    se = np.sqrt(1 / exposed - 1 / n_exp + 1 / control - 1 / n_ctrl)
    lo = np.exp(np.log(rr) - 1.96 * se)
    hi = np.exp(np.log(rr) + 1.96 * se)
    return rr, lo, hi

COVS = [
    "age_years",
    "sex_group",
    "indication",
    "patient_source_group",
    "hp_model",
    "endoscopist_volume_group",
    "platform_group",
]

A_COVS = [
    "age_years",
    "sex_group",
    "indication",
    "patient_source_group",
    "hp_model",
    "endoscopist_volume_group",
    "platform_group",
    "study_phase",
]


def design_matrix(work, covs, extra=None):
    parts = [pd.Series(1.0, index=work.index, name="const"), work[["g"]].astype(float)]
    for col in covs:
        if str(work[col].dtype) in ("object", "category", "bool"):
            d = pd.get_dummies(work[col].astype("category"), prefix=col, drop_first=True)
            parts.append(d.astype(float))
        else:
            parts.append(work[[col]].astype(float))
    if extra is not None:
        parts.append(pd.get_dummies(work[extra].astype("category"), prefix=extra, drop_first=True).astype(float))
    return pd.concat(parts, axis=1)


def poisson_cluster(frame, outcome, covs, extra=None):
    cols = [outcome, "g", "cluster"] + covs
    if extra is not None:
        cols.append(extra)
    work = frame[cols].copy().dropna()
    x = design_matrix(work, covs, extra)
    res = sm.GLM(work[outcome].astype(float), x, family=sm.families.Poisson()).fit(
        cov_type="cluster",
        cov_kwds={"groups": work["cluster"], "use_correction": True, "df_correction": True},
    )
    coef = float(res.params["g"])
    se = float(res.bse["g"])
    return (float(np.exp(coef)), float(np.exp(coef - 1.96 * se)), float(np.exp(coef + 1.96 * se)), float(res.pvalues["g"]))

def linear_cluster(frame, outcome, covs, extra=None):
    cols = [outcome, "g", "cluster"] + covs
    if extra is not None:
        cols.append(extra)
    work = frame[cols].copy().dropna()
    x = design_matrix(work, covs, extra)
    res = sm.OLS(work[outcome].astype(float), x).fit(
        cov_type="cluster",
        cov_kwds={"groups": work["cluster"], "use_correction": True, "df_correction": True},
    )
    coef = float(res.params["g"])
    se = float(res.bse["g"])
    return (coef, coef - 1.96 * se, coef + 1.96 * se, float(res.pvalues["g"]))

file_name = r""

df = pd.read_csv(file_name)
df["cluster"] = df["platform_session_cluster"]
df["g"] = (df["cade"] == "A").astype(float)
dates = sorted(pd.to_datetime(df["report_date"]).dropna().unique())
blocks = np.array_split(np.array(dates), 10)
mapping = {}
for i, block in enumerate(blocks, start=1):
    for day in block:
        mapping[pd.Timestamp(day)] = f"Time block {i}"
df["time_block"] = pd.to_datetime(df["report_date"]).map(mapping)

print("时间块分布")
for label, sub in df.groupby("time_block", sort=False):
    days = pd.to_datetime(sub["report_date"]).sort_values().unique()
    print(label, "日期", pd.Timestamp(days[0]).date(), "到", pd.Timestamp(days[-1]).date(), "总人数", len(sub), "A组", int((sub["cade"] == "A").sum()), "M组", int((sub["cade"] == "M").sum()))
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

count_list = [
    ("biopsy_count_any", "平均活检次数"),
    ("biopsy_count_cardia_fundus", "贲门-胃底平均活检次数"),
    ("biopsy_count_body", "胃体平均活检次数"),
    ("biopsy_count_angulus", "胃角平均活检次数"),
    ("biopsy_count_antrum_pylorus", "胃窦-幽门平均活检次数"),
    ("biopsy_count_remnant", "残胃平均活检次数"),
]

tt = df[df["turnover_valid"] == 1]

print("时间块调整模型（聚类稳健）")
# ============ 率终点 ============
for col, label in rate_list:
    rr, lo, hi, p = poisson_cluster(df, col, COVS, "time_block")
    print(label)
    print("调整RR", rr)
    print("调整RR置信区间", lo, hi)
    print("调整P值", p)
    print()

# ============ 计数终点 ============
for col, label in count_list:
    d, lo, hi, p = linear_cluster(df, col, COVS, "time_block")
    print(label)
    print("调整差值", d)
    print("调整差值置信区间", lo, hi)
    print("调整P值", p)
    print()

# ============ 周转时间 ============
d, lo, hi, p = linear_cluster(tt, "turnover_time", COVS, "time_block")
print("周转时间")
print("调整差值", d)
print("调整差值置信区间", lo, hi)
print("调整P值", p)
print()

# ============ 逐个排除时间块 ============
def loo_block(frame, outcome):
    for label in list(dict.fromkeys(frame["time_block"].tolist())):
        sub = frame[frame["time_block"] != label].copy()
        ai = sub[sub["cade"] == "A"]
        ctrl = sub[sub["cade"] == "M"]
        ac = int(ai[outcome].sum())
        bc = int(ctrl[outcome].sum())
        rr, lo, hi = calc_rr(ac, len(ai), bc, len(ctrl))
        p = stats.chi2_contingency([[ac, len(ai) - ac], [bc, len(ctrl) - bc]], correction=False)[1]
        crr, clo, chi, cp = poisson_cluster(sub, outcome, A_COVS)
        print("排除" + label)
        print("RR", rr)
        print("RR置信区间", lo, hi)
        print("P值", p)
        print("聚类RR", crr)
        print("聚类RR置信区间", clo, chi)
        print("聚类P值", cp)
        print()


print("逐个排除时间块：胃肿瘤检出率")
loo_block(df, "outcome_neoplasm")
print("逐个排除时间块：胃息肉检出率")
loo_block(df, "outcome_polyp")