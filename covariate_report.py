import numpy as np
import pandas as pd
import statsmodels.api as sm

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


# 拼设计矩阵
def design_matrix(work):
    parts = [pd.Series(1, index=work.index, name="const"), work[["g"]].astype(float)]
    for col in COVS:
        if str(work[col].dtype) in ("object", "category", "bool"):
            d = pd.get_dummies(work[col].astype("category"), prefix=col, drop_first=True)
            parts.append(d.astype(float))
        else:
            parts.append(work[[col]].astype(float))
    return pd.concat(parts, axis=1)


# 二分类终点用修正泊松
def poisson_adj(frame, outcome):
    cols = [outcome, "g"] + COVS
    work = frame[cols].copy().dropna()
    x = design_matrix(work)
    res = sm.GLM(work[outcome].astype(float), x, family=sm.families.Poisson()).fit(cov_type="HC0")
    coef = float(res.params["g"])
    se = float(res.bse["g"])
    return (float(np.exp(coef)), float(np.exp(coef - 1.96 * se)), float(np.exp(coef + 1.96 * se)), float(res.pvalues["g"]))


# 连续终点用 OLS
def ols_adj(frame, outcome):
    cols = [outcome, "g"] + COVS
    work = frame[cols].copy().dropna()
    x = design_matrix(work)
    res = sm.OLS(work[outcome].astype(float), x).fit(cov_type="HC0")
    coef = float(res.params["g"])
    se = float(res.bse["g"])
    return (coef, coef - 1.96 * se, coef + 1.96 * se, float(res.pvalues["g"]))

file_name = r""

df = pd.read_csv(file_name)
a = df[df["cade"] == "A"]
b = df[df["cade"] == "M"]
df["g"] = (df["cade"] == "A").astype(float)
print("A组人数", len(a))
print("B组人数", len(b))
print()

# 胃肿瘤检出率
ac = int(a["outcome_neoplasm"].sum())
bc = int(b["outcome_neoplasm"].sum())
rr, lo, hi, p = poisson_adj(df, "outcome_neoplasm")
print("胃肿瘤检出率")
print("A组事件数", ac)
print("A组事件百分比", ac / len(a) * 100)
print("B组事件数", bc)
print("B组事件百分比", bc / len(b) * 100)
print("调整RR", rr)
print("调整RR置信区间", lo, hi)
print("调整P值", p)
print()

# 肠化或萎缩性胃炎检出率
ac = int(a["outcome_atrophy_im"].sum())
bc = int(b["outcome_atrophy_im"].sum())
rr, lo, hi, p = poisson_adj(df, "outcome_atrophy_im")
print("肠化或萎缩性胃炎检出率")
print("A组事件数", ac)
print("A组事件百分比", ac / len(a) * 100)
print("B组事件数", bc)
print("B组事件百分比", bc / len(b) * 100)
print("调整RR", rr)
print("调整RR置信区间", lo, hi)
print("调整P值", p)
print()

# 胃息肉检出率
ac = int(a["outcome_polyp"].sum())
bc = int(b["outcome_polyp"].sum())
rr, lo, hi, p = poisson_adj(df, "outcome_polyp")
print("胃息肉检出率")
print("A组事件数", ac)
print("A组事件百分比", ac / len(a) * 100)
print("B组事件数", bc)
print("B组事件百分比", bc / len(b) * 100)
print("调整RR", rr)
print("调整RR置信区间", lo, hi)
print("调整P值", p)
print()

# 胃炎检出率
ac = int(a["outcome_gastritis"].sum())
bc = int(b["outcome_gastritis"].sum())
rr, lo, hi, p = poisson_adj(df, "outcome_gastritis")
print("胃炎检出率")
print("A组事件数", ac)
print("A组事件百分比", ac / len(a) * 100)
print("B组事件数", bc)
print("B组事件百分比", bc / len(b) * 100)
print("调整RR", rr)
print("调整RR置信区间", lo, hi)
print("调整P值", p)
print()

# 周转时间
ta = a[a["turnover_valid"] == 1]["turnover_time"]
tb = b[b["turnover_valid"] == 1]["turnover_time"]
use = df[df["turnover_valid"] == 1].copy()
mean, lo, hi, p = ols_adj(use, "turnover_time")
print("周转时间")
print("A组均值", ta.mean())
print("A组标准差", ta.std())
print("B组均值", tb.mean())
print("B组标准差", tb.std())
print("调整均值差", mean)
print("调整均值差置信区间", lo, hi)
print("调整P值", p)
print()

# 活检率
ac = int(a["biopsy_any"].sum())
bc = int(b["biopsy_any"].sum())
rr, lo, hi, p = poisson_adj(df, "biopsy_any")
print("活检率")
print("A组事件数", ac)
print("A组事件百分比", ac / len(a) * 100)
print("B组事件数", bc)
print("B组事件百分比", bc / len(b) * 100)
print("调整RR", rr)
print("调整RR置信区间", lo, hi)
print("调整P值", p)
print()

# 贲门-胃底活检率
ac = int(a["biopsy_rate_cardia_fundus"].sum())
bc = int(b["biopsy_rate_cardia_fundus"].sum())
rr, lo, hi, p = poisson_adj(df, "biopsy_rate_cardia_fundus")
print("贲门-胃底活检率")
print("A组事件数", ac)
print("A组事件百分比", ac / len(a) * 100)
print("B组事件数", bc)
print("B组事件百分比", bc / len(b) * 100)
print("调整RR", rr)
print("调整RR置信区间", lo, hi)
print("调整P值", p)
print()

# 胃体活检率
ac = int(a["biopsy_rate_body"].sum())
bc = int(b["biopsy_rate_body"].sum())
rr, lo, hi, p = poisson_adj(df, "biopsy_rate_body")
print("胃体活检率")
print("A组事件数", ac)
print("A组事件百分比", ac / len(a) * 100)
print("B组事件数", bc)
print("B组事件百分比", bc / len(b) * 100)
print("调整RR", rr)
print("调整RR置信区间", lo, hi)
print("调整P值", p)
print()

# 胃角活检率
ac = int(a["biopsy_rate_angulus"].sum())
bc = int(b["biopsy_rate_angulus"].sum())
rr, lo, hi, p = poisson_adj(df, "biopsy_rate_angulus")
print("胃角活检率")
print("A组事件数", ac)
print("A组事件百分比", ac / len(a) * 100)
print("B组事件数", bc)
print("B组事件百分比", bc / len(b) * 100)
print("调整RR", rr)
print("调整RR置信区间", lo, hi)
print("调整P值", p)
print()

# 胃窦-幽门活检率
ac = int(a["biopsy_rate_antrum_pylorus"].sum())
bc = int(b["biopsy_rate_antrum_pylorus"].sum())
rr, lo, hi, p = poisson_adj(df, "biopsy_rate_antrum_pylorus")
print("胃窦-幽门活检率")
print("A组事件数", ac)
print("A组事件百分比", ac / len(a) * 100)
print("B组事件数", bc)
print("B组事件百分比", bc / len(b) * 100)
print("调整RR", rr)
print("调整RR置信区间", lo, hi)
print("调整P值", p)
print()

# 残胃活检率
ac = int(a["biopsy_rate_remnant"].sum())
bc = int(b["biopsy_rate_remnant"].sum())
rr, lo, hi, p = poisson_adj(df, "biopsy_rate_remnant")
print("残胃活检率")
print("A组事件数", ac)
print("A组事件百分比", ac / len(a) * 100)
print("B组事件数", bc)
print("B组事件百分比", bc / len(b) * 100)
print("调整RR", rr)
print("调整RR置信区间", lo, hi)
print("调整P值", p)
print()

# 平均活检次数
av = a["biopsy_count_any"].fillna(0)
bv = b["biopsy_count_any"].fillna(0)
mean, lo, hi, p = ols_adj(df, "biopsy_count_any")
print("平均活检次数")
print("A组均值", av.mean())
print("B组均值", bv.mean())
print("调整均值差", mean)
print("调整均值差置信区间", lo, hi)
print("调整P值", p)
print()

# 贲门-胃底平均活检次数
av = a["biopsy_count_cardia_fundus"].fillna(0)
bv = b["biopsy_count_cardia_fundus"].fillna(0)
mean, lo, hi, p = ols_adj(df, "biopsy_count_cardia_fundus")
print("贲门-胃底平均活检次数")
print("A组均值", av.mean())
print("B组均值", bv.mean())
print("调整均值差", mean)
print("调整均值差置信区间", lo, hi)
print("调整P值", p)
print()

# 胃体平均活检次数
av = a["biopsy_count_body"].fillna(0)
bv = b["biopsy_count_body"].fillna(0)
mean, lo, hi, p = ols_adj(df, "biopsy_count_body")
print("胃体平均活检次数")
print("A组均值", av.mean())
print("B组均值", bv.mean())
print("调整均值差", mean)
print("调整均值差置信区间", lo, hi)
print("调整P值", p)
print()

# 胃角平均活检次数
av = a["biopsy_count_angulus"].fillna(0)
bv = b["biopsy_count_angulus"].fillna(0)
mean, lo, hi, p = ols_adj(df, "biopsy_count_angulus")
print("胃角平均活检次数")
print("A组均值", av.mean())
print("B组均值", bv.mean())
print("调整均值差", mean)
print("调整均值差置信区间", lo, hi)
print("调整P值", p)
print()

# 胃窦-幽门平均活检次数
av = a["biopsy_count_antrum_pylorus"].fillna(0)
bv = b["biopsy_count_antrum_pylorus"].fillna(0)
mean, lo, hi, p = ols_adj(df, "biopsy_count_antrum_pylorus")
print("胃窦-幽门平均活检次数")
print("A组均值", av.mean())
print("B组均值", bv.mean())
print("调整均值差", mean)
print("调整均值差置信区间", lo, hi)
print("调整P值", p)
print()

# 残胃平均活检次数
av = a["biopsy_count_remnant"].fillna(0)
bv = b["biopsy_count_remnant"].fillna(0)
mean, lo, hi, p = ols_adj(df, "biopsy_count_remnant")
print("残胃平均活检次数")
print("A组均值", av.mean())
print("B组均值", bv.mean())
print("调整均值差", mean)
print("调整均值差置信区间", lo, hi)
print("调整P值", p)