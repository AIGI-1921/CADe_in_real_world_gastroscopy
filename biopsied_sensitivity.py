import numpy as np
import pandas as pd
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

def mean_diff_ci(x, y):
    d = x.mean() - y.mean()
    v1 = x.var()
    v0 = y.var()
    n1 = len(x)
    n0 = len(y)
    se2 = v1 / n1 + v0 / n0
    se = np.sqrt(se2)
    dfnum = se2 ** 2
    dfden = (v1 / n1) ** 2 / (n1 - 1) + (v0 / n0) ** 2 / (n0 - 1)
    df = dfnum / dfden if dfden > 0 else np.inf
    t = stats.t.ppf(0.975, df) if np.isfinite(df) else 1.96
    return d, d - t * se, d + t * se

cohort_path = r""

df = pd.read_csv(cohort_path)
biopsied = df[df["biopsy_any"] == 1].copy()

a = biopsied[biopsied["cade"] == "A"]
b = biopsied[biopsied["cade"] == "M"]
print("活检患者A组人数", len(a))
print("活检患者M组人数", len(b))
print()

rate_list = [
    ("outcome_neoplasm", "胃肿瘤检出率"),
    ("outcome_atrophy_im", "肠化或萎缩性胃炎检出率"),
    ("outcome_polyp", "胃息肉检出率"),
    ("outcome_gastritis", "胃炎检出率"),
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

# ============ 率终点 ============
for col, label in rate_list:
    ac = int(a[col].sum())
    bc = int(b[col].sum())
    rr, lo, hi = calc_rr(ac, len(a), bc, len(b))
    p = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
    print(label)
    print("A组事件数", ac)
    print("A组人数", len(a))
    print("M组事件数", bc)
    print("M组人数", len(b))
    print("RR", rr)
    print("RR置信区间", lo, hi)
    print("P值", p)
    print()

# ============ 计数终点 ============
for col, label in count_list:
    av = a[col].fillna(0).astype(float)
    bv = b[col].fillna(0).astype(float)
    d, lo, hi = mean_diff_ci(av, bv)
    p = stats.mannwhitneyu(av, bv)[1]
    print(label)
    print("A组均值", av.mean())
    print("M组均值", bv.mean())
    print("差值", d)
    print("差值置信区间", lo, hi)
    print("P值", p)
    print()

# ============ 周转时间（只看有效） ============
tt_a = a[a["turnover_valid"] == 1]
tt_b = b[b["turnover_valid"] == 1]
print("周转时间")
print("A组均值", tt_a["turnover_time"].mean())
print("M组均值", tt_b["turnover_time"].mean())
d, lo, hi = mean_diff_ci(tt_a["turnover_time"], tt_b["turnover_time"])
p = stats.ttest_ind(tt_a["turnover_time"], tt_b["turnover_time"], equal_var=False)[1]
print("差值", d)
print("差值置信区间", lo, hi)
print("P值", p)