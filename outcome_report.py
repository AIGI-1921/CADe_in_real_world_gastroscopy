import numpy as np
import pandas as pd
from scipy import stats

# 算 RR
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

file_name = r""

df = pd.read_csv(file_name)
a = df[df["cade"] == "A"]
b = df[df["cade"] == "M"]
print("A组人数", len(a))
print("B组人数", len(b))
print()

# 胃肿瘤检出率
ac = int(a["outcome_neoplasm"].sum())
bc = int(b["outcome_neoplasm"].sum())
rr, lo, hi = calc_rr(ac, len(a), bc, len(b))
p = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
print("胃肿瘤检出率")
print("A组事件数", ac)
print("A组事件百分比", ac / len(a) * 100)
print("B组事件数", bc)
print("B组事件百分比", bc / len(b) * 100)
print("RR", rr)
print("RR置信区间", lo, hi)
print("P值", p)
print()

# 肠化或萎缩性胃炎检出率
ac = int(a["outcome_atrophy_im"].sum())
bc = int(b["outcome_atrophy_im"].sum())
rr, lo, hi = calc_rr(ac, len(a), bc, len(b))
p = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
print("肠化或萎缩性胃炎检出率")
print("A组事件数", ac)
print("A组事件百分比", ac / len(a) * 100)
print("B组事件数", bc)
print("B组事件百分比", bc / len(b) * 100)
print("RR", rr)
print("RR置信区间", lo, hi)
print("P值", p)
print()

# 胃息肉检出率
ac = int(a["outcome_polyp"].sum())
bc = int(b["outcome_polyp"].sum())
rr, lo, hi = calc_rr(ac, len(a), bc, len(b))
p = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
print("胃息肉检出率")
print("A组事件数", ac)
print("A组事件百分比", ac / len(a) * 100)
print("B组事件数", bc)
print("B组事件百分比", bc / len(b) * 100)
print("RR", rr)
print("RR置信区间", lo, hi)
print("P值", p)
print()

# 胃炎检出率
ac = int(a["outcome_gastritis"].sum())
bc = int(b["outcome_gastritis"].sum())
rr, lo, hi = calc_rr(ac, len(a), bc, len(b))
p = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
print("胃炎检出率")
print("A组事件数", ac)
print("A组事件百分比", ac / len(a) * 100)
print("B组事件数", bc)
print("B组事件百分比", bc / len(b) * 100)
print("RR", rr)
print("RR置信区间", lo, hi)
print("P值", p)
print()

# 周转时间
ta = a[a["turnover_valid"] == 1]["turnover_time"]
tb = b[b["turnover_valid"] == 1]["turnover_time"]
p = stats.ttest_ind(ta, tb, equal_var=False)[1]
d, lo, hi = mean_diff_ci(ta, tb)
print("转诊时间")
print("A组均值", ta.mean())
print("A组标准差", ta.std())
print("B组均值", tb.mean())
print("B组标准差", tb.std())
print("差值", d)
print("差值置信区间", lo, hi)
print("P值", p)
print()

# 活检率
ac = int(a["biopsy_any"].sum())
bc = int(b["biopsy_any"].sum())
rr, lo, hi = calc_rr(ac, len(a), bc, len(b))
p = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
print("活检率")
print("A组事件数", ac)
print("A组事件百分比", ac / len(a) * 100)
print("B组事件数", bc)
print("B组事件百分比", bc / len(b) * 100)
print("RR", rr)
print("RR置信区间", lo, hi)
print("P值", p)
print()

# 贲门-胃底活检率
ac = int(a["biopsy_rate_cardia_fundus"].sum())
bc = int(b["biopsy_rate_cardia_fundus"].sum())
rr, lo, hi = calc_rr(ac, len(a), bc, len(b))
p = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
print("贲门-胃底活检率")
print("A组事件数", ac)
print("A组事件百分比", ac / len(a) * 100)
print("B组事件数", bc)
print("B组事件百分比", bc / len(b) * 100)
print("RR", rr)
print("RR置信区间", lo, hi)
print("P值", p)
print()

# 胃体活检率
ac = int(a["biopsy_rate_body"].sum())
bc = int(b["biopsy_rate_body"].sum())
rr, lo, hi = calc_rr(ac, len(a), bc, len(b))
p = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
print("胃体活检率")
print("A组事件数", ac)
print("A组事件百分比", ac / len(a) * 100)
print("B组事件数", bc)
print("B组事件百分比", bc / len(b) * 100)
print("RR", rr)
print("RR置信区间", lo, hi)
print("P值", p)
print()

# 胃角活检率
ac = int(a["biopsy_rate_angulus"].sum())
bc = int(b["biopsy_rate_angulus"].sum())
rr, lo, hi = calc_rr(ac, len(a), bc, len(b))
p = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
print("胃角活检率")
print("A组事件数", ac)
print("A组事件百分比", ac / len(a) * 100)
print("B组事件数", bc)
print("B组事件百分比", bc / len(b) * 100)
print("RR", rr)
print("RR置信区间", lo, hi)
print("P值", p)
print()

# 胃窦-幽门活检率
ac = int(a["biopsy_rate_antrum_pylorus"].sum())
bc = int(b["biopsy_rate_antrum_pylorus"].sum())
rr, lo, hi = calc_rr(ac, len(a), bc, len(b))
p = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
print("胃窦-幽门活检率")
print("A组事件数", ac)
print("A组事件百分比", ac / len(a) * 100)
print("B组事件数", bc)
print("B组事件百分比", bc / len(b) * 100)
print("RR", rr)
print("RR置信区间", lo, hi)
print("P值", p)
print()

# 残胃活检率
ac = int(a["biopsy_rate_remnant"].sum())
bc = int(b["biopsy_rate_remnant"].sum())
rr, lo, hi = calc_rr(ac, len(a), bc, len(b))
p = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
print("残胃活检率")
print("A组事件数", ac)
print("A组事件百分比", ac / len(a) * 100)
print("B组事件数", bc)
print("B组事件百分比", bc / len(b) * 100)
print("RR", rr)
print("RR置信区间", lo, hi)
print("P值", p)
print()

# 平均活检次数
av = a["biopsy_count_any"].fillna(0)
bv = b["biopsy_count_any"].fillna(0)
p = stats.mannwhitneyu(av, bv)[1]
d, lo, hi = mean_diff_ci(av, bv)
print("平均活检次数")
print("A组均值", av.mean())
print("B组均值", bv.mean())
print("差值", d)
print("差值置信区间", lo, hi)
print("P值", p)
print()

# 贲门-胃底平均活检次数
av = a["biopsy_count_cardia_fundus"].fillna(0)
bv = b["biopsy_count_cardia_fundus"].fillna(0)
p = stats.mannwhitneyu(av, bv)[1]
d, lo, hi = mean_diff_ci(av, bv)
print("贲门-胃底平均活检次数")
print("A组均值", av.mean())
print("B组均值", bv.mean())
print("差值", d)
print("差值置信区间", lo, hi)
print("P值", p)
print()

# 胃体平均活检次数
av = a["biopsy_count_body"].fillna(0)
bv = b["biopsy_count_body"].fillna(0)
p = stats.mannwhitneyu(av, bv)[1]
d, lo, hi = mean_diff_ci(av, bv)
print("胃体平均活检次数")
print("A组均值", av.mean())
print("B组均值", bv.mean())
print("差值", d)
print("差值置信区间", lo, hi)
print("P值", p)
print()

# 胃角平均活检次数
av = a["biopsy_count_angulus"].fillna(0)
bv = b["biopsy_count_angulus"].fillna(0)
p = stats.mannwhitneyu(av, bv)[1]
d, lo, hi = mean_diff_ci(av, bv)
print("胃角平均活检次数")
print("A组均值", av.mean())
print("B组均值", bv.mean())
print("差值", d)
print("差值置信区间", lo, hi)
print("P值", p)
print()

# 胃窦-幽门平均活检次数
av = a["biopsy_count_antrum_pylorus"].fillna(0)
bv = b["biopsy_count_antrum_pylorus"].fillna(0)
p = stats.mannwhitneyu(av, bv)[1]
d, lo, hi = mean_diff_ci(av, bv)
print("胃窦-幽门平均活检次数")
print("A组均值", av.mean())
print("B组均值", bv.mean())
print("差值", d)
print("差值置信区间", lo, hi)
print("P值", p)
print()

# 残胃平均活检次数
av = a["biopsy_count_remnant"].fillna(0)
bv = b["biopsy_count_remnant"].fillna(0)
p = stats.mannwhitneyu(av, bv)[1]
d, lo, hi = mean_diff_ci(av, bv)
print("残胃平均活检次数")
print("A组均值", av.mean())
print("B组均值", bv.mean())
print("差值", d)
print("差值置信区间", lo, hi)
print("P值", p)
