import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

# 路径
cohort_path = r""
screen_path = r""

df = pd.read_csv(cohort_path)
sub = df[df["outcome_neoplasm"] == 0].copy()

xls = pd.ExcelFile(screen_path)
body_barcodes = set(pd.read_excel(xls, sheet_name="weiti")["Barcode"].dropna().astype(str))
target_body_barcodes = set(pd.read_excel(xls, sheet_name="weiti_baxiangbingbian")["Barcode"].dropna().astype(str))

sub["body_doc"] = sub["Barcode"].astype(str).isin(body_barcodes)
sub["body_target_doc"] = sub["Barcode"].astype(str).isin(target_body_barcodes)

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

def interaction_p_binary(flag):
    work = sub[["biopsy_rate_body", "cade", flag]].dropna().copy()
    work["cade"] = (work["cade"] == "A").astype(float)
    work["subgroup01"] = work[flag].astype(int)
    x = pd.DataFrame(
        {
            "const": 1.0,
            "cade": work["cade"],
            "subgroup01": work["subgroup01"].astype(float),
        }
    )
    x["interaction"] = x["cade"] * x["subgroup01"]
    result = sm.GLM(work["biopsy_rate_body"].astype(float), x, family=sm.families.Poisson()).fit(cov_type="HC1")
    return float(result.pvalues["interaction"])

def interaction_p_count(flag):
    work = sub[["biopsy_count_body", "cade", flag]].dropna().copy()
    work["cade"] = (work["cade"] == "A").astype(float)
    work["subgroup01"] = work[flag].astype(int)
    x = pd.DataFrame(
        {
            "const": 1.0,
            "cade": work["cade"],
            "subgroup01": work["subgroup01"].astype(float),
        }
    )
    x["interaction"] = x["cade"] * x["subgroup01"]
    result = sm.OLS(work["biopsy_count_body"].astype(float), x).fit(cov_type="HC1")
    return float(result.pvalues["interaction"])

def rate_block(use):
    a = use[use["cade"] == "A"]
    b = use[use["cade"] == "M"]
    ac = int(a["biopsy_rate_body"].sum())
    bc = int(b["biopsy_rate_body"].sum())
    rr, lo, hi = calc_rr(ac, len(a), bc, len(b))
    p = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
    print("A组事件数", ac)
    print("A组人数", len(a))
    print("M组事件数", bc)
    print("M组人数", len(b))
    print("RR", rr)
    print("RR置信区间", lo, hi)
    print("P值", p)
    print()

def count_block(use):
    av = use[use["cade"] == "A"]["biopsy_count_body"].fillna(0).astype(float)
    bv = use[use["cade"] == "M"]["biopsy_count_body"].fillna(0).astype(float)
    d, lo, hi = mean_diff_ci(av, bv)
    p = stats.mannwhitneyu(av, bv)[1]
    print("A组均值", av.mean())
    print("M组均值", bv.mean())
    print("差值", d)
    print("差值置信区间", lo, hi)
    print("P值", p)
    print()


print("胃体活检率")
print("整体")
rate_block(sub)
for flag, label in [("body_doc", "累及胃体"), ("body_target_doc", "目标病变累及胃体")]:
    print(label)
    print("交互P值", interaction_p_binary(flag))
    for level, level_label in [(True, "Yes"), (False, "No")]:
        print(level_label)
        rate_block(sub[sub[flag] == level])

print("胃体活检次数")
print("整体")
count_block(sub)
for flag, label in [("body_doc", "累及胃体"), ("body_target_doc", "目标病变累及胃体")]:
    print(label)
    print("交互P值", interaction_p_count(flag))
    for level, level_label in [(True, "Yes"), (False, "No")]:
        print(level_label)
        count_block(sub[sub[flag] == level])