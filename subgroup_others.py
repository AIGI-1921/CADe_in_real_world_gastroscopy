import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

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


# 算二分类亚组的交互P值
def interaction_p_binary(frame, outcome, subgroup_var):
    work = frame[[outcome, "cade", subgroup_var]].dropna().copy()
    work["cade"] = pd.Categorical(work["cade"]).codes
    work["subgroup01"] = pd.Categorical(work[subgroup_var]).codes
    x = pd.DataFrame(
        {
            "const": 1.0,
            "cade": work["cade"].astype(float),
            "subgroup01": work["subgroup01"].astype(float),
        }
    )
    x["interaction"] = x["cade"] * x["subgroup01"]
    result = sm.GLM(work[outcome].astype(float), x, family=sm.families.Poisson()).fit(cov_type="HC1")
    return float(result.pvalues["interaction"])


# 算连续终点亚组的交互P值
def interaction_p_continuous(frame, outcome, subgroup_var):
    work = frame[[outcome, "cade", subgroup_var]].dropna().copy()
    work["cade"] = pd.Categorical(work["cade"]).codes
    work["subgroup01"] = pd.Categorical(work[subgroup_var]).codes
    x = pd.DataFrame(
        {
            "const": 1.0,
            "cade": work["cade"].astype(float),
            "subgroup01": work["subgroup01"].astype(float),
        }
    )
    x["interaction"] = x["cade"] * x["subgroup01"]
    result = sm.OLS(work[outcome].astype(float), x).fit(cov_type="HC1")
    return float(result.pvalues["interaction"])

file_name = r""

df = pd.read_csv(file_name)

df["age_group"] = np.where(df["age_years"] >= 40, "≥40", "<40")
df["fatigue_subgroup"] = np.where(df["fatigue_time"], "Yes", "No")
df["platform_subgroup"] = np.where(df["platform_code"].astype(str) == "1", "Platform 1", "Platform 2")

# 亚组定义
subgroups = [
    ("Age", "age_group", ["≥40", "<40"]),
    ("Sex", "sex_group", ["Male", "Female"]),
    ("Source", "patient_source_group", ["Gastroenterology", "Non-gastroenterology"]),
    ("Fatigue period", "fatigue_subgroup", ["Yes", "No"]),
    ("Platform", "platform_subgroup", ["Platform 1", "Platform 2"]),
    ("Study phase", "study_phase", ["Phase 1", "Phase 2"]),
]

a = df[df["cade"] == "A"]
b = df[df["cade"] == "M"]

# 周转时间
tt_valid = df[df["turnover_valid"] == 1]
tt_a = a[a["turnover_valid"] == 1]
tt_b = b[b["turnover_valid"] == 1]

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

print("a组人数", len(a))
print("b组人数", len(b))
print()

# ============ 率终点 ============
for outcome, label in rate_list:
    for sub_label, sub_var, levels in subgroups:
        print("{} - {}亚组".format(label, sub_label))
        for level in levels:
            a_sub = a[a[sub_var] == level]
            b_sub = b[b[sub_var] == level]
            ac = int(a_sub[outcome].sum())
            bc = int(b_sub[outcome].sum())
            rr, lo, hi = calc_rr(ac, len(a_sub), bc, len(b_sub))
            print(level)
            print("a组事件数", ac)
            print("a组人数", len(a_sub))
            print("b组事件数", bc)
            print("b组人数", len(b_sub))
            print("RR", rr)
            print("RR置信区间", lo, hi)
            print()
        print("交互P值", interaction_p_binary(df, outcome, sub_var))
        print()

# ============ 计数终点 ============
for outcome, label in count_list:
    for sub_label, sub_var, levels in subgroups:
        print("{} - {}亚组".format(label, sub_label))
        for level in levels:
            a_sub = a[a[sub_var] == level][outcome].fillna(0).astype(float)
            b_sub = b[b[sub_var] == level][outcome].fillna(0).astype(float)
            print(level)
            print("a组均值", a_sub.mean())
            print("b组均值", b_sub.mean())
            print()
        print("交互P值", interaction_p_continuous(df, outcome, sub_var))
        print()

# ============ 周转时间（只看有效） ============
for sub_label, sub_var, levels in subgroups:
    print("周转时间 - {}亚组".format(sub_label))
    for level in levels:
        a_sub = tt_a[tt_a[sub_var] == level]["turnover_time"].astype(float)
        b_sub = tt_b[tt_b[sub_var] == level]["turnover_time"].astype(float)
        print(level)
        print("a组均值", a_sub.mean())
        print("b组均值", b_sub.mean())
        print()
    print("交互P值", interaction_p_continuous(tt_valid, "turnover_time", sub_var))
    print()