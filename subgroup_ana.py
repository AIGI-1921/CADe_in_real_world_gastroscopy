import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

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

a = df[df["cade"] == "A"]
b = df[df["cade"] == "M"]

hp_base = df[df["hp_model"] != "Missing"]
tt_valid = df[df["turnover_valid"] == 1]
tt_a = a[a["turnover_valid"] == 1]
tt_b = b[b["turnover_valid"] == 1]

print("a组人数", len(a))
print("b组人数", len(b))
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

# ============ 率终点 ============
for outcome, label in rate_list:
    print("{} - HP亚组".format(label))
    for sub_label in ["Negative_or_indeterminate", "Positive"]:
        a_sub = a[a["hp_model"] == sub_label]
        b_sub = b[b["hp_model"] == sub_label]
        ac = int(a_sub[outcome].sum())
        bc = int(b_sub[outcome].sum())
        rr, lo, hi = calc_rr(ac, len(a_sub), bc, len(b_sub))
        print(sub_label)
        print("a组事件数", ac)
        print("a组人数", len(a_sub))
        print("b组事件数", bc)
        print("b组人数", len(b_sub))
        print("RR", rr)
        print("RR置信区间", lo, hi)
        print()
    print("交互P值", interaction_p_binary(hp_base, outcome, "hp_model"))
    print()

    print("{} - 内镜医师量亚组".format(label))
    for sub_label in ["<1000", "≥1000"]:
        a_sub = a[a["endoscopist_volume_group"] == sub_label]
        b_sub = b[b["endoscopist_volume_group"] == sub_label]
        ac = int(a_sub[outcome].sum())
        bc = int(b_sub[outcome].sum())
        rr, lo, hi = calc_rr(ac, len(a_sub), bc, len(b_sub))
        print(sub_label)
        print("a组事件数", ac)
        print("a组人数", len(a_sub))
        print("b组事件数", bc)
        print("b组人数", len(b_sub))
        print("RR", rr)
        print("RR置信区间", lo, hi)
        print()
    print("交互P值", interaction_p_binary(df, outcome, "endoscopist_volume_group"))
    print()

# ============ 计数终点 ============
for outcome, label in count_list:
    print("{} - HP亚组".format(label))
    for sub_label in ["Negative_or_indeterminate", "Positive"]:
        a_sub = a[a["hp_model"] == sub_label][outcome].fillna(0).astype(float)
        b_sub = b[b["hp_model"] == sub_label][outcome].fillna(0).astype(float)
        print(sub_label)
        print("a组均值", a_sub.mean())
        print("b组均值", b_sub.mean())
        print()
    print("交互P值", interaction_p_continuous(hp_base, outcome, "hp_model"))
    print()

    print("{} - 内镜医师量亚组".format(label))
    for sub_label in ["<1000", "≥1000"]:
        a_sub = a[a["endoscopist_volume_group"] == sub_label][outcome].fillna(0).astype(float)
        b_sub = b[b["endoscopist_volume_group"] == sub_label][outcome].fillna(0).astype(float)
        print(sub_label)
        print("a组均值", a_sub.mean())
        print("b组均值", b_sub.mean())
        print()
    print("交互P值", interaction_p_continuous(df, outcome, "endoscopist_volume_group"))
    print()

# ============ 周转时间 ============
print("周转时间 - HP亚组")
for sub_label in ["Negative_or_indeterminate", "Positive"]:
    a_sub = tt_a[tt_a["hp_model"] == sub_label]["turnover_time"].astype(float)
    b_sub = tt_b[tt_b["hp_model"] == sub_label]["turnover_time"].astype(float)
    print(sub_label)
    print("a组均值", a_sub.mean())
    print("b组均值", b_sub.mean())
    print()

tt_hp = tt_valid[tt_valid["hp_model"] != "Missing"]
print("交互P值", interaction_p_continuous(tt_hp, "turnover_time", "hp_model"))
print()

print("周转时间 - 内镜医师量亚组")
for sub_label in ["<1000", "≥1000"]:
    a_sub = tt_a[tt_a["endoscopist_volume_group"] == sub_label]["turnover_time"].astype(float)
    b_sub = tt_b[tt_b["endoscopist_volume_group"] == sub_label]["turnover_time"].astype(float)
    print(sub_label)
    print("a组均值", a_sub.mean())
    print("b组均值", b_sub.mean())
    print()
print("交互P值", interaction_p_continuous(tt_valid, "turnover_time", "endoscopist_volume_group"))