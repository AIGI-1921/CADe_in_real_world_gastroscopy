import numpy as np
import pandas as pd
import statsmodels.api as sm
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

def poisson_cluster(frame, outcome):
    cols = [outcome, "g", "platform_session_cluster"] + COVS
    work = frame[cols].copy().dropna()
    x = design_matrix(work)
    res = sm.GLM(work[outcome].astype(float), x, family=sm.families.Poisson()).fit(
        cov_type="cluster",
        cov_kwds={"groups": work["platform_session_cluster"], "use_correction": True, "df_correction": True},
    )
    coef = float(res.params["g"])
    se = float(res.bse["g"])
    return (float(np.exp(coef)), float(np.exp(coef - 1.96 * se)), float(np.exp(coef + 1.96 * se)), float(res.pvalues["g"]))

cohort_path = r""
screen_path = r""

df = pd.read_csv(cohort_path)
df["g"] = (df["cade"] == "A").astype(float)
sub = df[df["outcome_neoplasm"] == 0].copy()

print("A组人数", len(sub[sub["cade"] == "A"]))
print("M组人数", len(sub[sub["cade"] == "M"]))
print()

xls = pd.ExcelFile(screen_path)

def barcode_set(sheet_name):
    return set(pd.read_excel(xls, sheet_name=sheet_name)["Barcode"].dropna().astype(str))

order = [
    ("总体病变记录率", "S3_screen_by_type"),
    ("贲门-胃底", "benmen_weidi"),
    ("胃体", "weiti"),
    ("胃角", "weijiao"),
    ("胃窦-幽门", "weidou_youmen"),
    ("残胃或胸胃", "canwei"),
    ("溃疡凹陷型", "kuiyang_aoxian"),
    ("糜烂充血型", "milan_fahong"),
    ("息肉隆起型", "xirou_longqi"),
    ("其他", "qita"),
]

for label, sheet in order:
    col = f"doc_{sheet}"
    sub[col] = sub["Barcode"].astype(str).isin(barcode_set(sheet)).astype(int)
    a = sub[sub["cade"] == "A"]
    b = sub[sub["cade"] == "M"]
    ac = int(a[col].sum())
    bc = int(b[col].sum())
    rr, lo, hi = calc_rr(ac, len(a), bc, len(b))
    p = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
    crr, clo, chi, cp = poisson_cluster(sub, col)
    print(label)
    print("A组记录数", ac)
    print("A组人数", len(a))
    print("M组记录数", bc)
    print("M组人数", len(b))
    print("RR", rr)
    print("RR置信区间", lo, hi)
    print("P值", p)
    print("聚类调整RR", crr)
    print("聚类调整RR置信区间", clo, chi)
    print("聚类P值", cp)
    print()