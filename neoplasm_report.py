import numpy as np
import pandas as pd
from scipy import stats

cancer_path = r""
cohort_path = r""

# 病理类型映射
PATHOLOGY_MAP = {
    "腺癌": "Adenocarcinoma",
    "低级别异型增生": "LGIN",
    "高级别异型增生": "HGIN",
    "腺瘤": "Adenoma",
    "神经内分泌瘤": "NET",
    "淋巴瘤": "Lymphoma",
    "间质瘤": "GIST",
    "转移癌": "Metastasis",
}


# 解析病变大小成类别
def parse_size(value):
    if pd.isna(value) or str(value).strip() == "":
        return "NA"
    text = str(value).strip()
    if text == "大":
        return ">20 mm"
    try:
        size_mm = float(text)
    except ValueError:
        return "NA"
    if size_mm <= 5:
        return "≤5 mm"
    if size_mm <= 10:
        return "5-10 mm"
    if size_mm <= 20:
        return "10-20 mm"
    return ">20 mm"


# 二分类比较
def chi_or_fisher(table):
    try:
        return float(stats.chi2_contingency(table, correction=False)[1])
    except ValueError:
        if table.shape == (2, 2):
            return float(stats.fisher_exact(table)[1])
        return float("nan")

def ffh_p(table, n_sim=200000, seed=20260725):
    row_sums = table.sum(axis=1)
    col_sums = table.sum(axis=0)
    keep_r = row_sums > 0
    keep_c = col_sums > 0
    table = table[np.ix_(keep_r, keep_c)]
    row_sums = table.sum(axis=1)
    col_sums = table.sum(axis=0)
    expected = np.outer(row_sums, col_sums) / table.sum()
    observed = float((((table - expected) ** 2) / expected).sum())
    rng = np.random.default_rng(seed)
    sims = stats.random_table(row_sums, col_sums).rvs(size=n_sim, random_state=rng)
    sim_stats = ((sims - expected) ** 2 / expected).sum(axis=(1, 2))
    return float((np.sum(sim_stats >= observed) + 1) / (n_sim + 1))

def sparse_multicat(table):
    row_sums = table.sum(axis=1)
    col_sums = table.sum(axis=0)
    table = table[np.ix_(row_sums > 0, col_sums > 0)]
    try:
        _, p, _, expected = stats.chi2_contingency(table, correction=False)
    except ValueError:
        return ffh_p(table)
    if float(expected.min()) < 5:
        return ffh_p(table)
    return float(p)

cancer = pd.read_excel(cancer_path)

cohort = pd.read_csv(cohort_path, usecols=["Barcode", "cade"])
cancer = cancer.merge(cohort, on="Barcode", how="left")

lesions = pd.DataFrame(
    {
        "group": cancer["cade"],
        "pathology": cancer["病变1病理类型"].map(lambda x: PATHOLOGY_MAP.get(str(x).strip(), np.nan)),
        "site_raw": cancer["病变1受累部位"],
        "size_category": cancer["病变1大小"].map(parse_size),
        "elevated": cancer["病变1形态"].map(lambda x: {"隆起": "Yes", "平坦": "No"}.get(str(x).strip(), np.nan)),
        "ulcerative": cancer["病变1是否有溃疡"].map(lambda x: {"有": "Yes", "无": "No"}.get(str(x).strip(), np.nan)),
    }
)
lesions["size_category"] = lesions["size_category"].replace({">20 mm": "＞20 mm"})

exp = lesions[lesions["group"] == "A"]
ctrl = lesions[lesions["group"] == "M"]

print("A组病灶数", len(exp))
print("M组病灶数", len(ctrl))
print()

for label in ["Adenocarcinoma", "HGIN", "LGIN", "Adenoma", "NET", "Lymphoma", "GIST", "Metastasis"]:
    print(label)
    print("A组", int((exp["pathology"] == label).sum()))
    print("M组", int((ctrl["pathology"] == label).sum()))
print("病理亚型检验P值", sparse_multicat(pd.crosstab(lesions["pathology"], lesions["group"]).reindex(
    index=["Adenocarcinoma", "HGIN", "LGIN", "Adenoma", "NET", "Lymphoma", "GIST", "Metastasis"],
    columns=["A", "M"], fill_value=0).to_numpy()))
print()

# 按部位
print("按部位")
for label, kw in [
    ("Cardia-Fundus", "贲门|胃底"),
    ("Body", "胃体"),
    ("Angulus", "胃角"),
    ("Antrum-Pylorus", "胃窦|幽门"),
    ("Gastric remnant or thoracic stomach", "残胃|胸胃"),
]:
    exp_count = int(exp["site_raw"].astype(str).str.contains(kw, regex=True, na=False).sum())
    ctrl_count = int(ctrl["site_raw"].astype(str).str.contains(kw, regex=True, na=False).sum())
    p = chi_or_fisher(np.array([[exp_count, len(exp) - exp_count], [ctrl_count, len(ctrl) - ctrl_count]]))
    print(label)
    print("A组", exp_count)
    print("M组", ctrl_count)
    print("P值", p)
print()

# 按大小
print("按大小")
for label in ["≤5 mm", "5-10 mm", "10-20 mm", "＞20 mm", "NA"]:
    print(label)
    print("A组", int((exp["size_category"] == label).sum()))
    print("M组", int((ctrl["size_category"] == label).sum()))
print("大小检验P值", sparse_multicat(pd.crosstab(lesions["size_category"], lesions["group"]).reindex(
    index=["≤5 mm", "5-10 mm", "10-20 mm", "＞20 mm", "NA"],
    columns=["A", "M"], fill_value=0).to_numpy()))
print()

# 隆起
print("隆起")
for label in ["Yes", "No"]:
    print(label)
    print("A组", int((exp["elevated"] == label).sum()))
    print("M组", int((ctrl["elevated"] == label).sum()))
print("隆起检验P值", chi_or_fisher(pd.crosstab(lesions["elevated"], lesions["group"]).reindex(
    index=["Yes", "No"], columns=["A", "M"], fill_value=0).to_numpy()))
print()

# 溃疡
print("溃疡")
for label in ["Yes", "No"]:
    print(label)
    print("A组", int((exp["ulcerative"] == label).sum()))
    print("M组", int((ctrl["ulcerative"] == label).sum()))
print("溃疡检验P值", chi_or_fisher(pd.crosstab(lesions["ulcerative"], lesions["group"]).reindex(
    index=["Yes", "No"], columns=["A", "M"], fill_value=0).to_numpy()))