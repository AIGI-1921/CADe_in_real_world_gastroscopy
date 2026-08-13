import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

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
    cols = [outcome, "g", "cluster"] + COVS
    work = frame[cols].copy().dropna()
    x = design_matrix(work)
    res = sm.GLM(work[outcome].astype(float), x, family=sm.families.Poisson()).fit(
        cov_type="cluster",
        cov_kwds={"groups": work["cluster"], "use_correction": True, "df_correction": True},
    )
    return float(res.pvalues["g"])


def ols_cluster(frame, outcome):
    cols = [outcome, "g", "cluster"] + COVS
    work = frame[cols].copy().dropna()
    x = design_matrix(work)
    res = sm.OLS(work[outcome].astype(float), x).fit(
        cov_type="cluster",
        cov_kwds={"groups": work["cluster"], "use_correction": True, "df_correction": True},
    )
    return float(res.pvalues["g"])

file_name = r""

df = pd.read_csv(file_name)
df["cluster"] = df["platform_session_cluster"]
a = df[df["cade"] == "A"]
b = df[df["cade"] == "M"]
df["g"] = (df["cade"] == "A").astype(float)
print("A组人数", len(a))
print("B组人数", len(b))
print()

labels = []
raw_pvals = []
cluster_pvals = []

# 胃肿瘤检出率
ac = int(a["outcome_neoplasm"].sum())
bc = int(b["outcome_neoplasm"].sum())
rawp = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
clup = poisson_cluster(df, "outcome_neoplasm")
labels.append("胃肿瘤检出率")
raw_pvals.append(rawp)
cluster_pvals.append(clup)

# 肠化或萎缩性胃炎检出率
ac = int(a["outcome_atrophy_im"].sum())
bc = int(b["outcome_atrophy_im"].sum())
rawp = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
clup = poisson_cluster(df, "outcome_atrophy_im")
labels.append("肠化或萎缩性胃炎检出率")
raw_pvals.append(rawp)
cluster_pvals.append(clup)

# 胃息肉检出率
ac = int(a["outcome_polyp"].sum())
bc = int(b["outcome_polyp"].sum())
rawp = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
clup = poisson_cluster(df, "outcome_polyp")
labels.append("胃息肉检出率")
raw_pvals.append(rawp)
cluster_pvals.append(clup)

# 胃炎检出率
ac = int(a["outcome_gastritis"].sum())
bc = int(b["outcome_gastritis"].sum())
rawp = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
clup = poisson_cluster(df, "outcome_gastritis")
labels.append("胃炎检出率")
raw_pvals.append(rawp)
cluster_pvals.append(clup)

# 周转时间
ta = a[a["turnover_valid"] == 1]["turnover_time"]
tb = b[b["turnover_valid"] == 1]["turnover_time"]
rawp = stats.ttest_ind(ta, tb, equal_var=False)[1]
use = df[df["turnover_valid"] == 1].copy()
clup = ols_cluster(use, "turnover_time")
labels.append("周转时间")
raw_pvals.append(rawp)
cluster_pvals.append(clup)

# 活检率
ac = int(a["biopsy_any"].sum())
bc = int(b["biopsy_any"].sum())
rawp = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
clup = poisson_cluster(df, "biopsy_any")
labels.append("活检率")
raw_pvals.append(rawp)
cluster_pvals.append(clup)

# 贲门-胃底活检率
ac = int(a["biopsy_rate_cardia_fundus"].sum())
bc = int(b["biopsy_rate_cardia_fundus"].sum())
rawp = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
clup = poisson_cluster(df, "biopsy_rate_cardia_fundus")
labels.append("贲门-胃底活检率")
raw_pvals.append(rawp)
cluster_pvals.append(clup)

# 胃体活检率
ac = int(a["biopsy_rate_body"].sum())
bc = int(b["biopsy_rate_body"].sum())
rawp = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
clup = poisson_cluster(df, "biopsy_rate_body")
labels.append("胃体活检率")
raw_pvals.append(rawp)
cluster_pvals.append(clup)

# 胃角活检率
ac = int(a["biopsy_rate_angulus"].sum())
bc = int(b["biopsy_rate_angulus"].sum())
rawp = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
clup = poisson_cluster(df, "biopsy_rate_angulus")
labels.append("胃角活检率")
raw_pvals.append(rawp)
cluster_pvals.append(clup)

# 胃窦-幽门活检率
ac = int(a["biopsy_rate_antrum_pylorus"].sum())
bc = int(b["biopsy_rate_antrum_pylorus"].sum())
rawp = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
clup = poisson_cluster(df, "biopsy_rate_antrum_pylorus")
labels.append("胃窦-幽门活检率")
raw_pvals.append(rawp)
cluster_pvals.append(clup)

# 残胃活检率
ac = int(a["biopsy_rate_remnant"].sum())
bc = int(b["biopsy_rate_remnant"].sum())
rawp = stats.chi2_contingency([[ac, len(a) - ac], [bc, len(b) - bc]], correction=False)[1]
clup = poisson_cluster(df, "biopsy_rate_remnant")
labels.append("残胃活检率")
raw_pvals.append(rawp)
cluster_pvals.append(clup)

# 平均活检次数
av = a["biopsy_count_any"].fillna(0)
bv = b["biopsy_count_any"].fillna(0)
rawp = stats.mannwhitneyu(av, bv)[1]
clup = ols_cluster(df, "biopsy_count_any")
labels.append("平均活检次数")
raw_pvals.append(rawp)
cluster_pvals.append(clup)

# 贲门-胃底平均活检次数
av = a["biopsy_count_cardia_fundus"].fillna(0)
bv = b["biopsy_count_cardia_fundus"].fillna(0)
rawp = stats.mannwhitneyu(av, bv)[1]
clup = ols_cluster(df, "biopsy_count_cardia_fundus")
labels.append("贲门-胃底平均活检次数")
raw_pvals.append(rawp)
cluster_pvals.append(clup)

# 胃体平均活检次数
av = a["biopsy_count_body"].fillna(0)
bv = b["biopsy_count_body"].fillna(0)
rawp = stats.mannwhitneyu(av, bv)[1]
clup = ols_cluster(df, "biopsy_count_body")
labels.append("胃体平均活检次数")
raw_pvals.append(rawp)
cluster_pvals.append(clup)

# 胃角平均活检次数
av = a["biopsy_count_angulus"].fillna(0)
bv = b["biopsy_count_angulus"].fillna(0)
rawp = stats.mannwhitneyu(av, bv)[1]
clup = ols_cluster(df, "biopsy_count_angulus")
labels.append("胃角平均活检次数")
raw_pvals.append(rawp)
cluster_pvals.append(clup)

# 胃窦-幽门平均活检次数
av = a["biopsy_count_antrum_pylorus"].fillna(0)
bv = b["biopsy_count_antrum_pylorus"].fillna(0)
rawp = stats.mannwhitneyu(av, bv)[1]
clup = ols_cluster(df, "biopsy_count_antrum_pylorus")
labels.append("胃窦-幽门平均活检次数")
raw_pvals.append(rawp)
cluster_pvals.append(clup)

# 残胃平均活检次数
av = a["biopsy_count_remnant"].fillna(0)
bv = b["biopsy_count_remnant"].fillna(0)
rawp = stats.mannwhitneyu(av, bv)[1]
clup = ols_cluster(df, "biopsy_count_remnant")
labels.append("残胃平均活检次数")
raw_pvals.append(rawp)
cluster_pvals.append(clup)

raw_reject, raw_q = fdrcorrection(raw_pvals)
cluster_reject, cluster_q = fdrcorrection(cluster_pvals)

print("各终点原始P值及FDR校正后Q值")
for i in range(17):
    print(labels[i])
    print("原始P值", raw_pvals[i])
    print("原始校正Q值", raw_q[i])
    print("聚类P值", cluster_pvals[i])
    print("聚类校正Q值", cluster_q[i])
    print()