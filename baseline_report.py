import pandas as pd
from scipy import stats

df = pd.read_csv(r"")

# 输入实际分组标签
a = df[df["cade"] == "A"]
b = df[df["cade"] == "M"]

# 样本量
print("A组", len(a))
print("B组", len(b))

# 性别
for s in df["sex_group"].unique():
    print("A组性别", s, (a["sex_group"] == s).sum())
    print("B组性别", s, (b["sex_group"] == s).sum())
print("性别检验", stats.chi2_contingency(pd.crosstab(df["sex_group"], df["cade"]))[1])

# 年龄
print("A组年龄均值", a["age_years"].mean())
print("B组年龄均值", b["age_years"].mean())
print("A组年龄标准差", a["age_years"].std()) # 默认 ddot = 1
print("B组年龄标准差", b["age_years"].std())
print("年龄t检验", stats.ttest_ind(a["age_years"], b["age_years"], equal_var=False)[1]) # Welch检验

# 适应症
for v in df["indication"].unique():
    print("A组适应症", v, (a["indication"] == v).sum())
    print("B组适应症", v, (b["indication"] == v).sum())
print("适应症检验", stats.chi2_contingency(pd.crosstab(df["indication"], df["cade"]))[1])

# 患者来源
print("A组消化内科", (a["patient_source_group"] == "Gastroenterology").sum())
print("B组消化内科", (b["patient_source_group"] == "Gastroenterology").sum())
print("A组非消化", (a["patient_source_group"] == "Non-gastroenterology").sum())
print("B组非消化", (b["patient_source_group"] == "Non-gastroenterology").sum())
print("来源检验", stats.chi2_contingency(pd.crosstab(df["patient_source_group"], df["cade"]))[1])

# HP状态
for raw in df["hp_detail"].unique():
    print("A组HP", raw, (a["hp_detail"] == raw).sum())
    print("B组HP", raw, (b["hp_detail"] == raw).sum())
print("HP检验", stats.chi2_contingency(pd.crosstab(df["hp_detail"], df["cade"]))[1])

# 内镜量
for vol in df["endoscopist_volume_group"].unique():
    print("A组内镜量", vol, (a["endoscopist_volume_group"] == vol).sum())
    print("B组内镜量", vol, (b["endoscopist_volume_group"] == vol).sum())
print("内镜量检验", stats.chi2_contingency(pd.crosstab(df["endoscopist_volume_group"], df["cade"]))[1])

# 平台
for code in df["platform_code"].unique():
    print("A组平台", code, (a["platform_code"] == code).sum())
    print("B组平台", code, (b["platform_code"] == code).sum())
print("平台检验", stats.chi2_contingency(pd.crosstab(df["platform_code"], df["cade"]))[1])

# 研究阶段
for phase in df["study_phase"].unique():
    print("A组", phase, (a["study_phase"] == phase).sum())
    print("B组", phase, (b["study_phase"] == phase).sum())
print("阶段检验", stats.chi2_contingency(pd.crosstab(df["study_phase"], df["cade"]))[1])