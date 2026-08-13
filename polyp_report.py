import numpy as np
import pandas as pd
from scipy import stats

# 路径
cohort_path = r""
polyp_path = r""

# 息肉亚型/部位映射
SUB_MAP = {
    1: "Fundic gland polyp",
    2: "Hyperplastic polyp",
    3: "Hamartomatous polyp",
    4: "Inflammatory polyp",
    5: "Hamartomatous polyp",
}
SITE_MAP = {
    1: "Cardia-Fundus",
    2: "Body",
    3: "Angulus",
    4: "Antrum-Pylorus",
    5: "Gastric remnant or thoracic stomach",
}

SUB_ORDER = ["Fundic gland polyp", "Hyperplastic polyp", "Hamartomatous polyp", "Inflammatory polyp"]
SITE_ORDER = ["Cardia-Fundus", "Body", "Angulus", "Antrum-Pylorus", "Gastric remnant or thoracic stomach", "NA"]
SIZE_ORDER = ["≤5 mm", "5-10 mm", "10-20 mm", "＞20 mm", "NA"]
THREE_ORDER = ["Yes", "No", "NA"]


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

def sparse_multicategory_p(table):
    keep_rows = table.sum(axis=1) > 0
    keep_cols = table.sum(axis=0) > 0
    table = table[np.ix_(keep_rows, keep_cols)]
    row_sums = table.sum(axis=1)
    col_sums = table.sum(axis=0)
    try:
        _, p, _, expected = stats.chi2_contingency(table, correction=False)
    except ValueError:
        p = None
        expected = np.full(table.shape, np.inf)
    if p is None or float(expected.min()) < 5:
        # 蒙特卡洛 Fisher-Freeman-Halton
        expected = np.outer(row_sums, col_sums) / table.sum()
        stat = float((((table - expected) ** 2) / expected).sum())
        rng = np.random.default_rng(20260725)
        sims = stats.random_table(row_sums, col_sums).rvs(size=200000, random_state=rng)
        sim_stats = ((sims - expected) ** 2 / expected).sum(axis=(1, 2))
        return float((np.sum(sim_stats >= stat) + 1) / (200000 + 1))
    return float(p)

source = pd.read_excel(polyp_path, sheet_name="Sheet1")
cohort = pd.read_csv(cohort_path)[["Barcode", "cade"]].drop_duplicates()

records = []
for _, row in source.iterrows():
    barcode = str(row["Barcode"]).strip()
    for slot in (1, 2, 3):
        if slot == 1:
            path_col = "息肉病变1病理（1-胃底腺息肉；2-增生性息肉；3-错构瘤性息肉；4-炎症性息肉；5-P-J息肉）"
            site_col = "息肉病变位置1 (贲门-胃底-1/胃体-2/胃角-3/胃窦-幽门-4/残胃-5)"
            size_col = "息肉病变大小1 (mm)"
            morph_col = "息肉病变形态1 (非扁平-1/扁平-2/多发非扁平-3/多发扁平-4)"
        else:
            path_col = f"息肉病变{slot}病理"
            site_col = f"息肉病变位置{slot}"
            size_col = f"息肉病变大小{slot}"
            morph_col = f"息肉病变形态{slot}"
        path_value = row.get(path_col)
        site_value = row.get(site_col)
        size_value = row.get(size_col)
        morph_value = row.get(morph_col)
        if all(pd.isna(v) or str(v).strip() == "" for v in [path_value, site_value, size_value, morph_value]):
            continue
        subtype_code = pd.to_numeric(path_value, errors="coerce")
        site_code = pd.to_numeric(site_value, errors="coerce")
        morph_code = pd.to_numeric(morph_value, errors="coerce")
        size_num = pd.to_numeric(size_value, errors="coerce")
        size_cat = "NA"
        if pd.notna(size_num):
            if size_num <= 5:
                size_cat = "≤5 mm"
            elif size_num <= 10:
                size_cat = "5-10 mm"
            elif size_num <= 20:
                size_cat = "10-20 mm"
            else:
                size_cat = ">20 mm"
        records.append(
            {
                "Barcode": barcode,
                "subtype": SUB_MAP.get(int(subtype_code), f"Other-{int(subtype_code)}") if pd.notna(subtype_code) else "NA",
                "site": SITE_MAP.get(int(site_code), "NA") if pd.notna(site_code) else "NA",
                "size_mm": size_num,
                "size_cat": size_cat,
                "elevated": "Yes" if morph_code in {1, 3} else "No" if morph_code in {2, 4} else "NA",
                "multifocal": "Yes" if morph_code in {3, 4} else "No" if morph_code in {1, 2} else "NA",
            }
        )

polyp_records = pd.DataFrame(records).merge(cohort, on="Barcode", how="left")
polyp_records["size_cat"] = polyp_records["size_cat"].replace({">20 mm": "＞20 mm"})
polyp_records = polyp_records[polyp_records["cade"].isin(["A", "M"])].copy()

analysis = polyp_records[polyp_records["subtype"].isin(SUB_ORDER)]
a = analysis[analysis["cade"] == "A"]
b = analysis[analysis["cade"] == "M"]
print("A组息肉记录数", len(a))
print("M组息肉记录数", len(b))
print()

def block_p(column, order):
    table = pd.crosstab(analysis[column], analysis["cade"]).reindex(index=order, columns=["A", "M"], fill_value=0)
    return sparse_multicategory_p(table.to_numpy())

for title, column, order in [
    ("病理亚型", "subtype", SUB_ORDER),
    ("部位", "site", SITE_ORDER),
    ("大小", "size_cat", SIZE_ORDER),
    ("隆起", "elevated", THREE_ORDER),
    ("多发", "multifocal", THREE_ORDER),
]:
    p = block_p(column, order)
    print(title)
    print("整体P值", p)
    for label in order:
        print(label)
        print("A组", int((a[column] == label).sum()))
        print("M组", int((b[column] == label).sum()))
    print()

hyper_barcodes = set(polyp_records.loc[polyp_records["subtype"] == "Hyperplastic polyp", "Barcode"].astype(str))
large_barcodes = set(polyp_records.loc[pd.to_numeric(polyp_records["size_mm"], errors="coerce") > 10, "Barcode"].astype(str))

cohort_all = pd.read_csv(cohort_path)[["Barcode", "cade"]].drop_duplicates()
cohort_all["hyper"] = cohort_all["Barcode"].astype(str).isin(hyper_barcodes).astype(int)
cohort_all["large"] = cohort_all["Barcode"].astype(str).isin(large_barcodes).astype(int)

ca = cohort_all[cohort_all["cade"] == "A"]
cb = cohort_all[cohort_all["cade"] == "M"]
for col, label in [
    ("hyper", "增生性息肉检出率"),
    ("large", ">10mm息肉样病变检出率"),
]:
    ac = int(ca[col].sum())
    bc = int(cb[col].sum())
    rr, lo, hi = calc_rr(ac, len(ca), bc, len(cb))
    p = stats.chi2_contingency([[ac, len(ca) - ac], [bc, len(cb) - bc]], correction=False)[1]
    print(label)
    print("A组事件数", ac)
    print("A组人数", len(ca))
    print("M组事件数", bc)
    print("M组人数", len(cb))
    print("RR", rr)
    print("RR置信区间", lo, hi)
    print("P值", p)
    print()