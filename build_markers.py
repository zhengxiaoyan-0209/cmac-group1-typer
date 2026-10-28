# -*- coding: utf-8 -*-
"""
build_markers.py
================
从原始分析产物构建 "C. macginleyi Group 1 快速鉴定" 网页所需的数据文件。

输入
----
1) 参考基因组 fasta          : CP069516.1 (= GCA_016889465.1 = FDAARGOS_1195)
2) Group1 特异 SNP 数据库    : mcg_snp_database_FDAARGOS_1195.json
3) SNP 基因注释表            : Group1特异SNP基因定位表.xlsx
4) 菌株分组表                : group.csv

注：组特异/富集基因序列（group_634 等）已不在网页里使用，故不再打进 markers。

输出
----
assets/markers.js  —— 浏览器可直接 <script src> 加载（避免 file:// 下 fetch 被 CORS 拦截）
data/markers.json  —— 同一份数据的 JSON 版本，便于脚本复用 / 后续更新

设计要点（为什么不存"整条探针"，而存"侧翼锚"）
------------------------------------
若存 51-mer 完整探针并要求完全匹配，查询基因组在该窗口内只要多出任意一个
其他变异就会漏判。改为：只锚定 SNP 紧邻的一侧 25 bp，用子串定位成功后读取
其相邻的一个碱基。这样侧翼（25 bp 以外、甚至缺口以外）的其他变异不影响判定，
召回率显著更高；同时 25 bp 在 2.5 Mb 基因组中期望出现 << 1 次，特异性足够。
"""

import json
import os
import re
import zipfile
from xml.etree import ElementTree as ET

# ---------------------------------------------------------------- 路径配置
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DESK = r"C:/Users/15963/Desktop/麦氏棒杆菌"

REF_FNA = os.path.join(BASE, "data", "reference_CP069516.1.fna")
SNP_JSON = os.path.join(DESK, "组内SNP", "mcg_snp_database_FDAARGOS_1195.json")
SNP_XLSX = os.path.join(DESK, "组内SNP", "Group1特异SNP基因定位表.xlsx")
GROUP_CSV = os.path.join(DESK, "组内SNP", "group.csv")

OUT_JS = os.path.join(BASE, "assets", "markers.js")
OUT_JSON = os.path.join(BASE, "data", "markers.json")

# 锚定长度梯度：先用长锚保证唯一性，找不到再降级以提高召回
ANCHOR_LEN = 25
FLANK_LEN = 25          # 实际截取的侧翼长度（>= ANCHOR_LEN）
ANCHOR_FALLBACK = [25, 20, 16]

# xlsx 里 11 个位点是基因间区，翻成英文
FEATURE_EN = {"基因间区": "intergenic", "CDS": "CDS"}

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


# ---------------------------------------------------------------- 基础工具
def read_single_fasta(path):
    """读取单条序列，返回 (header, seq)。序列统一大写、非法字符转 N。"""
    head, chunks = None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                head = line[1:]
            else:
                chunks.append(line)
    seq = "".join(chunks).upper()
    seq = re.sub(r"[^ACGTN]", "N", seq)
    return head, seq


def read_multi_fasta(path):
    """读取多序列 fasta，返回 [(header, seq), ...]"""
    out, head, chunks = [], None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if head is not None:
                    out.append((head, "".join(chunks)))
                head, chunks = line[1:], []
            else:
                chunks.append(line)
    if head is not None:
        out.append((head, "".join(chunks)))
    res = []
    for h, s in out:
        s = re.sub(r"[^ACGTN]", "N", s.upper())
        res.append((h, s))
    return res


def read_xlsx_sheet(zf, idx):
    """最小依赖地读取 xlsx 的某个 sheet（inlineStr 格式）。"""
    root = ET.fromstring(zf.read("xl/worksheets/sheet%d.xml" % idx))
    rows = []
    for row in root.iter(NS + "row"):
        d = {}
        for c in row:
            key = "".join(ch for ch in c.get("r") if ch.isalpha())
            t = c.find(NS + "is/" + NS + "t")
            v = c.find(NS + "v")
            d[key] = t.text if t is not None else (v.text if v is not None else "")
        rows.append(d)
    return rows


def complement(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))


def revcomp(s):
    return complement(s)[::-1]


# ---------------------------------------------------------------- 1. 参考基因组
ref_head, ref = read_single_fasta(REF_FNA)
print("参考基因组: %s" % ref_head)
print("长度      : %d bp" % len(ref))

# ---------------------------------------------------------------- 2. SNP 位点
snp_db = json.load(open(SNP_JSON, encoding="utf-8"))
snps = snp_db["group_specific_snps"]["1"]
groups = snp_db["groups"]
print("Group1 特异 SNP 数: %d" % len(snps))

# ---------------------------------------------------------------- 3. 注释表
zf = zipfile.ZipFile(SNP_XLSX)
anno = {}
for r in read_xlsx_sheet(zf, 1)[1:]:          # sheet1 = SNP明细
    if not r.get("B"):
        continue
    pos = int(float(r["B"]))
    anno[pos] = {
        "feature": FEATURE_EN.get(r.get("E", ""), r.get("E", "")),
        "gene": r.get("F", "") or "-",
        "locus_tag": r.get("G", "") or "-",
        "product": r.get("H", "") or "-",
        "part": r.get("I", ""),
        "region_kb": r.get("J", ""),
    }

# xlsx 里的功能归类是中文，翻成英文供英文界面使用
CATEGORY_EN = {
    "脂质/水解酶": "Lipid / hydrolase",
    "未注释": "Uncharacterised",
    "转运系统": "Transport",
    "代谢": "Metabolism",
    "应激/能量调控": "Stress / energy regulation",
    "其他": "Other",
}

gene_summary = []
for r in read_xlsx_sheet(zf, 2)[1:]:          # sheet2 = 基因命中汇总
    if not r.get("A"):
        continue
    cat = r["E"]
    gene_summary.append({
        "locus": r["A"], "product": r["B"], "n_snp": int(float(r["C"])),
        "span": r["D"], "category": CATEGORY_EN.get(cat, cat),
    })

# ---------------------------------------------------------------- 4. 分组信息
groups_map = {}
with open(GROUP_CSV, encoding="utf-8-sig") as fh:
    next(fh)
    for line in fh:
        line = line.strip()
        if not line:
            continue
        s, g = line.split(",")
        groups_map[s] = g
print("分组标本数: Group1=%d Group2=%d Group3=%d" % (
    sum(1 for v in groups_map.values() if v == "1"),
    sum(1 for v in groups_map.values() if v == "2"),
    sum(1 for v in groups_map.values() if v == "3"),
))

# ---------------------------------------------------------------- 5. 构建位点探针
markers = []
bad = 0
for key, alt in sorted(snps.items(), key=lambda kv: int(kv[0].split("|")[1])):
    acc, spos, ref_base = key.split("|")
    pos = int(spos)                       # 1-based
    if acc != "CP069516.1":
        raise SystemExit("未知参考序列编号: %s" % acc)
    obs = ref[pos - 1]
    if obs != ref_base:
        bad += 1
        print("  !! 坐标校验失败 %s 参考碱基=%s 记录=%s" % (key, obs, ref_base))

    left = ref[max(0, pos - 1 - FLANK_LEN): pos - 1]      # SNP 紧邻上游
    right = ref[pos: pos + FLANK_LEN]                     # SNP 紧邻下游
    a = anno.get(pos, {})
    markers.append({
        "pos": pos,
        "ref": ref_base,
        "alt": alt,
        "gene": a.get("gene", "-"),
        "locus_tag": a.get("locus_tag", "-"),
        "product": a.get("product", "-"),
        "feature": a.get("feature", ""),
        "region_kb": a.get("region_kb", ""),
        "L": left,
        "R": right,
    })

if bad:
    raise SystemExit("有 %d 个位点与参考基因组不符，请检查坐标系统" % bad)
print("坐标校验: %d/%d 通过（1-based）" % (len(markers) - bad, len(markers)))

# ---- 单倍型感知侧翼 ----
# 85 个位点里存在密集簇（最密处 25 bp 窗口内有 6 个位点，如 2390 kb 区）。
# 若侧翼锚只按参考序列取，锚本身就覆盖了邻近的突变位点，真实 Group 1 菌株
# 在这些邻位是 alt 等位基因 -> 锚失配 -> 该位点漏判。
# 因此每个侧翼额外生成一份 "Group 1 单倍型版本"：把落在窗口内的其它位点
# 一并替换成 alt。判定引擎两套都用，避免只偏向某一种单倍型。
altmap = {m["pos"]: m["alt"] for m in markers}
for m in markers:
    pos = m["pos"]
    L = m["L"]
    R = m["R"]
    # L 覆盖的参考坐标区间(1-based): pos-FLANK_LEN .. pos-1
    start = max(1, pos - FLANK_LEN)
    la = "".join(altmap.get(pp, ref[pp - 1]) for pp in range(start, pos))
    # R 覆盖: pos+1 .. pos+FLANK_LEN
    end = min(len(ref), pos + FLANK_LEN)
    ra = "".join(altmap.get(pp, ref[pp - 1]) for pp in range(pos + 1, end + 1))
    m["La"] = la if la != L else None          # None 表示该窗口内无邻近位点
    m["Ra"] = ra if ra != R else None

# 基因组范围内的寡核苷酸自检：确认侧翼锚在全基因组足够唯一
uniq_ok = 0
for m in markers:
    nL = ref.count(m["L"]) + ref.count(revcomp(m["L"]))
    nR = ref.count(m["R"]) + ref.count(revcomp(m["R"]))
    m["uniq_L"] = nL
    m["uniq_R"] = nR
    if nL == 1 or nR == 1:
        uniq_ok += 1
print("侧翼锚在参考基因组内唯一: %d/%d" % (uniq_ok, len(markers)))

# ---------------------------------------------------------------- 7. 输出
payload = {
    "meta": {
        "reference": "CP069516.1 (GCA_016889465.1, C. macginleyi FDAARGOS_1195)",
        "ref_header": ref_head,
        "ref_length": len(ref),
        "n_markers": len(markers),
        "anchor_fallback": ANCHOR_FALLBACK,
        "groups": groups,
        "n_group_members": {
            "1": sum(1 for v in groups_map.values() if v == "1"),
            "2": sum(1 for v in groups_map.values() if v == "2"),
            "3": sum(1 for v in groups_map.values() if v == "3"),
        },
        "note": "The 85 loci are strictly Group 1 specific SNPs (present in 100% of "
                "Group 1 genomes, absent from all others); they are not 85 genes.",
    },
    "snps": markers,
    "gene_summary": gene_summary,
}

os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
with open(OUT_JSON, "w", encoding="utf-8") as fh:
    json.dump(payload, fh, ensure_ascii=False, indent=1)

with open(OUT_JS, "w", encoding="utf-8") as fh:
    fh.write("/* Generated file - do not edit by hand. Rebuild: scripts/build_markers.py */\n")
    fh.write("window.CMAC_MARKERS = ")
    json.dump(payload, fh, ensure_ascii=False, separators=(",", ":"))
    fh.write(";\n")

print("\n输出:")
print("  %s  (%.1f KB)" % (OUT_JS, os.path.getsize(OUT_JS) / 1024))
print("  %s  (%.1f KB)" % (OUT_JSON, os.path.getsize(OUT_JSON) / 1024))
