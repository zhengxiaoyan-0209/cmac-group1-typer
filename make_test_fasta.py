# -*- coding: utf-8 -*-
"""生成用于自测的小型 FASTA（每个 SNP 位点上下游各 500 bp 的窗口）。

- non_group1_windows.fa : 完全取参考序列 —— 应判为「非 Group 1」，85/85
- group1_windows.fa     : 同样窗口，但把所有组特异位点换成 Group 1 等位基因 —— 应判为「Group 1」，85/85
"""
import json
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WIN = 500

ref = "".join(
    l.strip() for l in open(os.path.join(BASE, "data", "reference_CP069516.1.fna"), encoding="utf-8")
    if not l.startswith(">")
).upper()

db = json.load(open(os.path.join(BASE, "data", "markers.json"), encoding="utf-8"))
snps = sorted(db["snps"], key=lambda m: m["pos"])
altmap = {m["pos"]: m["alt"] for m in snps}

# 相邻窗口合并，避免重复片段
merged = []
for m in snps:
    s, e = max(1, m["pos"] - WIN), min(len(ref), m["pos"] + WIN)
    if merged and s <= merged[-1][1] + 1:
        merged[-1][1] = max(merged[-1][1], e)
    else:
        merged.append([s, e])

out_dir = os.path.join(BASE, "test")
os.makedirs(out_dir, exist_ok=True)


def wrap(seq, width=70):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


for tag, use_alt in (("non_group1", False), ("group1", True)):
    path = os.path.join(out_dir, "%s_windows.fa" % tag)
    with open(path, "w", encoding="utf-8") as fh:
        for i, (s, e) in enumerate(merged, 1):
            piece = list(ref[s - 1:e])
            for p in range(s, e + 1):
                if p in altmap and use_alt:
                    piece[p - s] = altmap[p]
            fh.write(">window_%02d CP069516.1:%d-%d (%s)\n" % (i, s, e, "Group1 单倍型" if use_alt else "参考单倍型"))
            fh.write(wrap("".join(piece)) + "\n")
    print("%-24s %d 条片段  %.0f KB" % (os.path.basename(path), len(merged), os.path.getsize(path) / 1024))

print("\n期望结果：")
print("  non_group1_windows.fa -> 非 Group 1（匹配率 0.0%%，85/85 可判定）")
print("  group1_windows.fa     -> Group 1（匹配率 100.0%%，85/85 可判定）")
