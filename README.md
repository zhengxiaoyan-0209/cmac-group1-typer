# C. macginleyi Group 1 Typer

眼部分离 **麦氏棒状杆菌（Corynebacterium macginleyi）** 高危谱系 Group 1 的在线鉴定工具。
拖入组装好的基因组，浏览器本地几秒内给出判定，**序列不上传、不需要服务器**。
界面为英文（面向论文读者），此 README 面向工具的维护者（你）。

> 论文：Zheng X, Xie J, Zheng H, et al. *Population structure and antimicrobial resistance of
> ocular Corynebacterium macginleyi isolates.*（IOVS 投稿中）
> 65 株基因组分为 3 个地理结构明显的谱系；其中 Group 1（8 株，均来自瑞典）
> 与更严重的临床表现相关。

---

## 一个容易混淆的概念

原始分析里有两组数字，容易混为一谈：

| 维度 | 数量 | 位置 / 名称 | 论文中的描述 |
|---|---|---|---|
| **Group 1 组特异 SNP** | **85 个位点** | 落在 18 个基因 + 若干基因间区 | "85 group-specific SNPs" |
| Group 1 组特异基因 | 1 个 | `group_634`（hypothetical protein, 600 bp） | "the unique gene group_634" |
| Group 1 **富集**基因（非绝对特异） | 5 个 | `cas9` `cas1` `cas2` `adh` `group_780`（7/8，87.5%） | "five genes were carried by 87.5%" |

**「85 个」是 SNP 位点，不是 85 个基因。** 本工具只使用这 85 个 SNP 位点作为判定依据；
`group_634` 与富集基因不在判定逻辑里（早期版本曾作为辅助旁证展示，现已移除）。

---

## 在线 vs 本地

直接用 GitHub Pages 托管即可；也可以**双击 `index.html` 本地离线使用**
（数据以 `<script>` 方式内联在 `assets/markers.js`，不走 `fetch`，所以 `file://` 协议下不会被 CORS 拦截）。

---

## 目录结构

```
cmac-group1-typer/
├── index.html                 单页应用
├── assets/
│   ├── markers.js             85 个位点标记数据（探针）——自动生成，勿手改
│   ├── core.js                鉴定引擎（浏览器 / Node 共用）
│   ├── app.js                 界面逻辑
│   └── style.css
├── data/
│   ├── markers.json           同一份数据的 JSON 版本
│   └── reference_CP069516.1.fna  参考基因组（2.4 MB）
├── scripts/
│   ├── build_markers.py       从原始分析产物重建 markers.js
│   ├── make_test_fasta.py     生成自测用 FASTA
│   ├── typer.js               命令行批处理版（与网页共用同一引擎）
│   └── validate.js            离线自洽验证（Node）
├── test/
│   ├── non_group1_windows.fa  自测：应判为非 Group 1
│   └── group1_windows.fa      自测：应判为 Group 1
└── README.md
```

> `data/reference_CP069516.1.fna` 从 NCBI 下载，仅由构建脚本读取，
> **页面上不会用到**，如果不介意无法重新生成探针，可以从仓库里删掉以节省体积。

---

## 使用方法

1. 打开页面
2. 拖入（或点击选择）**组装好的基因组 FASTA**，支持 `.fa / .fna / .fasta / .txt`，也支持 `.gz`；
   **可以一次拖入整个文件夹**（比如 65 株的 assembly 目录），每个文件算一个样本
3. 点【Run typing】
4. 多个样本时先出**结果汇总表**，点任意一行可查看该样本的 85 个位点明细

参数（可调，改完自动重算）：

| 参数 | 默认 | 含义 |
|---|---|---|
| 最小位点覆盖 | 70% | 可判定位点不足这个比例时不出结论（证据不足） |
| 判为 Group 1 阈值 | 90% | Group 1 匹配率 ≥ 此值 → Group 1 |
| 判为非 Group 1 阈值 | 10% | Group 1 匹配率 ≤ 此值 → 非 Group 1 |

**Group 1 匹配率** = Group 1 型位点数 ÷ 可判定位点数。

可导出两种 CSV：**汇总表**（每样本一行，适合直接进 R/Excel 做图）和**位点明细**（当前选中样本的 85 行）。

### 自测

仓库自带两个 25 KB 的小文件，拖进去就应该得到：

| 文件 | 期望结论 | 期望匹配率 |
|---|---|---|
| `test/non_group1_windows.fa` | 非 Group 1 | 0.0%（85/85 可判定） |
| `test/group1_windows.fa` | Group 1 | 100.0%（85/85 可判定） |

---

## 判定原理

论文流程是：把待检基因组比对到参考基因组 `CP069516.1`
（= `GCA_016889465.1` = 菌株 FDAARGOS_1195，2,493,537 bp），
读取 85 个 Group 1 特异位点的等位基因，按"匹配比例"判定归属。

本工具用**「锚定法」**在浏览器里复现同一逻辑，无需全基因组比对：

1. 取 SNP 紧邻的一侧序列（默认 25 bp）作为锚
2. 在待检基因组**及其反向互补序列**中定位该锚
3. 锚定成功后，读取与锚相邻的那 1 个碱基 → 该位点的实测等位基因
4. 25 bp 锚找不到就降级到 20 / 16 bp
5. 锚多次命中且读出不同碱基 → 判为重复区，记为「未判定」（**宁可漏判，不误判**）

### 为什么每个侧翼要存两套版本

85 个位点里有密集簇，最密处 25 bp 窗口内挤了 6 个位点（2390 kb 区）。
若侧翼只按参考序列取，锚本身就覆盖邻近突变位点，而真实 Group 1 菌株在邻位是 alt 等位基因，
结果锚不上 → **漏判**。实测发现这一步不做的话，体外 Group 1 基因组只有 73/85 能判出来。

因此每个侧翼额外保存一份「Group 1 单倍型版本」（把窗口内落着的其它组特异位点一并换成 alt），
两套都搜索。加上这一步后，体外 Group 1 基因组恢复到 **85/85 命中**。

### 与"整条探针完全匹配"方案的对比

| 方案 | 体外 Group 1 召回 | 说明 |
|---|---|---|
| 51-mer 探针要求完全匹配 | 偏低 | 窗口内有任何额外变异即失效 |
| **侧翼锚定 + 单碱基读取**（本工具） | **85/85** | 容忍锚之外的任何变异 |

---

## 验证结果

```bash
node scripts/validate.js
```

**27 项断言全部通过**（详见下方）。同一份 `assets/core.js` 同时被网页和命令行调用，不存在"两套实现"。

| 用例 | 输入 | 结果 |
|---|---|---|
| T1 | 参考基因组自身（FDAARGOS_1195，属 Group 2） | 非 Group 1，85/85 可判定，匹配率 0.0% |
| T2 | 体外 Group 1 基因组（85 个位点替换成 alt） | **Group 1**，85/85 可判定，匹配率 100.0% |
| T3 | T2 + 0.1% / 0.5% / 1.0% 随机噪声 | 均判为 Group 1；1% 时仍 84/85、97.6% |
| T4 | 打断成平均 30 kb / 8 kb 的 contig 并随机翻转方向 | Group 1，85/85，匹配率 100.0% |
| T5 | 只有一半位点呈 Group 1 型（模拟重组） | **不确定**（50.6%），不误报 |
| T6 | 参考基因组整条反向互补 | 非 Group 1，85/85，匹配率 0.0%（不受 contig 方向影响） |

> 注：T1 / T2 属于自洽性检验（构造数据即判定依据），**不能替代对 65 株真实基因组的回顾性验证**。
> 正式投稿建议补充这一步 —— 见下方「建议补做的验证」。

---

## 已知边界（重要）

- **输入必须是组装后的基因组**（草图 contig 完全可以），不能直接丢原始 fastq。
- 位点恰好落在 contig 末端、或侧翼 16 bp 内还有别的变异时，该位点判不出来（记为「未判定」），不影响其他位点。
- 本工具是比对流程的**近似实现**。用于发表的正式结果，请按论文流程跑
  `bwa` / `snippy` + `snp-sites`；本工具适合快速筛查和分享给同行。
- 判定只基于 85 个 SNP 位点，不含组特异基因 `group_634` 的有无检测；
  该基因的存在/缺失请按论文流程用 Roary + BLAST 判定。
- 因为 Group 2、Group 3 在「组内 100% 且组外 0%」的标准下**没有**可用位点，
  本工具只能判定 `"是否属于 Group 1"`，不能反过来做三分类。

---

## GitHub Pages 部署步骤

### 方式 A：命令行推送（常规网络）

```bash
# 1) 在本机 初始化仓库
cd cmac-group1-typer
git init
git add .
git commit -m "feat: C. macginleyi Group 1 快速鉴定工具"

# 2) 在 GitHub 上新建一个空仓库（不要勾选 Initialize with README）
#    仓库名建议：cmac-group1-typer

# 3) 关联并推送
git remote add origin https://github.com/<你的用户名>/cmac-group1-typer.git
git branch -M main
git push -u origin main
```

然后在 GitHub 仓库页面：

1. **Settings** → 左侧 **Pages**
2. **Source** 选 `Deploy from a branch`
3. **Branch** 选 `main`、目录选 `/ (root)`
4. 点 **Save**，等 1–2 分钟

### 方式 B：走 REST API（`github.com:443` 不通时用）

有些网络下 `github.com` 的 443 端口连不上（`git push` 会一直超时），但
`api.github.com` 是通的。这时用仓库自带的发布脚本，它会创建仓库、上传全部文件、
把 `main` 设为默认分支并直接打开 Pages：

```bash
# token 在 https://github.com/settings/tokens 生成，勾选 repo 权限即可
export GITHUB_TOKEN=ghp_xxxxxxxxxxxxxxxxxxxx
node scripts/publish_github.js
```

脚本幂等，重复运行只会在原有历史上追加一个新提交。

访问地址就是（脚本结束时也会打印）：

```
https://<你的用户名>.github.io/cmac-group1-typer/
```

> 不需要 GitHub Actions，不需要构建步骤 —— 页面是纯静态的 HTML/CSS/JS，推上去就能用。
> 所有资源都用相对路径引用，部署在子目录下同样正常。

### 想绑定的 doi / 补充材料

工具本身可以直接作为论文的 supplementary material 引用，建议同时把仓库打一个 tag 并
关联 Zenodo（`git tag v1.0.0 && git push --tags`），就能拿到一个长期可用的 DOI。

---

## 重新生成标记数据

换了 SNP 列表或参考基因组后：

```bash
python scripts/build_markers.py    # 重建 assets/markers.js + data/markers.json
python scripts/make_test_fasta.py  # 重建 test/ 下的自测文件
node scripts/validate.js           # 回归验证
```

`build_markers.py` 会同时做坐标校验：如果发现 85 个位点的参考碱基与参考基因组对不上，会直接报错退出，**不会静默产出错误数据**。

---

## 建议补做的验证

工具目前的自检是构造数据上的自洽性检验。要让审稿人信服，建议再补两步：

1. **回顾性验证**：拿原来那 65 株的组装序列逐个跑一遍，确认 8 株 Group 1 全部 100%、
   其余 57 株全部 0%。可以直接用命令行版批处理：

   ```bash
   node scripts/typer.js /path/to/assemblies/*.fna -o 65株回顾验证.csv --detail
   # 支持gz；--min-cov / --hi / --lo 可覆盖默认阈值
   ```

   也可以直接把整个文件夹拖进网页。

2. **与正统流程对照**：对其中若干株同时跑 `snippy` + `snp-sites`，核对 85 个位点的等位基因
   与本工具读取结果完全一致。`assets/core.js` 里的 `callLocus()` 在 Node 下可直接调用，
   便于写对照脚本。

这两步跑完，README 里的验证表可以从「自洽性检验」升级为「真实数据回顾性验证」。
