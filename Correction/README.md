# Module 2：抗体 VH 序列 QC → 修正 → 审核

上游（Module 1 / 序列生成）交来一批候选 VH 氨基酸序列，本模块负责把它们收拾干净再往下传。每条序列走三步：

| 步骤 | 做什么 |
| --- | --- |
| **QC** | 找问题——非标氨基酸、5 个保守位点是否被突变、长度是否正常、CDR 里有没有质量隐患 motif |
| **修正** | 修可修的——非标氨基酸替换、被突变的保守位点回退 |
| **审核** | 给修完的序列打分，判定 **PASS** 或 **FAIL** |

最终产物是 `manifest.csv`，交给 stage 3。规则依据见 [docs/02_规则说明_中文.md](docs/02_规则说明_中文.md)。

---

## 一、环境准备

需要 Python 3.14（在 3.14.6 / Windows 11 上验证过）。

```bash
pip install -r requirements.txt
```

依赖只有两个：`anarci` 和 `pyhmmer`。选 pyhmmer 而不是原版 ANARCI 的 `hmmscan`，是因为 Windows 上没有 hmmscan 可执行文件。

**以下所有命令都在项目根目录（本文件所在目录）执行。**

---

## 二、跑一批序列

### 1. 把输入放到位

默认读这个路径：

```
data/01_generated/<BatchID>/Abseqs_<BatchID>.fa
```

比如 BatchID 是 `BATCH_HER2_20260809_01`，输入就放在
`data/01_generated/BATCH_HER2_20260809_01/Abseqs_BATCH_HER2_20260809_01.fa`。

### 2. 运行

```bash
python scripts/pipeline.py -b BATCH_HER2_20260809_01
```

### 3. 取产物

跑完屏幕会打印：

```
============================================================
Batch Processing Complete: BATCH_HER2_20260809_01
============================================================
  Total:   25
  PASS:    21
  FAIL:    4

  Manifest: C:\Users\Lenovo\Desktop\GCLSC\data\02_correction\BATCH_HER2_20260809_01\manifest.csv
  Summary:  C:\Users\Lenovo\Desktop\GCLSC\data\02_correction\BATCH_HER2_20260809_01\batch_correction_summary/batch_correction_summary.md
```

**交给 stage 3 的是 `manifest.csv`**，别的都是过程记录，自己看。

### 输入不在默认路径怎么办

上游直接甩来一个文件（比如 `D2_candidates_50.fa`）时，用 `-i` 指路径就行，`-b` 可以省：

```bash
python scripts/pipeline.py -i D2_candidates_50.fa
```

BatchID 自动从文件名推出来（`Abseqs_` 前缀和小写扩展名会被去掉）：`D2_candidates_50.fa` → 批目录 `data/02_correction/D2_candidates_50/`。

两个都给也行，`-b` 优先：

```bash
python scripts/pipeline.py -b D2_candidates_50 -i D2_candidates_50.fa
```

---

## 三、输入格式

标准 FASTA，一次一批：

```
>BATCH_EXAMPLE_001_SEQ_00001|Score=0.92|Length=111
EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTY
YADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDYWGQGTLVTVSS
```

- 序列可以折行，也可以一行到底，两种都认。
- `Score=`、`Length=` 这些字段本模块不读，保留原样即可。

序列 ID 按下面的顺序取，取到就停：

1. header 里带 `Global_Seq_ID=` 字段 → 用它
2. 否则找**键名里含 `SEQ`** 的字段 → 用它
3. 都没有 → 取 header 里第一个 `|` 之前的内容

上面例子走的是第 3 条。建议直接用第 3 种的写法，header 里不重复放 `Global_Seq_ID=`。

---

## 四、输出

```
data/02_correction/<BatchID>/
├── manifest.csv                        ← 交付给 stage 3
├── batch_correction_summary/
│   └── batch_correction_summary.md     ← 整批汇总，人看的
└── <Seq_ID>/
    ├── raw.fa                          ← 原始序列
    ├── processed/
    │   ├── v01.fa                      ← 修正后的序列（有修正才生成）
    │   └── final.fa                    ← 最终序列（判 PASS 才生成）
    ├── logs/
    │   ├── correction_log.json         ← 改了哪些位点、改前改后是什么
    │   └── audit_log.json              ← 7 项得分明细、总分、判定
    └── reports/
        └── structure_report.md         ← 单条序列的完整报告
```

`processed/` 里有什么取决于这条序列的情况：

| 情况 | 生成的文件 |
| --- | --- |
| 本来就干净，PASS | `final.fa` |
| 有修正，PASS | `v01.fa` + `final.fa` |
| 有修正，FAIL | `v01.fa`（没有 `final.fa`） |
| 无修正，FAIL | 都不生成 |

**`manifest.csv`** 三列，一行一条，含 FAIL 的：

| 列 | 内容 |
| --- | --- |
| `Global_Seq_ID` | 序列 ID |
| `Sequence` | **修正后**的序列（不是原始序列） |
| `Status` | `PASS` 或 `FAIL` |

FAIL 的行保留在表里，不删——出问题的序列要能追溯到，不然 stage 3 只会看到"少了 1 条"。

---

## 五、PASS / FAIL 怎么判的

先看**硬性失败条件**，中任一条直接 FAIL，不进打分：

| # | 条件 | 原因 |
| --- | --- | --- |
| 1 | G119 突变 | 紧转角，改了没法预测，不可修正 |
| 2 | 长度 < 90 或 > 140 | 不是有效的 V 结构域 |
| 3 | 非抗体特征——FR2 的 W-motif 和 FR4 的 J-motif 都找不到 | 输入压根不像抗体序列 |
| 4 | C23 和 C104 同时缺失 | 二硫键没了，Ig 折叠起不来 |
| 5 | C23 或 C104 任一仍缺失 | 单条二硫键断裂 |

没中硬性条件就按 7 项加权打分（满分 1.0），**≥ 0.94 判 PASS**：

| 项 | 权重 |
| --- | --- |
| C23 完整 | 25% |
| C104 完整 | 25% |
| 无非标氨基酸 | 15% |
| W41 完整 | 10% |
| F118 完整 | 10% |
| 长度正常（95–135 满分，90–94 或 136–140 减半） | 10% |
| 无不可修正问题 | 5% |

**修正后结构合规的位点按满分计**——修正与否只记进报告，不影响得分。所以一条被修好的序列和一条本来就干净的序列，得分是一样的。

阈值 0.94 是从数据里校准出来的——正负样本扫描后，**误杀为 0 的最高阈值**（正样本最低分 0.950）。

本模块**只有 PASS 和 FAIL 两种状态**，没有中间档、没有 WARN 状态写进 manifest。

---

## 六、附加工具

### `scripts/annotate_imgt.py` — 生成 IMGT 编号注释版

把序列摆到固定的 IMGT 1–128 框架上，缺位用 `.` 占位。用于人工核对保守位点读数。

```bash
python scripts/annotate_imgt.py -i D2_candidates_50.fa \
    -o data/02_correction/D2_candidates_50/D2_candidates_50_imgt.fa
```

输出每条一行、固定 128 列，header 里带 `IMGT=`、`Nterm_missing=`、`C23=`、`C104=` 等字段。格式说明和 ANARCI 编号偏移问题见 [docs/02_IMGT注释说明.md](docs/02_IMGT注释说明.md)。

### `scripts/generate_mock.py` — 造示例输入

新克隆下来想试跑，但手上没有序列时用这个：生成 25 条覆盖各类缺陷的序列，写到默认输入路径，接着就能直接跑 `pipeline.py -b BATCH_HER2_20260809_01`。

```bash
python scripts/generate_mock.py
```

### `scripts/validation/` — 验证脚本（一般用不到）

v1–v4 四套验证脚本和正样本存档脚本，产物是 `docs/validation_report_v1–v4.md`。

**这些脚本会重新生成对应的报告文档。** 报告是定稿的，别随手重跑。

---

## 七、目录结构

```
<项目根目录>/
├── README.md                  ← 本文件
├── requirements.txt
├── scripts/
│   ├── pipeline.py            ← 主流程，就这一个入口
│   ├── anarci_numbering.py    ← ANARCI + pyhmmer 的 IMGT 编号封装
│   ├── annotate_imgt.py       ← IMGT 注释版生成
│   ├── generate_mock.py       ← 示例输入生成
│   └── validation/            ← v1–v4 验证脚本，一般不动
├── docs/
│   ├── 02_规则说明_中文.md     ← 规则依据（中文，主文档）
│   ├── 02_conservation_rules.md
│   ├── 02_IMGT注释说明.md
│   └── validation_report_v1–v4.md
└── data/                      ← 本地目录，整个不入库
    ├── 01_generated/          ← 输入：generate_mock.py 生成，或上游给
    └── 02_correction/         ← 输出
```

**`data/` 整个目录不入库**（仓库 `.gitignore` 里已屏蔽）。要跑主流程，用 `generate_mock.py` 造输入，或者直接把上游的 FASTA 丢进来即可。

`scripts/validation/` 里的脚本额外需要 `data/*.json`——从 PDB、IEDB 抓下来的公开抗体序列，是当初做验证时缓存下来的。这部分同样不入库，所以**验证脚本在新克隆的仓库里跑不起来**，需要那份缓存才能复现 `docs/validation_report_*.md`。主流程 `pipeline.py` 不受影响。

---

## 八、参考文档

- [docs/02_规则说明_中文.md](docs/02_规则说明_中文.md) — 规则、修正策略、阈值校准过程
- [docs/02_IMGT注释说明.md](docs/02_IMGT注释说明.md) — IMGT 注释格式，ANARCI 编号框架偏移与修复
- [docs/validation_report_v1.md](docs/validation_report_v1.md) — 端到端验证
- [docs/validation_report_v4.md](docs/validation_report_v4.md) — 阈值扫描

**文献：**

- Lefranc, M. P., et al. (2003). IMGT unique numbering for immunoglobulin and T cell receptor variable domains and Ig superfamily V-like domains. *Developmental & Comparative Immunology*, 27(1), 55–77.
- Dunbar, J., & Deane, C. M. (2016). ANARCI: antigen receptor numbering and receptor classification. *Bioinformatics*, 32(2), 298–300.
- Chothia, C., & Lesk, A. M. (1987). Canonical structures for the hypervariable regions of immunoglobulins. *Journal of Molecular Biology*, 196(4), 901–917.
