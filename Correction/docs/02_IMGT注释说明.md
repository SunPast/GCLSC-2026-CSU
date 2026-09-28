# 上游候选序列的 IMGT 编号注释

对上游交来的 `D2_candidates_50.fa` 做 IMGT 编号注释。

| 文件 | 内容 |
| --- | --- |
| `data/02_correction/D2_candidates_50/D2_candidates_50_imgt.fa` | 注释版序列，50 条，每条一行，固定 128 列 |
| `data/02_correction/D2_candidates_50/D2_candidates_50_imgt.ruler.txt` | 列号标尺，用于把列位置换算成 IMGT 位号 |
| `scripts/annotate_imgt.py` | 注释脚本 |

复现命令：

```bash
python scripts/annotate_imgt.py -i D2_candidates_50.fa \
    -o data/02_correction/D2_candidates_50/D2_candidates_50_imgt.fa
```

编号方法：ANARCI 通过 pyhmmer 后端做精确 IMGT 编号，再加一步框架校正，见第四节。

**Reference:** Dunbar, J., & Deane, C. M. (2016). "ANARCI: antigen receptor numbering and receptor classification." *Bioinformatics*, 32(2), 298–300.

---

## 一、格式约定

每条序列一行，长度 128 列，第 1 列对应 IMGT 位号 1，第 128 列对应 IMGT 位号 128。

| 符号 | 含义 |
| --- | --- |
| 字母 | 该 IMGT 位号上的残基 |
| `.` | 该 IMGT 位号上没有残基 |

`.` 覆盖两种情况，读的时候要分开看：

1. **短环空位**——IMGT 位号存在，但这条序列的环比较短，本来就没有残基。例如 CDR1 长度不到 12 个残基时，位号 31–34 就是空的。这是正常的，完整抗体也会有，通常只有 1–4 个连续的 `.`。
2. **整段缺失**——某个区域整段不在序列里，表现为开头或中间一长串连续的 `.`。这是缺陷。

表头的 `IMGT=` 字段直接给出实际覆盖的位号范围，不用自己数。

表头字段：

| 字段 | 含义 |
| --- | --- |
| `IMGT=1-128` | 这条序列实际覆盖的 IMGT 位号范围 |
| `Nterm_missing=0` | N 端缺失的位号数，即 IMGT 1 到该位号全空 |
| `Cterm_extra=0` | 超出 IMGT 128 的残余残基数，附在第 128 列之后 |
| `frame_shift=+0` | 编号框架校正量，`+0` 表示 ANARCI 原编号未动，见第四节 |
| `C23=C` / `C104=C` | 两个保守半胱氨酸在 IMGT 23 / 104 上实际读到的残基 |

---

## 二、当前批次（2026-09-22 收到的版本）

50 条，来自两个上游子批次 `BATCH_HER2_20260921_01` 27 条、`BATCH_HER2_20260921_02` 23 条。

| 覆盖范围 | 条数 |
| --- | --- |
| 1–128 | 50 |

**50 条全部从 IMGT 位号 1 起始，N 端完整。** 序列长度 117–120，起始残基 `EVQLVE` 27 条、`QVQLVE` 21 条，与天然 VH 的 N 端一致。

2 条在 IMGT 128 之后多出 1 个残基，表头记 `Cterm_extra=1`。

### 上一版（已作废）

2026-09-20 收到的同名文件里，50 条中有 38 条 N 端被砍掉 14–16 个位号，缺的是 IMGT 1–16 这一段，即 FR1 的第一条 β 折叠股。这段是 Ig 结构域 ABED β 片层的边，缺了域折不起来。该版本的分析结果保留在 `data/02_correction/D2_candidates_50_v1_superseded/`，作废原因是上游已重新生成。

对比上一批 `data/01_generated/BATCH_HER2_20260809_01/Abseqs_BATCH_HER2_20260809_01.fa` 的 25 条，也全部从 IMGT 位号 1 起始。所以那次截断是 09-20 那一版的个别退化，本版已恢复。

**Reference:**

- Lefranc, M. P., et al. (2003). "IMGT unique numbering for immunoglobulin and T cell receptor variable domains and Ig superfamily V-like domains." *Developmental & Comparative Immunology*, 27(1), 55–77.
- Chothia, C., & Lesk, A. M. (1987). "Canonical structures for the hypervariable regions of immunoglobulins." *Journal of Molecular Biology*, 196(4), 901–917.

---

## 三、当前批次的保守位点读数

从注释版直接读出 50 条在 5 个保守位号上的残基：

| IMGT 23 | IMGT 41 | IMGT 104 | IMGT 118 | IMGT 119 | 条数 |
| --- | --- | --- | --- | --- | --- |
| C | W | C | W | G | 48 |
| C | W | C | W | **W** | 1 |
| C | W | C | **Y** | G | 1 |

- **C23 与 C104：50 条全部是半胱氨酸**，没有真正的二硫键位点突变。
- **G119 异常 1 条**：`BATCH_HER2_20260921_02_SEQ_00126`，J 区读数 IMGT 118–122 为 `W-W-Q-G-T`，位号 119 是色氨酸。按规则不自动修正，判不通过。
- **F118 异常 1 条**：`BATCH_HER2_20260921_01_SEQ_00167`，位号 118 是酪氨酸，按 J-motif 规则回退为色氨酸后通过。

这两个读数与 `pipeline.py` 同批产出的 `manifest.csv` 完全一致，可互相印证。

本批没有出现游离半胱氨酸（上一版有 3 条在 IMGT 24 多一个 Cys）。

---

## 四、ANARCI 编号框架偏移与修复

### 现象

ANARCI 在部分序列上会把整条序列的位号整体偏移一格。偏移后 IMGT 23 读成 `S`、IMGT 104 读成 `Y`，表面看像是两个保守半胱氨酸同时突变，而半胱氨酸其实都在序列里。

### 受控实验

取一条完整 VH，分别删掉开头 0、14、15、16、17、18 个残基，各跑一次 ANARCI：

| 删除残基数 | ANARCI 起始位号 | ANARCI 读到的 IMGT 23 | 判定 |
| --- | --- | --- | --- |
| 0 / 14 / 16 / 17 / 18 | 正确 | C | 编号正确 |
| **15** | 偏 1 格 | **S** | **编号错位** |

删 15 个残基时 HMM 比对滑一格，整条序列位号全部偏移一位。IMGT 在位号 10 上有一个空位，删掉 `QVQLVQSGAEVKKPG` 这 15 个残基后，前端空位与内部空位叠加，是滑格的直接原因。

这个问题不限于人为截断。在框架区突变较多的序列上也会出现：本批 `BATCH_HER2_20260921_01_SEQ_00088`（FR1 起始为 `ESQLLESGG`，FR2、FR3 多处替换）同样触发了偏移。

### 修复

用两个保守半胱氨酸当锚点：一个能折叠的 V 结构域，IMGT 23 和 IMGT 104 必须是半胱氨酸。在 ±3 格范围内搜位移，取使两个位置同时落在 `C` 上的那一格。找不到就不动 ANARCI 的原编号（说明半胱氨酸是真的坏了，属于另一类问题）。

脚本里对应 `annotate_imgt.py` 的 `_cys_anchored()`。表头 `frame_shift` 记录实际用了多少格。

### 修复后校验

对当前批次 50 条逐条核对：

| 校验项 | 结果 |
| --- | --- |
| 残基无丢失、无乱序（对齐后拼回的字符串与原序列完全一致） | 50/50 通过 |
| IMGT 23 落在 Cys 上 | 50/50 |
| IMGT 104 落在 Cys 上 | 50/50 |
| 触发框架校正 | 1 条（`SEQ_00088`） |

上一版 50 条中触发校正的是 28 条，作废前同样逐条通过上述校验。

### 对 pipeline 的影响

`scripts/pipeline.py` 的 `find_conserved_positions()` 先取 ANARCI 编号，取不到再退回 motif。偏移发生时 ANARCI 给出的 IMGT 23/104 不是 Cys，于是 motif 兜底接管——按【前 40 位里第一个 Cys】和【40 位之后最后一个 Cys】定位。

**结论：最终判定不受影响。** motif 兜底落到的正是真实位点，两批共 100 条的通过/不通过与手工核对一致。

但 `02_规则说明_中文.md` 第 2.3 节写的【先用 ANARCI 做精确 IMGT 编号】，在这类序列上不成立——实际生效的是 motif 兜底。该例外已补入 `02_规则说明_中文.md` 第 2.3 节和 `02_conservation_rules.md` 第 3 节。

---

## 五、许可证

| 工具 | 许可证 | 商用 |
| --- | --- | --- |
| ANARCI | BSD-3-Clause | ✅ |
| pyhmmer | MIT | ✅ |
