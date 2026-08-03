# AI4S Drug Design Collaboration

> 第一届全球大学生生命科学挑战赛（CSU 2024级生信）AI 大分子与多肽药物设计赛道协作仓库

本项目聚焦 **AI 驱动的抗体序列设计与优化**，以 **ERBB2（HER2）** 为当前研究靶点，探索生成式建模、潜空间扩散与强化学习在抗体设计中的应用，构建从序列表征、可控生成到计算验证的研究流程。

---

## 1. 项目目标

本项目拟构建一个基于 **潜空间表征学习** 与 **潜空间条件扩散 + 强化学习** 的抗体设计平台，实现：

- 高质量抗体序列表征学习
- 可用于下游生成的潜空间建模
- 抗体序列可控生成
- 基于结合相关反馈的强化学习优化
- 面向 HER2 靶点的计算设计与验证

整体遵循 **大规模无监督预训练 → 条件生成 → 任务导向微调** 的分阶段范式：先在大规模抗体序列上学习表征与生成先验，再引入属性条件与有限标注信号，逐步面向靶点优化。

---

## 2. 研究概览

### 2.1 🧬 当前研究对象

- **Target Gene:** ERBB2
- **Protein:** HER2 (Human Epidermal Growth Factor Receptor 2)

HER2 是经典肿瘤治疗靶点，在乳腺癌、胃癌等多种恶性肿瘤中过表达，具有丰富的抗体结构、亲和力及实验数据，为 AI 模型训练与计算验证提供了可靠的数据基础。

### 2.2 🏗️ 技术路线

整体框架如下：

```text
Antibody Sequence
        │
        ▼
Sequence VAE
(Latent Representation)
        │
        ▼
Latent Diffusion
(Conditional Generation)
        │
        ▼
Reinforcement Learning
(Policy Fine-tuning)
        │
        ▼
Antibody Optimization
        │
        ▼
Docking / Structure Prediction
```

整个模型采用三阶段训练策略。

**Stage I：抗体序列 VAE 预训练**

在大规模无标签抗体序列（当前以 OAS 为主，并补充公开 BCR 数据）上训练序列 VAE，学习可供重建与采样的潜空间表示。训练完成后冻结编码器–解码器，作为后续阶段的表征与解码模块。

潜空间质量验证包括（但不限于）：

- 重建质量与抗体框架保守性检查；
- 潜空间分布与利用率分析；
- 先验采样序列的基本生物学合理性评估。

**Stage II：潜空间条件扩散**

冻结 Stage I 的 VAE，在潜空间上训练条件扩散模型，学习潜变量分布，并按连续生物物理 / 可开发性相关属性进行条件生成。属性信号用于引导采样，而非在本阶段重新做大规模序列表征预训练。

本阶段完成后，应能够从噪声出发采样潜变量，并经冻结解码器得到候选抗体序列。

**Stage III：强化学习微调**

在 Stage II 扩散策略基础上，以强化学习（计划采用扩散策略优化一类方法）对生成过程做任务导向微调，面向 HER2 等相关目标提升候选质量。

奖励函数拟综合考虑：

- Binding Affinity（或可代理的结合相关分数）
- Sequence Diversity
- Developability
- Sequence Validity

实现针对靶点场景的抗体序列优化。

### 2.3 🚀 奖励预测器设计

强化学习若每一步直接调用分子对接计算奖励，计算成本过高。

因此计划引入轻量级 **Reward Predictor**（架构候选包括 Mamba 等序列模型）：

1. 真实结合 / 亲和力相关数据
2. Reward Predictor
3. 预测奖励
4. RL 更新

当预测器认为候选序列具有较高潜力时，可调用物理引擎（Docking）等进行最终校验，以控制训练成本。

### 2.4 ✨ 方法特点

- **分阶段解耦**：表征、条件生成与策略优化分开训练，降低联调难度与峰值算力需求。
- **潜空间生成**：在压缩潜空间上做扩散，避免直接在离散序列空间上承担全部生成负担。
- **条件可控**：属性条件用于引导生成，服务于可开发性与设计约束。
- **反馈可扩展**：以可学习奖励近似昂贵物理评估，并保留对接 / 结构预测作为后验校验。



### 2.5 📚 数据集

**Stage I 预训练**

- OAS（Observed Antibody Space；当前主力）
- 公开 BCR / 单细胞抗体序列（补充）

**Stage II 条件信息**

- 由序列计算或标注得到的连续性质（如等电点、疏水性、稳定性相关指标等；具体集合随实验迭代）

**Stage III 与奖励学习**

- AB-Bind
- SabDab-Affinity
- HER2 相关公开抗体数据（可补充）



### 2.6 📊 模型评估

**阶段一：表征与重建评估**

- 重建质量
- 抗体框架 / 保守位点合理性
- 潜空间结构与采样质量
- 必要的可视化与探针分析（如 UMAP / t-SNE）

**阶段二：生成质量评估**

与公开生成式基线进行比较（视可复现性选用），例如：

- RFdiffusion
- AbDiffuser
- ESM3

评价指标包括：

- Binding Energy（或等价代理指标）
- Diversity
- Novelty
- Developability

并结合 AlphaFold3 或 tFold 等工具评估抗原–抗体复合物结构（以实际可运行工具为准）。

**阶段三：闭环验证**

针对已有 HER2 抗体：

- 保持抗原不变；
- 重新设计抗体；
- 比较设计结果与天然抗体的结合能力与结构稳定性。

若条件允许，将进一步开展湿实验验证；否则采用成熟计算流程完成干实验验证。

---



### 3. 协作理念

- 兴趣驱动，自主选择研究模块。
- 鼓励合理使用 AI 工具提升研发效率。
- 重视学习过程与团队成长。
- 尊重成员学习与生活安排，关键节点集中协作。

---



### 4. 项目进度


| 节点      | 时间        |
| ------- | --------- |
| 报名与摘要提交 | 9 月 5 日前  |
| 项目成果提交  | 10 月 8 日前 |
| 区赛评审    | 10 月 15 日 |
| 全国总决赛   | 10 月下旬    |


---



### 5. 联系方式

项目讨论以群聊为主，也欢迎通过 issue 或 PR 参与GitHub仓库建设。

也欢迎邮件联系：

- 📧 [2210240103@csu.edu.cn](mailto:2210240103@csu.edu.cn)
- 📧 [jiading682@qq.com](mailto:jiading682@qq.com)

---



### 6. 快速开始

权重文件体积较大，**不包含在本仓库中**。请先通过上述邮箱联系获取，并放置到指定槽位：

```text
DiffEboAb/checkpoints/stage1_vae/set.pt      # 主模型（set-latent）
DiffEboAb/checkpoints/stage1_vae/model.pt    # 指向 set.pt 的默认加载名（可为同名文件或符号链接）
DiffEboAb/checkpoints/stage1_vae/vec.pt      # 可选：vector-latent 对照
```

依赖：Python ≥ 3.10。安装：

```bash
pip install -r requirements.txt
# GPU 建议按 https://pytorch.org 选择对应 CUDA 的 torch 轮子
```

在 `DiffEboAb/` 目录下执行无条件先验采样：

```bash
cd DiffEboAb

# -o / --out 为必填；默认加载 checkpoints/stage1_vae/model.pt
python scripts/stage1_vae/generate.py -o outputs/uncond.fasta --n 1000

# 指定对照权重或设备
python scripts/stage1_vae/generate.py \
  -o outputs/uncond_vec.fasta \
  --ckpt checkpoints/stage1_vae/vec.pt \
  --n 100 \
  --device cuda
```

运行单测（可选）：

```bash
cd DiffEboAb
python -m pytest tests/ -q
```

---

> **Disclaimer**
>
> 本项目仅用于算法研究与学术交流，不能直接作为任何临床医疗建议。所有生成结果均需经过进一步的计算验证及实验验证。

