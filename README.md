# AI4S Drug Design Collaboration

> 第一届全球大学生生命科学挑战赛（CSU 2024级生信）AI 大分子与多肽药物设计赛道协作仓库

本项目聚焦 **AI 驱动的抗体序列设计与优化**，以 **ERBB2（HER2）** 为当前研究靶点，探索生成式人工智能、潜空间建模与强化学习在抗体药物设计中的应用，构建从序列表征、性质学习、序列生成到计算验证的完整研究流程。

-----

## 1. 项目目标

本项目拟构建一个基于隐空间解耦（Latent-space Disentanglement）与潜空间扩散强化学习（Latent Diffusion + Reinforcement Learning）的抗体设计平台，实现：

- 高质量抗体序列表征学习
- 可解释的潜空间建模
- 抗体序列可控生成
- 基于结合能反馈的强化学习优化
- 面向 HER2 靶点的计算设计与验证

整体采用目前蛋白质基础模型广泛使用的 **Pretrain → Finetune** 范式，充分利用大规模无标签数据与有限标注数据，实现高效迁移学习。

-----

## 2. 研究概览

### 2.1 🧬 当前研究对象

- **Target Gene:** ERBB2
- **Protein:** HER2 (Human Epidermal Growth Factor Receptor 2)

HER2 是经典肿瘤治疗靶点，在乳腺癌、胃癌等多种恶性肿瘤中过表达，具有丰富的抗体结构、亲和力及实验数据，为 AI 模型训练与计算验证提供了可靠的数据基础。



### 2.2 🏗️ 技术路线

整体框架如下：

```text
Protein Sequence
        │
        ▼
ESM Protein Encoder
        │
        ▼
Transformer-VAE
(Latent Disentanglement)
        │
        ▼
Latent Representation
        │
        ├────────► Property Prediction
        │
        ▼
Latent Diffusion Transformer
        │
        ▼
Reinforcement Learning
        │
        ▼
Antibody Optimization
        │
        ▼
Docking / Structure Prediction
```


整个模型采用三阶段训练策略。

**Stage I：无监督预训练**

利用 ESM 提取蛋白序列表示，并通过 Transformer-VAE 学习具有可解释性的潜空间。

模型设计半显式潜空间：

- 显式潜变量负责学习具有明确生物学意义的属性；
- 剩余潜变量保留复杂序列信息，提高生成能力。

完成训练后，通过以下方式验证潜空间质量：

- UMAP / t-SNE 聚类可视化；
- 随机森林分类及特征重要性分析；
- 下游分类器性能评估。



**Stage II：连续性质继续预训练**

引入带标签数据进行继续预训练，预测多个连续生物物理性质，包括：

- Protein Fitness
- Thermal Stability
- GB1 Fitness

采用多任务学习，使潜空间同时学习：

- Sequence Representation
- Physicochemical Property
- Energy-aware Representation

进一步提升表示能力。



**Stage III：监督微调**

冻结 VAE，仅训练潜空间扩散模型（Latent Diffusion Transformer）。

扩散过程引入预训练条件表示进行条件生成，并利用强化学习不断优化生成策略。

奖励函数综合考虑：

- Binding Affinity
- Sequence Diversity
- Developability
- Sequence Validity

实现针对 HER2 的抗体序列优化。



### 2.3 🚀 奖励预测器设计

强化学习若直接调用分子对接计算奖励，将产生极高计算成本。

因此设计轻量级 **Mamba Reward Predictor**：

1. 真实结合数据
2. Reward Predictor
3. 预测奖励
4. RL 更新

仅当预测器认为候选序列具有较高潜力时，再调用物理引擎（Docking）进行最终验证。

这一设计能够显著降低强化学习训练成本，提高有限算力条件下的训练效率。



### 2.4 ✨ 方法特点

- **高可解释性**：显式潜变量提高模型可解释性。
- **高表达能力**：潜空间扩散突破传统 VAE 高斯先验限制。
- **高生成质量**：扩散模型结合强化学习，实现抗体序列持续优化。
- **工程可行性**：分阶段训练降低峰值算力需求，更适合学生团队实施。



### 2.5 📚 数据集

**预训练**

- FLIP Benchmark
- UniProt（可扩展）
- OAS（可扩展）

**连续性质学习**

- Protein Fitness
- Thermal Stability
- GB1

**监督微调**

- AB-Bind
- SabDab-Affinity
- HER2 相关公开抗体数据（可补充）



### 2.6 📊 模型评估

**阶段一：潜空间评估**

- UMAP / t-SNE 聚类
- Batch Effect 分析
- 随机森林分类
- 特征重要性分析

**阶段二：生成质量评估**

与以下模型进行比较：

- RFdiffusion
- AbDiffuser
- ESM3

评价指标包括：

- Binding Energy
- Diversity
- Novelty
- Developability

并结合 AlphaFold3 或 tFold 预测抗原-抗体复合物结构。

**阶段三：闭环验证**

针对已有 HER2 抗体：

- 保持抗原不变；
- 重新设计抗体；
- 比较设计结果与天然抗体的结合能力和结构稳定性。

若条件允许，将进一步开展湿实验验证；否则采用成熟计算流程完成干实验验证。

-----



### 3. 协作理念

- 兴趣驱动，自主选择研究模块。
- 鼓励合理使用 AI 工具提升研发效率。
- 重视学习过程与团队成长。
- 尊重成员学习与生活安排，关键节点集中协作。

-----



### 4. 项目进度

| 节点 | 时间 |
| :--- | :--- |
| 报名与摘要提交 | 9 月 5 日前 |
| 项目成果提交 | 10 月 8 日前 |
| 区赛评审 | 10 月 15 日 |
| 全国总决赛 | 10 月下旬 |


-----



### 5. 联系方式

项目讨论以群聊为主，也欢迎通过 issue 或 PR 参与GitHub仓库建设。

也欢迎邮件联系：

- 📧 2210240103@csu.edu.cn

- 📧 jiading682@qq.com

------

> **Disclaimer**
>
> 本项目仅用于算法研究与学术交流，不作为任何临床医疗建议。所有生成结果均需经过进一步的计算验证及实验验证。
