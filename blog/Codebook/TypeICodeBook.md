# NR Type I 码本设计原理

> 本文把 MIMO-OFDM 信道、信道估计、双极化阵列和 NR Type I 码本放在同一个数学框架中，重点推导 Type I 码本如何用有限索引表示下行预编码矩阵。背景和设计原则只保留必要内容，重点放在单层/多层 Type I 预编码矩阵、宽带/子带 PMI 获取以及量化误差建模。

## 目录

- [1 背景、模型与问题定义](#1-背景模型与问题定义)
  - [1.1 MIMO-OFDM 预编码模型](#11-mimo-ofdm-预编码模型)
  - [1.2 为什么需要码本反馈](#12-为什么需要码本反馈)
  - [1.3 理想预编码与量化目标](#13-理想预编码与量化目标)
- [2 Notation](#2-notation)
- [3 Type I 码本设计原则](#3-type-i-码本设计原则)
- [4 单层 Type I 预编码矩阵表达](#4-单层-type-i-预编码矩阵表达)
  - [4.1 空间波束选择](#41-空间波束选择)
  - [4.2 双极化合成](#42-双极化合成)
  - [4.3 宽带和子带相位](#43-宽带和子带相位)
  - [4.4 单层码本量化误差](#44-单层码本量化误差)
- [5 多层 Type I 预编码矩阵表达](#5-多层-type-i-预编码矩阵表达)
  - [5.1 多层传输与秩指示](#51-多层传输与秩指示)
  - [5.2 多波束/多极化层的矩阵构造](#52-多波束多极化层的矩阵构造)
  - [5.3 归一化与功率约束](#53-归一化与功率约束)
- [6 宽带与子带预编码矩阵获取](#6-宽带与子带预编码矩阵获取)
  - [6.1 宽带预编码矩阵获取](#61-宽带预编码矩阵获取)
  - [6.2 子带预编码矩阵获取](#62-子带预编码矩阵获取)
  - [6.3 宽带方向与子带补充的组合](#63-宽带方向与子带补充的组合)
  - [6.4 gNB 侧如何使用反馈得到的预编码矩阵](#64-gnb-侧如何使用反馈得到的预编码矩阵)
  - [6.5 一个端到端获取示例](#65-一个端到端获取示例)
- [7 Type I 码本与信道建模的关系](#7-type-i-码本与信道建模的关系)
- [8 小结](#8-小结)
- [9 参考文献](#9-参考文献)

# 1 背景、模型与问题定义

NR 下行 MIMO-OFDM 预编码的核心问题是：gNB 如何根据 UE 侧看到的下行信道，选择一个发射预编码矩阵，使多天线发送信号落在有利的空间方向、极化方向和频域子带上。Type I 码本不是为了完整反馈信道矩阵，而是为了用少量 RI/PMI/CQI 比特表示一个可用于下行传输的预编码矩阵。

## 1.1 MIMO-OFDM 预编码模型

NR Type I 码本的预编码矩阵定义在 CSI-RS antenna ports 维度上。也就是说，码本矩阵的行数对应 CSI-RS antenna ports 数，而不是直接对应基站实现中的每一个物理天线单元。设用于码本反馈的 CSI-RS 端口数为：

$$
P_{\mathrm{CSI-RS}}=N_{\mathrm{AP}}.
$$

固定第 $m$ 个 OFDM 符号和第 $k$ 个子载波，UE 从 CSI-RS 端口看到的频域信道为：

$$
\mathbf H_{\mathrm{CSI}}[m,k]\in\mathbb C^{N_r\times P_{\mathrm{CSI-RS}}}.
$$

若 PDSCH 传输 $\nu$ 层数据，层符号为 $\mathbf s[m,k]\in\mathbb C^{\nu\times1}$，Type I 码本生成的预编码矩阵为：

$$
\mathbf W[m,k]\in\mathbb C^{P_{\mathrm{CSI-RS}}\times \nu}.
$$

它把 PDSCH 层，也等价地把对应的 PDSCH DM-RS antenna ports，映射到 CSI-RS antenna ports：

$$
\mathbf x[m,k]=\mathbf W[m,k]\mathbf s[m,k],
$$

$$
\mathbf y[m,k]
=
\mathbf H_{\mathrm{CSI}}[m,k]\mathbf W[m,k]\mathbf s[m,k]
+\mathbf w[m,k].
$$

其中 $\nu\leq\min(P_{\mathrm{CSI-RS}},N_r)$。预编码目标是让等效信道 $\mathbf H_{\mathrm{CSI}}[m,k]\mathbf W[m,k]$ 具有更高增益、更低层间干扰和更好的条件数。

从参考信号关系看，CSI-RS 用于让 UE 测量 $\mathbf H_{\mathrm{CSI}}[m,k]$ 并完成 PMI/RI/CQI 选择；PDSCH DM-RS 与 PDSCH 数据经历相同的 $\mathbf W[m,k]$，因此 UE 解调 PDSCH 时通常估计的是等效信道 $\mathbf H_{\mathrm{CSI}}[m,k]\mathbf W[m,k]$，不需要把 $\mathbf H_{\mathrm{CSI}}$ 和 $\mathbf W$ 分别解出来。

如果把一个 OFDM 符号上的所有子载波合并起来，整体频域信道是块对角矩阵：

$$
\mathbf H_{\mathrm{OFDM}}[m]
=
\mathrm{blkdiag}
\left(
\mathbf H_{\mathrm{CSI}}[m,0],
\mathbf H_{\mathrm{CSI}}[m,1],
\cdots,
\mathbf H_{\mathrm{CSI}}[m,N-1]
\right)
\in
\mathbb C^{NN_r\times NP_{\mathrm{CSI-RS}}}.
$$

因此，本文后续的 $\mathbf H_{\mathrm{CSI}}[m,k]$ 始终指某个资源粒子上从 CSI-RS antenna ports 到 UE 接收天线的 MIMO 信道矩阵；宽带和子带预编码是在多个 $k$ 上共同选择一个或一组码本矩阵。

## 1.2 为什么需要码本反馈

在 FDD 下行中，上下行频点不同，gNB 不能直接用上行信道推断完整下行信道。在 TDD 中，虽然可以利用上下行互易性，但当 UE 上行覆盖受限、UE 接收天线数多于发送天线数、SRS 开销较大或上下行校准误差不可忽略时，仍可能需要下行 CSI-RS 辅助获取下行预编码方向。

典型 Type I 反馈闭环为：

```mermaid
flowchart LR
    A[gNB 发送 CSI-RS] --> B[UE 测量 CSI-RS]
    B --> C[UE 估计 CSI-RS 端口信道 H_CSI]
    C --> D[UE 在 Type I 码本中搜索候选预编码矩阵]
    D --> E[UE 选择 RI、PMI 和 CQI]
    E --> F[UE 通过上行控制信道反馈 CSI]
    F --> G[gNB 根据 PMI 重建预编码矩阵 W]
    G --> H[gNB 使用 W 进行 PDSCH 下行预编码]
```

也就是说，UE 不反馈完整 $\hat{\mathbf H}_{\mathrm{CSI}}[m,k]$，而是在标准码本集合中选择 RI、PMI 和 CQI；gNB 根据 PMI 使用同一套码本规则重建 $\mathbf W$。Type I 码本适合低开销 SU-MIMO 反馈，也可以作为 MU-MIMO 调度输入之一，但它本身不是多用户零迫或块对角化预编码算法。gNB 可以参考 UE 推荐的 PMI，但调度器并不被强制使用 UE 推荐的预编码矩阵。

## 1.3 理想预编码与量化目标

如果没有码本限制，单用户 MIMO 的理想发送方向来自信道右奇异子空间。对某个子载波或窄带等效信道：

$$
\mathbf H
=
\mathbf U\mathbf \Sigma\mathbf V^H,
\qquad
\mathbf W_{\mathrm{opt}}
=
\mathbf V_\nu,
$$

其中 $\mathbf V_\nu$ 是 $\mathbf V$ 的前 $\nu$ 列。此时：

$$
\mathbf H\mathbf W_{\mathrm{opt}}
=
\mathbf U\mathbf \Sigma\mathbf V^H\mathbf V_\nu
=
\mathbf U_\nu\mathbf \Sigma_\nu.
$$

宽带或子带预编码不是为每个子载波独立反馈一个 $\mathbf W[m,k]$，而是在一个子载波集合上选择代表性的发射端子空间。例如第 $b$ 个子带可构造：

$$
\mathbf R_b
=
\sum_{k\in\mathcal K_b}
\hat{\mathbf H}_{\mathrm{CSI}}^H[m,k]\hat{\mathbf H}_{\mathrm{CSI}}[m,k],
\qquad
\mathbf W_{\mathrm{opt},b}
=
\mathrm{eig}_\nu(\mathbf R_b).
$$

Type I 码本把连续的 $\mathbf W_{\mathrm{opt}}$ 或 $\mathbf W_{\mathrm{opt},b}$ 限制到有限集合 $\mathcal W$ 中：

$$
\boxed{
\mathbf W^\star
=
\arg\max_{\mathbf W_i\in\mathcal W}
\mathcal M(\hat{\mathbf H}_{\mathrm{CSI}},\mathbf W_i)
}
$$

常见度量包括等效信道增益 $\|\hat{\mathbf H}_{\mathrm{CSI}}\mathbf W_i\|_F^2$、子空间距离或互信息近似：

$$
\mathcal M(\hat{\mathbf H}_{\mathrm{CSI}},\mathbf W_i)
=
\log_2\det
\left(
\mathbf I+
\frac{\rho}{L}
\hat{\mathbf H}_{\mathrm{CSI}}\mathbf W_i\mathbf W_i^H\hat{\mathbf H}_{\mathrm{CSI}}^H
\right).
$$

因此，本文后面的重点不是再展开一般 MIMO 预编码理论，而是推导 Type I 码本如何构造 $\mathcal W$，以及 UE 如何在 $\mathcal W$ 中选择宽带/子带 PMI。

# 2 Notation

| 符号 | 含义 |
|---|---|
| $m$ | OFDM 符号索引 |
| $k$ | OFDM 子载波索引 |
| $b$ | 子带索引 |
| $N$ | OFDM FFT 点数或子载波数 |
| $N_t$ | gNB 实现侧发射维度；本文只在需要区分物理阵列实现时使用 |
| $N_r$ | UE 接收天线端口数 |
| $P_{\mathrm{CSI-RS}}$ | 用于码本反馈的 CSI-RS antenna ports 数 |
| $N_{\mathrm{AP}}$ | 论文中使用的 antenna port 数，本文取 $N_{\mathrm{AP}}=P_{\mathrm{CSI-RS}}$ |
| $N_1,N_2$ | 发射端二维天线面板在两个空间维度上的端口数 |
| $O_1,O_2$ | 2D DFT 码本在两个维度上的过采样因子 |
| $N_t^{\mathrm{pos}}$ | 同一个极化方向的阵列位置数 |
| $N_t^{\mathrm{pol}}$ | 发射端极化方向数，双极化时通常为 2 |
| $\nu$ | 传输层数，也称 rank，由 RI 指示或约束 |
| $L$ | 码本中报告/组合的 beam 数；若讨论 rank 时本文优先使用 $\nu$ |
| $\mathbf H_{\mathrm{CSI}}[m,k]$ | CSI-RS antenna ports 到 UE 接收端口的 MIMO CFR，$\mathbb C^{N_r\times P_{\mathrm{CSI-RS}}}$ |
| $\hat{\mathbf H}_{\mathrm{CSI}}[m,k]$ | UE 基于 CSI-RS 得到的下行信道估计 |
| $\mathbf x[m,k]$ | 预编码后的 CSI-RS antenna ports 维度发射向量，$\mathbb C^{P_{\mathrm{CSI-RS}}\times 1}$ |
| $\mathbf s[m,k]$ | 预编码前的第 $m$ 个符号、第 $k$ 个子载波的层符号向量，$\mathbb C^{\nu\times 1}$ |
| $\mathbf y[m,k]$ | 接收向量，$\mathbb C^{N_r\times 1}$ |
| $\mathbf W[m,k]$ | 码本预编码矩阵，$\mathbb C^{P_{\mathrm{CSI-RS}}\times \nu}$ |
| $\mathcal W$ | 码本集合 |
| $\mathcal W_{\mathrm{WB}}$ | 宽带 Type I 码本候选集合 |
| $\mathcal W_{\mathrm{SB}}(i_{\mathrm{WB}})$ | 给定宽带 PMI 后的子带补充候选集合 |
| $\mathbf W_i$ | 码本中的第 $i$ 个候选预编码矩阵 |
| $\mathbf W_{\mathrm{WB}}$ | 宽带预编码矩阵 |
| $\mathbf W_b$ | 第 $b$ 个子带的预编码矩阵 |
| $i_{\mathrm{WB}}$ | 宽带 PMI |
| $i_{\mathrm{SB},b}$ | 第 $b$ 个子带的子带 PMI |
| $N_{\mathrm{cb}}$ | 码本候选预编码矩阵数量 |
| $\epsilon_{\mathrm{q}}$ | 预编码方向量化误差 |
| RI | Rank Indicator，秩指示 |
| PMI | Precoding Matrix Indicator，预编码矩阵指示 |
| CQI | Channel Quality Indicator，信道质量指示 |
| $\mathbf g_{m_1},\mathbf u_{m_2}$ | 水平和垂直方向的 DFT beam 向量 |
| $\mathbf w_{m_1,m_2}$ | 2D DFT beam，通常由 $\mathbf u_{m_2}\otimes\mathbf g_{m_1}$ 得到 |
| $i_{1,1},i_{1,2}$ | Type I PMI 中的宽带 beam indicators，用于确定 $m_1,m_2$ |
| $i_2$ | Type I PMI 中的 phase/co-phasing indicator |
| $\phi_n$ | 双极化 co-phasing 因子，常见形式为 $e^{j\pi n/2}$ |
| $\mathbf U,\mathbf \Sigma,\mathbf V$ | 信道矩阵 SVD 分解中的左奇异矩阵、奇异值矩阵、右奇异矩阵 |
| $\sigma_i$ | 第 $i$ 个奇异值 |
| $\mathbf V_\nu$ | 最优 $\nu$ 层预编码方向，由前 $\nu$ 个右奇异向量组成 |
| $\|\cdot\|_F$ | Frobenius 范数 |
| $(\cdot)^H$ | 共轭转置 |
| $\otimes$ | Kronecker 积 |

# 3 Type I 码本设计原则

Type I 码本的设计原则可以压缩成一句话：

$$
\boxed{
\text{用有限索引表示 DFT 空间波束、双极化合成、rank 和宽带/子带补充信息}
}
$$

它不是直接量化 $\mathbf H_{\mathrm{CSI}}[m,k]$，而是量化预编码方向：

$$
\hat{\mathbf H}_{\mathrm{CSI}}
\rightarrow
\mathbf W_{\mathrm{opt}}
\approx
\mathbf W(i_{1,1},i_{1,2},i_2,\dots)
\rightarrow
\mathrm{PMI}.
$$

PMI 本质上是码本生成函数的参数：

$$
\boxed{
\mathrm{PMI}
\longleftrightarrow
(i_{1,1},i_{1,2},i_2,\mathrm{RI},\text{subband index},\dots)
}
$$

对 Type I single-panel，PMI 可以理解为 $I=[I_1,I_2]$：$I_1$ 主要报告 2D DFT beam selection，$I_2$ 主要报告 dual-polarization co-phasing。gNB 和 UE 都知道标准定义的码本生成规则，因此 UE 只需反馈索引，gNB 即可重建预编码矩阵：

$$
\boxed{
\mathbf W
=
f_{\mathrm{codebook}}(\mathrm{PMI})
}
$$

主要参数及作用如下。

| 参数 | 作用 | 对性能和开销的影响 |
|---|---|---|
| $N_1,N_2$ | 面板两个维度的端口数 | 决定阵列孔径和空间分辨率 |
| $O_1,O_2$ | DFT 过采样因子 | 增大候选波束数量，降低角度量化误差 |
| $\Phi$ 或 $i_2$ | co-phasing 候选集合/索引 | 增强两个极化端口的合成能力 |
| $\nu_{\max}$ | 最大支持 rank | 支持更高层数，但需要更多 RI/PMI 搜索 |
| 子带数量 $N_{\mathrm{SB}}$ | 频域反馈粒度 | 子带越细，频选适配越好，但反馈越大 |

# 4 单层 Type I 预编码矩阵表达

## 4.1 空间波束选择

对单层传输，预编码矩阵退化为向量：

$$
\mathbf W=\mathbf w\in\mathbb C^{P_{\mathrm{CSI-RS}}\times1}.
$$

假设单极化 UPA 的空间位置为 $N_1\times N_2$，两个维度上的过采样因子分别为 $O_1,O_2$。第一维和第二维的 DFT 向量定义为：

$$
\mathbf g_{m_1}^{(N_1O_1)}
=
\begin{bmatrix}
1 &
e^{j2\pi\frac{m_1}{N_1O_1}} &
\cdots &
e^{j2\pi\frac{(N_1-1)m_1}{N_1O_1}}
\end{bmatrix}^T,
$$

$$
\mathbf u_{m_2}^{(N_2O_2)}
=
\begin{bmatrix}
1 &
e^{j2\pi\frac{m_2}{N_2O_2}} &
\cdots &
e^{j2\pi\frac{(N_2-1)m_2}{N_2O_2}}
\end{bmatrix}^T.
$$

二维 DFT 空间波束为：

$$
\mathbf w_{m_1,m_2}
=
\mathbf u_{m_2}^{(N_2O_2)}
\otimes
\mathbf g_{m_1}^{(N_1O_1)},
\qquad
\bar{\mathbf w}_{m_1,m_2}
=
\frac{1}{\sqrt{N_1N_2}}\mathbf w_{m_1,m_2}.
$$

若不考虑双极化，候选向量可以写成：

$$
\mathbf w_{m_1,m_2}
=
\bar{\mathbf w}_{m_1,m_2}.
$$

UE 可根据估计信道选择最大增益波束：

$$
(m_1^\star,m_2^\star)
=
\arg\max_{m_1,m_2}
\sum_{k\in\mathcal K}
\left\|
\hat{\mathbf H}_{\mathrm{CSI}}[m,k]\bar{\mathbf w}_{m_1,m_2}
\right\|_2^2.
$$

这里 $\mathcal K$ 可以是宽带 CSI 带宽内的子载波集合，也可以是某个子带。标准中的 $i_{1,1}$ 和 $i_{1,2}$ 用于指示 $m_1,m_2$。

## 4.2 双极化合成

对 Type I single-panel，双极化结构通常使用相同的 2D beam，并通过 co-phasing 因子合成两个极化端口。单层候选向量可写成：

$$
\mathbf w_{m_1,m_2,n}
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
\bar{\mathbf w}_{m_1,m_2}\\
\phi_n\bar{\mathbf w}_{m_1,m_2}
\end{bmatrix}.
$$

其中 co-phasing 因子可写成：

$$
\phi_n=e^{j\pi n/2},
\qquad
n\in\{0,1,2,3\}.
$$

这相当于用有限比特量化两个极化端口之间的相对相位。标准中的 $i_2$ 用于指示这类 co-phasing 信息；在多层和不同表格下，$i_2$ 的取值范围和解释会随 rank、端口数和码本模式变化。

UE 选择：

$$
(m_1^\star,m_2^\star,n^\star)
=
\arg\max_{m_1,m_2,n}
\sum_{k\in\mathcal K}
\left\|
\hat{\mathbf H}_{\mathrm{CSI}}[m,k]\mathbf w_{m_1,m_2,n}
\right\|_2^2.
$$

这个表达也说明了 Type I 码本和双极化信道建模的关系：如果真实信道中两个极化端口之间存在显著的相对相位，码本需要用 $i_2$ 或等价 co-phasing 索引捕获这种极化合成方向。

## 4.3 宽带和子带相位

宽带信道的主导空间方向通常随频率变化较慢，而相位和等效增益可能随频率变化更明显。因此 Type I 码本可采用宽带方向 + 子带补充的思想：

$$
\mathbf w_b
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
\bar{\mathbf w}_{m_1,m_2}\\
\phi_{n_b}\bar{\mathbf w}_{m_1,m_2}
\end{bmatrix}.
$$

其中：

- $(m_1,m_2)$ 由宽带 beam indicators 给出；
- $\phi_{n_b}$ 或对应 co-phasing matrix 可随子带 $b$ 变化；
- CQI 也可按宽带或子带反馈。

论文中也强调，Type I single-panel 的 2D beam 通常按宽带方式报告，而 co-phasing matrix 可以是频率相关的，需要按子带报告。这种结构降低了反馈开销：不需要每个子带重新反馈完整空间波束，只需反馈少量 co-phasing 或补充索引。

## 4.4 单层码本量化误差

为了评价 Type I 码本是否能准确表征理想预编码方向，可以定义单层预编码方向量化误差。设第 $k$ 个子载波的理想 SVD 预编码向量为：

$$
\mathbf w_{\mathrm{svd}}[k],
$$

Type I 码本选择得到的向量为：

$$
\mathbf w_{\mathrm{cb}}[k].
$$

由于预编码向量的公共相位不影响等效信道功率，比较时应使用子空间距离，而不是直接比较元素差值。常用量化误差为：

$$
\boxed{
\epsilon_{\mathrm{q}}[k]
=
1-
\frac{
\left|
\mathbf w_{\mathrm{svd}}^H[k]
\mathbf w_{\mathrm{cb}}[k]
\right|^2
}{
\|\mathbf w_{\mathrm{svd}}[k]\|_2^2
\|\mathbf w_{\mathrm{cb}}[k]\|_2^2
}
}
$$

若两者方向完全一致，则：

$$
\epsilon_{\mathrm{q}}[k]=0.
$$

若两者正交，则：

$$
\epsilon_{\mathrm{q}}[k]=1.
$$

对宽带 Type I PMI，可以统计整个 CSI 带宽上的平均误差：

$$
\bar{\epsilon}_{\mathrm{WB}}
=
\frac{1}{|\mathcal K_{\mathrm{WB}}|}
\sum_{k\in\mathcal K_{\mathrm{WB}}}
\epsilon_{\mathrm{q}}[k].
$$

对子带 PMI，可以对每个子带分别统计：

$$
\bar{\epsilon}_{b}
=
\frac{1}{|\mathcal K_b|}
\sum_{k\in\mathcal K_b}
\epsilon_{\mathrm{q}}[k].
$$

通常情况下，若信道主导方向随频率变化明显，子带 Type I 预编码的量化误差会小于宽带 Type I 预编码；若信道角度结构稳定、频率选择性主要体现在幅度相位上，则宽带 PMI 已经可以获得较好的预编码方向。

# 5 多层 Type I 预编码矩阵表达

## 5.1 多层传输与秩指示

多层传输时：

$$
\mathbf W
=
\begin{bmatrix}
\mathbf w_1 & \mathbf w_2 & \cdots & \mathbf w_\nu
\end{bmatrix}
\in
\mathbb C^{P_{\mathrm{CSI-RS}}\times \nu}.
$$

RI 表示 UE 推荐的层数 $\nu$。从信道矩阵角度看，$\nu$ 应和信道有效秩、奇异值分布、SINR 和接收机能力相关。

若奇异值满足：

$$
\sigma_1\gg\sigma_2\gg\cdots,
$$

则适合较低 rank；若多个奇异值都较大，则可支持更高 rank。

## 5.2 多波束/多极化层的矩阵构造

Type I 多层码本可以抽象理解为选择多个结构化列向量：

$$
\mathbf W
=
\left[
\mathbf w_{i_1},
\mathbf w_{i_2},
\dots,
\mathbf w_{i_\nu}
\right],
$$

其中每个 $\mathbf w_{i_\ell}$ 来自 Type I 单层候选集合或标准定义的规则化组合。

例如 rank-2 时，一个简单的双极化结构可以写成：

$$
\mathbf W
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
\bar{\mathbf w}_{m_{1,1},m_{2,1}} & \bar{\mathbf w}_{m_{1,2},m_{2,2}}\\
\phi_{n_1}\bar{\mathbf w}_{m_{1,1},m_{2,1}} &
\phi_{n_2}\bar{\mathbf w}_{m_{1,2},m_{2,2}}
\end{bmatrix}.
$$

若两个层使用相同空间波束但不同极化组合，则：

$$
\mathbf W
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
\bar{\mathbf w}_{m_1,m_2} & \bar{\mathbf w}_{m_1,m_2}\\
\phi_{n_1}\bar{\mathbf w}_{m_1,m_2} &
\phi_{n_2}\bar{\mathbf w}_{m_1,m_2}
\end{bmatrix}.
$$

但实际码本需要保证列之间具有合适的正交性或低相关性，否则层间干扰会增大。

## 5.3 归一化与功率约束

预编码矩阵通常需要满足总功率约束：

$$
\mathrm{tr}(\mathbf W^H\mathbf W)=\nu.
$$

如果每层单位范数：

$$
\|\mathbf w_\ell\|_2^2=1,
$$

且不同层近似正交：

$$
\mathbf w_i^H\mathbf w_j\approx0,\qquad i\ne j,
$$

则：

$$
\mathbf W^H\mathbf W\approx\mathbf I_\nu.
$$

这有助于降低发射侧层间相关性，也便于接收侧均衡。

# 6 宽带与子带预编码矩阵获取

Type I 码本并不是给出连续预编码优化问题的任意解，而是规定一组离散、结构化的候选预编码矩阵。实际获取过程可以概括为：

$$
\boxed{
\hat{\mathbf H}_{\mathrm{CSI}}[m,k]
\rightarrow
\text{UE 选择 RI/PMI/CQI}
\rightarrow
\text{gNB 根据 PMI 重建 } \mathbf W
}
$$

其中 UE 负责根据下行 CSI-RS 估计信道并选择推荐码本索引，gNB 负责根据标准码本规则和调度结果生成实际用于 PDSCH 的预编码矩阵。

## 6.1 宽带预编码矩阵获取

宽带预编码矩阵用于描述整个 CSI 带宽上相对稳定的主导空间方向。设 CSI 带宽内的子载波集合为 $\mathcal K_{\mathrm{WB}}$，UE 基于 CSI-RS 得到估计信道：

$$
\hat{\mathbf H}_{\mathrm{CSI}}[m,k],\qquad k\in\mathcal K_{\mathrm{WB}}.
$$

一种宽带选择方式是先构造宽带发射端等效协方差：

$$
\mathbf R_{\mathrm{WB}}
=
\sum_{k\in\mathcal K_{\mathrm{WB}}}
\hat{\mathbf H}_{\mathrm{CSI}}^{H}[m,k]\hat{\mathbf H}_{\mathrm{CSI}}[m,k].
$$

若不考虑码本量化，理想 $\nu$ 层宽带预编码方向可由 $\mathbf R_{\mathrm{WB}}$ 的前 $\nu$ 个主特征向量给出：

$$
\mathbf R_{\mathrm{WB}}\mathbf v_i
=
\lambda_i\mathbf v_i,\qquad i=1,\dots,\nu,
$$

$$
\mathbf W_{\mathrm{opt,WB}}
=
\begin{bmatrix}
\mathbf v_1 & \mathbf v_2 & \cdots & \mathbf v_\nu
\end{bmatrix}.
$$

Type I 码本不是反馈 $\mathbf W_{\mathrm{opt,WB}}$ 的所有复数元素，而是在宽带候选集合中选择最接近它的结构化矩阵：

$$
\boxed{
i_{\mathrm{WB}}^\star
=
\arg\max_{i_{\mathrm{WB}}}
\sum_{k\in\mathcal K_{\mathrm{WB}}}
\mathcal M
\left(
\hat{\mathbf H}_{\mathrm{CSI}}[m,k],
\mathbf W_{\mathrm{WB}}(i_{\mathrm{WB}})
\right)
}
$$

其中 $\mathbf W_{\mathrm{WB}}(i_{\mathrm{WB}})\in\mathcal W_{\mathrm{WB}}$ 由 Type I 码本规则生成。常用度量可以是等效信道增益：

$$
\mathcal M
\left(
\hat{\mathbf H}_{\mathrm{CSI}}[m,k],
\mathbf W
\right)
=
\left\|
\hat{\mathbf H}_{\mathrm{CSI}}[m,k]\mathbf W
\right\|_F^2,
$$

也可以是吞吐率近似：

$$
\mathcal M
\left(
\hat{\mathbf H}_{\mathrm{CSI}}[m,k],
\mathbf W
\right)
=
\log_2
\det
\left(
\mathbf I
+
\frac{\rho}{\nu}
\hat{\mathbf H}_{\mathrm{CSI}}[m,k]\mathbf W\mathbf W^H\hat{\mathbf H}_{\mathrm{CSI}}^{H}[m,k]
\right).
$$

选定宽带 PMI 后，宽带预编码矩阵为：

$$
\boxed{
\mathbf W_{\mathrm{WB}}
=
f_{\mathrm{WB}}(i_{\mathrm{WB}}^\star)
}
$$

其中 $f_{\mathrm{WB}}(\cdot)$ 是标准码本规则中的 PMI 到预编码矩阵的映射函数。对双极化面板，$i_{\mathrm{WB}}^\star$ 通常对应 2D DFT 波束方向、极化相位和 rank 相关结构。

## 6.2 子带预编码矩阵获取

宽带 PMI 假设主导空间方向在整个 CSI 带宽内变化较慢。但在频率选择性信道中，不同子带可能存在不同的相位、增益和局部最优方向。因此 NR 可以在宽带 PMI 的基础上引入子带 PMI。

设第 $b$ 个子带的子载波集合为 $\mathcal K_b$，其子带等效协方差为：

$$
\mathbf R_b
=
\sum_{k\in\mathcal K_b}
\hat{\mathbf H}_{\mathrm{CSI}}^{H}[m,k]\hat{\mathbf H}_{\mathrm{CSI}}[m,k].
$$

理想子带预编码方向为：

$$
\mathbf W_{\mathrm{opt},b}
=
\mathrm{eig}_{\nu}
\left(
\mathbf R_b
\right),
$$

其中 $\mathrm{eig}_{\nu}(\cdot)$ 表示取前 $\nu$ 个主特征向量组成的矩阵。

在 Type I 反馈中，子带预编码矩阵通常不是完全独立重新选择一个任意矩阵，而是在宽带 PMI 约束下选择补充索引：

$$
i_{\mathrm{SB},b}^\star
=
\arg\max_{i_{\mathrm{SB},b}}
\sum_{k\in\mathcal K_b}
\mathcal M
\left(
\hat{\mathbf H}_{\mathrm{CSI}}[m,k],
\mathbf W_b(i_{\mathrm{WB}}^\star,i_{\mathrm{SB},b})
\right).
$$

第 $b$ 个子带的预编码矩阵为：

$$
\boxed{
\mathbf W_b
=
f_{\mathrm{SB}}
\left(
i_{\mathrm{WB}}^\star,
i_{\mathrm{SB},b}^\star
\right)
}
$$

其中 $f_{\mathrm{SB}}(\cdot)$ 表示宽带 PMI 和子带 PMI 共同决定子带预编码矩阵。这样做的原因是：宽带 PMI 反馈主导空间方向，子带 PMI 只反馈频率选择性导致的局部修正，从而在性能和反馈开销之间折中。

## 6.3 宽带方向与子带补充的组合

对单层双极化 Type I 码本，一个直观的宽带 + 子带组合可写成：

$$
\mathbf w_b
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
\bar{\mathbf w}_{m_1^\star,m_2^\star}\\
\phi_{n_b^\star}
\bar{\mathbf w}_{m_1^\star,m_2^\star}
\end{bmatrix}.
$$

其中：

- $(m_1^\star,m_2^\star)$ 来自宽带 PMI 中的 beam indicators，表示宽带主导 2D DFT beam；
- $\phi_{n_b^\star}$ 来自第 $b$ 个子带的 co-phasing indicator，表示该子带上的双极化相位合成；
- $\mathbf w_b$ 是第 $b$ 个子带的单层预编码向量。

更一般地，对 $\nu$ 层传输：

$$
\mathbf W_b
=
\begin{bmatrix}
\mathbf w_{b,1} &
\mathbf w_{b,2} &
\cdots &
\mathbf w_{b,\nu}
\end{bmatrix},
$$

每一列由 Type I 码本的空间波束、极化相位和层间结构共同决定。

从矩阵分解角度看，这个过程可以概括为：

$$
\boxed{
\mathbf V_\nu[m,k]
\approx
\mathbf W_b
=
f_{\mathrm{codebook}}
\left(
i_{\mathrm{WB}}^\star,
i_{\mathrm{SB},b}^\star,
\nu^\star
\right)
}
$$

其中 $\mathbf V_\nu[m,k]$ 是理想信道右奇异子空间，$\mathbf W_b$ 是 Type I 码本可表达的量化近似。

## 6.4 gNB 侧如何使用反馈得到的预编码矩阵

UE 反馈 RI/PMI/CQI 后，gNB 并不是接收一个显式矩阵，而是接收码本索引。gNB 使用相同的 Type I 码本生成规则得到：

$$
\mathbf W_{\mathrm{WB}}
=
f_{\mathrm{WB}}(i_{\mathrm{WB}}^\star),
$$

或对每个子带得到：

$$
\mathbf W_b
=
f_{\mathrm{SB}}
\left(
i_{\mathrm{WB}}^\star,
i_{\mathrm{SB},b}^\star
\right).
$$

实际 PDSCH 传输时，gNB 还会结合调度、MCS、DM-RS 端口、功率分配、层映射和实现侧约束。若采用宽带预编码，则同一个 $\mathbf W_{\mathrm{WB}}$ 可用于整个调度带宽：

$$
\mathbf x[m,k]
=
\mathbf W_{\mathrm{WB}}\mathbf s[m,k],
\qquad
k\in\mathcal K_{\mathrm{sched}}.
$$

若采用子带预编码，则不同子带使用不同的 $\mathbf W_b$：

$$
\mathbf x[m,k]
=
\mathbf W_b\mathbf s[m,k],
\qquad
k\in\mathcal K_b.
$$

因此，Type I 码本下的预编码矩阵获取可以分成两个层次：

1. UE 侧：基于信道估计结果选择最优或近似最优的 RI/PMI/CQI；
2. gNB 侧：基于 PMI 和标准码本规则重建宽带或子带预编码矩阵，并将其用于下行 PDSCH。

## 6.5 一个端到端获取示例

为了把上述过程串起来，考虑一个简化的单用户、单层、双极化 UPA 场景：

$$
N_1=4,\quad N_2=2,\quad P_{\mathrm{CSI-RS}}=2N_1N_2=16.
$$

设 Type I 码本使用：

$$
O_1=4,\quad O_2=4,\quad
\phi_n=e^{j\pi n/2},\quad n\in\{0,1,2,3\}.
$$

**Step 1：UE 基于 CSI-RS 估计信道**

UE 得到：

$$
\hat{\mathbf H}_{\mathrm{CSI}}[m,k]\in\mathbb C^{N_r\times16},
\qquad
k\in\mathcal K_{\mathrm{CSI}}.
$$

**Step 2：UE 生成候选 Type I 预编码向量**

对每个候选 $(m_1,m_2,n)$，生成：

$$
\mathbf w_{m_1,m_2,n}
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
\bar{\mathbf w}_{m_1,m_2}\\
\phi_n\bar{\mathbf w}_{m_1,m_2}
\end{bmatrix}.
$$

**Step 3：UE 做宽带 PMI 搜索**

UE 计算：

$$
(m_1^\star,m_2^\star,n^\star)
=
\arg\max_{m_1,m_2,n}
\sum_{k\in\mathcal K_{\mathrm{WB}}}
\left\|
\hat{\mathbf H}_{\mathrm{CSI}}[m,k]\mathbf w_{m_1,m_2,n}
\right\|_2^2.
$$

这一步得到宽带主导方向：

$$
i_{\mathrm{WB}}^\star
\leftrightarrow
(m_1^\star,m_2^\star,n^\star).
$$

**Step 4：如果配置子带 PMI，UE 做子带补充搜索**

对每个子带 $b$，在固定或部分固定宽带方向的基础上搜索补充相位或补充索引：

$$
i_{\mathrm{SB},b}^\star
=
\arg\max_{i_{\mathrm{SB},b}}
\sum_{k\in\mathcal K_b}
\left\|
\hat{\mathbf H}_{\mathrm{CSI}}[m,k]
\mathbf w_b(i_{\mathrm{WB}}^\star,i_{\mathrm{SB},b})
\right\|_2^2.
$$

**Step 5：UE 反馈 CSI**

UE 反馈：

$$
\mathrm{CSI}
=
\{\mathrm{RI},\mathrm{PMI},\mathrm{CQI}\}.
$$

在单层例子中：

$$
\mathrm{RI}=1,
$$

PMI 包括宽带索引和可选的子带索引：

$$
\mathrm{PMI}
=
\left(
i_{\mathrm{WB}}^\star,
\{i_{\mathrm{SB},b}^\star\}_{b=1}^{N_{\mathrm{SB}}}
\right).
$$

**Step 6：gNB 重建预编码向量并用于 PDSCH**

若使用宽带预编码：

$$
\mathbf W[m,k]
=
f_{\mathrm{WB}}(i_{\mathrm{WB}}^\star),
\qquad
k\in\mathcal K_{\mathrm{sched}}.
$$

若使用子带预编码：

$$
\mathbf W[m,k]
=
f_{\mathrm{SB}}
\left(
i_{\mathrm{WB}}^\star,
i_{\mathrm{SB},b}^\star
\right),
\qquad
k\in\mathcal K_b.
$$

最终下行发送信号为：

$$
\mathbf x[m,k]
=
\mathbf W[m,k]\mathbf s[m,k].
$$

这个例子说明，Type I 码本反馈的核心不是让 UE 报告完整信道矩阵，而是让 UE 和 gNB 基于相同码本规则，对一个有限集合中的预编码方向达成一致。

# 7 Type I 码本与信道建模的关系

Type I 码本设计和 MIMO 信道建模密切相关。

**1. 与角度域信道建模的关系**

Type I 使用 DFT 波束，本质上利用了 ULA/UPA 阵列的角度域结构。若信道具有有限个主导 AoD：

$$
\mathbf H_{\mathrm{CSI}}[m,k]
=
\sum_{\ell=1}^{L_p}
\alpha_\ell[k]
\mathbf a_r(\Omega_{r,\ell})
\mathbf a_t^H(\Omega_{t,\ell}),
$$

则最优预编码方向往往接近某个发射端阵列导向矢量：

$$
\mathbf w_{\mathrm{opt}}
\approx
\mathbf a_t(\Omega_{t,\ell^\star}).
$$

DFT 码本就是用离散波束集合近似连续角度方向。

**2. 与双极化信道建模的关系**

双极化信道可写成：

$$
\mathbf H_{\mathrm{CSI}}[m,k]
=
\sum_{\ell}
\sqrt{P_\ell}
\mathbf A_r(\Omega_{r,\ell})
\mathbf P_\ell
\mathbf A_t^H(\Omega_{t,\ell})
e^{-j\frac{2\pi}{N}k\tau_\ell}.
$$

其中 $\mathbf P_\ell$ 是 $2\times2$ 极化耦合矩阵。Type I 码本中的 co-phasing 因子和双极化端口组合，正是为了用低开销方式适配这种极化耦合结构。

若 XPR 很高，交叉极化泄漏弱，两个极化之间相对独立，码本中 co-phasing 的影响可能较小。若 XPR 较低，交叉极化较强，合适的 co-phasing 对性能影响更明显。

**3. 与 Kronecker 相关模型的关系**

若信道采用空间-极化 Kronecker 相关模型：

$$
\mathbf R_t
=
\mathbf R_{t,\mathrm{space}}
\otimes
\mathbf R_{t,\mathrm{pol}},
$$

则 Type I 码本可以理解为对 $\mathbf R_t$ 主导特征方向的结构化量化：

- DFT 波束量化 $\mathbf R_{t,\mathrm{space}}$ 的主导空间方向；
- co-phasing 量化 $\mathbf R_{t,\mathrm{pol}}$ 或极化耦合导致的主导极化方向；
- rank/层结构量化多个主导特征方向。

**4. 与信道估计误差的关系**

UE 搜索码本依赖 $\hat{\mathbf H}_{\mathrm{CSI}}$。若信道估计误差为：

$$
\hat{\mathbf H}_{\mathrm{CSI}}[m,k]
=
\mathbf H_{\mathrm{CSI}}[m,k]+\mathbf E[m,k],
$$

则 PMI 选择可能发生错误：

$$
\arg\max_i
\mathcal M(\hat{\mathbf H}_{\mathrm{CSI}},\mathbf W_i)
\ne
\arg\max_i
\mathcal M(\mathbf H_{\mathrm{CSI}},\mathbf W_i).
$$

因此，CSI-RS 密度、信道插值、LMMSE/Wiener 估计、双选择性衰落跟踪都会影响 Type I 码本反馈质量。

**5. 信道参数对 Type I 码本性能的影响**

Type I 码本的性能和信道建模参数直接相关。可以从以下几个维度理解。

| 信道参数 | 对 Type I 码本的影响 |
|---|---|
| AoD 角度扩展 | 角度扩展越小，主导方向越接近单个 DFT 波束；角度扩展越大，单波束 Type I 量化误差越大 |
| 路径数 $L_p$ | 路径越少，信道越稀疏，DFT 波束越有效；路径越多，主导子空间可能需要多波束组合 |
| XPR | XPR 越低，交叉极化越强，co-phasing 选择越重要 |
| 极化相关 $\rho_{\mathrm{pol}}$ | 极化相关越高，两个极化端口越不独立，可支持 rank 可能下降 |
| RMS 时延扩展 | 时延扩展越大，频率选择性越强，子带 PMI 比宽带 PMI 更有价值 |
| 多普勒 | 多普勒越大，PMI 老化越明显，需要更密集 CSI-RS 或更短反馈周期 |
| 信道估计 NMSE | NMSE 越高，PMI 误选概率越高，码本反馈增益下降 |

以前文的角度域和双极化模型为例，若 AoD 集中，即 $\Omega_{t,\ell}$ 分布很窄，则 $\mathbf A_t(\Omega_{t,\ell})$ 的主导方向接近某个 DFT beam，Type I 码本容易选到接近最优的 $(m_1,m_2)$。若 AoD 分布很宽，则一个 Type I 单波束方向难以同时匹配所有路径，量化误差增大。

若 XPR 较低，$\mathbf P_\ell$ 的非对角项较强：

$$
\left|\alpha_{\ell,VH}\right|^2,\,
\left|\alpha_{\ell,HV}\right|^2
\not\ll
\left|\alpha_{\ell,VV}\right|^2,\,
\left|\alpha_{\ell,HH}\right|^2,
$$

则两个极化端口之间的合成相位会显著影响等效信道：

$$
\mathbf H_{\mathrm{CSI}}[m,k]\mathbf w_{m_1,m_2,n}.
$$

这时增加 co-phasing 候选集合或使用更精细的子带 co-phasing，可能降低量化误差。反之，若 XPR 很高、交叉极化很弱，co-phasing 的收益可能小于空间 beam selection 的收益。

对 Kronecker 空间-极化相关模型，如果 $\mathbf R_{t,\mathrm{space}}$ 的最大特征值明显大于其他特征值，Type I 的主导 DFT 波束量化通常有效；如果 $\mathbf R_{t,\mathrm{pol}}$ 的两个极化高度相关，则双极化端口提供的独立空间维度有限，高 rank 反馈不一定带来明显增益。

因此，在仿真 Type I 码本时，建议至少扫描：

- AoD/AoA 角度扩展；
- XPR；
- 极化相关系数；
- RMS 时延扩展；
- 信道估计 NMSE；
- 宽带 PMI 与子带 PMI 的反馈粒度。

# 8 小结

NR Type I 码本是面向低开销 CSI 反馈的结构化预编码码本。它的核心不是反馈完整信道，而是用有限索引表示一个适合下行传输的预编码矩阵。按照 NR 术语，这个码本矩阵定义在 CSI-RS antenna ports 维度上，行数对应 $P_{\mathrm{CSI-RS}}$，列数对应传输 layers/DM-RS antenna ports。

从数学上看，理想预编码来自信道矩阵的右奇异向量：

$$
\mathbf H
=
\mathbf U\mathbf \Sigma\mathbf V^H,
\qquad
\mathbf W_{\mathrm{opt}}
=
\mathbf V_\nu.
$$

Type I 码本用标准化结构近似 $\mathbf V_\nu$：

$$
\mathbf W_{\mathrm{opt}}
\approx
\mathbf W(\mathrm{PMI}).
$$

其主要设计思想包括：

- 用 2D DFT 波束表示空间方向；
- 用 co-phasing 因子表示双极化端口之间的相位合成；
- 用 RI 表示传输层数；
- 用宽带/子带 PMI 兼顾反馈开销和频率选择性；
- 用 CQI 反映该预编码假设下的链路质量。

在 NR FDD 系统中，Type I 码本是下行信道获取和预编码反馈的重要机制。在 TDD 系统中，当上行覆盖受限、SRS 不可靠或校准误差较大时，Type I 码本也可以作为下行 CSI 获取的补充手段。

与信道建模的关系可以概括为：

- MIMO-OFDM 信道建模给出 $\mathbf H_{\mathrm{CSI}}[m,k]$；
- CSI-RS 信道估计给出 $\hat{\mathbf H}_{\mathrm{CSI}}[m,k]$；
- SVD 或互信息准则给出理想预编码方向；
- Type I 码本用 DFT beam、co-phasing 和层结构对理想预编码方向做有限比特量化。

因此，理解 Type I 码本需要同时理解 MIMO-OFDM 信道矩阵、双极化阵列结构、信道估计误差和有限反馈预编码之间的关系。

# 9 参考文献

[1] 3GPP TS 38.214, "NR; Physical layer procedures for data," section 5.2.2, CSI reporting and Type I/Type II codebook related procedures.

[2] 3GPP TS 38.211, "NR; Physical channels and modulation," CSI-RS, antenna port, OFDM, and downlink physical channel definitions.

[3] 3GPP TS 38.212, "NR; Multiplexing and channel coding," CSI report payload and uplink control information encoding.

[4] 3GPP TS 38.331, "NR; Radio Resource Control (RRC) protocol specification," CSI-MeasConfig, CSI-ResourceConfig, CSI-ReportConfig, and codebook related RRC configuration.

[5] 3GPP TR 38.901, "Study on channel model for frequencies from 0.5 to 100 GHz," MIMO, polarization, antenna array, cluster delay line, and spatial channel modeling.

[6] D. J. Love, R. W. Heath Jr., V. K. N. Lau, D. Gesbert, B. D. Rao, and M. Andrews, "An overview of limited feedback in wireless communication systems," IEEE Journal on Selected Areas in Communications, vol. 26, no. 8, pp. 1341-1365, Oct. 2008.

[7] D. J. Love and R. W. Heath Jr., "Limited feedback unitary precoding for spatial multiplexing systems," IEEE Transactions on Information Theory, vol. 51, no. 8, pp. 2967-2976, Aug. 2005.

[8] A. Paulraj, R. Nabar, and D. Gore, Introduction to Space-Time Wireless Communications, Cambridge University Press, 2003.

[9] D. Tse and P. Viswanath, Fundamentals of Wireless Communication, Cambridge University Press, 2005.

[10] T. L. Marzetta, E. G. Larsson, H. Yang, and H. Q. Ngo, Fundamentals of Massive MIMO, Cambridge University Press, 2016.

[11] E. Björnson, J. Hoydis, and L. Sanguinetti, Massive MIMO Networks: Spectral, Energy, and Hardware Efficiency, Foundations and Trends in Signal Processing, 2017.

[12] Z. Qin and H. Yin, "A Review of Codebooks for CSI Feedback in 5G New Radio and Beyond," arXiv:2302.09222v2, Jun. 2023.

[13] ShareTechnote, "5G/NR - CSI RS Codebook," https://www.sharetechnote.com/html/5G/5G_CSI_RS_Codebook.html.
