# NR Type I 码本设计原理

> 本文基于 MIMO-OFDM 信道建模、双极化收发天线建模以及 MIMO-OFDM 信道估计的背景，介绍 NR Type I 码本的设计动机、使用场景和数学原理。重点说明 Type I 码本如何用有限反馈比特表征下行预编码矩阵，并从信道矩阵分解过渡到码本预编码矩阵的构造。

## 目录

- [1 背景与动机](#1-背景与动机)
  - [1.1 从 MIMO-OFDM 信道到预编码](#11-从-mimo-ofdm-信道到预编码)
  - [1.2 为什么需要码本](#12-为什么需要码本)
  - [1.3 Type I 码本的定位](#13-type-i-码本的定位)
- [2 使用场景](#2-使用场景)
  - [2.1 NR FDD 系统的下行信道获取](#21-nr-fdd-系统的下行信道获取)
  - [2.2 TDD 系统中上行覆盖受限时的信道获取](#22-tdd-系统中上行覆盖受限时的信道获取)
  - [2.3 与非码本反馈和 Type II 码本的关系](#23-与非码本反馈和-type-ii-码本的关系)
- [3 Notation](#3-notation)
- [4 MIMO-OFDM 信道、信道估计与预编码目标](#4-mimo-ofdm-信道信道估计与预编码目标)
  - [4.1 频域 MIMO-OFDM 系统模型](#41-频域-mimo-ofdm-系统模型)
  - [4.2 信道估计与 CSI 获取](#42-信道估计与-csi-获取)
  - [4.3 最优预编码的基本形式](#43-最优预编码的基本形式)
- [5 从信道矩阵分解到预编码矩阵](#5-从信道矩阵分解到预编码矩阵)
  - [5.1 单子载波窄带信道的 SVD 预编码](#51-单子载波窄带信道的-svd-预编码)
  - [5.2 宽带 OFDM 中的子带预编码](#52-宽带-ofdm-中的子带预编码)
  - [5.3 预编码矩阵量化问题](#53-预编码矩阵量化问题)
- [6 NR Type I 码本设计原理](#6-nr-type-i-码本设计原理)
  - [6.1 Type I 码本的核心思想](#61-type-i-码本的核心思想)
  - [6.2 2D DFT 波束和端口索引](#62-2d-dft-波束和端口索引)
  - [6.3 双极化端口结构](#63-双极化端口结构)
  - [6.4 PMI 如何表征预编码矩阵](#64-pmi-如何表征预编码矩阵)
- [7 单层 Type I 预编码矩阵表达](#7-单层-type-i-预编码矩阵表达)
  - [7.1 空间波束选择](#71-空间波束选择)
  - [7.2 双极化合成](#72-双极化合成)
  - [7.3 宽带和子带相位](#73-宽带和子带相位)
- [8 多层 Type I 预编码矩阵表达](#8-多层-type-i-预编码矩阵表达)
  - [8.1 多层传输与秩指示](#81-多层传输与秩指示)
  - [8.2 多波束/多极化层的矩阵构造](#82-多波束多极化层的矩阵构造)
  - [8.3 归一化与功率约束](#83-归一化与功率约束)
- [9 预编码矩阵的获取方法](#9-预编码矩阵的获取方法)
  - [9.1 Type I 码本的 SU-MIMO 定位](#91-type-i-码本的-su-mimo-定位)
  - [9.2 宽带预编码矩阵获取](#92-宽带预编码矩阵获取)
  - [9.3 子带预编码矩阵获取](#93-子带预编码矩阵获取)
  - [9.4 宽带方向与子带补充的组合](#94-宽带方向与子带补充的组合)
  - [9.5 gNB 侧如何使用反馈得到的预编码矩阵](#95-gnb-侧如何使用反馈得到的预编码矩阵)
- [10 Type I 码本反馈流程](#10-type-i-码本反馈流程)
- [11 Type I 码本与信道建模的关系](#11-type-i-码本与信道建模的关系)
- [12 小结](#12-小结)
- [13 参考文献](#13-参考文献)

# 1 背景与动机

## 1.1 从 MIMO-OFDM 信道到预编码

在 MIMO-OFDM 系统中，第 $m$ 个 OFDM 符号、第 $k$ 个子载波上的下行频域模型可写为：

$$
\mathbf y[m,k]
=
\mathbf H[m,k]\mathbf x[m,k]+\mathbf w[m,k],
$$

其中：

$$
\mathbf H[m,k]\in\mathbb C^{N_r\times N_t}.
$$

$N_t$ 是 gNB 发射端口数，$N_r$ 是 UE 接收端口数。若发射端进行多层传输，发送符号向量为：

$$
\mathbf s[m,k]\in\mathbb C^{L\times 1},
$$

其中 $L$ 是传输层数。预编码后：

$$
\boxed{
\mathbf x[m,k]
=
\mathbf W[m,k]\mathbf s[m,k]
}
$$

其中：

$$
\mathbf W[m,k]\in\mathbb C^{N_t\times L}
$$

是下行预编码矩阵。接收端看到的等效信道为：

$$
\mathbf y[m,k]
=
\mathbf H[m,k]\mathbf W[m,k]\mathbf s[m,k]
+\mathbf w[m,k].
$$

预编码的目的，是让 $\mathbf H[m,k]\mathbf W[m,k]$ 具有更好的信号增益、更低的层间干扰、更好的条件数，最终提升吞吐率、覆盖和可靠性。

## 1.2 为什么需要码本

如果 gNB 完全知道下行信道 $\mathbf H[m,k]$，可以直接根据 SVD、MMSE、ZF、SLNR 或其他准则计算预编码矩阵。但在 NR FDD 下行中，gNB 不能直接利用上行信道推断下行信道，因为上下行频率不同，信道不满足严格互易性。

因此 FDD 下行 CSI 获取通常采用如下流程：

1. gNB 发送 CSI-RS；
2. UE 基于 CSI-RS 估计下行信道；
3. UE 从标准定义的码本集合中选择一个或多个预编码候选；
4. UE 向 gNB 反馈 CSI，包括 RI、PMI、CQI；
5. gNB 根据反馈的 PMI 选择预编码矩阵。

码本的核心作用是：

$$
\boxed{
\text{用有限反馈比特描述一个接近最优的预编码矩阵}
}
$$

如果没有码本，UE 需要反馈完整的复数信道矩阵或完整预编码矩阵，其开销很高。例如一个 $N_t=32$、$L=4$ 的预编码矩阵有 $128$ 个复数元素。如果每个复数用有限精度量化，反馈开销远高于 NR 实际可承受范围。

码本把连续空间中的预编码矩阵量化为有限集合：

$$
\mathcal W
=
\left\{
\mathbf W_1,\mathbf W_2,\dots,\mathbf W_{|\mathcal W|}
\right\}.
$$

UE 只需反馈索引：

$$
i^\star
=
\arg\max_i
\mathcal M(\mathbf H,\mathbf W_i),
$$

其中 $\mathcal M(\cdot)$ 是某种性能度量，例如互信息、SINR、接收功率或与最优右奇异向量的夹角。

## 1.3 Type I 码本的定位

NR 下行码本反馈中，Type I 码本主要面向较低反馈开销的 beam-based CSI 反馈。它的设计目标不是完整恢复下行信道，而是让 UE 用少量比特报告一个适合下行传输的预编码方向。

Type I 码本具有以下特点：

- 基于 2D DFT 波束，适配 UPA 和双极化天线面板；
- 反馈开销相对较低；
- 主要表征一个或少数主导空间方向；
- 适合常规 SU-MIMO 或低复杂度 CSI 反馈；
- 对频率选择性的刻画较粗，通常以宽带/子带 PMI 方式表达。

相比 Type II 码本，Type I 更像是“选择一个标准化波束方向及极化合成方式”；Type II 则更强调用多个波束线性组合更精细地逼近信道特征向量。

因此，Type I 码本首先应理解为 SU-MIMO 下行链路的预编码矩阵量化机制。它可以作为 MU-MIMO 调度和配对的输入信息之一，但 Type I 码本本身不直接解决多用户之间的零干扰约束。多用户场景下，gNB 还需要结合多个 UE 的反馈、调度策略和多用户预编码算法进一步处理。

# 2 使用场景

## 2.1 NR FDD 系统的下行信道获取

FDD 系统上下行使用不同频点。由于传播信道的频率响应随频率变化，尤其在大带宽和高频段下，上行信道不能直接代表下行信道。因此 gNB 需要依赖 UE 反馈下行 CSI。

Type I 码本在 FDD 下行中的典型流程为：

$$
\text{CSI-RS}
\rightarrow
\hat{\mathbf H}_{\mathrm{DL}}
\rightarrow
\text{码本搜索}
\rightarrow
\text{RI/PMI/CQI}
\rightarrow
\text{gNB 下行预编码}.
$$

其中：

- RI 表示推荐传输层数；
- PMI 表示推荐预编码矩阵索引；
- CQI 表示在该预编码和秩假设下的可支持调制编码等级。

Type I 码本的作用是降低 FDD 下行 CSI 反馈开销。UE 不反馈完整信道矩阵，而是反馈对下行传输最有用的方向信息。

## 2.2 TDD 系统中上行覆盖受限时的信道获取

TDD 系统理论上可以利用上下行信道互易性。gNB 可通过上行 SRS 估计上行信道，再经射频校准得到下行预编码方向。

但实际系统中，TDD 仍可能需要下行 CSI 反馈或码本辅助：

- UE 发射功率有限，上行 SRS 覆盖不足；
- 小区边缘 UE 上行链路质量差，gNB 无法可靠估计高维上行信道；
- 大规模阵列中，全端口 SRS 开销较高；
- 多 TRP、多面板和高频场景中，上下行校准误差不可忽略；
- gNB 可能需要 UE 侧基于下行 CSI-RS 的波束/PMI 辅助。

因此，在 TDD 上行覆盖受限时，可以使用类似 FDD 的下行 CSI-RS + PMI 反馈方式，辅助 gNB 获取下行预编码方向。

这种场景下，Type I 码本提供的是一种低开销的下行方向反馈机制。它不能完全替代高精度信道估计，但能在上行不可用或不可靠时提供可用的预编码方向。

## 2.3 与非码本反馈和 Type II 码本的关系

CSI 反馈大致可以分为三类：

| 方式 | 反馈内容 | 优点 | 缺点 |
|---|---|---|---|
| Type I 码本 | 主导波束、极化合成、秩和 CQI | 开销低、标准化、实现简单 | 精度有限 |
| Type II 码本 | 多波束线性组合系数 | 精度高，适合 MU-MIMO | 开销更高，复杂度更高 |
| 非码本/显式反馈 | 信道或预编码矩阵量化值 | 最灵活 | 开销很大，标准化和实现复杂 |

Type I 码本适合常规 FDD CSI 反馈和低开销波束选择。若需要更高空间分辨率、更强 MU-MIMO 分离能力，通常需要 Type II 或增强型码本。

# 3 Notation

| 符号 | 含义 |
|---|---|
| $m$ | OFDM 符号索引 |
| $k$ | OFDM 子载波索引 |
| $b$ | 子带索引 |
| $N$ | OFDM FFT 点数或子载波数 |
| $N_t$ | gNB 发射天线端口数 |
| $N_r$ | UE 接收天线端口数 |
| $N_1,N_2$ | 发射端二维天线面板在两个空间维度上的端口数 |
| $O_1,O_2$ | 2D DFT 码本在两个维度上的过采样因子 |
| $N_t^{\mathrm{pos}}$ | 发射端物理天线位置数 |
| $N_t^{\mathrm{pol}}$ | 发射端每个位置的极化端口数，双极化时通常为 2 |
| $L$ | 传输层数，也称 rank |
| $\mathbf H[m,k]$ | 第 $m$ 个 OFDM 符号、第 $k$ 个子载波上的 MIMO CFR |
| $\hat{\mathbf H}[m,k]$ | UE 基于 CSI-RS 得到的下行信道估计 |
| $\mathbf x[m,k]$ | 预编码后的发射向量，$\mathbb C^{N_t\times1}$ |
| $\mathbf s[m,k]$ | 层符号向量，$\mathbb C^{L\times1}$ |
| $\mathbf y[m,k]$ | 接收向量，$\mathbb C^{N_r\times1}$ |
| $\mathbf W[m,k]$ | 预编码矩阵，$\mathbb C^{N_t\times L}$ |
| $\mathcal W$ | 码本集合 |
| $\mathcal W_{\mathrm{WB}}$ | 宽带 Type I 码本候选集合 |
| $\mathcal W_{\mathrm{SB}}(i_{\mathrm{WB}})$ | 给定宽带 PMI 后的子带补充候选集合 |
| $\mathbf W_i$ | 码本中的第 $i$ 个预编码候选 |
| $\mathbf W_{\mathrm{WB}}$ | 宽带预编码矩阵 |
| $\mathbf W_b$ | 第 $b$ 个子带的预编码矩阵 |
| $i_{\mathrm{WB}}$ | 宽带 PMI |
| $i_{\mathrm{SB},b}$ | 第 $b$ 个子带的子带 PMI |
| RI | Rank Indicator，秩指示 |
| PMI | Precoding Matrix Indicator，预编码矩阵指示 |
| CQI | Channel Quality Indicator，信道质量指示 |
| $\mathbf v_{l,m}$ | 2D DFT 空间波束向量 |
| $l,m$ | 2D DFT 波束在两个空间维度上的索引 |
| $\phi$ | 双极化端口之间的相对相位 |
| $\mathbf U,\mathbf \Sigma,\mathbf V$ | 信道矩阵 SVD 分解中的左奇异矩阵、奇异值矩阵、右奇异矩阵 |
| $\sigma_i$ | 第 $i$ 个奇异值 |
| $\mathbf V_L$ | 最优 $L$ 层预编码方向，由前 $L$ 个右奇异向量组成 |
| $\|\cdot\|_F$ | Frobenius 范数 |
| $(\cdot)^H$ | 共轭转置 |
| $\otimes$ | Kronecker 积 |

# 4 MIMO-OFDM 信道、信道估计与预编码目标

## 4.1 频域 MIMO-OFDM 系统模型

考虑下行 MIMO-OFDM。gNB 有 $N_t$ 个发射端口，UE 有 $N_r$ 个接收端口。若传输 $L$ 层数据，则：

$$
\mathbf x[m,k]
=
\mathbf W[m,k]\mathbf s[m,k],
$$

$$
\mathbf y[m,k]
=
\mathbf H[m,k]\mathbf W[m,k]\mathbf s[m,k]
+\mathbf w[m,k].
$$

其中通常要求：

$$
L\le \min(N_t,N_r).
$$

预编码矩阵需要满足功率归一化约束。常见写法为：

$$
\mathrm{tr}\left(\mathbf W^H\mathbf W\right)=L
$$

或对每一列单位范数：

$$
\|\mathbf w_i\|_2^2=1.
$$

具体归一化方式取决于标准定义和发射功率分配方式。

## 4.2 信道估计与 CSI 获取

在下行 CSI 获取中，UE 首先基于 CSI-RS 估计下行信道。对某个 CSI-RS 资源，简化模型为：

$$
\mathbf Y_p
=
\mathbf H_p\mathbf X_p+\mathbf W_p.
$$

若 CSI-RS 端口正交，UE 可得到各端口到接收天线的信道估计：

$$
\hat{\mathbf H}[m,k].
$$

在 MIMO-OFDM 信道估计中，UE 可使用 LS、MMSE、LMMSE、DFT 插值、Wiener 插值等方法获得完整时频网格或子带上的信道估计。码本反馈并不要求 UE 把完整信道反馈给 gNB，而是基于 $\hat{\mathbf H}$ 选择一个预编码矩阵索引。

因此信道估计和码本反馈之间的关系是：

$$
\boxed{
\hat{\mathbf H}
\rightarrow
\text{码本搜索}
\rightarrow
\text{PMI}
}
$$

信道估计质量会直接影响 PMI 的可靠性。若 $\hat{\mathbf H}$ 噪声较大或插值误差较大，UE 选择的预编码方向可能偏离真实最优方向。

## 4.3 最优预编码的基本形式

如果不考虑码本限制，单用户 MIMO 的理想预编码通常与信道矩阵的右奇异向量有关。

对某个子载波或子带等效信道：

$$
\mathbf H
=
\mathbf U\mathbf \Sigma\mathbf V^H.
$$

若传输 $L$ 层，理想预编码矩阵可取：

$$
\boxed{
\mathbf W_{\mathrm{opt}}
=
\mathbf V_L
}
$$

其中 $\mathbf V_L$ 由 $\mathbf V$ 的前 $L$ 列构成：

$$
\mathbf V_L
=
\begin{bmatrix}
\mathbf v_1 & \mathbf v_2 & \cdots & \mathbf v_L
\end{bmatrix}.
$$

此时等效信道为：

$$
\mathbf H\mathbf W_{\mathrm{opt}}
=
\mathbf U\mathbf \Sigma\mathbf V^H\mathbf V_L
=
\mathbf U_L\mathbf \Sigma_L.
$$

理想 SVD 预编码能把 MIMO 信道分解成若干近似并行的空间子信道。

码本设计的核心问题就是：如何用有限集合 $\mathcal W$ 逼近这些最优右奇异向量。

# 5 从信道矩阵分解到预编码矩阵

## 5.1 单子载波窄带信道的 SVD 预编码

先考虑单子载波窄带 MIMO 信道：

$$
\mathbf y=\mathbf H\mathbf W\mathbf s+\mathbf w.
$$

对信道做 SVD：

$$
\mathbf H
=
\mathbf U
\begin{bmatrix}
\sigma_1 & & \\
& \sigma_2 & \\
& & \ddots
\end{bmatrix}
\mathbf V^H,
\qquad
\sigma_1\ge\sigma_2\ge\cdots.
$$

若传输 1 层，最优方向是最大奇异值对应的右奇异向量：

$$
\mathbf w_{\mathrm{opt}}
=
\mathbf v_1.
$$

若传输 $L$ 层：

$$
\mathbf W_{\mathrm{opt}}
=
\left[
\mathbf v_1,\dots,\mathbf v_L
\right].
$$

这说明预编码矩阵的列空间应尽量对齐信道的主导右奇异子空间。

## 5.2 宽带 OFDM 中的子带预编码

在 OFDM 系统中，信道随子载波变化：

$$
\mathbf H[k]
=
\sum_{n=0}^{L_h-1}
\mathbf H[n]e^{-j\frac{2\pi}{N}kn}.
$$

如果每个子载波都独立反馈 PMI，开销会很高。因此 NR 通常按宽带或子带进行反馈：

- 宽带 PMI：反映整个 CSI 带宽上的平均主导方向；
- 子带 PMI：在宽带方向基础上补充频率选择性变化；
- CQI 可宽带或子带反馈；
- RI 通常是宽带或半静态变化较慢的量。

对第 $b$ 个子带，可定义等效信道协方差：

$$
\mathbf R_b
=
\sum_{k\in\mathcal K_b}
\hat{\mathbf H}^H[k]\hat{\mathbf H}[k].
$$

理想子带预编码方向可由 $\mathbf R_b$ 的主特征向量给出：

$$
\mathbf R_b\mathbf v_i
=
\lambda_i\mathbf v_i.
$$

Type I 码本相当于用预定义 DFT 波束集合去量化这些主特征方向。

## 5.3 预编码矩阵量化问题

设理想预编码矩阵为 $\mathbf W_{\mathrm{opt}}$，码本集合为 $\mathcal W$。UE 需要选择：

$$
\boxed{
\mathbf W^\star
=
\arg\max_{\mathbf W_i\in\mathcal W}
\mathcal M(\hat{\mathbf H},\mathbf W_i)
}
$$

常见度量包括：

1. **最大等效信道增益**

   单层时：

   $$
   \mathbf w^\star
   =
   \arg\max_{\mathbf w_i\in\mathcal W}
   \|\hat{\mathbf H}\mathbf w_i\|_2^2.
   $$

2. **最小子空间距离**

   $$
   \mathbf W^\star
   =
   \arg\min_{\mathbf W_i\in\mathcal W}
   \left\|
   \mathbf W_{\mathrm{opt}}\mathbf W_{\mathrm{opt}}^H
   -
   \mathbf W_i\mathbf W_i^H
   \right\|_F^2.
   $$

3. **最大互信息或吞吐率指标**

   $$
   \mathbf W^\star
   =
   \arg\max_{\mathbf W_i\in\mathcal W}
   \log_2
   \det
   \left(
   \mathbf I
   +
   \frac{\rho}{L}
   \hat{\mathbf H}\mathbf W_i\mathbf W_i^H\hat{\mathbf H}^H
   \right).
   $$

实际 UE 实现不一定显式做 SVD 或穷举全部精确吞吐率，但这些表达揭示了码本选择的数学目标。

# 6 NR Type I 码本设计原理

## 6.1 Type I 码本的核心思想

Type I 码本的核心思想是：

$$
\boxed{
\text{用 DFT 空间波束 + 极化合成 + 层间结构来表示预编码矩阵}
}
$$

这和 MIMO 信道的角度域结构直接相关。对 ULA/UPA 阵列，远场平面波的阵列响应天然接近 DFT 波束。若信道主导路径集中在某个 AoD 区域，则最优预编码方向往往接近某个 DFT 波束或少数相邻 DFT 波束的组合。

Type I 码本不是直接量化 $\mathbf H$，而是量化预编码方向：

$$
\hat{\mathbf H}
\rightarrow
\mathbf W_{\mathrm{opt}}
\approx
\mathbf W(i_1,i_2,\dots)
\rightarrow
\mathrm{PMI}.
$$

其中 PMI 的不同部分通常表示：

- 选择哪个 2D DFT 空间波束；
- 选择哪个极化相位；
- 选择哪个层数和层间组合；
- 对频率选择性信道，选择宽带或子带相关的补充索引。

## 6.2 2D DFT 波束和端口索引

NR 大规模 MIMO 常用二维天线面板。假设单极化空间位置为 $N_1\times N_2$，其中：

- $N_1$：第一维端口数；
- $N_2$：第二维端口数；
- $N_t^{\mathrm{pos}}=N_1N_2$。

对第一维定义 DFT 向量：

$$
\mathbf u_{l}^{(N_1O_1)}
=
\begin{bmatrix}
1 &
e^{j2\pi\frac{l}{N_1O_1}} &
\cdots &
e^{j2\pi\frac{(N_1-1)l}{N_1O_1}}
\end{bmatrix}^T.
$$

对第二维：

$$
\mathbf u_{m}^{(N_2O_2)}
=
\begin{bmatrix}
1 &
e^{j2\pi\frac{m}{N_2O_2}} &
\cdots &
e^{j2\pi\frac{(N_2-1)m}{N_2O_2}}
\end{bmatrix}^T.
$$

二维 DFT 空间波束可写成 Kronecker 积：

$$
\boxed{
\mathbf v_{l,m}
=
\mathbf u_l^{(N_1O_1)}
\otimes
\mathbf u_m^{(N_2O_2)}
}
$$

为了满足单位功率约束，通常归一化为：

$$
\bar{\mathbf v}_{l,m}
=
\frac{1}{\sqrt{N_1N_2}}
\mathbf v_{l,m}.
$$

这里 $O_1,O_2$ 是过采样因子。过采样的作用是让候选波束方向更密，从而降低真实 AoD 和离散 DFT 网格不对齐造成的量化误差。

## 6.3 双极化端口结构

NR 基站常使用双极化天线面板。如果每个空间位置有两个极化端口，则总端口数为：

$$
N_t=2N_1N_2.
$$

一个自然的端口排列方式是：

$$
\mathbf w
=
\begin{bmatrix}
\mathbf w^{(0)}\\
\mathbf w^{(1)}
\end{bmatrix},
$$

其中：

- $\mathbf w^{(0)}$ 是第一组极化端口上的权重；
- $\mathbf w^{(1)}$ 是第二组极化端口上的权重；
- 两者维度均为 $N_1N_2$。

若两个极化使用相同空间波束，只差一个相对相位，则单层预编码向量可写成：

$$
\boxed{
\mathbf w_{l,m,\phi}
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
\bar{\mathbf v}_{l,m}\\
e^{j\phi}\bar{\mathbf v}_{l,m}
\end{bmatrix}
}
$$

其中 $\phi$ 是两个极化端口之间的相对相位。

这个结构清楚反映了 Type I 码本对双极化阵列的建模方式：

- 用 $\bar{\mathbf v}_{l,m}$ 表示空间方向；
- 用 $e^{j\phi}$ 表示两个极化端口之间的相位合成；
- 用有限索引 $(l,m,\phi)$ 表示预编码向量。

## 6.4 PMI 如何表征预编码矩阵

Type I 码本中的 PMI 并不是直接反馈矩阵元素，而是反馈生成矩阵所需的离散索引。

抽象地，可写为：

$$
\boxed{
\mathrm{PMI}
\longleftrightarrow
(l,m,\phi,\text{rank},\text{subband index},\dots)
}
$$

gNB 和 UE 都知道标准定义的码本生成规则。因此 UE 只需反馈索引，gNB 即可重建预编码矩阵：

$$
\boxed{
\mathbf W
=
f_{\mathrm{codebook}}(\mathrm{PMI})
}
$$

对单层传输，PMI 主要表征一个波束方向和极化相位；对多层传输，PMI 还需要表征多个列向量的组合方式。

# 7 单层 Type I 预编码矩阵表达

## 7.1 空间波束选择

对单层传输，预编码矩阵退化为向量：

$$
\mathbf W=\mathbf w\in\mathbb C^{N_t\times1}.
$$

若不考虑双极化，候选向量可以写成：

$$
\mathbf w_{l,m}
=
\bar{\mathbf v}_{l,m}.
$$

UE 可根据估计信道选择最大增益波束：

$$
(l^\star,m^\star)
=
\arg\max_{l,m}
\sum_{k\in\mathcal K}
\left\|
\hat{\mathbf H}[k]\bar{\mathbf v}_{l,m}
\right\|_2^2.
$$

这里 $\mathcal K$ 可以是宽带 CSI 带宽内的子载波集合，也可以是某个子带。

## 7.2 双极化合成

对双极化阵列，单层候选向量可写成：

$$
\mathbf w_{l,m,\phi}
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
\bar{\mathbf v}_{l,m}\\
e^{j\phi}\bar{\mathbf v}_{l,m}
\end{bmatrix}.
$$

候选相位 $\phi$ 取自有限集合，例如：

$$
\phi\in
\left\{
0,\frac{\pi}{2},\pi,\frac{3\pi}{2}
\right\}.
$$

这相当于用有限比特量化两个极化之间的相对相位。

UE 选择：

$$
(l^\star,m^\star,\phi^\star)
=
\arg\max_{l,m,\phi}
\sum_{k\in\mathcal K}
\left\|
\hat{\mathbf H}[k]\mathbf w_{l,m,\phi}
\right\|_2^2.
$$

这个表达也说明了 Type I 码本和双极化信道建模的关系：如果真实信道中两个极化端口之间存在显著的相对相位和功率差，码本需要用 $\phi$ 等索引捕获这种极化合成方向。

## 7.3 宽带和子带相位

宽带信道的主导空间方向通常随频率变化较慢，而相位和等效增益可能随频率变化更明显。因此 Type I 码本可采用宽带方向 + 子带补充的思想：

$$
\mathbf w_b
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
\bar{\mathbf v}_{l,m}\\
e^{j\phi_b}\bar{\mathbf v}_{l,m}
\end{bmatrix}.
$$

其中：

- $(l,m)$ 可作为宽带 PMI；
- $\phi_b$ 可随子带 $b$ 变化；
- CQI 也可按宽带或子带反馈。

这种结构降低了反馈开销：不需要每个子带重新反馈完整空间波束，只需反馈少量相位或补充索引。

# 8 多层 Type I 预编码矩阵表达

## 8.1 多层传输与秩指示

多层传输时：

$$
\mathbf W
=
\begin{bmatrix}
\mathbf w_1 & \mathbf w_2 & \cdots & \mathbf w_L
\end{bmatrix}
\in
\mathbb C^{N_t\times L}.
$$

RI 表示 UE 推荐的层数 $L$。从信道矩阵角度看，$L$ 应和信道有效秩、奇异值分布、SINR 和接收机能力相关。

若奇异值满足：

$$
\sigma_1\gg\sigma_2\gg\cdots,
$$

则适合较低 rank；若多个奇异值都较大，则可支持更高 rank。

## 8.2 多波束/多极化层的矩阵构造

Type I 多层码本可以抽象理解为选择多个结构化列向量：

$$
\mathbf W
=
\left[
\mathbf w_{i_1},
\mathbf w_{i_2},
\dots,
\mathbf w_{i_L}
\right],
$$

其中每个 $\mathbf w_{i_\ell}$ 来自 Type I 单层候选集合或其规则化组合。

例如 rank-2 时，一个简单的双极化结构可以写成：

$$
\mathbf W
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
\bar{\mathbf v}_{l_1,m_1} & \bar{\mathbf v}_{l_2,m_2}\\
e^{j\phi_1}\bar{\mathbf v}_{l_1,m_1} &
e^{j\phi_2}\bar{\mathbf v}_{l_2,m_2}
\end{bmatrix}.
$$

若两个层使用相同空间波束但不同极化组合，则：

$$
\mathbf W
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
\bar{\mathbf v}_{l,m} & \bar{\mathbf v}_{l,m}\\
e^{j\phi_1}\bar{\mathbf v}_{l,m} &
e^{j\phi_2}\bar{\mathbf v}_{l,m}
\end{bmatrix}.
$$

但实际码本需要保证列之间具有合适的正交性或低相关性，否则层间干扰会增大。

## 8.3 归一化与功率约束

预编码矩阵通常需要满足总功率约束：

$$
\mathrm{tr}(\mathbf W^H\mathbf W)=L.
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
\mathbf W^H\mathbf W\approx\mathbf I_L.
$$

这有助于降低发射侧层间相关性，也便于接收侧均衡。

# 9 预编码矩阵的获取方法

Type I 码本并不是给出一个连续优化问题的任意解，而是规定一组离散、结构化的候选预编码矩阵。预编码矩阵的获取过程可以理解为：

$$
\boxed{
\hat{\mathbf H}[k]
\rightarrow
\text{UE 选择 RI/PMI/CQI}
\rightarrow
\text{gNB 根据 PMI 重建 } \mathbf W
}
$$

其中 UE 负责根据下行 CSI-RS 估计信道并选择推荐码本索引，gNB 负责根据标准码本规则和调度结果生成实际用于 PDSCH 的预编码矩阵。

## 9.1 Type I 码本的 SU-MIMO 定位

NR Type I 码本的基本目标是单用户下行预编码反馈。对一个 UE，gNB 希望选择：

$$
\mathbf W[k]\in\mathbb C^{N_t\times L}
$$

使该 UE 的等效信道

$$
\mathbf H[k]\mathbf W[k]
$$

在所选 rank $L$ 下具有较高的有效增益和较低的层间干扰。Type I 码本将这个问题限制在有限集合内：

$$
\mathbf W[k]\in\mathcal W.
$$

因此，Type I 的核心输出是“该 UE 推荐的单用户预编码方向”。在 MU-MIMO 中，gNB 可以把多个 UE 的 PMI、CQI、RI 作为调度参考，但最终多用户预编码还需要考虑用户间干扰，例如用户配对、波束正交性、功率分配和调度约束。换句话说，Type I 码本本身是 SU-MIMO 码本反馈机制，不是完整的 MU-MIMO 零迫或块对角化预编码算法。

## 9.2 宽带预编码矩阵获取

宽带预编码矩阵用于描述整个 CSI 带宽上相对稳定的主导空间方向。设 CSI 带宽内的子载波集合为 $\mathcal K_{\mathrm{WB}}$，UE 基于 CSI-RS 得到估计信道：

$$
\hat{\mathbf H}[k],\qquad k\in\mathcal K_{\mathrm{WB}}.
$$

一种宽带选择方式是先构造宽带发射端等效协方差：

$$
\mathbf R_{\mathrm{WB}}
=
\sum_{k\in\mathcal K_{\mathrm{WB}}}
\hat{\mathbf H}^{H}[k]\hat{\mathbf H}[k].
$$

若不考虑码本量化，理想 $L$ 层宽带预编码方向可由 $\mathbf R_{\mathrm{WB}}$ 的前 $L$ 个主特征向量给出：

$$
\mathbf R_{\mathrm{WB}}\mathbf v_i
=
\lambda_i\mathbf v_i,\qquad i=1,\dots,L,
$$

$$
\mathbf W_{\mathrm{opt,WB}}
=
\begin{bmatrix}
\mathbf v_1 & \mathbf v_2 & \cdots & \mathbf v_L
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
\hat{\mathbf H}[k],
\mathbf W_{\mathrm{WB}}(i_{\mathrm{WB}})
\right)
}
$$

其中 $\mathbf W_{\mathrm{WB}}(i_{\mathrm{WB}})\in\mathcal W_{\mathrm{WB}}$ 由 Type I 码本规则生成。常用度量可以是等效信道增益：

$$
\mathcal M
\left(
\hat{\mathbf H}[k],
\mathbf W
\right)
=
\left\|
\hat{\mathbf H}[k]\mathbf W
\right\|_F^2,
$$

也可以是吞吐率近似：

$$
\mathcal M
\left(
\hat{\mathbf H}[k],
\mathbf W
\right)
=
\log_2
\det
\left(
\mathbf I
+
\frac{\rho}{L}
\hat{\mathbf H}[k]\mathbf W\mathbf W^H\hat{\mathbf H}^{H}[k]
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

## 9.3 子带预编码矩阵获取

宽带 PMI 假设主导空间方向在整个 CSI 带宽内变化较慢。但在频率选择性信道中，不同子带可能存在不同的相位、增益和局部最优方向。因此 NR 可以在宽带 PMI 的基础上引入子带 PMI。

设第 $b$ 个子带的子载波集合为 $\mathcal K_b$，其子带等效协方差为：

$$
\mathbf R_b
=
\sum_{k\in\mathcal K_b}
\hat{\mathbf H}^{H}[k]\hat{\mathbf H}[k].
$$

理想子带预编码方向为：

$$
\mathbf W_{\mathrm{opt},b}
=
\mathrm{eig}_{L}
\left(
\mathbf R_b
\right),
$$

其中 $\mathrm{eig}_{L}(\cdot)$ 表示取前 $L$ 个主特征向量组成的矩阵。

在 Type I 反馈中，子带预编码矩阵通常不是完全独立重新选择一个任意矩阵，而是在宽带 PMI 约束下选择补充索引：

$$
i_{\mathrm{SB},b}^\star
=
\arg\max_{i_{\mathrm{SB},b}}
\sum_{k\in\mathcal K_b}
\mathcal M
\left(
\hat{\mathbf H}[k],
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

## 9.4 宽带方向与子带补充的组合

对单层双极化 Type I 码本，一个直观的宽带 + 子带组合可写成：

$$
\mathbf w_b
=
\frac{1}{\sqrt{2}}
\begin{bmatrix}
\bar{\mathbf v}_{l^\star,m^\star}\\
e^{j\phi_b^\star}
\bar{\mathbf v}_{l^\star,m^\star}
\end{bmatrix}.
$$

其中：

- $(l^\star,m^\star)$ 来自宽带 PMI，表示宽带主导 2D DFT 波束；
- $\phi_b^\star$ 来自第 $b$ 个子带的子带 PMI，表示该子带上的极化相位或补充相位；
- $\mathbf w_b$ 是第 $b$ 个子带的单层预编码向量。

更一般地，对 $L$ 层传输：

$$
\mathbf W_b
=
\begin{bmatrix}
\mathbf w_{b,1} &
\mathbf w_{b,2} &
\cdots &
\mathbf w_{b,L}
\end{bmatrix},
$$

每一列由 Type I 码本的空间波束、极化相位和层间结构共同决定。

从矩阵分解角度看，这个过程可以概括为：

$$
\boxed{
\mathbf V_L[k]
\approx
\mathbf W_b
=
f_{\mathrm{codebook}}
\left(
i_{\mathrm{WB}}^\star,
i_{\mathrm{SB},b}^\star,
L^\star
\right)
}
$$

其中 $\mathbf V_L[k]$ 是理想信道右奇异子空间，$\mathbf W_b$ 是 Type I 码本可表达的量化近似。

## 9.5 gNB 侧如何使用反馈得到的预编码矩阵

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
\mathbf x[k]
=
\mathbf W_{\mathrm{WB}}\mathbf s[k],
\qquad
k\in\mathcal K_{\mathrm{sched}}.
$$

若采用子带预编码，则不同子带使用不同的 $\mathbf W_b$：

$$
\mathbf x[k]
=
\mathbf W_b\mathbf s[k],
\qquad
k\in\mathcal K_b.
$$

因此，Type I 码本下的预编码矩阵获取可以分成两个层次：

1. UE 侧：基于信道估计结果选择最优或近似最优的 RI/PMI/CQI；
2. gNB 侧：基于 PMI 和标准码本规则重建宽带或子带预编码矩阵，并将其用于下行 PDSCH。

# 10 Type I 码本反馈流程

Type I 码本反馈可以抽象为以下步骤。

**Step 1：gNB 配置 CSI-RS 和 CSI report**

gNB 配置 CSI-RS 端口、带宽、报告周期、码本类型、可支持 rank、PMI/CQI 粒度等。

**Step 2：UE 估计下行信道**

UE 根据 CSI-RS 得到：

$$
\hat{\mathbf H}[k].
$$

**Step 3：UE 对候选 rank 做评估**

对 $L=1,2,\dots,L_{\max}$，UE 评估不同 rank 下的码本候选。

**Step 4：UE 搜索 Type I 码本**

对每个候选 $\mathbf W_i$，计算度量：

$$
\mathcal M_i
=
\sum_{k\in\mathcal K}
\log_2
\det
\left(
\mathbf I
+
\frac{\rho}{L}
\hat{\mathbf H}[k]\mathbf W_i
\mathbf W_i^H\hat{\mathbf H}^H[k]
\right).
$$

选择最大度量对应的 rank 和 PMI：

$$
(L^\star,i^\star)
=
\arg\max_{L,i}\mathcal M_i.
$$

**Step 5：UE 反馈 RI/PMI/CQI**

反馈内容包括：

- RI：推荐 rank；
- PMI：推荐预编码矩阵索引；
- CQI：在该 rank/PMI 下可支持的链路质量。

**Step 6：gNB 重建预编码矩阵**

gNB 根据 PMI 和标准码本规则得到：

$$
\mathbf W^\star
=
f_{\mathrm{codebook}}(i^\star).
$$

然后用于 PDSCH 下行传输。

# 11 Type I 码本与信道建模的关系

Type I 码本设计和 MIMO 信道建模密切相关。

**1. 与角度域信道建模的关系**

Type I 使用 DFT 波束，本质上利用了 ULA/UPA 阵列的角度域结构。若信道具有有限个主导 AoD：

$$
\mathbf H[k]
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
\mathbf H[k]
=
\sum_{\ell}
\sqrt{P_\ell}
\mathbf A_r(\Omega_{r,\ell})
\mathbf P_\ell
\mathbf A_t^H(\Omega_{t,\ell})
e^{-j\frac{2\pi}{N}k\tau_\ell}.
$$

其中 $\mathbf P_\ell$ 是 $2\times2$ 极化耦合矩阵。Type I 码本中的极化相位和双极化端口组合，正是为了用低开销方式适配这种极化耦合结构。

若 XPR 很高，交叉极化泄漏弱，两个极化之间相对独立，码本中极化相位的影响可能较小。若 XPR 较低，交叉极化较强，合适的极化合成对性能影响更明显。

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
- 极化相位量化 $\mathbf R_{t,\mathrm{pol}}$ 或极化耦合导致的主导极化方向；
- rank/层结构量化多个主导特征方向。

**4. 与信道估计误差的关系**

UE 搜索码本依赖 $\hat{\mathbf H}$。若信道估计误差为：

$$
\hat{\mathbf H}[k]
=
\mathbf H[k]+\mathbf E[k],
$$

则 PMI 选择可能发生错误：

$$
\arg\max_i
\mathcal M(\hat{\mathbf H},\mathbf W_i)
\ne
\arg\max_i
\mathcal M(\mathbf H,\mathbf W_i).
$$

因此，CSI-RS 密度、信道插值、LMMSE/Wiener 估计、双选择性衰落跟踪都会影响 Type I 码本反馈质量。

# 12 小结

NR Type I 码本是面向低开销 CSI 反馈的结构化预编码码本。它的核心不是反馈完整信道，而是用有限索引表示一个适合下行传输的预编码矩阵。

从数学上看，理想预编码来自信道矩阵的右奇异向量：

$$
\mathbf H
=
\mathbf U\mathbf \Sigma\mathbf V^H,
\qquad
\mathbf W_{\mathrm{opt}}
=
\mathbf V_L.
$$

Type I 码本用标准化结构近似 $\mathbf V_L$：

$$
\mathbf W_{\mathrm{opt}}
\approx
\mathbf W(\mathrm{PMI}).
$$

其主要设计思想包括：

- 用 2D DFT 波束表示空间方向；
- 用有限相位集合表示双极化端口之间的合成；
- 用 RI 表示传输层数；
- 用宽带/子带 PMI 兼顾反馈开销和频率选择性；
- 用 CQI 反映该预编码假设下的链路质量。

在 NR FDD 系统中，Type I 码本是下行信道获取和预编码反馈的重要机制。在 TDD 系统中，当上行覆盖受限、SRS 不可靠或校准误差较大时，Type I 码本也可以作为下行 CSI 获取的补充手段。

与信道建模的关系可以概括为：

- MIMO-OFDM 信道建模给出 $\mathbf H[m,k]$；
- 信道估计给出 $\hat{\mathbf H}[m,k]$；
- SVD 或互信息准则给出理想预编码方向；
- Type I 码本用 DFT 波束、极化相位和层结构对理想预编码方向做有限比特量化。

因此，理解 Type I 码本需要同时理解 MIMO-OFDM 信道矩阵、双极化阵列结构、信道估计误差和有限反馈预编码之间的关系。

# 13 参考文献

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
