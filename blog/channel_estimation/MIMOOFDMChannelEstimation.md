# MIMO-OFDM 导频设计、信道估计算法与信道插值

> 本文聚焦 MIMO-OFDM 系统中的导频设计、信道估计算法和时频空/角度选择性信道插值。信道建模部分可参考 [MIMO OFDM channel modeling.md](MIMO%20OFDM%20channel%20modeling.md)。

## 目录

- [1 导频设计](#1-导频设计)
  - [1.1 系统模型](#11-系统模型)
  - [1.2 为什么 MIMO 需要可分离导频](#12-为什么-mimo-需要可分离导频)
  - [1.3 常见导频结构](#13-常见导频结构)
  - [1.4 导频间隔设计](#14-导频间隔设计)
- [2 信道估计算法](#2-信道估计算法)
  - [2.1 基于 CIR/CFR 的联合观测模型](#21-基于-circfr-的联合观测模型)
  - [2.2 ZF/LS 导频 CFR 初估计](#22-zfls-导频-cfr-初估计)
  - [2.3 MMSE 信道估计](#23-mmse-信道估计)
  - [2.4 LMMSE 信道估计](#24-lmmse-信道估计)
  - [2.5 OMP 角度-时延稀疏信道估计](#25-omp-角度-时延稀疏信道估计)
- [3 选择性衰落信道重构与插值](#3-选择性衰落信道重构与插值)
  - [3.1 从导频 RE 到数据 RE 的模型化重构](#31-从导频-re-到数据-re-的模型化重构)
  - [3.2 基于 IDFT/DFT 的频率选择性重构](#32-基于-idftdft-的频率选择性重构)
  - [3.3 基于时间选择性的符号级插值](#33-基于时间选择性的符号级插值)
  - [3.4 时频空 Wiener 预测](#34-时频空-wiener-预测)
  - [3.5 角度-时延域重构与去噪](#35-角度-时延域重构与去噪)
  - [3.6 块衰落信道下的简化](#36-块衰落信道下的简化)
- [4 附录](#4-附录)

# 1 导频设计

导频设计的目标是在有限时频资源内，使接收端能够区分不同发射天线的信道，并且在时间、频率和空间维度上获得足够采样。MIMO-OFDM 的导频设计比 SISO 更复杂，因为同一接收天线上会叠加来自多根发射天线的信号。

## 1.1 系统模型

先给出后续推导默认使用的仿真参数和基本假设。

| 参数 | 含义 | 典型约束或关系 |
|---|---|---|
| $f_c$ | 载波频率 | 决定波长 $\lambda=c/f_c$ 和最大多普勒 |
| $B$ | 系统带宽 | 近似 $B=N_{\mathrm{act}}\Delta f$ |
| $\Delta f$ | 子载波间隔 | 有效 OFDM 符号时长 $T_u=1/\Delta f$ |
| $N_{\mathrm{FFT}}$ | FFT 点数 | 本文简记为 $N$ |
| $N_{\mathrm{act}}$ | 有效子载波数 | $N_{\mathrm{act}}\le N_{\mathrm{FFT}}$ |
| $T_{\mathrm{cp}}$ | 循环前缀长度 | 应覆盖主要多径时延扩展 |
| $T_{\mathrm{sym}}$ | 含 CP 的 OFDM 符号时长 | $T_{\mathrm{sym}}=T_u+T_{\mathrm{cp}}$ |
| $N_t,N_r$ | 发射、接收天线数 | 决定 MIMO 信道矩阵维度 |
| $N_{\mathrm{sy}}$ | 一个训练块内的连续 OFDM 符号数 | 块内信道近似不变 |
| $L_h$ | 离散时延 tap 数 | 通常由最大时延和采样率决定 |
| $\tau_{\max}$ | 最大多径时延 | 影响频域导频间隔 |
| $v$ | 移动速度 | 与 $f_c$ 共同决定多普勒 |
| $f_D$ | 最大多普勒频移 | $f_D=vf_c/c$ |
| $\sigma_w^2$ | 噪声功率 | 决定导频接收信号的 SNR |

本文默认以下条件成立：

- CP 足够长，使 OFDM 子载波间干扰可忽略；
- 一个训练块内满足块衰落近似，即 $N_{\mathrm{sy}}T_{\mathrm{sym}}$ 小于信道相干时间；
- 宽带信道由有限个离散时延 tap 描述，不同子载波由同一组离散时延 tap 经 DFT 得到；
- 若使用统计估计，PDP、多普勒谱和空间相关矩阵可由信道模型或测量获得；
- 若使用角度-时延稀疏估计，阵列响应和角度字典与信道建模文档保持一致。

考虑 $N_t$ 发射天线、$N_r$ 接收天线、$N=N_{\mathrm{FFT}}$ 个 OFDM 子载波。实际估计通常以一个训练块或相干块为单位：假设连续 $N_{\mathrm{sy}}$ 个 OFDM 符号内信道近似不变。记块索引为 $b$，块内符号索引为 $q=0,1,\dots,N_{\mathrm{sy}}-1$，第 $k$ 个子载波上的频域模型为：

$$
\mathbf y_b[q,k]
=
\mathbf H_b[k]\mathbf x_b[q,k]+\mathbf w_b[q,k],
$$

其中：

- $\mathbf y_b[q,k]\in\mathbb C^{N_r\times 1}$ 是接收向量；
- $\mathbf x_b[q,k]\in\mathbb C^{N_t\times 1}$ 是发射向量；
- $\mathbf H_b[k]\in\mathbb C^{N_r\times N_t}$ 是第 $b$ 个训练块内共享的 MIMO 频域信道；
- $\mathbf w_b[q,k]\sim\mathcal{CN}(\mathbf 0,\sigma_w^2\mathbf I)$ 是噪声。

块内共享信道的假设为：

$$
\mathbf H[b,0,k]
=
\mathbf H[b,1,k]
=
\cdots
=
\mathbf H[b,N_{\mathrm{sy}}-1,k]
\triangleq
\mathbf H_b[k].
$$

这不是说每个 OFDM 符号的信道独立变化，而是说在 $N_{\mathrm{sy}}T_{\mathrm{sym}}$ 时间内信道相干，可用多个连续符号上的导频共同估计同一个宽带 MIMO 信道。

若只考虑某一组导频资源，且信道在这些资源内近似不变，可写成局部块训练模型：

$$
\boxed{
\mathbf Y=\mathbf H\mathbf X+\mathbf W
}
$$

其中：

- $\mathbf Y\in\mathbb C^{N_r\times N_p}$；
- $\mathbf H\in\mathbb C^{N_r\times N_t}$；
- $\mathbf X\in\mathbb C^{N_t\times N_p}$；
- $N_p$ 是导频接收信号数量。

这个模型是窄带或局部块初始估计的特例。后续信道估计算法会进一步采用宽带 MIMO-OFDM 建模，把同一训练块内不同符号、不同子载波上的导频都作为同一组离散时延 tap 或角度-时延路径的观测，而不是把每个符号、每个子载波看成互相独立的未知信道。

## 1.2 为什么 MIMO 需要可分离导频

在 SISO-OFDM 中，导频 RE 上的 CFR 初估计可以写成：

$$
Y_p=X_pH_p+W_p.
$$

若 $X_p\ne0$，直接相除即可得到 LS 估计：

$$
\hat H_p=\frac{Y_p}{X_p}.
$$

但在 MIMO 中，某个导频 RE 上的接收信号是多根发射天线或多个端口的叠加：

$$
\mathbf y_p
=
\sum_{t=1}^{N_t}\mathbf h_t x_t+\mathbf w_p,
$$

其中 $\mathbf h_t$ 是当前导频子载波 CFR 矩阵 $\mathbf H_b[k_p]$ 的第 $t$ 列。如果多个发射天线使用不可区分的导频，接收端无法判断接收能量分别来自哪根天线或端口。

因此在导频 CFR 初估计阶段，局部导频矩阵至少需要满行秩：

$$
\boxed{
\operatorname{rank}(\mathbf X)=N_t
}
$$

这要求局部导频接收信号数量满足：

$$
N_p\ge N_t.
$$

更常用、更稳健的设计是正交训练：

$$
\boxed{
\mathbf X\mathbf X^H=E_p\mathbf I_{N_t}
}
$$

此时不同发射天线的导频互不干扰，ZF/LS 可以先得到导频 RE 上的 CFR 初估计。

对宽带信道估计而言，仅让不同发射天线的导频相互正交还不够。正交性只能保证接收端可以先分离不同发射端口的导频 CFR；频域上还必须为每个发射端口放置足够数量的导频子载波，才能从这些 CFR 初估计恢复有限时延信道。

对某个固定的接收天线 $r$ 和发射端口 $t$，设该发射端口在一个训练块内占用的导频子载波集合为：

$$
\mathcal P_t=\{k_1,k_2,\dots,k_{N_{p,t}}\}.
$$

在这些导频子载波上，ZF/LS 初估计得到的 CFR 向量定义为：

$$
\hat{\mathbf h}_{p}
=
\begin{bmatrix}
\hat H_{r,t}[k_1] &
\hat H_{r,t}[k_2] &
\cdots &
\hat H_{r,t}[k_{N_{p,t}}]
\end{bmatrix}^T
=
\mathbf F_{\mathcal P,L_h}\mathbf h_{\tau}
+\mathbf e_p,
$$

其中：

- $\hat{\mathbf h}_{p}\in\mathbb C^{N_{p,t}\times 1}$ 是导频位置上的 CFR 初估计；
- $\mathbf h_{\tau}\in\mathbb C^{L_h\times 1}$ 是该天线对的离散时延 tap 向量；
- $\mathbf F_{\mathcal P,L_h}\in\mathbb C^{N_{p,t}\times L_h}$ 是由导频子载波索引抽取出的部分 DFT 矩阵，其第 $i$ 行、第 $n$ 列为 $e^{-j2\pi k_i n/N}$；
- $\mathbf e_p$ 是导频 CFR 初估计误差。

若不引入稀疏先验，只做普通有限离散时延 tap 的 LS/DFT 重构，则至少需要：

$$
\boxed{
N_{p,t}\ge L_h
}
$$

并且部分 DFT 矩阵必须满列秩：

$$
\boxed{
\operatorname{rank}(\mathbf F_{\mathcal P,L_h})=L_h
}
$$

这里的 $L_h$ 是有效离散时延 tap 数，来自最大时延扩展和采样率，不是物理路径数。多个物理路径可能落在同一个离散时延 tap，off-grid 路径也可能泄漏到多个离散时延 tap。因此更准确的说法是：每个发射端口的频域导频数量应不少于需要估计的离散时延 tap 数；导频数量越多，噪声平均和 LS 稳定性越好。

若采用角度-时延稀疏估计，才需要额外假设有效路径数或稀疏度。此时导频接收信号数量不一定要大于完整的 $L_h$，但必须足以恢复目标稀疏支撑，并且测量矩阵 $\mathbf \Phi$ 不能让不同角度-时延原子高度混叠。

## 1.3 常见导频结构

MIMO-OFDM 中常见导频结构主要可以分为时分、频分和码分。它们对应导频在时间、频率和码域上的正交化方式。

**时分导频**  
不同发射天线在不同 OFDM 符号发送导频：

| OFDM 符号 | Tx1 | Tx2 | Tx3 |
|---|---|---|---|
| $m_1$ | pilot | 0 | 0 |
| $m_2$ | 0 | pilot | 0 |
| $m_3$ | 0 | 0 | pilot |

优点是简单，缺点是时间开销较大，且高速移动时不同天线的估计不在同一时刻。

**频分导频**  
不同发射天线占用不同子载波集合：

$$
\mathcal P_i\cap\mathcal P_j=\varnothing,\qquad i\ne j.
$$

优点是可在同一 OFDM 符号内估计多根天线，缺点是每根天线的频域采样密度下降。

**码分导频**  
不同发射天线在同一组资源上使用正交序列：

$$
\sum_{p\in\mathcal P}
x_i[p]x_j^*[p]
=0,\qquad i\ne j.
$$

优点是资源利用率较高，缺点是对频偏、同步误差和强频率选择性更敏感。

## 1.4 导频间隔设计

导频间隔由信道的时间选择性和频率选择性决定。

**频域间隔**  
令 $S_f$ 表示频域上相邻导频之间间隔的子载波数。实际导频频率间隔为：

$$
\Delta f_{\mathrm{pilot}}
=
S_f\Delta f.
$$

其中 $\Delta f$ 是子载波间隔，$\Delta f_{\mathrm{pilot}}$ 的单位是 Hz。为跟踪频率选择性信道，导频频率间隔应小于相干带宽量级：

$$
\boxed{
\Delta f_{\mathrm{pilot}}
=
S_f\Delta f
\lesssim
\frac{1}{\tau_{\max}}
}
$$

其中 $\tau_{\max}$ 是最大多径时延。

**时域间隔**  
一个训练块包含连续 $N_{\mathrm{sy}}$ 个 OFDM 符号。若相邻导频训练块的间隔为 $S_b$ 个训练块，则实际导频时间间隔为：

$$
\Delta t_{\mathrm{pilot}}
=
S_bN_{\mathrm{sy}}T_{\mathrm{sym}}.
$$

其中 $T_{\mathrm{sym}}$ 是含 CP 的 OFDM 符号时长。为跟踪时间选择性信道，导频时间间隔应小于相干时间量级：

$$
\boxed{
\Delta t_{\mathrm{pilot}}
=
S_bN_{\mathrm{sy}}T_{\mathrm{sym}}
\lesssim
\frac{1}{2f_D}
}
$$

其中 $f_D$ 是最大多普勒频移。

**角度域影响**  
角度域选择性衰落不直接决定 $\Delta f_{\mathrm{pilot}}$ 或 $\Delta t_{\mathrm{pilot}}$。频域间隔主要由时延扩展决定，时域间隔主要由多普勒扩展决定。角度域影响的是空间维度上的导频测量资源，即需要多少彼此线性独立的发射端口、空间层、预编码方向或接收合并方向来估计并利用 MIMO 信道的空间自由度。

若角度扩展较大，AoA/AoD 分布更分散，信道空间秩和角度域有效支撑通常更大。此时如果系统希望利用这些空间自由度，就需要配置足够维度的导频端口或预编码测量资源。若角度扩展很小，信道可能低秩或高度相关，增加端口数的边际收益较小。

因此不能简单说“角度扩展大就必须增加导频端口数”。端口数由系统要估计的 CSI 维度、天线端口配置和空间层数决定；角度域选择性决定这些空间维度是否有效、是否容易分离，以及角度-时延稀疏估计所需的导频接收信号数量。

MIMO 不改变时间和频率采样的基本原则，但会增加空间维度上的导频正交和测量需求。若每根发射天线或每个端口都需要独立导频，则导频开销通常随待估计的端口维度增加。

# 2 信道估计算法

本节以实际宽带 MIMO-OFDM 接收机的处理链路为主线：先在导频 RE 上估计 CFR 初值，再利用有限时延、PDP/多普勒/空间相关或角度-时延稀疏结构，对这些 CFR 初值做宽带滤波、预测和重构。窄带模型是退化特例；宽带处理中，不同子载波的 CFR 由同一组离散时延 tap 或同一组物理路径共同决定。

## 2.1 基于 CIR/CFR 的联合观测模型

离散 MIMO CIR 写为：

$$
\mathbf H_b[n]\in\mathbb C^{N_r\times N_t},
\qquad n=0,1,\dots,L_h-1.
$$

第 $b$ 个训练块、第 $k$ 个子载波上的 CFR 是沿离散时延维度做 DFT 得到的：

$$
\boxed{
\mathbf H_b[k]
=
\sum_{n=0}^{L_h-1}
\mathbf H_b[n]e^{-j2\pi kn/N}
}
$$

因此同一训练块内，不同符号共享同一个信道，不同子载波也不是独立未知量，而是由同一组 $\{\mathbf H_b[n]\}$ 共同决定。对块内某个导频 RE $(q_i,k_i)$，接收模型为：

$$
\mathbf y_i
=
\mathbf H_b[k_i]\mathbf x_i+\mathbf w_i.
$$

其中：

- $i$ 是训练块内第 $i$ 个导频 RE 的索引，对应位置 $(q_i,k_i)$；
- $\mathbf y_i\in\mathbb C^{N_r\times 1}$ 是该导频 RE 上的接收频域信号；
- $\mathbf x_i\in\mathbb C^{N_t\times 1}$ 是该导频 RE 上的发射导频向量；
- $\mathbf H_b[k_i]\in\mathbb C^{N_r\times N_t}$ 是第 $b$ 个训练块、第 $k_i$ 个子载波上的 MIMO CFR；
- $\mathbf w_i\in\mathbb C^{N_r\times 1}$ 是该导频 RE 上的噪声向量。

定义离散时延 tap 向量：

$$
\mathbf h_b
=
\operatorname{vec}
\left(
\mathbf H_b[0],\mathbf H_b[1],\dots,\mathbf H_b[L_h-1]
\right)
\in\mathbb C^{N_rN_tL_h\times 1}.
$$

令：

$$
\mathbf a(k)
=
\left[
1,e^{-j2\pi k/N},\dots,e^{-j2\pi k(L_h-1)/N}
\right]^T.
$$

则：

$$
\operatorname{vec}(\mathbf H_b[k])
=
(\mathbf a^T(k)\otimes\mathbf I_{N_rN_t})\mathbf h_b.
$$

代入 $\mathbf y_i=(\mathbf x_i^T\otimes\mathbf I_{N_r})\operatorname{vec}(\mathbf H_b[k_i])+\mathbf w_i$，得到单个导频接收信号：

$$
\mathbf y_i
=
\underbrace{
(\mathbf x_i^T\otimes\mathbf I_{N_r})
(\mathbf a^T(k_i)\otimes\mathbf I_{N_rN_t})
}_{\mathbf A_i}
\mathbf h_b
+\mathbf w_i.
$$

其中 $\mathbf A_i\in\mathbb C^{N_r\times N_rN_tL_h}$ 是第 $i$ 个导频 RE 对宽带离散时延 tap 向量 $\mathbf h_b$ 的观测矩阵。

设第 $b$ 个训练块内共有 $N_{\mathrm{obs}}$ 个导频 RE。这里的 $N_{\mathrm{obs}}$ 是时间-频率导频位置数量，例如所有 $(q_i,k_i)$ 的个数；它不包含接收天线维度。每个导频 RE 给出一个 $\mathbf y_i\in\mathbb C^{N_r\times 1}$。将这些接收向量按列方向纵向堆叠：

$$
\mathbf y
=
\begin{bmatrix}
\mathbf y_1\\
\mathbf y_2\\
\vdots\\
\mathbf y_{N_{\mathrm{obs}}}
\end{bmatrix}
\in\mathbb C^{N_rN_{\mathrm{obs}}\times 1}.
$$

对应地，将每个 $\mathbf A_i\in\mathbb C^{N_r\times N_rN_tL_h}$ 纵向堆叠：

$$
\mathbf A
=
\begin{bmatrix}
\mathbf A_1\\
\mathbf A_2\\
\vdots\\
\mathbf A_{N_{\mathrm{obs}}}
\end{bmatrix}
\in
\mathbb C^{N_rN_{\mathrm{obs}}\times N_rN_tL_h}.
$$

噪声向量同样堆叠为：

$$
\mathbf w
=
\begin{bmatrix}
\mathbf w_1\\
\mathbf w_2\\
\vdots\\
\mathbf w_{N_{\mathrm{obs}}}
\end{bmatrix}
\in\mathbb C^{N_rN_{\mathrm{obs}}\times 1}.
$$

于是所有导频接收信号可统一写成：

$$
\boxed{
\mathbf y=\mathbf A\mathbf h+\mathbf w
}
$$

在本节的天线-时延域模型中，$\mathbf h=\mathbf h_b\in\mathbb C^{N_rN_tL_h\times 1}$。更一般地，$\mathbf h$ 是用于宽带建模的信道参数，可以是：

- 天线-离散时延域（tap 域）向量；
- 时频空联合信道向量；
- 角度-时延域稀疏向量。

实际实现中，ZF/LS 常先给出导频子载波上的 CFR 初估计；LMMSE/Wiener、IDFT/DFT 重构和 OMP 再把这些导频接收信号或 CFR 初估计放入宽带模型中，生成目标 RE 上的 $\hat{\mathbf H}[q,k]$。若采用块衰落简化，$\hat{\mathbf H}[q,k]$ 可退化为块内共享的 $\hat{\mathbf H}_b[k]$。

## 2.2 ZF/LS 导频 CFR 初估计

在实际 OFDM 接收机中，ZF/LS 最常见的用途是在导频 RE 上先估计 CFR。对某个导频子载波 $k_p$，若经过时分、频分或码分设计后，第 $t$ 个发射端口的导频可被分离，则有：

$$
\mathbf y_b[q,k_p]
=
\mathbf H_b[k_p]\mathbf x_b[q,k_p]+\mathbf w_b[q,k_p].
$$

对固定训练块 $b$ 和固定导频子载波 $k_p$，将由正交导频产生的 $N_p$ 次接收频域信号表示为矩阵形式：

$$
\mathbf Y_b[k_p]
=
\mathbf H_b[k_p]\mathbf X_b[k_p]
+\mathbf W_b[k_p].
$$

其中：

- $\mathbf Y_b[k_p]\in\mathbb C^{N_r\times N_p}$ 是第 $b$ 个训练块、第 $k_p$ 个导频子载波上收集到的频域接收信号矩阵，每一列对应一次导频接收信号；
- $\mathbf H_b[k_p]\in\mathbb C^{N_r\times N_t}$ 是第 $b$ 个训练块、第 $k_p$ 个子载波上的 MIMO CFR；
- $\mathbf X_b[k_p]\in\mathbb C^{N_t\times N_p}$ 是第 $b$ 个训练块、第 $k_p$ 个子载波上的发射导频矩阵，每一列是一次导频发送向量；
- $\mathbf W_b[k_p]\in\mathbb C^{N_r\times N_p}$ 是对应的噪声矩阵；
- $N_p$ 是该子载波上用于分离 $N_t$ 个发射端口的导频观测次数，需要满足 $N_p\ge N_t$。

若 $\mathbf X_b[k_p]$ 满行秩，即 $\operatorname{rank}(\mathbf X_b[k_p]) = N_t$，则导频位置上的 ZF/LS CFR 初估计为：

$$
\boxed{
\hat{\mathbf H}_{\mathrm{LS},b}[k_p]
=
\mathbf Y_b[k_p]\mathbf X_b[k_p]^H
\left(\mathbf X_b[k_p]\mathbf X_b[k_p]^H\right)^{-1}
}
$$

若导频正交，$\mathbf X_b[k_p]\mathbf X_b[k_p]^H=E_p\mathbf I$，则：

$$
\boxed{
\hat{\mathbf H}_{\mathrm{LS},b}[k_p]
=
\frac{1}{E_p}\mathbf Y_b[k_p]\mathbf X_b[k_p]^H
}
$$

这一步的估计对象是导频子载波上的 CFR，而不是物理路径数，也不是必然直接估计 CIR。它只完成“导频位置初估计”。宽带结构的利用属于后续重构/插值步骤：可以把这些导频 CFR 初估计送入第 3 章的 IDFT/DFT 有限时延重构、Wiener/LMMSE 时频空预测或 OMP 角度-时延重构。

因此，ZF/LS 本身不需要假定物理路径数。若后续使用 IDFT/DFT 有限时延重构，需要假定的是最大有效离散时延长度 $L_h$，不是物理路径数；物理路径数或稀疏度假设主要出现在 OMP 一类稀疏估计算法中。

**算法伪代码**

```text
Input:
  Received pilot matrices Y_b[k_p] on pilot subcarriers
  Pilot matrices X_b[k_p]

Steps:
  1. Separate pilots by time/frequency/code orthogonality
  2. For each pilot subcarrier k_p:
       estimate CFR H_LS_b[k_p] = Y_b[k_p] X_b[k_p]^H (X_b[k_p] X_b[k_p]^H)^(-1)
  3. Pass pilot CFR estimates to the selected wideband reconstruction method if needed

Output:
  Pilot-position CFR estimates H_LS_b[k_p]
```

ZF/LS 的优点是实现简单、与实际导频处理一致；缺点是导频 CFR 初估计本身不利用 PDP、多普勒和空间相关，低 SNR 或导频稀疏时需要后续宽带重构或时频空插值。

## 2.3 MMSE 信道估计

MMSE 估计的目标是在完整先验分布下求后验均值：

$$
\boxed{
\hat{\mathbf h}_{\mathrm{MMSE}}
=
\mathbb E[\mathbf h\mid \mathbf y]
}
$$

这里的 $\mathbf h$ 不再是单个子载波的信道矩阵向量，而是宽带 MIMO 信道参数。它的先验可来自信道建模文档中的模型：

- i.i.d. Rayleigh/Rician 离散时延 tap 模型：先验由 PDP 和 K 因子决定；
- Kronecker 空间相关模型：先验由 $\mathbf R_r,\mathbf R_t$ 和 PDP 决定；
- 几何角度-时延模型：先验由路径数、AoA/AoD、时延、路径增益分布决定；
- CDL/TDL 模型：先验由标准化簇、射线、时延、角度和多普勒参数决定。

沿用 2.1 的堆叠模型时，各变量维度为：

- $\mathbf y\in\mathbb C^{N_rN_{\mathrm{obs}}\times 1}$：训练块内所有导频接收信号的堆叠向量；
- $\mathbf A\in\mathbb C^{N_rN_{\mathrm{obs}}\times N_h}$：联合观测矩阵；
- $\mathbf h\in\mathbb C^{N_h\times 1}$：待估计的宽带信道参数向量；
- $\mathbf w\in\mathbb C^{N_rN_{\mathrm{obs}}\times 1}$：噪声向量；
- $\mathbf R_w\in\mathbb C^{N_rN_{\mathrm{obs}}\times N_rN_{\mathrm{obs}}}$：噪声协方差矩阵。

若 $\mathbf h$ 采用天线-离散时延域（tap 域）表示，则 $N_h=N_rN_tL_h$。

MMSE 由贝叶斯公式给出：

$$
p(\mathbf h\mid\mathbf y)
=
\frac{p(\mathbf y\mid\mathbf h)p(\mathbf h)}{p(\mathbf y)}.
$$

其中：

$$
p(\mathbf y\mid\mathbf h)
\propto
\exp\left(
-\|\mathbf y-\mathbf A\mathbf h\|_{\mathbf R_w^{-1}}^2
\right).
$$

所以：

$$
\boxed{
\hat{\mathbf h}_{\mathrm{MMSE}}
=
\int \mathbf h\,p(\mathbf h\mid\mathbf y)\,d\mathbf h
}
$$

一般几何信道或 CDL 信道的先验并不一定是简单高斯分布，此时 MMSE 可能没有闭式解，需要通过网格积分、粒子滤波、EM、消息传递或神经网络近似。若 $\mathbf h$ 和 $\mathbf w$ 联合高斯，MMSE 退化为 LMMSE。

**算法伪代码**

```text
Input:
  Pilot observation vector y
  Joint observation matrix A
  Prior distribution p(h) from the MIMO-OFDM channel model
  Noise covariance Rw

Steps:
  1. Build likelihood p(y | h) from y = A h + w
  2. Build posterior p(h | y) proportional to p(y | h) p(h)
  3. Compute posterior mean h_hat = E[h | y]
  4. Reshape h_hat into CIR taps or model parameters
  5. Generate H_hat_b[k] from the estimated model

Output:
  H_hat_b[k] on the required subcarriers in the block
```

MMSE 是理论最优均方误差估计；它是否可实现，取决于信道模型先验是否足够简单。

## 2.4 LMMSE 信道估计

LMMSE 只使用二阶统计量，把估计器限制为线性形式：

$$
\hat{\mathbf h}=\mathbf B\mathbf y.
$$

沿用 2.3 的记号，$\mathbf y\in\mathbb C^{N_rN_{\mathrm{obs}}\times 1}$，$\mathbf h\in\mathbb C^{N_h\times 1}$。因此：

- $\mathbf B\in\mathbb C^{N_h\times N_rN_{\mathrm{obs}}}$ 是线性估计矩阵；
- $\mathbf R_h\in\mathbb C^{N_h\times N_h}$ 是宽带信道参数协方差；
- $\mathbf R_w\in\mathbb C^{N_rN_{\mathrm{obs}}\times N_rN_{\mathrm{obs}}}$ 是噪声协方差；
- $\mathbf R_{hy}\in\mathbb C^{N_h\times N_rN_{\mathrm{obs}}}$；
- $\mathbf R_{yy}\in\mathbb C^{N_rN_{\mathrm{obs}}\times N_rN_{\mathrm{obs}}}$。

对联合观测模型 $\mathbf y=\mathbf A\mathbf h+\mathbf w$，正交性原理给出：

$$
\mathbb E[
(\mathbf h-\mathbf B\mathbf y)\mathbf y^H
]
=\mathbf 0.
$$

因此：

$$
\mathbf B
=
\mathbf R_{hy}\mathbf R_{yy}^{-1}.
$$

若 $\mathbf h$ 与 $\mathbf w$ 不相关：

$$
\mathbf R_{hy}
=
\mathbf R_h\mathbf A^H,
\qquad
\mathbf R_{yy}
=
\mathbf A\mathbf R_h\mathbf A^H+\mathbf R_w.
$$

得到联合 LMMSE：

$$
\boxed{
\hat{\mathbf h}_{\mathrm{LMMSE}}
=
\mathbf R_h\mathbf A^H
(\mathbf A\mathbf R_h\mathbf A^H+\mathbf R_w)^{-1}
\mathbf y
}
$$

关键在于 $\mathbf R_h$ 应来自 MIMO-OFDM 信道建模，而不是假设各子载波独立。若采用 WSSUS + Kronecker 空间相关模型，离散时延域（tap 域）协方差可写为：

$$
\boxed{
\mathbf R_h
=
\mathbf R_{\mathrm{time}}
\otimes
\mathbf R_{\mathrm{delay}}
\otimes
\mathbf R_t^T
\otimes
\mathbf R_r
}
$$

其中：

- $\mathbf R_{\mathrm{time}}$ 描述训练块之间的时间相关；若只估计单个训练块，可省略该项；
- $\mathbf R_{\mathrm{delay}}\in\mathbb C^{L_h\times L_h}$ 通常由 PDP 决定；若不同离散时延 tap 独立，则 $\mathbf R_{\mathrm{delay}}=\operatorname{diag}(P_0,\dots,P_{L_h-1})$；
- $\mathbf R_t\in\mathbb C^{N_t\times N_t}$ 是发射端空间相关矩阵；
- $\mathbf R_r\in\mathbb C^{N_r\times N_r}$ 是接收端空间相关矩阵。

若只估计单个训练块且采用天线-离散时延域（tap 域）向量，常用简化形式为：

$$
\mathbf R_h
=
\mathbf R_{\mathrm{delay}}
\otimes
\mathbf R_t^T
\otimes
\mathbf R_r
\in\mathbb C^{N_rN_tL_h\times N_rN_tL_h}.
$$

这些相关矩阵可以通过模型参数或测量估计获得：

**时间相关矩阵**  
若已知最大多普勒 $f_D$，可用 Jakes/Clarke 模型生成训练块之间的时间相关。第 $b_1$、$b_2$ 个训练块之间的相关系数可写为：

$$
[\mathbf R_{\mathrm{time}}]_{b_1,b_2}
=
J_0\!\left(
2\pi f_DN_{\mathrm{sy}}T_{\mathrm{sym}}|b_1-b_2|
\right).
$$

如果只在单个训练块内估计，且假设块内信道不变，则通常不显式构造 $\mathbf R_{\mathrm{time}}$。

**时延相关矩阵**  
$\mathbf R_{\mathrm{delay}}$ 由 PDP 决定。若离散时延 tap 之间按 WSSUS 假设相互独立，则：

$$
\mathbf R_{\mathrm{delay}}
=
\operatorname{diag}(P_0,P_1,\dots,P_{L_h-1}),
$$

其中 $P_n=\mathbb E\|\mathbf H_b[n]\|_F^2$ 表示第 $n$ 个离散时延 tap 的平均功率，并通常归一化为 $\sum_n P_n=1$。PDP 可来自标准模型、测量或仿真设定。

**空间相关矩阵**  
$\mathbf R_t$ 和 $\mathbf R_r$ 可由角度功率谱、指数相关模型、样本协方差或 CDL/几何路径参数生成。以 ULA 为例，若归一化角度为 $\Omega=(d/\lambda)\sin\theta$，角度功率谱为 $p(\Omega)$，则空间相关矩阵元素可写为：

$$
[\mathbf R]_{p,q}
=
\int p(\Omega)e^{-j2\pi(p-q)\Omega}\,d\Omega.
$$

若采用离散路径/射线模型，则：

$$
[\mathbf R]_{p,q}
=
\frac{\sum_{\ell}P_\ell e^{-j2\pi(p-q)\Omega_\ell}}
{\sum_{\ell}P_\ell}.
$$

也可以用指数相关近似：

$$
[\mathbf R]_{p,q}=\rho^{|p-q|},
\qquad 0\le|\rho|<1.
$$

在实测或链路仿真中，也可以由多个信道快照估计样本协方差：

$$
\hat{\mathbf R}
=
\frac{1}{N_s}\sum_{s=1}^{N_s}
\mathbf h_s\mathbf h_s^H.
$$

估计误差协方差为：

$$
\boxed{
\mathbf R_e
=
\mathbf R_h
-\mathbf R_h\mathbf A^H
(\mathbf A\mathbf R_h\mathbf A^H+\mathbf R_w)^{-1}
\mathbf A\mathbf R_h
}
$$

得到 $\hat{\mathbf h}_{\mathrm{LMMSE}}$ 后，同样先恢复 $\hat{\mathbf H}_b[n]$，再沿离散时延维度做 DFT 生成 $\hat{\mathbf H}_b[k]$。

**算法伪代码**

```text
Input:
  Pilot observation vector y
  Joint observation matrix A
  PDP or delay covariance R_delay
  Spatial covariances Rt and Rr
  Optional time covariance R_time
  Noise covariance Rw

Steps:
  1. Build Rh from the channel model:
       Rh = R_time kron R_delay kron Rt^T kron Rr
  2. Compute Cyy = A Rh A^H + Rw
  3. Compute h_hat = Rh A^H Cyy^(-1) y
  4. Reshape h_hat into H_hat_b[n]
  5. Generate H_hat_b[k] = sum_n H_hat_b[n] exp(-j 2 pi k n / N)

Output:
  H_hat_b[k] on the required subcarriers in the block
```

LMMSE 是工程上最常见的“基于信道模型”的联合估计方法。它显式利用 PDP、空间相关和时间相关，因此比逐子载波 LS 更符合宽带 MIMO-OFDM 信道结构。

## 2.5 OMP 角度-时延稀疏信道估计

对第 $b$ 个训练块，信道建模文档中的角度域表示为：

$$
\mathbf H_b^a[k]
=
\mathbf U_r^H\mathbf H_b[k]\mathbf U_t.
$$

再沿子载波维度做 IDFT，可得到角度-时延域信道：

$$
\mathbf H_b^a[n]
=
\mathbf U_r^H\mathbf H_b[n]\mathbf U_t.
$$

在大规模 MIMO、毫米波或路径数较少的场景中，$\mathbf H_b^a[n]$ 在三维网格：

$$
(\text{AoA bin},\text{AoD bin},\text{delay tap})
$$

上近似稀疏。令：

$$
\mathbf s
=
\operatorname{vec}\{\mathbf H_b^a[n]\}
$$

为角度-时延稀疏向量。由：

$$
\mathbf H_b[n]
=
\mathbf U_r\mathbf H_b^a[n]\mathbf U_t^H
$$

和沿离散时延维度做 DFT，可得到：

$$
\mathbf H_b[k_i]
=
\sum_{n=0}^{L_h-1}
\mathbf U_r\mathbf H_b^a[n]\mathbf U_t^H
e^{-j2\pi k_i n/N}.
$$

对第 $i$ 个导频 RE，有：

$$
\mathbf y_i
=
\mathbf H_b[k_i]\mathbf x_i+\mathbf w_i.
$$

利用恒等式 $\mathbf A\mathbf B\mathbf C=(\mathbf C^T\otimes\mathbf A)\operatorname{vec}(\mathbf B)$，得到：

$$
\mathbf y_i
=
\mathbf \Phi_i\mathbf s+\mathbf w_i.
$$

若定义：

$$
\mathbf s
=
\begin{bmatrix}
\operatorname{vec}(\mathbf H_b^a[0])\\
\operatorname{vec}(\mathbf H_b^a[1])\\
\vdots\\
\operatorname{vec}(\mathbf H_b^a[L_h-1])
\end{bmatrix}
\in\mathbb C^{N_rN_tL_h\times 1},
$$

则第 $i$ 个导频 RE 的测量矩阵为：

$$
\boxed{
\mathbf \Phi_i
=
\begin{bmatrix}
e^{-j2\pi k_i0/N}
(\mathbf x_i^T\mathbf U_t^*\otimes\mathbf U_r)
&
e^{-j2\pi k_i1/N}
(\mathbf x_i^T\mathbf U_t^*\otimes\mathbf U_r)
&
\cdots
&
e^{-j2\pi k_i(L_h-1)/N}
(\mathbf x_i^T\mathbf U_t^*\otimes\mathbf U_r)
\end{bmatrix}
}
$$

其中：

- $\mathbf \Phi_i\in\mathbb C^{N_r\times N_rN_tL_h}$；
- $\mathbf x_i^T\mathbf U_t^*\in\mathbb C^{1\times N_t}$ 表示发射导频向量在发射角度字典上的等效权重；
- $\mathbf U_r\in\mathbb C^{N_r\times N_r}$ 是接收角度字典；
- 每个横向 block 对应一个离散时延 tap。

堆叠所有导频：

$$
\boxed{
\mathbf y=\mathbf \Phi\mathbf s+\mathbf w
}
$$

其中：

$$
\mathbf \Phi
=
\begin{bmatrix}
\mathbf \Phi_1\\
\mathbf \Phi_2\\
\vdots\\
\mathbf \Phi_{N_{\mathrm{obs}}}
\end{bmatrix}
\in
\mathbb C^{N_rN_{\mathrm{obs}}\times N_rN_tL_h}.
$$

因此 $\mathbf \Phi$ 同时包含：

- 发射导频 $\mathbf x_i$；
- 接收/发射角度字典 $\mathbf U_r,\mathbf U_t$；
- 子载波到离散时延 tap 的相位 $e^{-j2\pi k_in/N}$；
- 可能存在的接收合并或其它前端线性处理矩阵。

OMP 求解：

$$
\min_{\mathbf s}\|\mathbf s\|_0,
\qquad
\text{s.t.}
\quad
\|\mathbf y-\mathbf \Phi\mathbf s\|_2^2\le \epsilon.
$$

每次迭代选择与残差最相关的角度-时延原子：

$$
j_i
=
\arg\max_j
\left|
\boldsymbol\phi_j^H\mathbf r^{(i-1)}
\right|.
$$

在已选支撑集 $\mathcal S^{(i)}$ 上做 LS：

$$
\boxed{
\hat{\mathbf s}_{\mathcal S^{(i)}}
=
(\mathbf \Phi_{\mathcal S^{(i)}}^H
\mathbf \Phi_{\mathcal S^{(i)}})^{-1}
\mathbf \Phi_{\mathcal S^{(i)}}^H\mathbf y
}
$$

残差更新：

$$
\mathbf r^{(i)}
=
\mathbf y
-\mathbf \Phi_{\mathcal S^{(i)}}
\hat{\mathbf s}_{\mathcal S^{(i)}}.
$$

恢复 $\hat{\mathbf s}$ 后，先得到角度-时延域 $\hat{\mathbf H}_b^a[n]$，再回到天线-时延域：

$$
\hat{\mathbf H}_b[n]
=
\mathbf U_r\hat{\mathbf H}_b^a[n]\mathbf U_t^H,
$$

最后生成 CFR：

$$
\hat{\mathbf H}_b[k]
=
\sum_n
\hat{\mathbf H}_b[n]e^{-j2\pi kn/N}.
$$

**算法伪代码**

```text
Input:
  Pilot observations y_i on positions (q_i, k_i) within one N_sy-symbol block
  Pilot vectors x_i
  Receive and transmit DFT dictionaries Ur, Ut
  Delay grid length Lh
  Maximum path count L or residual threshold epsilon

Steps:
  1. Build the angle-delay measurement matrix Phi:
       combine pilot vectors, Ur, Ut, and delay phase exp(-j 2 pi k_i n / N)
  2. Initialize residual r = y, support S = empty set
  3. For iter = 1 to L:
       c = Phi^H r
       j = argmax_j |c_j|
       S = S union {j}
       s_S = (Phi_S^H Phi_S)^(-1) Phi_S^H y
       r = y - Phi_S s_S
       if ||r||_2^2 <= epsilon, stop
  4. Put s_S back into the full sparse vector s_hat
  5. Reshape s_hat into H_a_hat_b[n]
  6. Transform back:
       H_hat_b[n] = Ur H_a_hat_b[n] Ut^H
       H_hat_b[k] = sum_n H_hat_b[n] exp(-j 2 pi k n / N)

Output:
  H_hat_b[k] on the required subcarriers in the block
```

OMP 的核心优势是它直接匹配角度-时延稀疏信道模型；缺点是依赖字典精度。若真实 AoA/AoD 或时延不落在网格上，会出现 off-grid 泄漏。

# 3 选择性衰落信道重构与插值

本章讨论时间、频率和角度选择性衰落信道下的插值。此时信道不是只随子载波变化，也不是在一个训练块内严格不变，而是同时具有：

- 频率选择性：由有限时延扩展或 PDP 决定；
- 时间选择性：由多普勒扩展决定，不同 OFDM 符号上的信道缓慢变化且相关；
- 空间/角度选择性：由 AoA/AoD 角度功率谱、阵列响应和空间相关决定。

因此，“插值”不应理解为对彼此独立的子载波做曲线拟合，而应理解为模型化重构：先得到导频 RE 上的 CFR 初估计，再根据有限时延、时间相关、空间/角度相关或角度-时延稀疏结构，把导频 RE 上的信息外推到数据 RE。核心思想是从导频 CFR 初估计中恢复共享的选择性衰落结构，再生成目标时频资源上的 $\hat{\mathbf H}[q,k]$。

## 3.1 从导频 RE 到数据 RE 的模型化重构

设一个插值窗口内的导频位置集合为 $\mathcal P$，数据位置集合为 $\mathcal D$。每个位置由 OFDM 符号和子载波二元组表示：

$$
(q,k)\in\mathcal P\cup\mathcal D.
$$

ZF/LS 首先给出导频位置 CFR 初估计：

$$
\{\hat{\mathbf H}_{p}[q,k]:(q,k)\in\mathcal P\}.
$$

选择性衰落信道插值/重构进一步把这些 CFR 初估计转化为某种模型参数：

$$
\hat{\boldsymbol\theta}
\in
\{\hat{\mathbf H}[q,n],\hat{\mathbf h},\hat{\mathbf H}^a[q,n],\hat{\mathbf s}\}.
$$

信道插值就是利用两级映射：

$$
\boxed{
\{\hat{\mathbf H}_{p}[q,k]:(q,k)\in\mathcal P\}
\longrightarrow
\hat{\boldsymbol\theta}
\longrightarrow
\hat{\mathbf H}[q,k],
\qquad
(q,k)\in\mathcal D
}
$$

不同方法的差别在于这个映射使用了多少选择性衰落信道模型信息：

- IDFT/DFT 重构：利用有限时延扩展；
- Wiener/LMMSE 预测：利用 PDP、多普勒和空间相关；
- 角度-时延/角度-时延-多普勒重构：利用路径在角度、时延和多普勒域的结构或稀疏性。

## 3.2 基于 IDFT/DFT 的频率选择性重构

IDFT/DFT 重构只处理频率选择性。它的核心步骤是：对某个固定 OFDM 符号 $q$，先把导频子载波上的 CFR 初估计转换到离散时延域，得到一个等效 CIR；再只保留有效长度 $L_h$ 内的离散时延 tap；最后通过 DFT 重构该符号上的全频域 CFR。

若已经有一个全频域粗 CFR 向量：

$$
\hat{\mathbf h}_{f,\mathrm{raw}}
=
\begin{bmatrix}
\hat H[q,0] & \hat H[q,1] & \cdots & \hat H[q,N-1]
\end{bmatrix}^T,
$$

则先做 IDFT：

$$
\hat{\mathbf h}_{\tau,\mathrm{raw}}
=
\mathbf F_N^H\hat{\mathbf h}_{f,\mathrm{raw}}.
$$

只保留前 $L_h$ 个有效离散时延 tap：

$$
\hat h_{\tau}[n]
=
\begin{cases}
\hat h_{\tau,\mathrm{raw}}[n], & 0\le n<L_h,\\
0, & L_h\le n<N.
\end{cases}
$$

再做 DFT 得到等效 CFR：

$$
\boxed{
\hat{\mathbf h}_{f}
=
\mathbf F_N\hat{\mathbf h}_{\tau}
}
$$

如果已经从导频 CFR 初估计中得到第 $q$ 个 OFDM 符号的等效离散时延 tap：

$$
\hat{\mathbf H}[q,n],
\qquad n=0,\dots,L_h-1,
$$

则任意子载波上的等效 CFR 可直接由 DFT 得到：

$$
\boxed{
\hat{\mathbf H}[q,k]
=
\sum_{n=0}^{L_h-1}
\hat{\mathbf H}[q,n]e^{-j2\pi kn/N}
}
$$

这个公式自动保证同一个 OFDM 符号内不同子载波共享同一组离散时延 tap。它是宽带 OFDM 中最基本的频域重构方式。若信道随时间变化，则需要对每个目标符号 $q$ 分别得到或预测 $\hat{\mathbf H}[q,n]$；若某个符号上的频域导频不足，则应结合 3.3 或 3.4 的时间相关模型联合恢复。

若系统先在导频子载波上做了逐 RE LS，得到粗估计 $\hat{\mathbf H}_{p}[q,k]$，则对每个有足够频域导频的符号 $q$，可以对每个天线对或联合向量执行 IDFT/DFT 重构。定义该符号上的导频子载波集合为：

$$
\mathcal P_q=\{k:(q,k)\in\mathcal P\}.
$$

若导频只覆盖部分子载波，等价做法是用部分 DFT 拟合长度为 $L_h$ 的 CIR：

$$
\hat{\mathbf h}_{\tau}
=
(\mathbf F_{\mathcal P_q,L_h}^H\mathbf F_{\mathcal P_q,L_h})^{-1}
\mathbf F_{\mathcal P_q,L_h}^H
\hat{\mathbf h}_{p}.
$$

再由 DFT 生成所有子载波上的等效 CFR：

$$
\hat{\mathbf h}_{f}
=
\mathbf F_{N,L_h}\hat{\mathbf h}_{\tau}.
$$

**算法伪代码**

```text
Input:
  Pilot-subcarrier estimates H_p_hat[q,k] on one OFDM symbol q
  Pilot subcarrier set P
  FFT size N
  Channel length Lh

Steps:
  1. Select pilot CFR estimates on the target symbol q
  2. Estimate an equivalent length-Lh CIR for symbol q:
       h_tau_hat = pinv(F_P_Lh) h_p_hat
  3. Apply DFT to reconstruct full-band equivalent CFR:
       h_f_hat = F_N_Lh h_tau_hat
  4. Reshape h_f_hat into H_hat[q,k] for all k

Output:
  Full-band H_hat[q,k] on the target symbol
```

## 3.3 基于时间选择性的符号级插值

时间选择性来自多普勒。对双选择性衰落信道，同一个训练窗口内不同 OFDM 符号上的信道通常不是严格相同，而是缓慢变化并具有时间相关性。因此时域插值应该在导频符号和数据符号之间进行，插值对象优先选择离散时延 tap、角度-时延系数或联合信道向量，而不是每个子载波独立插值。

若第 $q_1$ 和第 $q_2$ 个导频符号的离散时延 tap 估计为：

$$
\hat{\mathbf H}[q_1,n],
\qquad
\hat{\mathbf H}[q_2,n].
$$

则对两个导频符号之间的数据符号 $q$，低复杂度做法是对每个离散时延 tap 做符号级线性插值：

$$
\boxed{
\hat{\mathbf H}[q,n]
=
\frac{q_2-q}{q_2-q_1}\hat{\mathbf H}[q_1,n]
+
\frac{q-q_1}{q_2-q_1}\hat{\mathbf H}[q_2,n]
}
$$

然后再沿离散时延维度做 DFT 生成 $\hat{\mathbf H}[q,k]$。若使用多普勒统计，可用 Jakes/Clarke 相关函数。两个 OFDM 符号间隔为 $\Delta q$ 时，近似时间间隔为 $\Delta qT_{\mathrm{sym}}$：

$$
R_{\mathrm{time}}(\Delta q)
=
J_0(2\pi f_D T_{\mathrm{sym}}\Delta q).
$$

导频符号间隔仍需满足相干时间约束：

$$
S_tT_{\mathrm{sym}}\lesssim \frac{1}{2f_D}.
$$

## 3.4 时频空 Wiener 预测

Wiener 预测是 LMMSE 在数据 RE 上的预测形式，也是时间、频率、空间/角度选择性衰落信道中最统一的插值表达。令 $\mathbf h_p$ 表示导频位置的联合信道向量，$\mathbf h_d$ 表示数据位置的联合信道向量。导频位置估计为：

$$
\hat{\mathbf h}_p=\mathbf h_p+\mathbf e_p.
$$

线性预测器：

$$
\hat{\mathbf h}_d=\mathbf B\hat{\mathbf h}_p.
$$

最小化：

$$
\mathbb E\|\mathbf h_d-\mathbf B\hat{\mathbf h}_p\|_2^2.
$$

由正交性原理：

$$
\boxed{
\hat{\mathbf h}_d
=
\mathbf R_{d,p}
(\mathbf R_{p,p}+\mathbf R_{e,p})^{-1}
\hat{\mathbf h}_p
}
$$

其中：

- $\mathbf R_{p,p}$：导频 RE 之间的时间、频率、空间/角度协方差；
- $\mathbf R_{d,p}$：数据 RE 与导频 RE 之间的时间、频率、空间/角度交叉协方差；
- $\mathbf R_{e,p}$：导频估计误差协方差。

这些协方差来自信道建模。WSSUS + Kronecker 模型下可写成：

$$
\mathbf R
=
\mathbf R_{\mathrm{time}}
\otimes
\mathbf R_{\mathrm{freq}}
\otimes
\mathbf R_{\mathrm{tx}}^T
\otimes
\mathbf R_{\mathrm{rx}}.
$$

频域相关由 PDP 决定：

$$
\boxed{
R_f(\Delta k)
=
\sum_{n=0}^{L_h-1}
P_n e^{-j2\pi \Delta k n/N}
}
$$

时间相关可由多普勒模型决定。若以 OFDM 符号为时间采样单位，则：

$$
R_{\mathrm{time}}(\Delta q)
=
J_0(2\pi f_DT_{\mathrm{sym}}\Delta q).
$$

空间相关 $\mathbf R_{\mathrm{tx}},\mathbf R_{\mathrm{rx}}$ 可由角度功率谱、指数相关模型或几何路径生成。角度选择性衰落对插值的影响体现在空间/角度维协方差或角度字典中：若不同天线端口或波束方向看到的信道相关性不同，Wiener 预测需要利用 $\mathbf R_{\mathrm{tx}},\mathbf R_{\mathrm{rx}}$ 或角度域表示把导频端口上的信息映射到目标端口/目标波束上的信道。

**算法伪代码**

```text
Input:
  Pilot-position channel estimates h_p_hat
  Pilot positions P and data positions D
  PDP P_n
  Doppler model or time covariance
  Spatial covariances Rtx and Rrx
  Pilot estimation error covariance Re_p

Steps:
  1. Build symbol-level time correlation R_time from Doppler
  2. Build frequency correlation R_freq from PDP
  3. Build spatial correlation Rtx^T kron Rrx
  4. Extract R_pp and R_dp for pilot and data RE positions
  5. Compute h_d_hat = R_dp (R_pp + Re_p)^(-1) h_p_hat
  6. Map h_d_hat back to H_hat[q,k] on data REs

Output:
  H_hat[q,k] on data positions
```

## 3.5 角度-时延域重构与去噪

若估计算法输出角度-时延稀疏系数 $\hat{\mathbf s}$，则数据 RE 的信道由物理结构直接重构。对时间选择性信道，路径增益或角度-时延系数应允许随 OFDM 符号 $q$ 缓慢变化。先恢复：

$$
\hat{\mathbf H}^a[q,n],
$$

再回到天线-时延域：

$$
\hat{\mathbf H}[q,n]
=
\mathbf U_r\hat{\mathbf H}^a[q,n]\mathbf U_t^H.
$$

最后沿离散时延维度做 DFT：

$$
\boxed{
\hat{\mathbf H}[q,k]
=
\sum_{n=0}^{L_h-1}
\mathbf U_r\hat{\mathbf H}^a[q,n]\mathbf U_t^H
e^{-j2\pi kn/N}
}
$$

这类方法本质上不是普通插值，而是“路径参数重构”：只要角度、时延和随时间变化的路径增益估计准确，就可以在目标 OFDM 符号和所有子载波上生成一致的 MIMO 信道。若进一步采用角度-时延-多普勒字典，则稀疏系数可以直接描述不同多普勒 bin 上的路径变化。

**算法伪代码**

```text
Input:
  Estimated sparse vector s_hat
  Angle dictionaries Ur, Ut
  Delay grid length Lh
  Target OFDM symbol q
  FFT size N

Steps:
  1. Reshape or predict s_hat into angle-delay matrices H_a_hat[q,n]
  2. For each discrete-delay tap n:
       H_hat[q,n] = Ur H_a_hat[q,n] Ut^H
  3. For each desired subcarrier k:
       H_hat[q,k] = sum_n H_hat[q,n] exp(-j 2 pi k n / N)

Output:
  Model-consistent H_hat[q,k] on all required REs
```

角度-时延重构适合大规模 MIMO 和毫米波系统。若角度或时延网格过粗，需要过采样字典、off-grid 修正或连续参数估计。

## 3.6 块衰落信道下的简化

如果系统满足块衰落近似，即在连续 $N_{\mathrm{sy}}$ 个 OFDM 符号内：

$$
\mathbf H[0,k]\approx \mathbf H[1,k]\approx \cdots \approx \mathbf H[N_{\mathrm{sy}}-1,k],
$$

则时间选择性可以在该块内忽略。此时第 3 章的选择性衰落插值会简化为：

- 时间维：不需要对块内不同符号做时域插值，可先合并同一子载波、不同导频符号上的 CFR 初估计以提高 SNR；
- 频率维：仍需根据有限时延扩展或 PDP 做 IDFT/DFT、LMMSE 或 Wiener 频域重构；
- 空间/角度维：仍可利用空间相关、角度功率谱或角度-时延稀疏性做多天线联合估计。

若相邻训练块之间还需要预测，则可把每个训练块视为一个时间采样点，对块索引做低速率时间插值。此时 Jakes 相关函数中的时间间隔应使用训练块中心间隔，例如 $\Delta t=\Delta bN_{\mathrm{sy}}T_{\mathrm{sym}}$，而不是符号间隔 $\Delta qT_{\mathrm{sym}}$。

# 4 附录

## A.1 符号说明

| 符号 | 含义 |
|---|---|
| $N_t$ | 发射天线数 |
| $N_r$ | 接收天线数 |
| $f_c$ | 载波频率 |
| $B$ | 系统带宽 |
| $\Delta f$ | 子载波间隔 |
| $N_{\mathrm{FFT}}$ | FFT 点数，正文中常简记为 $N$ |
| $N_{\mathrm{act}}$ | 有效子载波数 |
| $T_{\mathrm{cp}}$ | 循环前缀长度 |
| $T_u$ | 不含 CP 的有效 OFDM 符号时长，$T_u=1/\Delta f$ |
| $T_{\mathrm{sym}}$ | 含 CP 的 OFDM 符号时长 |
| $N_p$ | 导频接收信号数量 |
| $N_{\mathrm{obs}}$ | 一个训练块内的导频 RE 数量，不包含接收天线维度 |
| $N_h$ | 宽带信道参数向量 $\mathbf h$ 的维度 |
| $\mathbf X$ | 局部导频矩阵 |
| $\mathbf Y$ | 局部接收导频矩阵 |
| $N_{\mathrm{sy}}$ | 一个训练窗口或块内包含的连续 OFDM 符号数 |
| $b$ | 训练块或相干块索引，主要用于块衰落简化模型 |
| $q$ | OFDM 符号索引 |
| $\mathbf H[q,k]$ | 第 $q$ 个 OFDM 符号、第 $k$ 个子载波上的 MIMO CFR |
| $\mathbf H[q,n]$ | 第 $q$ 个 OFDM 符号、第 $n$ 个离散时延 tap 上的 MIMO CIR |
| $\mathbf H_b[k]$ | 块衰落简化下，第 $b$ 个训练块、第 $k$ 个子载波上的 MIMO CFR |
| $\mathbf H_b[n]$ | 块衰落简化下，第 $b$ 个训练块、第 $n$ 个离散时延 tap 上的 MIMO CIR |
| $\mathbf A$ | 宽带联合观测矩阵 |
| $\mathbf h$ | 向量化宽带信道参数，可表示离散时延域（tap 域）信道或联合时频空信道 |
| $\mathbf R_h$ | 信道协方差矩阵 |
| $\mathbf R_w$ | 噪声协方差矩阵 |
| $\mathbf R_t$ | 发射端空间相关矩阵 |
| $\mathbf R_r$ | 接收端空间相关矩阵 |
| $\mathbf \Phi$ | 稀疏测量矩阵 |
| $\mathbf s$ | 角度-时延域稀疏向量 |
| $\mathbf U_t,\mathbf U_r$ | 发射端和接收端 DFT 角度字典 |
| $v$ | 移动速度 |
| $f_D$ | 最大多普勒频移 |
| $\tau_{\max}$ | 最大多径时延 |
| $L_h$ | 最大离散时延 tap 数 |
| $\sigma_w^2$ | 噪声功率 |

## A.2 方法对比

| 方法 | 估计对象 | 是否利用宽带模型 | 主要优点 | 主要缺点 | 适用场景 |
|---|---|---|---|---|---|
| ZF/LS | 导频 RE 上的 CFR 初估计 | 基础 LS 不利用；后续可接 IDFT/DFT、Wiener/LMMSE 或 OMP 重构 | 简单、贴近实际接收机 | 初估计抗噪声差，导频稀疏时需后续重构/插值 | 导频初估计 |
| MMSE | 宽带信道参数后验均值 | 是，通过完整先验分布 | MSE 最优 | 先验和积分通常难获得 | 理论分析、贝叶斯估计 |
| LMMSE | CIR/时频空联合信道向量 | 是，通过 PDP、空间相关和多普勒协方差 | 抗噪声、能利用相关性 | 协方差估计和矩阵求逆复杂 | 工程信道估计、Wiener 插值 |
| OMP | 角度-时延稀疏系数 | 是，通过角度字典和 delay DFT | 导频开销低、匹配稀疏物理路径 | 依赖稀疏性和网格精度 | 大规模 MIMO、毫米波 |

## A.3 参考文献

1. D. Tse and P. Viswanath, *Fundamentals of Wireless Communication*, Cambridge University Press, 2005. https://web.stanford.edu/~dntse/wireless_book.html
2. A. M. Sayeed, "Deconstructing Multiantenna Fading Channels," *IEEE Transactions on Signal Processing*, vol. 50, no. 10, pp. 2563-2579, 2002. https://minds.wisconsin.edu/handle/1793/9386
3. J. Jo and I. Sohn, "On the optimality of training signals for MMSE channel estimation in MIMO-OFDM systems," *EURASIP Journal on Wireless Communications and Networking*, 2015. https://link.springer.com/article/10.1186/s13638-015-0345-y
4. J. P. Nair and R. V. Raja Kumar, "Optimal Superimposed Training Sequences for Channel Estimation in MIMO-OFDM Systems," *EURASIP Journal on Advances in Signal Processing*, 2010. https://link.springer.com/article/10.1155/2010/140506
5. D. Katselis et al., "Training sequence design for MIMO channels: an application-oriented approach," *EURASIP Journal on Wireless Communications and Networking*, 2013. https://link.springer.com/article/10.1186/1687-1499-2013-245
6. G. Li and G. Liao, "A Pilot-Pattern Based Algorithm for MIMO-OFDM Channel Estimation," *Algorithms*, 2017. https://www.mdpi.com/1999-4893/10/1/3
7. MathWorks, "Estimate the wireless channel in a MIMO OFDM system using pilot signals." https://www.mathworks.com/help/comm/ref/ofdmchannelestimate.html
