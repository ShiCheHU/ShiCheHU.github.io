# MIMO OFDM 信道建模

> **导频设计与信道估计** 部分已拆分为独立文档：[MIMOOFDMChannelEstimation.md](MIMOOFDMChannelEstimation.md)。

## 目录

- [1 概述](#1-概述)
- [2 MIMO-OFDM 系统模型](#2-mimo-ofdm-系统模型)
  - [2.1 参数与符号](#21-参数与符号)
  - [2.2 连续域 MIMO CIR](#22-连续域-mimo-cir)
    - [2.2.1 注意：频域矩阵乘法 vs 时域卷积](#221-注意频域矩阵乘法-vs-时域卷积)
  - [2.3 离散域 MIMO CIR/CFR](#23-离散域-mimo-circfr)
  - [2.4 OFDM 频域输入输出模型](#24-ofdm-频域输入输出模型)
- [3 MIMO 信道建模方法](#3-mimo-信道建模方法)
  - [3.1 独立 Rayleigh/Rician tap 模型](#31-独立-rayleighrician-tap-模型)
  - [3.2 空间相关 Kronecker 模型](#32-空间相关-kronecker-模型)
    - [3.2.1 Kronecker 模型的基本假设](#321-kronecker-模型的基本假设)
    - [3.2.2 Kronecker 模型下的离散域 MIMO CIR](#322-kronecker-模型下的离散域-mimo-cir)
    - [3.2.3 Kronecker 模型下的 MIMO CFR](#323-kronecker-模型下的-mimo-cfr)
    - [3.2.4 空间相关矩阵的构造与推导](#324-空间相关矩阵的构造与推导)
      - [3.2.4.1 从角度功率谱推导相关矩阵](#3241-从角度功率谱推导相关矩阵)
      - [3.2.4.1b 角度扩展（Angular Spread）的生成](#3241b-角度扩展angular-spread的生成)
      - [3.2.4.2 指数相关模型](#3242-指数相关模型)
      - [3.2.4.3 利用阵列导向矢量构造相关矩阵](#3243-利用阵列导向矢量构造相关矩阵)
    - [3.2.5 模型优缺点](#325-模型优缺点)
  - [3.3 角度域信道建模与稀疏表示](#33-角度域信道建模与稀疏表示)
    - [3.3.1 连续角度域——从阵列相位差到导向矢量](#331-连续角度域从阵列相位差到导向矢量)
    - [3.3.2 离散角度域——DFT 基向量与虚拟信道表示](#332-离散角度域dft-基向量与虚拟信道表示)
    - [3.3.3 角度-时延域联合稀疏性](#333-角度-时延域联合稀疏性)
    - [3.3.4 扩展：UPA（均匀平面阵）的角度域表示](#334-扩展upa均匀平面阵的角度域表示)
  - [3.4 3GPP TR 38.901 CDL/TDL 模型](#34-3gpp-tr-38901-cdltdl-模型)
- [4 仿真建模流程](#4-仿真建模流程)
- [5 常见模型选择建议](#5-常见模型选择建议)
- [6 References](#6-references)

# 1 概述

MIMO（Multiple-Input Multiple-Output）是多天线输入输出技术。利用发送端和接收端部署的多根天线，提供空域传输的自由度，提升无线通信系统的容量、频谱效率、覆盖范围。尤其的，在高频段，MIMO技术尤为重要。因为毫米波、太赫兹频段覆盖距离短，MIMO技术可以提升覆盖。

MIMO主要应用空间分集、空分复用、波束赋形三种技术实现自身优势。其中，空间分集、空分复用可以归并到空时块编码技术中。这3种技术介绍如下：

- 空间分集：将同一组数据在不同的信道发送。这样，同一组数据有了多个副本，接收端可以对多个信道的数据合并，等效提示信噪比。
- 空分复用：将同一组数据分段，在不同信道发送不同数据部分。这样，即使部分数据经过的信道质量差，但其他的数据传输质量好。等效降低了信道的影响。
- 波束赋形：通过调整天线阵列，让多根天线发送的信号叠加后形成有向信号。从而提升给定方向信号的覆盖范围。在超大规模MIMO系统中得到应用。

MIMO技术的应用在移动蜂窝系统、Wi-Fi技术中：

- 在5G NR中，混合波束赋形用在FR 2频段，数字波束赋形用在sub-6G 频段。
- 在Wi-Fi系统中，主要采用空时块编码。

MIMO在无线技术演进是贯穿5G标准化，天线数增加、不同天线面板、不同传输点之间的协作是发展趋势。到了6G，MIMO依然是讨论的主要方向，包括：信道的设计、波束管理、CSI获取、上行覆盖、初始接入、移动性管理等。MIMO技术等演进面临如下挑战：

- 天线数增加导致的信道获取开销过高。上下行资源中不得不分配更多的资源给导频，导频资源和数据是相互正交的，这使得用于获取信道的导频资源过高。
- 硬件成本提升。要支持更多的数据流、更多的用户数，要求配置的射频链路数量也线性增加，而这导致成本和功耗低增加。
- 数据处理的复杂度提升。天线数增加，采用混合波束赋形，需要波束管理，避免波束不对齐。数据流增加，接收端做信道估计和实时数据接收的复杂度提升。
- 多用户干扰消除和覆盖。多用户在5G标准讨论很多，也在标准中得到支持，但在实际中没有得到应用。因为多用户完全可以通过空分、时分和频分复用来支持。支持多用户预编码成本很高，而且只能用在移动慢（移动快实时更新预编码成本高）、距离近（距离远空间上就分开了）的设备。在6G阶段，很多公司提的更多用户的支持就很难支持。

这片文章主要是介绍MIMO OFDM信道的建模。之前已经介绍过SISO OFDM信道建模，可扩展到MIMO OFDM是完全不同的。相比SISO OFDM系统的信道建模，MIMO OFDM信道建模有如下新的特点：

- SISO OFDM系统主要考虑频率和时间选择性，MIMO OFDM系统还需要额外考虑空间选择性，通常用AoA/AoD进行建模。
- MIMO OFDM系统需要考虑天线间耦合的影响的建模。例如，天线间距的影响。
- 需要考虑信道矩阵的秩和条件数的变化，做秩自适应。
- 导频设计也需要考虑不同天线端口之间的正交。

在 SISO OFDM 中，信道估计的核心是恢复每个 OFDM 符号、每个子载波上的标量信道频率响应：

$$
H[m,k].
$$

在 MIMO-OFDM 中，标量信道扩展为矩阵信道：

$$
\boxed{
\mathbf H[m,k]\in \mathbb C^{N_r\times N_t}
}
$$

其中 $N_t$ 为发射天线数，$N_r$ 为接收天线数。矩阵元素 $H_{r,t}[m,k]$ 表示第 $t$ 根发射天线到第 $r$ 根接收天线的 SISO 子信道。

因此 MIMO-OFDM 信道建模需要同时描述四类相关性：

| 维度 | 物理来源 | 数学对象 |
|------|----------|----------|
| 时延相关 | 多径时延扩展 | PDP、CIR tap |
| 频率相关 | 有限时延扩展 | CFR 子载波相关 |
| 时间相关 | 多普勒扩展 | Jakes/AR/BEM 过程 |
| 空间相关 | 多天线阵列、AoA/AoD、角度扩展 | 阵列响应、空间协方差、角度域稀疏性 |

MIMO 相比 SISO 的关键新增问题是：不同发射天线和接收天线之间的子信道一般不是独立的。空间相关性会影响：

- 信道矩阵秩与条件数；
- 空间复用层数；
- 分集增益；
- 预编码和检测性能；
- 信道估计所需导频开销。

# 2 MIMO-OFDM 系统模型

## 2.1 参数与符号

记：

- 发射天线数：$N_t$
- 接收天线数：$N_r$
- OFDM 子载波数：$N$
- OFDM 符号数：$M$
- 子载波间隔：$\Delta f$
- 有效 OFDM 符号时长：$T_u=1/\Delta f$
- 循环前缀长度：$T_{cp}$
- 总 OFDM 符号时长：$T_{\mathrm{sym}}=T_u+T_{cp}$
- 采样间隔：$T_s=1/(N\Delta f)$
- 路径数或离散 tap 数：$L$
- 第 $m$ 个 OFDM 符号对应时刻：$t_m=mT_{\mathrm{sym}}$

通常要求：

$$
\boxed{
T_{cp}>\tau_{\max}
}
$$

这样去 CP 后，线性卷积可等效为循环卷积，每个子载波上的 MIMO-OFDM 信道可近似解耦为窄带 MIMO 平坦衰落信道。

## 2.2 连续域 MIMO CIR

第 $t$ 根发射天线到第 $r$ 根接收天线的连续域时变冲激响应为：

$$
h_{r,t}(t,\tau)=\sum_{\ell=0}^{L-1}\alpha_{\ell}^{(r,t)}(t)\delta(\tau-\tau_\ell).
$$

堆叠所有天线对，可得矩阵形式：

$$
\boxed{
\mathbf H(t,\tau)=\sum_{\ell=0}^{L-1}\mathbf G_\ell(t)\delta(\tau-\tau_\ell)
}
$$

其中：

$$
\mathbf G_\ell(t)=
\begin{bmatrix}
\alpha_\ell^{(1,1)}(t) & \cdots & \alpha_\ell^{(1,N_t)}(t)\\
\vdots & \ddots & \vdots\\
\alpha_\ell^{(N_r,1)}(t) & \cdots & \alpha_\ell^{(N_r,N_t)}(t)
\end{bmatrix}
\in \mathbb C^{N_r\times N_t}.
$$

若第 $\ell$ 条路径的功率为 $P_\ell$，也常写为：

$$
\mathbf H(t,\tau)=\sum_{\ell=0}^{L-1}\sqrt{P_\ell}\mathbf G_\ell(t)\delta(\tau-\tau_\ell),
\qquad \sum_{\ell=0}^{L-1}P_\ell=1.
$$

这里 $\mathbf G_\ell(t)$ 描述空间衰落和时间选择性，$P_\ell$ 描述功率时延谱。

### 2.2.1 信道的时频域关系

时域中MIMO 信道的输入-输出是**矩阵卷积**关系，频域MIMO的输入输出时是矩阵乘法。下面从 LTI 系统角度给出一个严格的推导。

**时域卷积形式**：对于 $N_t$ 输入 $N_r$ 输出的 MIMO LTI 系统，信道是矩阵值冲激响应 $\mathbf H(t, \tau)$。时域输入 $\mathbf x(t) \in \mathbb{C}^{N_t}$ 和输出 $\mathbf y(t) \in \mathbb{C}^{N_r}$ 的关系为矩阵卷积：

$
\boxed{
\mathbf y(t) = \int_{-\infty}^{\infty} \mathbf H(t, \tau)\,\mathbf x(t - \tau)\,d\tau
}
$

将上式写成分量形式：

$
y_r(t) = \sum_{t=1}^{N_t} \int_{-\infty}^{\infty} h_{r,t}(t, \tau)\,x_t(t - \tau)\,d\tau,
\quad r = 1, \dots, N_r.
$

**关键观察**：在分量形式下，每一个 $(r, t)$ 对都是标量值函数的普通卷积。内层积分与有限天线求和可交换，因此可以对每个分量**分别**做傅里叶变换：

$
\begin{aligned}
Y_r(f) &= \sum_{t=1}^{N_t} \mathcal{F}\!\left\{ \int h_{r,t}(t, \tau)\,x_t(t - \tau)\,d\tau \right\} \\[4pt]
&= \sum_{t=1}^{N_t} H_{r,t}(f)\,X_t(f).
\end{aligned}
$

将所有 $r$ 和 $t$ 的标量关系按矩阵乘法重组，即得频域矩阵形式：

$
\boxed{
\mathbf Y(f) = \mathbf H(f)\,\mathbf X(f)
}
$

**为什么不能直接对向量值函数做傅里叶变换？**

从数学上讲，标量值函数（绝对可积/平方可积/施瓦茨类）构成无穷维向量空间，傅里叶变换是该空间上的线性算子。但值域为 $\mathbb{C}^{N_r}$ 的向量值函数所构成的集合并非向量空间（它是函数环 $\mathcal{A}$ 上的模，$\mathcal{A}$-module），向量空间中的正交分解理论不能直接套用到模上。因此"对向量值函数逐分量做 FT，然后把各分量拼成向量"这一操作依赖于 MIMO **卷积系统**的结构（即各 SISO 子信道可交换 FT 与求和），而非 FT 的一般性质。

**MIMO-OFDM 中的结论**：在 OFDM 系统中，每个子载波上的信道是窄带平坦的，因此频域矩阵乘法 $\mathbf y[m,k] = \mathbf H[m,k]\,\mathbf x[m,k]$ 在每个时频资源上成立。这是后续一切 MIMO 信道建模和估计算法的基础。

## 2.3 离散域 MIMO CIR/CFR

在整数采样 tap 模型下，离散 MIMO CIR 为：

$$
\boxed{
\mathbf H[m,n]\in \mathbb C^{N_r\times N_t},\qquad n=0,1,\dots,N-1
}
$$

其中只有有限个 tap 非零。第 $m$ 个 OFDM 符号、第 $k$ 个子载波上的 MIMO CFR 为离散 CIR 沿时延 tap 维度的 DFT。由于 $\mathbf H[m,n]$ 是矩阵，严格写法应先固定收发天线对 $(r,t)$：

$$
\boxed{
H_{r,t}[m,k]
=
\sum_{n=0}^{N-1}
H_{r,t}[m,n]\ e^{-j\frac{2\pi}{N}kn},
\qquad r=1,\dots,N_r,\ t=1,\dots,N_t.
}
$$

把所有 $(r,t)$ 的结果重新排列成矩阵，才得到 $\mathbf H[m,k]\in\mathbb C^{N_r\times N_t}$。矩阵紧凑写法

$$
\mathbf H[m,k]
=
\sum_{n=0}^{N-1}\mathbf H[m,n]e^{-j\frac{2\pi}{N}kn}
$$

表示对矩阵 $\mathbf H[m,n]$ 的每个元素乘以同一个 DFT 相位标量，而不是矩阵与 DFT 向量相乘。

逆变换同样逐天线对定义：

$$
\boxed{
H_{r,t}[m,n]
=
\frac{1}{N}\sum_{k=0}^{N-1}
H_{r,t}[m,k]\ e^{j\frac{2\pi}{N}kn}
}
$$

若采用连续时延路径模型，则对每个天线对：

$$
\boxed{
H_{r,t}[m,k]
=
\sum_{\ell=0}^{L-1}
\sqrt{P_\ell}\,[\mathbf G_\ell[m]]_{r,t}\ e^{-j2\pi k\Delta f\tau_\ell}
}
$$

该式是 SISO CFR 对所有天线对的并行应用。若写成矩阵形式，应理解为逐元素叠加路径相位。

## 2.4 OFDM 频域输入输出模型

第 $m$ 个 OFDM 符号、第 $k$ 个子载波上，发射向量为：

$$
\mathbf x[m,k]\in \mathbb C^{N_t\times 1}.
$$

接收向量为：

$$
\mathbf y[m,k]\in \mathbb C^{N_r\times 1}.
$$

若 CP 足够长且一个 OFDM 符号内信道近似不变，则频域模型为：

$$
\boxed{
\mathbf y[m,k]=\mathbf H[m,k]\mathbf x[m,k]+\mathbf w[m,k]
}
$$

其中 $\mathbf w[m,k]\sim\mathcal{CN}(\mathbf 0,\sigma_w^2\mathbf I_{N_r})$。

展开后：

$$
y_r[m,k]=\sum_{t=1}^{N_t}H_{r,t}[m,k]x_t[m,k]+w_r[m,k].
$$

如果一个子载波上同时发送多个空间流，则 $\mathbf H[m,k]$ 还会进入 MIMO 检测或均衡：

$$
\hat{\mathbf x}[m,k]=\mathbf W[m,k]\mathbf y[m,k].
$$

常见线性均衡器为：

$$
\mathbf W_{\mathrm{ZF}}[m,k]
=
\left(\mathbf H^H[m,k]\mathbf H[m,k]\right)^{-1}\mathbf H^H[m,k],
$$

$$
\mathbf W_{\mathrm{MMSE}}[m,k]
=
\left(\mathbf H^H[m,k]\mathbf H[m,k]+\sigma_w^2\mathbf I\right)^{-1}\mathbf H^H[m,k].
$$

因此信道估计误差不仅影响相位补偿，也会影响空间流分离。

# 3 MIMO 信道建模方法

## 3.1 独立 Rayleigh/Rician tap 模型

独立 tap 模型是 MIMO 信道建模中最基础的方法，其核心假设是：

- 不同时延 tap 之间统计独立（WSSUS 假设）；
- 同一 tap 内，不同天线对 $(r, t)$ 之间也统计独立——即无视阵列几何和空间相关。

该模型的思想可追溯到 Bello (1963) 的 WSSUS 信道理论框架，并在 COST 207 (1986) 等标准化工作中作为默认假设广泛使用，至今仍是算法验证阶段的首选基线。

### 3.1.1 多径信道的构造

系统先设定一条功率时延谱（PDP），给出一组时延和功率对：

$
\{(\tau_\ell, P_\ell)\}_{\ell=0}^{L-1}, \qquad \sum_{\ell=0}^{L-1} P_\ell = 1.
$

以采样间隔 $T_s$ 离散化为整数 tap 索引 $n_\ell = \mathrm{round}(\tau_\ell / T_s)$。

每条路径 $\ell$ 的 MIMO 空间信道矩阵为：

$
\boxed{
\mathbf G_\ell[m] \in \mathbb C^{N_r \times N_t}
}
$

其元素 $[\mathbf G_\ell[m]]_{r,t}$ 独立地服从零均值复高斯分布，方差为 1。第 $m$ 个 OFDM 符号、第 $n$ 个采样 tap 上的 MIMO CIR 为：

$
\boxed{
\mathbf H[m, n] = \sqrt{P_n}\ \mathbf G_n[m],
\qquad \mathbf G_n[m] = \mathbf G_{\sigma(n)}[m]
}
$

其中 $n$ 与 $\ell$ 通过 PDP 中的 tap 索引映射联系。若某个 tap 上叠加了多条路径，则其功率 $P_n$ 为对应路径功率之和。

### 3.1.2 LOS 径的概率化生成

物理上，首径（最小传播时延）最可能包含 LOS 分量；后续径经过多次反射和散射，LOS 的概率极小。因此工程上通常**仅对前几个 tap 按一定概率施加 Rician 衰落，其余 tap 一律按 Rayleigh 处理**。

具体做法（以首径 LOS、其余 NLOS 为例）：

- 生成一个伯努利随机变量 $b_{\mathrm{LOS}} \in \{0, 1\}$，$\Pr(b_{\mathrm{LOS}} = 1) = p_{\mathrm{LOS}}$。
- 若 $b_{\mathrm{LOS}} = 1$，首径（$\ell = 0$）采用 **Rician** 衰落，后续径采用 **Rayleigh**；
- 若 $b_{\mathrm{LOS}} = 0$，所有路径均采用 **Rayleigh**。

$p_{\mathrm{LOS}}$ 可以是固定值（如 20%），也可以按距离相关函数取值（见 3GPP TR 38.901 各场景定义）。

### 3.1.3 单条路径的 Rician 空间矩阵

对一条被标记为 LOS 的路径 $\ell$，其空间矩阵 $\mathbf G_\ell[m]$ 的第 $(r, t)$ 个元素按 Rician 分布生成：

$$
\boxed{
\left[\mathbf G_\ell[m]\right]_{r,t}
=
\sqrt{\frac{K_\ell}{K_\ell + 1}}\ \alpha_{\mathrm{LOS}}
+
\sqrt{\frac{1}{K_\ell + 1}}\ \alpha_{\mathrm{NLOS}},
\qquad
\alpha_{\mathrm{NLOS}} \sim \mathcal{CN}(0, 1)
}
$$

其中：

| 符号 | 含义 |
|------|------|
| $K_\ell$ | 第 $\ell$ 条路径的 Rician K 因子（线性），$K_\ell \ge 0$ |
| $\alpha_{\mathrm{LOS}}$ | 确定性 LOS 分量，通常取 $\lvert\alpha_{\mathrm{LOS}}\rvert = 1$，相位取决于 AoA/AoD 几何；在简化模型中可直接取常数 |
| $\alpha_{\mathrm{NLOS}}$ | 随机散射分量，每根天线对独立采样 |

当 $K_\ell = 0$ 时，Rician 退化为 Rayleigh——此时散射分量方差为 1。因此 **Rayleigh 可视为 Rician 在 $K=0$ 时的特例**，实现上可统一处理。

K 因子的典型取值：

| 场景 | $K_\ell$ (dB) | 含义 |
|------|---------------|------|
| 城市宏小区 LOS | 4~10 dB | 弱 LOS |
| 农村开阔地 LOS | 10~20 dB | 强 LOS |
| 室内 LOS | 5~10 dB | 中等 LOS |
| 任何 NLOS 径 | $-\infty$ dB (即 0) | Rayleigh |

### 3.1.4 MIMO CFR 的组装

将各路径贡献在频域叠加，得到第 $m$ 个 OFDM 符号、第 $k$ 个子载波上的 MIMO CFR。对固定收发天线对 $(r,t)$：

$$
\boxed{
H_{r,t}[m, k]
=
\sum_{\ell=0}^{L-1}
\sqrt{P_\ell}\ [\mathbf G_\ell[m]]_{r,t}\ e^{-j 2\pi k \Delta f \tau_\ell}
}
$$

所有 $H_{r,t}[m,k]$ 组成矩阵 $\mathbf H[m,k]$。也可按离散 tap 形式逐元素写出：

$$
[\mathbf H[m, k]]_{r,t}
=
\sum_{n=0}^{N-1}
[\mathbf H[m, n]]_{r,t}\ e^{-j\frac{2\pi}{N}kn}.
$$

### 3.1.5 完整生成流程

对每个 OFDM 符号 $m$：

1. **设定 PDP**：确定路径数 $L$、时延 $\tau_\ell$、功率 $P_\ell$。
2. **决定 LOS/NLOS**：按 $p_{\mathrm{LOS}}$ 抽样，决定首径（或前几条径）是否为 LOS。
3. **逐径生成空间矩阵 $\mathbf G_\ell[m]$**：
   - LOS 径：按 Rician 公式生成，K 因子按场景取定值；
   - NLOS 径：所有元素独立 $\sim \mathcal{CN}(0, 1)$（Rayleigh）。
4. **频域叠加**：按 CFR 公式对各子载波求和。
5. *可选*：若需模拟时间选择性，可通过 Jakes 或 AR(1) 模型在 OFDM 符号间驱动 $\alpha_{\mathrm{NLOS}}$ 和 $\alpha_{\mathrm{LOS}}$ 的相位旋转。

### 3.1.6 模型优缺点

**优点**：

- 实现简单，参数少（仅 PDP + K 因子 + $p_{\mathrm{LOS}}$），适合快速算法验证；
- 各 tap 独立，可直接沿用 SISO 多径的理论框架；
- 逐天线对独立使其成为容量和 BER 理论上界的基准。

**缺点**：

- 完全忽略阵列几何和空间相关性，当天线数较多或间距较小时偏乐观；
- 不同发射-接收天线对之间独立，不能描述 AoD-AoA 耦合、波束效应和预编码增益的空域结构；
- 当 $N_t, N_r$ 较大时，独立假设可导致信道秩显著高估，空间复用层数预测失真。

## 3.2 空间相关 Kronecker 模型

Kronecker 模型是 MIMO 信道建模中最经典的解析模型之一，由 Shiu、Foschini、Gans 和 Kahn 等人在 2000 年前后提出并系统分析。其核心思想是用**可分离的空间相关结构**来刻画多天线之间的统计依赖关系，在理论推导的简洁性和物理准确性之间取得了良好平衡。

### 3.2.1 Kronecker 模型的基本假设

Kronecker 模型的核心假设可以概括为三条：

**假设 1：发射端与接收端空间相关可分离**

这是 Kronecker 模型最本质的假设。它要求任意两根发射天线 $(t_1, t_2)$ 和任意两根接收天线 $(r_1, r_2)$ 之间的互相关可以分解为发射端相关和接收端相关的乘积：

$$
\mathbb E\left[ H_{r_1,t_1} H_{r_2,t_2}^* \right]
=
[\mathbf R_t]_{t_1,t_2} \cdot [\mathbf R_r]_{r_1,r_2}.
$$

换言之，一对天线对 $(r_1, t_1)$ 和 $(r_2, t_2)$ 的信道相关性**完全由发射端两天线之间的空间相关性**和**接收端两天线之间的空间相关性**的乘积决定，两者互不干扰。

该假设等价于：信道矩阵向量化 $\mathbf h = \mathrm{vec}(\mathbf H)$ 的协方差矩阵具有 Kronecker 积结构：

$$
\mathbf R_h
=
\mathbb E\left[ \mathbf h \mathbf h^H \right]
=
\mathbf R_t^T \otimes \mathbf R_r.
$$

**假设 2：发射端散射环境与接收端散射环境统计独立**

物理上，这意味着发射端周围的散射体分布和接收端周围的散射体分布是彼此独立的。数学上，这保证了 $\mathbf R_t$ 和 $\mathbf R_r$ 可以分别独立定义。如果发射端和接收端之间存在公共散射体（例如街道峡谷中两侧建筑同时影响 AoD 和 AoA），则该假设会被破坏。

**假设 3：每个 tap 独立满足上述可分离结构**

对于宽带多径信道，Kronecker 模型通常逐 tap 应用。即对每个时延 tap $n$，可以有不同的相关矩阵 $\mathbf R_{t,n}$ 和 $\mathbf R_{r,n}$。这允许不同时延簇具有不同的角度扩展特性（例如首达簇角度扩展较小，后续簇角度扩展较大）。但如果不同 tap 之间也存在相关性，则 Kronecker 模型无法直接刻画。

### 3.2.2 Kronecker 模型下的离散域 MIMO CIR

在第 2.3 节中，离散 MIMO CIR 定义为 $\mathbf H[m, n] \in \mathbb C^{N_r \times N_t}$，表示第 $m$ 个 OFDM 符号在第 $n$ 个采样 tap 上的信道矩阵。在 Kronecker 模型下，该矩阵按下式生成：

$$
\boxed{
\mathbf H[m,n]
=
\sqrt{P_n}\ \mathbf R_{r,n}^{1/2}\ \mathbf W_n[m]\ \mathbf R_{t,n}^{1/2}
}
$$

其中各符号含义为：

| 符号 | 维度 | 含义 |
|------|------|------|
| $P_n$ | 标量 | 第 $n$ 个 tap 的平均功率，由 PDP 给出 |
| $\mathbf W_n[m]$ | $\mathbb C^{N_r \times N_t}$ | i.i.d. 零均值单位方差复高斯矩阵，$\mathrm{vec}(\mathbf W_n[m]) \sim \mathcal{CN}(\mathbf 0, \mathbf I_{N_r N_t})$ |
| $\mathbf R_{t,n}$ | $\mathbb C^{N_t \times N_t}$ | 发射端空间相关矩阵，描述 $N_t$ 根发射天线之间的统计相关性 |
| $\mathbf R_{r,n}$ | $\mathbb C^{N_r \times N_r}$ | 接收端空间相关矩阵，描述 $N_r$ 根接收天线之间的统计相关性 |

这里 $\mathbf R^{1/2}$ 表示矩阵的 Hermitian 平方根，满足 $\mathbf R^{1/2} (\mathbf R^{1/2})^H = \mathbf R$。

将 $\mathbf H[m,n]$ 逐列向量化 $\mathbf h_n[m] = \mathrm{vec}\{\mathbf H[m,n]\}$，利用恒等式 $\mathrm{vec}(\mathbf A \mathbf B \mathbf C) = (\mathbf C^T \otimes \mathbf A)\,\mathrm{vec}(\mathbf B)$，可得：

$$
\begin{aligned}
\mathbf h_n[m]
&=
\sqrt{P_n}\ \mathrm{vec}\left( \mathbf R_{r,n}^{1/2}\ \mathbf W_n[m]\ \mathbf R_{t,n}^{1/2} \right) \\[4pt]
&=
\sqrt{P_n}\ \left( (\mathbf R_{t,n}^{1/2})^T \otimes \mathbf R_{r,n}^{1/2} \right)\ \mathrm{vec}(\mathbf W_n[m]).
\end{aligned}
$$

令 $\mathbf w_n[m] = \mathrm{vec}(\mathbf W_n[m]) \sim \mathcal{CN}(\mathbf 0, \mathbf I)$，则 $\mathbf h_n[m]$ 是 $\mathbf w_n[m]$ 的线性变换。由复高斯向量的线性变换性质，$\mathbf h_n[m]$ 仍服从零均值复高斯分布，其协方差矩阵为：

$$
\begin{aligned}
\mathbb E\left[ \mathbf h_n[m] \mathbf h_n^H[m] \right]
&=
P_n\ \left( (\mathbf R_{t,n}^{1/2})^T \otimes \mathbf R_{r,n}^{1/2} \right)
\ \mathbb E\left[ \mathbf w_n[m] \mathbf w_n^H[m] \right]
\ \left( (\mathbf R_{t,n}^{1/2})^T \otimes \mathbf R_{r,n}^{1/2} \right)^H \\[4pt]
&=
P_n\ \left( (\mathbf R_{t,n}^{1/2})^T \otimes \mathbf R_{r,n}^{1/2} \right)
\left( (\mathbf R_{t,n}^{1/2})^* \otimes (\mathbf R_{r,n}^{1/2})^H \right) \\[4pt]
&=
P_n\ \left( (\mathbf R_{t,n}^{1/2})^T (\mathbf R_{t,n}^{1/2})^* \right) \otimes
\left( \mathbf R_{r,n}^{1/2} (\mathbf R_{r,n}^{1/2})^H \right).
\end{aligned}
$$

利用 $(\mathbf R_{t,n}^{1/2})^T (\mathbf R_{t,n}^{1/2})^* = ( \mathbf R_{t,n}^{1/2} (\mathbf R_{t,n}^{1/2})^H )^T = \mathbf R_{t,n}^T$，最终得到向量化信道的高斯分布形式：

$$
\boxed{
\mathbf h_n[m]
=
\mathrm{vec}\{\mathbf H[m,n]\}
\sim
\mathcal{CN}\left(
\mathbf 0,\;
P_n\left(\mathbf R_{t,n}^{T}\otimes \mathbf R_{r,n}\right)
\right)
}
$$

该式是 Kronecker 模型的向量化表示，也是后续 LMMSE 信道估计中构造 $\mathbf R_h$ 的基础。

**完整 MIMO CIR 的生成流程**（每个 OFDM 符号 $m$、每个 tap $n$）：

1. 生成 i.i.d. 矩阵 $\mathbf W_n[m] \in \mathbb C^{N_r \times N_t}$，每个元素独立 $\sim \mathcal{CN}(0, 1)$；
2. 左乘 $\mathbf R_{r,n}^{1/2}$，施加接收端空间相关；
3. 右乘 $\mathbf R_{t,n}^{1/2}$，施加发射端空间相关；
4. 乘以 $\sqrt{P_n}$，赋予该 tap 的功率。

### 3.2.3 Kronecker 模型下的 MIMO CFR

由第 2.3 节，第 $m$ 个 OFDM 符号、第 $k$ 个子载波上的 MIMO CFR 为 CIR 的 DFT：

$$
\mathbf H[m, k] = \sum_{n=0}^{N-1} \mathbf H[m, n]\ e^{-j\frac{2\pi}{N}kn}.
$$

将 Kronecker 模型生成的 CIR 代入，并假设各 tap 之间独立（不同时延的散射体统计独立），则：

$$
\boxed{
\mathbf H[m, k]
=
\sum_{n=0}^{N-1}
\sqrt{P_n}\ \mathbf R_{r,n}^{1/2}\ \mathbf W_n[m]\ \mathbf R_{t,n}^{1/2}
\ e^{-j\frac{2\pi}{N}kn}
}
$$

若各 tap 的空间相关特性相同（即所有 tap 共享同一组相关矩阵 $\mathbf R_t, \mathbf R_r$），则可进一步简化为：

$$
\mathbf H[m, k]
=
\mathbf R_r^{1/2}
\left(
\sum_{n=0}^{N-1} \sqrt{P_n}\ \mathbf W_n[m]\ e^{-j\frac{2\pi}{N}kn}
\right)
\mathbf R_t^{1/2}.
$$

该形式清晰地表明：在 Kronecker 假设下，空间相关结构与频率选择性是**解耦**的——相关矩阵不随子载波变化，仅内部的 i.i.d. 衰落部分受 DFT 相位旋转的影响。

对 CFR 进行向量化，同样可以得到其统计分布。由于线性变换保持高斯性，且各 tap 独立，向量化 CFR 的协方差为各 tap 协方差的加权和：

$$
\mathrm{vec}\{\mathbf H[m, k]\}
\sim
\mathcal{CN}\left(
\mathbf 0,\;
\sum_{n=0}^{N-1} P_n \left( \mathbf R_{t,n}^T \otimes \mathbf R_{r,n} \right)
\right).
$$

若所有 tap 共享相同相关矩阵，则协方差简化为 $(\mathbf R_t^T \otimes \mathbf R_r)$（因 $\sum_n P_n = 1$），此时 CFR 的统计特性与 CIR 的单个 tap 一致——这也是 Kronecker 模型的一个便利性质。

### 3.2.4 空间相关矩阵的构造与推导

Kronecker 模型的实用性取决于能否合理地构造出 $\mathbf R_t$ 和 $\mathbf R_r$。下面给出从物理参数到相关矩阵的完整推导链路。

#### 3.2.4.1 从角度功率谱推导相关矩阵

考虑一维均匀线阵（ULA），天线间距为 $d$，波长为 $\lambda$。定义归一化方向余弦（或角频率）：

$$
\Omega = \frac{d}{\lambda} \sin\theta,
$$

其中 $\theta \in [-\pi/2, \pi/2]$ 为相对于阵列法线的物理到达角（或出发角）。对于间距为半波长 $d = \lambda/2$ 的阵列，$\Omega \in [-0.5, 0.5]$。

在该阵列上，角度为 $\theta$ 的窄带平面波的阵列导向矢量为：

$$
\mathbf a(\Omega)
=
\frac{1}{\sqrt{N}}
\begin{bmatrix}
1 \\ e^{-j2\pi\Omega} \\ e^{-j2\pi\cdot 2\Omega} \\ \vdots \\ e^{-j2\pi(N-1)\Omega}
\end{bmatrix}.
$$

物理上，第 $\ell$ 条路径（或散射簇中的一条射线）以角度 $\theta_\ell$ 到达阵列，对第 $p$ 根天线的贡献携带相位 $e^{-j2\pi p\,\Omega_\ell}$。接收信号的协方差矩阵是所有路径贡献的统计平均。

假设散射体在角度域上连续分布，其功率密度由**角度功率谱（Angular Power Spectrum, APS）** $p(\Omega)$ 描述，满足 $\int p(\Omega)\,d\Omega = 1$。则接收空间相关矩阵的第 $(p, q)$ 个元素（第 $p$ 根天线与第 $q$ 根天线之间的复相关）为：

$$
\boxed{
[\mathbf R]_{p,q}
=
\mathbb E\left[ h_p h_q^* \right]
=
\int e^{-j2\pi(p-q)\Omega}\ p(\Omega)\ d\Omega
}
$$

**推导过程**：

考虑第 $p$ 根天线的信道系数为所有角度路径贡献的叠加：

$$
h_p = \int \alpha(\Omega)\ e^{-j2\pi p\,\Omega}\ d\Omega,
$$

其中 $\alpha(\Omega)$ 是复增益密度函数，满足 $\mathbb E[\alpha(\Omega)\alpha^*(\Omega')] = p(\Omega)\,\delta(\Omega - \Omega')$（不同角度不相关）。则：

$$
\begin{aligned}
[\mathbf R]_{p,q}
&= \mathbb E[h_p h_q^*] \\[4pt]
&= \mathbb E\left[
\int \alpha(\Omega) e^{-j2\pi p\Omega} d\Omega
\int \alpha^*(\Omega') e^{j2\pi q\Omega'} d\Omega'
\right] \\[4pt]
&= \iint
\mathbb E[\alpha(\Omega)\alpha^*(\Omega')]\
e^{-j2\pi(p\Omega - q\Omega')}\ d\Omega\,d\Omega' \\[4pt]
&= \iint p(\Omega)\,\delta(\Omega - \Omega')\
e^{-j2\pi(p\Omega - q\Omega')}\ d\Omega\,d\Omega' \\[4pt]
&= \int p(\Omega)\ e^{-j2\pi(p-q)\Omega}\ d\Omega.
\end{aligned}
$$

**关键观察**：

- $[\mathbf R]_{p,q}$ 仅依赖于天线索引差 $\Delta = p - q$，这是 ULA 的**Toeplitz** 结构，也是**广义平稳（WSS）**假设在空间域的体现。
- 相关矩阵的元素本质上是 $p(\Omega)$ 的逆 Fourier 变换（在变量 $(p-q)$ 处取值），因此 $p(\Omega)$ 完全决定了空间相关结构。
- 该公式同时适用于发射端和接收端：对发射端，$p(\Omega)$ 为出发角功率谱（APS at transmitter, APS$_{\text{Tx}}$）；对接收端，$p(\Omega)$ 为到达角功率谱（APS at receiver, APS$_{\text{Rx}}$）。

**两种常见的角度功率谱及其对应的相关矩阵**：

| 角度分布 | $p(\Omega)$ | $[\mathbf R]_{p,q}$ |
|----------|-------------|---------------------|
| 均匀分布（全向散射） | $p(\Omega) = 1$（归一化区间内） | $\mathrm{sinc}(2(p-q)\Omega_{\max})$ |
| 截断高斯分布 | $p(\Omega) \propto \exp(-\frac{(\Omega-\bar\Omega)^2}{2\sigma_\Omega^2})$ | 无闭式，需数值积分 |
| 拉普拉斯分布 | $p(\Omega) \propto \exp\!\left(-\frac{\sqrt{2}\,\lvert\Omega-\bar{\Omega}\rvert}{\sigma_\Omega}\right)$ | 无闭式，需数值积分 |
| 单点源（纯 LOS） | $p(\Omega) = \delta(\Omega - \Omega_0)$ | $e^{-j2\pi(p-q)\Omega_0}$（秩 1） |

拉普拉斯分布在宏小区场景中与实测数据吻合较好；截断高斯分布在微小区中更为常见。3GPP TR 38.901 中，簇内射线的角度偏移通常建模为拉普拉斯分布或包裹高斯分布。

#### 3.2.4.1b 角度扩展（Angular Spread）的生成

角度功率谱 $p(\Omega)$ 中的参数 $\sigma_\Omega$ 即为**角度扩展**——它量化了散射能量在角度域上的分散程度，是 MIMO 信道建模中最核心的大尺度参数之一。

**1) 角度扩展的严格定义**

角度扩展（Angular Spread, AS）定义为角度功率谱的二阶中心矩的平方根（RMS angular spread）：

$
\boxed{
\sigma_{\mathrm{AS}}
\triangleq
\sqrt{
\frac{\int (\Omega - \bar\Omega)^2\, p(\Omega)\, d\Omega}
{\int p(\Omega)\, d\Omega}
},
\qquad
\bar\Omega \triangleq
\frac{\int \Omega\, p(\Omega)\, d\Omega}
{\int p(\Omega)\, d\Omega}
}
$

其中 $\bar\Omega$ 是**平均角度**（mean angle），即角度域功率分布的一阶矩。以度为单位时，定义相同，只需将 $\Omega$ 换算为物理角度 $\theta$。

**物理直觉**：

- $\sigma_{\mathrm{AS}}$ 小 → 散射体集中在窄角度范围内 → 多根天线接收/发射的信号高度相关 → 信道矩阵趋于低秩 → 空间复用能力弱，但波束赋形增益高；
- $\sigma_{\mathrm{AS}}$ 大 → 散射体分布广泛（丰富散射环境）→ 不同天线的信道趋于独立 → 信道矩阵趋于满秩 → 空间复用层数多。

在常见角度分布下，$\sigma_\Omega$ 与 AS 的关系如下：

| 角度分布 | $\sigma_{\mathrm{AS}}$ 与参数 $\sigma_\Omega$ 的关系 |
|----------|-----------------------------------------------------|
| 截断高斯 | AS ≈ $\sigma_\Omega$（当截断区间远大于 $\sigma_\Omega$ 时精确成立） |
| 拉普拉斯 | AS = $\sigma_\Omega$（严格相等，$\sigma_\Omega$ 即为 RMS 角度扩展） |
| 均匀分布（半宽 $\Omega_{\max}$） | AS = $\Omega_{\max}/\sqrt{3}$ |

因此在实际仿真中，指定 AS 等价于指定 $p(\Omega)$ 的宽度参数。

**2) 如何用平均角度和角度扩展生成单条路径的角度**

在仿真中，对每条路径（或簇、射线）的角度取值可按如下流程生成：

**Step 1：确定大尺度参数**
- 选择场景（UMa/UMi/RMa/InH 等）和传播条件（LOS/NLOS）；
- 查表或根据经验公式获得该场景下的平均角度 $\bar\Omega$（或 $\bar\theta$）和角度扩展 $\sigma_{\mathrm{AS}}$。

**Step 2：选择角度分布类型**
- 最常见选择：拉普拉斯分布（宏小区）或包裹高斯分布（微小区）；
- 纯 LOS 场景可对 LOS 径取 $\sigma_{\mathrm{AS}} \to 0$，此时角度退化为确定性值 $\bar\Omega$。

**Step 3：为每条路径/射线生成具体角度**

假设选定拉普拉斯分布，第 $\ell$ 条路径（或簇中第 $q$ 条射线）的角度为：

$
\Omega_\ell = \bar\Omega + \Delta\Omega_\ell,
$

其中 $\Delta\Omega_\ell$ 是角度偏移量，按拉普拉斯分布生成：

$
p(\Delta\Omega) = \frac{1}{\sqrt{2}\,\sigma_{\mathrm{AS}}} \exp\!\left( -\frac{\sqrt{2}\,|\Delta\Omega|}{\sigma_{\mathrm{AS}}} \right).
$

若使用高斯分布：

$
p(\Delta\Omega) = \frac{1}{\sqrt{2\pi}\,\sigma_{\mathrm{AS}}} \exp\!\left( -\frac{(\Delta\Omega)^2}{2\sigma_{\mathrm{AS}}^2} \right),
\qquad \text{需截断至合理角度范围}.
$

**Step 4：将归一化角度转换为物理角度**

$\Omega = \frac{d}{\lambda}\sin\theta$ → $\theta = \arcsin\!\left(\frac{\lambda}{d}\,\Omega\right)$。对于 $d = \lambda/2$，$\theta = \arcsin(2\Omega)$。

**完整生成示例**（发射端 Azimuth AoD，拉普拉斯分布）：

$
\begin{aligned}
&\text{给定: } \bar\theta_{\mathrm{AoD}} = 30^\circ,\quad \sigma_{\mathrm{ASD}} = 5^\circ,\quad L = 3\ \text{条路径}. \\[4pt]
&\text{step 1: } \bar\Omega_t = \frac{d}{\lambda}\sin\bar\theta_{\mathrm{AoD}} = 0.5 \times \sin 30^\circ = 0.25. \\[4pt]
&\text{step 2: } \sigma_\Omega = \frac{d}{\lambda}\sin(\sigma_{\mathrm{ASD}}) \approx \frac{d}{\lambda}\,\sigma_{\mathrm{ASD}}\ \text{(小角度近似)}. \\[4pt]
&\text{step 3: 对每条路径 }\ell,\ \Omega_{t,\ell} = \bar\Omega_t + \Delta\Omega_\ell,\quad \Delta\Omega_\ell \sim \mathrm{Laplace}(0, \sigma_\Omega/\sqrt{2}). \\[4pt]
&\text{step 4: } \theta_{t,\ell} = \arcsin\!\left(\frac{\lambda}{d}\,\Omega_{t,\ell}\right).
\end{aligned}
$

**3) 多径/多簇场景下 AS 的含义**

在宽带 MIMO 信道中，不同时延簇通常具有不同的 AS。TR 38.901 定义了两级角度扩展：

| 级别 | 参数 | 含义 |
|------|------|------|
| **簇级角度扩展** | $c_{\mathrm{ASD}}$, $c_{\mathrm{ASA}}$, $c_{\mathrm{ZSD}}$, $c_{\mathrm{ZSA}}$ | 簇内 20 条射线的角度分布宽度 |
| **全局角度扩展** | ASD, ASA, ZSD, ZSA | 所有簇在全空间上的总体角度扩展（大尺度参数） |

生成流程为：
1. 先按场景获得**全局 AS**（如 UMa NLOS 的 ASD ≈ 22°）；
2. 用全局 AS 缩放 CDL 表中的归一化簇级角度偏移 $c_{\mathrm{ASD}}$ 等；
3. 每个簇内射线以簇中心角为均值、$c_{\mathrm{ASD}} \cdot \mathrm{ASD}_{\text{desired}}$ 为扩展量生成具体角度。

**4) 典型场景的角度扩展参考值**

以下数值来自 3GPP TR 38.901（UMa 场景，方位角 AS，单位：度）：

| 场景 | LOS/NLOS | ASD | ASA | ZSD | ZSA |
|------|----------|-----|-----|-----|-----|
| UMa | LOS | $\lg \mathrm{ASD} \sim \mathcal{N}(1.06, 0.28)$ | $\lg \mathrm{ASA} \sim \mathcal{N}(1.60, 0.18)$ | $\lg \mathrm{ZSD} \sim \mathcal{N}(0.64, 0.32)$ | $\lg \mathrm{ZSA} \sim \mathcal{N}(0.90, 0.23)$ |
| UMa | NLOS | $\lg \mathrm{ASD} \sim \mathcal{N}(1.46, 0.28)$ | $\lg \mathrm{ASA} \sim \mathcal{N}(1.79, 0.20)$ | $\lg \mathrm{ZSD} \sim \mathcal{N}(0.81, 0.23)$ | $\lg \mathrm{ZSA} \sim \mathcal{N}(0.96, 0.17)$ |

注：此处 AS 值服从对数正态分布，仿真时先抽样 $\lg\mathrm{AS}$ 再取 $10^{\lg\mathrm{AS}}$。表中给出的是 $\lg\mathrm{AS}$ 的均值 $\mu$ 和标准差 $\sigma$。例如 UMa NLOS 的 ASD 对数均值为 1.46，即几何平均 ASD ≈ $10^{1.46} \approx 28.8°$（该值因具体参数集而有波动，与 CDL 表中的角度扩展区间一致）。

**5) 小结**

| 问题 | 答案 |
|------|------|
| AS 的物理含义 | 角度域散射能量的 RMS 宽度 |
| 如何参数化 $p(\Omega)$ | 指定分布类型 + 平均角度 $\bar\Omega$ + 角度扩展 $\sigma_{\mathrm{AS}}$ |
| 如何为单条路径赋角度 | $\Omega_\ell = \bar\Omega + \Delta\Omega$，$\Delta\Omega$ 按指定分布随机生成 |
| AS 与空间相关的联系 | $[\mathbf R]_{p,q} = \int e^{-j2\pi(p-q)\Omega}p(\Omega)d\Omega$，AS 决定 $p(\Omega)$ 宽度 → 决定相关矩阵衰减速度 |
| AS 与信道秩的关系 | AS 越大 → 相关越弱 → 信道秩越高 |

#### 3.2.4.2 指数相关模型

当不需要精确匹配特定角度分布、或仅需调节相关强度时，可使用**指数相关模型**（Exponential Correlation Model）简化：

$$
\boxed{
[\mathbf R_t]_{p,q} = \rho_t^{|p-q|},
\qquad
[\mathbf R_r]_{p,q} = \rho_r^{|p-q|}
}
$$

其中 $0 \le \rho_t, \rho_r < 1$ 为相关系数。

**模型性质**：

- $\rho_t = 0$：发射天线间完全独立（i.i.d.）；
- $\rho_t \to 1$：发射天线间完全相关（信道矩阵秩退化为 1）；
- 相关随天线索引差 $|p-q|$ 呈指数衰减——物理直觉是天间距越远，信道越独立；
- 相关矩阵为 Toeplitz、Hermitian、正定矩阵。

该模型的合理性在于：对于一个角度扩展有限、中心角为零的均匀散射环境，$p(\Omega)$ 近似为拉普拉斯分布时，其 Fourier 变换在整数采样点上的衰减近似为指数形式。但需注意，指数相关模型是一种**唯象近似**，并不严格对应某个具体的物理角度分布。

**$\rho$ 与有效秩的关系**：

定义相关矩阵的**有效秩**（Effective Rank）：

$$
r_{\mathrm{eff}}(\mathbf R) = \exp\left( -\sum_i \frac{\lambda_i}{\sum_j \lambda_j} \ln \frac{\lambda_i}{\sum_j \lambda_j} \right),
$$

其中 $\lambda_i$ 是 $\mathbf R$ 的特征值。当 $\rho = 0$ 时，所有特征值相等，$r_{\mathrm{eff}} = N$；当 $\rho \to 1$ 时，最大特征值主导，$r_{\mathrm{eff}} \to 1$。信道矩阵 $\mathbf H = \mathbf R_r^{1/2} \mathbf W \mathbf R_t^{1/2}$ 的秩受 $\min(r_{\mathrm{eff}}(\mathbf R_t), r_{\mathrm{eff}}(\mathbf R_r))$ 约束，因此 $\rho$ 越大，可支持的空间复用层数越少。

#### 3.2.4.3 利用阵列导向矢量构造相关矩阵

除了从角度功率谱积分，也可以从离散路径参数直接构造相关矩阵。这在几何模型中更为自然。这里的 $p(\Omega)$ 不是 CIR；它不描述路径到达的**时延位置**，也不包含每条路径的瞬时复增益相位，而是描述信道功率在**角度方向**上的统计分布。可以把它理解为角度域的功率谱（APS），类似 PDP 描述功率随时延的分布，但变量从 $\tau$ 换成了归一化角度 $\Omega$。

若角度域只有 $L$ 个离散方向，则可将 $p(\Omega)$ 写成离散角度分布：

$$
p(\Omega) = \sum_{\ell=1}^{L} \gamma_\ell\ \delta(\Omega - \Omega_\ell),
\qquad \sum_\ell \gamma_\ell = 1,
$$

其中 $\gamma_\ell$ 为第 $\ell$ 条路径的归一化功率，$\Omega_\ell$ 为其归一化方向余弦。代入积分公式：

$$
[\mathbf R]_{p,q}
=
\sum_{\ell=1}^{L} \gamma_\ell\ e^{-j2\pi(p-q)\Omega_\ell}
=
\sum_{\ell=1}^{L} \gamma_\ell\
[\mathbf a(\Omega_\ell)]_p\ [\mathbf a(\Omega_\ell)]_q^*.
$$

因此相关矩阵也可写为各路径导向矢量的加权外积：

$$
\boxed{
\mathbf R = \sum_{\ell=1}^{L} \gamma_\ell\ \mathbf a(\Omega_\ell)\ \mathbf a^H(\Omega_\ell)
}
$$

这是相关矩阵的**谱表示**，清晰地表明 $\mathbf R$ 的秩不超过路径数 $L$。该形式便于在仿真中直接通过给定的 AoA/AoD 集合生成精确的相关矩阵，而不需要数值积分。

##### 3.2.4.3.1 实际中如何获得角度功率谱

实际系统中，$p(\Omega)$ 通常不是直接已知的解析函数，而是通过测量、标准化信道模型或经验假设得到。

**1) 由阵列测量估计**

若接收端有天线阵列，可在多个时间快照或多个导频资源上测量信道向量：

$$
\mathbf h(t)
=
\sum_{\ell=1}^{L}
\alpha_\ell(t)\ \mathbf a(\Omega_\ell).
$$

先估计空间协方差：

$$
\hat{\mathbf R}
=
\frac{1}{T}\sum_{t=1}^{T}\mathbf h(t)\mathbf h^H(t).
$$

再由 $\hat{\mathbf R}$ 估计角度功率谱。对 ULA，相关矩阵元素与角度功率谱近似满足 Fourier 对偶关系：

$$
[\mathbf R]_{p,q}
=
\int p(\Omega)e^{-j2\pi(p-q)\Omega}d\Omega.
$$

因此可用波束扫描、Bartlett、Capon/MVDR、MUSIC/ESPRIT、SAGE 或稀疏恢复方法估计角度谱。

**2) 由标准化信道模型生成**

在链路级仿真中，更常见的做法是从标准模型生成角度参数。例如 3GPP TR 38.901 根据场景、LOS/NLOS 状态、延迟扩展和角度扩展生成簇/射线的 AoA、AoD、功率等参数。得到路径集合

$$
\{(\Omega_\ell,\gamma_\ell)\}_{\ell=1}^{L}
$$

后，可直接写成：

$$
p(\Omega)=\sum_{\ell=1}^{L}\gamma_\ell\delta(\Omega-\Omega_\ell),
$$

或构造相关矩阵：

$$
\mathbf R
=
\sum_{\ell=1}^{L}
\gamma_\ell\mathbf a(\Omega_\ell)\mathbf a^H(\Omega_\ell).
$$

**3) 由经验分布假设**

如果只是算法验证，也可直接假设角度分布。例如：

- 均匀分布：表示散射较丰富、角度覆盖较宽；
- 高斯分布：表示围绕某个主方向有有限角度扩展；
- 拉普拉斯分布：宏小区建模中常见；
- 单点 delta 分布：表示纯 LOS 或强主径。

典型拉普拉斯角度谱可写为：

$$
p(\Omega)
\propto
\exp\left(
-\frac{\sqrt{2}\lvert\Omega-\bar{\Omega}\rvert}{\sigma_\Omega}
\right),
$$

其中 $\bar{\Omega}$ 是平均方向，$\sigma_\Omega$ 描述角度扩展。

因此实际建模中常见链路是：

$$
\text{测量信道}
\rightarrow
\hat{\mathbf R}
\rightarrow
\text{角度谱估计}
\rightarrow
p(\Omega)
$$

或：

$$
\text{选择场景}
\rightarrow
\text{生成 AoA/AoD 与角度扩展}
\rightarrow
p(\Omega)\ \text{或}\ \mathbf R.
$$

### 3.2.5 模型优缺点

**优点**：

- **理论简洁**：向量化协方差为 Kronecker 积形式，便于推导容量、估计误差的解析表达式；
- **参数少**：仅需两个 $N_t \times N_t$ 和 $N_r \times N_r$ 的相关矩阵，而非 $N_t N_r \times N_t N_r$ 的全相关矩阵；
- **可分性**：发射端与接收端相关可以独立调节，适合研究不对称场景；
- **闭合形式**：相关矩阵可通过角度功率谱或指数模型快速生成，无需射线追踪。

**缺点**：

- **不能描述 AoD-AoA 的联合耦合**：物理上，某个特定发射角度往往只对应特定接收角度（如街道峡谷的波导效应），这种耦合在 Kronecker 分离假设下被平均化，窄化了信道的角度-角度双选性结构；
- **对某些传播场景预测偏乐观/偏悲观**：
  - 当实际存在强耦合时，Kronecker 模型会**高估**信道容量（因为分离假设降低了有效条件数）；
  - 在某些散射丰富的场景中，又可能**低估**容量（因为分离假设限制了互信息的上界）；
- **不适于描述"针孔"（keyhole/pinhole）信道**：当发射端和接收端都处于丰富散射环境，但两者之间仅通过一个狭窄开口耦合时，信道矩阵秩为 1 但相关矩阵满秩，Kronecker 模型无法复现这一现象；
- **大规模阵列下偏差增大**：随着天线数增加，近场效应和球面波前不可忽略，平面波假设下的相关矩阵构造精度下降。

针对联合耦合问题，Weichselberger 模型通过特征分解引入了发射-接收联合特征模式，是 Kronecker 模型的直接推广。实际工程仿真中，如果对 AoD-AoA 耦合有明确要求，建议改用 3.2 节之外的几何模型（§3.3.1）或 3GPP CDL 模型（§3.4）。

## 3.3 角度域信道建模与稀疏表示

前两节（§3.1、§3.2）从统计相关矩阵的角度建模 MIMO 信道，本节则从**阵列信号处理**的角度出发：利用天线阵列的相位结构，将连续的空间角度（AoA/AoD）映射为信道矩阵的代数结构。核心逻辑链为：

$
\text{ULA 相位差}\ \longrightarrow\ \text{导向矢量}\ \longrightarrow\ \text{连续角度多径叠加}\ \longrightarrow\ \text{角度离散化}\ \longrightarrow\ \text{DFT 基矩阵}\ \longrightarrow\ \text{虚拟信道稀疏表示}
$

这一框架不仅是信道建模的工具，更是压缩感知信道估计的理论基础。

### 3.3.1 连续角度域——从阵列相位差到导向矢量

**1) 单个平面波的阵列响应**

考虑接收端为 $N$ 元均匀线阵（ULA），阵元间距 $d$。一远场窄带平面波以角度 $\theta$（相对于阵列法线）入射，波长为 $\lambda$。

相邻阵元的**波程差**为 $d\sin\theta$，对应的**相位差**为：

$$
\Delta\phi = 2\pi \cdot \frac{d\sin\theta}{\lambda}
$$

定义**方向余弦**（direction cosine）：

$$
\boxed{
\Omega \triangleq \frac{d}{\lambda}\sin\theta,
\qquad \Omega \in \left[-\frac{d}{\lambda},\ \frac{d}{\lambda}\right]
}
$$

以第 $0$ 号阵元为参考（相位 = 0），第 $p$ 号阵元相对于参考的累积相位差为 $p \cdot \Delta\phi = 2\pi p\Omega$。因此 $N$ 元阵列在该入射方向上的响应——即**阵列导向矢量**（steering vector）为：

$$
\boxed{
\mathbf a(\Omega)
=
\frac{1}{\sqrt{N}}
\begin{bmatrix}
1 \\
e^{-j2\pi\Omega} \\
e^{-j2\pi\cdot 2\Omega} \\
\vdots \\
e^{-j2\pi(N-1)\Omega}
\end{bmatrix}
}
$$

其中 $1/\sqrt{N}$ 为归一化因子，保证 $\|\mathbf a(\Omega)\|^2 = 1$。

**关键观察**：导向矢量的第 $p$ 个元素为 $e^{-j2\pi p\Omega}$——这正是以天线索引 $p$ 为"离散时间"、$\Omega$ 为"数字频率"的**复指数序列**。这一形式是后续 DFT 离散化的根源。

**2) 多径叠加——几何多径信道**

设信道有 $L$ 条可分辨物理路径，第 $\ell$ 条路径具有：
- 发射端离开角方向余弦 $\Omega_{t,\ell}$，接收端到达角方向余弦 $\Omega_{r,\ell}$；
- 复路径增益 $\alpha_\ell[m]$；
- 传播时延 $\tau_\ell$。

则 MIMO 频域信道矩阵为各路径秩-1 分量的叠加：

$$
\boxed{
\mathbf H[m,k]
=
\sum_{\ell=0}^{L-1}
\alpha_\ell[m]\,
e^{-j2\pi k\Delta f\tau_\ell}\,
\mathbf a_r(\Omega_{r,\ell})\,
\mathbf a_t^H(\Omega_{t,\ell})
}
$$

若需考虑多普勒：$\alpha_\ell[m] = \beta_\ell\, e^{j2\pi\nu_\ell mT_{\mathrm{sym}}}$。

**结构与秩**：每一条路径的空间分量 $\mathbf a_r(\Omega_{r,\ell})\mathbf a_t^H(\Omega_{t,\ell})$ 是**秩-1** 矩阵——它在 $N_r \times N_t$ 维空间中仅占据一个方向。整个信道矩阵是 $L$ 个秩-1 矩阵的加权和，因此

$
\operatorname{rank}(\mathbf H) \leq \min(N_t, N_r, L).
$

路径角度越分散 → 各秩-1 分量的方向差异越大 → 信道矩阵秩越高，空间复用潜力越大。反之，角度扩展越小（如 LOS 主导场景），信道越接近低秩。

**连续角度**：上述模型中，$\Omega_{t,\ell}$ 和 $\Omega_{r,\ell}$ 可取 $[-d/\lambda,\ d/\lambda]$ 内的**任意连续实数值**，这是"连续角度域"的含义。

**3) 波束赋形与空间标签的内积**

导向矢量的另一个名字是**空间标签**（Spatial Signature）——它是对特定入射方向的"数字指纹"。理解导向矢量之间内积的含义，是连接连续角度与离散 DFT 基的关键桥梁，也是从阵列物理结构自然引出角度离散化的动机所在。

**波束赋形问题**：假设我们希望将发射功率集中在方向 $\Omega_0$ 上。发射信号为标量 $s$，各天线上的复加权系数组成向量 $\mathbf w \in \mathbb{C}^N$，满足功率约束 $\|\mathbf w\|^2 = 1$。发送出去的信号向量为 $\mathbf w s$，在远场方向 $\Omega$ 处的接收信号（不计路径损耗）为

$
y(\Omega) = \mathbf a^H(\Omega)\,\mathbf w\,s.
$

若目标方向为 $\Omega_0$，则接收功率正比于 $\bigl|\mathbf a^H(\Omega_0)\mathbf w\bigr|^2$。由 **Cauchy-Schwarz 不等式**，在 $\|\mathbf a(\Omega_0)\|^2 = 1$ 且 $\|\mathbf w\|^2 = 1$ 的条件下

$
\bigl|\mathbf a^H(\Omega_0)\mathbf w\bigr|^2 \le \|\mathbf a(\Omega_0)\|^2\ \|\mathbf w\|^2 = 1,
$

等号成立当且仅当 $\mathbf w$ 与 $\mathbf a(\Omega_0)$ 共线。因此**最优波束赋形向量即导向矢量本身**：

$$
\boxed{\mathbf w_{\mathrm{opt}} = \mathbf a(\Omega_0)}.
$$

这称为**匹配滤波波束赋形**（Matched Filter Beamforming）：每根天线以导向矢量的共轭作为复权重，使各天线信号在目标方向上同相叠加。

**波束图样**：若以 $\mathbf a(\Omega_0)$ 为波束赋形向量，在其他方向 $\Omega$ 上的增益为两个空间标签的内积模平方，即**阵列的波束图样**（Beam Pattern）：

$$
\boxed{
G(\Omega_0, \Omega)
\triangleq \bigl|\mathbf a^H(\Omega_0)\,\mathbf a(\Omega)\bigr|^2
= \left|\frac{1}{N}\sum_{p=0}^{N-1} e^{-j2\pi p(\Omega - \Omega_0)}\right|^2
}.
$$

展开有限等比级数求和，得到其闭式表达——**Dirichlet 核的平方**：

$
G(\Omega_0, \Omega) = \left|
\frac{\sin\!\bigl[\pi N(\Omega - \Omega_0)\bigr]}
{N\sin\!\bigl[\pi(\Omega - \Omega_0)\bigr]}
\right|^2.
$

波束图样具有如下关键性质：

- 在 $\Omega = \Omega_0$ 处取最大值 $1$（主瓣峰值）；
- 在 $\Omega = \Omega_0 \pm \frac{k}{N}$（$k = 1, 2, \dots, N-1$）处**严格为零**——天线数 $N$ 越大，零点越密集，波束越尖锐；
- 主瓣零点间距（即主瓣宽度）约为 $2/N$（以 $\Omega$ 计），因此天线越多指向性越强、空间角度分辨能力越高。

**关键洞察——正交性的来源与离散化的动机**：取离散角度网格 $\Omega_i = i/N$（$i = 0, 1, \dots, N-1$），则对任意两个网格点 $\Omega_i \neq \Omega_j$：

$
\mathbf a^H(\Omega_i)\,\mathbf a(\Omega_j)
= \frac{1}{N}\sum_{p=0}^{N-1} e^{-j2\pi p(j-i)/N}
= \begin{cases}
1, & i = j, \\[2pt]
0, & i \neq j,
\end{cases}
$

即**不同离散网格点上的空间标签彼此正交**。全部 $N$ 个离散角度网格点 $\{\Omega_i = i/N\}_{i=0}^{N-1}$ 上的空间标签 $\{\mathbf a(\Omega_i)\}$ 构成 $\mathbb{C}^N$ 的一组**完备标准正交基**。

这一结论揭示了角度离散化的深层道理：天线数 $N$ 有限意味着阵列只能在 $N$ 个正交方向上提供独立的空间自由度——这 $N$ 个正交方向恰好由 DFT 基给出，不是人为选择，而是阵列物理结构的内禀性质。在 $d = \lambda/2$ 的常用配置下，$|\Omega_i| \le 1/2$ 对应的方向是物理可实现的真实空间指向；$|\Omega_i| > 1/2$ 对应"不可见区"，但作为完备基必须保留（类比数字信号处理中超出 Nyquist 频率的 DFT 分量）。

下面将这一正交基形式化地写为 DFT 基矩阵，并完成离散角度域的系统构建。


### 3.3.2 离散角度域——DFT 基向量与虚拟信道表示

前一小节中，导向矢量 $\mathbf a(\Omega)$ 的参数 $\Omega$（方向余弦）在 $[-d/\lambda,\ d/\lambda]$ 内**连续取值**。本节回答一个核心问题：**如何将连续的角度参数离散化为一组固定的正交基，使得任何一个导向矢量都能用这组基的线性组合来表示？**

答案在于导向矢量的代数形式与离散时间傅里叶变换（DTFT）之间的深刻联系。

**1) 导向矢量即空间序列的 DTFT 核**

回顾 3.3.1 中导出的导向矢量：

$$
\mathbf a(\Omega)
= \frac{1}{\sqrt{N}}
\begin{bmatrix} 1 \\ e^{-j2\pi\Omega} \\ e^{-j2\pi\cdot 2\Omega} \\ \vdots \\ e^{-j2\pi(N-1)\Omega} \end{bmatrix}
$$

其第 $p$ 个元素（$p = 0, 1, \dots, N-1$）为：

$$
[\mathbf a(\Omega)]_p = \frac{1}{\sqrt{N}}\ e^{-j2\pi\Omega p}
$$

将天线索引 $p$ 视为离散"时间"变量、$\Omega$ 视为"数字频率"，这正是长度为 $N$ 的矩形窗序列 $\frac{1}{\sqrt{N}}\operatorname{rect}_N[p]$ 的**DTFT 核**在频率点 $\Omega$ 处的取值。换言之，$\mathbf a(\Omega)$ 不是别的，正是**在 $N$ 个等间距采样点 $p=0,1,\dots,N-1$ 上对复指数信号 $e^{-j2\pi\Omega t}$ 的离散观测**。

**2) 均匀频率采样 $\to$ DFT 基**

现在对频率 $\Omega$ 进行**均匀采样**。将区间 $[0, 1)$ 划分为 $N$ 个等间距网格点：

$$
\boxed{
\Omega_i = \frac{i}{N},\qquad i = 0, 1, \dots, N-1
}
$$

取整是因为复指数满足 $\mathbf a(\Omega + 1) = \mathbf a(\Omega)$（周期为 $1$），$[0,1)$ 内的 $N$ 个采样点覆盖了一个完整周期。

**为什么恰好是 $N$ 个网格点？** 能不能用 $M > N$ 个更密的网格？

- **完备正交基的约束**：$\mathbf a(\Omega_i)$ 都是 $\mathbb{C}^N$ 中的向量，而 $\mathbb{C}^N$ 是 $N$ 维空间。正交基最多只有 $N$ 个向量。取 $\Omega_i = i/N$ 时，这 $N$ 个导向矢量恰好两两正交，构成标准正交基。若用 $M > N$ 个网格点，则这些导向矢量必然线性相关——它们构成的是**冗余框架（overcomplete frame）**，不再是一组基。$N$ 是正交基的「临界采样」上限。
- **过采样的实际用途**：实际系统（如 5G NR）常用 $O \times N$ 个波束（过采样因子 $O = 2$ 或 $4$）。这些过采样波束之间依然正交，且能提供更密的角度候选方向，减少 off-grid 偏差。空间自由度仍是 $N$——过采样是让候选方向更密，不是让你同时用更多正交波束。

当 $d = \lambda/2$（常用半波长间距）时，$\Omega = \frac{1}{2}\sin\theta \in [-1/2,\ 1/2]$。此时 $\Omega_i$ 中落在 $[-1/2,\ 1/2]$ 内的部分对应物理可实现的空间角度；落在区间外的部分（$|\Omega_i| > 1/2$）在物理上不对应任何实角度，**但作为 DFT 谱的完整基必须保留**——好比数字信号处理的频域分析中，整个 $[0, 2\pi)$ 区间必须全部考虑。

将 $\Omega_i = i/N$ 逐一代入 $\mathbf a(\cdot)$：

$$
\begin{aligned}
\mathbf a(\Omega_0) &= \frac{1}{\sqrt{N}}
\begin{bmatrix} 1 & 1 & 1 & \cdots & 1 \end{bmatrix}^T \\[4pt]
\mathbf a(\Omega_1) &= \frac{1}{\sqrt{N}}
\begin{bmatrix} 1 & e^{-j2\pi/N} & e^{-j2\pi\cdot 2/N} & \cdots & e^{-j2\pi(N-1)/N} \end{bmatrix}^T \\[4pt]
\mathbf a(\Omega_i) &= \frac{1}{\sqrt{N}}
\begin{bmatrix} 1 & e^{-j2\pi i/N} & e^{-j2\pi\cdot 2i/N} & \cdots & e^{-j2\pi(N-1)i/N} \end{bmatrix}^T
\end{aligned}
$$

每一个 $\mathbf a(\Omega_i)$ 的第 $p$ 个元素为 $\frac{1}{\sqrt{N}} e^{-j2\pi p i/N}$——**这正是 $N$ 点 DFT 矩阵的共轭转置的第 $i$ 列**。

**3) 构造酉 DFT 基矩阵**

将所有 $N$ 个采样的导向矢量按列排列，得到**酉 DFT 基矩阵**：

$$
\boxed{
\mathbf U
\triangleq
\begin{bmatrix}
\mid & \mid & & \mid \\
\mathbf a(\Omega_0) & \mathbf a(\Omega_1) & \cdots & \mathbf a(\Omega_{N-1}) \\
\mid & \mid & & \mid
\end{bmatrix}
=
\frac{1}{\sqrt{N}}
\begin{bmatrix}
1 & 1 & 1 & \cdots & 1 \\
1 & \omega^{-1} & \omega^{-2} & \cdots & \omega^{-(N-1)} \\
1 & \omega^{-2} & \omega^{-4} & \cdots & \omega^{-2(N-1)} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
1 & \omega^{-(N-1)} & \omega^{-2(N-1)} & \cdots & \omega^{-(N-1)^2}
\end{bmatrix}
}
$$

其中 $\omega \triangleq e^{j2\pi/N}$。$\mathbf U$ 是酉矩阵：$\mathbf U^H \mathbf U = \mathbf I_N$。其列向量构成 $\mathbb{C}^N$ 的一组标准正交基——这就是 **DFT 空间基**。

发射端和接收端分别构造 $\mathbf U_t \in \mathbb{C}^{N_t \times N_t}$ 和 $\mathbf U_r \in \mathbb{C}^{N_r \times N_r}$。

**4) 关键引理：任意导向矢量的 DFT 基展开**

这是连接连续角度与离散基的核心桥梁。对于**任意连续的** $\Omega$（不要求恰好落在网格上），其导向矢量 $\mathbf a(\Omega)$ 可用 DFT 基的线性组合**精确表示**：

$$
\boxed{
\mathbf a(\Omega) = \mathbf U\ \mathbf d(\Omega)
}
$$

其中 $\mathbf d(\Omega) \in \mathbb{C}^N$ 为 **DFT 系数向量**，满足 $\mathbf d(\Omega) = \mathbf U^H \mathbf a(\Omega)$。其第 $i$ 个分量具体为：

$$
\begin{aligned}
d_i(\Omega)
&= \bigl[\mathbf U^H \mathbf a(\Omega)\bigr]_i
= \mathbf a^H(\Omega_i)\,\mathbf a(\Omega) \\[4pt]
&= \frac{1}{N} \sum_{p=0}^{N-1}
e^{j2\pi p i/N}\cdot e^{-j2\pi p\Omega} \\[4pt]
&= \frac{1}{N} \sum_{p=0}^{N-1} e^{-j2\pi p(\Omega - i/N)} \\[4pt]
&= \frac{1}{N} \cdot \frac{1 - e^{-j2\pi N(\Omega - i/N)}}{1 - e^{-j2\pi(\Omega - i/N)}} \\[4pt]
&= e^{-j\pi(N-1)(\Omega - i/N)} \cdot
\frac{\sin\!\bigl[\pi N(\Omega - i/N)\bigr]}{N\sin\!\bigl[\pi(\Omega - i/N)\bigr]}
\end{aligned}
$$

这个系数是**Dirichlet 核**（周期 sinc 函数），以 $\Delta\Omega = \Omega - i/N$ 为自变量：

- 当 $\Omega$ **恰好落在**第 $i$ 个网格点时（$\Delta\Omega = 0$）：$d_i = 1$，其余所有 $d_j = 0$——导向矢量退化为 $\mathbf U$ 的单个列向量；
- 当 $\Omega$ **不在任何网格点上**时（off-grid）：$d_i(\Omega)$ 的幅值按 $\left|\frac{\sin(\pi N\Delta\Omega)}{N\sin(\pi\Delta\Omega)}\right|$ 衰减——能量**泄漏**到多个相邻的 DFT 基向量上。

这回答了本节的核心问题：**DFT 基是 $\mathbb{C}^N$ 中的一组完备正交基，任意连续角度的导向矢量都可以精确展开为 $N$ 个基向量的线性组合。** 当角度恰好对齐网格时，展开仅需一项（严格稀疏）；不对齐时则需要多项叠加（近似稀疏）。这就是角度域稀疏表示的本质。

**5) 几何多径信道的 DFT 基展开**

将上述展开代入 3.3.1 中的几何多径信道：

$$
\begin{aligned}
\mathbf H[m,k]
&= \sum_{\ell=0}^{L-1} \alpha_\ell e^{-j2\pi k\Delta f\tau_\ell}\,
\underbrace{\Bigl(\mathbf U_r\, \mathbf d_r(\Omega_{r,\ell})\Bigr)}_{\mathbf a_r(\Omega_{r,\ell})}
\underbrace{\Bigl(\mathbf U_t\, \mathbf d_t(\Omega_{t,\ell})\Bigr)^H}_{\mathbf a_t^H(\Omega_{t,\ell})} \\[4pt]
&= \mathbf U_r \left(
\sum_{\ell=0}^{L-1} \alpha_\ell e^{-j2\pi k\Delta f\tau_\ell}\,
\mathbf d_r(\Omega_{r,\ell})\,
\mathbf d_t^H(\Omega_{t,\ell})
\right) \mathbf U_t^H
\end{aligned}
$$

**信道矩阵完全可以用 DFT 基矩阵和对应的离散系数矩阵来表示。**

**6) 虚拟信道表示**

定义**角度域虚拟信道矩阵**（virtual channel representation）：

$$
\boxed{
\mathbf H^a[m,k] \triangleq \mathbf U_r^H\ \mathbf H[m,k]\ \mathbf U_t
}
$$

这是一个酉变换：$\|\mathbf H^a\|_F = \|\mathbf H\|_F$，且两者具有相同的互信息、条件数等。$\mathbf H^a$ 的第 $(i,j)$ 个元素的物理含义为：

$$
[\mathbf H^a]_{i,j}
= \mathbf a_r^H\!\left(\frac{i}{N_r}\right)\ \mathbf H[m,k]\ \mathbf a_t\!\left(\frac{j}{N_t}\right)
$$

即从**第 $j$ 个发射角度 bin**（$\Omega_t = j/N_t$）到**第 $i$ 个接收角度 bin**（$\Omega_r = i/N_r$）的复合信道增益。

从展开式看出，$\mathbf H^a$ 是各路径 $\mathbf d_r(\Omega_{r,\ell})\mathbf d_t^H(\Omega_{t,\ell})$ 的加权叠加。每条物理路径在 $\mathbf H^a$ 中产生一个以 $(i_\ell^*, j_\ell^*)$ 为中心、由 Dirichlet 核形状决定的**能量斑块**（blob）。

**向量化形式**：

$$
\mathbf h^a = \operatorname{vec}(\mathbf H^a)
= \bigl(\mathbf U_t^T \otimes \mathbf U_r^H\bigr)\,
\operatorname{vec}(\mathbf H)
$$

**7) on-grid 与 off-grid 的稀疏性对比**

| 情况 | AoA/AoD 位置 | $\mathbf H^a$ 的结构 |
|------|-------------|---------------------|
| 严格 on-grid | $\Omega_{t,\ell} \in \bigl\{\frac{j}{N_t}\bigr\},\ \Omega_{r,\ell} \in \bigl\{\frac{i}{N_r}\bigr\}$ | 严格稀疏：仅 $L$ 个非零元 |
| off-grid（实际） | $\Omega$ 不等于任何 $i/N$ | 近似稀疏：每路径能量按 Dirichlet 核泄漏到若干相邻 bin |

on-grid 假设下，每路径仅贡献一个非零元：

$$
\mathbf H^a = \sum_{\ell=0}^{L-1} \alpha_\ell\, e^{-j2\pi k\Delta f\tau_\ell}\,
\mathbf e_{i_\ell} \mathbf e_{j_\ell}^T
$$

其中 $\mathbf e_i$ 是第 $i$ 个标准基向量。这是压缩感知信道估计的理想工作模型。

实际系统永远是 off-grid，但：
- 大规模阵列（$N$ 大）$\to$ 角度网格分辨率更高 $\to$ 泄漏更集中在少数 bin 中 $\to$ 近似稀疏性更强；
- 在估计算法中可通过网格细化（grid refinement）或无网格（gridless）方法（如原子范数最小化）来克服 off-grid 偏差。

### 3.3.3 角度-时延域联合稀疏性

OFDM 系统的子载波维度提供了**时延分辨能力**。对 $\mathbf H^a[m,k]$ 沿子载波维度作 IDFT，得到**角度-时延域信道**：

$$
\boxed{
\mathbf H^a[m,n]
= \frac{1}{N_c}\sum_{k=0}^{N_c-1} \mathbf H^a[m,k]\ e^{j2\pi kn/N_c}
}
$$

或等价地，从空-时-频三维联合变换：

$$
\mathbf H^a[m,n] = \mathbf U_r^H\ \mathbf H[m,n]\ \mathbf U_t
$$

其中 $\mathbf H[m,n]$ 是时延 tap $n$ 处的 MIMO CIR 矩阵。

在角度-时延域中：
- 每一组 $(\text{AoA bin},\ \text{AoD bin},\ \text{delay tap})$ 对应一个候选传播路径；
- 物理路径数 $L$ 远小于角度-时延域的总网格点数 $N_r \times N_t \times N_{\mathrm{tap}}$；
- 因此 $\mathbf H^a[m,n]$ 在该三维空间中具有**显著的稀疏性**。

将所有维度堆叠为超向量：

$$
\mathbf h_{\mathrm{ad}} = \operatorname{vec}\bigl\{\mathbf H^a[m,n]\bigr\}
$$

$\mathbf h_{\mathrm{ad}}$ 的稀疏支撑集直接对应物理路径的角度和时延参数。这一性质是以下信道估计方法的理论基础：

- **压缩感知（Compressive Sensing）**：利用稀疏先验，从远少于网格点数的导频中恢复高维信道；
- **OMP / SOMP**：贪婪地逐路径估计 AoA/AoD/时延/增益；
- **LASSO / 稀疏贝叶斯学习**：通过 $\ell_1$ 正则化或分层贝叶斯先验促进稀疏解；
- **深度学习**：利用神经网络隐式学习角度-时延域的稀疏结构。


### 3.3.4 扩展：UPA（均匀平面阵）的角度域表示

实际 MIMO 系统（尤其是 5G NR 大规模 MIMO）常使用 UPA（Uniform Planar Array）而非 ULA。本节将 ULA 的 DFT 基分析推广到 UPA，核心思想是：**UPA 的导向矢量和 DFT 基矩阵是 ULA 对应部分的 Kronecker 积**。

**1) UPA 的阵列几何与导向矢量**

考虑 $N_h \times N_v$ 的 UPA（水平 $N_h$ 元，垂直 $N_v$ 元），总阵元数 $N = N_h N_v$。天线阵元的索引为 $(p, q)$，其中 $p = 0, \dots, N_h-1$ 为水平索引，$q = 0, \dots, N_v-1$ 为垂直索引。

在 3D 传播环境中，平面波由两个角度参数描述：
- $\theta$：俯仰角（elevation，相对于阵列法线 $\mathbf z$ 轴），$0 \le \theta \le \pi$；
- $\phi$：方位角（azimuth，在 $\mathbf{xy}$ 平面内的旋转角），$0 \le \phi < 2\pi$。

水平方向和垂直方向的波数分量分别为：

$$
k_x = \frac{2\pi}{\lambda} \sin\theta \cos\phi,\qquad
k_y = \frac{2\pi}{\lambda} \sin\theta \sin\phi
$$

阵元 $(p, q)$ 相对于参考阵元 $(0, 0)$ 的相位差为 $p \cdot k_x d_h + q \cdot k_y d_v$，其中 $d_h, d_v$ 分别为水平和垂直阵元间距。定义两个方向余弦：

$$
\boxed{
\Omega_h \triangleq \frac{d_h}{\lambda} \sin\theta \cos\phi,\qquad
\Omega_v \triangleq \frac{d_v}{\lambda} \sin\theta \sin\phi
}
$$

则导向矢量的 $(p, q)$ 分量为 $e^{-j2\pi(p\Omega_h + q\Omega_v)}$。将所有 $N_h N_v$ 个元素按列主序（column-major）堆叠为 $N \times 1$ 向量，UPA 的导向矢量可写为：

$$
\boxed{
\mathbf a_{\mathrm{UPA}}(\Omega_h, \Omega_v)
=
\mathbf a_h(\Omega_h) \otimes \mathbf a_v(\Omega_v)
}
$$

其中：
- $\mathbf a_h(\Omega_h) \in \mathbb{C}^{N_h}$ 是水平维 ULA 导向矢量：$[\mathbf a_h]_p = \frac{1}{\sqrt{N_h}} e^{-j2\pi p\Omega_h}$；
- $\mathbf a_v(\Omega_v) \in \mathbb{C}^{N_v}$ 是垂直维 ULA 导向矢量：$[\mathbf a_v]_q = \frac{1}{\sqrt{N_v}} e^{-j2\pi q\Omega_v}$；
- 归一化因子 $\frac{1}{\sqrt{N_h N_v}}$ 来自 $\otimes$ 积的结合。

**为什么需要 UPA？** ULA 只能在方位角方向（一维）上操控波束；当通信系统需要在**俯仰角和方位角两个维度上同时进行波束赋形**（例如高层建筑覆盖、无人机通信、卫星通信），ULA 力不从心。UPA 的天线分布在二维平面上，可同时控制水平与垂直两个维度的波束方向，实现**全空间 3D 波束操控**。这就是 5G NR 大规模 MIMO 大量采用 UPA 的根本原因。

**1b) UPA 的 2D 波束图样**

借由导向矢量的 Kronecker 积结构，UPA 的波束图样自然分解为水平和垂直两个一维波束图样的乘积。若以导向矢量 $\mathbf a_{\mathrm{UPA}}(\Omega_{h,0}, \Omega_{v,0})$ 做匹配滤波波束赋形，则在方向 $(\Omega_h, \Omega_v)$ 上的增益为

$
\boxed{
\begin{aligned}
G_{\mathrm{UPA}}\bigl((\Omega_{h,0}, \Omega_{v,0}),\ (\Omega_h, \Omega_v)\bigr)
&= \bigl|\mathbf a_{\mathrm{UPA}}^H(\Omega_{h,0}, \Omega_{v,0})\ \mathbf a_{\mathrm{UPA}}(\Omega_h, \Omega_v)\bigr|^2 \\[4pt]
&= \bigl|\mathbf a_h^H(\Omega_{h,0})\,\mathbf a_h(\Omega_h)\bigr|^2\ \cdot\
\bigl|\mathbf a_v^H(\Omega_{v,0})\,\mathbf a_v(\Omega_v)\bigr|^2 \\[4pt]
&= G_h(\Omega_{h,0}, \Omega_h)\ \cdot\ G_v(\Omega_{v,0}, \Omega_v),
\end{aligned}
}
$

其中 $G_h$ 和 $G_v$ 即为上一节定义的 1D Dirichlet 核平方波束图样（分别以 $N_h$ 和 $N_v$ 为天线数）。这种**可分离性**（separability）是 UPA 角度域分析的核心便利：2D 波束是水平和垂直两个 1D 波束的张量积，所有关于主瓣宽度、零点位置、off-grid 泄漏的结论均可逐维沿用。

2D 波束图样的关键特征：

- **主瓣**：位于 $(\Omega_h, \Omega_v) = (\Omega_{h,0}, \Omega_{v,0})$，峰值增益为 $1$；
- **零点网格**：在 $\Omega_h = \Omega_{h,0} \pm k/N_h$ 或 $\Omega_v = \Omega_{v,0} \pm \ell/N_v$（$k = 1,\dots,N_h-1$，$\ell = 1,\dots,N_v-1$）处出现零点，形成规则的"零点栅格"；
- **旁瓣**：第一旁瓣电平约为 $-13$ dB（与 1D 一致，由 Dirichlet 核的旁瓣结构决定）；
- **角度分辨率**：水平维分辨率 $\propto 1/N_h$，垂直维分辨率 $\propto 1/N_v$。因此 $32 \times 32$ 的 UPA 不仅在总阵元数上等价于 $1024$ 元 ULA，还能在两个维度上同时提供高分辨率。

**2) UPA 的 DFT 基矩阵**

对水平维：将 $\Omega_h$ 均匀离散化为 $N_h$ 个网格点 $\Omega_{h,i} = i/N_h\ (i = 0, \dots, N_h-1)$，构造 $N_h \times N_h$ 酉 DFT 矩阵 $\mathbf U_h$。

对垂直维：将 $\Omega_v$ 均匀离散化为 $N_v$ 个网格点 $\Omega_{v,j} = j/N_v\ (j = 0, \dots, N_v-1)$，构造 $N_v \times N_v$ 酉 DFT 矩阵 $\mathbf U_v$。

**UPA 的 DFT 基矩阵**为两个维度的 Kronecker 积：

$$
\boxed{
\mathbf U_{\mathrm{UPA}} = \mathbf U_h \otimes \mathbf U_v \ \in \mathbb{C}^{N_h N_v \times N_h N_v}
}
$$

这实质上是 **2D-DFT 矩阵**：

- $\mathbf U_{\mathrm{UPA}}$ 的每一列对应一对离散角度 bin $(\Omega_{h,i}, \Omega_{v,j})$，即 2D 角度域的一个网格点；
- 两个维度共 $N_h \times N_v$ 个角度 bin，和阵列总阵元数一致。

**3) UPA 的虚拟信道表示**

MIMO 系统中发射端和接收端均可以是 UPA。设发射端为 $N_{t,h} \times N_{t,v}$ 的 UPA，接收端为 $N_{r,h} \times N_{r,v}$ 的 UPA。信道矩阵 $\mathbf H \in \mathbb{C}^{N_r N_{r,h} N_{r,v} \times N_t N_{t,h} N_{t,v}}$（这里假设为平坦衰落或单子载波）。角度域虚拟信道为：

$$
\boxed{
\mathbf H^a = \bigl(\mathbf U_{r,h} \otimes \mathbf U_{r,v}\bigr)^H\ \mathbf H\ \bigl(\mathbf U_{t,h} \otimes \mathbf U_{t,v}\bigr)
}
$$

注意利用 Kronecker 积的性质 $(\mathbf A \otimes \mathbf B)^H = \mathbf A^H \otimes \mathbf B^H$：

$$
\mathbf H^a = \bigl(\mathbf U_{r,h}^H \otimes \mathbf U_{r,v}^H\bigr)\ \mathbf H\ \bigl(\mathbf U_{t,h} \otimes \mathbf U_{t,v}\bigr)
$$

$\mathbf H^a$ 的维度也是 $N_{r,h}N_{r,v} \times N_{t,h}N_{t,v}$，每个元素 $(I, J)$ 对应一组 4D 角度 bin：

$$
(\Omega_{r,h,i}, \Omega_{r,v,j}) \to (\Omega_{t,h,i'}, \Omega_{t,v,j'})
$$

**4) on-grid 下的 2D 稀疏性**

若所有物理路径的 AoA/AoD 恰好落在 2D 角度网格上，则 $\mathbf H^a$ 中仅有 $L$ 个非零元——每条路径贡献一个 4D 角度 bin 上的复增益。路径数 $L$ 远小于总角度 bin 数 $N_{r,h}N_{r,v}N_{t,h}N_{t,v}$，因此 $\mathbf H^a$ 是**严格稀疏的 2D 信道张量**。

**5) 与 ULA 对比**

| 特性 | ULA | UPA |
|------|-----|-----|
| 角度参数 | 1D: $\Omega$（$\sin\theta$） | 2D: $(\Omega_h, \Omega_v)$ = $(\sin\theta\cos\phi,\ \sin\theta\sin\phi)$ |
| 导向矢量 | $\mathbf a(\Omega)$ | $\mathbf a_h(\Omega_h) \otimes \mathbf a_v(\Omega_v)$ |
| DFT 基矩阵 | $\mathbf U$（1D-DFT） | $\mathbf U_h \otimes \mathbf U_v$（2D-DFT） |
| 角度 bin 总数 | $N$ | $N_h \times N_v$ |
| 波束赋形维度 | 仅方位角方向 | 方位角 + 俯仰角（全空间 3D） |
| 波束图样结构 | 1D Dirichlet 核平方 | 两 1D Dirichlet 核平方的乘积 |
| 稀疏域示意图 | 1D 角度谱 | 2D 角度图（俯仰-方位） |
| 典型应用 | 3G/4G/Wi-Fi | 5G NR 大规模 MIMO |

**6) 张量视角——爱因斯坦求和约定**

UPA 的信道具有天然的**多索引结构**：发射端天线由水平维索引 $h$ 和垂直维索引 $v$ 共同标记，接收端同理。与其将所有天线强制展平为一维向量（丢失行列对应关系），不如直接用**张量（Tensor）**语言来书写。

**上/下指标记号**：采用爱因斯坦求和约定（Einstein Summation Convention）——重复出现的一对上/下指标自动隐含求和。

将发射信号记为 2-指标张量 $X^{hv}$（第 $h$ 列、第 $v$ 行天线上的复信号），接收信号为 $Y^{h'v'}$，则空域 MIMO 输入-输出关系可写为 4-指标形式：

$
Y^{h'v'} = H^{h'v'}_{\ \ \ hv}\ X^{hv},
$

其中：
- $H^{h'v'}_{\ \ \ hv}$ 是从发射天线 $(h, v)$ 到接收天线 $(h', v')$ 的复信道增益；
- 按爱因斯坦约定，上标 $hv$ 重复出现，对 $h = 0,\dots,N_{t,h}-1$ 和 $v = 0,\dots,N_{t,v}-1$ 自动求和；
- $H^{h'v'}_{\ \ \ hv}$ 既可看作 $N_{r,h}N_{r,v} \times N_{t,h}N_{t,v}$ 的矩阵，也可看作 $N_{r,h} \times N_{r,v} \times N_{t,h} \times N_{t,v}$ 的四维数组。

**空间标签的张量形式**：水平维的空间标签 $\mathbf e_h(\Omega_h)$ 和垂直维的空间标签 $\mathbf e_v(\Omega_v)$ 的张量积构成 UPA 的空间标签：

$
e^{hv}(\Omega_h, \Omega_v) = e_h^h(\Omega_h)\ e_v^v(\Omega_v),
$

其中 $e^{hv}$ 表示具体天线 $(h, v)$ 处的空间标签值（即 $\mathbf a_{\mathrm{UPA}}$ 的第 $(h, v)$ 个分量）。

两个空间标签的内积在爱因斯坦约定下极为简洁——共轭转置对应上下标颠倒：

$
\langle \mathbf a_1, \mathbf a_2 \rangle = \bar e_1^{hv}\ e_2^{hv},
$

这里 $\bar e_1^{hv}$ 已将共轭隐含（上标颠倒记号）。

**角度域的张量变换**：空域向量 $X^{hv}$ 到角度域向量 $x^{ij}$ 的变换为 2D-IDFT（在张量记号下同样优美）：

$
x^{ij} = X^{hv}\ e_h^h(\Omega_{h,i})\ e_v^v(\Omega_{v,j}),
$

其中 $\Omega_{h,i} = i/N_h$，$\Omega_{v,j} = j/N_v$。这正是将空域信号按水平和垂直两个维度的 DFT 基分别展开。

相应地，角度域信道张量为：

$
\boxed{
h^{i'j'}_{\ \ \ ij} = \bar e_{r,h}^{i'}\ \bar e_{r,v}^{j'}\ H^{h'v'}_{\ \ \ hv}\ e_{t,h}^{i}\ e_{t,v}^{j}
}
$

$h^{i'j'}_{\ \ \ ij}$ 的物理意义为：从发射端水平角度 bin $i$、垂直角度 bin $j$ 发出，在接收端水平角度 bin $i'$、垂直角度 bin $j'$ 被捕获的路径的等效复增益。展开后即为 2D-DFT：

$
h^{i'j'}_{\ \ \ ij} = \frac{1}{N_{t,h}N_{t,v}N_{r,h}N_{r,v}}
\sum_{h,v,h',v'}
H^{h'v'}_{\ \ \ hv}\
e^{j2\pi(h'i'/N_{r,h} + v'j'/N_{r,v})}\
e^{-j2\pi(hi/N_{t,h} + vj/N_{t,v})}.
$

**张量记号的便利**：完整的 MIMO 空域到角域变换只需一行：

$
Y^{h'v'} = H^{h'v'}_{\ \ \ hv}\ X^{hv}
\quad\longleftrightarrow\quad
y^{i'j'} = h^{i'j'}_{\ \ \ ij}\ x^{ij}.
$

这种写法天然反映了 UPA 的可分离结构——水平和垂直维度在指标上独立演变，无需人为展平和重组矩阵，是所有 Kronecker 积分解的自然语言。

UPA 的 Kronecker 结构使其角度域分析完全分解为水平和垂直两个独立维度，UPA 的 2D-DFT 基仍是酉矩阵。本节 3.3.2 中关于 Dirichlet 核展开、off-grid 泄漏等所有结论均可逐维推广到 UPA。



## 3.4 3GPP TR 38.901 CDL/TDL 模型

3GPP TR 38.901 是 5G/5G-Advanced 链路级和系统级仿真的标准化信道模型，适用频率 0.5–100 GHz，支持 0–1000 km/h 的移动速度。它整合了 3GPP 此前多个信道模型（3D-UMa/UMi SCM、IMT-Advanced 等），基本上在每个协议版本都有更新，不断扩充。

### 3.4.1 模型架构

38.901 将信道划分为两层：

| 层次 | 内容 | 更新频率 |
|------|------|----------|
| **大尺度（Large-Scale）** | 路径损耗、阴影衰落、穿透损耗、LOS/NLOS 状态判定 | 每 drop（通常数十毫秒） |
| **小尺度（Small-Scale）** | 簇、射线、时延、功率、角参数、极化、多普勒 | 每 snapshot（OFDM 符号级） |

仿真流程可概括为三个步骤：

1. **Drop 初始化**：按场景（UMa/UMi/RMa/InH 等）确定大尺度参数，判定 LOS/NLOS；
2. **生成小尺度参数**：从表中查簇数 $N$、每簇射线数 $Q$，生成簇的时延、功率、到达角（AoA/ZoA）、离开角（AoD/ZoD）、交叉极化比（XPR）等；
3. **按时间步进**：在每个采样时刻，由射线角度和 UE 速度计算多普勒相位，乘以天线方向图和阵列相位，叠加所有簇和射线。

### 3.4.2 CDL 信道冲激响应

CDL（Clustered Delay Line）显式建模每条射线的角度参数，因此能精确刻画 MIMO 阵列的空间特性。第 $n$ 个簇（$n = 1,\dots,N$）的第 $q$ 条射线（$q = 1,\dots,Q$）对发射天线 $s$ 到接收天线 $u$ 的贡献为：

$
\begin{aligned}
h_{u,s,n}^{\mathrm{NLOS}}(t)
&=
\sqrt{\frac{P_n}{Q}}
\sum_{q=1}^{Q}
\mathbf F_{r,u}^T(\theta_{n,q,\mathrm{ZoA}},\phi_{n,q,\mathrm{AoA}})
\ \mathbf \Phi_{n,q}\ 
\mathbf F_{t,s}(\theta_{n,q,\mathrm{ZoD}},\phi_{n,q,\mathrm{AoD}}) \\
&\quad \times
\exp\!\Bigl(j\frac{2\pi}{\lambda_0}\,\hat{\mathbf r}_{r,n,q}^T\mathbf d_{r,u}\Bigr)
\exp\!\Bigl(j\frac{2\pi}{\lambda_0}\,\hat{\mathbf r}_{t,n,q}^T\mathbf d_{t,s}\Bigr)
\exp\!\Bigl(j\frac{2\pi}{\lambda_0}\,\hat{\mathbf r}_{r,n,q}^T\mathbf v\,t\Bigr).
\end{aligned}
$

若为 LOS 场景，还需叠加直射径分量：

$
h_{u,s}^{\mathrm{LOS}}(t) =
\sqrt{\frac{1}{K_R+1}}\,
\mathbf F_{r,u}^T(\theta_{\mathrm{LOS,ZoA}},\phi_{\mathrm{LOS,AoA}})
\begin{bmatrix} 1 & 0 \\ 0 & -1 \end{bmatrix}
\mathbf F_{t,s}(\theta_{\mathrm{LOS,ZoD}},\phi_{\mathrm{LOS,AoD}}) \\
\times
\exp\!\Bigl(j\frac{2\pi}{\lambda_0}\,\hat{\mathbf r}_{r,\mathrm{LOS}}^T\mathbf d_{r,u}\Bigr)
\exp\!\Bigl(j\frac{2\pi}{\lambda_0}\,\hat{\mathbf r}_{t,\mathrm{LOS}}^T\mathbf d_{t,s}\Bigr)
\exp\!\Bigl(j\frac{2\pi}{\lambda_0}\,\hat{\mathbf r}_{r,\mathrm{LOS}}^T\mathbf v\,t\Bigr),
$

完整冲激响应为 LOS 分量与小尺度 NLOS 分量之和；$K_R$ 为 Rician K 因子（线性）。

### 3.4.3 参数逐一解释

下表逐一说明 CDL 公式中各符号的物理含义与取值范围：

| 符号 | 含义 | 说明 |
|------|------|------|
| $N$ | 簇数 | CDL-A/B/C 含 23–24 个簇，CDL-D/E 含 13–14 个簇（含 LOS 簇） |
| $Q$ | 每簇射线数 | 标准取 $Q = 20$。这 20 条射线在簇内偏离簇中心角 $\pm c_{\mathrm{ASD/ASA/ZSD/ZSA}}$ 范围内分布 |
| $P_n$ | 第 $n$ 个簇的归一化功率 | $\sum_n P_n = 1$；首簇功率最大，后续按指数衰减。LOS 时 LOS 簇功率由 K 因子单独分配 |
| $\tau_n$ | 第 $n$ 个簇的时延 | 由归一化时延表给出，仿真时按期望时延扩展（DS）缩放 |
| $\theta, \phi$ | 俯仰角（ZoD/ZoA）、方位角（AoD/AoA） | Zo = Zenith（天顶角），A = Azimuth。四个角度各由簇中心角 + 射线偏移构成 |
| $c_{\mathrm{ASD}},\ c_{\mathrm{ASA}}$ | 离开/到达方位角扩展因子 | 控制簇内射线在方位角上的散布，由标准预定义 |
| $c_{\mathrm{ZSD}},\ c_{\mathrm{ZSA}}$ | 离开/到达天顶角扩展因子 | 同上，用于俯仰维度 |
| $\mathbf F_{t,s},\ \mathbf F_{r,u}$ | 发射天线 $s$ / 接收天线 $u$ 的复场方向图 | 2×1 或 1×2 向量，包含极化（$\theta$ 和 $\phi$ 分量）。若不考虑极化天线，可取为 $[1, 1]^T$ 或标量 1 |
| $\mathbf \Phi_{n,q}$ | 簇 $n$ 射线 $q$ 的极化耦合矩阵 | $2 \times 2$ 矩阵，含交叉极化比 $\kappa$（XPR），描述 $\theta/\phi$ 分量间的功率转移 |
| $\lambda_0$ | 载波波长 | $\lambda_0 = c / f_c$ |
| $\hat{\mathbf r}_{t,n,q}$ | 发射方向单位矢量 | $[\sin\theta\cos\phi,\ \sin\theta\sin\phi,\ \cos\theta]^T$，由 ZoD 和 AoD 确定 |
| $\hat{\mathbf r}_{r,n,q}$ | 接收方向单位矢量 | 同上，由 ZoA 和 AoA 确定 |
| $\mathbf d_{t,s}$ | 发射天线 $s$ 的位置矢量 | 相对于阵列参考点（通常为原点），在天线阵列坐标系中 |
| $\mathbf d_{r,u}$ | 接收天线 $u$ 的位置矢量 | 同上 |
| $\mathbf v$ | UE 速度矢量 | $\mathbf v = v \cdot [\sin\theta_v\cos\phi_v,\ \sin\theta_v\sin\phi_v,\ \cos\theta_v]^T$，其中 $(\theta_v, \phi_v)$ 为速度方向角 |
| $K_R$ | Rician K 因子（线性） | LOS 场景专用；NLOS 场景 $K_R = 0$（即无直射径项） |

**核心要点**：CDL 公式本质上是 §3.3 的几何多径模型的**标准化实现**——它同样用导向矢量（阵列相位项）叠加各路径，但所有路径参数（角度、时延、功率、XPR）都来自 TR 38.901 的预定义表格，而非用户自行设定。

### 3.4.4 五个预定义 CDL 模型

TR 38.901 表 7.7.1-1 至 7.7.1-5 定义了五种 CDL 模型，覆盖从散射丰富的 NLOS 到强 LOS 的典型场景。以下给出关键参数，更完整的每簇时延/功率/角度表请参考规范原文。

| 参数 | CDL-A | CDL-B | CDL-C | CDL-D | CDL-E |
|------|-------|-------|-------|-------|-------|
| **传播条件** | NLOS | NLOS | NLOS | LOS | LOS |
| **簇数 $N$** | 23 | 23 | 24 | 13 | 14 |
| **每簇射线 $Q$** | 20 | 20 | 20 | 20 | 20 |
| **归一化时延扩展** | ~0.44 | ~1.07 | ~2.46 | ~0.29 | ~0.30 |
| **角度扩展** | **大** | **中等** | **中-大** | **较小** | **最小** |
| ASD（离开方位角）/° | ~22–26 | ~10–11 | ~15–18 | ~4–6 | ~3–4 |
| ASA（到达方位角）/° | ~37–44 | ~19–22 | ~26–29 | ~7–9 | ~5–7 |
| ZSD（离开天顶角）/° | ~6 | ~6 | ~6 | ~3 | ~3 |
| ZSA（到达天顶角）/° | ~8 | ~8 | ~8 | ~3–4 | ~3–4 |
| **典型 K 因子 / dB** | — | — | — | ~11 | ~20 |
| **XPR（交叉极化比）/ dB** | 11 | 9 | 7 | 11 | 8 |
| **典型场景** | 城区宏小区 NLOS | 城区微小区 NLOS | 郊区/农村 NLOS | 城区宏小区 LOS | 郊区/开阔地 LOS |

**如何解读**：

- **CDL-A**：角度扩展最大、散射最丰富，MIMO 信道矩阵趋向满秩，空间复用潜力高；
- **CDL-B**：中等角度扩展，是大量仿真中默认的 NLOS 基线；
- **CDL-C**：时延扩展最大（多径更长），频率选择性最强，对 OFDM 导频密度要求更高；
- **CDL-D**：弱 LOS（$K_R \approx 10$–$11$ dB），角度扩展小，信道矩阵趋向低秩，波束赋形增益显著；
- **CDL-E**：强 LOS（$K_R \approx 20$ dB），角度扩展极小，信道接近秩 1，空间复用层数远少于天线数。

### 3.4.5 时延与角度的缩放

表中的归一化参数需要按期望的时延扩展（DS, Delay Spread）和角度扩展（AS, Angular Spread）进行缩放，以匹配具体场景：

- **时延缩放**：$\tau_{n,\mathrm{scaled}} = \tau_{n,\mathrm{norm}} \cdot \mathrm{DS}_{\mathrm{desired}}$，其中 $\mathrm{DS}_{\mathrm{desired}}$ 是大尺度参数（如 UMa NLOS 场景典型 DS = 129 ns）；
- **角度缩放**：四个角度（ASD、ASA、ZSD、ZSA）各自按期望值缩放簇内射线偏移 $c_{\mathrm{ASD}}$ 等，确保最终角度扩展匹配场景。

这种"归一化表 + 按场景缩放"的设计使同一组 CDL 模型可适配 UMa、UMi、RMa、InH 等多种场景，只需改变大尺度参数而非更换整个模型。

### 3.4.6 TDL 模型

TDL（Tapped Delay Line）是 CDL 的简化版。它**丢弃所有角度和极化信息**，仅保留时延和平均功率：

| 模型 | 传播条件 | 特点 |
|------|----------|------|
| TDL-A | NLOS | 23 taps，大时延扩展 |
| TDL-B | NLOS | 23 taps，中等 |
| TDL-C | NLOS | 24 taps，最大时延扩展 |
| TDL-D | LOS | 13 taps + LOS tap |
| TDL-E | LOS | 14 taps + LOS tap，强 LOS |

TDL 适合以下场景：
- SISO / SIMO 仿真，不需要天线阵列几何；
- 只关心频率选择性，不关心空间选择性；
- 快速原型验证，省去天线方向图和阵列相位计算。

TDL 的每个 tap 独立按 Rayleigh（NLOS tap）或 Rician（LOS tap）衰落，本质上退化为 §3.1 的独立 tap 模型，但参数来自标准化表格。

### 3.4.7 LOS/NLOS 状态判定

TR 38.901 中，LOS 概率是 UE-BS 距离的函数，不同场景有不同的经验公式。以 UMa（城市宏小区）为例：

$
P_{\mathrm{LOS}}(d_{2D}) =
\begin{cases}
1, & d_{2D} \le 18\ \mathrm{m}, \\[4pt]
\displaystyle\Bigl(\frac{18}{d_{2D}} + \exp\!\bigl(-\frac{d_{2D}}{63}\bigr)\bigl(1 - \frac{18}{d_{2D}}\bigr)\Bigr), & d_{2D} > 18\ \mathrm{m}.
\end{cases}
$

仿真时，每个 drop 按 $P_{\mathrm{LOS}}$ 抽样一次，决定使用 CDL-D/E（LOS）还是 CDL-A/B/C（NLOS）。

### 3.4.8 与前几节模型的对比

| 特性 | §3.1 i.i.d. tap | §3.2 Kronecker | §3.3 几何 AoA/AoD | §3.4 CDL/TDL |
|------|-----------------|----------------|-------------------|--------------|
| 角度信息 | 无 | 统计（相关矩阵） | 任意（用户指定） | 标准化表格 |
| 极化 | 无 | 无 | 可选 | 有（XPR + 极化场方向图） |
| 天线方向图 | 无 | 无 | 可选 | 有（$\mathbf F_{t,s}, \mathbf F_{r,u}$） |
| 标准化程度 | 无 | 低 | 无 | 高（3GPP 规范） |
| 仿真复杂度 | 最低 | 低-中 | 中 | 中-高 |
| 适用场景 | 算法验证 | 相关分析 | 波束/稀疏研究 | 标准对齐仿真 |

**选择建议**：若目标是和 5G NR 链路级评估对齐，或需要完整的天线方向图、极化效应、3D 空间信息，使用 CDL；若只需快速验证信道估计/预编码算法正确性，§3.1–§3.3 的简化模型已足够。

> 导频设计与信道估计的内容已移至 [MIMOOFDMChannelEstimation.md](MIMOOFDMChannelEstimation.md)。

# 4 仿真建模流程

一个可复现实验的 MIMO-OFDM 信道仿真通常按以下流程组织。

### Step 1：设置系统参数

给定：

$$
N_t,N_r,N,M,\Delta f,T_{cp},f_c,v_{\max}.
$$

计算：

$$
T_s=\frac{1}{N\Delta f},\qquad
T_{\mathrm{sym}}=\frac{1}{\Delta f}+T_{cp},\qquad
f_D=\frac{v_{\max}}{\lambda}.
$$

### Step 2：设置 PDP 和 tap

给定：

$$
\{P_\ell,\tau_\ell\}_{\ell=0}^{L-1}.
$$

若使用离散 tap，令：

$$
n_\ell=\mathrm{round}(\tau_\ell/T_s).
$$

### Step 3：生成空间信道

可选模型：

1. i.i.d. Rayleigh/Rician；
2. Kronecker 空间相关；
3. 几何 AoA/AoD 模型；
4. 3GPP CDL/TDL。

以 Kronecker 为例：

$$
\mathbf G_\ell[m]=
\mathbf R_{r,\ell}^{1/2}\mathbf W_\ell[m]\mathbf R_{t,\ell}^{1/2}.
$$

### Step 4：生成时间相关

可用：

- Jakes/Clarke sum-of-sinusoids；
- AR(1) 近似；
- BEM；
- 3GPP 中由射线方向和速度产生的多普勒相位。

### Step 5：构造 CFR

$$
\mathbf H[m,k]=
\sum_{\ell=0}^{L-1}
\sqrt{P_\ell}\mathbf G_\ell[m]
e^{-j2\pi k\Delta f\tau_\ell}.
$$

### Step 6：发送导频和数据

对每个 $(m,k)$：

$$
\mathbf y[m,k]=\mathbf H[m,k]\mathbf x[m,k]+\mathbf w[m,k].
$$

### Step 7：估计、插值和均衡

1. 在导频位置做 LS 或 LMMSE；
2. 在时频空维度插值；
3. 对每个子载波做 ZF/MMSE/ML/MAP 检测；
4. 统计 NMSE、BER、EVM、容量或吞吐量。

常用 NMSE：

$$
\boxed{
\mathrm{NMSE}
=
\frac{
\sum_{m,k}\|\hat{\mathbf H}[m,k]-\mathbf H[m,k]\|_F^2
}{
\sum_{m,k}\|\mathbf H[m,k]\|_F^2
}
}
$$

# 5 常见模型选择建议

| 目标 | 推荐模型 | 说明 |
|------|----------|------|
| 验证 LS/MMSE 公式 | i.i.d. Rayleigh tap | 简单、可控、便于复现 |
| 研究空间相关影响 | Kronecker/Weichselberger | 可显式调节相关矩阵 |
| 研究阵列、角度、波束 | 几何 AoA/AoD 模型 | 能解释秩、角度稀疏和波束训练 |
| 对齐 5G NR 链路级仿真 | 3GPP TR 38.901 CDL/TDL | 标准化、参数完整 |
| 大规模 MIMO/mmWave 估计 | 角度-时延稀疏模型 | 适合压缩感知和深度学习估计 |
| 高速移动双选择性 | Jakes/AR/BEM + MIMO 空间模型 | 需要显式时间相关 |

经验上：

- 若只是调通信道估计算法，先用 i.i.d. Rayleigh；
- 若要验证 MIMO 特有现象，至少加入空间相关或几何角度；
- 若要写接近通信标准的仿真，使用 3GPP CDL；
- 若阵列规模很大，不要只在天线域估计，应考虑角度域或低秩/稀疏先验。

# 6 References

1. D. Tse and P. Viswanath, *Fundamentals of Wireless Communication*, Cambridge University Press, 2005. Chapter 7: MIMO I, Spatial Multiplexing and Channel Modeling. https://web.stanford.edu/~dntse/wireless_book.html
2. A. M. Sayeed, "Deconstructing Multiantenna Fading Channels," *IEEE Transactions on Signal Processing*, vol. 50, no. 10, pp. 2563-2579, 2002. https://minds.wisconsin.edu/handle/1793/9386
3. D. Shiu, G. J. Foschini, M. J. Gans, and J. M. Kahn, "Fading correlation and its effect on the capacity of multielement antenna systems," *IEEE Transactions on Communications*, vol. 48, no. 3, pp. 502-513, 2000.
4. 3GPP TR 38.901, *Study on channel model for frequencies from 0.5 to 100 GHz*. ETSI TR 138 901 V19.2.0, 2026-02. https://www.etsi.org/deliver/etsi_tr/138900_138999/138901/19.02.00_60/tr_138901v190200p.pdf
5. 3D MIMO channel modeling introduction. https://zhuanlan.zhihu.com/p/664298097
6. The spatial correlation of MIMO channels. https://zhuanlan.zhihu.com/p/721013608
7. 逸风晴, "MIMO 通信的角域表示," 知乎专栏「杂谈通信中的线性代数」, 2021. https://zhuanlan.zhihu.com/p/363018716
