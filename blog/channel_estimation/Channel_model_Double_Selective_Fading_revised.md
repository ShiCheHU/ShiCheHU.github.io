# Channel model -- Doubly Selective Fading Channel

## 目录

- [1 概述](#1-概述)
- [2 系统参数与一般形式](#2-系统参数与一般形式)
  - [2.1 基本参数](#21-基本参数)
  - [2.2 连续域 CIR](#22-连续域-cir)
  - [2.3 连续域 CFR](#23-连续域-cfr)
  - [2.4 离散域表示](#24-离散域表示)
- [3 双选择性衰落信道建模](#3-双选择性衰落信道建模)
  - [3.1 小尺度衰落分类](#31-小尺度衰落分类)
  - [3.2 双选择性信道的统一模型](#32-双选择性信道的统一模型)
  - [3.3 块衰落与相干时间内慢变](#33-块衰落与相干时间内慢变)
- [4 路径增益过程的建模](#4-路径增益过程的建模)
  - [4.1 Jakes/Clarke 建模](#41-jakesclarke-建模)
  - [4.2 自回归建模](#42-自回归建模)
  - [4.3 基展开模型（BEM）](#43-基展开模型bem)
- [5 Rayleigh / Rician 双选择性信道](#5-rayleigh--rician-双选择性信道)
  - [5.1 双选择性 Rayleigh 衰落信道](#51-双选择性-rayleigh-衰落信道)
  - [5.2 双选择性 Rician 衰落信道](#52-双选择性-rician-衰落信道)
  - [5.3 生成步骤](#53-生成步骤)
- [6 特殊配置与退化情形](#6-特殊配置与退化情形)
- [7 参考](#7-参考)

# 1 概述

双选择性衰落信道（doubly selective fading channel）同时体现：

- **时延扩展**导致的频率选择性；
- **多普勒扩展**导致的时间选择性。

对于 OFDM 系统而言，双选择性信道是最一般的建模情形。平坦衰落、频率选择性慢衰落、时间选择性平坦衰落等常见场景，都可以看作该模型在特定参数设定下的退化情形。

本文首先给出连续域和离散域的统一表示，然后介绍双选择性信道的常见建模方法，以及 Rayleigh / Rician 两类典型统计模型。

# 2 系统参数与一般形式

## 2.1 基本参数

- 路径数：$L$
- 最大时延：$\tau_{\max}$
- 载波频率：$f_c$
- 最大移动速度：$v_{\max}$
- 最大多普勒频移：
  $$
  \nu_{\max}=\frac{v_{\max}}{\lambda},\qquad \lambda=\frac{c}{f_c}
  $$
- 子载波间隔：$\Delta f$
- 每个 OFDM 符号的子载波数：$N$
- 有效 OFDM 符号持续时间：
  $$
  T_u=\frac{1}{\Delta f}
  $$
- 循环前缀长度：$T_{cp}$
- OFDM 符号总时长：
  $$
  T_{\mathrm{sym}}=T_u+T_{cp}
  $$
- OFDM 符号数：$M$，其中 $m=0,1,\dots,M-1$

为避免 ISI，通常要求：

$$
T_{cp}>\tau_{\max}.
$$

记第 $m$ 个 OFDM 符号对应的时刻为：

$$
t_m=mT_{\mathrm{sym}}.
$$

## 2.2 连续域 CIR

一般将信道冲激响应（CIR）表示为：

$$
\boxed{
h(t,\tau)=\sum_{\ell=0}^{L-1}\alpha_\ell(t)\,\delta\bigl(\tau-\tau_\ell(t)\bigr)
}
$$

其中：

- $t$ 表示观测时间，用于描述多普勒引起的时间变化；
- $\tau$ 表示时延变量，用于描述多径传播时延；
- $\alpha_\ell(t)$ 为第 $\ell$ 条路径的时变复增益；
- $\tau_\ell(t)$ 为第 $\ell$ 条路径的时变传播时延。

通常在一个 OFDM 符号的持续时间内，可认为路径时延变化较慢，因此近似写成：

$$
\tau_\ell(t)\approx \tau_\ell.
$$

于是得到：

$$
\boxed{
h(t,\tau)=\sum_{\ell=0}^{L-1}\alpha_\ell(t)\,\delta(\tau-\tau_\ell)
}
$$

如果进一步采用单频多普勒近似，可令

$$
\alpha_\ell(t)=\beta_\ell e^{j2\pi \nu_\ell t},
$$

其中 $\nu_\ell$ 为第 $\ell$ 条路径的多普勒频移，则有

$$
\boxed{
h(t,\tau)=\sum_{\ell=0}^{L-1}\beta_\ell e^{j2\pi\nu_\ell t}\,\delta(\tau-\tau_\ell)
}
$$

该表达式适合描述少量镜面路径的相位旋转。对于散射丰富的信道，更常见的做法是将 $\alpha_\ell(t)$ 建模为一个复随机过程，后文将统一写成 $\sqrt{P_\ell}g_\ell(t)$。

## 2.3 连续域 CFR

信道频率响应（CFR）是 CIR 关于 $\tau$ 的傅里叶变换：

$$
\begin{aligned}
H(t,f)
&=\int h(t,\tau)e^{-j2\pi f\tau}d\tau \\
&=\sum_{\ell=0}^{L-1}\alpha_\ell(t)e^{-j2\pi f\tau_\ell}.
\end{aligned}
$$

因此，时延扩展决定了频域选择性，而 $\alpha_\ell(t)$ 的时间变化决定了时间选择性。

## 2.4 离散域表示

在 OFDM 中，通常把连续模型离散为：

- OFDM 符号索引 $m$；
- 离散时延 tap 索引 $n$；
- 子载波索引 $k$。

离散 CIR 记为

$$
h[m,n],
$$

表示第 $m$ 个 OFDM 符号时刻、第 $n$ 个离散时延 tap 上的信道。

对应的二维 CFR 为

$$
H[m,k],
$$

表示第 $m$ 个 OFDM 符号时刻、第 $k$ 个子载波上的信道。两者满足 DFT / IDFT 关系：

$$
\boxed{
H[m,k]=\sum_{n=0}^{N-1}h[m,n]e^{-j\frac{2\pi}{N}kn}
}
$$

$$
\boxed{
h[m,n]=\frac{1}{N}\sum_{k=0}^{N-1}H[m,k]e^{j\frac{2\pi}{N}kn}
}
$$

# 3 双选择性衰落信道建模

## 3.1 小尺度衰落分类

小尺度衰落通常可按如下两个维度分类：

| 分类依据 | 类型 | 物理成因 |
|---|---|---|
| 时延扩展 vs. 信号带宽 | 平坦衰落 / 频率选择性衰落 | 多径时延 |
| 相干时间 vs. 符号周期 | 慢衰落 / 快衰落 | 多普勒扩展 |

这里需要特别注意：

- **慢衰落**：$T_c \gg T_{\mathrm{sym}}$，信道在一个 OFDM 符号内变化很小；
- **快衰落**：$T_c \lesssim T_{\mathrm{sym}}$，信道在一个 OFDM 符号内变化不可忽略，可能引入 ICI。

即为：

$$
\boxed{
T_c \gg T_{\mathrm{sym}} \Rightarrow \text{慢衰落},
\qquad
T_c \lesssim T_{\mathrm{sym}} \Rightarrow \text{快衰落}
}
$$

## 3.2 双选择性信道的统一模型

对第 $m$ 个 OFDM 符号，若采用准静态近似，则其 CIR 可写为

$$
\begin{aligned}
h_m(\tau)
&\triangleq h(t_m,\tau) \\
&=\sum_{\ell=0}^{L-1}\alpha_{\ell,m}\delta(\tau-\tau_\ell).
\end{aligned}
$$

为了把平均功率与时间变化分离，通常进一步写成

$$
\boxed{
h_m(\tau)=\sum_{\ell=0}^{L-1}\sqrt{P_\ell}\,g_{\ell,m}\,\delta(\tau-\tau_\ell)
}
$$

其中：

- $P_\ell$ 为第 $\ell$ 条路径的平均功率，由 PDP（power delay profile）给出；
- $g_{\ell,m}$ 为第 $\ell$ 条路径在第 $m$ 个 OFDM 符号上的归一化复衰落过程；
- 一般约定
  $$
  \mathbb E|g_{\ell,m}|^2=1.
  $$

于是离散域 CFR 可写为

$$
\boxed{
H[m,k]=\sum_{\ell=0}^{L-1}\sqrt{P_\ell}\,g_{\ell,m}
 e^{-j2\pi k\Delta f\tau_\ell}
}
$$

这是一种非常常用的双选择性 OFDM 信道表达。

## 3.3 块衰落与相干时间内慢变

在上述统一模型下，可以得到几种常见简化：

### 1) 块衰落

设一个块内包含 $M_b$ 个 OFDM 符号，若在块内可认为

$$
g_{\ell,m}=g_{\ell,b},\qquad m\in \mathcal B_b,
$$

其中$\mathcal B_b$ 表示 $b$ 个 OFDM 符号的索引集合，其信道CIR在每个符号相同，则得到块衰落模型。通常

$$
M_b\approx \left\lfloor \frac{T_c}{T_{\mathrm{sym}}}\right\rfloor.
$$

### 2) 相干时间内慢变

块衰落只是粗粒度近似。若希望更准确描述相干时间内的连续慢变，则应保留 $g_{\ell,m}$ 随 $m$ 的缓慢变化。这类模型通常通过：

- Jakes / Clarke 随机过程采样；
- 低阶自回归（AR）近似；
- 基展开模型（BEM）；

来实现。

### 3) 快衰落

若信道在一个 OFDM 符号内的变化不可忽略，则仅用 $h_m(\tau)$ 描述每个符号已不充分，需要考虑符号内连续时间变化及其引起的 ICI。
正教时频调制（OTFS）是OFDM的扩展，把时域快变也考虑在建模中，在5G和6G的初期讨论波形方向。不过存在复杂度过高的问题，所以标准化困难。
本文不进一步展开该快衰弱的详细离散化，涉及到ISI的影响。

# 4 路径增益过程的建模

## 4.1 Jakes/Clarke 建模

在各向同性散射假设下，可将第 $\ell$ 条路径的归一化衰落过程表示为一组谐波的叠加：

$$
\boxed{
g_{\ell,m}=\frac{1}{\sqrt{N_c}}
\sum_{n=1}^{N_c}
\exp\left(j\bigl(2\pi \nu_{\max}\cos\theta_n\, mT_{\mathrm{sym}}+\phi_{n,\ell}\bigr)\right)
}
$$

其中：

- $N_c$ 为谐波数，常取 8、16、32 等；
- $\theta_n$ 为离散角度；
- $\phi_{n,\ell}\sim U[0,2\pi)$ 为与路径相关的独立初始相位。

一种常用的角度选取方式为

$$
\theta_n=\frac{2\pi n-\pi+\theta}{4N_c},
\qquad \theta\sim U[0,2\pi).
$$

该模型的时间自相关函数近似为 0 阶第一类贝塞尔函数：

$$
\boxed{
R_g(\Delta t)=J_0(2\pi\nu_{\max}\Delta t)
}
$$

离散到 OFDM 符号时刻，则有

$$
\boxed{
R_g(\Delta m)=J_0(2\pi\nu_{\max}\Delta m T_{\mathrm{sym}})
}
$$

其中 $\Delta m=m_1-m_2$。

## 4.2 自回归建模

为了降低复杂度，常用低阶自回归过程拟合时间相关性。以一阶 AR 模型为例：

$$
\boxed{
g_{\ell,m}=\rho g_{\ell,m-1}+\sqrt{1-|\rho|^2}\,\gamma_{\ell,m}
}
$$

其中：

- $\gamma_{\ell,m}\sim \mathcal{CN}(0,1)$ 为独立创新项；
- $\rho$ 为相邻 OFDM 符号之间的相关系数。

若希望该模型逼近 Jakes 自相关，可取

$$
\boxed{
\rho=J_0(2\pi\nu_{\max}T_{\mathrm{sym}})
}
$$

此时 AR(1) 模型在实现简洁性与时间相关性之间取得了良好折中。

## 4.3 基展开模型（BEM）

另一类常见做法是用一组基函数展开时间变化过程。例如：

$$
g_{\ell,m}=\sum_{q=0}^{Q-1} c_{\ell,q}\,\psi_q[m]
$$

其中：

- $\psi_q[m]$ 为预设基函数；
- $c_{\ell,q}$ 为待定系数。

BEM 适合用于推导时变信道估计、跟踪和预测算法，在 OFDM 双选择性信道处理中尤其常见。

# 5 Rayleigh / Rician 双选择性信道

Rayleigh 和 Rician 分布用于刻画路径增益过程 $g_{\ell,m}$ 的统计性质。

## 5.1 双选择性 Rayleigh 衰落信道

### 5.1.1 定义

对于所有路径 $\ell$，若 $g_{\ell,m}$ 是零均值、单位方差的复高斯随机过程，且不同路径相互独立，则称该信道为双选择性 Rayleigh 衰落信道：

$$
\boxed{
g_{\ell,m}\sim \mathcal{CN}(0,1)
}
$$

此时每条路径的实际增益为 $\sqrt{P_\ell}g_{\ell,m}$。

### 5.1.2 幅度分布

对于任意固定 $m$，幅度

$$
|\sqrt{P_\ell}g_{\ell,m}|
$$

服从 Rayleigh 分布，相位服从均匀分布。

## 5.2 双选择性 Rician 衰落信道

### 5.2.1 物理情形

Rician 信道通常用于描述存在明显 LOS 分量的场景。常见建模方式是：

- 第 0 条路径包含确定性直射分量；
- 其余路径仍为散射分量。

### 5.2.2 建模方式

设第 0 条路径为 LOS 路径，引入 Rician 因子 $K$：

$$
K=\frac{\text{LOS 功率}}{\text{散射功率}}.
$$

则第 0 条路径的归一化增益可以写成

$$
\boxed{
g_{0,m}
=
\sqrt{\frac{K}{K+1}}e^{j(2\pi \nu_{\mathrm{LOS}}mT_{\mathrm{sym}}+\phi_0)}
+
\sqrt{\frac{1}{K+1}}\,\tilde g_{0,m}
}
$$

其中：

- $\tilde g_{0,m}\sim\mathcal{CN}(0,1)$ 为散射分量；
- $\nu_{\mathrm{LOS}}$ 为 LOS 分量的多普勒频移；
- $\phi_0$ 为初始相位。

其余路径仍为 Rayleigh：

$$
g_{\ell,m}=\tilde g_{\ell,m}\sim\mathcal{CN}(0,1),
\qquad \ell=1,2,\dots,L-1.
$$

因此，总 CIR 为

$$
\boxed{
h_m(\tau)=\sqrt{P_0}g_{0,m}\delta(\tau-\tau_0)
+
\sum_{\ell=1}^{L-1}\sqrt{P_\ell}g_{\ell,m}\delta(\tau-\tau_\ell)
}
$$

### 5.2.3 幅度分布

此时第一径的幅度 $|\sqrt{P_0}g_{0,m}|$ 服从 Rician 分布，其余路径幅度仍服从 Rayleigh 分布。

## 5.3 生成步骤

一个常用的双选择性 Rayleigh / Rician 信道生成流程如下：

1. **设定时延与 PDP**
   - 给定路径时延 $\{\tau_\ell\}$；
   - 给定功率延迟谱 $\{P_\ell\}$，满足
     $$
     \sum_{\ell=0}^{L-1}P_\ell=1.
     $$
2. **设定时间相关模型**
   - 根据 Jakes / Clarke、AR 或 BEM 生成 $g_{\ell,m}$。
3. **选择统计模型**
   - Rayleigh：对全部路径使用零均值复高斯过程；
   - Rician：对 LOS 路径加入确定性直射分量。
4. **组合 CIR**
   $$
   h_m(\tau)=\sum_{\ell=0}^{L-1}\sqrt{P_\ell}g_{\ell,m}\delta(\tau-\tau_\ell).
   $$
5. **转换到频域**
   $$
   H[m,k]=\sum_{\ell=0}^{L-1}\sqrt{P_\ell}g_{\ell,m}e^{-j2\pi k\Delta f\tau_\ell}.
   $$

# 6 特殊配置与退化情形

在统一模型

$$
h_m(\tau)=\sum_{\ell=0}^{L-1}\sqrt{P_\ell}g_{\ell,m}\delta(\tau-\tau_\ell)
$$

下，可以得到若干特殊场景：

1. **$L=1$ 且 $g_{0,m}$ 恒定**
   - 平坦衰落 + 慢衰落。
2. **$L\ge 2$ 且 $g_{\ell,m}$ 恒定**
   - 频率选择性衰落 + 慢衰落。
3. **$L=1$ 且 $g_{0,m}$ 变化显著**
   - 平坦衰落 + 时间选择性衰落。
4. **$L\ge 2$ 且 $g_{\ell,m}$ 变化显著**
   - 频率选择性衰落 + 时间选择性衰落，即双选择性衰落。

由此可见，双选择性模型是统一的上层描述，其余常见模型均可视为其特殊情形。

# 7 参考

1. Andrea Goldsmith, *Wireless Communications*.
2. Mingduo Liao, PhD Thesis, 2022, NYU, *Channel Estimation for Massive MIMO Systems*.
3. 3GPP TR 38.901, *Study on channel model for frequencies from 0.5 to 100 GHz*.
