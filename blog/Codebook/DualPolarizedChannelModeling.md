# 双极化收发天线下的 MIMO-OFDM 信道建模

> 本文讨论在 MIMO-OFDM 信道建模中引入双极化收发天线后的系统模型、建模动机、与传统 MIMO-OFDM 信道模型的区别，以及常用的统计相关模型和角度-时延域几何模型。

## 目录

- [1 背景与动机](#1-背景与动机)
  - [1.1 传统 MIMO-OFDM 信道建模](#11-传统-mimo-ofdm-信道建模)
  - [1.2 为什么需要双极化建模](#12-为什么需要双极化建模)
  - [1.3 双极化与普通多天线端口的区别](#13-双极化与普通多天线端口的区别)
- [2 Notation](#2-notation)
- [3 双极化 MIMO-OFDM 系统模型](#3-双极化-mimo-ofdm-系统模型)
  - [3.1 端口、物理天线位置与极化维度](#31-端口物理天线位置与极化维度)
  - [3.2 频域输入输出模型](#32-频域输入输出模型)
  - [3.3 极化子信道分块表示](#33-极化子信道分块表示)
- [4 与传统 MIMO-OFDM 建模的区别](#4-与传统-mimo-ofdm-建模的区别)
- [5 Kronecker 空间-极化相关模型](#5-kronecker-空间-极化相关模型)
  - [5.1 传统 Kronecker 模型](#51-传统-kronecker-模型)
  - [5.2 空间-极化联合相关矩阵](#52-空间-极化联合相关矩阵)
  - [5.3 XPR/XPD 与极化耦合](#53-xprxpd-与极化耦合)
  - [5.4 双极化 Kronecker TDL/CFR 生成](#54-双极化-kronecker-tdlcfr-生成)
- [6 角度-时延域双极化几何模型](#6-角度-时延域双极化几何模型)
  - [6.1 普通角度-时延域 MIMO 模型](#61-普通角度-时延域-mimo-模型)
  - [6.2 路径级 2x2 极化耦合矩阵](#62-路径级-2x2-极化耦合矩阵)
  - [6.3 带极化响应的阵列导向矢量](#63-带极化响应的阵列导向矢量)
  - [6.4 双极化 OFDM CFR](#64-双极化-ofdm-cfr)
- [7 仿真建模流程](#7-仿真建模流程)
- [8 模型选择建议](#8-模型选择建议)
- [9 小结](#9-小结)

# 1 背景与动机

## 1.1 传统 MIMO-OFDM 信道建模

在传统 MIMO-OFDM 系统中，第 $m$ 个 OFDM 符号、第 $k$ 个子载波上的频域输入输出模型通常写为：

$$
\mathbf y[m,k]
=
\mathbf H[m,k]\mathbf x[m,k]+\mathbf w[m,k],
$$

其中：

$$
\mathbf H[m,k]\in\mathbb C^{N_r\times N_t}.
$$

$N_t$ 是发射天线端口数，$N_r$ 是接收天线端口数。矩阵元素 $H_{r,t}[m,k]$ 表示第 $t$ 个发射端口到第 $r$ 个接收端口的子信道。

传统 MIMO-OFDM 信道建模主要关注四类选择性：

| 维度 | 物理来源 | 常用数学对象 |
|---|---|---|
| 时延选择性 | 多径传播 | CIR、PDP、TDL/CDL |
| 频率选择性 | 有限时延扩展 | 子载波 CFR、频域相关矩阵 |
| 时间选择性 | 多普勒扩展 | Jakes、AR、BEM |
| 空间选择性 | 多天线阵列、AoA/AoD | 空间相关矩阵、阵列导向矢量、角度域稀疏性 |

若不考虑极化，每个天线端口通常被视为一个标量端口。信道矩阵只描述不同收发端口之间的复增益、相位、相关性和频率选择性。

## 1.2 为什么需要双极化建模

实际蜂窝基站和大规模 MIMO 阵列中，常用双极化天线单元。例如：

- 垂直/水平极化：$V/H$；
- 正负 45 度交叉极化：$+45^\circ/-45^\circ$；
- 每个物理天线位置上集成两个极化端口。

双极化天线的工程动机包括：

1. **在有限阵列尺寸内增加端口数**

   若每个物理位置只有一个极化端口，要增加端口数通常需要增加阵列尺寸。双极化可以在同一物理位置提供两个端口，从而提升端口密度。

2. **利用极化分集**

   不同极化方向对传播环境的响应不同。即使空间位置相同，两个极化端口的衰落也不完全相同，因此可以提供额外分集。

3. **支持空分复用和多用户传输**

   双极化端口可以作为独立 MIMO 端口参与预编码、信道估计和数据传输。

4. **贴近 4G/5G/6G 阵列形态**

   现代基站阵列通常不是单极化 ULA，而是双极化 UPA 或多面板阵列。若仿真中忽略极化，可能高估或低估真实系统中的信道秩、端口相关性和波束性能。

## 1.3 双极化与普通多天线端口的区别

从基带信号处理角度看，双极化端口可以被并入普通 MIMO 端口数：

$$
N_t^{\mathrm{port}} = N_t^{\mathrm{pos}}N_t^{\mathrm{pol}},
\qquad
N_r^{\mathrm{port}} = N_r^{\mathrm{pos}}N_r^{\mathrm{pol}}.
$$

如果每个物理天线位置有两个极化端口，则：

$$
N_t^{\mathrm{pol}}=2,\qquad N_r^{\mathrm{pol}}=2.
$$

因此，信道矩阵仍然可以写成：

$$
\mathbf H[m,k]\in
\mathbb C^{N_r^{\mathrm{port}}\times N_t^{\mathrm{port}}}.
$$

但是，双极化端口不能简单理解为“多了几根独立天线”。它额外引入：

- 同极化传播：$V\to V$、$H\to H$；
- 交叉极化传播：$V\to H$、$H\to V$；
- 极化隔离和交叉极化泄漏；
- 同一物理位置两个极化端口之间的相关性；
- 路径级极化旋转、反射和散射导致的极化混合。

所以双极化建模的关键不是只增加端口数，而是要描述**空间维度和极化维度的联合统计结构或联合几何结构**。

# 2 Notation

| 符号 | 含义 |
|---|---|
| $m$ | OFDM 符号索引 |
| $k$ | OFDM 子载波索引 |
| $n$ | 离散时延 tap 索引 |
| $\ell$ | 物理路径、簇或射线索引 |
| $N$ | OFDM FFT 点数或子载波数 |
| $N_t^{\mathrm{pos}}$ | 发射端物理天线位置数 |
| $N_r^{\mathrm{pos}}$ | 接收端物理天线位置数 |
| $N_t^{\mathrm{pol}}$ | 发射端每个位置的极化端口数，双极化时通常为 2 |
| $N_r^{\mathrm{pol}}$ | 接收端每个位置的极化端口数，双极化时通常为 2 |
| $N_t^{\mathrm{port}}$ | 发射端总端口数，$N_t^{\mathrm{pos}}N_t^{\mathrm{pol}}$ |
| $N_r^{\mathrm{port}}$ | 接收端总端口数，$N_r^{\mathrm{pos}}N_r^{\mathrm{pol}}$ |
| $\mathcal P_t,\mathcal P_r$ | 发射端和接收端极化集合，如 $\{V,H\}$ 或 $\{+45^\circ,-45^\circ\}$ |
| $\mathbf x[m,k]$ | 发射向量，$\mathbb C^{N_t^{\mathrm{port}}\times 1}$ |
| $\mathbf y[m,k]$ | 接收向量，$\mathbb C^{N_r^{\mathrm{port}}\times 1}$ |
| $\mathbf H[m,k]$ | 双极化 MIMO-OFDM CFR 矩阵 |
| $\mathbf H[m,n]$ | 双极化 MIMO-OFDM CIR tap 矩阵 |
| $\mathbf R_t,\mathbf R_r$ | 发射端、接收端联合空间-极化相关矩阵 |
| $\mathbf R_{t,\mathrm{space}},\mathbf R_{r,\mathrm{space}}$ | 发射端、接收端空间位置相关矩阵 |
| $\mathbf R_{t,\mathrm{pol}},\mathbf R_{r,\mathrm{pol}}$ | 发射端、接收端极化相关矩阵 |
| $\rho_{\mathrm{pol}}$ | 两个极化端口之间的相关系数 |
| $\mathrm{XPR}$ | Cross-Polarization Power Ratio，交叉极化功率比 |
| $\mathrm{XPD}$ | Cross-Polarization Discrimination，交叉极化鉴别度 |
| $\mathbf P_\ell$ | 第 $\ell$ 条路径的 $2\times2$ 极化耦合矩阵 |
| $\mathbf a_t(\Omega_{t,\ell})$ | 发射端空间阵列导向矢量 |
| $\mathbf a_r(\Omega_{r,\ell})$ | 接收端空间阵列导向矢量 |
| $\mathbf b_t^{(p)}(\Omega)$ | 发射端第 $p$ 个极化的方向图或极化响应 |
| $\mathbf b_r^{(p)}(\Omega)$ | 接收端第 $p$ 个极化的方向图或极化响应 |
| $P_\ell$ | 第 $\ell$ 条路径平均功率 |
| $\tau_\ell$ | 第 $\ell$ 条路径时延 |
| $\alpha_\ell$ | 第 $\ell$ 条路径复增益 |
| $\sigma_w^2$ | AWGN 噪声方差 |
| $\otimes$ | Kronecker 积 |
| $\mathrm{vec}(\cdot)$ | 按列向量化 |

# 3 双极化 MIMO-OFDM 系统模型

## 3.1 端口、物理天线位置与极化维度

双极化阵列中，需要区分两个概念：

- **物理天线位置**：阵列中不同空间位置的天线单元；
- **极化端口**：同一位置上不同极化方向的端口。

发射端端口可以用二元索引表示：

$$
(i_t,p_t),
$$

其中：

- $i_t=0,\dots,N_t^{\mathrm{pos}}-1$ 是发射端物理位置索引；
- $p_t\in\mathcal P_t$ 是发射端极化索引。

接收端类似：

$$
(i_r,p_r).
$$

为了写成普通 MIMO 矩阵，可以把二元索引展平为端口索引：

$$
t = i_t N_t^{\mathrm{pol}} + p_t,
\qquad
r = i_r N_r^{\mathrm{pol}} + p_r.
$$

展平后，双极化系统仍然是 MIMO 系统，但端口相关结构不再只是空间相关。

## 3.2 频域输入输出模型

第 $m$ 个 OFDM 符号、第 $k$ 个子载波上，双极化 MIMO-OFDM 的输入输出模型为：

$$
\boxed{
\mathbf y[m,k]
=
\mathbf H[m,k]\mathbf x[m,k]+\mathbf w[m,k]
}
$$

其中：

$$
\mathbf x[m,k]\in\mathbb C^{N_t^{\mathrm{port}}\times1},
\quad
\mathbf y[m,k]\in\mathbb C^{N_r^{\mathrm{port}}\times1},
$$

$$
\mathbf H[m,k]\in
\mathbb C^{N_r^{\mathrm{port}}\times N_t^{\mathrm{port}}}.
$$

离散时延域和频域的关系仍然是 OFDM 信道建模中的 DFT 关系：

$$
\boxed{
\mathbf H[m,k]
=
\sum_{n=0}^{L_h-1}
\mathbf H[m,n]e^{-j\frac{2\pi}{N}kn}
}
$$

区别在于 $\mathbf H[m,n]$ 的行列索引现在同时包含空间位置和极化端口。

## 3.3 极化子信道分块表示

以 $V/H$ 双极化为例，如果发射端和接收端都采用双极化端口，则可把信道矩阵按极化分块：

$$
\mathbf H[m,k]
=
\begin{bmatrix}
\mathbf H_{VV}[m,k] & \mathbf H_{VH}[m,k]\\
\mathbf H_{HV}[m,k] & \mathbf H_{HH}[m,k]
\end{bmatrix}.
$$

这里每个子块的维度为：

$$
\mathbf H_{p_rp_t}[m,k]
\in
\mathbb C^{N_r^{\mathrm{pos}}\times N_t^{\mathrm{pos}}},
\qquad
p_r,p_t\in\{V,H\}.
$$

各子块含义为：

| 子块 | 含义 |
|---|---|
| $\mathbf H_{VV}$ | 发射 V 极化到接收 V 极化的同极化信道 |
| $\mathbf H_{HH}$ | 发射 H 极化到接收 H 极化的同极化信道 |
| $\mathbf H_{VH}$ | 发射 H 极化到接收 V 极化的交叉极化信道 |
| $\mathbf H_{HV}$ | 发射 V 极化到接收 H 极化的交叉极化信道 |

理想极化隔离下，交叉极化子块很小：

$$
\|\mathbf H_{VH}\|_F^2,\ \|\mathbf H_{HV}\|_F^2
\ll
\|\mathbf H_{VV}\|_F^2,\ \|\mathbf H_{HH}\|_F^2.
$$

实际传播中，反射、绕射、散射和天线方向图都会导致极化旋转，因此交叉极化项通常不能直接忽略。

# 4 与传统 MIMO-OFDM 建模的区别

双极化 MIMO-OFDM 建模和传统 MIMO-OFDM 建模的关系可以概括为：

> 双极化模型仍然是 MIMO-OFDM 模型，但端口维度从“空间端口”扩展为“空间位置 $\times$ 极化端口”的联合维度。

主要区别如下。

| 对比项 | 传统 MIMO-OFDM | 双极化 MIMO-OFDM |
|---|---|---|
| 端口含义 | 通常只表示不同天线位置或不同阵元 | 同时包含空间位置和极化方向 |
| 信道矩阵维度 | $N_r\times N_t$ | $N_r^{\mathrm{pos}}N_r^{\mathrm{pol}}\times N_t^{\mathrm{pos}}N_t^{\mathrm{pol}}$ |
| 相关性 | 主要考虑空间相关 | 需要考虑空间相关、极化相关、空间-极化耦合 |
| 子信道结构 | 每个天线对一个标量信道 | 可分为同极化和交叉极化子信道 |
| 几何模型 | 路径增益多为标量 $\alpha_\ell$ | 路径增益扩展为 $2\times2$ 极化耦合矩阵 |
| 关键参数 | PDP、AoA/AoD、Doppler、空间相关 | 额外需要 XPR/XPD、极化方向图、极化相关 |
| 适用性 | 普通多天线链路级仿真 | 更贴近双极化基站阵列和实际 5G/6G 天线面板 |

因此，双极化建模不是简单给普通 MIMO 信道矩阵加一个更大的相关矩阵，而是要明确端口结构和极化传播机制。

# 5 Kronecker 空间-极化相关模型

## 5.1 传统 Kronecker 模型

传统平坦 MIMO Kronecker 相关模型常写为：

$$
\boxed{
\mathbf H
=
\mathbf R_r^{1/2}
\mathbf H_w
\left(\mathbf R_t^{1/2}\right)^T
}
$$

其中：

- $\mathbf H_w$ 是 i.i.d. 复高斯矩阵；
- $\mathbf R_t$ 是发射端空间相关矩阵；
- $\mathbf R_r$ 是接收端空间相关矩阵。

向量化后：

$$
\mathrm{vec}(\mathbf H)
\sim
\mathcal{CN}
\left(
\mathbf 0,\,
\mathbf R_t^T\otimes\mathbf R_r
\right).
$$

在宽带 TDL 模型中，可对每个时延 tap 使用类似结构：

$$
\mathbf H[m,n]
=
\sqrt{P_n}
\mathbf R_r^{1/2}
\mathbf W[m,n]
\left(\mathbf R_t^{1/2}\right)^T.
$$

## 5.2 空间-极化联合相关矩阵

双极化下，$\mathbf R_t$ 和 $\mathbf R_r$ 不应只表示空间位置相关，而应表示**空间-极化联合相关**。

一种常用简化方式是假设空间相关和极化相关可分离：

$$
\boxed{
\mathbf R_t
=
\mathbf R_{t,\mathrm{space}}
\otimes
\mathbf R_{t,\mathrm{pol}}
}
$$

$$
\boxed{
\mathbf R_r
=
\mathbf R_{r,\mathrm{space}}
\otimes
\mathbf R_{r,\mathrm{pol}}
}
$$

其中：

$$
\mathbf R_{t,\mathrm{space}}
\in
\mathbb C^{N_t^{\mathrm{pos}}\times N_t^{\mathrm{pos}}},
\quad
\mathbf R_{t,\mathrm{pol}}
\in
\mathbb C^{N_t^{\mathrm{pol}}\times N_t^{\mathrm{pol}}},
$$

$$
\mathbf R_{r,\mathrm{space}}
\in
\mathbb C^{N_r^{\mathrm{pos}}\times N_r^{\mathrm{pos}}},
\quad
\mathbf R_{r,\mathrm{pol}}
\in
\mathbb C^{N_r^{\mathrm{pol}}\times N_r^{\mathrm{pol}}}.
$$

双极化时，极化相关矩阵常取 $2\times2$ 形式。例如：

$$
\boxed{
\mathbf R_{\mathrm{pol}}
=
\begin{bmatrix}
1 & \rho_{\mathrm{pol}}\\
\rho_{\mathrm{pol}}^* & 1
\end{bmatrix}
}
$$

$\rho_{\mathrm{pol}}$ 描述同一位置上两个极化端口之间的统计相关性。

若只用该相关矩阵，模型能够描述两个极化端口“相关或不相关”，但还不能完整描述交叉极化功率泄漏。因此通常还需要 XPR/XPD 参数。

## 5.3 XPR/XPD 与极化耦合

XPR（Cross-Polarization Power Ratio）定义为同极化功率与交叉极化功率之比：

$$
\boxed{
\mathrm{XPR}
=
\frac{P_{\mathrm{co}}}{P_{\mathrm{cross}}}
}
$$

以线性值表示时，若 $\mathrm{XPR}=10$，表示交叉极化功率约为同极化功率的 $1/10$。以 dB 表示：

$$
\mathrm{XPR}_{\mathrm{dB}}
=
10\log_{10}(\mathrm{XPR}).
$$

对 $V/H$ 双极化信道，一个简单的极化功率缩放矩阵可以写成：

$$
\mathbf S_{\mathrm{pol}}
=
\begin{bmatrix}
1 & \frac{1}{\sqrt{\kappa}}\\
\frac{1}{\sqrt{\kappa}} & 1
\end{bmatrix},
\qquad
\kappa=\mathrm{XPR}.
$$

这里 $1/\sqrt{\kappa}$ 作用在幅度上，因此交叉极化功率约为 $1/\kappa$。

更一般地，对第 $\ell$ 条路径，可以定义路径级极化耦合矩阵：

$$
\boxed{
\mathbf P_\ell
=
\begin{bmatrix}
\alpha_{\ell,VV} & \alpha_{\ell,VH}\\
\alpha_{\ell,HV} & \alpha_{\ell,HH}
\end{bmatrix}
}
$$

其中常见功率关系为：

$$
\mathbb E|\alpha_{\ell,VV}|^2
\approx
\mathbb E|\alpha_{\ell,HH}|^2
=
1,
$$

$$
\mathbb E|\alpha_{\ell,VH}|^2
\approx
\mathbb E|\alpha_{\ell,HV}|^2
=
\frac{1}{\kappa_\ell}.
$$

$\kappa_\ell$ 可以是固定值，也可以随路径、簇、场景随机变化。

## 5.4 双极化 Kronecker TDL/CFR 生成

在双极化 Kronecker TDL 模型中，每个 tap 的 CIR 可以写成：

$$
\boxed{
\mathbf H[m,n]
=
\sqrt{P_n}
\mathbf R_r^{1/2}
\mathbf W[m,n]
\left(\mathbf R_t^{1/2}\right)^T
}
$$

其中：

$$
\mathbf R_t
=
\mathbf R_{t,\mathrm{space}}\otimes\mathbf R_{t,\mathrm{pol}},
\qquad
\mathbf R_r
=
\mathbf R_{r,\mathrm{space}}\otimes\mathbf R_{r,\mathrm{pol}}.
$$

$\mathbf W[m,n]$ 是维度为：

$$
N_r^{\mathrm{port}}\times N_t^{\mathrm{port}}
$$

的 i.i.d. 复高斯矩阵。若需要显式 XPR，可在生成 $\mathbf W[m,n]$ 后按极化子块缩放，或者直接构造带极化功率结构的 $\mathbf W_{\mathrm{pol}}[m,n]$。

然后通过 DFT 得到 CFR：

$$
\boxed{
\mathbf H[m,k]
=
\sum_{n=0}^{L_h-1}
\mathbf H[m,n]e^{-j\frac{2\pi}{N}kn}
}
$$

这种模型的优点是实现简单、参数少、便于分析信道估计和 BER/NMSE。缺点是物理解释较弱，尤其不能准确描述每条路径的 AoD-AoA-极化联合耦合。

# 6 角度-时延域双极化几何模型

## 6.1 普通角度-时延域 MIMO 模型

传统几何 MIMO 信道通常写成路径叠加：

$$
\mathbf H[m,k]
=
\sum_{\ell=1}^{L}
\sqrt{P_\ell}
\alpha_\ell[m]
\mathbf a_r(\Omega_{r,\ell})
\mathbf a_t^H(\Omega_{t,\ell})
e^{-j\frac{2\pi}{N}k\tau_\ell}.
$$

其中每条路径贡献一个秩一矩阵：

$$
\mathbf a_r(\Omega_{r,\ell})
\mathbf a_t^H(\Omega_{t,\ell}).
$$

该模型显式描述：

- AoA：$\Omega_{r,\ell}$；
- AoD：$\Omega_{t,\ell}$；
- 时延：$\tau_\ell$；
- 多普勒或时间变化：$\alpha_\ell[m]$；
- 路径功率：$P_\ell$。

但普通模型中路径增益 $\alpha_\ell[m]$ 是标量，无法描述极化旋转和交叉极化泄漏。

## 6.2 路径级 2x2 极化耦合矩阵

双极化几何模型中，第 $\ell$ 条路径的极化传播由 $2\times2$ 矩阵描述：

$$
\boxed{
\mathbf P_\ell[m]
=
\begin{bmatrix}
\alpha_{\ell,VV}[m] & \alpha_{\ell,VH}[m]\\
\alpha_{\ell,HV}[m] & \alpha_{\ell,HH}[m]
\end{bmatrix}
}
$$

矩阵元素表示：

| 元素 | 含义 |
|---|---|
| $\alpha_{\ell,VV}$ | 发射 V 到接收 V 的同极化路径增益 |
| $\alpha_{\ell,HH}$ | 发射 H 到接收 H 的同极化路径增益 |
| $\alpha_{\ell,VH}$ | 发射 H 到接收 V 的交叉极化路径增益 |
| $\alpha_{\ell,HV}$ | 发射 V 到接收 H 的交叉极化路径增益 |

一个常用随机参数化为：

$$
\mathbf P_\ell[m]
=
\begin{bmatrix}
e^{j\phi_{\ell,VV}[m]} &
\sqrt{\kappa_\ell^{-1}}e^{j\phi_{\ell,VH}[m]}\\
\sqrt{\kappa_\ell^{-1}}e^{j\phi_{\ell,HV}[m]} &
e^{j\phi_{\ell,HH}[m]}
\end{bmatrix},
$$

其中 $\kappa_\ell=\mathrm{XPR}_\ell$，$\phi_{\ell,\cdot}$ 是随机相位或随时间演化的相位。

如果需要包含路径复幅度，可写成：

$$
\mathbf P_\ell[m]
=
\beta_\ell[m]
\begin{bmatrix}
e^{j\phi_{\ell,VV}[m]} &
\sqrt{\kappa_\ell^{-1}}e^{j\phi_{\ell,VH}[m]}\\
\sqrt{\kappa_\ell^{-1}}e^{j\phi_{\ell,HV}[m]} &
e^{j\phi_{\ell,HH}[m]}
\end{bmatrix}.
$$

其中 $\beta_\ell[m]$ 描述该路径整体的时间选择性衰落。

## 6.3 带极化响应的阵列导向矢量

如果每个位置有两个极化端口，阵列响应不只是空间相位，还包含极化方向图响应。

对第 $\ell$ 条路径，发射端可定义双极化阵列响应矩阵：

$$
\mathbf A_t(\Omega_{t,\ell})
=
\begin{bmatrix}
\mathbf a_{t,V}(\Omega_{t,\ell}) &
\mathbf a_{t,H}(\Omega_{t,\ell})
\end{bmatrix}
\in
\mathbb C^{N_t^{\mathrm{port}}\times 2}.
$$

接收端类似：

$$
\mathbf A_r(\Omega_{r,\ell})
=
\begin{bmatrix}
\mathbf a_{r,V}(\Omega_{r,\ell}) &
\mathbf a_{r,H}(\Omega_{r,\ell})
\end{bmatrix}
\in
\mathbb C^{N_r^{\mathrm{port}}\times 2}.
$$

其中 $\mathbf a_{t,V}$ 和 $\mathbf a_{t,H}$ 不仅包含阵列空间相位，还包含对应极化端口的方向图、增益和相位响应。

若忽略方向图差异，只考虑每个位置的空间相位，则可以近似写成：

$$
\mathbf A_t(\Omega)
\approx
\mathbf a_{t,\mathrm{space}}(\Omega)\otimes\mathbf I_2,
$$

$$
\mathbf A_r(\Omega)
\approx
\mathbf a_{r,\mathrm{space}}(\Omega)\otimes\mathbf I_2.
$$

该近似适合算法验证，但不适合分析真实天线方向图、俯仰角、方位角和极化方向图耦合。

## 6.4 双极化 OFDM CFR

引入极化耦合矩阵后，第 $\ell$ 条路径对 CFR 的贡献为：

$$
\boxed{
\mathbf H_\ell[m,k]
=
\sqrt{P_\ell}
\mathbf A_r(\Omega_{r,\ell})
\mathbf P_\ell[m]
\mathbf A_t^H(\Omega_{t,\ell})
e^{-j\frac{2\pi}{N}k\tau_\ell}
}
$$

因此完整双极化几何 MIMO-OFDM 信道为：

$$
\boxed{
\mathbf H[m,k]
=
\sum_{\ell=1}^{L}
\sqrt{P_\ell}
\mathbf A_r(\Omega_{r,\ell})
\mathbf P_\ell[m]
\mathbf A_t^H(\Omega_{t,\ell})
e^{-j\frac{2\pi}{N}k\tau_\ell}
}
$$

该式和普通角度-时延域 MIMO 模型的区别在于：

- 普通模型中路径增益是标量 $\alpha_\ell[m]$；
- 双极化模型中路径增益扩展为 $2\times2$ 极化耦合矩阵 $\mathbf P_\ell[m]$；
- 普通阵列导向矢量扩展为带极化响应的阵列响应矩阵 $\mathbf A_t,\mathbf A_r$。

如果进一步把 CFR 沿角度和时延离散化，可以得到角度-时延-极化域的稀疏表示。此时每条路径不仅对应 AoA/AoD/时延支撑，还对应一个极化耦合矩阵。

# 7 仿真建模流程

双极化 MIMO-OFDM 信道仿真可按以下流程实现。

**Step 1：设置系统和阵列参数**

- $N_t^{\mathrm{pos}},N_r^{\mathrm{pos}}$；
- 极化数 $N_t^{\mathrm{pol}}=N_r^{\mathrm{pol}}=2$；
- OFDM 参数 $N,\Delta f,T_{\mathrm{cp}}$；
- 载频、速度、多普勒；
- 天线间距和阵列类型。

**Step 2：设置路径或 tap 参数**

- 路径数 $L$ 或 tap 数 $L_h$；
- PDP $P_\ell$；
- 时延 $\tau_\ell$；
- AoA/AoD；
- 多普勒频移。

**Step 3：设置极化参数**

- 极化基：$V/H$ 或 $+45^\circ/-45^\circ$；
- XPR/XPD；
- 极化相关系数 $\rho_{\mathrm{pol}}$；
- 是否使用极化方向图。

**Step 4A：若采用 Kronecker 模型**

1. 构造 $\mathbf R_{t,\mathrm{space}}$ 和 $\mathbf R_{r,\mathrm{space}}$；
2. 构造 $\mathbf R_{t,\mathrm{pol}}$ 和 $\mathbf R_{r,\mathrm{pol}}$；
3. 得到联合相关矩阵：

   $$
   \mathbf R_t
   =
   \mathbf R_{t,\mathrm{space}}\otimes\mathbf R_{t,\mathrm{pol}},
   \quad
   \mathbf R_r
   =
   \mathbf R_{r,\mathrm{space}}\otimes\mathbf R_{r,\mathrm{pol}}.
   $$

4. 生成每个 tap 的 i.i.d. 矩阵并施加相关；
5. 根据 XPR 缩放交叉极化子块；
6. 沿时延维 DFT 得到 CFR。

**Step 4B：若采用角度-时延域模型**

1. 为每条路径生成 AoA/AoD、时延和功率；
2. 为每条路径生成 XPR 和极化耦合矩阵 $\mathbf P_\ell$；
3. 构造带极化响应的阵列矩阵 $\mathbf A_t,\mathbf A_r$；
4. 按路径叠加公式生成 $\mathbf H[m,k]$；
5. 若需要时变信道，对 $\mathbf P_\ell[m]$ 或 $\beta_\ell[m]$ 引入多普勒相位。

# 8 模型选择建议

| 研究目标 | 推荐模型 | 原因 |
|---|---|---|
| 链路级 BER/NMSE 快速仿真 | 双极化 Kronecker TDL | 参数少、实现简单、容易控制相关强度 |
| 研究空间相关和极化相关影响 | 空间-极化 Kronecker | 可直接调节 $\rho_{\mathrm{space}}$ 和 $\rho_{\mathrm{pol}}$ |
| 研究 XPR 对秩和容量的影响 | Kronecker + 极化子块缩放 | 便于单独扫描 XPR |
| 研究波束赋形和角度分辨率 | 双极化角度-时延域模型 | 显式包含 AoA/AoD 和阵列响应 |
| 研究真实阵列方向图 | 几何模型 + 极化方向图 | 可描述方向图、俯仰角、方位角和极化响应 |
| 对齐 3GPP 标准仿真 | CDL/TDL 双极化模型 | 3GPP CDL 显式包含角度、XPR 和极化 |

一般建议：

- 如果目标是验证信道估计、插值、均衡算法，先使用双极化 Kronecker TDL；
- 如果目标是研究波束、CSI 反馈、角度域稀疏估计，应使用角度-时延域模型；
- 如果目标是贴近 5G NR 系统级或链路级评估，应参考 3GPP TR 38.901 CDL 模型；
- 如果只把两个极化端口当成独立天线，而不引入 XPR/极化相关，仿真结果通常会过于理想。

# 9 小结

双极化收发天线建模可以看成 MIMO-OFDM 信道建模的扩展，但不是简单增加天线数。它要求把端口结构从单一空间维度扩展为空间-极化联合维度。

核心区别是：

- 传统 MIMO-OFDM 主要描述天线位置之间的空间相关和角度选择性；
- 双极化 MIMO-OFDM 还需要描述同极化、交叉极化、XPR/XPD、极化相关和路径级极化耦合。

在建模方法上，传统 MIMO-OFDM 中常用的两类方式仍然适用：

1. **Kronecker 相关模型**

   将空间相关矩阵扩展为空间-极化联合相关矩阵：

   $$
   \mathbf R
   =
   \mathbf R_{\mathrm{space}}
   \otimes
   \mathbf R_{\mathrm{pol}}.
   $$

2. **角度-时延域几何模型**

   将每条路径的标量复增益扩展为 $2\times2$ 极化耦合矩阵：

   $$
   \mathbf P_\ell
   =
   \begin{bmatrix}
   \alpha_{\ell,VV} & \alpha_{\ell,VH}\\
   \alpha_{\ell,HV} & \alpha_{\ell,HH}
   \end{bmatrix}.
   $$

如果要在仿真中兼顾复杂度和物理解释，可以先用 Kronecker 空间-极化模型建立基线，再用角度-时延域双极化模型验证波束和几何路径相关结论。
