# Block Diagonalization Precoding

## 目录

- [1 概述](#1-概述)
- [2 系统模型](#2-系统模型)
- [3 Block Diagonalization 的基本思想](#3-block-diagonalization-的基本思想)
- [4 发送端预编码矩阵设计](#4-发送端预编码矩阵设计)
  - [4.1 构造其他用户的聚合信道](#41-构造其他用户的聚合信道)
  - [4.2 第一阶段：消除多用户干扰](#42-第一阶段消除多用户干扰)
  - [4.3 第二阶段：对等效单用户信道并行化](#43-第二阶段对等效单用户信道并行化)
- [5 接收端线性处理](#5-接收端线性处理)
- [6 功率分配](#6-功率分配)
- [7 Regularized BD](#7-regularized-bd)
- [8 参考](#8-参考)

# 1 概述

Block Diagonalization（BD）是一类经典的多用户 MIMO 下行线性预编码方法。其核心思想是：

> **为每个用户设计发送端预编码矩阵，使其信号落在其他用户聚合信道的零空间中，从而在发送端消除多用户干扰。**

与单用户 MIMO 中的 SVD 预编码相比，BD 的关键难点在于：

- 不仅要提升目标用户的有效信道质量；
- 还必须保证对其他用户不产生干扰。

因此，BD 通常分两步完成：

1. 在发送端构造零空间，消除多用户干扰；
2. 在该零空间内进一步对目标用户的等效信道做 SVD，以获得并行空间子信道。

# 2 系统模型

考虑下行多用户 MIMO 广播信道。设：

- 基站发送天线数：$N_t$
- 用户数：$K$
- 每个用户接收天线数：$N_r$
- 每个用户发送流数：$S$

记第 $k$ 个用户的信道矩阵为

$$
\mathbf H_k\in\mathbb C^{N_r\times N_t},\qquad k=1,2,\dots,K.
$$

第 $k$ 个用户的发送符号向量为

$$
\mathbf d_k\in\mathbb C^{S\times 1}.
$$

第 $k$ 个用户的发送端预编码矩阵记为

$$
\mathbf F_k\in\mathbb C^{N_t\times S}.
$$

则总发送信号为

$$
\mathbf x=\sum_{k=1}^{K}\mathbf F_k\mathbf d_k.
$$

第 $k$ 个用户接收到的信号为

$$
\mathbf y_k=\mathbf H_k\mathbf x+\mathbf n_k
=\mathbf H_k\mathbf F_k\mathbf d_k+\sum_{j\ne k}\mathbf H_k\mathbf F_j\mathbf d_j+\mathbf n_k,
$$

其中 $\mathbf n_k\in\mathbb C^{N_r\times 1}$ 为接收噪声。

若进一步引入接收端线性处理矩阵

$$
\mathbf W_k\in\mathbb C^{S\times N_r},
$$

则判决变量为

$$
\hat{\mathbf d}_k=\mathbf W_k\mathbf y_k.
$$

为了使 BD 至少可行，需要零空间维度足够大。一个常见充分条件是

$$
\boxed{N_t\ge K N_r}
$$

更一般地，若每个用户发送 $S$ 条流，则需要目标用户在消除其他用户后仍保留至少 $S$ 维有效空间。

# 3 Block Diagonalization 的基本思想

对于每个用户 $k$，希望预编码矩阵 $\mathbf F_k$ 满足：

$$
\boxed{\mathbf H_j\mathbf F_k=\mathbf 0,\qquad j\ne k}
$$

也就是说，用户 $k$ 的发送信号不对其他任何用户造成干扰。

若该条件成立，则第 $k$ 个用户接收信号化简为

$$
\mathbf y_k=\mathbf H_k\mathbf F_k\mathbf d_k+\mathbf n_k.
$$

这样，多用户系统就被分解成若干个相互无干扰的单用户 MIMO 子问题。

# 4 发送端预编码矩阵设计

## 4.1 构造其他用户的聚合信道

对第 $k$ 个用户，定义除其自身外其他用户的聚合信道矩阵：

$$
\tilde{\mathbf H}_k=
\begin{bmatrix}
\mathbf H_1^T & \cdots & \mathbf H_{k-1}^T & \mathbf H_{k+1}^T & \cdots & \mathbf H_K^T
\end{bmatrix}^T
\in\mathbb C^{(K-1)N_r\times N_t}.
$$

那么要满足

$$
\mathbf H_j\mathbf F_k=\mathbf 0,\quad j\ne k,
$$

等价于要求

$$
\boxed{\tilde{\mathbf H}_k\mathbf F_k=\mathbf 0.}
$$

因此，$\mathbf F_k$ 的列空间必须位于 $\tilde{\mathbf H}_k$ 的右零空间中。

## 4.2 第一阶段：消除多用户干扰

对 $\tilde{\mathbf H}_k$ 做 SVD：

$$
\tilde{\mathbf H}_k
=\mathbf U_k\boldsymbol\Sigma_k
\begin{bmatrix}
\mathbf V_k^{(1)} & \mathbf V_k^{(0)}
\end{bmatrix}^H.
$$

其中：

- $\mathbf U_k\in\mathbb C^{(K-1)N_r\times (K-1)N_r}$；
- $\boldsymbol\Sigma_k\in\mathbb C^{(K-1)N_r\times N_t}$；
- $\mathbf V_k^{(1)}$ 对应非零奇异值；
- $\mathbf V_k^{(0)}$ 对应零奇异值。

由于 $\mathbf V_k^{(0)}$ 张成 $\tilde{\mathbf H}_k$ 的右零空间，因此有

$$
\tilde{\mathbf H}_k\mathbf V_k^{(0)}=\mathbf 0.
$$

于是，只要令

$$
\mathbf F_k=\mathbf V_k^{(0)}\mathbf T_k,
$$

其中 $\mathbf T_k$ 为任意维度匹配矩阵，就自动满足对其他用户零干扰。

## 4.3 第二阶段：对等效单用户信道并行化

为了进一步优化用户 $k$ 自身的有效信道，引入等效信道

$$
\mathbf H_{\mathrm{eff},k}=\mathbf H_k\mathbf V_k^{(0)}.
$$

注意其维度为

$$
\mathbf H_{\mathrm{eff},k}\in\mathbb C^{N_r\times d_k},
$$

其中 $d_k=\dim\mathcal N(\tilde{\mathbf H}_k)$ 为零空间维数。

再对等效信道做 SVD：

$$
\mathbf H_{\mathrm{eff},k}
=\bar{\mathbf U}_k\bar{\boldsymbol\Sigma}_k\bar{\mathbf V}_k^H.
$$

若发送 $S$ 条流，则取 $\bar{\mathbf V}_k$ 的前 $S$ 列，记为

$$
\bar{\mathbf V}_{k,s}\in\mathbb C^{d_k\times S}.
$$

则 BD 预编码矩阵可写为

$$
\boxed{
\mathbf F_k=\mathbf V_k^{(0)}\bar{\mathbf V}_{k,s}
}
$$

其维度为

$$
(N_t\times d_k)(d_k\times S)=N_t\times S.
$$

这是标准 BD 预编码结构。

将其代回可得：

$$
\mathbf H_k\mathbf F_k
=\mathbf H_k\mathbf V_k^{(0)}\bar{\mathbf V}_{k,s}
=\bar{\mathbf U}_k\bar{\boldsymbol\Sigma}_k
\begin{bmatrix}
\mathbf I_S \\
\mathbf 0
\end{bmatrix},
$$

即用户 $k$ 的有效信道被进一步并行化为若干独立空间子信道。

# 5 接收端线性处理

由上式可见，对用户 $k$ 的等效单用户信道，最自然的接收端线性处理矩阵为左奇异矩阵的前 $S$ 行，即

$$
\boxed{
\mathbf W_k=\bar{\mathbf U}_{k,s}^H
}
$$

其中 $\bar{\mathbf U}_{k,s}\in\mathbb C^{N_r\times S}$ 为 $\bar{\mathbf U}_k$ 的前 $S$ 列。

于是判决变量为

$$
\hat{\mathbf d}_k
=\mathbf W_k\mathbf y_k
=\bar{\mathbf U}_{k,s}^H\mathbf H_k\mathbf F_k\mathbf d_k
+\bar{\mathbf U}_{k,s}^H\mathbf n_k.
$$

由于

$$
\bar{\mathbf U}_{k,s}^H\mathbf H_k\mathbf F_k
=\bar{\boldsymbol\Sigma}_{k,s},
$$

其中 $\bar{\boldsymbol\Sigma}_{k,s}$ 为前 $S$ 个非零奇异值形成的对角矩阵，因此最终得到一组并行子信道。

# 6 功率分配

在得到并行子信道之后，可以对每个用户做功率分配。记第 $k$ 个用户分配矩阵为

$$
\boldsymbol\Lambda_k=\operatorname{diag}(p_{k,1},\dots,p_{k,S}),
$$

则实际发送预编码矩阵可写为

$$
\boxed{
\mathbf F_k=\mathbf V_k^{(0)}\bar{\mathbf V}_{k,s}\boldsymbol\Lambda_k^{1/2}
}
$$

对应的功率约束为

$$
\operatorname{tr}(\mathbf F_k\mathbf F_k^H)\le P_k.
$$

由于经 SVD 并行化后，每个用户得到 $S$ 个近似独立子信道，因此在每个用户内部可采用注水（water-filling）或等功率分配。若系统总和速率最大化问题在 BD 约束下分解到各用户并行子信道上，则注水分配是常见做法。

需要说明的是：

- 原始的多用户和速率最大化问题在一般形式下并不简单地“直接是凸问题”；
- BD 的价值在于先通过零空间构造消除多用户干扰，再把问题转化为若干单用户并行子信道上的功率分配问题。

# 7 Regularized BD

经典 BD 完全消除多用户干扰，但在低信噪比下可能过度牺牲阵列增益。为兼顾干扰抑制与噪声增强，常引入 **Regularized BD (RBD)**。

其思想是在第一阶段的零空间设计中，不再要求严格满足

$$
\tilde{\mathbf H}_k\mathbf F_k=\mathbf 0,
$$

而是通过正则化最小化

$$
\left\|\tilde{\mathbf H}_k\mathbf F_k\right\|_F^2
+
\gamma_k\left\|\mathbf F_k\right\|_F^2,
$$

其中 $\gamma_k$ 为正则化系数。这样得到的预编码方向不再严格落在零空间中，但通常在中低 SNR 下可获得更好的性能折中。

因此可以将 RBD 视为：

- 高 SNR 时逼近 BD；
- 低 SNR 时适当保留主信道增益、减少噪声放大。

# 8 参考

1. *Low-Complexity Design of Generalized Block Diagonalization Precoding Algorithms for Multiuser MIMO Systems*, 2014.
2. Q. H. Spencer, A. L. Swindlehurst, and M. Haardt, “Zero-forcing methods for downlink spatial multiplexing in multiuser MIMO channels.”
3. D. Tse and P. Viswanath, *Fundamentals of Wireless Communication*.
