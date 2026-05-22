# SISO OFDM channel estimation

## 目录

- [1 1概述](#1-1概述)

- [2 2 系统模型](#2-2-系统模型)
  - [2.1 2.1 信道模型](#21-21-信道模型)
    - [2.1.1 离散域表示](#211-离散域表示)

  - [2.2 2.2 接收信号模型](#22-22-接收信号模型)

- [3 3 信道估计与插值](#3-3-信道估计与插值)
  - [3.1 3.1 导频位置上的信道估计](#31-31-导频位置上的信道估计)
    - [3.1.1 3.1.1 Least Squares (LS)](#311-311-Least-Squares-LS)

    - [3.1.2 3.1.2 Minimum Mean Square Error (MMSE)](#312-312-Minimum-Mean-Square-Error-MMSE)

    - [3.1.3 3.1.3 Linear MMSE (LMMSE)](#313-313-Linear-MMSE-LMMSE)

  - [3.2 3.2 插值](#32-32-插值)
    - [3.2.1 3.2.1 线性插值](#321-321-线性插值)
      - [3.2.1.1 一维线性插值](#3211-一维线性插值)

      - [3.2.1.2 二维线性插值](#3212-二维线性插值)

    - [3.2.2 3.2.2 二阶插值](#322-322-二阶插值)
      - [3.2.2.1 一维二阶插值](#3221-一维二阶插值)

      - [3.2.2.2 二维二阶插值](#3222-二维二阶插值)

    - [3.2.3 3.2.3 Wiener (LMMSE) 插值](#323-323-Wiener-LMMSE-插值)
      - [3.2.3.1 一维插值](#3231-一维插值)

      - [3.2.3.2 导频估计噪声的影响](#3232-导频估计噪声的影响)

      - [3.2.3.3 二维插值](#3233-二维插值)

    - [3.2.4 3.2.4 DFT 插值](#324-324-DFT-插值)
      - [3.2.4.1 3.2.4.1 一维频域 DFT 插值](#3241-3241-一维频域-DFT-插值)
        - [3.2.4.1.1 Step 1：由导频 CFR 恢复低维 CIR](#32411-Step-1由导频-CFR-恢复低维-CIR)

        - [3.2.4.1.2 Step 2：时域截断 / 去噪](#32412-Step-2时域截断--去噪)

        - [3.2.4.1.3 Step 3：零填充并恢复全部子载波的 CFR](#32413-Step-3零填充并恢复全部子载波的-CFR)

      - [3.2.4.2 3.2.4.2 二维 DFT 插值](#3242-3242-二维-DFT-插值)

    - [3.2.5 3.2.5 DFT 插值的成立条件](#325-325-DFT-插值的成立条件)

- [4 4 References](#4-4-References)


# 1 1概述

信道估计与信道模型、导频图样和接收机复杂度密切相关。通常需要根据：

- 信道的时间选择性与频率选择性；
- 导频的时频分布；
- 期望复杂度与性能折中；

来设计信道估计与插值算法。

在 OFDM 系统中，信道处理通常分为两步：

1. **导频位置上的信道估计**：在导频所在的时频位置估计信道；
2. **插值/外推**：利用导频估计结果恢复数据子载波处的信道。

常见导频图样主要有两类：

- **块状导频（block-type pilots）**：某些 OFDM 符号上的全部子载波均承载导频，适用于**时间选择性较弱、频率选择性较强**的场景；
- **梳状导频（comb-type pilots）**：每个或部分 OFDM 符号上仅部分子载波承载导频，适用于**时间选择性较强、频率选择性较弱**的场景。

对于多符号场景，通常还需要利用时间相关性进行二维插值，甚至进一步做信道跟踪或预测。本文仅讨论**导频辅助的信道估计与插值**。

# 2 2 系统模型

## 2.1 2.1 信道模型

考虑单输入单输出（SISO）OFDM 系统。记：

- 路径数：$L$
- 最大时延：$\tau_{\max}$
- 最大移动速度：$v_{\max}$
- 最大多普勒频移：$\nu_{\max} = \dfrac{v_{\max}}{\lambda},\ \lambda = \dfrac{c}{f_c}$
- 载波频率：$f_c$
- 子载波间隔：$\Delta f$
- 每个 OFDM 符号的子载波数：$N$
- 有效 OFDM 符号时长：$T_u = 1/\Delta f$
- 循环前缀长度：$T_{cp}$
- OFDM 总符号时长：$T_{\mathrm{sym}} = T_u + T_{cp}$
- OFDM 符号数：$M$，其中 $m=0,1,\dots,M-1$

为避免 ISI，通常要求：$T _{cp} > \tau_{\max}$.

记第 $m$ 个 OFDM 符号对应的时刻为：$t_m = mT_{\mathrm{sym}}.$

若在第 $m$ 个 OFDM 符号内采用准静态近似，则该符号经历的连续域 CIR 可表示为：

$$
\begin{aligned}
h_m(\tau)
&\triangleq h(t_m,\tau) 
&= \sum _{\ell=0}^{L-1} \beta_\ell e^{j2\pi \nu _\ell t_m}\delta(\tau-\tau _\ell) 
&= \sum_{\ell=0}^{L-1} \alpha _{\ell,m}\delta(\tau-\tau_\ell),
\end{aligned}
$$

其中：

- $\tau_\ell$ 为第 $\ell$ 条路径时延；
- $\beta_\ell$ 为第 $\ell$ 条路径的初始复增益；
- $\nu_\ell$ 为第 $\ell$ 条路径的多普勒频移；
- $\alpha_{\ell,m} \triangleq \beta_\ell e^{j2\pi \nu_\ell t_m}$ 为第 $m$ 个 OFDM 符号上第 $\ell$ 条路径的等效复增益。

### 2.1.1 离散域表示

采样间隔为：$T _{\mathrm{samp}} = \frac{1}{N\Delta f}$.

对应的离散 CIR 记为：$h[m,n],\qquad n=0,1,\dots,N-1$.

在整数时延 tap 模型下，$h[m,n]$ 仅在有限个 tap 上非零。

相应的 CFR 定义为对离散 CIR 做 **DFT**：$H[m,k]=\sum _{n=0}^{N-1} h[m,n] e^{-j\frac{2\pi}{N}kn},\qquad k=0,1,\dots,N-1.$

其逆变换为 **IDFT**：$h[m,n]=\frac{1}{N}\sum _{k=0}^{N-1} H[m,k] e^{j\frac{2\pi}{N}kn},\qquad n=0,1,\dots,N-1.$

## 2.2 2.2 接收信号模型

设第 $m$ 个 OFDM 符号在频域上的发送符号为：$X_m[k],\qquad k=0,1,\dots,N-1,$ 其中 $X_m[k]\in \mathcal Q$，$\mathcal Q$ 为调制星座集合，如 QPSK、16-QAM、64-QAM 等。

经 IDFT 后得到时域发送信号 $x_m[n]$。假设符号功率已归一化，即：$\mathbb E{|X_m[k]|^2}=1.$

若 CP 足够长，去 CP 后可将线性卷积等效为循环卷积，则时域接收信号为: $y_m[n]=h_m[n]\circledast x_m[n]+w_m[n],\qquad n=0,1,\dots,N-1,$ 其中：

- $\circledast$ 表示循环卷积；
- $w_m[n]$ 为复高斯白噪声，噪声功率为 $\sigma^2$。

矩阵形式可写为：
$\mathbf y_m=\mathrm{circ}(\mathbf h_m)\mathbf x_m+\mathbf w_m$

或等价写为
$\mathbf y_m=\mathrm{circ}(\mathbf x_m)\mathbf h_m+\mathbf w_m,$

二者只是循环卷积矩阵由哪个向量生成的记法不同。

对接收信号做 DFT，得到频域输入输出关系：


$$
Y_m[k]=H_m[k]X_m[k]+W_m[k],\qquad k=0,1,\dots,N-1.
$$

矩阵形式为：

$$
\mathbf Y_m=\operatorname{diag}(\mathbf X_m)\mathbf H_m+\mathbf W_m.
$$

将频域位置划分为导频集合 $\mathcal P$ 与数据集合 $\mathcal D$，则有：

$$
\begin{aligned} \mathbf Y _{m,p}&=\operatorname{diag}(\mathbf X_{m,p})\mathbf H _{m,p}+\mathbf W_{m,p},\\
\mathbf Y _{m,d}&=\operatorname{diag}(\mathbf X_{m,d})\mathbf H _{m,d}+\mathbf W_{m,d}. \end{aligned} 
$$

若把多个 OFDM 符号上的导频统一向量化，则可进一步写成：

$$
\begin{aligned}
\mathbf Y_p&=\operatorname{diag}(\mathbf X_p)\mathbf H_p+\mathbf W_p, \\
\mathbf Y_d&=\operatorname{diag}(\mathbf X_d)\mathbf H_d+\mathbf W_d,
\end{aligned} 
$$

其中：

- $\mathbf X_p,\mathbf H_p,\mathbf Y_p,\mathbf W_p\in\mathbb C^{N_p\times 1}$；
- $\mathbf X_d,\mathbf H_d,\mathbf Y_d,\mathbf W_d\in\mathbb C^{N_d\times 1}$。

这里 $N_p$ 为所考虑时频区域内的导频总数，$N_d$ 为数据位置总数。

# 3 3 信道估计与插值

信道恢复通常分为两部分：

1. **导频位置上的信道估计**
2. **从导频位置到数据位置的插值**

本文介绍三类常见导频估计算法：

- LS
- MMSE
- LMMSE

以及四类常见插值算法：

- 线性插值
- 二阶插值
- Wiener (LMMSE) 插值
- DFT 插值

## 3.1 3.1 导频位置上的信道估计

### 3.1.1 3.1.1 Least Squares (LS)

在导频位置上，LS 估计通过最小化残差平方和得到：
$\hat{\mathbf H}_p =  \arg\min _{\mathbf H_p\in\mathbb C^{N_p\times 1}}
\left|
\mathbf Y_p-\operatorname{diag}(\mathbf X_p)\mathbf H_p
\right|_2^2.$

其闭式解为：
$\hat{\mathbf H}_p = \left(
\operatorname{diag}(\mathbf X_p)^H
\operatorname{diag}(\mathbf X_p)
\right)^{-1}
\operatorname{diag}(\mathbf X_p)^H\mathbf Y_p$.

由于 $\operatorname{diag}(\mathbf X_p)$ 为对角矩阵，只要导频符号均非零，则可进一步化简为：


$$
\boxed{
\hat{\mathbf H}_p=\operatorname{diag}(\mathbf X_p)^{-1}\mathbf Y_p
}
$$

即逐导频点相除。

若导频功率统一满足 $|X_p[i]|^2=\sigma_x^2$，则 LS 估计误差均方值与 SNR 成反比。LS 的优点是实现简单；缺点是未利用信道统计信息，抗噪性能较弱。

估计误差:

$$
\begin{aligned}
\epsilon &= \mathbb{E}\left\{ \lVert \hat{\mathbf{H}}_p - \mathbf{H}_p \rVert^2_F \right\} \\
& = \mathbb{E}\left\{ \lVert \left( \mathrm{diag}(\mathbf{X}_p) \mathrm{diag}(\mathbf{X}_p)^H\right)^{-1} \mathrm{diag}(\mathbf{X}_p)^H \mathbf{Y}_p - \mathbf{H}_p \rVert^2_F \right\} \\
& = \mathbb{E}\left\{ \lVert \left( \mathrm{diag}(\mathbf{X}_p) \mathrm{diag}(\mathbf{X}_p)^H\right)^{-1} \mathrm{diag}(\mathbf{X}_p)^H (\mathrm{diag}(\mathbf{X}_p) \mathbf{H}_p + \mathbf{W}_p) - \mathbf{H}_p \rVert^2_F \right\} \\
& = \frac{\sigma^2}{\sigma^2_x} = \frac{1}{SNR}.

\end{aligned} 
$$

***

### 3.1.2 3.1.2 Minimum Mean Square Error (MMSE)

MMSE 估计试图最小化导频位置上的均方误差。设估计器为线性形式：


$$
\hat{\mathbf H}_p=\mathbf Q_p \mathbf Y_p,
$$

其中 $\mathbf Q_p\in\mathbb C^{N_p\times N_p}$，优化目标为：

$$
\mathbf Q_p^\star =  \arg\min {\mathbf Q_p}
\mathbb E\left\{
\left|
\mathbf H_p-\mathbf Q_p\mathbf Y_p
\right| 2^2
\right\}.
$$

由正交原理可得最优解：


$$
\mathbf Q_p^\star = \mathbb E[\mathbf H_p\mathbf Y_p^H]
\left(
\mathbb E[\mathbf Y_p\mathbf Y_p^H]
\right)^{-1}.
$$

根据模型

$$
\mathbf Y_p=\operatorname{diag}(\mathbf X_p)\mathbf H_p+\mathbf W_p,
$$

有：

$$
\begin{aligned}
\mathbb E[\mathbf H_p\mathbf Y_p^H] &= \mathbf R _{H_p}\operatorname{diag}(\mathbf X_p)^H,\\
\mathbb E[\mathbf Y_p\mathbf Y_p^H] &=  \operatorname{diag}(\mathbf X_p)\mathbf R_{H_p}\operatorname{diag}(\mathbf X_p)^H
+\sigma^2\mathbf I,
\end{aligned} 
$$

其中


$$
\mathbf R _{H_p}=\mathbb E[\mathbf H_p\mathbf H_p^H]
$$

为导频位置上的信道自相关矩阵。

因此：

$$
\boxed{
\hat{\mathbf H}_p = \mathbf R _{H_p}\operatorname{diag}(\mathbf X_p)^H
\left(
\operatorname{diag}(\mathbf X_p)\mathbf R _{H_p}\operatorname{diag}(\mathbf X_p)^H+\sigma^2\mathbf I
\right)^{-1}
\mathbf Y_p
}
$$

MMSE 估计利用了二阶统计信息，因此性能通常优于 LS，但需要已知或估计 $\mathbf R_{H_p}$ 和噪声功率。

估计误差：

$$
\begin{align}
\epsilon &= \mathrm{E}\left\{ \lVert \hat{\mathbf{H}}_p - \mathbf{H}_p \rVert^2_F \right\} \\
& = \sigma^2 Tr(R_h \left(\operatorname{diag}(\mathbf{X}_p)^H \operatorname{diag}(\mathbf{X}_p) R_h + \sigma^2 \mathbf{I}_{N_p} \right)^{-1})
\end{align} 
$$

***

### 3.1.3 3.1.3 Linear MMSE (LMMSE)

在 OFDM 文献中，LMMSE 一般指在**线性估计器**类中求最优 MMSE 解。对于当前模型，上式本身已经是 LMMSE 解。

若进一步假设导频符号模值相同：


$$
\operatorname{diag}(\mathbf X_p)\operatorname{diag}(\mathbf X_p)^H = \sigma_x^2\mathbf I,

$$

则上式可化简为：

$$
\boxed{
\hat{\mathbf H}_p = \mathbf R _{H_p}
\left(
\mathbf R _{H_p}+\frac{\sigma^2}{\sigma_x^2}\mathbf I
\right)^{-1}
\operatorname{diag}(\mathbf X_p)^{-1}\mathbf Y_p
}
$$

这就是常见的 LMMSE 信道估计表达式。

## 3.2 3.2 插值

插值用于根据导频位置上的信道估计恢复数据位置上的信道。
对于单个 OFDM 符号，通常只需要**频域插值**；对于多个符号上的散布导频，还需要**时间-频率二维插值**。

设：

- 频域导频间隔为 $S_f$
- 时间域导频间隔为 $S_t$

则导频位于：

- 频域索引：$0,S_f,2S_f,\dots$
- 时间索引：$0,S_t,2S_t,\dots$

***

### 3.2.1 3.2.1 线性插值

#### 3.2.1.1 一维线性插值

设 $gS_f\le k < (g+1)S_f$，定义


$$
\beta = \frac{k-gS_f}{S_f},\qquad 0\le \beta <1.
$$

则频域线性插值为：


$$
\boxed{ \hat H[m,k] =  (1-\beta)H_p[m,g]+\beta H_p[m,g+1]
}
$$

这里 $H_p[m,g]$ 表示第 $m$ 个 OFDM 符号上第 $g$ 个导频频率位置处的信道估计。

#### 3.2.1.2 二维线性插值

设：


$$
qS_t \le m < (q+1)S_t,\qquad
gS_f \le k < (g+1)S_f
$$

并定义

$$
\alpha = \frac{m-qS_t}{S_t},\qquad
\beta = \frac{k-gS_f}{S_f}.
$$

则双线性插值为：

$$
\boxed{
\begin{aligned}
\hat H[m,k]
=&\ (1-\alpha)(1-\beta)H_p[q,g]
+(1-\alpha)\beta H_p[q,g+1] \\
&+\alpha(1-\beta)H_p[q+1,g]
+\alpha\beta H_p[q+1,g+1]
\end{aligned}
}
$$

***

### 3.2.2 3.2.2 二阶插值

#### 3.2.2.1 一维二阶插值

设 $gS_f < k \le (g+1)S_f$，定义


$$
a=\frac{k-gS_f}{S_f},\qquad 0<a\le 1.
$$

采用三个导频点 $g-1,g,g+1$ 做二阶 Lagrange 插值，则：

$$
\boxed{
\hat H[m,k] = c _{-1}^{(f)}H_p[m,g-1]
+
c_0^{(f)}H_p[m,g]
+
c _{+1}^{(f)}H_p[m,g+1]
}
$$

其中

$$
c _{-1}^{(f)}=\frac{a(a-1)}{2},
\qquad
c_0^{(f)}=1-a^2=-(a-1)(a+1),
\qquad
c _{+1}^{(f)}=\frac{a(a+1)}{2}.
$$

并且

$$
c _{-1}^{(f)}+c_0^{(f)}+c _{+1}^{(f)}=1.
$$

> 注：其他权重也可以，只要让权重之和为1”。

在频域边界（如 $g=0$）无法取得 $g-1$ 时，通常退化为线性插值或采用单边二阶插值。

#### 3.2.2.2 二维二阶插值

定义频域插值系数：

$$
c_{-1}^{(f)}(a)=\frac{a(a-1)}{2},
\qquad
c_{0}^{(f)}(a)=1-a^2,
\qquad
c_{+1}^{(f)}(a)=\frac{a(a+1)}{2}.
$$

定义时间域插值系数：

$$
c_{-1}^{(t)}(b)=\frac{b(b-1)}{2},
\qquad
c_{0}^{(t)}(b)=1-b^2,
\qquad
c _{+1}^{(t)}(b)=\frac{b(b+1)}{2}.
$$

其中

$$
a=\frac{k-gS_f}{S_f},\qquad
b=\frac{m-qS_t}{S_t}.
$$

则二维二阶插值可写成张量积形式：


$$
\boxed{
\hat H[m,k]=\sum_{u=-1}^{1}\sum_{v=-1}^{1}
c_u^{(f)}(a),c_v^{(t)}(b),H_p[q+v,g+u]
}
$$

这相当于用周围 $3\times 3$ 个导频点联合进行插值。

***

### 3.2.3 3.2.3 Wiener (LMMSE) 插值

#### 3.2.3.1 一维插值

设导频位置上的估计为 $\hat{\mathbf H}_p$，希望恢复数据位置上的信道 $\hat{\mathbf H}_d$。假设采用线性估计器：


$$
\hat{\mathbf H}d=\mathbf A\hat{\mathbf H} p,
\qquad
\mathbf A\in\mathbb C^{N_d\times N_p}.
$$

目标是最小化均方误差：

$$
\mathbf A^\star = \arg\min _{\mathbf A}
\mathbb E\left[
\left|
\mathbf H_d-\mathbf A\hat{\mathbf H}p
\right| 2^2
\right].
$$

利用正交原理可得 Wiener-Hopf 方程：


$$
\mathbf A = \mathbf R_{\hat p\hat p}  \mathbf R_{d\hat p},
$$

其中

$$
\mathbf{R}_{d\hat p}=\mathbb E[\mathbf{H}_d \hat{\mathbf{H}}_p^H],\qquad
\mathbf R_{\hat p\hat p}=\mathbb E[\hat{\mathbf H}_p\hat{\mathbf H}_p^H].
$$

因此 Wiener 插值矩阵为：

$$
\boxed{
\mathbf A^\star=\mathbf R {d\hat p}\mathbf R{\hat p\hat p}^{-1}
}
$$

代回可得：

$$
\boxed{
\hat{\mathbf H}_d = \mathbf R_{d\hat p}\mathbf R_{\hat p\hat p}^{-1}\hat{\mathbf H}_p
}
$$

#### 3.2.3.2 导频估计噪声的影响

若导频位置上的 LS 估计可表示为


$$
\hat{\mathbf H}_p=\mathbf H_p+\mathbf e_p,
$$

其中噪声与真实信道不相关，且


$$
\mathbb E[\mathbf e_p\mathbf e_p^H]=\sigma_e^2\mathbf I,
$$

则有：

$$
\mathbf R_{d\hat p}=\mathbf R_{dp},\qquad \mathbf R _{\hat p\hat p} = \mathbf R _{pp}+\sigma_e^2\mathbf I.

$$

因此更常用的 Wiener 插值公式为：

$$
\boxed{
\hat{\mathbf H}_d = \mathbf R _{dp}
\left(
\mathbf R_{pp}+\sigma_e^2\mathbf I
\right)^{-1}
\hat{\mathbf H}_ p
}
$$

#### 3.2.3.3 二维插值

二维 Wiener 插值在形式上与一维并无本质区别。关键区别在于：

- 一维情况下，$\mathbf H_p$ 和 $\mathbf H_d$ 通常只包含单个维度上的位置；
- 二维情况下，需要把时频网格上的导频点和数据点统一向量化，再构造对应的二维相关矩阵。

因此，二维 Wiener 插值依然可写为：

$$
\boxed{
\hat{\mathbf H}_d = \mathbf R _{dp}
\left(
\mathbf R_{pp}+\sigma_e^2\mathbf I
\right)^{-1}
\hat{\mathbf H}_ p
}
$$

只是这里的 $\mathbf R_{pp}$ 和 $\mathbf R_{dp}$ 是基于二维时频坐标构造的相关矩阵。

***

### 3.2.4 3.2.4 DFT 插值

DFT 插值基于一个关键先验：

> **时域 CIR 长度有限。**

若导频频域采样足够密，则可以先由导频位置上的频域响应恢复 CIR，再由恢复的 CIR 重建全部子载波上的 CFR。

#### 3.2.4.1 3.2.4.1 一维频域 DFT 插值

考虑单个 OFDM 符号(省略符号索引m)、梳状导频。设导频间隔为 $S_f$，导频数为$N_p = \frac{N}{S_f}$（为简洁起见，假设可整除）。

设导频位置为：$k = gS_f,\qquad g=0,1,\dots,N_p-1.$

导频位置上的 LS 信道估计记为：$\hat H_p[g].$

##### 3.2.4.1.1 Step 1：由导频 CFR 恢复低维 CIR

由频域模型

$$
H[gS_f] =  \sum _{n=0}^{L_h-1} h[n] e^{-j\frac{2\pi}{N}(gS_f)n} = \sum _{n=0}^{L_h-1} h[n] e^{-j\frac{2\pi}{N_p}gn},
$$

其中 $L_h$ 为离散 CIR 的有效长度。只要$L_h \le N_p$, 则导频位置上的频域样本可唯一确定长度不超过 $N_p$ 的 CIR。因此可对导频估计做 $N_p$ 点 IDFT：

$$
\hat h_p[n] = \frac{1}{N_p}
\sum_{g=0}^{N_p-1}
\hat H_p[g] e^{j\frac{2\pi}{N_p}gn},
\qquad n=0,1,\dots,N_p-1.
$$

##### 3.2.4.1.2 Step 2：时域截断 / 去噪

由于噪声存在，$\hat h_p[n]$ 在所有 tap 上都会有扰动。若已知 CP 长度为 $N_{cp}$，并且有效 CIR 长度满足

$L_h \le N_{cp} \le N_p$, 则可进行时域截断：

$$
\tilde h[n] =  \begin{cases}
\hat h_p[n], & 0\le n\le N_{cp}-1
0, & N _{cp}\le n\le N_p-1
\end{cases}
$$

这一步利用了“有效 CIR 长度有限”的先验，能够抑制高时延区域的噪声。

##### 3.2.4.1.3 Step 3：零填充并恢复全部子载波的 CFR

将截断后的 CIR 零填充到长度 $N$：

$$
\bar h[n] =  \begin{cases}
\tilde h[n], & 0\le n\le N_p-1\
0, & N_p\le n\le N-1
\end{cases}
$$

然后做 $N$ 点 DFT：

$$
\boxed{
\hat H[k] =  \sum _{n=0}^{N-1}\bar h[n]e^{-j\frac{2\pi}{N}kn},
\qquad k=0,1,\dots,N-1
}
$$

从而得到所有子载波上的信道估计。

***

#### 3.2.4.2 3.2.4.2 二维 DFT 插值

对于多个 OFDM 符号上的导频，最常见的二维 DFT 插值并不是直接在整个时频平面上做二维 DFT，而是采用**分离式处理**：

1. **频率方向**：在每个导频 OFDM 符号上，对梳状导频先做一维 DFT 插值，恢复该符号上全部子载波的 CFR；
2. **时间方向**：在导频 OFDM 符号之间，对各子载波上的信道做线性、二阶或 Wiener 时间插值。

也就是说，常见实现是：

$$
\text{DFT (频率方向)} + \text{时间域插值}
$$

这种处理方式与前面线性插值、二阶插值、Wiener 插值的逻辑保持一致，也最符合工程实现习惯。

更一般地，如果信道在时延-时间域或时延-多普勒域具有更强的低维结构，也可以进一步构造二维 DFT/BEM 类方法，但那已经超出了本文的基础插值范畴。

***

### 3.2.5 3.2.5 DFT 插值的成立条件

DFT 插值成立的关键条件是：

$$
\boxed{
L_h \le N_p
}
$$

实际实现时，常用更保守、也更工程化的写法：

$$
\boxed{
N _{cp} \le N_p
}
$$

因为若有效 CIR 超过 $N_p$，则对 $N_p$ 点 IDFT 得到的时域响应会发生混叠，无法唯一恢复原始 CIR。此时 DFT 插值性能会明显下降。

# 4 4 References

1. [https://pmc.ncbi.nlm.nih.gov/articles/PMC8309705/](https://pmc.ncbi.nlm.nih.gov/articles/PMC8309705/ "https://pmc.ncbi.nlm.nih.gov/articles/PMC8309705/")
2. LS / MMSE 推导参考： [https://www.cnblogs.com/hjd21/p/16634313.html](https://www.cnblogs.com/hjd21/p/16634313.html "https://www.cnblogs.com/hjd21/p/16634313.html")
3. A. Hjorungnes and D. Gesbert, "Complex-Valued Matrix Differentiation: Techniques and Key Results," IEEE Transactions on Signal Processing, vol. 55, no. 6, pp. 2740-2746, June 2007.
