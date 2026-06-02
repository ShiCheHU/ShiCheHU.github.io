# SISO OFDM Channel Estimation

## Table of Contents

- [1 Overview](#1-overview)

- [2 System Model](#2-system-model)
  - [2.1 Channel Model](#21-channel-model)
    - [2.1.1 Discrete-Domain Representation](#211-discrete-domain-representation)

  - [2.2 Received Signal Model](#22-received-signal-model)

- [3 Channel Estimation and Interpolation](#3-channel-estimation-and-interpolation)
  - [3.1 Channel Estimation at Pilot Positions](#31-channel-estimation-at-pilot-positions)
    - [3.1.1 Least Squares (LS)](#311-least-squares-ls)

    - [3.1.2 Minimum Mean Square Error (MMSE)](#312-minimum-mean-square-error-mmse)

    - [3.1.3 Linear MMSE (LMMSE)](#313-linear-mmse-lmmse)

  - [3.2 Interpolation](#32-interpolation)
    - [3.2.1 Linear Interpolation](#321-linear-interpolation)
      - [3.2.1.1 One-Dimensional Linear Interpolation](#3211-one-dimensional-linear-interpolation)

      - [3.2.1.2 Two-Dimensional Linear Interpolation](#3212-two-dimensional-linear-interpolation)

    - [3.2.2 Second-Order Interpolation](#322-second-order-interpolation)
      - [3.2.2.1 One-Dimensional Second-Order Interpolation](#3221-one-dimensional-second-order-interpolation)

      - [3.2.2.2 Two-Dimensional Second-Order Interpolation](#3222-two-dimensional-second-order-interpolation)

    - [3.2.3 Wiener (LMMSE) Interpolation](#323-wiener-lmmse-interpolation)
      - [3.2.3.1 One-Dimensional Interpolation](#3231-one-dimensional-interpolation)

      - [3.2.3.2 Impact of Pilot Estimation Noise](#3232-impact-of-pilot-estimation-noise)

      - [3.2.3.3 Two-Dimensional Interpolation](#3233-two-dimensional-interpolation)

    - [3.2.4 DFT Interpolation](#324-dft-interpolation)
      - [3.2.4.1 One-Dimensional Frequency-Domain DFT Interpolation](#3241-one-dimensional-frequency-domain-dft-interpolation)
        - [3.2.4.1.1 Step 1: Recover a Low-Dimensional CIR from Pilot CFR Samples](#32411-step-1-recover-a-low-dimensional-cir-from-pilot-cfr-samples)

        - [3.2.4.1.2 Step 2: Time-Domain Truncation / Denoising](#32412-step-2-time-domain-truncation--denoising)

        - [3.2.4.1.3 Step 3: Zero Padding and Full-Subcarrier CFR Reconstruction](#32413-step-3-zero-padding-and-full-subcarrier-cfr-reconstruction)

      - [3.2.4.2 Two-Dimensional DFT Interpolation](#3242-two-dimensional-dft-interpolation)

    - [3.2.5 Conditions for DFT Interpolation](#325-conditions-for-dft-interpolation)

- [4 References](#4-references)


# 1 Overview

Channel estimation is closely related to the channel model, pilot pattern, and receiver complexity. In practice, the design of channel estimation and interpolation algorithms usually depends on:

- time selectivity and frequency selectivity of the channel;
- the time-frequency distribution of pilots;
- the desired complexity-performance tradeoff.

In an OFDM system, channel recovery is usually divided into two steps:

1. **Channel estimation at pilot positions**: estimate the channel at the time-frequency resource elements that carry pilots.
2. **Interpolation / extrapolation**: recover the channel on data subcarriers from the pilot estimates.

Common pilot patterns mainly include:

- **Block-type pilots**: all subcarriers in some OFDM symbols carry pilots. This is suitable for channels with weak time selectivity and strong frequency selectivity.
- **Comb-type pilots**: only a subset of subcarriers in each or some OFDM symbols carry pilots. This is suitable for channels with stronger time selectivity and weaker frequency selectivity.

For multi-symbol scenarios, time-domain correlation is often used for two-dimensional interpolation, channel tracking, or prediction. This note focuses only on **pilot-aided channel estimation and interpolation**.

# 2 System Model

## 2.1 Channel Model

Consider a single-input single-output (SISO) OFDM system. Define:

- number of paths: $L$
- maximum delay: $\tau_{\max}$
- maximum velocity: $v_{\max}$
- maximum Doppler shift: $\nu_{\max} = \dfrac{v_{\max}}{\lambda},\ \lambda = \dfrac{c}{f_c}$
- carrier frequency: $f_c$
- subcarrier spacing: $\Delta f$
- number of subcarriers per OFDM symbol: $N$
- useful OFDM symbol duration: $T_u = 1/\Delta f$
- cyclic prefix length: $T_{cp}$
- total OFDM symbol duration: $T_{\mathrm{sym}} = T_u + T_{cp}$
- number of OFDM symbols: $M$, where $m=0,1,\dots,M-1$

To avoid ISI, one usually requires $T_{cp} > \tau_{\max}$.

Let the time corresponding to the $m$-th OFDM symbol be $t_m = mT_{\mathrm{sym}}$.

Under the quasi-static approximation within the $m$-th OFDM symbol, the continuous-domain CIR can be written as:

$$
\begin{aligned}
h_m(\tau)
&\triangleq h(t_m,\tau) \\
&= \sum _{\ell=0}^{L-1} \beta_\ell e^{j2\pi \nu _\ell t_m}\delta(\tau-\tau _\ell) \\
&= \sum_{\ell=0}^{L-1} \alpha _{\ell,m}\delta(\tau-\tau_\ell),
\end{aligned}
$$

where:

- $\tau_\ell$ is the delay of the $\ell$-th path;
- $\beta_\ell$ is the initial complex gain of the $\ell$-th path;
- $\nu_\ell$ is the Doppler shift of the $\ell$-th path;
- $\alpha_{\ell,m} \triangleq \beta_\ell e^{j2\pi \nu_\ell t_m}$ is the equivalent complex gain of the $\ell$-th path in the $m$-th OFDM symbol.

### 2.1.1 Discrete-Domain Representation

The sampling interval is $T_{\mathrm{samp}} = \frac{1}{N\Delta f}$.

The corresponding discrete CIR is denoted by $h[m,n],\qquad n=0,1,\dots,N-1$.

In an integer-delay tap model, $h[m,n]$ is nonzero only on a finite number of taps.

The corresponding CFR is defined as the **DFT** of the discrete CIR:

$$
H[m,k]=\sum _{n=0}^{N-1} h[m,n] e^{-j\frac{2\pi}{N}kn},\qquad k=0,1,\dots,N-1.
$$

The inverse transform is the **IDFT**:

$$
h[m,n]=\frac{1}{N}\sum _{k=0}^{N-1} H[m,k] e^{j\frac{2\pi}{N}kn},\qquad n=0,1,\dots,N-1.
$$

## 2.2 Received Signal Model

Let the transmitted frequency-domain symbol on the $m$-th OFDM symbol be $X_m[k],\qquad k=0,1,\dots,N-1$, where $X_m[k]\in \mathcal Q$ and $\mathcal Q$ is the modulation constellation, such as QPSK, 16-QAM, or 64-QAM.

After IDFT, the time-domain transmitted signal is $x_m[n]$. Assume the symbol power is normalized, i.e., $\mathbb E{|X_m[k]|^2}=1$.

If the CP is long enough, after CP removal the linear convolution can be treated as circular convolution. The received time-domain signal is:

$$
y_m[n]=h_m[n]\circledast x_m[n]+w_m[n],\qquad n=0,1,\dots,N-1,
$$

where:

- $\circledast$ denotes circular convolution;
- $w_m[n]$ is complex Gaussian white noise with power $\sigma^2$.

In matrix form:

$$
\mathbf y_m=\mathrm{circ}(\mathbf h_m)\mathbf x_m+\mathbf w_m
$$

or equivalently:

$$
\mathbf y_m=\mathrm{circ}(\mathbf x_m)\mathbf h_m+\mathbf w_m.
$$

The two expressions only differ in which vector is used to generate the circular convolution matrix.

After applying the DFT to the received signal, the frequency-domain input-output relation is:

$$
Y_m[k]=H_m[k]X_m[k]+W_m[k],\qquad k=0,1,\dots,N-1.
$$

In matrix form:

$$
\mathbf Y_m=\operatorname{diag}(\mathbf X_m)\mathbf H_m+\mathbf W_m.
$$

Divide the frequency-domain positions into a pilot set $\mathcal P$ and a data set $\mathcal D$:

$$
\begin{aligned}
\mathbf Y _{m,p}&=\operatorname{diag}(\mathbf X_{m,p})\mathbf H _{m,p}+\mathbf W_{m,p},\\
\mathbf Y _{m,d}&=\operatorname{diag}(\mathbf X_{m,d})\mathbf H _{m,d}+\mathbf W_{m,d}.
\end{aligned}
$$

If pilots across multiple OFDM symbols are vectorized, this can be written as:

$$
\begin{aligned}
\mathbf Y_p&=\operatorname{diag}(\mathbf X_p)\mathbf H_p+\mathbf W_p, \\
\mathbf Y_d&=\operatorname{diag}(\mathbf X_d)\mathbf H_d+\mathbf W_d,
\end{aligned}
$$

where:

- $\mathbf X_p,\mathbf H_p,\mathbf Y_p,\mathbf W_p\in\mathbb C^{N_p\times 1}$;
- $\mathbf X_d,\mathbf H_d,\mathbf Y_d,\mathbf W_d\in\mathbb C^{N_d\times 1}$.

Here $N_p$ is the total number of pilot positions in the considered time-frequency region, and $N_d$ is the total number of data positions.

# 3 Channel Estimation and Interpolation

Channel recovery is usually divided into:

1. **channel estimation at pilot positions**
2. **interpolation from pilot positions to data positions**

This note introduces three common pilot-domain channel estimation algorithms:

- LS
- MMSE
- LMMSE

and four common interpolation algorithms:

- linear interpolation
- second-order interpolation
- Wiener (LMMSE) interpolation
- DFT interpolation

## 3.1 Channel Estimation at Pilot Positions

### 3.1.1 Least Squares (LS)

At pilot positions, LS estimation minimizes the squared residual:

$$
\hat{\mathbf H}_p =  \arg\min _{\mathbf H_p\in\mathbb C^{N_p\times 1}}
\left\|
\mathbf Y_p-\operatorname{diag}(\mathbf X_p)\mathbf H_p
\right\|_2^2.
$$

The closed-form solution is:

$$
\hat{\mathbf H}_p = \left(
\operatorname{diag}(\mathbf X_p)^H
\operatorname{diag}(\mathbf X_p)
\right)^{-1}
\operatorname{diag}(\mathbf X_p)^H\mathbf Y_p.
$$

Since $\operatorname{diag}(\mathbf X_p)$ is diagonal, as long as all pilot symbols are nonzero, the solution simplifies to:

$$
\boxed{
\hat{\mathbf H}_p=\operatorname{diag}(\mathbf X_p)^{-1}\mathbf Y_p
}
$$

That is, the received pilot symbols are divided by the known transmitted pilot symbols element by element.

If all pilots have the same power, $|X_p[i]|^2=\sigma_x^2$, the LS estimation MSE is inversely proportional to SNR. LS is simple to implement, but it does not use channel statistics and therefore has weaker noise robustness.

The estimation error is:

$$
\begin{aligned}
\epsilon &= \mathbb{E}\left\{ \lVert \hat{\mathbf{H}}_p - \mathbf{H}_p \rVert^2_2 \right\} \\
&= \mathbb{E}\left\{ \lVert \operatorname{diag}(\mathbf{X}_p)^{-1}\mathbf W_p \rVert^2_2 \right\} \\
&= \frac{\sigma^2}{\sigma^2_x}
= \frac{1}{\mathrm{SNR}}.
\end{aligned}
$$

***

### 3.1.2 Minimum Mean Square Error (MMSE)

MMSE estimation minimizes the mean square error at pilot positions. Assume a linear estimator:

$$
\hat{\mathbf H}_p=\mathbf Q_p \mathbf Y_p,
$$

where $\mathbf Q_p\in\mathbb C^{N_p\times N_p}$. The optimization objective is:

$$
\mathbf Q_p^\star =  \arg\min _{\mathbf Q_p}
\mathbb E\left\{
\left\|
\mathbf H_p-\mathbf Q_p\mathbf Y_p
\right\|_2^2
\right\}.
$$

By the orthogonality principle, the optimal solution is:

$$
\mathbf Q_p^\star = \mathbb E[\mathbf H_p\mathbf Y_p^H]
\left(
\mathbb E[\mathbf Y_p\mathbf Y_p^H]
\right)^{-1}.
$$

From the model

$$
\mathbf Y_p=\operatorname{diag}(\mathbf X_p)\mathbf H_p+\mathbf W_p,
$$

we have:

$$
\begin{aligned}
\mathbb E[\mathbf H_p\mathbf Y_p^H] &= \mathbf R _{H_p}\operatorname{diag}(\mathbf X_p)^H,\\
\mathbb E[\mathbf Y_p\mathbf Y_p^H] &=  \operatorname{diag}(\mathbf X_p)\mathbf R_{H_p}\operatorname{diag}(\mathbf X_p)^H
+\sigma^2\mathbf I,
\end{aligned}
$$

where

$$
\mathbf R _{H_p}=\mathbb E[\mathbf H_p\mathbf H_p^H]
$$

is the channel autocorrelation matrix at pilot positions.

Therefore:

$$
\boxed{
\hat{\mathbf H}_p = \mathbf R _{H_p}\operatorname{diag}(\mathbf X_p)^H
\left(
\operatorname{diag}(\mathbf X_p)\mathbf R _{H_p}\operatorname{diag}(\mathbf X_p)^H+\sigma^2\mathbf I
\right)^{-1}
\mathbf Y_p
}
$$

MMSE estimation uses second-order statistics, so it usually outperforms LS, but it requires known or estimated $\mathbf R_{H_p}$ and noise power.

The estimation error is:

$$
\begin{aligned}
\epsilon
&= \mathbb E\left\{ \lVert \hat{\mathbf H}_p-\mathbf H_p\rVert_2^2 \right\} \\
&= \operatorname{Tr}\left(
\mathbf R_{H_p}
-\mathbf R_{H_p}\operatorname{diag}(\mathbf X_p)^H
\left(
\operatorname{diag}(\mathbf X_p)\mathbf R_{H_p}\operatorname{diag}(\mathbf X_p)^H
+\sigma^2\mathbf I_{N_p}
\right)^{-1}
\operatorname{diag}(\mathbf X_p)\mathbf R_{H_p}
\right).
\end{aligned}
$$

If $\mathbf R_{H_p}$ is nonsingular, the same expression can be written via the matrix inversion lemma as:

$$
\boxed{
\epsilon
=\operatorname{Tr}\left[
\sigma^2\mathbf R_{H_p}
\left(
\operatorname{diag}(\mathbf X_p)^H\operatorname{diag}(\mathbf X_p)\mathbf R_{H_p}
+\sigma^2\mathbf I_{N_p}
\right)^{-1}
\right]
}
$$

***

### 3.1.3 Linear MMSE (LMMSE)

In OFDM literature, LMMSE usually refers to the MMSE solution within the class of linear estimators. For the current model, the expression above is already the LMMSE solution.

If the pilot symbols further have the same magnitude:

$$
\operatorname{diag}(\mathbf X_p)\operatorname{diag}(\mathbf X_p)^H = \sigma_x^2\mathbf I,
$$

then the estimator can be simplified to:

$$
\boxed{
\hat{\mathbf H}_p = \mathbf R _{H_p}
\left(
\mathbf R _{H_p}+\frac{\sigma^2}{\sigma_x^2}\mathbf I
\right)^{-1}
\operatorname{diag}(\mathbf X_p)^{-1}\mathbf Y_p
}
$$

This is the commonly used LMMSE channel estimation expression.

## 3.2 Interpolation

Interpolation is used to recover the channel at data positions from the channel estimates at pilot positions.
For a single OFDM symbol, frequency-domain interpolation is usually sufficient. For scattered pilots across multiple symbols, time-frequency two-dimensional interpolation is also required.

Let:

- frequency-domain pilot spacing be $S_f$
- time-domain pilot spacing be $S_t$

Then the pilots are located at:

- frequency indices: $0,S_f,2S_f,\dots$
- time indices: $0,S_t,2S_t,\dots$

***

### 3.2.1 Linear Interpolation

#### 3.2.1.1 One-Dimensional Linear Interpolation

Assume $gS_f\le k < (g+1)S_f$, and define:

$$
\beta = \frac{k-gS_f}{S_f},\qquad 0\le \beta <1.
$$

The frequency-domain linear interpolation is:

$$
\boxed{ \hat H[m,k] =  (1-\beta)H_p[m,g]+\beta H_p[m,g+1]
}
$$

Here $H_p[m,g]$ denotes the channel estimate at the $g$-th pilot frequency position in the $m$-th OFDM symbol.

#### 3.2.1.2 Two-Dimensional Linear Interpolation

Assume:

$$
qS_t \le m < (q+1)S_t,\qquad
gS_f \le k < (g+1)S_f
$$

and define:

$$
\alpha = \frac{m-qS_t}{S_t},\qquad
\beta = \frac{k-gS_f}{S_f}.
$$

The bilinear interpolation is:

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

### 3.2.2 Second-Order Interpolation

#### 3.2.2.1 One-Dimensional Second-Order Interpolation

Assume $gS_f < k \le (g+1)S_f$, and define:

$$
a=\frac{k-gS_f}{S_f},\qquad 0<a\le 1.
$$

Using three pilot points $g-1,g,g+1$ for second-order Lagrange interpolation gives:

$$
\boxed{
\hat H[m,k] = c _{-1}^{(f)}H_p[m,g-1]
+
c_0^{(f)}H_p[m,g]
+
c _{+1}^{(f)}H_p[m,g+1]
}
$$

where

$$
c _{-1}^{(f)}=\frac{a(a-1)}{2},
\qquad
c_0^{(f)}=1-a^2=-(a-1)(a+1),
\qquad
c _{+1}^{(f)}=\frac{a(a+1)}{2}.
$$

and

$$
c _{-1}^{(f)}+c_0^{(f)}+c _{+1}^{(f)}=1.
$$

> Other weights can also be used as long as their sum is 1.

At frequency-domain boundaries, such as $g=0$, $g-1$ is unavailable. One usually falls back to linear interpolation or uses a one-sided second-order interpolation.

#### 3.2.2.2 Two-Dimensional Second-Order Interpolation

Define frequency-domain interpolation coefficients:

$$
c_{-1}^{(f)}(a)=\frac{a(a-1)}{2},
\qquad
c_{0}^{(f)}(a)=1-a^2,
\qquad
c_{+1}^{(f)}(a)=\frac{a(a+1)}{2}.
$$

Define time-domain interpolation coefficients:

$$
c_{-1}^{(t)}(b)=\frac{b(b-1)}{2},
\qquad
c_{0}^{(t)}(b)=1-b^2,
\qquad
c _{+1}^{(t)}(b)=\frac{b(b+1)}{2}.
$$

where

$$
a=\frac{k-gS_f}{S_f},\qquad
b=\frac{m-qS_t}{S_t}.
$$

The two-dimensional second-order interpolation can be written as a tensor product:

$$
\boxed{
\hat H[m,k]=\sum_{u=-1}^{1}\sum_{v=-1}^{1}
c_u^{(f)}(a)c_v^{(t)}(b)H_p[q+v,g+u]
}
$$

This is equivalent to jointly interpolating with the surrounding $3\times 3$ pilot points.

***

### 3.2.3 Wiener (LMMSE) Interpolation

#### 3.2.3.1 One-Dimensional Interpolation

Let the channel estimate at pilot positions be $\hat{\mathbf H}_p$, and suppose we want to recover the channel at data positions $\hat{\mathbf H}_d$. Assume a linear estimator:

$$
\hat{\mathbf H}_d=\mathbf A\hat{\mathbf H}_p,
\qquad
\mathbf A\in\mathbb C^{N_d\times N_p}.
$$

The objective is to minimize the mean square error:

$$
\mathbf A^\star = \arg\min _{\mathbf A}
\mathbb E\left[
\left\|
\mathbf H_d-\mathbf A\hat{\mathbf H}_p
\right\|_2^2
\right].
$$

By the orthogonality principle, the Wiener-Hopf solution is:

$$
\mathbf A^\star = \mathbf R_{d\hat p}\mathbf R_{\hat p\hat p}^{-1},
$$

where

$$
\mathbf{R}_{d\hat p}=\mathbb E[\mathbf{H}_d \hat{\mathbf{H}}_p^H],\qquad
\mathbf R_{\hat p\hat p}=\mathbb E[\hat{\mathbf H}_p\hat{\mathbf H}_p^H].
$$

Therefore, the Wiener interpolation matrix is:

$$
\boxed{
\mathbf A^\star=\mathbf R_{d\hat p}\mathbf R_{\hat p\hat p}^{-1}
}
$$

Substituting it back gives:

$$
\boxed{
\hat{\mathbf H}_d = \mathbf R_{d\hat p}\mathbf R_{\hat p\hat p}^{-1}\hat{\mathbf H}_p
}
$$

#### 3.2.3.2 Impact of Pilot Estimation Noise

If the LS estimate at pilot positions can be written as:

$$
\hat{\mathbf H}_p=\mathbf H_p+\mathbf e_p,
$$

where the estimation noise is uncorrelated with the true channel and

$$
\mathbb E[\mathbf e_p\mathbf e_p^H]=\sigma_e^2\mathbf I,
$$

then:

$$
\mathbf R_{d\hat p}=\mathbf R_{dp},\qquad
\mathbf R _{\hat p\hat p} = \mathbf R _{pp}+\sigma_e^2\mathbf I.
$$

A more commonly used Wiener interpolation formula is therefore:

$$
\boxed{
\hat{\mathbf H}_d = \mathbf R _{dp}
\left(
\mathbf R_{pp}+\sigma_e^2\mathbf I
\right)^{-1}
\hat{\mathbf H}_p
}
$$

#### 3.2.3.3 Two-Dimensional Interpolation

Two-dimensional Wiener interpolation has the same form as the one-dimensional case. The main difference is:

- in the one-dimensional case, $\mathbf H_p$ and $\mathbf H_d$ usually contain positions along only one dimension;
- in the two-dimensional case, the pilot and data positions on the time-frequency grid must be vectorized, and the corresponding two-dimensional correlation matrices must be constructed.

Thus, two-dimensional Wiener interpolation can still be written as:

$$
\boxed{
\hat{\mathbf H}_d = \mathbf R _{dp}
\left(
\mathbf R_{pp}+\sigma_e^2\mathbf I
\right)^{-1}
\hat{\mathbf H}_p
}
$$

Here $\mathbf R_{pp}$ and $\mathbf R_{dp}$ are correlation matrices built from two-dimensional time-frequency coordinates.

***

### 3.2.4 DFT Interpolation

DFT interpolation is based on a key prior:

> **The time-domain CIR has finite length.**

If the pilot frequency-domain sampling is dense enough, one can first recover the CIR from the CFR samples at pilot positions, and then reconstruct the CFR on all subcarriers from the recovered CIR.

#### 3.2.4.1 One-Dimensional Frequency-Domain DFT Interpolation

Consider a single OFDM symbol and omit the symbol index $m$. Assume comb-type pilots. Let the pilot spacing be $S_f$, and let the number of pilots be $N_p = \frac{N}{S_f}$, assuming divisibility for simplicity.

The pilot positions are:

$$
k = gS_f,\qquad g=0,1,\dots,N_p-1.
$$

The LS channel estimate at pilot positions is denoted by $\hat H_p[g]$.

##### 3.2.4.1.1 Step 1: Recover a Low-Dimensional CIR from Pilot CFR Samples

From the frequency-domain model:

$$
H[gS_f] =  \sum _{n=0}^{L_h-1} h[n] e^{-j\frac{2\pi}{N}(gS_f)n}
= \sum _{n=0}^{L_h-1} h[n] e^{-j\frac{2\pi}{N_p}gn},
$$

where $L_h$ is the effective length of the discrete CIR. As long as $L_h \le N_p$, the frequency-domain samples at pilot positions uniquely determine a CIR with length no greater than $N_p$. Therefore, an $N_p$-point IDFT can be applied to the pilot estimates:

$$
\hat h_p[n] = \frac{1}{N_p}
\sum_{g=0}^{N_p-1}
\hat H_p[g] e^{j\frac{2\pi}{N_p}gn},
\qquad n=0,1,\dots,N_p-1.
$$

##### 3.2.4.1.2 Step 2: Time-Domain Truncation / Denoising

Because noise is present, $\hat h_p[n]$ is perturbed on all taps. If the CP length $N_{cp}$ is known and the effective CIR length satisfies

$$
L_h \le N_{cp} \le N_p,
$$

then time-domain truncation can be applied:

$$
\tilde h[n] =
\begin{cases}
\hat h_p[n], & 0\le n\le N_{cp}-1,\\
0, & N_{cp}\le n\le N_p-1.
\end{cases}
$$

This step uses the finite effective CIR length prior and suppresses noise in high-delay taps.

##### 3.2.4.1.3 Step 3: Zero Padding and Full-Subcarrier CFR Reconstruction

Zero-pad the truncated CIR to length $N$:

$$
\bar h[n] =
\begin{cases}
\tilde h[n], & 0\le n\le N_p-1,\\
0, & N_p\le n\le N-1.
\end{cases}
$$

Then apply an $N$-point DFT:

$$
\boxed{
\hat H[k] =  \sum _{n=0}^{N-1}\bar h[n]e^{-j\frac{2\pi}{N}kn},
\qquad k=0,1,\dots,N-1
}
$$

This gives the channel estimates on all subcarriers.

***

#### 3.2.4.2 Two-Dimensional DFT Interpolation

For pilots across multiple OFDM symbols, the most common two-dimensional DFT interpolation is not to directly apply a two-dimensional DFT over the whole time-frequency plane. Instead, it usually uses a separable procedure:

1. **Frequency direction**: for each pilot OFDM symbol, first perform one-dimensional DFT interpolation on the comb pilots to recover the CFR on all subcarriers in that symbol.
2. **Time direction**: between pilot OFDM symbols, perform linear, second-order, or Wiener time interpolation for each subcarrier.

In other words, a common implementation is:

$$
\text{DFT (frequency direction)} + \text{time-domain interpolation}
$$

This is consistent with the logic of linear interpolation, second-order interpolation, and Wiener interpolation above, and it is also the most common engineering implementation.

More generally, if the channel has a stronger low-dimensional structure in the delay-time or delay-Doppler domain, one can further construct two-dimensional DFT/BEM-type methods, but that is beyond the scope of this introductory interpolation note.

***

### 3.2.5 Conditions for DFT Interpolation

The key condition for DFT interpolation is:

$$
\boxed{
L_h \le N_p
}
$$

In practical implementations, a more conservative engineering condition is often written as:

$$
\boxed{
N_{cp} \le N_p
}
$$

If the effective CIR length exceeds $N_p$, the $N_p$-point IDFT will produce time-domain aliasing, and the original CIR cannot be uniquely recovered. In this case, the performance of DFT interpolation degrades significantly.

# 4 References

1. [https://pmc.ncbi.nlm.nih.gov/articles/PMC8309705/](https://pmc.ncbi.nlm.nih.gov/articles/PMC8309705/ "https://pmc.ncbi.nlm.nih.gov/articles/PMC8309705/")
2. LS / MMSE derivation reference: [https://www.cnblogs.com/hjd21/p/16634313.html](https://www.cnblogs.com/hjd21/p/16634313.html "https://www.cnblogs.com/hjd21/p/16634313.html")
3. A. Hjorungnes and D. Gesbert, "Complex-Valued Matrix Differentiation: Techniques and Key Results," IEEE Transactions on Signal Processing, vol. 55, no. 6, pp. 2740-2746, June 2007.
