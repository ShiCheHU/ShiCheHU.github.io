# MIMO-OFDM Pilot Design, Channel Estimation, and Interpolation

> This article focuses on pilot design, channel estimation algorithms, and time-frequency-space/angular-selective channel interpolation in MIMO-OFDM systems. See [MIMO OFDM channel modeling.md](MIMO%20OFDM%20channel%20modeling.md) for channel modeling.

## Table of Contents

- [1 Pilot Design](#1-pilot-design)
  - [1.1 System Model](#11-system-model)
  - [1.2 Why MIMO Requires Separable Pilots](#12-why-mimo-requires-separable-pilots)
  - [1.3 Common Pilot Structures](#13-common-pilot-structures)
  - [1.4 Pilot Spacing Design](#14-pilot-spacing-design)
- [2 Channel Estimation Algorithms](#2-channel-estimation-algorithms)
  - [2.1 Joint Observation Model Based on CIR/CFR](#21-joint-observation-model-based-on-circfr)
  - [2.2 ZF/LS Pilot CFR Initial Estimation](#22-zfls-pilot-cfr-initial-estimation)
  - [2.3 MMSE Channel Estimation](#23-mmse-channel-estimation)
  - [2.4 LMMSE Channel Estimation](#24-lmmse-channel-estimation)
  - [2.5 OMP Angle-Delay Sparse Channel Estimation](#25-omp-angle-delay-sparse-channel-estimation)
- [3 Selective-Fading Channel Reconstruction and Interpolation](#3-selective-fading-channel-reconstruction-and-interpolation)
  - [3.1 Model-Based Reconstruction from Pilot REs to Data REs](#31-model-based-reconstruction-from-pilot-res-to-data-res)
  - [3.2 Frequency-Selective Reconstruction via IDFT/DFT](#32-frequency-selective-reconstruction-via-idftdft)
  - [3.3 Symbol-Level Interpolation Based on Time Selectivity](#33-symbol-level-interpolation-based-on-time-selectivity)
  - [3.4 Time-Frequency-Space Wiener Prediction](#34-time-frequency-space-wiener-prediction)
  - [3.5 Angle-Delay Domain Reconstruction and Denoising](#35-angle-delay-domain-reconstruction-and-denoising)
  - [3.6 Simplification Under Block-Fading Channels](#36-simplification-under-block-fading-channels)
- [4 Appendix](#4-appendix)

# 1 Pilot Design

The goal of pilot design is to enable the receiver to distinguish the channels from different transmit antennas within a limited time-frequency resource, while obtaining sufficient samples in the time, frequency, and space dimensions. Pilot design for MIMO-OFDM is more complex than for SISO because multiple transmit antennas' signals superimpose at each receive antenna.

## 1.1 System Model

We first lay out the simulation parameters and basic assumptions used throughout the subsequent derivations.

| Parameter | Meaning | Typical Constraint or Relationship |
|---|---|---|
| $f_c$ | Carrier frequency | Determines wavelength $\lambda=c/f_c$ and maximum Doppler shift |
| $B$ | System bandwidth | Approximately $B=N_{\mathrm{act}}\Delta f$ |
| $\Delta f$ | Subcarrier spacing | Effective OFDM symbol duration $T_u=1/\Delta f$ |
| $N_{\mathrm{FFT}}$ | FFT size | Denoted simply as $N$ throughout this article |
| $N_{\mathrm{act}}$ | Number of active subcarriers | $N_{\mathrm{act}}\le N_{\mathrm{FFT}}$ |
| $T_{\mathrm{cp}}$ | Cyclic prefix length | Should cover the dominant multipath delay spread |
| $T_{\mathrm{sym}}$ | OFDM symbol duration (including CP) | $T_{\mathrm{sym}}=T_u+T_{\mathrm{cp}}$ |
| $N_t,N_r$ | Number of transmit and receive antennas | Determines MIMO channel matrix dimensions |
| $N_{\mathrm{sy}}$ | Consecutive OFDM symbols per training block | Channel is approximately constant within a block |
| $L_h$ | Number of discrete delay taps | Typically determined by maximum delay and sampling rate |
| $\tau_{\max}$ | Maximum multipath delay | Affects frequency-domain pilot spacing |
| $v$ | Mobile velocity | Together with $f_c$ determines Doppler shift |
| $f_D$ | Maximum Doppler shift | $f_D=vf_c/c$ |
| $\sigma_w^2$ | Noise power | Determines the SNR of received pilot signals |

This article assumes the following conditions hold:

- The CP is sufficiently long so that OFDM inter-subcarrier interference is negligible;
- The block-fading approximation holds within each training block, i.e., $N_{\mathrm{sy}}T_{\mathrm{sym}}$ is shorter than the channel coherence time;
- The wideband channel is described by a finite number of discrete delay taps, and different subcarriers are related to the same set of discrete delay taps via DFT;
- If statistical estimation is used, the PDP, Doppler spectrum, and spatial correlation matrices can be obtained from channel models or measurements;
- If angle-delay sparse estimation is used, the array response and angle dictionaries are consistent with the channel modeling document.

Consider $N_t$ transmit antennas, $N_r$ receive antennas, and $N=N_{\mathrm{FFT}}$ OFDM subcarriers. In practice, estimation is typically performed per training block (or coherence block): it is assumed that the channel is approximately constant over $N_{\mathrm{sy}}$ consecutive OFDM symbols. Let $b$ denote the block index and $q=0,1,\dots,N_{\mathrm{sy}}-1$ denote the symbol index within the block. The frequency-domain model on the $k$-th subcarrier is:

$$
\mathbf y_b[q,k]
=
\mathbf H_b[k]\mathbf x_b[q,k]+\mathbf w_b[q,k],
$$

where:

- $\mathbf y_b[q,k]\in\mathbb C^{N_r\times 1}$ is the received vector;
- $\mathbf x_b[q,k]\in\mathbb C^{N_t\times 1}$ is the transmitted vector;
- $\mathbf H_b[k]\in\mathbb C^{N_r\times N_t}$ is the MIMO frequency-domain channel shared across the $b$-th training block;
- $\mathbf w_b[q,k]\sim\mathcal{CN}(\mathbf 0,\sigma_w^2\mathbf I)$ is the noise.

The block-shared channel assumption is:

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

This does not mean the channel varies independently across OFDM symbols. Rather, it means the channel remains coherent over the interval $N_{\mathrm{sy}}T_{\mathrm{sym}}$, so pilots on multiple consecutive symbols can jointly estimate the same wideband MIMO channel.

When considering only a particular set of pilot resources, with the channel approximately constant over those resources, we can write a local block training model:

$$
\boxed{
\mathbf Y=\mathbf H\mathbf X+\mathbf W
}
$$

where:

- $\mathbf Y\in\mathbb C^{N_r\times N_p}$;
- $\mathbf H\in\mathbb C^{N_r\times N_t}$;
- $\mathbf X\in\mathbb C^{N_t\times N_p}$;
- $N_p$ is the number of pilot received signals.

This model is a special case for narrowband or local block initial estimation. Later channel estimation algorithms will adopt wideband MIMO-OFDM modeling, treating pilots on different symbols and different subcarriers within the same training block as observations of the same set of discrete delay taps or angle-delay paths, rather than regarding each symbol and each subcarrier as having an independently unknown channel.

## 1.2 Why MIMO Requires Separable Pilots

In SISO-OFDM, the initial CFR estimate on a pilot RE can be expressed as:

$$
Y_p=X_pH_p+W_p.
$$

As long as $X_p\ne0$, a simple division yields the LS estimate:

$$
\hat H_p=\frac{Y_p}{X_p}.
$$

In MIMO, however, the received signal on a given pilot RE is a superposition from multiple transmit antennas or ports:

$$
\mathbf y_p
=
\sum_{t=1}^{N_t}\mathbf h_t x_t+\mathbf w_p,
$$

where $\mathbf h_t$ is the $t$-th column of the CFR matrix $\mathbf H_b[k_p]$ at that pilot subcarrier. If multiple transmit antennas use indistinguishable pilots, the receiver cannot tell which antenna or port contributed which received energy.

Therefore, at the pilot CFR initial estimation stage, the local pilot matrix must at least have full row rank:

$$
\boxed{
\operatorname{rank}(\mathbf X)=N_t
}
$$

This requires the number of local pilot received signals to satisfy:

$$
N_p\ge N_t.
$$

A more common and robust design is orthogonal training:

$$
\boxed{
\mathbf X\mathbf X^H=E_p\mathbf I_{N_t}
}
$$

Under this design, pilots from different transmit antennas do not interfere with each other, and ZF/LS can first obtain the CFR initial estimate on the pilot REs.

For wideband channel estimation, merely making the pilots of different transmit antennas orthogonal is not enough. Orthogonality only guarantees that the receiver can first separate the pilot CFRs of different transmit ports; we must also place a sufficient number of pilot subcarriers in the frequency domain for each transmit port in order to recover the finite-delay channel from these CFR initial estimates.

For a fixed receive antenna $r$ and transmit port $t$, let the set of pilot subcarriers occupied by that port within one training block be:

$$
\mathcal P_t=\{k_1,k_2,\dots,k_{N_{p,t}}\}.
$$

On these pilot subcarriers, the CFR vector obtained by ZF/LS initial estimation is defined as:

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

where:

- $\hat{\mathbf h}_{p}\in\mathbb C^{N_{p,t}\times 1}$ is the CFR initial estimate at the pilot positions;
- $\mathbf h_{\tau}\in\mathbb C^{L_h\times 1}$ is the discrete delay tap vector for this antenna pair;
- $\mathbf F_{\mathcal P,L_h}\in\mathbb C^{N_{p,t}\times L_h}$ is a partial DFT matrix extracted from the pilot subcarrier indices, with its $i$-th row, $n$-th column element being $e^{-j2\pi k_i n/N}$;
- $\mathbf e_p$ is the pilot CFR initial estimation error.

Without introducing a sparsity prior, if we only perform ordinary LS/DFT reconstruction with finite discrete delay taps, we need at least:

$$
\boxed{
N_{p,t}\ge L_h
}
$$

and the partial DFT matrix must have full column rank:

$$
\boxed{
\operatorname{rank}(\mathbf F_{\mathcal P,L_h})=L_h
}
$$

Here $L_h$ is the effective number of discrete delay taps, determined by the maximum delay spread and sampling rate — not the number of physical paths. Multiple physical paths may fall into the same discrete delay tap, and off-grid paths may leak into several discrete delay taps. A more precise statement is therefore: the number of frequency-domain pilots per transmit port should be no fewer than the number of discrete delay taps to be estimated; the more pilots, the better the noise averaging and LS stability.

Only when angle-delay sparse estimation is employed is an additional assumption on the effective number of paths (or the sparsity level) required. In that case the number of pilot received signals need not exceed the full $L_h$, but it must be sufficient to recover the target sparse support, and the measurement matrix $\mathbf \Phi$ must not cause severe aliasing among different angle-delay atoms.

## 1.3 Common Pilot Structures

Common pilot structures in MIMO-OFDM can be broadly classified into TDM, FDM, and CDM. They correspond to orthogonalizing pilots in the time, frequency, and code domains, respectively.

**Time-Division Pilots**  
Different transmit antennas send pilots on different OFDM symbols:

| OFDM Symbol | Tx1 | Tx2 | Tx3 |
|---|---|---|---|
| $m_1$ | pilot | 0 | 0 |
| $m_2$ | 0 | pilot | 0 |
| $m_3$ | 0 | 0 | pilot |

Advantage: simplicity. Disadvantages: large time overhead, and estimates for different antennas are not captured at the same time instant under high mobility.

**Frequency-Division Pilots**  
Different transmit antennas occupy disjoint sets of subcarriers:

$$
\mathcal P_i\cap\mathcal P_j=\varnothing,\qquad i\ne j.
$$

Advantage: multiple antennas can be estimated within the same OFDM symbol. Disadvantage: reduced frequency-domain sampling density per antenna.

**Code-Division Pilots**  
Different transmit antennas use orthogonal sequences on the same set of resources:

$$
\sum_{p\in\mathcal P}
x_i[p]x_j^*[p]
=0,\qquad i\ne j.
$$

Advantage: higher resource utilization. Disadvantage: more sensitive to frequency offsets, synchronization errors, and strong frequency selectivity.

## 1.4 Pilot Spacing Design

Pilot spacing is determined by the channel's time selectivity and frequency selectivity.

**Frequency-Domain Spacing**  
Let $S_f$ denote the number of subcarriers between adjacent pilots in frequency. The actual pilot frequency spacing is:

$$
\Delta f_{\mathrm{pilot}}
=
S_f\Delta f.
$$

Here $\Delta f$ is the subcarrier spacing and $\Delta f_{\mathrm{pilot}}$ has units of Hz. To track a frequency-selective channel, the pilot frequency spacing should be on the order of, or smaller than, the coherence bandwidth:

$$
\boxed{
\Delta f_{\mathrm{pilot}}
=
S_f\Delta f
\lesssim
\frac{1}{\tau_{\max}}
}
$$

where $\tau_{\max}$ is the maximum multipath delay.

**Time-Domain Spacing**  
A training block contains $N_{\mathrm{sy}}$ consecutive OFDM symbols. If adjacent pilot training blocks are spaced $S_b$ training blocks apart, the actual pilot time spacing is:

$$
\Delta t_{\mathrm{pilot}}
=
S_bN_{\mathrm{sy}}T_{\mathrm{sym}}.
$$

Here $T_{\mathrm{sym}}$ is the OFDM symbol duration including CP. To track a time-selective channel, the pilot time spacing should be on the order of, or smaller than, the coherence time:

$$
\boxed{
\Delta t_{\mathrm{pilot}}
=
S_bN_{\mathrm{sy}}T_{\mathrm{sym}}
\lesssim
\frac{1}{2f_D}
}
$$

where $f_D$ is the maximum Doppler shift.

**Influence of the Angle Domain**  
Angle-domain selective fading does not directly determine $\Delta f_{\mathrm{pilot}}$ or $\Delta t_{\mathrm{pilot}}$. Frequency-domain spacing is primarily determined by delay spread, and time-domain spacing by Doppler spread. The angle domain affects the pilot measurement resources in the spatial dimension — that is, how many mutually linearly independent transmit ports, spatial layers, precoding directions, or receive combining directions are needed to estimate and exploit the spatial degrees of freedom of the MIMO channel.

When the angular spread is large, the AoA/AoD distribution is more dispersed, and the channel spatial rank and angle-domain effective support are typically larger. In this case, if the system wishes to exploit these spatial degrees of freedom, it must configure pilot ports or precoding measurement resources of sufficient dimensionality. When the angular spread is small, the channel may be low-rank or highly correlated, and the marginal benefit of increasing the number of ports is limited.

Thus, one cannot simply assert that "a larger angular spread necessarily requires more pilot ports." The number of ports is determined by the CSI dimensionality the system needs to estimate, the antenna port configuration, and the number of spatial layers; angle-domain selectivity determines whether these spatial dimensions are effective, whether they are easy to separate, and how many pilot received signals are needed for angle-delay sparse estimation.

MIMO does not alter the fundamental principles of time and frequency sampling, but it does add pilot orthogonality and measurement requirements in the spatial dimension. If each transmit antenna or port requires independent pilots, the pilot overhead generally increases with the dimensionality of the ports to be estimated.

# 2 Channel Estimation Algorithms

This section follows the processing chain of a practical wideband MIMO-OFDM receiver: first estimate CFR initial values on pilot REs, then use finite-delay, PDP/Doppler/spatial-correlation, or angle-delay sparse structure to perform wideband filtering, prediction, and reconstruction on these CFR initial estimates. The narrowband model is a degenerate special case; in wideband processing, the CFRs on different subcarriers are jointly determined by the same set of discrete delay taps or the same set of physical paths.

## 2.1 Joint Observation Model Based on CIR/CFR

The discrete MIMO CIR is written as:

$$
\mathbf H_b[n]\in\mathbb C^{N_r\times N_t},
\qquad n=0,1,\dots,L_h-1.
$$

The CFR on the $k$-th subcarrier of the $b$-th training block is obtained by DFT along the discrete delay dimension:

$$
\boxed{
\mathbf H_b[k]
=
\sum_{n=0}^{L_h-1}
\mathbf H_b[n]e^{-j2\pi kn/N}
}
$$

Thus, within the same training block, different symbols share the same channel, and different subcarriers are not independent unknowns — they are jointly determined by the same set $\{\mathbf H_b[n]\}$. For a particular pilot RE $(q_i,k_i)$ within the block, the received signal model is:

$$
\mathbf y_i
=
\mathbf H_b[k_i]\mathbf x_i+\mathbf w_i.
$$

where:

- $i$ is the index of the $i$-th pilot RE within the training block, corresponding to position $(q_i,k_i)$;
- $\mathbf y_i\in\mathbb C^{N_r\times 1}$ is the received frequency-domain signal on this pilot RE;
- $\mathbf x_i\in\mathbb C^{N_t\times 1}$ is the transmitted pilot vector on this pilot RE;
- $\mathbf H_b[k_i]\in\mathbb C^{N_r\times N_t}$ is the MIMO CFR on the $k_i$-th subcarrier of the $b$-th training block;
- $\mathbf w_i\in\mathbb C^{N_r\times 1}$ is the noise vector on this pilot RE.

Define the discrete delay tap vector:

$$
\mathbf h_b
=
\operatorname{vec}
\left(
\mathbf H_b[0],\mathbf H_b[1],\dots,\mathbf H_b[L_h-1]
\right)
\in\mathbb C^{N_rN_tL_h\times 1}.
$$

Let:

$$
\mathbf a(k)
=
\left[
1,e^{-j2\pi k/N},\dots,e^{-j2\pi k(L_h-1)/N}
\right]^T.
$$

Then:

$$
\operatorname{vec}(\mathbf H_b[k])
=
(\mathbf a^T(k)\otimes\mathbf I_{N_rN_t})\mathbf h_b.
$$

Substituting into $\mathbf y_i=(\mathbf x_i^T\otimes\mathbf I_{N_r})\operatorname{vec}(\mathbf H_b[k_i])+\mathbf w_i$, the single-pilot received signal becomes:

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

Here $\mathbf A_i\in\mathbb C^{N_r\times N_rN_tL_h}$ is the observation matrix of the $i$-th pilot RE onto the wideband discrete delay tap vector $\mathbf h_b$.

Suppose the $b$-th training block contains $N_{\mathrm{obs}}$ pilot REs in total. Here $N_{\mathrm{obs}}$ is the number of time-frequency pilot positions — that is, the count of all $(q_i,k_i)$; it does not include the receive-antenna dimension. Each pilot RE yields one $\mathbf y_i\in\mathbb C^{N_r\times 1}$. Stack these received vectors vertically:

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

Correspondingly, stack the $\mathbf A_i\in\mathbb C^{N_r\times N_rN_tL_h}$ vertically:

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

The noise vectors are stacked likewise:

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

All pilot received signals can then be written uniformly as:

$$
\boxed{
\mathbf y=\mathbf A\mathbf h+\mathbf w
}
$$

In the antenna-delay domain model of this section, $\mathbf h=\mathbf h_b\in\mathbb C^{N_rN_tL_h\times 1}$. More generally, $\mathbf h$ is the channel parameter used for wideband modeling, which may be:

- An antenna-discrete-delay domain (tap domain) vector;
- A joint time-frequency-space channel vector;
- An angle-delay domain sparse vector.

In practice, ZF/LS often first provides CFR initial estimates on pilot subcarriers; LMMSE/Wiener, IDFT/DFT reconstruction, and OMP then feed these pilot received signals or CFR initial estimates into the wideband model to generate $\hat{\mathbf H}[q,k]$ on the target REs. Under the block-fading simplification, $\hat{\mathbf H}[q,k]$ can degenerate to the block-shared $\hat{\mathbf H}_b[k]$.

## 2.2 ZF/LS Pilot CFR Initial Estimation

In a practical OFDM receiver, the most common use of ZF/LS is to first estimate the CFR on pilot REs. For a given pilot subcarrier $k_p$, if pilots of the $t$-th transmit port can be separated via TDM, FDM, or CDM design, we have:

$$
\mathbf y_b[q,k_p]
=
\mathbf H_b[k_p]\mathbf x_b[q,k_p]+\mathbf w_b[q,k_p].
$$

For a fixed training block $b$ and fixed pilot subcarrier $k_p$, the $N_p$ frequency-domain received signals produced by orthogonal pilots are expressed in matrix form:

$$
\mathbf Y_b[k_p]
=
\mathbf H_b[k_p]\mathbf X_b[k_p]
+\mathbf W_b[k_p].
$$

where:

- $\mathbf Y_b[k_p]\in\mathbb C^{N_r\times N_p}$ is the frequency-domain received signal matrix collected on the $k_p$-th pilot subcarrier of the $b$-th training block; each column corresponds to one pilot receive observation;
- $\mathbf H_b[k_p]\in\mathbb C^{N_r\times N_t}$ is the MIMO CFR on the $k_p$-th subcarrier of the $b$-th training block;
- $\mathbf X_b[k_p]\in\mathbb C^{N_t\times N_p}$ is the transmitted pilot matrix on the $k_p$-th subcarrier of the $b$-th training block; each column is one pilot transmit vector;
- $\mathbf W_b[k_p]\in\mathbb C^{N_r\times N_p}$ is the corresponding noise matrix;
- $N_p$ is the number of pilot observations on this subcarrier used to separate the $N_t$ transmit ports, and must satisfy $N_p\ge N_t$.

If $\mathbf X_b[k_p]$ has full row rank, i.e., $\operatorname{rank}(\mathbf X_b[k_p]) = N_t$, the ZF/LS CFR initial estimate at the pilot position is:

$$
\boxed{
\hat{\mathbf H}_{\mathrm{LS},b}[k_p]
=
\mathbf Y_b[k_p]\mathbf X_b[k_p]^H
\left(\mathbf X_b[k_p]\mathbf X_b[k_p]^H\right)^{-1}
}
$$

If the pilots are orthogonal, $\mathbf X_b[k_p]\mathbf X_b[k_p]^H=E_p\mathbf I$, then:

$$
\boxed{
\hat{\mathbf H}_{\mathrm{LS},b}[k_p]
=
\frac{1}{E_p}\mathbf Y_b[k_p]\mathbf X_b[k_p]^H
}
$$

This step estimates the CFR on pilot subcarriers — it does not estimate the number of physical paths, nor does it necessarily estimate the CIR directly. It only completes the "initial estimation at pilot positions." Exploiting the wideband structure belongs to the subsequent reconstruction/interpolation steps: one can feed these pilot CFR estimates into the IDFT/DFT finite-delay reconstruction, Wiener/LMMSE time-frequency-space prediction, or OMP angle-delay reconstruction described in Chapter 3.

Thus, ZF/LS itself does not require any assumption on the number of physical paths. If IDFT/DFT finite-delay reconstruction is used subsequently, what must be assumed is the maximum effective discrete delay length $L_h$, not the number of physical paths; the assumption on the number of physical paths or the sparsity level appears primarily in sparse estimation algorithms such as OMP.

**Algorithm Pseudocode**

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

The advantage of ZF/LS is its simplicity and alignment with practical pilot processing. The disadvantage is that the pilot CFR initial estimate itself does not exploit PDP, Doppler, or spatial correlation; at low SNR or with sparse pilots, subsequent wideband reconstruction or time-frequency-space interpolation is needed.

## 2.3 MMSE Channel Estimation

The goal of MMSE estimation is to find the posterior mean under the full prior distribution:

$$
\boxed{
\hat{\mathbf h}_{\mathrm{MMSE}}
=
\mathbb E[\mathbf h\mid \mathbf y]
}
$$

Here $\mathbf h$ is no longer the channel matrix vector of a single subcarrier, but a wideband MIMO channel parameter. Its prior can be derived from the models in the channel modeling document:

- i.i.d. Rayleigh/Rician discrete delay tap model: the prior is determined by the PDP and K-factor;
- Kronecker spatial correlation model: the prior is determined by $\mathbf R_r,\mathbf R_t$ and the PDP;
- Geometric angle-delay model: the prior is determined by the number of paths, AoA/AoD, delays, and path-gain distributions;
- CDL/TDL models: the prior is determined by standardized clusters, rays, delays, angles, and Doppler parameters.

Using the stacked model from 2.1, the variable dimensions are:

- $\mathbf y\in\mathbb C^{N_rN_{\mathrm{obs}}\times 1}$: the stacked vector of all pilot received signals within the training block;
- $\mathbf A\in\mathbb C^{N_rN_{\mathrm{obs}}\times N_h}$: the joint observation matrix;
- $\mathbf h\in\mathbb C^{N_h\times 1}$: the wideband channel parameter vector to be estimated;
- $\mathbf w\in\mathbb C^{N_rN_{\mathrm{obs}}\times 1}$: the noise vector;
- $\mathbf R_w\in\mathbb C^{N_rN_{\mathrm{obs}}\times N_rN_{\mathrm{obs}}}$: the noise covariance matrix.

If $\mathbf h$ is expressed in the antenna-discrete-delay domain (tap domain), then $N_h=N_rN_tL_h$.

MMSE is derived from Bayes' formula:

$$
p(\mathbf h\mid\mathbf y)
=
\frac{p(\mathbf y\mid\mathbf h)p(\mathbf h)}{p(\mathbf y)}.
$$

where:

$$
p(\mathbf y\mid\mathbf h)
\propto
\exp\left(
-\|\mathbf y-\mathbf A\mathbf h\|_{\mathbf R_w^{-1}}^2
\right).
$$

Hence:

$$
\boxed{
\hat{\mathbf h}_{\mathrm{MMSE}}
=
\int \mathbf h\,p(\mathbf h\mid\mathbf y)\,d\mathbf h
}
$$

The prior of a general geometric channel or CDL channel is not necessarily a simple Gaussian distribution. In such cases MMSE may not have a closed-form solution and requires approximation via grid integration, particle filtering, EM, message passing, or neural networks. If $\mathbf h$ and $\mathbf w$ are jointly Gaussian, MMSE reduces to LMMSE.

**Algorithm Pseudocode**

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

MMSE is the theoretically optimal mean-squared-error estimator; whether it is realizable depends on whether the channel model prior is sufficiently simple.

## 2.4 LMMSE Channel Estimation

LMMSE uses only second-order statistics, restricting the estimator to a linear form:

$$
\hat{\mathbf h}=\mathbf B\mathbf y.
$$

Using the notation from 2.3, $\mathbf y\in\mathbb C^{N_rN_{\mathrm{obs}}\times 1}$, $\mathbf h\in\mathbb C^{N_h\times 1}$. Therefore:

- $\mathbf B\in\mathbb C^{N_h\times N_rN_{\mathrm{obs}}}$ is the linear estimation matrix;
- $\mathbf R_h\in\mathbb C^{N_h\times N_h}$ is the wideband channel parameter covariance;
- $\mathbf R_w\in\mathbb C^{N_rN_{\mathrm{obs}}\times N_rN_{\mathrm{obs}}}$ is the noise covariance;
- $\mathbf R_{hy}\in\mathbb C^{N_h\times N_rN_{\mathrm{obs}}}$;
- $\mathbf R_{yy}\in\mathbb C^{N_rN_{\mathrm{obs}}\times N_rN_{\mathrm{obs}}}$.

For the joint observation model $\mathbf y=\mathbf A\mathbf h+\mathbf w$, the orthogonality principle gives:

$$
\mathbb E[
(\mathbf h-\mathbf B\mathbf y)\mathbf y^H
]
=\mathbf 0.
$$

Therefore:

$$
\mathbf B
=
\mathbf R_{hy}\mathbf R_{yy}^{-1}.
$$

If $\mathbf h$ and $\mathbf w$ are uncorrelated:

$$
\mathbf R_{hy}
=
\mathbf R_h\mathbf A^H,
\qquad
\mathbf R_{yy}
=
\mathbf A\mathbf R_h\mathbf A^H+\mathbf R_w.
$$

This yields the joint LMMSE:

$$
\boxed{
\hat{\mathbf h}_{\mathrm{LMMSE}}
=
\mathbf R_h\mathbf A^H
(\mathbf A\mathbf R_h\mathbf A^H+\mathbf R_w)^{-1}
\mathbf y
}
$$

The key point is that $\mathbf R_h$ should come from MIMO-OFDM channel modeling, not from assuming the subcarriers are independent. Under the WSSUS + Kronecker spatial correlation model, the discrete delay domain (tap domain) covariance can be written as:

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

where:

- $\mathbf R_{\mathrm{time}}$ describes the time correlation between training blocks; it can be omitted if only a single training block is estimated;
- $\mathbf R_{\mathrm{delay}}\in\mathbb C^{L_h\times L_h}$ is typically determined by the PDP; if different discrete delay taps are independent, then $\mathbf R_{\mathrm{delay}}=\operatorname{diag}(P_0,\dots,P_{L_h-1})$;
- $\mathbf R_t\in\mathbb C^{N_t\times N_t}$ is the transmit-side spatial correlation matrix;
- $\mathbf R_r\in\mathbb C^{N_r\times N_r}$ is the receive-side spatial correlation matrix.

When estimating only a single training block with the antenna-discrete-delay domain (tap domain) vector, a common simplified form is:

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

These correlation matrices can be obtained from model parameters or measurement estimates:

**Time Correlation Matrix**  
If the maximum Doppler shift $f_D$ is known, the Jakes/Clarke model can be used to generate the time correlation between training blocks. The correlation coefficient between the $b_1$-th and $b_2$-th training blocks can be written as:

$$
[\mathbf R_{\mathrm{time}}]_{b_1,b_2}
=
J_0\!\left(
2\pi f_DN_{\mathrm{sy}}T_{\mathrm{sym}}|b_1-b_2|
\right).
$$

If estimation is only within a single training block and the channel is assumed constant within the block, $\mathbf R_{\mathrm{time}}$ is typically not constructed explicitly.

**Delay Correlation Matrix**  
$\mathbf R_{\mathrm{delay}}$ is determined by the PDP. If discrete delay taps are mutually independent under the WSSUS assumption, then:

$$
\mathbf R_{\mathrm{delay}}
=
\operatorname{diag}(P_0,P_1,\dots,P_{L_h-1}),
$$

where $P_n=\mathbb E\|\mathbf H_b[n]\|_F^2$ denotes the average power of the $n$-th discrete delay tap, usually normalized such that $\sum_n P_n=1$. The PDP can come from standard models, measurements, or simulation settings.

**Spatial Correlation Matrices**  
$\mathbf R_t$ and $\mathbf R_r$ can be generated from the angular power spectrum, exponential correlation model, sample covariance, or CDL/geometric path parameters. For a ULA, with normalized angle $\Omega=(d/\lambda)\sin\theta$ and angular power spectrum $p(\Omega)$, the entries of the spatial correlation matrix can be written as:

$$
[\mathbf R]_{p,q}
=
\int p(\Omega)e^{-j2\pi(p-q)\Omega}\,d\Omega.
$$

Under a discrete path/ray model:

$$
[\mathbf R]_{p,q}
=
\frac{\sum_{\ell}P_\ell e^{-j2\pi(p-q)\Omega_\ell}}
{\sum_{\ell}P_\ell}.
$$

An exponential correlation approximation is also possible:

$$
[\mathbf R]_{p,q}=\rho^{|p-q|},
\qquad 0\le|\rho|<1.
$$

In measurements or link simulations, the sample covariance can also be estimated from multiple channel snapshots:

$$
\hat{\mathbf R}
=
\frac{1}{N_s}\sum_{s=1}^{N_s}
\mathbf h_s\mathbf h_s^H.
$$

The estimation error covariance is:

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

After obtaining $\hat{\mathbf h}_{\mathrm{LMMSE}}$, one first recovers $\hat{\mathbf H}_b[n]$ and then performs a DFT along the discrete delay dimension to generate $\hat{\mathbf H}_b[k]$.

**Algorithm Pseudocode**

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

LMMSE is the most common "channel-model-based" joint estimation approach in engineering. It explicitly exploits the PDP, spatial correlation, and time correlation, making it more consistent with the wideband MIMO-OFDM channel structure than per-subcarrier LS.

## 2.5 OMP Angle-Delay Sparse Channel Estimation

For the $b$-th training block, the angle-domain representation from the channel modeling document is:

$$
\mathbf H_b^a[k]
=
\mathbf U_r^H\mathbf H_b[k]\mathbf U_t.
$$

Taking an IDFT along the subcarrier dimension yields the angle-delay domain channel:

$$
\mathbf H_b^a[n]
=
\mathbf U_r^H\mathbf H_b[n]\mathbf U_t.
$$

In massive MIMO, mmWave, or scenes with few paths, $\mathbf H_b^a[n]$ is approximately sparse over the three-dimensional grid:

$$
(\text{AoA bin},\text{AoD bin},\text{delay tap})
$$

Let:

$$
\mathbf s
=
\operatorname{vec}\{\mathbf H_b^a[n]\}
$$

be the angle-delay sparse vector. From:

$$
\mathbf H_b[n]
=
\mathbf U_r\mathbf H_b^a[n]\mathbf U_t^H
$$

and the DFT along the discrete delay dimension, we obtain:

$$
\mathbf H_b[k_i]
=
\sum_{n=0}^{L_h-1}
\mathbf U_r\mathbf H_b^a[n]\mathbf U_t^H
e^{-j2\pi k_i n/N}.
$$

For the $i$-th pilot RE:

$$
\mathbf y_i
=
\mathbf H_b[k_i]\mathbf x_i+\mathbf w_i.
$$

Using the identity $\mathbf A\mathbf B\mathbf C=(\mathbf C^T\otimes\mathbf A)\operatorname{vec}(\mathbf B)$, we get:

$$
\mathbf y_i
=
\mathbf \Phi_i\mathbf s+\mathbf w_i.
$$

If we define:

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

then the measurement matrix for the $i$-th pilot RE is:

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

where:

- $\mathbf \Phi_i\in\mathbb C^{N_r\times N_rN_tL_h}$;
- $\mathbf x_i^T\mathbf U_t^*\in\mathbb C^{1\times N_t}$ represents the equivalent weights of the transmit pilot vector on the transmit angle dictionary;
- $\mathbf U_r\in\mathbb C^{N_r\times N_r}$ is the receive angle dictionary;
- Each horizontal block corresponds to one discrete delay tap.

Stacking all pilots:

$$
\boxed{
\mathbf y=\mathbf \Phi\mathbf s+\mathbf w
}
$$

where:

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

Thus $\mathbf \Phi$ simultaneously encodes:

- The transmit pilots $\mathbf x_i$;
- The receive and transmit angle dictionaries $\mathbf U_r,\mathbf U_t$;
- The subcarrier-to-discrete-delay-tap phases $e^{-j2\pi k_in/N}$;
- Any receive combining or other front-end linear processing matrices that may be present.

OMP solves:

$$
\min_{\mathbf s}\|\mathbf s\|_0,
\qquad
\text{s.t.}
\quad
\|\mathbf y-\mathbf \Phi\mathbf s\|_2^2\le \epsilon.
$$

At each iteration, it selects the angle-delay atom most correlated with the residual:

$$
j_i
=
\arg\max_j
\left|
\boldsymbol\phi_j^H\mathbf r^{(i-1)}
\right|.
$$

LS is performed on the selected support $\mathcal S^{(i)}$:

$$
\boxed{
\hat{\mathbf s}_{\mathcal S^{(i)}}
=
(\mathbf \Phi_{\mathcal S^{(i)}}^H
\mathbf \Phi_{\mathcal S^{(i)}})^{-1}
\mathbf \Phi_{\mathcal S^{(i)}}^H\mathbf y
}
$$

Residual update:

$$
\mathbf r^{(i)}
=
\mathbf y
-\mathbf \Phi_{\mathcal S^{(i)}}
\hat{\mathbf s}_{\mathcal S^{(i)}}.
$$

After recovering $\hat{\mathbf s}$, first obtain the angle-delay domain $\hat{\mathbf H}_b^a[n]$, then transform back to the antenna-delay domain:

$$
\hat{\mathbf H}_b[n]
=
\mathbf U_r\hat{\mathbf H}_b^a[n]\mathbf U_t^H,
$$

and finally generate the CFR:

$$
\hat{\mathbf H}_b[k]
=
\sum_n
\hat{\mathbf H}_b[n]e^{-j2\pi kn/N}.
$$

**Algorithm Pseudocode**

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

The core advantage of OMP is that it directly matches the angle-delay sparse channel model; the downside is its reliance on dictionary accuracy. If the true AoA/AoD or delays fall off the grid, off-grid leakage occurs.

# 3 Selective-Fading Channel Reconstruction and Interpolation

This chapter discusses interpolation under time-, frequency-, and angle-selective fading channels. In this regime the channel does not simply vary across subcarriers, nor is it strictly constant within a training block. Instead, it simultaneously exhibits:

- Frequency selectivity: determined by the finite delay spread or PDP;
- Time selectivity: determined by the Doppler spread; the channel varies slowly across different OFDM symbols yet remains correlated;
- Space/angle selectivity: determined by the AoA/AoD angular power spectrum, array response, and spatial correlation.

Hence, "interpolation" should not be understood as fitting curves through mutually independent subcarriers. It should be understood as model-based reconstruction: first obtain CFR initial estimates on pilot REs, then, based on finite delay, time correlation, spatial/angular correlation, or angle-delay sparse structure, extrapolate the information from pilot REs to data REs. The core idea is to recover the shared selective-fading structure from the pilot CFR initial estimates and then generate $\hat{\mathbf H}[q,k]$ on the target time-frequency resources.

## 3.1 Model-Based Reconstruction from Pilot REs to Data REs

Let the set of pilot positions within an interpolation window be $\mathcal P$ and the set of data positions be $\mathcal D$. Each position is specified by an (OFDM symbol, subcarrier) pair:

$$
(q,k)\in\mathcal P\cup\mathcal D.
$$

ZF/LS first yields the CFR initial estimates at the pilot positions:

$$
\{\hat{\mathbf H}_{p}[q,k]:(q,k)\in\mathcal P\}.
$$

Selective-fading channel interpolation/reconstruction further transforms these CFR initial estimates into some model parameters:

$$
\hat{\boldsymbol\theta}
\in
\{\hat{\mathbf H}[q,n],\hat{\mathbf h},\hat{\mathbf H}^a[q,n],\hat{\mathbf s}\}.
$$

Channel interpolation is thus a two-stage mapping:

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

Different methods differ in how much selective-fading channel model information this mapping exploits:

- IDFT/DFT reconstruction: exploits the finite delay spread;
- Wiener/LMMSE prediction: exploits the PDP, Doppler, and spatial correlation;
- Angle-delay / angle-delay-Doppler reconstruction: exploits the structure or sparsity of paths in the angle, delay, and Doppler domains.

## 3.2 Frequency-Selective Reconstruction via IDFT/DFT

IDFT/DFT reconstruction only handles frequency selectivity. Its core steps are: for a given OFDM symbol $q$, first transform the CFR initial estimates on pilot subcarriers into the discrete delay domain to obtain an equivalent CIR; retain only the discrete delay taps within the effective length $L_h$; then reconstruct the full-band CFR for that symbol via DFT.

If a full-band coarse CFR vector is available:

$$
\hat{\mathbf h}_{f,\mathrm{raw}}
=
\begin{bmatrix}
\hat H[q,0] & \hat H[q,1] & \cdots & \hat H[q,N-1]
\end{bmatrix}^T,
$$

first apply IDFT:

$$
\hat{\mathbf h}_{\tau,\mathrm{raw}}
=
\mathbf F_N^H\hat{\mathbf h}_{f,\mathrm{raw}}.
$$

Retain only the first $L_h$ effective discrete delay taps:

$$
\hat h_{\tau}[n]
=
\begin{cases}
\hat h_{\tau,\mathrm{raw}}[n], & 0\le n<L_h,\\
0, & L_h\le n<N.
\end{cases}
$$

Then apply DFT to obtain the equivalent CFR:

$$
\boxed{
\hat{\mathbf h}_{f}
=
\mathbf F_N\hat{\mathbf h}_{\tau}
}
$$

If the equivalent discrete delay taps for the $q$-th OFDM symbol have already been obtained from the pilot CFR initial estimates:

$$
\hat{\mathbf H}[q,n],
\qquad n=0,\dots,L_h-1,
$$

then the equivalent CFR on any subcarrier can be obtained directly via DFT:

$$
\boxed{
\hat{\mathbf H}[q,k]
=
\sum_{n=0}^{L_h-1}
\hat{\mathbf H}[q,n]e^{-j2\pi kn/N}
}
$$

This formula automatically ensures that different subcarriers within the same OFDM symbol share the same set of discrete delay taps. It is the most basic frequency-domain reconstruction method in wideband OFDM. If the channel varies with time, $\hat{\mathbf H}[q,n]$ must be obtained or predicted separately for each target symbol $q$; if the frequency-domain pilots on a particular symbol are insufficient, the time-correlation models of sections 3.3 or 3.4 should be combined for joint recovery.

If the system first performs per-RE LS on pilot subcarriers to obtain coarse estimates $\hat{\mathbf H}_{p}[q,k]$, then for each symbol $q$ with sufficient frequency-domain pilots, IDFT/DFT reconstruction can be performed on each antenna pair or on the joint vector. Define the set of pilot subcarriers on that symbol as:

$$
\mathcal P_q=\{k:(q,k)\in\mathcal P\}.
$$

If pilots only cover a subset of subcarriers, an equivalent approach is to fit a length-$L_h$ CIR using a partial DFT:

$$
\hat{\mathbf h}_{\tau}
=
(\mathbf F_{\mathcal P_q,L_h}^H\mathbf F_{\mathcal P_q,L_h})^{-1}
\mathbf F_{\mathcal P_q,L_h}^H
\hat{\mathbf h}_{p}.
$$

Then generate the equivalent CFR on all subcarriers via DFT:

$$
\hat{\mathbf h}_{f}
=
\mathbf F_{N,L_h}\hat{\mathbf h}_{\tau}.
$$

**Algorithm Pseudocode**

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

## 3.3 Symbol-Level Interpolation Based on Time Selectivity

Time selectivity arises from Doppler. For doubly selective fading channels, the channel is typically not strictly identical across different OFDM symbols within the same observation window; rather, it varies slowly and exhibits time correlation. Hence, time-domain interpolation should be performed between pilot symbols and data symbols, and the preferred interpolation targets are the discrete delay taps, angle-delay coefficients, or the joint channel vector — not each subcarrier independently.

If the discrete delay tap estimates on the $q_1$-th and $q_2$-th pilot symbols are:

$$
\hat{\mathbf H}[q_1,n],
\qquad
\hat{\mathbf H}[q_2,n].
$$

Then for a data symbol $q$ lying between the two pilot symbols, a low-complexity approach is to perform symbol-level linear interpolation on each discrete delay tap:

$$
\boxed{
\hat{\mathbf H}[q,n]
=
\frac{q_2-q}{q_2-q_1}\hat{\mathbf H}[q_1,n]
+
\frac{q-q_1}{q_2-q_1}\hat{\mathbf H}[q_2,n]
}
$$

Then perform a DFT along the discrete delay dimension to generate $\hat{\mathbf H}[q,k]$. If Doppler statistics are used, the Jakes/Clarke correlation function can be employed. When the separation between two OFDM symbols is $\Delta q$, the approximate time interval is $\Delta qT_{\mathrm{sym}}$:

$$
R_{\mathrm{time}}(\Delta q)
=
J_0(2\pi f_D T_{\mathrm{sym}}\Delta q).
$$

The pilot symbol spacing must still satisfy the coherence time constraint:

$$
S_tT_{\mathrm{sym}}\lesssim \frac{1}{2f_D}.
$$

## 3.4 Time-Frequency-Space Wiener Prediction

Wiener prediction is the predictive form of LMMSE on data REs, and it is the most unified interpolation formulation for time-, frequency-, and space/angle-selective fading channels. Let $\mathbf h_p$ denote the joint channel vector at pilot positions and $\mathbf h_d$ denote the joint channel vector at data positions. The pilot-position estimate is:

$$
\hat{\mathbf h}_p=\mathbf h_p+\mathbf e_p.
$$

The linear predictor:

$$
\hat{\mathbf h}_d=\mathbf B\hat{\mathbf h}_p.
$$

Minimizing:

$$
\mathbb E\|\mathbf h_d-\mathbf B\hat{\mathbf h}_p\|_2^2.
$$

By the orthogonality principle:

$$
\boxed{
\hat{\mathbf h}_d
=
\mathbf R_{d,p}
(\mathbf R_{p,p}+\mathbf R_{e,p})^{-1}
\hat{\mathbf h}_p
}
$$

where:

- $\mathbf R_{p,p}$: the time, frequency, and space/angle covariance among pilot REs;
- $\mathbf R_{d,p}$: the time, frequency, and space/angle cross-covariance between data REs and pilot REs;
- $\mathbf R_{e,p}$: the pilot estimation error covariance.

These covariances come from the channel model. Under the WSSUS + Kronecker model they can be written as:

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

The frequency correlation is determined by the PDP:

$$
\boxed{
R_f(\Delta k)
=
\sum_{n=0}^{L_h-1}
P_n e^{-j2\pi \Delta k n/N}
}
$$

The time correlation can be determined by the Doppler model. Taking the OFDM symbol as the time sampling unit:

$$
R_{\mathrm{time}}(\Delta q)
=
J_0(2\pi f_DT_{\mathrm{sym}}\Delta q).
$$

The spatial correlations $\mathbf R_{\mathrm{tx}},\mathbf R_{\mathrm{rx}}$ can be generated from the angular power spectrum, exponential correlation model, or geometric paths. The influence of angle-selective fading on interpolation manifests in the spatial/angular covariance or the angle dictionaries: if the channel correlation differs across antenna ports or beam directions, Wiener prediction must use $\mathbf R_{\mathrm{tx}},\mathbf R_{\mathrm{rx}}$ or the angle-domain representation to map the information from pilot ports to the channel on the target ports or target beams.

**Algorithm Pseudocode**

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

## 3.5 Angle-Delay Domain Reconstruction and Denoising

If the estimation algorithm outputs angle-delay sparse coefficients $\hat{\mathbf s}$, the channel on data REs is reconstructed directly from the physical structure. For time-selective channels, path gains or angle-delay coefficients should be allowed to vary slowly across OFDM symbols $q$. First recover:

$$
\hat{\mathbf H}^a[q,n],
$$

then transform back to the antenna-delay domain:

$$
\hat{\mathbf H}[q,n]
=
\mathbf U_r\hat{\mathbf H}^a[q,n]\mathbf U_t^H.
$$

Finally, perform a DFT along the discrete delay dimension:

$$
\boxed{
\hat{\mathbf H}[q,k]
=
\sum_{n=0}^{L_h-1}
\mathbf U_r\hat{\mathbf H}^a[q,n]\mathbf U_t^H
e^{-j2\pi kn/N}
}
$$

This class of methods is essentially not ordinary interpolation, but "path-parameter reconstruction": as long as the angles, delays, and time-varying path gains are estimated accurately, a consistent MIMO channel can be generated on the target OFDM symbols and all subcarriers. If an angle-delay-Doppler dictionary is further adopted, the sparse coefficients can directly describe the path variations across different Doppler bins.

**Algorithm Pseudocode**

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

Angle-delay reconstruction is well suited to massive MIMO and mmWave systems. If the angle or delay grids are too coarse, oversampled dictionaries, off-grid corrections, or continuous parameter estimation are needed.

## 3.6 Simplification Under Block-Fading Channels

If the system satisfies the block-fading approximation — i.e., within $N_{\mathrm{sy}}$ consecutive OFDM symbols:

$$
\mathbf H[0,k]\approx \mathbf H[1,k]\approx \cdots \approx \mathbf H[N_{\mathrm{sy}}-1,k],
$$

then time selectivity can be neglected within that block. In this case, the selective-fading interpolation of Chapter 3 simplifies to:

- Time dimension: no symbol-level time-domain interpolation is needed within the block; CFR initial estimates on the same subcarrier but different pilot symbols can be combined to improve SNR;
- Frequency dimension: IDFT/DFT, LMMSE, or Wiener frequency-domain reconstruction is still required based on the finite delay spread or PDP;
- Space/angle dimension: multi-antenna joint estimation can still exploit spatial correlation, angular power spectra, or angle-delay sparsity.

If prediction is needed between adjacent training blocks, each training block can be treated as one time-sampling point, and low-rate time interpolation can be performed on the block index. In this case, the time interval in the Jakes correlation function should use the training-block center interval, e.g., $\Delta t=\Delta bN_{\mathrm{sy}}T_{\mathrm{sym}}$, rather than the symbol spacing $\Delta qT_{\mathrm{sym}}$.

# 4 Appendix

## A.1 Notation

| Symbol | Meaning |
|---|---|
| $N_t$ | Number of transmit antennas |
| $N_r$ | Number of receive antennas |
| $f_c$ | Carrier frequency |
| $B$ | System bandwidth |
| $\Delta f$ | Subcarrier spacing |
| $N_{\mathrm{FFT}}$ | FFT size, often abbreviated as $N$ in the text |
| $N_{\mathrm{act}}$ | Number of active subcarriers |
| $T_{\mathrm{cp}}$ | Cyclic prefix length |
| $T_u$ | Effective OFDM symbol duration (without CP), $T_u=1/\Delta f$ |
| $T_{\mathrm{sym}}$ | OFDM symbol duration including CP |
| $N_p$ | Number of pilot received signals |
| $N_{\mathrm{obs}}$ | Number of pilot REs within one training block, excluding the receive-antenna dimension |
| $N_h$ | Dimension of the wideband channel parameter vector $\mathbf h$ |
| $\mathbf X$ | Local pilot matrix |
| $\mathbf Y$ | Local received pilot matrix |
| $N_{\mathrm{sy}}$ | Number of consecutive OFDM symbols within one observation window or block |
| $b$ | Training block or coherence block index, mainly used in the block-fading simplification |
| $q$ | OFDM symbol index |
| $\mathbf H[q,k]$ | MIMO CFR on the $q$-th OFDM symbol, $k$-th subcarrier |
| $\mathbf H[q,n]$ | MIMO CIR on the $q$-th OFDM symbol, $n$-th discrete delay tap |
| $\mathbf H_b[k]$ | Under block-fading simplification: MIMO CFR on the $b$-th training block, $k$-th subcarrier |
| $\mathbf H_b[n]$ | Under block-fading simplification: MIMO CIR on the $b$-th training block, $n$-th discrete delay tap |
| $\mathbf A$ | Wideband joint observation matrix |
| $\mathbf h$ | Vectorized wideband channel parameter, can represent discrete-delay-domain (tap domain) channel or joint time-frequency-space channel |
| $\mathbf R_h$ | Channel covariance matrix |
| $\mathbf R_w$ | Noise covariance matrix |
| $\mathbf R_t$ | Transmit-side spatial correlation matrix |
| $\mathbf R_r$ | Receive-side spatial correlation matrix |
| $\mathbf \Phi$ | Sparse measurement matrix |
| $\mathbf s$ | Angle-delay domain sparse vector |
| $\mathbf U_t,\mathbf U_r$ | Transmit- and receive-side DFT angle dictionaries |
| $v$ | Mobile velocity |
| $f_D$ | Maximum Doppler shift |
| $\tau_{\max}$ | Maximum multipath delay |
| $L_h$ | Maximum number of discrete delay taps |
| $\sigma_w^2$ | Noise power |

## A.2 Method Comparison

| Method | Estimation Target | Uses Wideband Model? | Main Advantages | Main Disadvantages | Suitable Scenarios |
|---|---|---|---|---|---|
| ZF/LS | CFR initial estimate on pilot REs | Basic LS does not; can be followed by IDFT/DFT, Wiener/LMMSE, or OMP reconstruction | Simple, close to practical receiver | Initial estimate poor against noise; requires subsequent reconstruction/interpolation with sparse pilots | Pilot initial estimation |
| MMSE | Posterior mean of wideband channel parameters | Yes, via full prior distribution | MSE optimal | Prior and integration usually intractable | Theoretical analysis, Bayesian estimation |
| LMMSE | CIR or joint time-frequency-space channel vector | Yes, via PDP, spatial correlation, and Doppler covariance | Noise-robust, exploits correlations | Covariance estimation and matrix inversion are complex | Practical channel estimation, Wiener interpolation |
| OMP | Angle-delay sparse coefficients | Yes, via angle dictionaries and delay DFT | Low pilot overhead, matches sparse physical paths | Relies on sparsity and grid accuracy | Massive MIMO, mmWave |

## A.3 References

1. D. Tse and P. Viswanath, *Fundamentals of Wireless Communication*, Cambridge University Press, 2005. https://web.stanford.edu/~dntse/wireless_book.html
2. A. M. Sayeed, "Deconstructing Multiantenna Fading Channels," *IEEE Transactions on Signal Processing*, vol. 50, no. 10, pp. 2563-2579, 2002. https://minds.wisconsin.edu/handle/1793/9386
3. J. Jo and I. Sohn, "On the optimality of training signals for MMSE channel estimation in MIMO-OFDM systems," *EURASIP Journal on Wireless Communications and Networking*, 2015. https://link.springer.com/article/10.1186/s13638-015-0345-y
4. J. P. Nair and R. V. Raja Kumar, "Optimal Superimposed Training Sequences for Channel Estimation in MIMO-OFDM Systems," *EURASIP Journal on Advances in Signal Processing*, 2010. https://link.springer.com/article/10.1155/2010/140506
5. D. Katselis et al., "Training sequence design for MIMO channels: an application-oriented approach," *EURASIP Journal on Wireless Communications and Networking*, 2013. https://link.springer.com/article/10.1186/1687-1499-2013-245
6. G. Li and G. Liao, "A Pilot-Pattern Based Algorithm for MIMO-OFDM Channel Estimation," *Algorithms*, 2017. https://www.mdpi.com/1999-4893/10/1/3
7. MathWorks, "Estimate the wireless channel in a MIMO OFDM system using pilot signals." https://www.mathworks.com/help/comm/ref/ofdmchannelestimate.html
