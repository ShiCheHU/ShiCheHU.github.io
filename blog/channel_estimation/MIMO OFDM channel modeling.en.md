# MIMO OFDM Channel Modeling

> The **Pilot Design and Channel Estimation** section has been split into a standalone document: [MIMOOFDMChannelEstimation.md](MIMOOFDMChannelEstimation.md).

## Table of Contents

- [1 Overview](#1-overview)
- [2 MIMO-OFDM System Model](#2-mimo-ofdm-system-model)
  - [2.1 Parameters and Notation](#21-parameters-and-notation)
  - [2.2 Continuous-Domain MIMO CIR](#22-continuous-domain-mimo-cir)
    - [2.2.1 Time-Frequency Domain Relationship](#221-time-frequency-domain-relationship)
  - [2.3 Discrete-Domain MIMO CIR/CFR](#23-discrete-domain-mimo-circfr)
  - [2.4 OFDM Frequency-Domain Input-Output Model](#24-ofdm-frequency-domain-input-output-model)
- [3 MIMO Channel Modeling Methods](#3-mimo-channel-modeling-methods)
  - [3.1 Independent Rayleigh/Rician Tap Model](#31-independent-rayleighrician-tap-model)
  - [3.2 Spatially Correlated Kronecker Model](#32-spatially-correlated-kronecker-model)
    - [3.2.1 Basic Assumptions of the Kronecker Model](#321-basic-assumptions-of-the-kronecker-model)
    - [3.2.2 Discrete-Domain MIMO CIR Under the Kronecker Model](#322-discrete-domain-mimo-cir-under-the-kronecker-model)
    - [3.2.3 MIMO CFR Under the Kronecker Model](#323-mimo-cfr-under-the-kronecker-model)
    - [3.2.4 Construction and Derivation of Spatial Correlation Matrices](#324-construction-and-derivation-of-spatial-correlation-matrices)
      - [3.2.4.1 Deriving Correlation Matrices from the Angular Power Spectrum](#3241-deriving-correlation-matrices-from-the-angular-power-spectrum)
      - [3.2.4.1b Angular Spread — Definition and Generation](#3241b-angular-spread--definition-and-generation)
      - [3.2.4.2 Exponential Correlation Model](#3242-exponential-correlation-model)
      - [3.2.4.3 Constructing Correlation Matrices Using Array Steering Vectors](#3243-constructing-correlation-matrices-using-array-steering-vectors)
    - [3.2.5 Model Pros and Cons](#325-model-pros-and-cons)
  - [3.3 Angular-Domain Channel Modeling and Sparse Representation](#33-angular-domain-channel-modeling-and-sparse-representation)
    - [3.3.1 Continuous Angular Domain — From Array Phase Differences to Steering Vectors](#331-continuous-angular-domain--from-array-phase-differences-to-steering-vectors)
    - [3.3.2 Discrete Angular Domain — DFT Basis Vectors and Virtual Channel Representation](#332-discrete-angular-domain--dft-basis-vectors-and-virtual-channel-representation)
    - [3.3.3 Joint Angle-Delay Domain Sparsity](#333-joint-angle-delay-domain-sparsity)
    - [3.3.4 Extension: UPA (Uniform Planar Array) Angular Domain Representation](#334-extension-upa-uniform-planar-array-angular-domain-representation)
  - [3.4 3GPP TR 38.901 CDL/TDL Models](#34-3gpp-tr-38901-cdltdl-models)
- [4 Simulation Workflow](#4-simulation-workflow)
- [5 Model Selection Recommendations](#5-model-selection-recommendations)
- [6 References](#6-references)

# 1 Overview

MIMO (Multiple-Input Multiple-Output) is a multi-antenna transmission/reception technology. By exploiting multiple antennas deployed at the transmitter and receiver, it provides degrees of freedom in the spatial domain, improving the capacity, spectral efficiency, and coverage of wireless communication systems. In high-frequency bands in particular, MIMO technology is especially important, as millimeter-wave and terahertz bands suffer from short coverage distances and MIMO can enhance coverage.

MIMO primarily leverages three techniques—spatial diversity, spatial multiplexing, and beamforming—to achieve its advantages. Among these, spatial diversity and spatial multiplexing can be categorized under space-time block coding techniques. These three techniques are introduced below:

- **Spatial diversity**: The same data is transmitted over different channels. This way, multiple copies of the same data are available, and the receiver can combine data from multiple channels to effectively boost the signal-to-noise ratio.
- **Spatial multiplexing**: A single data stream is split into segments, and different segments are sent over different channels. Even if the channel quality for some data is poor, other data transmissions enjoy good quality. This effectively mitigates channel impairments.
- **Beamforming**: By adjusting the antenna array, the signals transmitted by multiple antennas add constructively to form a directional signal, thereby improving the coverage range in a given direction. This is employed in massive MIMO systems.

MIMO technology finds application in cellular systems and Wi-Fi:

- In 5G NR, hybrid beamforming is used in the FR2 band, while digital beamforming is used in the sub-6 GHz band.
- In Wi-Fi systems, space-time block coding is primarily adopted.

The evolution of MIMO in wireless technology spans 5G standardization; increasing antenna counts, different antenna panels, and coordination across transmission points are the development trends. Heading into 6G, MIMO remains a major topic of discussion, including: channel design, beam management, CSI acquisition, uplink coverage, initial access, mobility management, and more. The evolution of MIMO technology faces the following challenges:

- Excessively high channel acquisition overhead due to increasing antenna counts. More resources must be allocated to pilots in both uplink and downlink; since pilot resources and data are orthogonal, the pilot resources used for channel acquisition become too high.
- Increased hardware cost. Supporting more data streams and more users requires a linearly increasing number of RF chains, which in turn increases cost and power consumption.
- Increased data processing complexity. With more antennas and the use of hybrid beamforming, beam management is needed to avoid beam misalignment. With more data streams, the complexity of channel estimation and real-time data reception at the receiver increases.
- Multi-user interference cancellation and coverage. Multi-user operation was extensively discussed in 5G standards and is supported in the specifications, but has not been widely deployed in practice, since multi-user scenarios can be fully supported via spatial, time, and frequency division multiplexing. Multi-user precoding is expensive and can only be applied to slow-moving devices (real-time precoding updates are costly for fast mobility) and nearby devices (spatially separated at longer distances). In the 6G phase, supporting even larger numbers of users—as proposed by many companies—will be very difficult.

This article primarily introduces MIMO OFDM channel modeling. While SISO OFDM channel modeling has been covered previously, the extension to MIMO OFDM is fundamentally different. Compared with SISO channel modeling, MIMO OFDM channel modeling has the following new characteristics:

- SISO OFDM systems mainly consider frequency and time selectivity, whereas MIMO OFDM systems must additionally consider spatial selectivity, typically modeled via AoA/AoD.
- MIMO OFDM systems need to account for inter-antenna coupling effects, e.g., the influence of antenna spacing.
- Variations in the rank and condition number of the channel matrix must be considered for rank adaptation.
- Pilot design also needs to account for orthogonality across different antenna ports.

In SISO OFDM, the core of channel estimation is to recover the scalar channel frequency response for each OFDM symbol and each subcarrier:

$$
H[m,k].
$$

In MIMO-OFDM, the scalar channel extends to a matrix channel:

$$
\boxed{
\mathbf H[m,k]\in \mathbb C^{N_r\times N_t}
}
$$

where $N_t$ is the number of transmit antennas and $N_r$ is the number of receive antennas. The matrix element $H_{r,t}[m,k]$ denotes the SISO subchannel from the $t$-th transmit antenna to the $r$-th receive antenna.

Therefore, MIMO-OFDM channel modeling must simultaneously describe four types of correlation:

| Dimension | Physical Origin | Mathematical Object |
|-----------|----------------|---------------------|
| Delay correlation | Multipath delay spread | PDP, CIR tap |
| Frequency correlation | Finite delay spread | CFR subcarrier correlation |
| Time correlation | Doppler spread | Jakes/AR/BEM process |
| Spatial correlation | Multi-antenna array, AoA/AoD, angular spread | Array response, spatial covariance, angular-domain sparsity |

The key new issue introduced by MIMO compared with SISO is that subchannels between different transmit and receive antennas are generally not independent. Spatial correlation affects:

- Channel matrix rank and condition number;
- Number of spatial multiplexing layers;
- Diversity gain;
- Precoding and detection performance;
- Pilot overhead required for channel estimation.

# 2 MIMO-OFDM System Model

## 2.1 Parameters and Notation

Let:

- Number of transmit antennas: $N_t$
- Number of receive antennas: $N_r$
- Number of OFDM subcarriers: $N$
- Number of OFDM symbols: $M$
- Subcarrier spacing: $\Delta f$
- Useful OFDM symbol duration: $T_u=1/\Delta f$
- Cyclic prefix length: $T_{cp}$
- Total OFDM symbol duration: $T_{\mathrm{sym}}=T_u+T_{cp}$
- Sampling interval: $T_s=1/(N\Delta f)$
- Number of paths or discrete taps: $L$
- Time instant of the $m$-th OFDM symbol: $t_m=mT_{\mathrm{sym}}$

The following is typically required:

$$
\boxed{
T_{cp}>\tau_{\max}
}
$$

so that after CP removal, linear convolution reduces to circular convolution, and the MIMO-OFDM channel on each subcarrier can be approximately decoupled into a narrowband MIMO flat-fading channel.

## 2.2 Continuous-Domain MIMO CIR

The continuous-time, time-varying impulse response from the $t$-th transmit antenna to the $r$-th receive antenna is:

$$
h_{r,t}(t,\tau)=\sum_{\ell=0}^{L-1}\alpha_{\ell}^{(r,t)}(t)\delta(\tau-\tau_\ell).
$$

Stacking all antenna pairs yields matrix form:

$$
\boxed{
\mathbf H(t,\tau)=\sum_{\ell=0}^{L-1}\mathbf G_\ell(t)\delta(\tau-\tau_\ell)
}
$$

where:

$$
\mathbf G_\ell(t)=
\begin{bmatrix}
\alpha_\ell^{(1,1)}(t) & \cdots & \alpha_\ell^{(1,N_t)}(t)\\
\vdots & \ddots & \vdots\\
\alpha_\ell^{(N_r,1)}(t) & \cdots & \alpha_\ell^{(N_r,N_t)}(t)
\end{bmatrix}
\in \mathbb C^{N_r\times N_t}.
$$

If the $\ell$-th path has power $P_\ell$, one often writes:

$$
\mathbf H(t,\tau)=\sum_{\ell=0}^{L-1}\sqrt{P_\ell}\mathbf G_\ell(t)\delta(\tau-\tau_\ell),
\qquad \sum_{\ell=0}^{L-1}P_\ell=1.
$$

Here $\mathbf G_\ell(t)$ describes spatial fading and time selectivity, while $P_\ell$ describes the power delay profile.

### 2.2.1 Time-Frequency Domain Relationship

In the time domain, the MIMO channel input-output relation is a **matrix convolution**, while in the frequency domain it becomes matrix multiplication. We give a rigorous derivation from the LTI system perspective.

**Time-domain convolution form**: For an MIMO LTI system with $N_t$ inputs and $N_r$ outputs, the channel is a matrix-valued impulse response $\mathbf H(t, \tau)$. The time-domain input $\mathbf x(t) \in \mathbb{C}^{N_t}$ and output $\mathbf y(t) \in \mathbb{C}^{N_r}$ are related by matrix convolution:

$
\boxed{
\mathbf y(t) = \int_{-\infty}^{\infty} \mathbf H(t, \tau)\,\mathbf x(t - \tau)\,d\tau
}
$

Writing component-wise:

$
y_r(t) = \sum_{t=1}^{N_t} \int_{-\infty}^{\infty} h_{r,t}(t, \tau)\,x_t(t - \tau)\,d\tau,
\quad r = 1, \dots, N_r.
$

**Key observation**: In component form, each $(r, t)$ pair is an ordinary convolution of scalar-valued functions. The inner integral and the finite antenna summation commute; hence one can apply the Fourier transform to each component **separately**:

$
\begin{aligned}
Y_r(f) &= \sum_{t=1}^{N_t} \mathcal{F}\!\left\{ \int h_{r,t}(t, \tau)\,x_t(t - \tau)\,d\tau \right\} \\[4pt]
&= \sum_{t=1}^{N_t} H_{r,t}(f)\,X_t(f).
\end{aligned}
$

Recombining all scalar relations for all $r$ and $t$ as matrix multiplication gives the frequency-domain matrix form:

$
\boxed{
\mathbf Y(f) = \mathbf H(f)\,\mathbf X(f)
}
$

**Why cannot the Fourier transform be directly applied to a vector-valued function?**

Mathematically, scalar-valued functions (absolutely integrable / square-integrable / Schwartz class) form an infinite-dimensional vector space, and the Fourier transform is a linear operator on that space. However, the set of vector-valued functions with codomain $\mathbb{C}^{N_r}$ does not form a vector space; it is a module over the function ring $\mathcal{A}$ (an $\mathcal{A}$-module). The orthogonal decomposition theory of vector spaces cannot be directly applied to modules. Therefore, the operation of "applying FT component-wise to a vector-valued function and then assembling the components into a vector" relies on the structure of the MIMO **convolution** system (i.e., each SISO subchannel allows commutativity of FT and summation), not on a general property of the Fourier transform.

**Conclusion for MIMO-OFDM**: In an OFDM system, the channel on each subcarrier is narrowband flat; hence the frequency-domain matrix multiplication $\mathbf y[m,k] = \mathbf H[m,k]\,\mathbf x[m,k]$ holds at each time-frequency resource. This is the foundation for all subsequent MIMO channel modeling and estimation algorithms.

## 2.3 Discrete-Domain MIMO CIR/CFR

Under the integer-sampled tap model, the discrete MIMO CIR is:

$$
\boxed{
\mathbf H[m,n]\in \mathbb C^{N_r\times N_t},\qquad n=0,1,\dots,N-1
}
$$

where only finitely many taps are nonzero. The MIMO CFR at the $m$-th OFDM symbol and $k$-th subcarrier is the DFT of the discrete CIR along the delay tap dimension. Since $\mathbf H[m,n]$ is a matrix, the rigorous writing should first fix a transmit-receive antenna pair $(r,t)$:

$$
\boxed{
H_{r,t}[m,k]
=
\sum_{n=0}^{N-1}
H_{r,t}[m,n]\ e^{-j\frac{2\pi}{N}kn},
\qquad r=1,\dots,N_r,\ t=1,\dots,N_t.
}
$$

Rearranging all $(r,t)$ results into a matrix yields $\mathbf H[m,k]\in\mathbb C^{N_r\times N_t}$. The compact matrix notation

$$
\mathbf H[m,k]
=
\sum_{n=0}^{N-1}\mathbf H[m,n]e^{-j\frac{2\pi}{N}kn}
$$

means multiplying each element of $\mathbf H[m,n]$ by the same DFT phase scalar, rather than multiplying a matrix by a DFT vector.

The inverse transform is likewise defined per antenna pair:

$$
\boxed{
H_{r,t}[m,n]
=
\frac{1}{N}\sum_{k=0}^{N-1}
H_{r,t}[m,k]\ e^{j\frac{2\pi}{N}kn}
}
$$

If a continuous-delay path model is adopted, then for each antenna pair:

$$
\boxed{
H_{r,t}[m,k]
=
\sum_{\ell=0}^{L-1}
\sqrt{P_\ell}\,[\mathbf G_\ell[m]]_{r,t}\ e^{-j2\pi k\Delta f\tau_\ell}
}
$$

This equation is the parallel application of the SISO CFR to all antenna pairs. Writing in matrix form, it should be understood as element-wise superposition of path phases.

## 2.4 OFDM Frequency-Domain Input-Output Model

At the $m$-th OFDM symbol and $k$-th subcarrier, the transmit vector is:

$$
\mathbf x[m,k]\in \mathbb C^{N_t\times 1}.
$$

The receive vector is:

$$
\mathbf y[m,k]\in \mathbb C^{N_r\times 1}.
$$

If the CP is sufficiently long and the channel is approximately constant within one OFDM symbol, the frequency-domain model is:

$$
\boxed{
\mathbf y[m,k]=\mathbf H[m,k]\mathbf x[m,k]+\mathbf w[m,k]
}
$$

where $\mathbf w[m,k]\sim\mathcal{CN}(\mathbf 0,\sigma_w^2\mathbf I_{N_r})$.

Expanding:

$$
y_r[m,k]=\sum_{t=1}^{N_t}H_{r,t}[m,k]x_t[m,k]+w_r[m,k].
$$

If multiple spatial streams are sent simultaneously on one subcarrier, $\mathbf H[m,k]$ also enters MIMO detection or equalization:

$$
\hat{\mathbf x}[m,k]=\mathbf W[m,k]\mathbf y[m,k].
$$

Common linear equalizers are:

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

Thus channel estimation errors not only affect phase compensation but also spatial stream separation.

# 3 MIMO Channel Modeling Methods

## 3.1 Independent Rayleigh/Rician Tap Model

The independent tap model is the most fundamental approach in MIMO channel modeling. Its core assumptions are:

- Different delay taps are statistically independent (WSSUS assumption);
- Within the same tap, different antenna pairs $(r, t)$ are also statistically independent—i.e., array geometry and spatial correlation are ignored.

The idea of this model traces back to Bello's (1963) WSSUS channel framework and has been widely used as the default assumption in standardization efforts such as COST 207 (1986). It remains the baseline of choice for algorithm validation to this day.

### 3.1.1 Construction of the Multipath Channel

The system first specifies a Power Delay Profile (PDP), giving a set of delay–power pairs:

$
\{(\tau_\ell, P_\ell)\}_{\ell=0}^{L-1}, \qquad \sum_{\ell=0}^{L-1} P_\ell = 1.
$

These are discretized to integer tap indices $n_\ell = \mathrm{round}(\tau_\ell / T_s)$ at sampling interval $T_s$.

The MIMO spatial channel matrix for each path $\ell$ is:

$
\boxed{
\mathbf G_\ell[m] \in \mathbb C^{N_r \times N_t}
}
$

Its elements $[\mathbf G_\ell[m]]_{r,t}$ are independently drawn from a zero-mean complex Gaussian distribution with variance 1. The MIMO CIR at the $m$-th OFDM symbol and $n$-th sample tap is:

$
\boxed{
\mathbf H[m, n] = \sqrt{P_n}\ \mathbf G_n[m],
\qquad \mathbf G_n[m] = \mathbf G_{\sigma(n)}[m]
}
$

where $n$ is related to $\ell$ through the tap-index mapping in the PDP. If multiple paths overlap on the same tap, its power $P_n$ is the sum of the corresponding path powers.

### 3.1.2 Probabilistic Generation of LOS Paths

Physically, the first-arriving path (minimum propagation delay) is most likely to contain an LOS component; subsequent paths, after multiple reflections and scattering, have a very low probability of being LOS. Hence in engineering practice, Rician fading is typically applied **only to the first few taps with a certain probability, while all remaining taps are treated as Rayleigh**.

A concrete approach (LOS at first tap, NLOS otherwise) is:

- Draw a Bernoulli random variable $b_{\mathrm{LOS}} \in \{0, 1\}$, $\Pr(b_{\mathrm{LOS}} = 1) = p_{\mathrm{LOS}}$.
- If $b_{\mathrm{LOS}} = 1$, the first tap ($\ell = 0$) uses **Rician** fading and subsequent taps use **Rayleigh**;
- If $b_{\mathrm{LOS}} = 0$, all paths use **Rayleigh**.

$p_{\mathrm{LOS}}$ may be a fixed value (e.g., 20%) or a distance-dependent function (see 3GPP TR 38.901 scenario definitions).

### 3.1.3 Rician Spatial Matrix for a Single Path

For a path $\ell$ designated as LOS, the $(r, t)$-th element of its spatial matrix $\mathbf G_\ell[m]$ is generated according to the Rician distribution:

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

where:

| Symbol | Meaning |
|--------|---------|
| $K_\ell$ | Rician K-factor (linear) of the $\ell$-th path, $K_\ell \ge 0$ |
| $\alpha_{\mathrm{LOS}}$ | Deterministic LOS component; typically $|\alpha_{\mathrm{LOS}}| = 1$, with phase depending on AoA/AoD geometry; in simplified models it may be taken as a constant |
| $\alpha_{\mathrm{NLOS}}$ | Random scattering component, independently sampled per antenna pair |

When $K_\ell = 0$, Rician degenerates to Rayleigh—the scattering component variance is 1. Hence **Rayleigh can be regarded as the special case of Rician with $K=0$**, and the two can be handled uniformly in implementation.

Typical K-factor values:

| Scenario | $K_\ell$ (dB) | Meaning |
|----------|---------------|---------|
| Urban macrocell LOS | 4–10 dB | Weak LOS |
| Rural open area LOS | 10–20 dB | Strong LOS |
| Indoor LOS | 5–10 dB | Moderate LOS |
| Any NLOS path | $-\infty$ dB (i.e., 0) | Rayleigh |

### 3.1.4 Assembly of MIMO CFR

Superposing all path contributions in the frequency domain yields the MIMO CFR at the $m$-th OFDM symbol and $k$-th subcarrier. For a fixed antenna pair $(r,t)$:

$$
\boxed{
H_{r,t}[m, k]
=
\sum_{\ell=0}^{L-1}
\sqrt{P_\ell}\ [\mathbf G_\ell[m]]_{r,t}\ e^{-j 2\pi k \Delta f \tau_\ell}
}
$$

All $H_{r,t}[m,k]$ together form the matrix $\mathbf H[m,k]$. In discrete tap form, element-wise:

$$
[\mathbf H[m, k]]_{r,t}
=
\sum_{n=0}^{N-1}
[\mathbf H[m, n]]_{r,t}\ e^{-j\frac{2\pi}{N}kn}.
$$

### 3.1.5 Complete Generation Procedure

For each OFDM symbol $m$:

1. **Set PDP**: determine the number of paths $L$, delays $\tau_\ell$, and powers $P_\ell$.
2. **Determine LOS/NLOS**: draw according to $p_{\mathrm{LOS}}$ to decide whether the first tap (or first several taps) is LOS.
3. **Generate spatial matrix $\mathbf G_\ell[m]$ per path**:
   - LOS paths: generated according to the Rician formula, with K-factor taken as a scenario-dependent constant;
   - NLOS paths: all elements independent $\sim \mathcal{CN}(0, 1)$ (Rayleigh).
4. **Frequency-domain superposition**: sum over subcarriers according to the CFR formula.
5. *Optional*: if time selectivity is to be simulated, the phase rotation of $\alpha_{\mathrm{NLOS}}$ and $\alpha_{\mathrm{LOS}}$ across OFDM symbols can be driven via Jakes or AR(1) models.

### 3.1.6 Model Pros and Cons

**Pros**:

- Simple to implement, few parameters (only PDP + K-factor + $p_{\mathrm{LOS}}$), suitable for rapid algorithm validation;
- Independent taps allow direct reuse of the SISO multipath theoretical framework;
- Per-antenna-pair independence makes it the baseline for capacity and BER upper bounds.

**Cons**:

- Completely ignores array geometry and spatial correlation; optimistic when antenna count is large or spacing is small;
- Independence across different transmit-receive antenna pairs cannot describe AoD-AoA coupling, beam effects, or the spatial structure of precoding gains;
- When $N_t, N_r$ are large, the independence assumption can lead to significant overestimation of channel rank, distorting spatial multiplexing layer predictions.

## 3.2 Spatially Correlated Kronecker Model

The Kronecker model is one of the most classic analytical models in MIMO channel modeling, proposed and systematically analyzed around 2000 by Shiu, Foschini, Gans, Kahn, and others. Its core idea is to use a **separable spatial correlation structure** to characterize the statistical dependence among multiple antennas, striking a good balance between theoretical elegance and physical accuracy.

### 3.2.1 Basic Assumptions of the Kronecker Model

The core assumptions of the Kronecker model can be summarized as three items:

**Assumption 1: Transmit-side and receive-side spatial correlations are separable**

This is the most fundamental assumption of the Kronecker model. It requires that the cross-correlation between any two transmit antennas $(t_1, t_2)$ and any two receive antennas $(r_1, r_2)$ factorizes as the product of transmit correlation and receive correlation:

$$
\mathbb E\left[ H_{r_1,t_1} H_{r_2,t_2}^* \right]
=
[\mathbf R_t]_{t_1,t_2} \cdot [\mathbf R_r]_{r_1,r_2}.
$$

In other words, the channel correlation of an antenna pair $(r_1, t_1)$ and $(r_2, t_2)$ is **entirely determined by the product of the spatial correlation between the two transmit antennas** and **the spatial correlation between the two receive antennas**, with no mutual interference.

This assumption is equivalent to: the covariance matrix of the vectorized channel $\mathbf h = \mathrm{vec}(\mathbf H)$ has a Kronecker product structure:

$$
\mathbf R_h
=
\mathbb E\left[ \mathbf h \mathbf h^H \right]
=
\mathbf R_t^T \otimes \mathbf R_r.
$$

**Assumption 2: The scattering environments at transmitter and receiver are statistically independent**

Physically, this means the scatterer distribution around the transmitter and the scatterer distribution around the receiver are mutually independent. Mathematically, this ensures that $\mathbf R_t$ and $\mathbf R_r$ can be independently defined. If common scatterers exist between the transmitter and receiver (e.g., buildings on both sides of a street canyon simultaneously affecting AoD and AoA), this assumption is violated.

**Assumption 3: Each tap independently satisfies the above separable structure**

For wideband multipath channels, the Kronecker model is typically applied per tap. That is, each delay tap $n$ may have different correlation matrices $\mathbf R_{t,n}$ and $\mathbf R_{r,n}$. This allows different delay clusters to exhibit different angular spread characteristics (e.g., the first-arriving cluster has small angular spread while later clusters have larger angular spread). However, if correlations also exist between different taps, the Kronecker model cannot directly capture them.

### 3.2.2 Discrete-Domain MIMO CIR Under the Kronecker Model

In §2.3, the discrete MIMO CIR was defined as $\mathbf H[m, n] \in \mathbb C^{N_r \times N_t}$, the channel matrix at the $m$-th OFDM symbol and $n$-th sample tap. Under the Kronecker model, this matrix is generated as follows:

$$
\boxed{
\mathbf H[m,n]
=
\sqrt{P_n}\ \mathbf R_{r,n}^{1/2}\ \mathbf W_n[m]\ \mathbf R_{t,n}^{1/2}
}
$$

where the symbols are:

| Symbol | Dimensions | Meaning |
|--------|------------|---------|
| $P_n$ | scalar | Average power of the $n$-th tap, given by the PDP |
| $\mathbf W_n[m]$ | $\mathbb C^{N_r \times N_t}$ | i.i.d. zero-mean unit-variance complex Gaussian matrix, $\mathrm{vec}(\mathbf W_n[m]) \sim \mathcal{CN}(\mathbf 0, \mathbf I_{N_r N_t})$ |
| $\mathbf R_{t,n}$ | $\mathbb C^{N_t \times N_t}$ | Transmit-side spatial correlation matrix, describing statistical dependence among $N_t$ transmit antennas |
| $\mathbf R_{r,n}$ | $\mathbb C^{N_r \times N_r}$ | Receive-side spatial correlation matrix, describing statistical dependence among $N_r$ receive antennas |

Here $\mathbf R^{1/2}$ denotes the Hermitian square root, satisfying $\mathbf R^{1/2} (\mathbf R^{1/2})^H = \mathbf R$.

Column-wise vectorization $\mathbf h_n[m] = \mathrm{vec}\{\mathbf H[m,n]\}$, using the identity $\mathrm{vec}(\mathbf A \mathbf B \mathbf C) = (\mathbf C^T \otimes \mathbf A)\,\mathrm{vec}(\mathbf B)$, yields:

$$
\begin{aligned}
\mathbf h_n[m]
&=
\sqrt{P_n}\ \mathrm{vec}\left( \mathbf R_{r,n}^{1/2}\ \mathbf W_n[m]\ \mathbf R_{t,n}^{1/2} \right) \\[4pt]
&=
\sqrt{P_n}\ \left( (\mathbf R_{t,n}^{1/2})^T \otimes \mathbf R_{r,n}^{1/2} \right)\ \mathrm{vec}(\mathbf W_n[m]).
\end{aligned}
$$

Let $\mathbf w_n[m] = \mathrm{vec}(\mathbf W_n[m]) \sim \mathcal{CN}(\mathbf 0, \mathbf I)$. Then $\mathbf h_n[m]$ is a linear transformation of $\mathbf w_n[m]$. By the linear transformation property of complex Gaussian vectors, $\mathbf h_n[m]$ remains zero-mean complex Gaussian, with covariance:

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

Using $(\mathbf R_{t,n}^{1/2})^T (\mathbf R_{t,n}^{1/2})^* = ( \mathbf R_{t,n}^{1/2} (\mathbf R_{t,n}^{1/2})^H )^T = \mathbf R_{t,n}^T$, we finally obtain the Gaussian distribution of the vectorized channel:

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

This is the vectorized form of the Kronecker model and the basis for constructing $\mathbf R_h$ in subsequent LMMSE channel estimation.

**Complete MIMO CIR generation procedure** (for each OFDM symbol $m$ and each tap $n$):

1. Generate i.i.d. matrix $\mathbf W_n[m] \in \mathbb C^{N_r \times N_t}$, each element independent $\sim \mathcal{CN}(0, 1)$;
2. Left-multiply by $\mathbf R_{r,n}^{1/2}$ to impose receive-side spatial correlation;
3. Right-multiply by $\mathbf R_{t,n}^{1/2}$ to impose transmit-side spatial correlation;
4. Multiply by $\sqrt{P_n}$ to assign the tap power.

### 3.2.3 MIMO CFR Under the Kronecker Model

From §2.3, the MIMO CFR at the $m$-th OFDM symbol and $k$-th subcarrier is the DFT of the CIR:

$$
\mathbf H[m, k] = \sum_{n=0}^{N-1} \mathbf H[m, n]\ e^{-j\frac{2\pi}{N}kn}.
$$

Substituting the Kronecker-model CIR and assuming independence among taps (scatterers at different delays are statistically independent), we obtain:

$$
\boxed{
\mathbf H[m, k]
=
\sum_{n=0}^{N-1}
\sqrt{P_n}\ \mathbf R_{r,n}^{1/2}\ \mathbf W_n[m]\ \mathbf R_{t,n}^{1/2}
\ e^{-j\frac{2\pi}{N}kn}
}
$$

If all taps share the same spatial correlation characteristics (i.e., all taps share a single set of correlation matrices $\mathbf R_t, \mathbf R_r$), this simplifies to:

$$
\mathbf H[m, k]
=
\mathbf R_r^{1/2}
\left(
\sum_{n=0}^{N-1} \sqrt{P_n}\ \mathbf W_n[m]\ e^{-j\frac{2\pi}{N}kn}
\right)
\mathbf R_t^{1/2}.
$$

This form clearly shows that under the Kronecker assumption, the spatial correlation structure and the frequency selectivity are **decoupled**—the correlation matrices do not vary with subcarrier; only the internal i.i.d. fading part is affected by DFT phase rotation.

Vectorizing the CFR likewise yields its statistical distribution. Since linear transformations preserve Gaussianity and the taps are independent, the covariance of the vectorized CFR is the weighted sum of the per-tap covariances:

$$
\mathrm{vec}\{\mathbf H[m, k]\}
\sim
\mathcal{CN}\left(
\mathbf 0,\;
\sum_{n=0}^{N-1} P_n \left( \mathbf R_{t,n}^T \otimes \mathbf R_{r,n} \right)
\right).
$$

If all taps share the same correlation matrices, the covariance simplifies to $(\mathbf R_t^T \otimes \mathbf R_r)$ (since $\sum_n P_n = 1$), and the statistical properties of the CFR coincide with those of a single CIR tap—a convenient property of the Kronecker model.

### 3.2.4 Construction and Derivation of Spatial Correlation Matrices

The practicality of the Kronecker model depends on being able to construct $\mathbf R_t$ and $\mathbf R_r$ in a physically reasonable way. Below we give the complete derivation chain from physical parameters to correlation matrices.

#### 3.2.4.1 Deriving Correlation Matrices from the Angular Power Spectrum

Consider a 1D uniform linear array (ULA) with antenna spacing $d$ and wavelength $\lambda$. Define the normalized direction cosine (or angular frequency):

$$
\Omega = \frac{d}{\lambda} \sin\theta,
$$

where $\theta \in [-\pi/2, \pi/2]$ is the physical angle of arrival (or departure) relative to the array broadside. For a half-wavelength spacing $d = \lambda/2$, $\Omega \in [-0.5, 0.5]$.

On this array, the array steering vector for a narrowband plane wave at angle $\theta$ is:

$$
\mathbf a(\Omega)
=
\frac{1}{\sqrt{N}}
\begin{bmatrix}
1 \\ e^{-j2\pi\Omega} \\ e^{-j2\pi\cdot 2\Omega} \\ \vdots \\ e^{-j2\pi(N-1)\Omega}
\end{bmatrix}.
$$

Physically, the $\ell$-th path (or a ray within a scattering cluster) arrives at the array at angle $\theta_\ell$, contributing a phase $e^{-j2\pi p\,\Omega_\ell}$ to the $p$-th antenna. The covariance matrix of the received signal is the statistical average of all path contributions.

Assuming scatterers are continuously distributed in the angular domain, with their power density described by the **Angular Power Spectrum (APS)** $p(\Omega)$, normalized as $\int p(\Omega)\,d\Omega = 1$. Then the $(p, q)$-th element of the receive spatial correlation matrix (complex correlation between the $p$-th and $q$-th antennas) is:

$$
\boxed{
[\mathbf R]_{p,q}
=
\mathbb E\left[ h_p h_q^* \right]
=
\int e^{-j2\pi(p-q)\Omega}\ p(\Omega)\ d\Omega
}
$$

**Derivation**:

Consider the channel coefficient at the $p$-th antenna as the superposition of contributions from all angular paths:

$$
h_p = \int \alpha(\Omega)\ e^{-j2\pi p\,\Omega}\ d\Omega,
$$

where $\alpha(\Omega)$ is the complex gain density, satisfying $\mathbb E[\alpha(\Omega)\alpha^*(\Omega')] = p(\Omega)\,\delta(\Omega - \Omega')$ (different angles are uncorrelated). Then:

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

**Key observations**:

- $[\mathbf R]_{p,q}$ depends only on the antenna index difference $\Delta = p - q$; this is the **Toeplitz** structure of a ULA and the manifestation of the **wide-sense stationary (WSS)** assumption in the spatial domain.
- The elements of the correlation matrix are essentially the inverse Fourier transform of $p(\Omega)$ (evaluated at argument $(p-q)$); thus $p(\Omega)$ completely determines the spatial correlation structure.
- This formula applies to both transmitter and receiver: for the transmitter, $p(\Omega)$ is the angular power spectrum of departure (APS$_{\text{Tx}}$); for the receiver, $p(\Omega)$ is the angular power spectrum of arrival (APS$_{\text{Rx}}$).

**Two common APS shapes and their corresponding correlation matrices**:

| Angular Distribution | $p(\Omega)$ | $[\mathbf R]_{p,q}$ |
|----------------------|-------------|---------------------|
| Uniform (omnidirectional scattering) | $p(\Omega) = 1$ (over normalized interval) | $\mathrm{sinc}(2(p-q)\Omega_{\max})$ |
| Truncated Gaussian | $p(\Omega) \propto \exp(-\frac{(\Omega-\bar\Omega)^2}{2\sigma_\Omega^2})$ | No closed form; requires numerical integration |
| Laplace | $p(\Omega) \propto \exp\!\left(-\frac{\sqrt{2}\,|\Omega-\bar{\Omega}|}{\sigma_\Omega}\right)$ | No closed form; requires numerical integration |
| Single point source (pure LOS) | $p(\Omega) = \delta(\Omega - \Omega_0)$ | $e^{-j2\pi(p-q)\Omega_0}$ (rank 1) |

The Laplace distribution agrees well with measurements in macrocell scenarios; truncated Gaussian is more common in microcells. In 3GPP TR 38.901, the angular offsets of rays within a cluster are typically modeled as Laplace or wrapped Gaussian distributions.

#### 3.2.4.1b Angular Spread — Definition and Generation

The parameter $\sigma_\Omega$ in the angular power spectrum $p(\Omega)$ is precisely the **Angular Spread** (AS) — it quantifies the degree of dispersion of scattering energy in the angular domain and is one of the most fundamental large-scale parameters in MIMO channel modeling.

**1) Rigorous definition of angular spread**

The Angular Spread is defined as the square root of the second central moment (RMS angular spread) of the angular power spectrum:

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

where $\bar\Omega$ is the **mean angle**, i.e., the first moment of the angular power distribution. The definition is identical in degrees after converting $\Omega$ to the physical angle $\theta$.

**Physical intuition**:

- Small $\sigma_{\mathrm{AS}}$ → scatterers are concentrated within a narrow angular range → signals received/transmitted by multiple antennas are highly correlated → the channel matrix tends toward low rank → weak spatial multiplexing capability but high beamforming gain;
- Large $\sigma_{\mathrm{AS}}$ → scatterers are widely distributed (rich scattering environment) → channels across different antennas tend toward independence → the channel matrix tends toward full rank → many spatial multiplexing layers.

Under common angular distributions, the relationship between $\sigma_\Omega$ and AS is:

| Angular Distribution | Relationship between $\sigma_{\mathrm{AS}}$ and parameter $\sigma_\Omega$ |
|----------------------|---------------------------------------------------------------------------|
| Truncated Gaussian | AS ≈ $\sigma_\Omega$ (exact when the truncation interval is much larger than $\sigma_\Omega$) |
| Laplace | AS = $\sigma_\Omega$ (exact equality — $\sigma_\Omega$ is the RMS angular spread) |
| Uniform (half-width $\Omega_{\max}$) | AS = $\Omega_{\max}/\sqrt{3}$ |

Thus in practice, specifying AS is equivalent to specifying the width parameter of $p(\Omega)$.

**2) How to generate per-path angles from mean angle and angular spread**

In simulation, the angle for each path (or cluster, or ray) is generated according to the following workflow:

**Step 1: Determine large-scale parameters**
- Select the scenario (UMa/UMi/RMa/InH, etc.) and propagation condition (LOS/NLOS);
- Look up or derive from empirical formulas the mean angle $\bar\Omega$ (or $\bar\theta$) and angular spread $\sigma_{\mathrm{AS}}$ for that scenario.

**Step 2: Choose an angular distribution**
- Most common choices: Laplace (macrocells) or wrapped Gaussian (microcells);
- For pure-LOS scenarios, the LOS path may take $\sigma_{\mathrm{AS}} \to 0$, reducing the angle to the deterministic value $\bar\Omega$.

**Step 3: Generate a specific angle for each path/ray**

Assuming a Laplace distribution, the angle of the $\ell$-th path (or the $q$-th ray within a cluster) is:

$
\Omega_\ell = \bar\Omega + \Delta\Omega_\ell,
$

where $\Delta\Omega_\ell$ is the angular offset drawn from the Laplace distribution:

$
p(\Delta\Omega) = \frac{1}{\sqrt{2}\,\sigma_{\mathrm{AS}}} \exp\!\left( -\frac{\sqrt{2}\,|\Delta\Omega|}{\sigma_{\mathrm{AS}}} \right).
$

If a Gaussian distribution is used instead:

$
p(\Delta\Omega) = \frac{1}{\sqrt{2\pi}\,\sigma_{\mathrm{AS}}} \exp\!\left( -\frac{(\Delta\Omega)^2}{2\sigma_{\mathrm{AS}}^2} \right),
\qquad \text{truncated to a reasonable angular range}.
$

**Step 4: Convert normalized angle to physical angle**

$\Omega = \frac{d}{\lambda}\sin\theta$ → $\theta = \arcsin\!\left(\frac{\lambda}{d}\,\Omega\right)$. For $d = \lambda/2$, $\theta = \arcsin(2\Omega)$.

**Complete generation example** (transmit Azimuth AoD, Laplace distribution):

$
\begin{aligned}
&\text{Given: } \bar\theta_{\mathrm{AoD}} = 30^\circ,\quad \sigma_{\mathrm{ASD}} = 5^\circ,\quad L = 3\ \text{paths}. \\[4pt]
&\text{step 1: } \bar\Omega_t = \frac{d}{\lambda}\sin\bar\theta_{\mathrm{AoD}} = 0.5 \times \sin 30^\circ = 0.25. \\[4pt]
&\text{step 2: } \sigma_\Omega = \frac{d}{\lambda}\sin(\sigma_{\mathrm{ASD}}) \approx \frac{d}{\lambda}\,\sigma_{\mathrm{ASD}}\ \text{(small-angle approximation)}. \\[4pt]
&\text{step 3: for each path }\ell,\ \Omega_{t,\ell} = \bar\Omega_t + \Delta\Omega_\ell,\quad \Delta\Omega_\ell \sim \mathrm{Laplace}(0, \sigma_\Omega/\sqrt{2}). \\[4pt]
&\text{step 4: } \theta_{t,\ell} = \arcsin\!\left(\frac{\lambda}{d}\,\Omega_{t,\ell}\right).
\end{aligned}
$

**3) Meaning of AS in multipath/multi-cluster scenarios**

In wideband MIMO channels, different delay clusters typically exhibit different AS values. TR 38.901 defines two levels of angular spread:

| Level | Parameter | Meaning |
|-------|-----------|---------|
| **Per-cluster angular spread** | $c_{\mathrm{ASD}}$, $c_{\mathrm{ASA}}$, $c_{\mathrm{ZSD}}$, $c_{\mathrm{ZSA}}$ | The angular distribution width of the 20 rays within a cluster |
| **Global angular spread** | ASD, ASA, ZSD, ZSA | The overall angular spread of all clusters in full space (large-scale parameter) |

The generation workflow is:
1. Obtain the **global AS** according to the scenario (e.g., ASD ≈ 22° for UMa NLOS);
2. Scale the normalized per-cluster angular offsets $c_{\mathrm{ASD}}$, etc., from the CDL tables using the global AS;
3. For each cluster, generate ray-specific angles using the cluster center angle as the mean and $c_{\mathrm{ASD}} \cdot \mathrm{ASD}_{\text{desired}}$ as the spread.

**4) Reference angular spread values for typical scenarios**

The following values are from 3GPP TR 38.901 (UMa scenario, azimuth AS, in degrees):

| Scenario | LOS/NLOS | ASD | ASA | ZSD | ZSA |
|----------|----------|-----|-----|-----|-----|
| UMa | LOS | $\lg \mathrm{ASD} \sim \mathcal{N}(1.06, 0.28)$ | $\lg \mathrm{ASA} \sim \mathcal{N}(1.60, 0.18)$ | $\lg \mathrm{ZSD} \sim \mathcal{N}(0.64, 0.32)$ | $\lg \mathrm{ZSA} \sim \mathcal{N}(0.90, 0.23)$ |
| UMa | NLOS | $\lg \mathrm{ASD} \sim \mathcal{N}(1.46, 0.28)$ | $\lg \mathrm{ASA} \sim \mathcal{N}(1.79, 0.20)$ | $\lg \mathrm{ZSD} \sim \mathcal{N}(0.81, 0.23)$ | $\lg \mathrm{ZSA} \sim \mathcal{N}(0.96, 0.17)$ |

Note: AS values here follow a log-normal distribution; during simulation one first samples $\lg\mathrm{AS}$ and then takes $10^{\lg\mathrm{AS}}$. The table gives the mean $\mu$ and standard deviation $\sigma$ of $\lg\mathrm{AS}$. For example, the log-mean ASD of UMa NLOS is 1.46, i.e., geometric mean ASD ≈ $10^{1.46} \approx 28.8°$ (this value fluctuates depending on the specific parameter set, consistent with the AS ranges in the CDL tables).

**5) Summary**

| Question | Answer |
|----------|--------|
| Physical meaning of AS | RMS width of the scattering energy in the angular domain |
| How to parameterize $p(\Omega)$ | Specify distribution type + mean angle $\bar\Omega$ + angular spread $\sigma_{\mathrm{AS}}$ |
| How to assign an angle to a single path | $\Omega_\ell = \bar\Omega + \Delta\Omega$, where $\Delta\Omega$ is drawn randomly from the chosen distribution |
| Link between AS and spatial correlation | $[\mathbf R]_{p,q} = \int e^{-j2\pi(p-q)\Omega}p(\Omega)d\Omega$; AS determines the width of $p(\Omega)$ → determines the decay rate of the correlation matrix |
| Link between AS and channel rank | Larger AS → weaker correlation → higher channel rank |

#### 3.2.4.2 Exponential Correlation Model

When exact matching to a specific angular distribution is not required—or when only the correlation strength needs to be tuned—the **Exponential Correlation Model** can be used as a simplification:

$$
\boxed{
[\mathbf R_t]_{p,q} = \rho_t^{|p-q|},
\qquad
[\mathbf R_r]_{p,q} = \rho_r^{|p-q|}
}
$$

where $0 \le \rho_t, \rho_r < 1$ are the correlation coefficients.

**Model properties**:

- $\rho_t = 0$: transmit antennas completely independent (i.i.d.);
- $\rho_t \to 1$: transmit antennas completely correlated (channel matrix rank degenerates to 1);
- Correlation decays exponentially with antenna index difference $|p-q|$—the physical intuition is that the farther apart the antennas, the more independent the channels;
- The correlation matrix is Toeplitz, Hermitian, and positive definite.

The rationale for this model: for a uniform scattering environment with limited angular spread and zero central angle, when $p(\Omega)$ is approximately Laplace, the Fourier transform evaluated at integer sample points decays approximately exponentially. However, note that the exponential correlation model is a **phenomenological approximation** and does not strictly correspond to any specific physical angular distribution.

**Relationship between $\rho$ and effective rank**:

Define the **Effective Rank** of the correlation matrix:

$$
r_{\mathrm{eff}}(\mathbf R) = \exp\left( -\sum_i \frac{\lambda_i}{\sum_j \lambda_j} \ln \frac{\lambda_i}{\sum_j \lambda_j} \right),
$$

where $\lambda_i$ are the eigenvalues of $\mathbf R$. When $\rho = 0$, all eigenvalues are equal, $r_{\mathrm{eff}} = N$; when $\rho \to 1$, the largest eigenvalue dominates, $r_{\mathrm{eff}} \to 1$. The rank of the channel matrix $\mathbf H = \mathbf R_r^{1/2} \mathbf W \mathbf R_t^{1/2}$ is bounded by $\min(r_{\mathrm{eff}}(\mathbf R_t), r_{\mathrm{eff}}(\mathbf R_r))$; thus larger $\rho$ means fewer supportable spatial multiplexing layers.

#### 3.2.4.3 Constructing Correlation Matrices Using Array Steering Vectors

Besides integration from the angular power spectrum, correlation matrices can also be constructed directly from discrete path parameters. This is more natural in geometric models. The function $p(\Omega)$ here is not a CIR; it does not describe the **delay positions** of path arrivals, nor the instantaneous complex gain phase of each path; rather, it describes the statistical distribution of channel power in the **angular direction**. It can be understood as the angular-domain power spectrum (APS), analogous to how the PDP describes power distribution over delay, but with the variable changed from $\tau$ to the normalized angle $\Omega$.

If the angular domain contains only $L$ discrete directions, $p(\Omega)$ can be written as a discrete angular distribution:

$$
p(\Omega) = \sum_{\ell=1}^{L} \gamma_\ell\ \delta(\Omega - \Omega_\ell),
\qquad \sum_\ell \gamma_\ell = 1,
$$

where $\gamma_\ell$ is the normalized power of the $\ell$-th path and $\Omega_\ell$ its normalized direction cosine. Substituting into the integral formula:

$$
[\mathbf R]_{p,q}
=
\sum_{\ell=1}^{L} \gamma_\ell\ e^{-j2\pi(p-q)\Omega_\ell}
=
\sum_{\ell=1}^{L} \gamma_\ell\
[\mathbf a(\Omega_\ell)]_p\ [\mathbf a(\Omega_\ell)]_q^*.
$$

Thus the correlation matrix can also be written as the weighted outer product of the steering vectors of each path:

$$
\boxed{
\mathbf R = \sum_{\ell=1}^{L} \gamma_\ell\ \mathbf a(\Omega_\ell)\ \mathbf a^H(\Omega_\ell)
}
$$

This is the **spectral representation** of the correlation matrix, clearly showing that the rank of $\mathbf R$ does not exceed the number of paths $L$. This form facilitates generating exact correlation matrices in simulation directly from given AoA/AoD sets, without numerical integration.

##### 3.2.4.3.1 How to Obtain the Angular Power Spectrum in Practice

In real systems, $p(\Omega)$ is usually not a directly known analytic function but is obtained through measurements, standardized channel models, or empirical assumptions.

**1) Estimated from array measurements**

If the receiver has an antenna array, the channel vector can be measured over multiple time snapshots or pilot resources:

$$
\mathbf h(t)
=
\sum_{\ell=1}^{L}
\alpha_\ell(t)\ \mathbf a(\Omega_\ell).
$$

First, estimate the spatial covariance:

$$
\hat{\mathbf R}
=
\frac{1}{T}\sum_{t=1}^{T}\mathbf h(t)\mathbf h^H(t).
$$

Then estimate the angular power spectrum from $\hat{\mathbf R}$. For a ULA, the correlation matrix elements and the APS approximately satisfy a Fourier dual relationship:

$$
[\mathbf R]_{p,q}
=
\int p(\Omega)e^{-j2\pi(p-q)\Omega}d\Omega.
$$

Hence the angular spectrum can be estimated using beam scanning, Bartlett, Capon/MVDR, MUSIC/ESPRIT, SAGE, or sparse recovery methods.

**2) Generated from standardized channel models**

In link-level simulation, the more common approach is to generate angular parameters from a standard model. For example, 3GPP TR 38.901 generates cluster/ray parameters—AoA, AoD, power, etc.—according to the scenario, LOS/NLOS state, delay spread, and angular spread. Given the path set

$$
\{(\Omega_\ell,\gamma_\ell)\}_{\ell=1}^{L},
$$

one can directly write:

$$
p(\Omega)=\sum_{\ell=1}^{L}\gamma_\ell\delta(\Omega-\Omega_\ell),
$$

or construct the correlation matrix:

$$
\mathbf R
=
\sum_{\ell=1}^{L}
\gamma_\ell\mathbf a(\Omega_\ell)\mathbf a^H(\Omega_\ell).
$$

**3) Via empirical distribution assumptions**

For algorithm validation only, one may directly assume an angular distribution. Examples:

- Uniform: indicates rich scattering and wide angular coverage;
- Gaussian: indicates limited angular spread around a dominant direction;
- Laplace: common in macrocell modeling;
- Single-point delta: indicates pure LOS or a strong dominant path.

A typical Laplace angular spectrum can be written as:

$$
p(\Omega)
\propto
\exp\left(
-\frac{\sqrt{2}\,|\Omega-\bar{\Omega}|}{\sigma_\Omega}
\right),
$$

where $\bar{\Omega}$ is the mean direction and $\sigma_\Omega$ describes the angular spread.

Thus a common modeling chain in practice is:

$$
\text{Measure channel}
\rightarrow
\hat{\mathbf R}
\rightarrow
\text{Angular spectrum estimation}
\rightarrow
p(\Omega)
$$

or:

$$
\text{Select scenario}
\rightarrow
\text{Generate AoA/AoD and angular spread}
\rightarrow
p(\Omega)\ \text{or}\ \mathbf R.
$$

### 3.2.5 Model Pros and Cons

**Pros**:

- **Theoretically concise**: the vectorized covariance is in Kronecker product form, facilitating analytic expressions for capacity and estimation error;
- **Few parameters**: only two $N_t \times N_t$ and $N_r \times N_r$ correlation matrices are needed, rather than a full $N_t N_r \times N_t N_r$ correlation matrix;
- **Separability**: transmit and receive correlations can be adjusted independently, suitable for studying asymmetric scenarios;
- **Closed form**: correlation matrices can be quickly generated from APS or exponential models without ray tracing.

**Cons**:

- **Cannot describe joint AoD-AoA coupling**: physically, a particular transmit angle often corresponds only to a specific receive angle (e.g., waveguide effects in street canyons); such coupling is averaged out under the separable Kronecker assumption, narrowing the channel's angle-angle bi-selectivity structure;
- **Predictions may be over- or under-optimistic for certain propagation scenarios**:
  - When strong coupling is actually present, the Kronecker model can **overestimate** channel capacity (since the separable assumption lowers the effective condition number);
  - In some scattering-rich scenarios, it can also **underestimate** capacity (since the separable assumption limits the mutual information upper bound);
- **Unsuitable for describing "keyhole" (pinhole) channels**: when both transmitter and receiver are in rich scattering environments but are coupled only through a narrow aperture, the channel matrix has rank 1 while the correlation matrices are full rank; the Kronecker model cannot reproduce this phenomenon;
- **Increasing bias with large arrays**: as the number of antennas grows, near-field effects and spherical wavefronts become non-negligible, and the accuracy of correlation matrix construction under the plane-wave assumption degrades.

To address the joint coupling issue, the Weichselberger model introduces transmit-receive joint eigenmodes via eigendecomposition, serving as a direct generalization of the Kronecker model. In practical engineering simulation, if explicit AoD-AoA coupling is required, it is advisable to use geometric models (as in §3.3.1 outside §3.2) or the 3GPP CDL model (§3.4).


## 3.3 Angular-Domain Channel Modeling and Sparse Representation

The previous two sections (§3.1 and §3.2) modeled MIMO channels from the perspective of statistical correlation matrices. This section approaches the problem from the **array signal processing** angle: exploiting the phase structure of antenna arrays to map continuous spatial angles (AoA/AoD) into the algebraic structure of the channel matrix. The core logical chain is:

$$
\text{ULA phase differences}\ \longrightarrow\ \text{steering vectors}\ \longrightarrow\ \text{continuous-angle multipath superposition}\ \longrightarrow\ \text{angular discretization}\ \longrightarrow\ \text{DFT basis matrices}\ \longrightarrow\ \text{virtual channel sparse representation}
$$

This framework is not only a channel modeling tool but also the theoretical foundation of compressive-sensing channel estimation.

### 3.3.1 Continuous Angular Domain — From Array Phase Differences to Steering Vectors

**1) Array response to a single plane wave**

Consider a receive ULA of $N$ elements with spacing $d$. A far-field narrowband plane wave impinges at angle $\theta$ (relative to the array broadside), at wavelength $\lambda$.

The **path-length difference** between adjacent elements is $d\sin\theta$, and the corresponding **phase difference** is:

$$
\Delta\phi = 2\pi \cdot \frac{d\sin\theta}{\lambda}
$$

Define the **direction cosine**:

$$
\boxed{
\Omega \triangleq \frac{d}{\lambda}\sin\theta,
\qquad \Omega \in \left[-\frac{d}{\lambda},\ \frac{d}{\lambda}\right]
}
$$

Taking element $0$ as reference (phase = 0), the cumulative phase difference of element $p$ relative to the reference is $p \cdot \Delta\phi = 2\pi p\Omega$. Hence the response of the $N$-element array in this incident direction—i.e., the **array steering vector**—is:

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

where $1/\sqrt{N}$ is a normalization factor ensuring $\|\mathbf a(\Omega)\|^2 = 1$.

**Key observation**: the $p$-th element of the steering vector is $e^{-j2\pi p\Omega}$—precisely a **complex exponential sequence** with antenna index $p$ as "discrete time" and $\Omega$ as "digital frequency". This form is the root cause of the subsequent DFT discretization.

**2) Multipath superposition — geometric multipath channel**

Suppose the channel has $L$ resolvable physical paths, where the $\ell$-th path has:
- Transmit departure angle direction cosine $\Omega_{t,\ell}$ and receive arrival angle direction cosine $\Omega_{r,\ell}$;
- Complex path gain $\alpha_\ell[m]$;
- Propagation delay $\tau_\ell$.

Then the MIMO frequency-domain channel matrix is the superposition of rank-1 components from all paths:

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

If Doppler is to be considered: $\alpha_\ell[m] = \beta_\ell\, e^{j2\pi\nu_\ell mT_{\mathrm{sym}}}$.

**Structure and rank**: the spatial component of each path, $\mathbf a_r(\Omega_{r,\ell})\mathbf a_t^H(\Omega_{t,\ell})$, is a **rank-1** matrix—it occupies only one direction in the $N_r \times N_t$-dimensional space. The entire channel matrix is the weighted sum of $L$ rank-1 matrices, so

$$
\operatorname{rank}(\mathbf H) \leq \min(N_t, N_r, L).
$$

The more spread out the path angles → the more distinct the directions of the rank-1 components → the higher the channel matrix rank and the greater the spatial multiplexing potential. Conversely, smaller angular spread (as in LOS-dominated scenarios) leads to a more nearly low-rank channel.

**Continuous angle**: in the above model, $\Omega_{t,\ell}$ and $\Omega_{r,\ell}$ can take **any continuous real values** within $[-d/\lambda,\ d/\lambda]$; this is the meaning of "continuous angular domain".

**3) Beamforming and the inner product of spatial signatures**

Another name for the steering vector is **spatial signature**—it is the "digital fingerprint" of a specific incident direction. Understanding the meaning of the inner product between steering vectors is the key bridge connecting continuous angles to the discrete DFT basis, and the motivation for discretizing angles naturally emerges from the array's physical structure.

**The beamforming problem**: suppose we wish to concentrate transmitted power in the direction $\Omega_0$. The transmit signal is a scalar $s$, the complex weighting coefficients on each antenna form a vector $\mathbf w \in \mathbb{C}^N$ satisfying the power constraint $\|\mathbf w\|^2 = 1$. The transmitted signal vector is $\mathbf w s$, and the received signal in the far-field direction $\Omega$ (ignoring path loss) is

$$
y(\Omega) = \mathbf a^H(\Omega)\,\mathbf w\,s.
$$

If the target direction is $\Omega_0$, the received power is proportional to $\bigl|\mathbf a^H(\Omega_0)\mathbf w\bigr|^2$. By the **Cauchy-Schwarz inequality**, under $\|\mathbf a(\Omega_0)\|^2 = 1$ and $\|\mathbf w\|^2 = 1$,

$$
\bigl|\mathbf a^H(\Omega_0)\mathbf w\bigr|^2 \le \|\mathbf a(\Omega_0)\|^2\ \|\mathbf w\|^2 = 1,
$$

with equality if and only if $\mathbf w$ is collinear with $\mathbf a(\Omega_0)$. Therefore the **optimal beamforming vector is the steering vector itself**:

$$
\boxed{\mathbf w_{\mathrm{opt}} = \mathbf a(\Omega_0)}.
$$

This is called **matched filter beamforming**: each antenna uses the conjugate of the steering vector element as its complex weight, so that the signals from all antennas add constructively in the target direction.

**Beam pattern**: when beamforming with $\mathbf a(\Omega_0)$, the gain in another direction $\Omega$ is the squared modulus of the inner product of two spatial signatures, i.e., the **array beam pattern**:

$$
\boxed{
G(\Omega_0, \Omega)
\triangleq \bigl|\mathbf a^H(\Omega_0)\,\mathbf a(\Omega)\bigr|^2
= \left|\frac{1}{N}\sum_{p=0}^{N-1} e^{-j2\pi p(\Omega - \Omega_0)}\right|^2
}.
$$

Summing the finite geometric series yields its closed form—the **square of the Dirichlet kernel**:

$$
G(\Omega_0, \Omega) = \left|
\frac{\sin\!\bigl[\pi N(\Omega - \Omega_0)\bigr]}
{N\sin\!\bigl[\pi(\Omega - \Omega_0)\bigr]}
\right|^2.
$$

The beam pattern has the following key properties:

- Peaks at $\Omega = \Omega_0$ with value $1$ (mainlobe peak);
- Has exact nulls at $\Omega = \Omega_0 \pm \frac{k}{N}$ ($k = 1, 2, \dots, N-1$)—the larger $N$, the denser the nulls and the sharper the beam;
- The mainlobe null-to-null width is approximately $2/N$ (in $\Omega$), so more antennas give stronger directivity and higher spatial angular resolution.

**Key insight—origin of orthogonality and motivation for discretization**: take the discrete angular grid $\Omega_i = i/N$ ($i = 0, 1, \dots, N-1$). Then for any two grid points $\Omega_i \neq \Omega_j$:

$$
\mathbf a^H(\Omega_i)\,\mathbf a(\Omega_j)
= \frac{1}{N}\sum_{p=0}^{N-1} e^{-j2\pi p(j-i)/N}
= \begin{cases}
1, & i = j, \\[2pt]
0, & i \neq j,
\end{cases}
$$

i.e., **spatial signatures at different discrete grid points are mutually orthogonal**. The $N$ spatial signatures $\{\mathbf a(\Omega_i)\}$ at the $N$ discrete angular grid points $\{\Omega_i = i/N\}_{i=0}^{N-1}$ form a **complete orthonormal basis** of $\mathbb{C}^N$.

This conclusion reveals the deep rationale behind angular discretization: the finite number of antennas $N$ means the array can provide independent spatial degrees of freedom in only $N$ orthogonal directions—and these $N$ orthogonal directions are precisely given by the DFT basis, not an arbitrary choice, but an intrinsic property of the array's physical structure. For the commonly used $d = \lambda/2$ configuration, directions with $|\Omega_i| \le 1/2$ correspond to physically realizable spatial pointing directions; those with $|\Omega_i| > 1/2$ correspond to the "invisible region", but must be retained as part of a complete basis (analogous to DFT components beyond the Nyquist frequency in digital signal processing).

Below we formalize this orthonormal basis as DFT basis matrices and systematically construct the discrete angular domain.

### 3.3.2 Discrete Angular Domain — DFT Basis Vectors and Virtual Channel Representation

In the previous subsection, the parameter $\Omega$ (direction cosine) of the steering vector $\mathbf a(\Omega)$ **varies continuously** over $[-d/\lambda,\ d/\lambda]$. This subsection answers a core question: **how to discretize the continuous angular parameter into a fixed set of orthogonal basis vectors, such that any steering vector can be expressed as a linear combination of these basis vectors?**

The answer lies in the profound connection between the algebraic form of steering vectors and the Discrete-Time Fourier Transform (DTFT).

**1) The steering vector as the DTFT kernel of a spatial sequence**

Recall the steering vector derived in §3.3.1:

$$
\mathbf a(\Omega)
= \frac{1}{\sqrt{N}}
\begin{bmatrix} 1 \\ e^{-j2\pi\Omega} \\ e^{-j2\pi\cdot 2\Omega} \\ \vdots \\ e^{-j2\pi(N-1)\Omega} \end{bmatrix}
$$

Its $p$-th element ($p = 0, 1, \dots, N-1$) is:

$$
[\mathbf a(\Omega)]_p = \frac{1}{\sqrt{N}}\ e^{-j2\pi\Omega p}
$$

Viewing the antenna index $p$ as a discrete "time" variable and $\Omega$ as a "digital frequency", this is precisely the **DTFT kernel** of a length-$N$ rectangular-window sequence $\frac{1}{\sqrt{N}}\operatorname{rect}_N[p]$ evaluated at frequency $\Omega$. In other words, $\mathbf a(\Omega)$ is nothing other than the **discrete observation of the complex exponential signal $e^{-j2\pi\Omega t}$ at $N$ equally spaced sampling points $p=0,1,\dots,N-1$**.

**2) Uniform frequency sampling $\to$ DFT basis**

Now sample $\Omega$ **uniformly**. Partition the interval $[0, 1)$ into $N$ equally spaced grid points:

$$
\boxed{
\Omega_i = \frac{i}{N},\qquad i = 0, 1, \dots, N-1
}
$$

The interval $[0,1)$ is taken because complex exponentials satisfy $\mathbf a(\Omega + 1) = \mathbf a(\Omega)$ (period $1$), so $N$ sampling points within $[0,1)$ cover one full period.

**Why exactly $N$ grid points?** Could one use a denser grid of $M > N$ points?

- **Constraint of a complete orthonormal basis**: the $\mathbf a(\Omega_i)$ are all vectors in $\mathbb{C}^N$, and $\mathbb{C}^N$ is $N$-dimensional. An orthonormal basis can have at most $N$ vectors. With $\Omega_i = i/N$, these $N$ steering vectors happen to be pairwise orthogonal and form an orthonormal basis. If $M > N$ grid points were used, the steering vectors would necessarily be linearly dependent—they would form an **overcomplete frame**, not a basis. $N$ is the "critical sampling" upper bound for an orthonormal basis.
- **Practical use of oversampling**: real systems (e.g., 5G NR) commonly use $O \times N$ beams (oversampling factor $O = 2$ or $4$). These oversampled beams remain orthogonal and provide denser angular candidate directions, reducing off-grid bias. The number of spatial degrees of freedom remains $N$—oversampling densifies the candidate directions without enabling more simultaneous orthogonal beams.

When $d = \lambda/2$ (the common half-wavelength spacing), $\Omega = \frac{1}{2}\sin\theta \in [-1/2,\ 1/2]$. Among the $\Omega_i$, those lying inside $[-1/2,\ 1/2]$ correspond to physically realizable spatial angles; those outside ($|\Omega_i| > 1/2$) do not correspond to any real physical angle, **but must be retained as part of the complete DFT spectral basis**—just as in digital signal processing frequency-domain analysis, the entire $[0, 2\pi)$ interval must be considered.

Substituting $\Omega_i = i/N$ one by one into $\mathbf a(\cdot)$:

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

Each $\mathbf a(\Omega_i)$ has its $p$-th element equal to $\frac{1}{\sqrt{N}} e^{-j2\pi p i/N}$—**this is exactly the $i$-th column of the conjugate transpose of the $N$-point DFT matrix**.

**3) Constructing the unitary DFT basis matrix**

Arranging all $N$ sampled steering vectors column-wise yields the **unitary DFT basis matrix**:

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

where $\omega \triangleq e^{j2\pi/N}$. $\mathbf U$ is unitary: $\mathbf U^H \mathbf U = \mathbf I_N$. Its column vectors form an orthonormal basis of $\mathbb{C}^N$—this is the **DFT spatial basis**.

For the transmitter and receiver we construct $\mathbf U_t \in \mathbb{C}^{N_t \times N_t}$ and $\mathbf U_r \in \mathbb{C}^{N_r \times N_r}$ respectively.

**4) Key lemma: DFT basis expansion of an arbitrary steering vector**

This is the central bridge connecting continuous angles and the discrete basis. For **any continuous** $\Omega$ (not required to fall exactly on the grid), its steering vector $\mathbf a(\Omega)$ can be **exactly represented** as a linear combination of the DFT basis:

$$
\boxed{
\mathbf a(\Omega) = \mathbf U\ \mathbf d(\Omega)
}
$$

where $\mathbf d(\Omega) \in \mathbb{C}^N$ is the **DFT coefficient vector**, satisfying $\mathbf d(\Omega) = \mathbf U^H \mathbf a(\Omega)$. Its $i$-th component is specifically:

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

This coefficient is a **Dirichlet kernel** (periodic sinc function) with argument $\Delta\Omega = \Omega - i/N$:

- When $\Omega$ **exactly coincides** with the $i$-th grid point ($\Delta\Omega = 0$): $d_i = 1$, all other $d_j = 0$—the steering vector reduces to a single column of $\mathbf U$;
- When $\Omega$ **does not lie on any grid point** (off-grid): the magnitude of $d_i(\Omega)$ decays as $\left|\frac{\sin(\pi N\Delta\Omega)}{N\sin(\pi\Delta\Omega)}\right|$—energy **leaks** into multiple neighboring DFT basis vectors.

This answers the core question of this subsection: **the DFT basis is a complete orthonormal basis of $\mathbb{C}^N$, and the steering vector at any continuous angle can be exactly expanded as a linear combination of $N$ basis vectors.** When the angle is exactly aligned with the grid, the expansion requires only one term (strictly sparse); when misaligned, it requires a multi-term superposition (approximately sparse). This is the essence of angular-domain sparse representation.

**5) DFT basis expansion of the geometric multipath channel**

Substituting the above expansion into the geometric multipath channel of §3.3.1:

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

**The channel matrix can be fully represented in terms of DFT basis matrices and the corresponding discrete coefficient matrices.**

**6) Virtual channel representation**

Define the **angular-domain virtual channel matrix**:

$$
\boxed{
\mathbf H^a[m,k] \triangleq \mathbf U_r^H\ \mathbf H[m,k]\ \mathbf U_t
}
$$

This is a unitary transformation: $\|\mathbf H^a\|_F = \|\mathbf H\|_F$, and the two share the same mutual information, condition number, etc. The physical meaning of the $(i,j)$-th element of $\mathbf H^a$ is:

$$
[\mathbf H^a]_{i,j}
= \mathbf a_r^H\!\left(\frac{i}{N_r}\right)\ \mathbf H[m,k]\ \mathbf a_t\!\left(\frac{j}{N_t}\right)
$$

i.e., the composite channel gain from the **$j$-th transmit angular bin** ($\Omega_t = j/N_t$) to the **$i$-th receive angular bin** ($\Omega_r = i/N_r$).

From the expansion, $\mathbf H^a$ is the weighted superposition of $\mathbf d_r(\Omega_{r,\ell})\mathbf d_t^H(\Omega_{t,\ell})$ from each path. Each physical path produces an **energy blob** in $\mathbf H^a$ centered at $(i_\ell^*, j_\ell^*)$, shaped by the Dirichlet kernel.

**Vectorized form**:

$$
\mathbf h^a = \operatorname{vec}(\mathbf H^a)
= \bigl(\mathbf U_t^T \otimes \mathbf U_r^H\bigr)\,
\operatorname{vec}(\mathbf H)
$$

**7) On-grid vs. off-grid sparsity comparison**

| Case | AoA/AoD position | Structure of $\mathbf H^a$ |
|------|------------------|---------------------------|
| Strictly on-grid | $\Omega_{t,\ell} \in \bigl\{\frac{j}{N_t}\bigr\},\ \Omega_{r,\ell} \in \bigl\{\frac{i}{N_r}\bigr\}$ | Strictly sparse: only $L$ nonzero entries |
| Off-grid (realistic) | $\Omega$ not equal to any $i/N$ | Approximately sparse: each path's energy leaks into several neighboring bins according to the Dirichlet kernel |

Under the on-grid assumption, each path contributes only one nonzero entry:

$$
\mathbf H^a = \sum_{\ell=0}^{L-1} \alpha_\ell\, e^{-j2\pi k\Delta f\tau_\ell}\,
\mathbf e_{i_\ell} \mathbf e_{j_\ell}^T
$$

where $\mathbf e_i$ is the $i$-th standard basis vector. This is the ideal working model for compressive-sensing channel estimation.

Real-world systems are always off-grid, but:
- Large arrays (large $N$) → higher angular grid resolution → leakage is concentrated in fewer bins → approximate sparsity is stronger;
- In estimation algorithms, off-grid bias can be overcome via grid refinement or gridless methods (e.g., atomic norm minimization).

### 3.3.3 Joint Angle-Delay Domain Sparsity

OFDM systems' subcarrier dimension provides **delay resolution capability**. Applying an IDFT along the subcarrier dimension of $\mathbf H^a[m,k]$ yields the **angle-delay domain channel**:

$$
\boxed{
\mathbf H^a[m,n]
= \frac{1}{N_c}\sum_{k=0}^{N_c-1} \mathbf H^a[m,k]\ e^{j2\pi kn/N_c}
}
$$

or equivalently, via a joint space-time-frequency three-dimensional transform:

$$
\mathbf H^a[m,n] = \mathbf U_r^H\ \mathbf H[m,n]\ \mathbf U_t
$$

where $\mathbf H[m,n]$ is the MIMO CIR matrix at delay tap $n$.

In the angle-delay domain:
- Each triplet $(\text{AoA bin},\ \text{AoD bin},\ \text{delay tap})$ corresponds to a candidate propagation path;
- The number of physical paths $L$ is far smaller than the total number of angle-delay domain grid points $N_r \times N_t \times N_{\mathrm{tap}}$;
- Therefore $\mathbf H^a[m,n]$ exhibits **significant sparsity** in this three-dimensional space.

Stacking all dimensions into a super-vector:

$$
\mathbf h_{\mathrm{ad}} = \operatorname{vec}\bigl\{\mathbf H^a[m,n]\bigr\}
$$

The sparse support set of $\mathbf h_{\mathrm{ad}}$ directly corresponds to the angle and delay parameters of physical paths. This property is the theoretical foundation for the following channel estimation methods:

- **Compressive Sensing**: exploiting the sparsity prior to recover a high-dimensional channel from far fewer pilot resources than grid points;
- **OMP / SOMP**: greedily estimating AoA/AoD/delay/gain path by path;
- **LASSO / Sparse Bayesian Learning**: promoting sparse solutions via $\ell_1$ regularization or hierarchical Bayesian priors;
- **Deep Learning**: using neural networks to implicitly learn the sparse structure of the angle-delay domain.

### 3.3.4 Extension: UPA (Uniform Planar Array) Angular Domain Representation

Real MIMO systems (especially 5G NR massive MIMO) often employ UPA (Uniform Planar Array) rather than ULA. This section generalizes the DFT basis analysis from ULA to UPA; the core idea is: **the steering vectors and DFT basis matrices of a UPA are the Kronecker products of the corresponding ULA parts**.

**1) UPA array geometry and steering vector**

Consider an $N_h \times N_v$ UPA ($N_h$ horizontal elements, $N_v$ vertical elements), total $N = N_h N_v$ elements. Antenna elements are indexed by $(p, q)$, where $p = 0, \dots, N_h-1$ is the horizontal index and $q = 0, \dots, N_v-1$ the vertical index.

In 3D propagation, a plane wave is described by two angular parameters:
- $\theta$: elevation angle (relative to the array normal, the $\mathbf z$ axis), $0 \le \theta \le \pi$;
- $\phi$: azimuth angle (rotation in the $\mathbf{xy}$ plane), $0 \le \phi < 2\pi$.

The horizontal and vertical wavenumber components are:

$$
k_x = \frac{2\pi}{\lambda} \sin\theta \cos\phi,\qquad
k_y = \frac{2\pi}{\lambda} \sin\theta \sin\phi
$$

The phase difference of element $(p, q)$ relative to the reference element $(0, 0)$ is $p \cdot k_x d_h + q \cdot k_y d_v$, where $d_h, d_v$ are the horizontal and vertical element spacings. Define two direction cosines:

$$
\boxed{
\Omega_h \triangleq \frac{d_h}{\lambda} \sin\theta \cos\phi,\qquad
\Omega_v \triangleq \frac{d_v}{\lambda} \sin\theta \sin\phi
}
$$

Then the $(p, q)$ component of the steering vector is $e^{-j2\pi(p\Omega_h + q\Omega_v)}$. Stacking all $N_h N_v$ elements column-major into a vector of dimension $N \times 1$, the UPA steering vector can be written as:

$$
\boxed{
\mathbf a_{\mathrm{UPA}}(\Omega_h, \Omega_v)
=
\mathbf a_h(\Omega_h) \otimes \mathbf a_v(\Omega_v)
}
$$

where:
- $\mathbf a_h(\Omega_h) \in \mathbb{C}^{N_h}$ is the horizontal ULA steering vector: $[\mathbf a_h]_p = \frac{1}{\sqrt{N_h}} e^{-j2\pi p\Omega_h}$;
- $\mathbf a_v(\Omega_v) \in \mathbb{C}^{N_v}$ is the vertical ULA steering vector: $[\mathbf a_v]_q = \frac{1}{\sqrt{N_v}} e^{-j2\pi q\Omega_v}$;
- The normalization factor $\frac{1}{\sqrt{N_h N_v}}$ follows from the combination of $\otimes$ products.

**Why is UPA needed?** A ULA can steer beams only in the azimuth direction (one dimension); when communication systems need to **simultaneously steer beams in both elevation and azimuth dimensions** (e.g., high-rise building coverage, UAV communications, satellite communications), a ULA is insufficient. A UPA distributes antennas on a two-dimensional plane and can simultaneously control beam directions in both horizontal and vertical dimensions, achieving **full-space 3D beam steering**. This is the fundamental reason 5G NR massive MIMO extensively adopts UPA.

**1b) 2D beam pattern of a UPA**

Thanks to the Kronecker product structure of the steering vector, the UPA beam pattern naturally factorizes into the product of horizontal and vertical 1D beam patterns. When beamforming with the matched-filter vector $\mathbf a_{\mathrm{UPA}}(\Omega_{h,0}, \Omega_{v,0})$, the gain in direction $(\Omega_h, \Omega_v)$ is

$$
\boxed{
\begin{aligned}
G_{\mathrm{UPA}}\bigl((\Omega_{h,0}, \Omega_{v,0}),\ (\Omega_h, \Omega_v)\bigr)
&= \bigl|\mathbf a_{\mathrm{UPA}}^H(\Omega_{h,0}, \Omega_{v,0})\ \mathbf a_{\mathrm{UPA}}(\Omega_h, \Omega_v)\bigr|^2 \\[4pt]
&= \bigl|\mathbf a_h^H(\Omega_{h,0})\,\mathbf a_h(\Omega_h)\bigr|^2\ \cdot\
\bigl|\mathbf a_v^H(\Omega_{v,0})\,\mathbf a_v(\Omega_v)\bigr|^2 \\[4pt]
&= G_h(\Omega_{h,0}, \Omega_h)\ \cdot\ G_v(\Omega_{v,0}, \Omega_v),
\end{aligned}
}
$$

where $G_h$ and $G_v$ are the 1D Dirichlet-kernel-squared beam patterns defined in the previous section (with $N_h$ and $N_v$ antennas, respectively). This **separability** is the core convenience of UPA angular-domain analysis: the 2D beam is the tensor product of two 1D beams, and all conclusions regarding mainlobe width, null positions, and off-grid leakage carry over dimension by dimension.

Key features of the 2D beam pattern:

- **Mainlobe**: located at $(\Omega_h, \Omega_v) = (\Omega_{h,0}, \Omega_{v,0})$, peak gain $1$;
- **Null grid**: nulls appear at $\Omega_h = \Omega_{h,0} \pm k/N_h$ or $\Omega_v = \Omega_{v,0} \pm \ell/N_v$ ($k = 1,\dots,N_h-1$, $\ell = 1,\dots,N_v-1$), forming a regular "null lattice";
- **Sidelobes**: first sidelobe level approximately $-13$ dB (same as 1D, determined by the Dirichlet kernel sidelobe structure);
- **Angular resolution**: horizontal resolution $\propto 1/N_h$, vertical resolution $\propto 1/N_v$. Thus a $32 \times 32$ UPA not only equals a $1024$-element ULA in total element count but also provides high resolution simultaneously in both dimensions.

**2) UPA DFT basis matrix**

For the horizontal dimension: uniformly discretize $\Omega_h$ into $N_h$ grid points $\Omega_{h,i} = i/N_h\ (i = 0, \dots, N_h-1)$, and construct the $N_h \times N_h$ unitary DFT matrix $\mathbf U_h$.

For the vertical dimension: uniformly discretize $\Omega_v$ into $N_v$ grid points $\Omega_{v,j} = j/N_v\ (j = 0, \dots, N_v-1)$, and construct the $N_v \times N_v$ unitary DFT matrix $\mathbf U_v$.

The **UPA DFT basis matrix** is the Kronecker product of the two dimensions:

$$
\boxed{
\mathbf U_{\mathrm{UPA}} = \mathbf U_h \otimes \mathbf U_v \ \in \mathbb{C}^{N_h N_v \times N_h N_v}
}
$$

This is essentially a **2D-DFT matrix**:

- Each column of $\mathbf U_{\mathrm{UPA}}$ corresponds to a pair of discrete angular bins $(\Omega_{h,i}, \Omega_{v,j})$, i.e., one grid point of the 2D angular domain;
- The two dimensions together have $N_h \times N_v$ angular bins, matching the total number of array elements.

**3) UPA virtual channel representation**

Both the transmitter and receiver of a MIMO system can be UPAs. Suppose the transmitter is a $N_{t,h} \times N_{t,v}$ UPA and the receiver a $N_{r,h} \times N_{r,v}$ UPA. The channel matrix is $\mathbf H \in \mathbb{C}^{N_r N_{r,h} N_{r,v} \times N_t N_{t,h} N_{t,v}}$ (assuming flat fading or a single subcarrier). The angular-domain virtual channel is:

$$
\boxed{
\mathbf H^a = \bigl(\mathbf U_{r,h} \otimes \mathbf U_{r,v}\bigr)^H\ \mathbf H\ \bigl(\mathbf U_{t,h} \otimes \mathbf U_{t,v}\bigr)
}
$$

Using the Kronecker product property $(\mathbf A \otimes \mathbf B)^H = \mathbf A^H \otimes \mathbf B^H$:

$$
\mathbf H^a = \bigl(\mathbf U_{r,h}^H \otimes \mathbf U_{r,v}^H\bigr)\ \mathbf H\ \bigl(\mathbf U_{t,h} \otimes \mathbf U_{t,v}\bigr)
$$

The dimension of $\mathbf H^a$ is likewise $N_{r,h}N_{r,v} \times N_{t,h}N_{t,v}$, and each element $(I, J)$ corresponds to a set of 4D angular bins:

$$
(\Omega_{r,h,i}, \Omega_{r,v,j}) \to (\Omega_{t,h,i'}, \Omega_{t,v,j'})
$$

**4) On-grid 2D sparsity**

If the AoA/AoD of all physical paths happen to fall on the 2D angular grid, $\mathbf H^a$ will have only $L$ nonzero entries—each path contributes a complex gain in one 4D angular bin. The number of paths $L$ is far smaller than the total number of angular bins $N_{r,h}N_{r,v}N_{t,h}N_{t,v}$; hence $\mathbf H^a$ is a **strictly sparse 2D channel tensor**.

**5) Comparison with ULA**

| Property | ULA | UPA |
|----------|-----|-----|
| Angular parameters | 1D: $\Omega$ ($\sin\theta$) | 2D: $(\Omega_h, \Omega_v)$ = $(\sin\theta\cos\phi,\ \sin\theta\sin\phi)$ |
| Steering vector | $\mathbf a(\Omega)$ | $\mathbf a_h(\Omega_h) \otimes \mathbf a_v(\Omega_v)$ |
| DFT basis matrix | $\mathbf U$ (1D-DFT) | $\mathbf U_h \otimes \mathbf U_v$ (2D-DFT) |
| Total angular bins | $N$ | $N_h \times N_v$ |
| Beamforming dimension | Azimuth only | Azimuth + Elevation (full-space 3D) |
| Beam pattern structure | 1D Dirichlet-kernel-squared | Product of two 1D Dirichlet-kernel-squares |
| Sparse domain schematic | 1D angular spectrum | 2D angular map (elevation-azimuth) |
| Typical application | 3G/4G/Wi-Fi | 5G NR massive MIMO |

**6) Tensor perspective — Einstein summation convention**

UPA channels possess a natural **multi-index structure**: transmit antennas are jointly labeled by a horizontal index $h$ and a vertical index $v$, and similarly for the receiver. Rather than forcibly flattening all antennas into a 1D vector (losing row-column correspondence), it is more natural to use the language of **tensors**.

**Upper/lower index notation**: adopt the Einstein summation convention—a repeated upper-lower index pair implies automatic summation.

Denote the transmit signal as a 2-index tensor $X^{hv}$ (complex signal at column $h$, row $v$) and the receive signal as $Y^{h'v'}$. Then the spatial MIMO input-output relation can be written in 4-index form:

$$
Y^{h'v'} = H^{h'v'}_{\ \ \ hv}\ X^{hv},
$$

where:
- $H^{h'v'}_{\ \ \ hv}$ is the complex channel gain from transmit antenna $(h, v)$ to receive antenna $(h', v')$;
- By the Einstein convention, the repeated superscripts $hv$ are automatically summed over $h = 0,\dots,N_{t,h}-1$ and $v = 0,\dots,N_{t,v}-1$;
- $H^{h'v'}_{\ \ \ hv}$ can be viewed either as an $N_{r,h}N_{r,v} \times N_{t,h}N_{t,v}$ matrix, or as a four-dimensional array of size $N_{r,h} \times N_{r,v} \times N_{t,h} \times N_{t,v}$.

**Tensor form of spatial signatures**: the horizontal spatial signature $\mathbf e_h(\Omega_h)$ and the vertical spatial signature $\mathbf e_v(\Omega_v)$ combine via tensor product to form the UPA spatial signature:

$$
e^{hv}(\Omega_h, \Omega_v) = e_h^h(\Omega_h)\ e_v^v(\Omega_v),
$$

where $e^{hv}$ denotes the spatial signature value at a specific antenna $(h, v)$ (i.e., the $(h, v)$-th component of $\mathbf a_{\mathrm{UPA}}$).

The inner product of two spatial signatures is remarkably concise in Einstein notation—conjugate transpose corresponds to flipping upper/lower indices:

$$
\langle \mathbf a_1, \mathbf a_2 \rangle = \bar e_1^{hv}\ e_2^{hv},
$$

where $\bar e_1^{hv}$ already incorporates conjugation (implied by the upper-index-flip notation).

**Tensor transform to the angular domain**: the spatial-domain vector $X^{hv}$ transforms to the angular-domain vector $x^{ij}$ via a 2D-IDFT (equally elegant in tensor notation):

$$
x^{ij} = X^{hv}\ e_h^h(\Omega_{h,i})\ e_v^v(\Omega_{v,j}),
$$

where $\Omega_{h,i} = i/N_h$, $\Omega_{v,j} = j/N_v$. This is precisely expanding the spatial signal separately along the horizontal and vertical DFT bases.

Correspondingly, the angular-domain channel tensor is:

$$
\boxed{
h^{i'j'}_{\ \ \ ij} = \bar e_{r,h}^{i'}\ \bar e_{r,v}^{j'}\ H^{h'v'}_{\ \ \ hv}\ e_{t,h}^{i}\ e_{t,v}^{j}
}
$$

The physical meaning of $h^{i'j'}_{\ \ \ ij}$ is: the equivalent complex gain of the path emitted from transmit horizontal angular bin $i$, vertical angular bin $j$, and captured at receive horizontal angular bin $i'$, vertical angular bin $j'$. Expanded, this is a 2D-DFT:

$$
h^{i'j'}_{\ \ \ ij} = \frac{1}{N_{t,h}N_{t,v}N_{r,h}N_{r,v}}
\sum_{h,v,h',v'}
H^{h'v'}_{\ \ \ hv}\
e^{j2\pi(h'i'/N_{r,h} + v'j'/N_{r,v})}\
e^{-j2\pi(hi/N_{t,h} + vj/N_{t,v})}.
$$

**Convenience of tensor notation**: the full MIMO spatial-to-angular domain transform requires only one line:

$$
Y^{h'v'} = H^{h'v'}_{\ \ \ hv}\ X^{hv}
\quad\longleftrightarrow\quad
y^{i'j'} = h^{i'j'}_{\ \ \ ij}\ x^{ij}.
$$

This notation naturally reflects the separable structure of the UPA—the horizontal and vertical dimensions evolve independently in the indices, without the need to artificially flatten and reorganize matrices. It is the natural language for all Kronecker product decompositions.

The Kronecker structure of the UPA allows its angular-domain analysis to be completely decomposed into independent horizontal and vertical dimensions; the UPA 2D-DFT basis remains a unitary matrix. All conclusions from §3.3.2 concerning Dirichlet kernel expansion, off-grid leakage, etc., carry over to UPA dimension by dimension.

## 3.4 3GPP TR 38.901 CDL/TDL Models

3GPP TR 38.901 is the standardized channel model for 5G/5G-Advanced link-level and system-level simulation, covering frequencies from 0.5 to 100 GHz and supporting UE speeds from 0 to 1000 km/h. It integrates multiple previous 3GPP channel models (3D-UMa/UMi SCM, IMT-Advanced, etc.) and is updated in essentially every specification release with ongoing extensions.

### 3.4.1 Model Architecture

TR 38.901 divides the channel into two layers:

| Layer | Content | Update Rate |
|-------|---------|-------------|
| **Large-Scale** | Path loss, shadow fading, penetration loss, LOS/NLOS state determination | Per drop (typically tens of ms) |
| **Small-Scale** | Clusters, rays, delays, powers, angular parameters, polarization, Doppler | Per snapshot (OFDM symbol level) |

The simulation workflow can be summarized in three steps:

1. **Drop Initialization**: determine large-scale parameters according to the scenario (UMa/UMi/RMa/InH, etc.), and decide LOS/NLOS;
2. **Generate small-scale parameters**: look up the number of clusters $N$ and rays per cluster $Q$ from tables, and generate cluster delays, powers, arrival angles (AoA/ZoA), departure angles (AoD/ZoD), cross-polarization ratio (XPR), etc.;
3. **Time stepping**: at each sampling instant, compute the Doppler phase from ray angles and UE velocity, multiply by antenna patterns and array phases, and superpose all clusters and rays.

### 3.4.2 CDL Channel Impulse Response

CDL (Clustered Delay Line) explicitly models the angular parameters of each ray and can therefore accurately capture the spatial characteristics of MIMO arrays. The contribution of the $q$-th ray ($q = 1,\dots,Q$) of the $n$-th cluster ($n = 1,\dots,N$) from transmit antenna $s$ to receive antenna $u$ is:

$$
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
$$

For LOS scenarios, a direct-path component must be added:

$$
h_{u,s}^{\mathrm{LOS}}(t) =
\sqrt{\frac{1}{K_R+1}}\,
\mathbf F_{r,u}^T(\theta_{\mathrm{LOS,ZoA}},\phi_{\mathrm{LOS,AoA}})
\begin{bmatrix} 1 & 0 \\ 0 & -1 \end{bmatrix}
\mathbf F_{t,s}(\theta_{\mathrm{LOS,ZoD}},\phi_{\mathrm{LOS,AoD}}) \\
\times
\exp\!\Bigl(j\frac{2\pi}{\lambda_0}\,\hat{\mathbf r}_{r,\mathrm{LOS}}^T\mathbf d_{r,u}\Bigr)
\exp\!\Bigl(j\frac{2\pi}{\lambda_0}\,\hat{\mathbf r}_{t,\mathrm{LOS}}^T\mathbf d_{t,s}\Bigr)
\exp\!\Bigl(j\frac{2\pi}{\lambda_0}\,\hat{\mathbf r}_{r,\mathrm{LOS}}^T\mathbf v\,t\Bigr),
$$

The complete impulse response is the sum of the LOS component and the small-scale NLOS components; $K_R$ is the Rician K-factor (linear).

### 3.4.3 Parameter-by-Parameter Explanation

The following table explains each symbol in the CDL formula, its physical meaning, and its range of values:

| Symbol | Meaning | Remarks |
|--------|---------|---------|
| $N$ | Number of clusters | CDL-A/B/C contain 23–24 clusters; CDL-D/E contain 13–14 clusters (including LOS cluster) |
| $Q$ | Rays per cluster | Standard value $Q = 20$. These 20 rays are distributed within the cluster around the cluster center angle over the range $\pm c_{\mathrm{ASD/ASA/ZSD/ZSA}}$ |
| $P_n$ | Normalized power of the $n$-th cluster | $\sum_n P_n = 1$; the first cluster has the largest power, with exponential decay thereafter. Under LOS the LOS cluster power is separately allocated according to the K-factor |
| $\tau_n$ | Delay of the $n$-th cluster | Given by normalized delay tables; scaled by the desired delay spread (DS) in simulation |
| $\theta, \phi$ | Elevation angle (ZoD/ZoA), Azimuth angle (AoD/AoA) | Zo = Zenith, A = Azimuth. The four angles are each composed of the cluster center angle plus a ray offset |
| $c_{\mathrm{ASD}},\ c_{\mathrm{ASA}}$ | Azimuth spread factors for departure/arrival | Control the spread of rays within a cluster in azimuth; predefined by the standard |
| $c_{\mathrm{ZSD}},\ c_{\mathrm{ZSA}}$ | Zenith spread factors for departure/arrival | Same as above, for the elevation dimension |
| $\mathbf F_{t,s},\ \mathbf F_{r,u}$ | Complex field patterns of transmit antenna $s$ / receive antenna $u$ | $2 \times 1$ or $1 \times 2$ vectors containing polarization ($\theta$ and $\phi$ components). If polarization antennas are not considered, may be taken as $[1, 1]^T$ or scalar $1$ |
| $\mathbf \Phi_{n,q}$ | Polarization coupling matrix for ray $q$ of cluster $n$ | $2 \times 2$ matrix containing the cross-polarization ratio $\kappa$ (XPR), describing power transfer between $\theta/\phi$ components |
| $\lambda_0$ | Carrier wavelength | $\lambda_0 = c / f_c$ |
| $\hat{\mathbf r}_{t,n,q}$ | Unit vector of the transmit direction | $[\sin\theta\cos\phi,\ \sin\theta\sin\phi,\ \cos\theta]^T$, determined by ZoD and AoD |
| $\hat{\mathbf r}_{r,n,q}$ | Unit vector of the receive direction | Same as above, determined by ZoA and AoA |
| $\mathbf d_{t,s}$ | Position vector of transmit antenna $s$ | Relative to the array reference point (typically the origin), in the antenna array coordinate system |
| $\mathbf d_{r,u}$ | Position vector of receive antenna $u$ | Same as above |
| $\mathbf v$ | UE velocity vector | $\mathbf v = v \cdot [\sin\theta_v\cos\phi_v,\ \sin\theta_v\sin\phi_v,\ \cos\theta_v]^T$, where $(\theta_v, \phi_v)$ are the velocity direction angles |
| $K_R$ | Rician K-factor (linear) | Used only in LOS scenarios; $K_R = 0$ for NLOS (i.e., no LOS term) |

**Core insight**: the CDL formula is essentially a **standardized implementation** of the geometric multipath model of §3.3—it likewise superposes paths using steering vectors (array phase terms), but all path parameters (angles, delays, powers, XPR) come from the predefined tables of TR 38.901 rather than being user-specified.

### 3.4.4 Five Predefined CDL Models

Tables 7.7.1-1 through 7.7.1-5 of TR 38.901 define five CDL models, covering typical scenarios from scattering-rich NLOS to strong LOS. Key parameters are given below; for the full per-cluster delay/power/angle tables, please refer to the specification.

| Parameter | CDL-A | CDL-B | CDL-C | CDL-D | CDL-E |
|-----------|-------|-------|-------|-------|-------|
| **Propagation condition** | NLOS | NLOS | NLOS | LOS | LOS |
| **Clusters $N$** | 23 | 23 | 24 | 13 | 14 |
| **Rays per cluster $Q$** | 20 | 20 | 20 | 20 | 20 |
| **Normalized delay spread** | ~0.44 | ~1.07 | ~2.46 | ~0.29 | ~0.30 |
| **Angular spread** | **Large** | **Medium** | **Medium–large** | **Smaller** | **Smallest** |
| ASD (azimuth spread of departure) / ° | ~22–26 | ~10–11 | ~15–18 | ~4–6 | ~3–4 |
| ASA (azimuth spread of arrival) / ° | ~37–44 | ~19–22 | ~26–29 | ~7–9 | ~5–7 |
| ZSD (zenith spread of departure) / ° | ~6 | ~6 | ~6 | ~3 | ~3 |
| ZSA (zenith spread of arrival) / ° | ~8 | ~8 | ~8 | ~3–4 | ~3–4 |
| **Typical K-factor / dB** | — | — | — | ~11 | ~20 |
| **XPR (cross-polarization ratio) / dB** | 11 | 9 | 7 | 11 | 8 |
| **Typical scenario** | Urban macrocell NLOS | Urban microcell NLOS | Suburban/rural NLOS | Urban macrocell LOS | Suburban/open-area LOS |

**How to interpret**:

- **CDL-A**: largest angular spread, richest scattering; the MIMO channel matrix tends toward full rank, with high spatial multiplexing potential;
- **CDL-B**: medium angular spread; the default NLOS baseline in many simulations;
- **CDL-C**: largest delay spread (longer multipath); strongest frequency selectivity, requiring higher OFDM pilot density;
- **CDL-D**: weak LOS ($K_R \approx 10$–$11$ dB), small angular spread; the channel matrix tends toward low rank, and beamforming gain is significant;
- **CDL-E**: strong LOS ($K_R \approx 20$ dB), very small angular spread; the channel approaches rank 1, and the number of spatial multiplexing layers is far fewer than the number of antennas.

### 3.4.5 Delay and Angular Scaling

The normalized parameters in the tables need to be scaled according to the desired Delay Spread (DS) and Angular Spread (AS) to match a specific scenario:

- **Delay scaling**: $\tau_{n,\mathrm{scaled}} = \tau_{n,\mathrm{norm}} \cdot \mathrm{DS}_{\mathrm{desired}}$, where $\mathrm{DS}_{\mathrm{desired}}$ is a large-scale parameter (e.g., typical DS = 129 ns for UMa NLOS);
- **Angular scaling**: the four angles (ASD, ASA, ZSD, ZSA) each scale the intra-cluster ray offsets $c_{\mathrm{ASD}}$, etc., by the desired values, ensuring that the final angular spread matches the scenario.

This "normalized table + per-scenario scaling" design allows the same set of CDL models to be adapted to diverse scenarios such as UMa, UMi, RMa, and InH, simply by changing the large-scale parameters rather than replacing the entire model.

### 3.4.6 TDL Model

TDL (Tapped Delay Line) is a simplified version of CDL. It **discards all angular and polarization information**, retaining only delays and average powers:

| Model | Propagation condition | Characteristics |
|-------|----------------------|----------------|
| TDL-A | NLOS | 23 taps, large delay spread |
| TDL-B | NLOS | 23 taps, medium |
| TDL-C | NLOS | 24 taps, largest delay spread |
| TDL-D | LOS | 13 taps + LOS tap |
| TDL-E | LOS | 14 taps + LOS tap, strong LOS |

TDL is suitable for the following scenarios:
- SISO / SIMO simulations, where antenna array geometry is not needed;
- When only frequency selectivity matters, not spatial selectivity;
- Rapid prototyping, omitting antenna pattern and array phase computations.

Each TDL tap independently undergoes Rayleigh (NLOS tap) or Rician (LOS tap) fading, essentially reducing to the independent tap model of §3.1, but with parameters taken from standardized tables.

### 3.4.7 LOS/NLOS State Determination

In TR 38.901, the LOS probability is a function of the UE–BS distance, with different empirical formulas for different scenarios. Taking UMa (Urban Macro) as an example:

$$
P_{\mathrm{LOS}}(d_{2D}) =
\begin{cases}
1, & d_{2D} \le 18\ \mathrm{m}, \\[4pt]
\displaystyle\Bigl(\frac{18}{d_{2D}} + \exp\!\bigl(-\frac{d_{2D}}{63}\bigr)\bigl(1 - \frac{18}{d_{2D}}\bigr)\Bigr), & d_{2D} > 18\ \mathrm{m}.
\end{cases}
$$

During simulation, each drop draws from $P_{\mathrm{LOS}}$ once, deciding whether to use CDL-D/E (LOS) or CDL-A/B/C (NLOS).

### 3.4.8 Comparison with Models of Previous Sections

| Feature | §3.1 i.i.d. tap | §3.2 Kronecker | §3.3 Geometric AoA/AoD | §3.4 CDL/TDL |
|---------|-----------------|----------------|------------------------|--------------|
| Angular information | None | Statistical (correlation matrix) | Arbitrary (user-specified) | Standardized tables |
| Polarization | None | None | Optional | Yes (XPR + polarized field patterns) |
| Antenna patterns | None | None | Optional | Yes ($\mathbf F_{t,s}, \mathbf F_{r,u}$) |
| Standardization level | None | Low | None | High (3GPP specification) |
| Simulation complexity | Lowest | Low–Medium | Medium | Medium–High |
| Applicable scenarios | Algorithm validation | Correlation analysis | Beam/sparsity studies | Standards-aligned simulation |

**Selection advice**: if the goal is to align with 5G NR link-level evaluation, or if complete antenna patterns, polarization effects, and 3D spatial information are required, use CDL; if only rapid validation of channel estimation or precoding algorithm correctness is needed, the simplified models of §3.1–§3.3 suffice.

> The content on pilot design and channel estimation has been moved to [MIMOOFDMChannelEstimation.md](MIMOOFDMChannelEstimation.md).

# 4 Simulation Workflow

A reproducible MIMO-OFDM channel simulation is typically organized according to the following workflow.

### Step 1: Set System Parameters

Given:

$$
N_t,N_r,N,M,\Delta f,T_{cp},f_c,v_{\max}.
$$

Compute:

$$
T_s=\frac{1}{N\Delta f},\qquad
T_{\mathrm{sym}}=\frac{1}{\Delta f}+T_{cp},\qquad
f_D=\frac{v_{\max}}{\lambda}.
$$

### Step 2: Set PDP and Taps

Given:

$$
\{P_\ell,\tau_\ell\}_{\ell=0}^{L-1}.
$$

If using discrete taps, let:

$$
n_\ell=\mathrm{round}(\tau_\ell/T_s).
$$

### Step 3: Generate Spatial Channel

Available models:

1. i.i.d. Rayleigh/Rician;
2. Kronecker spatial correlation;
3. Geometric AoA/AoD model;
4. 3GPP CDL/TDL.

Using Kronecker as an example:

$$
\mathbf G_\ell[m]=
\mathbf R_{r,\ell}^{1/2}\mathbf W_\ell[m]\mathbf R_{t,\ell}^{1/2}.
$$

### Step 4: Generate Time Correlation

Options include:

- Jakes/Clarke sum-of-sinusoids;
- AR(1) approximation;
- BEM;
- 3GPP Doppler phase from ray directions and velocity.

### Step 5: Construct CFR

$$
\mathbf H[m,k]=
\sum_{\ell=0}^{L-1}
\sqrt{P_\ell}\mathbf G_\ell[m]
e^{-j2\pi k\Delta f\tau_\ell}.
$$

### Step 6: Transmit Pilots and Data

For each $(m,k)$:

$$
\mathbf y[m,k]=\mathbf H[m,k]\mathbf x[m,k]+\mathbf w[m,k].
$$

### Step 7: Estimation, Interpolation, and Equalization

1. Perform LS or LMMSE at pilot positions;
2. Interpolate across time-frequency-space dimensions;
3. Perform ZF/MMSE/ML/MAP detection on each subcarrier;
4. Compute statistics: NMSE, BER, EVM, capacity, or throughput.

Common NMSE:

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

# 5 Model Selection Recommendations

| Goal | Recommended Model | Notes |
|------|-------------------|-------|
| Verify LS/MMSE formulas | i.i.d. Rayleigh tap | Simple, controllable, easy to reproduce |
| Study spatial correlation effects | Kronecker / Weichselberger | Correlation matrices can be explicitly adjusted |
| Study arrays, angles, beams | Geometric AoA/AoD model | Can explain rank, angular sparsity, and beam training |
| Align with 5G NR link-level simulation | 3GPP TR 38.901 CDL/TDL | Standardized, complete parameters |
| Massive MIMO / mmWave estimation | Angle-delay sparse model | Suitable for compressive sensing and deep learning estimation |
| High-mobility doubly selective | Jakes/AR/BEM + MIMO spatial model | Requires explicit time correlation |

Empirical guidelines:

- If only tuning a channel estimation algorithm, start with i.i.d. Rayleigh;
- If validating MIMO-specific phenomena, at least add spatial correlation or geometric angles;
- If writing simulations close to communication standards, use 3GPP CDL;
- If the array size is very large, do not estimate solely in the antenna domain; consider angular-domain or low-rank/sparse priors.

# 6 References

1. D. Tse and P. Viswanath, *Fundamentals of Wireless Communication*, Cambridge University Press, 2005. Chapter 7: MIMO I, Spatial Multiplexing and Channel Modeling. https://web.stanford.edu/~dntse/wireless_book.html
2. A. M. Sayeed, "Deconstructing Multiantenna Fading Channels," *IEEE Transactions on Signal Processing*, vol. 50, no. 10, pp. 2563-2579, 2002. https://minds.wisconsin.edu/handle/1793/9386
3. D. Shiu, G. J. Foschini, M. J. Gans, and J. M. Kahn, "Fading correlation and its effect on the capacity of multielement antenna systems," *IEEE Transactions on Communications*, vol. 48, no. 3, pp. 502-513, 2000.
4. 3GPP TR 38.901, *Study on channel model for frequencies from 0.5 to 100 GHz*. ETSI TR 138 901 V19.2.0, 2026-02. https://www.etsi.org/deliver/etsi_tr/138900_138999/138901/19.02.00_60/tr_138901v190200p.pdf
5. 3D MIMO channel modeling introduction. https://zhuanlan.zhihu.com/p/664298097
6. The spatial correlation of MIMO channels. https://zhuanlan.zhihu.com/p/721013608
7. 逸风晴 (Yifengqing), "MIMO Communication Angular Domain Representation," Zhihu Column, 2021. https://zhuanlan.zhihu.com/p/363018716
