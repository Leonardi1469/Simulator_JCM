# Simulator JCM

This repository contains a Streamlit application for the numerical study of the Jaynes-Cummings model. The program solves the Schrodinger equation with QuTiP for a two-level atom interacting with a single cavity mode. The initial state and model parameters can be selected through the interface, and the calculated observables can be inspected as functions of time.

The application is an extended version of the simulator developed in connection with Emily Andrea Franco Escudero's undergraduate thesis, *Design and implementation of JC Simulator: an interactive simulator for the Jaynes-Cummings model*. The version associated with the thesis is kept in a [separate repository](https://github.com/Leonardi1469/JC-Simulator).

## Physical model

In units where $\hbar=1$, the Hamiltonian is

$$
H=\omega_c a^\dagger a+\frac{\omega_{eg}}{2}\sigma_z
  +g\left(a\sigma_+ + a^\dagger\sigma_-\right).
$$

The first term describes the cavity field, the second describes the two-level atom, and the last accounts for the exchange of one excitation between them. The parameters $\omega_c$, $\omega_{eg}$, and $g$ are the cavity frequency, atomic transition frequency, and coupling constant, respectively. The operators $a$ and $a^\dagger$ annihilate and create a cavity photon, while $\sigma_+$ and $\sigma_-$ raise and lower the atomic state. The detuning is $\Delta=\omega_{eg}-\omega_c$.

The Hamiltonian is written in the rotating-wave approximation (RWA), which neglects the counter-rotating terms $a\sigma_-$ and $a^\dagger\sigma_+$. The approximation is appropriate when the coupling is small compared with $\omega_c+\omega_{eg}$ for the photon numbers involved in the dynamics. Near resonance, the detuning should also be small compared with this sum. The interface reports the ratios $|\Delta|/(\omega_c+\omega_{eg})$ and $g\sqrt{\langle n(0)\rangle+1}/(\omega_c+\omega_{eg})$ as an initial guide. A warning is displayed when either ratio is at least 0.1; this value is indicative and does not define a strict validity limit. Since the photon distribution evolves, the user should also consider the photon numbers populated at later times.

The model assumes a closed system. Cavity losses, spontaneous emission, thermal effects, and counter-rotating interactions are not included. All frequencies and $g$ must be entered in the same angular-frequency units; time is expressed in the inverse unit.

## Initial states and numerical parameters

The atom can start in $\lvert e\rangle$, $\lvert g\rangle$, or the superposition

$$
\lvert\psi_a(0)\rangle=
\cos(\theta/2)\lvert e\rangle+
e^{i\phi}\sin(\theta/2)\lvert g\rangle.
$$

The initial field can be vacuum, a Fock state $\lvert n_0\rangle$, a coherent state $\lvert\alpha\rangle$, a squeezed vacuum $S(r)\lvert0\rangle$, or a displaced squeezed vacuum $D(\alpha)S(r)\lvert0\rangle$. In the current interface, $\alpha$ and $r$ are real. The order of $D(\alpha)$ and $S(r)$ is part of the definition of the last state.

The numerical parameters are the Fock-space dimension $N$, final time $t_{\max}$, and number of sampled times $N_t$. The field basis contains $\lvert0\rangle,\ldots,\lvert N-1\rangle$. When the highest retained level acquires a noticeable population, the program displays a warning. Results should be checked by repeating the calculation with a larger $N$, particularly for coherent and squeezed states. Increasing $N_t$ provides more output points but does not, by itself, establish the accuracy of the differential-equation solver.

## Calculated quantities

The main figure contains six time-dependent quantities:

1. Excited-state probability $P_e(t)$.
2. Ground-state probability $P_g(t)$.
3. Atomic inversion $W(t)=P_e(t)-P_g(t)$.
4. Mean photon number $\langle n(t)\rangle$.
5. Von Neumann entropy of the reduced atomic state $S_A(t)$, in bits.
6. Mean total number of excitations $\langle M(t)\rangle=\langle n(t)\rangle+P_e(t)$.

The photon-number distribution $P_n(t)$ is shown during the MP4 animation and recorded for every sampled time in the CSV file. The application also calculates the field statistics

$$
Q(t)=\frac{\langle n(n-1)\rangle-\langle n\rangle^2}{\langle n\rangle},
\qquad
g^{(2)}(0,t)=\frac{\langle n(n-1)\rangle}{\langle n\rangle^2}.
$$

These quantities are undefined when their denominators vanish. The displayed $g^{(2)}(0,t)$ is an equal-time field statistic, not a two-time correlation function. For the closed Jaynes-Cummings Hamiltonian, $\langle M(t)\rangle$ is conserved; the program reports its maximum numerical variation as a consistency check. This check does not replace a convergence study in $N$.

## Running the application

Install the dependencies with Python 3.11 or newer, then start Streamlit:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
streamlit run app.py
```

In Windows Command Prompt, replace the activation command with `.venv\Scripts\activate.bat`. In PowerShell, use `.\.venv\Scripts\Activate.ps1`.

Set the parameters and initial states, then select **Run simulation**. Select **Generate MP4 video** to view the evolving inversion, atomic populations, and photon distribution. The video can be paused or advanced with the player controls, and it can also be downloaded. It contains at most 90 frames selected from the calculated time points. MP4 encoding uses the executable provided by `imageio-ffmpeg`.

The six-panel figure is available in PDF and PNG, and the sampled observables and probabilities $P_n(t)$ are available in CSV. The exports correspond to the most recent run. The PDF contains the figure rather than a report of the input parameters; the CSV records those parameters in its header.

As a reference case, choose $\omega_c=\omega_{eg}=1$, $g=0.05$, an initially excited atom, and a vacuum field. The analytical result is $P_e(t)=\cos^2(gt)$ and $\langle M(t)\rangle=1$. A sufficiently long time interval is needed to observe an oscillation.

The application code is contained in `app.py`; dependencies are listed in `requirements.txt`. This repository currently has no license.
