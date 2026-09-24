# Simulator JCM

An interactive educational application for the closed Jaynes-Cummings model, built with [Streamlit](https://streamlit.io/) and [QuTiP](https://qutip.org/). Choose a two-level atom and a single-mode cavity field, inspect their time evolution, explore individual instants with a time slider, and generate a playable MP4 animation. Figures can be downloaded for teaching or scientific writing.

This expanded research version builds on the simulator associated with Emily Andrea Franco Escudero's undergraduate thesis, *Design and implementation of JC Simulator: an interactive simulator for the Jaynes-Cummings model*. The [thesis repository](https://github.com/Leonardi1469/JC-Simulator) remains separate; this repository is for the subsequent article.

## Model and physical assumptions

In units with $\hbar=1$, the Hamiltonian is

$$
H=\omega_c a^\dagger a+\frac{\omega_{eg}}{2}\sigma_z
  +g\left(a\sigma_+ + a^\dagger\sigma_-\right).
$$

The first term is the cavity-mode energy, the second is the two-level atom energy, and the last exchanges one excitation between the atom and field. $a$ annihilates a photon, $a^\dagger$ creates one, and $\sigma_+$ and $\sigma_-$ raise and lower the atom. The detuning is $\Delta=\omega_{eg}-\omega_c$. All frequencies and $g$ must use the same angular-frequency unit; time uses its inverse.

This Hamiltonian uses the **rotating-wave approximation** (RWA): the counter-rotating terms $a\sigma_-$ and $a^\dagger\sigma_+$ are omitted. A useful near-resonant regime has $|\Delta|/(\omega_c+\omega_{eg})\ll1$ and $g\sqrt{n+1}/(\omega_c+\omega_{eg})\ll1$ for photon numbers appreciably occupied during the evolution. The app displays estimates using the initial mean photon number and warns when an illustrative ratio reaches 0.1. This threshold is a teaching aid, not a universal boundary. Outside that regime, the program still solves the JC Hamiltonian, but its interpretation as an approximation to a particular physical system requires care.

The model contains **one atom, one field mode, and no cavity losses, spontaneous emission, or thermal bath**. It does not solve the full quantum Rabi Hamiltonian or a dissipative master equation.

## Inputs and observables

| Group | Choices or meaning |
| --- | --- |
| Physical parameters | Cavity frequency $\omega_c$, atomic transition frequency $\omega_{eg}$, coupling $g$. |
| Atomic state | Excited $\lvert e\rangle$, ground $\lvert g\rangle$, or $\cos(\theta/2)\lvert e\rangle+e^{i\phi}\sin(\theta/2)\lvert g\rangle$. |
| Field state | Vacuum, Fock $\lvert n_0\rangle$, coherent $\lvert\alpha\rangle$, squeezed vacuum $S(r)\lvert 0\rangle$, or displaced squeezed vacuum $D(\alpha)S(r)\lvert 0\rangle$. The interface uses real $\alpha$ and $r$. |
| Numerical parameters | Fock dimension $N$ (levels $0,\ldots,N-1$), final time $t_{\max}$, and $N_t$ sampled output times. |

The simulation calculates atomic populations $P_e(t)$ and $P_g(t)$, inversion $W(t)=P_e(t)-P_g(t)$, mean photon number $\langle n(t)\rangle$, atomic reduced-state von Neumann entropy in bits, and total-excitation mean $\langle M(t)\rangle=\langle n(t)\rangle+P_e(t)$. It also shows the photon distribution $P_n(t)$, Mandel parameter $Q(t)=[\langle n(n-1)\rangle-\langle n\rangle^2]/\langle n\rangle$, and equal-time $g^{(2)}(0,t)=\langle n(n-1)\rangle/\langle n\rangle^2$. Values of $Q$ and $g^{(2)}$ are undefined when their denominators vanish. The statistics describe the single-mode field at each time, not two-time correlations.

Conservation of $\langle M\rangle$ provides an internal check. The app warns if the highest retained Fock level becomes appreciably populated. Increase $N$ and compare results to check truncation convergence, particularly for squeezed and highly populated states. Increasing $N_t$ adds sampled times; it does **not** by itself control integrator error or demonstrate convergence.

## Use

1. Set the physical and numerical parameters, then select the initial states.
2. Read the RWA scale estimates and any truncation warning. Click **Run simulation**.
3. Inspect the curves and photon-distribution map. Drag **Inspect time** to compare atom populations with the instantaneous photon-number distribution.
4. Click **Generate MP4 video** for an animation of those populations and the developing inversion curve. The embedded video has playback controls and a download button. It uses at most 90 evenly spaced frames from the existing solution; it does not solve the equation again.
5. Download the six-panel figure as PDF or PNG, the photon map as PNG, or the sampled data as CSV. Exports correspond to the most recent run. The PDF is a figure, not a parameter report.

A validation case is $\omega_c=\omega_{eg}=1$, $g=0.05$, initially excited atom, and vacuum field. Analytically, $P_e(t)=\cos^2(gt)$ and $\langle M\rangle=1$. Use a sufficiently long interval to see an oscillation. For a coherent field, repeat with larger $N$ before interpreting collapse-and-revival behavior.

## Local installation

Use Python 3.11 or newer and install the pinned dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
streamlit run app.py
```

On Windows Command Prompt, use `.venv\Scripts\activate.bat` instead of `source .venv/bin/activate`. In PowerShell, use `.\.venv\Scripts\Activate.ps1`.

The video renderer uses `imageio-ffmpeg`, which bundles an executable; a separate system FFmpeg installation is not required. If video generation fails on a particular host, the static figures, data export, and time slider remain available.

All application logic is in **one `app.py`**. `requirements.txt` supplies dependencies and `.gitignore` excludes temporary files. This repository currently has **no license**: public visibility alone does not grant permission to redistribute or adapt its code.
