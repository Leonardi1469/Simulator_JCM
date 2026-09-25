"""Interactive, closed Jaynes–Cummings siulator (single-file Streamlit app)."""

from __future__ import annotations

from io import BytesIO, StringIO
from pathlib import Path
from tempfile import TemporaryDirectory

import imageio_ffmpeg
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
from matplotlib.animation import FFMpegWriter, FuncAnimation
from matplotlib.patches import Ellipse
from qutip import (
    basis,
    coherent,
    destroy,
    displace,
    fock,
    qeye,
    sesolve,
    sigmam,
    sigmap,
    sigmaz,
    squeeze,
    tensor,
)


def cavity_diagram():
    """A simple, original illustration of one two-level atom and one mode."""
    fig, ax = plt.subplots(figsize=(7.5, 2.6), dpi=130)
    ax.add_patch(Ellipse((0, 0), 6.1, 1.95, facecolor="#eff7ff", edgecolor="#37648c", lw=2))
    x = np.linspace(-2.45, 2.45, 350)
    ax.plot(x, 0.23 * np.sin(8 * x), color="#237dbe", lw=2, alpha=0.85)
    ax.plot([-0.35, 0.42], [0.69, 0.69], color="#b84d53", lw=2)
    ax.plot([-0.35, 0.42], [-0.63, -0.63], color="#b84d53", lw=2)
    ax.annotate("", xy=(0.05, 0.55), xytext=(0.05, -0.49),
                arrowprops={"arrowstyle": "<->", "color": "#b84d53", "lw": 1.6})
    ax.text(0.53, 0.73, r"$|e\rangle$", fontsize=11)
    ax.text(0.53, -0.65, r"$|g\rangle$", fontsize=11)
    ax.text(-2.6, 0.52, r"cavity mode $\omega_c$", fontsize=10)
    ax.text(1.12, -0.55, r"coupling $g$", fontsize=10)
    ax.text(0.58, 0.13, r"$\omega_{eg}$", color="#b84d53", fontsize=10)
    ax.set(xlim=(-3.25, 3.25), ylim=(-1.25, 1.25))
    ax.axis("off")
    fig.tight_layout(pad=0.2)
    return fig


def make_field(config):
    """The order of displacement and squeezing is defined explicitly."""
    n = config["N"]
    kind = config["field"]
    if kind == "Vacuum":
        return fock(n, 0)
    if kind == "Fock":
        return fock(n, config["n0"])
    if kind == "Coherent":
        return coherent(n, config["alpha"])
    squeezed = squeeze(n, config["r"]) * fock(n, 0)
    if kind == "Squeezed vacuum":
        return squeezed
    return displace(n, config["alpha"]) * squeezed


def make_atom(config):
    if config["atom"] == "Excited":
        return basis(2, 0)
    if config["atom"] == "Ground":
        return basis(2, 1)
    return (np.cos(config["theta"] / 2) * basis(2, 0)
            + np.exp(1j * config["phi"]) * np.sin(config["theta"] / 2) * basis(2, 1))


def rwa_guidance(config, mean_initial):
    total_frequency = config["wc"] + config["wa"]
    if total_frequency <= 0:
        return None
    detuning_ratio = abs(config["wa"] - config["wc"]) / total_frequency
    coupling_ratio = config["g"] * np.sqrt(mean_initial + 1) / total_frequency
    return detuning_ratio, coupling_ratio


def run_simulation(config):
    """Solve the Schrödinger equation and evaluate photon and atom statistics."""
    n = config["N"]
    times = np.linspace(0, config["tmax"], config["Nt"])
    field = make_field(config)
    atom = make_atom(config)
    psi0 = tensor(field, atom)
    a = tensor(destroy(n), qeye(2))
    sz = tensor(qeye(n), sigmaz())
    sp = tensor(qeye(n), sigmap())
    sm = tensor(qeye(n), sigmam())
    hamiltonian = (config["wc"] * a.dag() * a
                   + 0.5 * config["wa"] * sz
                   + config["g"] * (a * sp + a.dag() * sm))
    evolution = sesolve(hamiltonian, psi0, times)

    # Tensor ordering is |n> ⊗ |e/g>; the atomic excited state has index 0.
    amplitudes = np.stack([state.full().reshape(n, 2) for state in evolution.states])
    populations = np.sum(np.abs(amplitudes) ** 2, axis=2)
    excited = np.sum(np.abs(amplitudes[:, :, 0]) ** 2, axis=1)
    ground = np.sum(np.abs(amplitudes[:, :, 1]) ** 2, axis=1)
    inversion = excited - ground
    photon_numbers = np.arange(n)
    mean_photons = populations @ photon_numbers
    factorial_moment = populations @ (photon_numbers * (photon_numbers - 1))
    excitations = mean_photons + excited

    # The global state is pure; its reduced atomic density matrix has size 2×2.
    atomic_density = np.einsum("tni,tnj->tij", amplitudes, amplitudes.conj())
    eigenvalues = np.clip(np.linalg.eigvalsh(atomic_density), 0.0, 1.0)
    entropy_terms = np.zeros_like(eigenvalues)
    positive = eigenvalues > 0
    entropy_terms[positive] = -eigenvalues[positive] * np.log2(eigenvalues[positive])
    entropy = np.sum(entropy_terms, axis=1)

    mandel = np.full_like(mean_photons, np.nan)
    g2 = np.full_like(mean_photons, np.nan)
    nonzero = mean_photons > 1e-9
    mandel[nonzero] = (factorial_moment[nonzero]
                       - mean_photons[nonzero] ** 2) / mean_photons[nonzero]
    g2[nonzero] = factorial_moment[nonzero] / mean_photons[nonzero] ** 2

    return {
        "config": config.copy(), "times": times, "populations": populations,
        "excited": excited, "ground": ground, "inversion": inversion,
        "mean_photons": mean_photons, "excitations": excitations,
        "entropy": entropy, "mandel": mandel, "g2": g2,
        "mean_initial": float(populations[0] @ photon_numbers),
        "edge_max": float(np.max(populations[:, -1])),
    }


def style_axis(ax):
    ax.grid(alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=9)


def main_figure(data):
    t = data["times"]
    fig, axes = plt.subplots(3, 2, figsize=(8.1, 9.5), dpi=130)
    specifications = [
        (data["excited"], r"$P_e(t)$", "Excited-state probability", (0, 1.02)),
        (data["ground"], r"$P_g(t)$", "Ground-state probability", (0, 1.02)),
        (data["inversion"], r"$W(t)$", "Atomic inversion", (-1.05, 1.05)),
        (data["mean_photons"], r"$\langle n(t)\rangle$", "Mean photon number", None),
        (data["entropy"], r"$S_A(t)$ (bits)", "Atomic von Neumann entropy", (0, 1.02)),
        (data["excitations"], r"$\langle M(t)\rangle$", "Total excitations", None),
    ]
    for ax, (values, ylabel, title, ylim) in zip(axes.flat, specifications):
        ax.plot(t, values, lw=1.55, color="#176ba0")
        ax.set(xlabel=r"Time $t$", ylabel=ylabel, title=title, xlim=(t[0], t[-1]))
        if ylim is not None:
            ax.set_ylim(*ylim)
        elif title == "Total excitations":
            center = float(np.mean(values))
            ax.set_ylim(center - 0.35, center + 0.35)
        style_axis(ax)
    fig.suptitle("Jaynes–Cummings dynamics", fontsize=15)
    fig.tight_layout()
    return fig


def photon_map(data):
    populations = data["populations"]
    times = data["times"]
    n = populations.shape[1]
    occupied = np.flatnonzero(np.max(populations, axis=0) > 1e-3)
    upper = min(n - 1, max(4, int(occupied[-1]) + 2)) if occupied.size else min(n - 1, 4)
    fig, ax = plt.subplots(figsize=(8, 3.7), dpi=130)
    image = ax.imshow(populations.T, origin="lower", interpolation="nearest",
                      aspect="auto", extent=[times[0], times[-1], -0.5, n - 0.5],
                      vmin=0, cmap="viridis")
    ax.set(xlabel=r"Time $t$", ylabel="Photon number $n$",
           title=r"Photon-number probabilities $P_n(t)$", ylim=(-0.5, upper + 0.5))
    if upper == n - 1:
        ax.axhline(n - 1, color="white", ls="--", lw=1)
    fig.colorbar(image, ax=ax, label=r"$P_n(t)$")
    fig.tight_layout()
    return fig


def figure_bytes(fig, file_format, dpi=170):
    output = BytesIO()
    fig.savefig(output, format=file_format, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output.getvalue()


def csv_bytes(data):
    columns = ["time", "P_excited", "P_ground", "inversion", "mean_photons",
               "atom_entropy_bits", "mean_total_excitations", "Mandel_Q", "g2_equal_time"]
    values = np.column_stack([data[key] for key in [
        "times", "excited", "ground", "inversion", "mean_photons",
        "entropy", "excitations", "mandel", "g2",
    ]])
    output = StringIO()
    output.write("# Simulator JCM; model units: hbar=1; undefined statistics: nan\n")
    output.write("# parameters: " + "; ".join(f"{k}={v}" for k, v in data["config"].items()) + "\n")
    np.savetxt(output, values, delimiter=",", header=",".join(columns), comments="", fmt="%.10g")
    return output.getvalue().encode("utf-8")


def snapshot_figure(data, index):
    """A user-controlled snapshot of the already calculated evolution."""
    t = data["times"]
    probabilities = data["populations"][index]
    occupied = np.flatnonzero(probabilities > 0.001)
    max_n = min(len(probabilities) - 1, max(4, int(occupied[-1]) + 2)) if occupied.size else 4
    fig, (ax_atom, ax_field) = plt.subplots(1, 2, figsize=(8, 2.9), dpi=110)
    ax_atom.bar(["Excited", "Ground"], [data["excited"][index], data["ground"][index]],
                color=["#d56851", "#237dbe"])
    ax_atom.set(ylabel="Probability", ylim=(0, 1.05), title="Atom")
    ax_field.bar(np.arange(max_n + 1), probabilities[:max_n + 1], color="#237dbe")
    ax_field.set(xlabel="Photon number", ylabel=r"$P_n$", title="Cavity field",
                 xlim=(-0.6, max_n + 0.6), ylim=(0, max(0.2, probabilities.max() * 1.12)))
    ax_field.set_xticks(np.arange(max_n + 1))
    fig.suptitle(f"Time t = {t[index]:.3f}")
    fig.tight_layout()
    return fig


def make_video(data, max_frames=90, fps=12):
    """Build an MP4 for in-app playback from stored arrays; no new solve."""
    frame_indices = np.unique(np.linspace(0, len(data["times"]) - 1,
                                          min(max_frames, len(data["times"])), dtype=int))
    t = data["times"]
    probs = data["populations"]
    occupied = np.flatnonzero(np.max(probs, axis=0) > 0.001)
    max_n = min(probs.shape[1] - 1, max(4, int(occupied[-1]) + 2)) if occupied.size else 4
    fig, (ax_curve, ax_atom, ax_field) = plt.subplots(
        1, 3, figsize=(10, 3.2), dpi=90, gridspec_kw={"width_ratios": [1.6, 0.8, 1.4]}
    )
    ax_curve.set(xlim=(t[0], t[-1]), ylim=(-1.05, 1.05),
                 xlabel="Time", ylabel=r"$W(t)$", title="Atomic inversion")
    ax_curve.grid(alpha=0.25)
    curve, = ax_curve.plot([], [], lw=2, color="#176ba0")
    cursor = ax_curve.axvline(t[0], ls="--", color="#d56851")
    atom_bars = ax_atom.bar(["Excited", "Ground"], [0, 0], color=["#d56851", "#237dbe"])
    ax_atom.set(title="Atom", ylim=(0, 1.05), ylabel="Probability")
    photon_bars = ax_field.bar(np.arange(max_n + 1), np.zeros(max_n + 1), color="#237dbe")
    ax_field.set(title="Cavity field", xlabel="Photon number", ylabel=r"$P_n$",
                 xlim=(-0.6, max_n + 0.6), ylim=(0, 1.05))
    ax_field.set_xticks(np.arange(max_n + 1, max(1, (max_n + 1) // 8)))
    heading = fig.suptitle("")
    fig.tight_layout()

    def update(index):
        curve.set_data(t[:index + 1], data["inversion"][:index + 1])
        cursor.set_xdata([t[index], t[index]])
        for bar, value in zip(atom_bars, [data["excited"][index], data["ground"][index]]):
            bar.set_height(value)
        for bar, value in zip(photon_bars, probs[index, :max_n + 1]):
            bar.set_height(value)
        heading.set_text(f"Jaynes–Cummings evolution · t = {t[index]:.3f}")

    animation = FuncAnimation(fig, update, frames=frame_indices, blit=False)
    try:
        mpl.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
        with TemporaryDirectory() as directory:
            path = Path(directory) / "simulator_jcm.mp4"
            writer = FFMpegWriter(fps=fps, codec="libx264", bitrate=1200,
                                  extra_args=["-pix_fmt", "yuv420p", "-movflags", "+faststart"])
            animation.save(str(path), writer=writer, dpi=90)
            return path.read_bytes()
    finally:
        plt.close(fig)


st.set_page_config(page_title="Simulator JCM", layout="wide")
st.title("Jaynes–Cummings model simulator")
st.caption("Simulator JCM · Interactive teaching and research application · 2026")

with st.expander("Physical model and approximation", expanded=False):
    diagram = cavity_diagram()
    st.pyplot(diagram)
    plt.close(diagram)
    st.markdown(r"""
One **two-level atom** interacts with **one quantized cavity mode**. In units where $\hbar=1$,

$$H=\omega_c a^\dagger a+\frac{\omega_{eg}}{2}\sigma_z
       +g(a\sigma_+ + a^\dagger\sigma_-).$$

- $\omega_c a^\dagger a$ is the energy of the cavity field; $a$ and $a^\dagger$ annihilate and create a photon.
- $\omega_{eg}\sigma_z/2$ is the energy of the atom; $\omega_{eg}$ is its transition frequency.
- $g(a\sigma_+ + a^\dagger\sigma_-)$ exchanges one excitation between atom and field.

The **rotating-wave approximation** discards the counter-rotating terms $a\sigma_-$ and
$a^\dagger\sigma_+$. Check that $g\sqrt{n+1}$ is small relative to
$\omega_c+\omega_{eg}$ for photon numbers populated during evolution; the near-resonant
regime also has $|\omega_{eg}-\omega_c|\ll\omega_c+\omega_{eg}$. The displayed 0.1
warning threshold is illustrative, not a universal cutoff. The system is closed:
there are no cavity losses or spontaneous atomic emission.
""")

controls, output = st.columns([1, 2.1], gap="large")
with controls:
    with st.form("jcm_controls"):
        st.header("Physical parameters")
        st.caption("Angular frequencies and coupling share one unit; time uses its inverse.")
        wc = st.number_input("Cavity frequency ωc", min_value=0.0, value=1.0,
                             step=0.01, format="%.3f")
        wa = st.number_input("Atomic frequency ωeg", min_value=0.0, value=1.0,
                             step=0.01, format="%.3f")
        g = st.number_input("Coupling g", min_value=0.0, value=0.05,
                            step=0.01, format="%.3f")
        st.header("Numerical parameters")
        n = st.slider("Fock dimension N", min_value=2, max_value=100, value=20)
        tmax = st.number_input("Final time tmax", min_value=0.1, max_value=1000.0,
                               value=100.0, step=1.0)
        nt = st.number_input("Sampled time points Nt", min_value=100, max_value=5000,
                             value=800, step=100)
        st.caption("N retains |0⟩ through |N−1⟩. Nt changes sampling, not solver accuracy.")
        st.header("Initial states")
        atom = st.radio("Atom", ["Excited", "Ground", "Superposition"])
        theta, phi = 0.0, 0.0
        if atom == "Superposition":
            theta = st.slider("Polar angle θ / π", 0.0, 1.0, 0.5, 0.01) * np.pi
            phi = st.slider("Azimuthal angle φ / π", 0.0, 2.0, 0.0, 0.01) * np.pi
        field = st.radio("Cavity field", ["Vacuum", "Fock", "Coherent", "Squeezed vacuum",
                                             "Displaced squeezed vacuum"])
        n0, alpha, r = 0, 0.0, 0.0
        if field == "Fock":
            n0 = st.number_input("Fock photon number n₀", 0, n - 1, 1)
        if field in ("Coherent", "Displaced squeezed vacuum"):
            alpha = st.number_input("Real coherent amplitude α", value=2.0, step=0.1)
        if field in ("Squeezed vacuum", "Displaced squeezed vacuum"):
            r = st.number_input("Real squeezing parameter r", min_value=0.0,
                                max_value=3.0, value=0.5, step=0.1)
        submitted = st.form_submit_button("Run simulation", type="primary", use_container_width=True)

    config = {"wc": float(wc), "wa": float(wa), "g": float(g), "N": int(n),
              "tmax": float(tmax), "Nt": int(nt), "atom": atom, "theta": float(theta),
              "phi": float(phi), "field": field, "n0": int(n0),
              "alpha": float(alpha), "r": float(r)}

    # Exact mean for the finite-dimensional initial state; no approximation for
    # coherent or squeezed preparations is needed to compute the displayed ratio.
    prepared = make_field(config)
    initial_weights = np.abs(prepared.full().ravel()) ** 2
    mean_initial = float(initial_weights @ np.arange(n))
    ratios = rwa_guidance(config, mean_initial)
    if ratios is None:
        st.warning("Both frequencies are zero. The usual RWA scale test is undefined.")
    else:
        st.caption(f"RWA scale estimates: |Δ|/(ωc+ωeg) = {ratios[0]:.3g}; "
                   f"g√(⟨n⟩₀+1)/(ωc+ωeg) = {ratios[1]:.3g}.")
        if max(ratios) >= 0.1:
            st.warning("At least one RWA ratio is ≥ 0.1. Examine the physical "
                       "interpretation of this parameter set; 0.1 is an illustrative threshold.")
    if np.sum(initial_weights[max(0, n - 2):]) > 1e-3:
        st.warning("The initial field occupies the highest Fock levels. Increase N "
                   "and compare simulations to check convergence.")

with output:
    if submitted:
        try:
            with st.spinner("Solving the Schrödinger equation..."):
                result = run_simulation(config)
                figure = main_figure(result)
                result["main_pdf"] = figure_bytes(figure, "pdf")
                figure = main_figure(result)
                result["main_png"] = figure_bytes(figure, "png", dpi=200)
                figure = photon_map(result)
                result["map_png"] = figure_bytes(figure, "png", dpi=180)
                result["csv"] = csv_bytes(result)
            st.session_state["jcm_result"] = result
            st.session_state.pop("jcm_video", None)
        except Exception as exc:
            st.error(f"Simulation failed: {exc}")

    data = st.session_state.get("jcm_result")
    if data is None:
        st.info("Select parameters and initial states, then click **Run simulation**.")
    else:
        if data["config"] != config:
            st.info("The controls have changed. Displayed results still belong to the last run.")
        if data["edge_max"] > 1e-3:
            st.warning(f"The highest retained Fock level reaches {data['edge_max']:.2%} "
                       "probability. Increase N and check convergence.")
        conservation_error = float(np.max(np.abs(data["excitations"] - data["excitations"][0])))
        st.caption(f"Conservation check: maximum change in ⟨M⟩ = {conservation_error:.2e}. "
                   "This does not replace a convergence check in N.")
        st.subheader("Time-dependent observables")
        st.image(data["main_png"], width="stretch")
        st.subheader("Photon distribution")
        st.image(data["map_png"], width="stretch")

        st.subheader("Photon statistics")
        st.caption("Mandel Q and equal-time g² are undefined where the mean photon number is zero.")
        stat_fig, stat_axes = plt.subplots(1, 2, figsize=(8, 2.6), dpi=110)
        for ax, key, title, ylabel in zip(stat_axes, ["mandel", "g2"],
                                          ["Mandel parameter", "Equal-time correlation"],
                                          [r"$Q(t)$", r"$g^{(2)}(0,t)$"]):
            ax.plot(data["times"], data[key], lw=1.5, color="#176ba0")
            ax.set(xlabel="Time", ylabel=ylabel, title=title, xlim=(0, data["times"][-1]))
            style_axis(ax)
        stat_fig.tight_layout()
        st.pyplot(stat_fig)
        plt.close(stat_fig)

        st.subheader("Explore a selected instant")
        selected = st.slider("Inspect time (sample index)", 0, len(data["times"]) - 1,
                             0, key="jcm_time_index")
        snap = snapshot_figure(data, selected)
        st.pyplot(snap)
        plt.close(snap)
        st.caption(f"t = {data['times'][selected]:.3f}; "
                   f"Pe = {data['excited'][selected]:.4f}; "
                   f"⟨n⟩ = {data['mean_photons'][selected]:.4f}; "
                   f"atomic entropy = {data['entropy'][selected]:.4f} bits.")

        st.subheader("Download results")
        c1, c2, c3, c4 = st.columns(4)
        c1.download_button("Main PDF", data["main_pdf"], "jcm_observables.pdf",
                           "application/pdf", use_container_width=True)
        c2.download_button("Main PNG", data["main_png"], "jcm_observables.png",
                           "image/png", use_container_width=True)
        c3.download_button("Photon map", data["map_png"], "jcm_photon_map.png",
                           "image/png", use_container_width=True)
        c4.download_button("Data CSV", data["csv"], "jcm_data.csv",
                           "text/csv", use_container_width=True)

        st.subheader("Playable animation")
        st.caption("Render the evolving inversion, atomic populations and photon distribution "
                   "into a playable MP4. Playback controls allow pause, rewind and scrubbing. "
                   "Video generation uses at most 90 frames and does not repeat the simulation.")
        if st.button("Generate MP4 video", type="secondary"):
            try:
                with st.spinner("Rendering the MP4 video..."):
                    st.session_state["jcm_video"] = make_video(data)
            except Exception as exc:
                st.error(f"MP4 rendering failed: {exc}. Static figures remain available.")
        if "jcm_video" in st.session_state:
            st.video(st.session_state["jcm_video"], format="video/mp4")
            st.download_button("Download MP4", st.session_state["jcm_video"],
                               "jcm_animation.mp4", "video/mp4")
