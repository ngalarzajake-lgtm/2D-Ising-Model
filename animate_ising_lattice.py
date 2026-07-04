from pathlib import Path
import os

os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib-cache").resolve()))

import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter
from matplotlib.colors import ListedColormap
import numpy as np


FIGURE_DIR = Path("figures")
OUTPUT_FILE = FIGURE_DIR / "ising_lattice_evolution.gif"

LATTICE_SIZE = 72
FRAMES = 72
SWEEPS_PER_FRAME = 3
RANDOM_SEED = 50013

TEMPERATURES = [
    {"alpha": 1.600, "label": "cold: alpha = 1.6"},
    {"alpha": 2.269, "label": "critical-ish: alpha = 2.269"},
    {"alpha": 3.400, "label": "hot: alpha = 3.4"},
]


def checkerboard_masks(size):
    rows, cols = np.indices((size, size))
    even = (rows + cols) % 2 == 0
    return even, ~even


def metropolis_half_step(spins, alpha, mask, rng):
    neighbors = (
        np.roll(spins, 1, axis=0)
        + np.roll(spins, -1, axis=0)
        + np.roll(spins, 1, axis=1)
        + np.roll(spins, -1, axis=1)
    )
    delta_energy = 2 * spins * neighbors
    accept = (delta_energy <= 0) | (
        rng.random(spins.shape) < np.exp(-delta_energy / alpha)
    )
    spins[mask & accept] *= -1


def metropolis_sweep(spins, alpha, masks, rng):
    for mask in masks:
        metropolis_half_step(spins, alpha, mask, rng)


def simulate_histories():
    rng = np.random.default_rng(RANDOM_SEED)
    masks = checkerboard_masks(LATTICE_SIZE)
    histories = []

    for setup in TEMPERATURES:
        spins = rng.choice([-1, 1], size=(LATTICE_SIZE, LATTICE_SIZE)).astype(np.int8)
        frames = []

        for _ in range(FRAMES):
            for _ in range(SWEEPS_PER_FRAME):
                metropolis_sweep(spins, setup["alpha"], masks, rng)
            frames.append(spins.copy())

        histories.append(np.array(frames, dtype=np.int8))

    return histories


def save_animation(histories):
    FIGURE_DIR.mkdir(exist_ok=True)

    cmap = ListedColormap(["#1d4ed8", "#f8fafc"])
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.8), dpi=120)
    images = []

    for ax, setup, frames in zip(axes, TEMPERATURES, histories):
        image = ax.imshow(frames[0], cmap=cmap, vmin=-1, vmax=1, interpolation="nearest")
        ax.set_title(setup["label"], fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        images.append(image)

    sweep_label = fig.text(0.5, 0.04, "", ha="center", fontsize=11)
    fig.suptitle("2D Ising model lattice evolution", y=0.98, fontsize=14)
    fig.tight_layout(rect=(0, 0.07, 1, 0.93))

    writer = PillowWriter(fps=12)
    with writer.saving(fig, OUTPUT_FILE, dpi=120):
        for frame_index in range(FRAMES):
            sweep = (frame_index + 1) * SWEEPS_PER_FRAME
            sweep_label.set_text(f"{sweep} Monte Carlo sweeps")
            for image, frames in zip(images, histories):
                image.set_data(frames[frame_index])
            writer.grab_frame()

    plt.close(fig)


def main():
    histories = simulate_histories()
    save_animation(histories)
    print(f"Saved lattice animation to {OUTPUT_FILE.resolve()}")


if __name__ == "__main__":
    main()
