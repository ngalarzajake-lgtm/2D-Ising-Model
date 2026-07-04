import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt


FIGURE_DIR = Path("figures")
CRITICAL_ALPHA = 2.269


DATA_SETS = [
    {"N": 20, "file": "ising_data_N20.csv", "color": "#1f77b4"},
    {"N": 40, "file": "ising_data_N40.csv", "color": "#d62728"},
]


PLOTS = [
    {
        "column": "E_avg",
        "ylabel": r"$\langle E \rangle$",
        "title": "Average Energy",
        "name": "energy",
    },
    {
        "column": "M_avg_abs",
        "individual_column": "M_per_spin",
        "ylabel": r"$\langle M \rangle$",
        "individual_ylabel": r"$\langle |M| \rangle / N^2$",
        "title": "Magnetization",
        "name": "magnetization",
    },
    {
        "column": "Cv_per_spin",
        "ylabel": r"$C_p$",
        "title": "Specific Heat",
        "name": "specific_heat",
    },
    {
        "column": "chi_per_spin",
        "ylabel": r"$\chi$",
        "title": "Susceptibility",
        "name": "susceptibility",
    },
    {
        "column": "acceptance_ratio",
        "ylabel": "Acceptance ratio",
        "title": "Metropolis acceptance ratio",
        "name": "acceptance_ratio",
    },
]


def exact_magnetization(alpha):
    """Infinite-lattice spontaneous magnetization for the 2D Ising model."""
    if alpha >= CRITICAL_ALPHA:
        return 0.0

    value = 1.0 - math.sinh(2.0 / alpha) ** -4
    return value ** 0.125


def read_data(filename):
    """Read one Ising CSV file and convert every value to a float."""
    with open(filename, newline="") as file:
        rows = [
            {column: float(value) for column, value in row.items()}
            for row in csv.DictReader(file)
        ]

    return sorted(rows, key=lambda row: row["alpha"])


def load_all_data():
    all_data = []

    for data_set in DATA_SETS:
        rows = read_data(data_set["file"])
        all_data.append({**data_set, "rows": rows})

    return all_data


def format_axis(ax, ylabel):
    ax.set_xlabel("Temperature")
    ax.set_ylabel(ylabel)

    ax.axvline(
        CRITICAL_ALPHA,
        color="0.25",
        linestyle=":",
        linewidth=1.4,
        label=r"$\alpha_c \approx 2.269$",
    )

    ax.legend(frameon=False)


def plot_quantity(ax, all_data, plot_info, *, individual=False):
    column = plot_info["column"]
    ylabel = plot_info["ylabel"]

    if individual:
        column = plot_info.get("individual_column", column)
        ylabel = plot_info.get("individual_ylabel", ylabel)

    for data_set in all_data:
        rows = data_set["rows"]
        ax.plot(
            [row["alpha"] for row in rows],
            [row[column] for row in rows],
            color=data_set["color"],
            linewidth=1.8,
            label=f"N = {data_set['N']}",
        )

    if individual and plot_info["name"] == "magnetization":
        alpha_values = [row["alpha"] for row in all_data[0]["rows"]]
        ax.plot(
            alpha_values,
            [exact_magnetization(alpha) for alpha in alpha_values],
            color="black",
            linestyle="--",
            linewidth=1.5,
            label="Onsager-Yang exact",
        )

    ax.set_title(plot_info["title"])
    format_axis(ax, ylabel)


def save_summary_figure(all_data):
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.2), dpi=180)

    for ax, plot_info in zip(axes.ravel(), PLOTS[:4]):
        plot_quantity(ax, all_data, plot_info)

    fig.tight_layout()

    fig.savefig(FIGURE_DIR / "ising_observables_summary.png", bbox_inches="tight")
    fig.savefig(FIGURE_DIR / "ising_observables_summary.pdf", bbox_inches="tight")

    plt.close(fig)


def save_individual_figures(all_data):
    for plot_info in PLOTS:
        fig, ax = plt.subplots(figsize=(7, 4.8), dpi=180)

        plot_quantity(ax, all_data, plot_info, individual=True)
        fig.tight_layout()

        fig.savefig(FIGURE_DIR / f"ising_{plot_info['name']}.png", bbox_inches="tight")
        fig.savefig(FIGURE_DIR / f"ising_{plot_info['name']}.pdf", bbox_inches="tight")

        plt.close(fig)


def print_peak_locations(all_data):
    print("Peak estimates from the Monte Carlo data:")

    for data_set in all_data:
        rows = data_set["rows"]
        cv_peak = max(rows, key=lambda row: row["Cv_per_spin"])
        chi_peak = max(rows, key=lambda row: row["chi_per_spin"])
        print(
            f"N={data_set['N']}: "
            f"Cv peak at alpha={cv_peak['alpha']:.3f}, "
            f"chi peak at alpha={chi_peak['alpha']:.3f}"
        )


def main():
    FIGURE_DIR.mkdir(exist_ok=True)
    all_data = load_all_data()
    save_summary_figure(all_data)
    save_individual_figures(all_data)
    print_peak_locations(all_data)
    print(f"Saved graphs in {FIGURE_DIR.resolve()}")


if __name__ == "__main__":
    main()
