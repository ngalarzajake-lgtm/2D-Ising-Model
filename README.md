# 2D Ising Model Monte Carlo Simulation

Course project for **CH.50013 Statistical Thermodynamics**.

This repository contains a C++ implementation of the two-dimensional ferromagnetic Ising model using the Metropolis Monte Carlo algorithm, along with generated simulation data, plots, and the final report.

## Project Overview

The model simulates an \(N \times N\) square lattice of spins \(s_i = \pm 1\) with periodic boundary conditions and zero external field. The coupling constant is set to \(J = 1\), and temperature is represented by the dimensionless parameter

\[
\alpha = \frac{k_B T}{J}.
\]

The simulation compares finite lattices with \(N=20\) and \(N=40\), then evaluates how the observables behave near the exact infinite-lattice critical point

\[
\alpha_c = \frac{2}{\ln(1+\sqrt{2})} \approx 2.269.
\]

## What Is Included

| Path | Description |
| --- | --- |
| `Ising_Model.cpp` | Main C++ simulation code |
| `plot_ising_results.py` | Python plotting script for the generated CSV files |
| `ising_data_N20.csv` | Monte Carlo output for a 20 x 20 lattice |
| `ising_data_N40.csv` | Monte Carlo output for a 40 x 40 lattice |
| `figures/` | Generated PNG and PDF plots |
| `Ising_Model_Documentation.tex` | LaTeX source for the final report |
| `Ising_Model_Documentation.pdf` | Compiled final report |

## Methods

For each lattice size and temperature point:

- Spins are initialized randomly as \(+1\) or \(-1\).
- A trial spin flip is selected at a random lattice site.
- The energy change is computed from the four nearest neighbors.
- The move is accepted if it lowers the energy, or accepted probabilistically using \(e^{-\Delta E / \alpha}\).
- Measurements are collected after equilibration.

Simulation settings used for the included data:

| Setting | Value |
| --- | --- |
| Lattice sizes | \(N=20\), \(N=40\) |
| Temperature range | \(\alpha = 1.00\) to \(4.00\) |
| Temperature step | \(0.01\) |
| Equilibration | 100,000 sweeps per temperature |
| Sampling | 300,000 sweeps per temperature |
| One sweep | \(N^2\) attempted spin flips |

## Results Summary

The finite-size simulations reproduce the expected qualitative transition from an ordered low-temperature state to a disordered high-temperature state.

Key observations:

- The magnetization decreases sharply near the critical region.
- The specific heat peak occurs near \(\alpha \approx 2.31\) for \(N=20\) and \(\alpha \approx 2.29\) for \(N=40\).
- The susceptibility peak occurs near \(\alpha \approx 2.35\) for \(N=20\) and \(\alpha \approx 2.31\) for \(N=40\).
- The larger \(N=40\) lattice shows sharper finite-size behavior closer to the thermodynamic-limit critical point.

![Summary of Ising observables](figures/ising_observables_summary.png)

## Reproducing the Results

Compile the simulation:

```bash
g++ -std=c++17 -O2 Ising_Model.cpp -o ising_model
```

Run the simulation:

```bash
./ising_model
```

On Windows PowerShell, the executable may be run as:

```powershell
.\ising_model.exe
```

The full simulation is computationally expensive because it performs long equilibration and sampling runs for both lattice sizes across 301 temperature values. The included CSV files are the outputs used to generate the submitted figures and report.

Install the plotting dependency:

```bash
python -m pip install -r requirements.txt
```

Regenerate the figures from the included CSV files:

```bash
python plot_ising_results.py
```

## Report

The full write-up is available in [`Ising_Model_Documentation.pdf`](Ising_Model_Documentation.pdf). It includes the model background, algorithm description, finite-size comparison, plots, and additional discussion of Monte Carlo sampling in chemistry.

## Course Context

This project was completed by **Jacob Emmanuel Sadorra** for **CH.50013 Statistical Thermodynamics** in June 2026.
