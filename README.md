# 2D Ising Model

Tiny spins, big phase transition.

This was my **CH.50013 Statistical Thermodynamics** project: a Monte Carlo simulation of the 2D Ising model. The fun part is watching a random grid of up/down spins suddenly start behaving like a material with order, disorder, and a critical point.

## Lattice In Motion

Here is the lattice evolving over time at three temperatures:

![2D Ising lattice animation](figures/ising_lattice_evolution.gif)

What is happening:

- Cold lattice, `alpha = 1.6`: spins quickly form big ordered regions.
- Near critical, `alpha = 2.269`: the lattice gets dramatic, with large clusters constantly forming and breaking.
- Hot lattice, `alpha = 3.4`: thermal noise wins and the grid stays messy.

The critical temperature for the infinite 2D square-lattice Ising model is about `alpha_c = 2.269`, so the middle panel is the interesting one.

## The Bigger Run

The main C++ simulation runs the Metropolis algorithm for two lattice sizes:

- `N = 20`
- `N = 40`

It scans temperatures from `alpha = 1.00` to `alpha = 4.00` and records energy, magnetization, specific heat, susceptibility, and acceptance ratio.

![Summary of Ising observables](figures/ising_observables_summary.png)

The cool science bit:

- Magnetization drops as the system goes from ordered to disordered.
- Specific heat and susceptibility spike near the critical region.
- The `N = 40` lattice has sharper peaks than `N = 20`, which is what you expect as the system gets closer to the thermodynamic limit.

## Files

| File | What it is |
| --- | --- |
| `Ising_Model.cpp` | C++ Monte Carlo simulation for the full data run |
| `animate_ising_lattice.py` | Fast visual simulation that creates the lattice GIF |
| `plot_ising_results.py` | Makes the summary plots from the CSV files |
| `ising_data_N20.csv` | Data from the `N = 20` run |
| `ising_data_N40.csv` | Data from the `N = 40` run |
| `figures/` | Plots and the lattice animation |
| `Ising_Model_Documentation.pdf` | Full write-up/report |

## How To Make The Visuals Again

Regenerate the lattice animation:

```bash
python animate_ising_lattice.py
```

Regenerate the summary plots from the included CSV files:

```bash
python plot_ising_results.py
```

The full C++ data run is intentionally not the quick demo. It uses long equilibration and sampling runs, so the included CSV files are the saved results used for the plots.

## Built For

**Jacob Emmanuel Sadorra**  
CH.50013 Statistical Thermodynamics  
June 2026
