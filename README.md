# 2D Ising Model

Tiny spins, big phase transition.

I made this from a graduate-level statistical thermodynamics course: **CH50013_2026_1**. It is a Monte Carlo simulation of the 2D Ising model, which is one of those physics models that looks almost too simple at first: a grid of spins that can only point up or down. Then, as temperature changes, the grid starts acting like a real material with order, disorder, and a phase transition.

## Lattice In Motion

Here is the lattice evolving over time at three temperatures:

![2D Ising lattice animation](figures/ising_lattice_evolution.gif)

At low temperature, neighboring spins prefer to line up, so large ordered regions appear. At high temperature, thermal motion keeps breaking the order. Near the critical temperature, the lattice becomes the most interesting: clusters form at many length scales, and the system sits right between ordered and disordered behavior.

## Ising Model Background

The 2D Ising model places a spin on every site of a square lattice:

```text
s_i = +1 or -1
```

For the ferromagnetic model, neighboring spins want to point in the same direction. The energy is usually written as:

```text
E = -J * sum_<i,j>(s_i s_j)
```

`J` is the interaction strength, and the sum runs over nearest-neighbor pairs. If two neighboring spins are aligned, `s_i s_j = +1`, which lowers the energy. If they are opposite, `s_i s_j = -1`, which raises it.

The probability of seeing a spin configuration comes from the Boltzmann weight:

```text
P(configuration) = exp(-E / (k_B T)) / Z
```

where `Z` is the partition function: the sum of Boltzmann weights over every possible spin configuration. That is the beautiful but painful part: for an `N x N` lattice, there are `2^(N^2)` possible configurations. Monte Carlo lets us sample the important configurations instead of trying to enumerate all of them.

In this project I used the dimensionless temperature:

```text
alpha = k_B T / J
K = J / (k_B T) = 1 / alpha
```

The code uses `J = 1`, so `alpha` is the temperature scale.

## Metropolis Update

The simulation flips one randomly chosen spin at a time. Only the four nearest neighbors matter for the energy change:

```text
Delta E = 2 * J * s_ij * (s_up + s_down + s_left + s_right)
```

The Metropolis rule is:

```text
accept the flip if Delta E <= 0
otherwise accept with probability exp(-Delta E / (k_B T))
```

With `J = 1`, that becomes:

```text
accept with probability exp(-Delta E / alpha)
```

This is what lets the lattice sometimes make an energetically bad move. That randomness is important because real thermal systems fluctuate, especially near the critical point.

## Onsager's Critical Temperature

The 2D Ising model is famous because Lars Onsager solved the square-lattice model exactly in 1944. One of the biggest results is the exact critical temperature for the infinite lattice.

Using:

```text
K = J / (k_B T)
```

the critical point satisfies:

```text
sinh(2K_c) = 1
```

Solving that gives:

```text
K_c = 0.5 * ln(1 + sqrt(2))
```

Since `alpha = 1 / K`, the critical dimensionless temperature is:

```text
alpha_c = 2 / ln(1 + sqrt(2)) ~= 2.269
```

That number is the vertical reference line in the plots. My finite lattices do not produce a perfectly sharp singularity because `N = 20` and `N = 40` are still small compared with the infinite model. Instead, the transition shows up as rounded peaks in specific heat and susceptibility near Onsager's value.

The spontaneous magnetization also has a famous exact form, later derived by Yang:

```text
M = [1 - sinh(2K)^(-4)]^(1/8), for alpha < alpha_c
M = 0, for alpha >= alpha_c
```

That is why the magnetization curve drops so sharply near the critical region.

## The Bigger Run

The main C++ simulation runs the Metropolis algorithm for two lattice sizes:

- `N = 20`
- `N = 40`

It scans temperatures from `alpha = 1.00` to `alpha = 4.00` and records energy, magnetization, specific heat, susceptibility, and acceptance ratio.

![Summary of Ising observables](figures/ising_observables_summary.png)

## Built For

**Jacob Emmanuel Sadorra**  
Graduate-level course: **CH50013_2026_1**  
June 2026
