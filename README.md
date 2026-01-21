# Water-Melt-PT.jl

[![Julia](https://img.shields.io/badge/Julia-1.9%2B-blue.svg)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Path-dependent multi-system phase diagram modeling to calculate water-present partial melting reactions in a heterogeneous crust.**

## Overview

Water-Melt-PT.jl is a Julia package for modeling water release and melt production during crustal metamorphism. It uses the [MAGEMin](https://github.com/ComputationalThermodynamics/MAGEMin) thermodynamic engine to simulate:

- **Prograde metamorphism** of metasedimentary rocks (metapelites, metagraywackes)
- **Water release** from dehydration reactions along a P-T path
- **Garnet fractionation** during metamorphic crystallization
- **Water-fluxed melting** of orthogneiss from externally-derived water

The package is designed to investigate how water released from metasedimentary rocks can drive partial melting in adjacent orthogneiss, a key process in crustal anatexis.

## Installation

### Prerequisites

- Julia 1.9 or higher
- MAGEMin_C package

### Install from GitHub

```julia
using Pkg
Pkg.add(url="https://github.com/jvanardo/Water-Melt-PT.jl")
```

### Install dependencies manually

```julia
using Pkg
Pkg.add(["MAGEMin_C", "CairoMakie", "CSV", "DataFrames", "ProgressMeter", "Reexport"])
```

### Development installation

```bash
git clone https://github.com/jvanardo/Water-Melt-PT.jl.git
cd Water-Melt-PT.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Quick Start

```julia
using WaterMeltPT

# Run simulation with default parameters
results = run_water_melt_simulation()

# Generate plots
plot_results(results)

# Save results to CSV
save_results(results, "output/")
```

## Usage

### Basic Configuration

```julia
using WaterMeltPT

# Define P-T path configuration
config = PTPathConfig(
    Tmin = 600.0,      # Starting temperature (°C)
    Tmax = 750.0,      # Ending temperature (°C)
    T_step = 1.0,      # Temperature increment (°C)
    Pmin = 14.0,       # Starting pressure (kbar)
    Pmax = 9.5,        # Ending pressure (kbar)
    g_factor = 1.0     # Garnet fractionation factor (0-1)
)

# Run simulation
results = run_water_melt_simulation(config)
```

### Custom Bulk Compositions

```julia
# Define custom bulk compositions (wt%)
my_compositions = Dict(
    "metapelite" => BulkComposition(
        "my_metapelite",
        # SiO2, Al2O3, CaO, MgO, FeO, Fe2O3, K2O, Na2O, TiO2, MnO, H2O
        [59.0, 18.0, 0.6, 2.2, 6.3, 0.0, 3.7, 1.3, 0.8, 0.0, 7.5],
        "My reference"
    ),
    "metagraywacke" => BulkComposition(
        "my_metagraywacke",
        [70.5, 13.0, 1.7, 2.4, 4.9, 0.0, 2.4, 3.0, 0.7, 0.0, 7.5],
        "My reference"
    ),
    "orthogneiss" => BulkComposition(
        "my_orthogneiss", 
        [72.0, 15.0, 1.3, 0.8, 2.4, 0.0, 4.7, 3.0, 0.3, 0.0, 7.5],
        "My reference"
    )
)

results = run_water_melt_simulation(config; compositions=my_compositions)
```

### Accessing Results

```julia
# Access lithology-specific results
mp = results.lithologies["metapelite"]

# Solidus temperature
println("Metapelite solidus: $(mp.T_solidus)°C")

# Cumulative water release
println("Water released: $(mp.H2O_frac_total_cum[end] * 100)%")

# Phase assemblage at a specific step
step = 50
println("Phases at T=$(results.T_array[step])°C: $(mp.out[step].ph)")
```

## Default Bulk Compositions

| Lithology | Reference |
|-----------|-----------|
| Metapelite | Forshaw & Pattison (2024) |
| Metagraywacke | Villaros et al. (2018) |
| Orthogneiss | Weisbrod (1970) |

## Output

The package generates:

1. **H2O_frac.svg** - Cumulative water release plot
2. **Liq_frac.svg** - Melt fraction comparison plot
3. **Mode_Boxes.svg** - Phase mode diagrams for all lithologies
4. **simulation_results.csv** - Tabulated results

## Scientific Background

During prograde metamorphism, hydrous minerals in metasedimentary rocks break down through a series of dehydration reactions. The water released migrates upward and can flux adjacent rocks, lowering their solidus and enabling partial melting at lower temperatures than would occur under water-absent conditions.

This package models this process by:

1. Computing equilibrium phase assemblages along a P-T path
2. Extracting excess water at each step
3. Fractionating garnet to simulate its removal from the reactive bulk
4. Transferring released water to orthogneiss to simulate water-fluxed melting

## Citation

If you use this package in your research, please cite:

```bibtex
@software{WaterMeltPT,
  author = {Vanardois, J.},
  title = {Water-Melt-PT.jl: Path-dependent water-melt modeling in heterogeneous crust},
  year = {2026},
  url = {https://github.com/jvanardo/Water-Melt-PT.jl}
}
```

## References

- Forshaw, J.B. & Pattison, D.R.M. (2024) - Metapelite bulk compositions
- Villaros, A. et al. (2018) - Metagraywacke bulk compositions  
- Weisbrod, A. (1970) - Orthogneiss bulk compositions
- Lanari, P. & Tedeschi, M. (2024) - Phase color conventions
- Riel, N. et al. (2022) - MAGEMin thermodynamic engine

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
