"""
    WaterMeltPT

A Julia package for modeling water release and melt production during 
metamorphism using MAGEMin thermodynamic calculations.

This package simulates path-dependent multi-system phase diagram modeling 
to calculate water-present partial melting reactions in a heterogeneous crust.

# Main Features
- Calculate water release along P-T paths for multiple lithologies
- Model garnet fractionation during metamorphism  
- Simulate water-fluxed melting in orthogneiss
- Generate mode box diagrams and melt fraction plots

# References
- Forshaw, J.B. and Pattison, D.R.M. (2024) - Metapelite bulk compositions
- Villaros, A. et al. (2018) - Metagraywacke bulk compositions
- Weisbrod, A. (1970) - Orthogneiss bulk compositions
- Lanari, P. and Tedeschi, M. (2024) - Phase color conventions

# Example
```julia
using WaterMeltPT

# Define configuration
config = PTPathConfig(
    Tmin = 600.0, Tmax = 750.0,
    Pmin = 14.0, Pmax = 9.5,
    g_factor = 1.0
)

# Run simulation
results = run_water_melt_simulation(config)
```
"""
module WaterMeltPT

using Reexport: @reexport
@reexport using MAGEMin_C
using CSV
using DataFrames
using ProgressMeter
using CairoMakie

export PTPathConfig, BulkComposition, SimulationResults
export run_water_melt_simulation, plot_results, save_results
export extract_water_excess!, fractionate_garnet!, find_solidus_temperature

# ============================================================================
# Constants
# ============================================================================

"""Default oxide components for MAGEMin metapelite database."""
const XOXIDES = ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "MnO", "H2O"]

"""Input oxide components (with Fe2O3 separated)."""
const X_INIT_OX = ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "Fe2O3", "K2O", "Na2O", "TiO2", "MnO", "H2O"]

"""Pure H2O composition vector for water extraction calculations."""
const H2O_COMP = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]

"""
Phase color map following Lanari and Tedeschi (2024) conventions.
"""
const PHASE_COLOR_MAP = Dict{String, Any}(
    "liq"  => :red,
    "pl"   => (0.925, 0.863, 0.620),
    "q"    => (1.000, 1.000, 1.000),
    "bi"   => (0.463, 0.224, 0.129),
    "mu"   => (0.851, 0.780, 0.796),
    "afs"  => (0.914, 0.722, 0.800),
    "ilm"  => (0.878, 0.365, 0.165),
    "sill" => (0.000, 0.561, 0.737),
    "and"  => (0.498, 0.702, 0.757),
    "ky"   => (0.000, 0.357, 0.576),
    "crd"  => (0.533, 0.408, 0.659),
    "zo"   => (0.733, 0.706, 0.216),
    "H2O"  => :lightblue,
    "g"    => (0.710, 0.173, 0.122),
    "st"   => (0.886, 0.667, 0.000),
    "ttn"  => (0.745, 0.541, 0.263),
    "ru"   => (0.024, 0.259, 0.463),
    "chl"  => (0.506, 0.753, 0.443),
    "cld"  => (0.369, 0.467, 0.380),
    "opx"  => (0.882, 0.486, 0.565),
    "pat"  => (0.855, 0.514, 0.314)
)

# ============================================================================
# Configuration Types
# ============================================================================

"""
    PTPathConfig

Configuration for P-T path simulation.

# Fields
- `Tmin::Float64`: Minimum temperature (°C)
- `Tmax::Float64`: Maximum temperature (°C)  
- `T_step::Float64`: Temperature step size (°C)
- `Pmin::Float64`: Starting pressure (kbar)
- `Pmax::Float64`: Ending pressure (kbar)
- `g_factor::Float64`: Fraction of garnet to fractionate (0.0-1.0)
- `sys_in::String`: Input system type ("mol" or "wt")
"""
Base.@kwdef struct PTPathConfig
    Tmin::Float64 = 600.0
    Tmax::Float64 = 750.0
    T_step::Float64 = 1.0
    Pmin::Float64 = 14.0
    Pmax::Float64 = 9.5
    g_factor::Float64 = 1.0
    sys_in::String = "mol"
end

"""
    BulkComposition

Represents a bulk rock composition with metadata.

# Fields
- `name::String`: Name of the lithology
- `composition_wt::Vector{Float64}`: Composition in weight percent
- `reference::String`: Literature reference for the composition
"""
struct BulkComposition
    name::String
    composition_wt::Vector{Float64}
    reference::String
end

"""
    LithologyResults

Results for a single lithology along the P-T path.

# Fields
- `name::String`: Lithology name
- `out::Vector`: MAGEMin output structures
- `g_frac::Vector{Float64}`: Garnet fractions at each step
- `g_frac_total_cum::Vector{Float64}`: Cumulative garnet fraction
- `H2O_frac::Vector{Float64}`: Water fractions at each step
- `H2O_frac_total_cum::Vector{Float64}`: Cumulative water fraction
- `melt_frac_total_vol_cum::Vector{Float64}`: Cumulative melt volume fraction
- `T_solidus::Float64`: Solidus temperature
"""
mutable struct LithologyResults
    name::String
    out::Vector{MAGEMin_C.gmin_struct{Float64,Int64}}
    X::Vector{Float64}
    g_frac::Vector{Float64}
    g_frac_vol::Vector{Float64}
    g_frac_total::Vector{Float64}
    g_frac_total_vol::Vector{Float64}
    g_frac_total_cum::Vector{Float64}
    g_frac_total_vol_cum::Vector{Float64}
    H2O_frac::Vector{Float64}
    H2O_frac_vol::Vector{Float64}
    H2O_frac_total::Vector{Float64}
    H2O_frac_total_vol::Vector{Float64}
    H2O_frac_total_cum::Vector{Float64}
    H2O_frac_total_vol_cum::Vector{Float64}
    melt_frac_total_vol_cum::Vector{Float64}
    effective_fractions::Vector{Float64}
    effective_fractions_vol::Vector{Float64}
    effect_fract_total::Vector{Float64}
    effect_fract_total_vol::Vector{Float64}
    T_solidus::Float64
end

"""
    SimulationResults

Complete results from a water-melt simulation.

# Fields  
- `config::PTPathConfig`: Configuration used
- `T_array::Vector{Float64}`: Temperature array
- `P_array::Vector{Float64}`: Pressure array
- `lithologies::Dict{String, LithologyResults}`: Results by lithology name
"""
struct SimulationResults
    config::PTPathConfig
    T_array::Vector{Float64}
    P_array::Vector{Float64}
    lithologies::Dict{String, LithologyResults}
end

# ============================================================================
# Default Compositions
# ============================================================================

"""
    default_compositions()

Returns default bulk compositions for metapelite, metagraywacke, and orthogneiss.

# References
- Metapelite: Forshaw and Pattison (2024)
- Metagraywacke: Villaros et al. (2018)  
- Orthogneiss: Weisbrod (1970)
"""
function default_compositions()
    return Dict(
        "metapelite" => BulkComposition(
            "metapelite",
            # SiO2, Al2O3, CaO, MgO, FeO, Fe2O3, K2O, Na2O, TiO2, MnO, H2O
            [59.22, 18.12, 0.60, 2.22, 6.32, 0.0, 3.65, 1.27, 0.84, 0.0, 20.0],
            "Forshaw and Pattison (2024)"
        ),
        "metagraywacke" => BulkComposition(
            "metagraywacke", 
            [70.48, 13.05, 1.68, 2.38, 4.86, 0.0, 2.43, 2.97, 0.70, 0.0, 20.0],
            "Villaros et al. (2018)"
        ),
        "orthogneiss" => BulkComposition(
            "orthogneiss",
            [71.93, 14.91, 1.28, 0.79, 2.40, 0.0, 4.73, 2.95, 0.26, 0.0, 20.0],
            "Weisbrod (1970)"
        )
    )
end

# ============================================================================
# Helper Functions  
# ============================================================================

"""
    create_pt_arrays(config::PTPathConfig)

Create temperature and pressure arrays from configuration.

# Arguments
- `config::PTPathConfig`: P-T path configuration

# Returns
- `T_array::Vector{Float64}`: Temperature array
- `P_array::Vector{Float64}`: Pressure array
"""
function create_pt_arrays(config::PTPathConfig)
    T_array = collect(config.Tmin:config.T_step:config.Tmax)
    n_steps = length(T_array)
    
    P_resolution = (config.Pmax - config.Pmin) / (n_steps - 1)
    if P_resolution == 0.0
        P_array = fill(config.Pmin, n_steps)
    else
        P_array = collect(config.Pmin:P_resolution:config.Pmax)
    end
    
    return T_array, P_array
end

"""
    prepare_composition(comp::BulkComposition)

Convert weight percent composition to normalized molar composition for MAGEMin.

# Arguments
- `comp::BulkComposition`: Bulk composition in weight percent

# Returns
- `X::Vector{Float64}`: Normalized molar composition
"""
function prepare_composition(comp::BulkComposition)
    X, _ = convertBulk4MAGEMin(comp.composition_wt, X_INIT_OX, "wt", "mp")
    return X ./ sum(X)
end

"""
    extract_water_excess!(X::Vector{Float64}, out, H2O_comp::Vector{Float64})

Remove excess free water from bulk composition based on minimization results.

# Arguments
- `X::Vector{Float64}`: Bulk composition (modified in place)
- `out`: MAGEMin output structure
- `H2O_comp::Vector{Float64}`: Pure H2O composition vector

# Returns
- `H2O_frac::Float64`: Molar fraction of removed water
- `H2O_frac_vol::Float64`: Volume fraction of removed water
"""
function extract_water_excess!(X::Vector{Float64}, out, H2O_comp::Vector{Float64})
    if "H2O" in out.ph
        H2O_index = findfirst(==("H2O"), out.ph)
        H2O_frac = out.ph_frac[H2O_index]
        H2O_frac_vol = out.ph_frac_vol[H2O_index]
        X .= X .- (H2O_comp .* H2O_frac)
        return H2O_frac, H2O_frac_vol
    else
        return 0.0, 0.0
    end
end

"""
    fractionate_garnet!(X::Vector{Float64}, out, g_factor::Float64)

Remove fractionated garnet from bulk composition.

# Arguments
- `X::Vector{Float64}`: Bulk composition (modified in place)
- `out`: MAGEMin output structure  
- `g_factor::Float64`: Fraction of garnet to remove (0.0-1.0)

# Returns
- `g_frac::Float64`: Molar fraction of fractionated garnet
- `g_frac_vol::Float64`: Volume fraction of fractionated garnet
"""
function fractionate_garnet!(X::Vector{Float64}, out, g_factor::Float64)
    if "g" in out.ph
        g_index = findfirst(==("g"), out.ph)
        garnet_comp = out.SS_vec[g_index].Comp
        g_frac = out.ph_frac[g_index] * g_factor
        g_frac_vol = out.ph_frac_vol[g_index] * g_factor
        X .= X .- (garnet_comp .* g_frac)
        return g_frac, g_frac_vol
    else
        return 0.0, 0.0
    end
end

"""
    find_solidus_temperature(out_array::Vector, T_array::Vector{Float64})

Find the solidus temperature (first appearance of melt).

# Arguments
- `out_array::Vector`: Array of MAGEMin output structures
- `T_array::Vector{Float64}`: Temperature array

# Returns
- `T_solidus::Float64`: Solidus temperature in °C
"""
function find_solidus_temperature(out_array::Vector, T_array::Vector{Float64})
    idx = findfirst(x -> "liq" in out_array[x].ph, eachindex(T_array))
    if isnothing(idx)
        return NaN
    end
    return T_array[idx]
end

"""
    initialize_lithology_results(name::String, n_steps::Int, X_init::Vector{Float64})

Initialize a LithologyResults structure with zeroed arrays.
"""
function initialize_lithology_results(name::String, n_steps::Int, X_init::Vector{Float64})
    return LithologyResults(
        name,
        Vector{MAGEMin_C.gmin_struct{Float64,Int64}}(undef, n_steps),
        copy(X_init),
        zeros(n_steps), zeros(n_steps),  # g_frac, g_frac_vol
        zeros(n_steps), zeros(n_steps),  # g_frac_total, g_frac_total_vol
        zeros(n_steps), zeros(n_steps),  # g_frac_total_cum, g_frac_total_vol_cum
        zeros(n_steps), zeros(n_steps),  # H2O_frac, H2O_frac_vol
        zeros(n_steps), zeros(n_steps),  # H2O_frac_total, H2O_frac_total_vol
        zeros(n_steps), zeros(n_steps),  # H2O_frac_total_cum, H2O_frac_total_vol_cum
        zeros(n_steps),                   # melt_frac_total_vol_cum
        zeros(n_steps), zeros(n_steps),  # effective_fractions, effective_fractions_vol
        zeros(n_steps), zeros(n_steps),  # effect_fract_total, effect_fract_total_vol
        0.0                               # T_solidus
    )
end

# ============================================================================
# Core Simulation Functions
# ============================================================================

"""
    run_lithology_along_path!(results::LithologyResults, P_array, T_array, 
                              data, config::PTPathConfig)

Run minimization along P-T path for a single lithology with garnet fractionation.

# Arguments
- `results::LithologyResults`: Results structure (modified in place)
- `P_array::Vector{Float64}`: Pressure array
- `T_array::Vector{Float64}`: Temperature array
- `data`: MAGEMin data structure
- `config::PTPathConfig`: Configuration
"""
function run_lithology_along_path!(results::LithologyResults, P_array::Vector{Float64}, 
                                   T_array::Vector{Float64}, data, config::PTPathConfig)
    n_steps = length(T_array)
    
    for step in 1:n_steps
        T = T_array[step]
        P = P_array[step]
        
        # Calculate effective fractions from previous step
        if step > 1
            results.effective_fractions[step] = 1 - results.g_frac[step-1] - results.H2O_frac[step-1]
            results.effective_fractions_vol[step] = 1 - results.g_frac_vol[step-1] - results.H2O_frac_vol[step-1]
        else
            results.effective_fractions[step] = 1.0
            results.effective_fractions_vol[step] = 1.0
        end
        
        # Run minimization
        results.out[step] = deepcopy(single_point_minimization(
            P, T, data, 
            X=results.X, Xoxides=XOXIDES, 
            name_solvus=true, sys_in=config.sys_in
        ))
        
        # Extract water excess
        results.H2O_frac[step], results.H2O_frac_vol[step] = extract_water_excess!(
            results.X, results.out[step], H2O_COMP
        )
        
        # Fractionate garnet
        results.g_frac[step], results.g_frac_vol[step] = fractionate_garnet!(
            results.X, results.out[step], config.g_factor
        )
    end
    
    # Calculate cumulative fractions
    results.effect_fract_total .= accumulate(*, results.effective_fractions)
    results.effect_fract_total_vol .= accumulate(*, results.effective_fractions_vol)
    
    for i in 1:n_steps
        results.g_frac_total[i] = results.g_frac[i] * results.effect_fract_total[i]
        results.g_frac_total_vol[i] = results.g_frac_vol[i] * results.effect_fract_total_vol[i]
        results.H2O_frac_total[i] = results.H2O_frac[i] * results.effect_fract_total[i]
        results.H2O_frac_total_vol[i] = results.H2O_frac_vol[i] * results.effect_fract_total_vol[i]
        
        if "liq" in results.out[i].ph
            liq_idx = findfirst(==("liq"), results.out[i].ph)
            results.melt_frac_total_vol_cum[i] = results.out[i].ph_frac_vol[liq_idx] * results.effect_fract_total_vol[i]
        end
    end
    
    results.g_frac_total_cum .= accumulate(+, results.g_frac_total)
    results.g_frac_total_vol_cum .= accumulate(+, results.g_frac_total_vol)
    results.H2O_frac_total_cum .= accumulate(+, results.H2O_frac_total)
    results.H2O_frac_total_vol_cum .= accumulate(+, results.H2O_frac_total_vol)
    
    results.T_solidus = find_solidus_temperature(results.out, T_array)
end

"""
    run_water_fluxed_melting!(results::LithologyResults, water_source::LithologyResults,
                              og_init::LithologyResults, P_array, T_array, data, config;
                              water_ratio::Float64=1.0)

Run water-fluxed melting simulation for orthogneiss receiving water from another lithology.

# Arguments
- `results::LithologyResults`: Results structure for water-fluxed system (modified in place)
- `water_source::LithologyResults`: Source lithology providing water
- `og_init::LithologyResults`: Initial orthogneiss results (for solidus check)
- `P_array`, `T_array`: P-T arrays
- `data`: MAGEMin data structure
- `config::PTPathConfig`: Configuration
- `water_ratio::Float64`: Ratio of water source to orthogneiss (default 1.0)
"""
function run_water_fluxed_melting!(results::LithologyResults, water_source::LithologyResults,
                                   og_init::LithologyResults, P_array::Vector{Float64}, 
                                   T_array::Vector{Float64}, data, config::PTPathConfig;
                                   water_ratio::Float64=1.0)
    n_steps = length(T_array)
    Tmin = config.Tmin
    T_step = config.T_step
    
    for step in 1:n_steps
        T = T_array[step]
        P = P_array[step]
        
        # Add external water only when orthogneiss has melt
        if "liq" in og_init.out[step].ph
            external_H2O = water_ratio * water_source.H2O_frac_total[step]
            external_H2O_vol = water_ratio * water_source.H2O_frac_total_vol[step]
        else
            external_H2O = 0.0
            external_H2O_vol = 0.0
        end
        
        results.X .= results.X .+ external_H2O .* H2O_COMP
        
        # Calculate effective fractions
        if step > 1
            results.effective_fractions[step] = 1 - results.g_frac[step-1] + external_H2O
            results.effective_fractions_vol[step] = 1 - results.g_frac_vol[step-1] + external_H2O_vol
        else
            results.effective_fractions[step] = 1.0
            results.effective_fractions_vol[step] = 1.0
        end
        
        # Run minimization
        results.out[step] = deepcopy(single_point_minimization(
            P, T, data,
            X=results.X, Xoxides=XOXIDES,
            name_solvus=true, sys_in=config.sys_in
        ))
        
        # Fractionate garnet
        results.g_frac[step], results.g_frac_vol[step] = fractionate_garnet!(
            results.X, results.out[step], config.g_factor
        )
        
        # Extract water excess only before orthogneiss solidus
        solidus_step = Int((og_init.T_solidus - Tmin) / T_step) + 1
        if "H2O" in results.out[step].ph && step < solidus_step
            H2O_index = findfirst(==("H2O"), results.out[step].ph)
            results.H2O_frac[step] = results.out[step].ph_frac[H2O_index]
            results.H2O_frac_vol[step] = results.out[step].ph_frac_vol[H2O_index]
            results.X .= results.X .- (H2O_COMP .* results.out[step].ph_frac[H2O_index])
        else
            results.H2O_frac[step] = 0.0
            results.H2O_frac_vol[step] = 0.0
        end
    end
    
    # Calculate cumulative fractions
    results.effect_fract_total .= accumulate(*, results.effective_fractions)
    results.effect_fract_total_vol .= accumulate(*, results.effective_fractions_vol)
    
    for i in 1:n_steps
        results.g_frac_total[i] = results.g_frac[i] * results.effect_fract_total[i]
        results.g_frac_total_vol[i] = results.g_frac_vol[i] * results.effect_fract_total_vol[i]
        results.H2O_frac_total[i] = results.H2O_frac[i] * results.effect_fract_total[i]
        results.H2O_frac_total_vol[i] = results.H2O_frac_vol[i] * results.effect_fract_total_vol[i]
        
        if "liq" in results.out[i].ph
            liq_idx = findfirst(==("liq"), results.out[i].ph)
            results.melt_frac_total_vol_cum[i] = results.out[i].ph_frac_vol[liq_idx] * results.effect_fract_total_vol[i]
        end
    end
    
    results.g_frac_total_cum .= accumulate(+, results.g_frac_total)
    results.g_frac_total_vol_cum .= accumulate(+, results.g_frac_total_vol)
    results.H2O_frac_total_cum .= accumulate(+, results.H2O_frac_total)
    results.H2O_frac_total_vol_cum .= accumulate(+, results.H2O_frac_total_vol)
    
    results.T_solidus = find_solidus_temperature(results.out, T_array)
end

"""
    remove_initial_water_excess!(X::Vector{Float64}, P::Float64, T::Float64, data, config)

Remove excess water from initial composition at starting P-T conditions.

# Returns
- `X::Vector{Float64}`: Modified composition with water excess removed
"""
function remove_initial_water_excess!(X::Vector{Float64}, P::Float64, T::Float64, data, config::PTPathConfig)
    out_init = single_point_minimization(P, T, data, X=X, Xoxides=XOXIDES, name_solvus=true, sys_in=config.sys_in)
    
    if "H2O" in out_init.ph
        H2O_index = findfirst(==("H2O"), out_init.ph)
        X .= X .- (H2O_COMP .* out_init.ph_frac[H2O_index])
    else
        error("Initial composition must contain H2O to proceed with water-excess calculation.")
    end
    
    X .= X ./ sum(X)
    return X
end

# ============================================================================
# Main Simulation Function
# ============================================================================

"""
    run_water_melt_simulation(config::PTPathConfig=PTPathConfig(); 
                              compositions::Dict=default_compositions(),
                              mp_og_ratio::Float64=1.0,
                              mg_og_ratio::Float64=1.0,
                              verbose::Bool=true)

Run the complete water-melt simulation along a P-T path.

# Arguments
- `config::PTPathConfig`: P-T path configuration
- `compositions::Dict`: Dictionary of BulkComposition objects
- `mp_og_ratio::Float64`: Ratio of metapelite to orthogneiss for water flux (default 1.0)
- `mg_og_ratio::Float64`: Ratio of metagraywacke to orthogneiss for water flux (default 1.0)
- `verbose::Bool`: Print progress information (default true)

# Returns
- `SimulationResults`: Complete simulation results
"""
function run_water_melt_simulation(config::PTPathConfig=PTPathConfig();
                                   compositions::Dict=default_compositions(),
                                   mp_og_ratio::Float64=1.0,
                                   mg_og_ratio::Float64=1.0,
                                   verbose::Bool=true)
    
    # Initialize MAGEMin
    data = Initialize_MAGEMin("mp", verbose=false)
    
    try
        # Create P-T arrays
        T_array, P_array = create_pt_arrays(config)
        n_steps = length(T_array)
        
        verbose && println("Running simulation with $(n_steps) steps from $(config.Tmin)°C to $(config.Tmax)°C")
        
        # Prepare compositions
        X_mp = prepare_composition(compositions["metapelite"])
        X_mg = prepare_composition(compositions["metagraywacke"])
        X_og = prepare_composition(compositions["orthogneiss"])
        
        # Remove initial water excess
        remove_initial_water_excess!(X_mp, P_array[1], T_array[1], data, config)
        remove_initial_water_excess!(X_mg, P_array[1], T_array[1], data, config)
        remove_initial_water_excess!(X_og, P_array[1], T_array[1], data, config)
        
        # Initialize results structures
        results_mp = initialize_lithology_results("metapelite", n_steps, X_mp)
        results_mg = initialize_lithology_results("metagraywacke", n_steps, X_mg)
        results_og = initialize_lithology_results("orthogneiss", n_steps, X_og)
        
        # Run main lithologies
        verbose && println("Computing metapelite, metagraywacke, and orthogneiss...")
        @showprogress enabled=verbose for (results, name) in [(results_mp, "mp"), (results_mg, "mg"), (results_og, "og")]
            run_lithology_along_path!(results, P_array, T_array, data, config)
        end
        
        verbose && println("\nSolidus temperatures:")
        verbose && println("  Metapelite: $(results_mp.T_solidus)°C")
        verbose && println("  Metagraywacke: $(results_mg.T_solidus)°C")
        verbose && println("  Orthogneiss: $(results_og.T_solidus)°C")
        
        # Calculate water released between solidi
        if !isnan(results_mp.T_solidus) && !isnan(results_og.T_solidus)
            idx_mp = Int((results_mp.T_solidus - config.Tmin) / config.T_step) + 1
            idx_og = Int((results_og.T_solidus - config.Tmin) / config.T_step) + 1
            H2O_mp_released = results_mp.H2O_frac_total_cum[idx_mp] - results_mp.H2O_frac_total_cum[idx_og]
            verbose && println("\nH2O released (MP, between solidi): $(round(H2O_mp_released*100, digits=2))%")
        end
        
        if !isnan(results_mg.T_solidus) && !isnan(results_og.T_solidus)
            idx_mg = Int((results_mg.T_solidus - config.Tmin) / config.T_step) + 1
            idx_og = Int((results_og.T_solidus - config.Tmin) / config.T_step) + 1
            H2O_mg_released = results_mg.H2O_frac_total_cum[idx_mg] - results_mg.H2O_frac_total_cum[idx_og]
            verbose && println("H2O released (MG, between solidi): $(round(H2O_mg_released*100, digits=2))%")
        end
        
        # Water-fluxed melting simulations
        X_og_mp = prepare_composition(compositions["orthogneiss"])
        X_og_mg = prepare_composition(compositions["orthogneiss"])
        remove_initial_water_excess!(X_og_mp, P_array[1], T_array[1], data, config)
        remove_initial_water_excess!(X_og_mg, P_array[1], T_array[1], data, config)
        
        results_og_mp = initialize_lithology_results("orthogneiss+mp_water", n_steps, X_og_mp)
        results_og_mg = initialize_lithology_results("orthogneiss+mg_water", n_steps, X_og_mg)
        
        verbose && println("\nComputing water-fluxed orthogneiss (MP water source)...")
        run_water_fluxed_melting!(results_og_mp, results_mp, results_og, P_array, T_array, data, config; 
                                  water_ratio=mp_og_ratio)
        
        verbose && println("Computing water-fluxed orthogneiss (MG water source)...")
        run_water_fluxed_melting!(results_og_mg, results_mg, results_og, P_array, T_array, data, config;
                                  water_ratio=mg_og_ratio)
        
        # Compile results
        lithologies = Dict(
            "metapelite" => results_mp,
            "metagraywacke" => results_mg,
            "orthogneiss" => results_og,
            "orthogneiss+mp_water" => results_og_mp,
            "orthogneiss+mg_water" => results_og_mg
        )
        
        return SimulationResults(config, T_array, P_array, lithologies)
        
    finally
        Finalize_MAGEMin(data)
    end
end

# ============================================================================
# Plotting Functions
# ============================================================================

"""
    plot_water_release(results::SimulationResults; save_path::String="")

Plot cumulative water release for all lithologies.
"""
function plot_water_release(results::SimulationResults; save_path::String="")
    fig = Figure(size=(800, 700))
    ax = Axis(fig[1, 1], xlabel="T (°C)", ylabel="H₂O fraction", title="Cumulative Water Release")
    
    T = results.T_array
    
    lines!(ax, T, results.lithologies["metapelite"].H2O_frac_total_cum, color=:black, label="Metapelite")
    lines!(ax, T, results.lithologies["metagraywacke"].H2O_frac_total_cum, color=:blue, label="Metagraywacke")
    
    # Solidus lines
    vlines!(ax, results.lithologies["metapelite"].T_solidus, linestyle=:dash, color=:black, linewidth=1.5)
    vlines!(ax, results.lithologies["metagraywacke"].T_solidus, linestyle=:dash, color=:blue, linewidth=1.5)
    vlines!(ax, results.lithologies["orthogneiss"].T_solidus, linestyle=:dash, color=:red, linewidth=1.5, label="OG solidus")
    
    axislegend(ax, position=:lt)
    
    if !isempty(save_path)
        save(save_path, fig)
    end
    
    return fig
end

"""
    plot_melt_fractions(results::SimulationResults; save_path::String="")

Plot melt fractions for all lithologies including water-fluxed systems.
"""
function plot_melt_fractions(results::SimulationResults; save_path::String="")
    fig = Figure(size=(800, 700))
    ax = Axis(fig[1, 1], xlabel="T (°C)", ylabel="Melt fraction (vol)", title="Melt Production")
    
    T = results.T_array
    
    lines!(ax, T, results.lithologies["orthogneiss+mp_water"].melt_frac_total_vol_cum, color=:red, label="OG + H₂O (MP)")
    lines!(ax, T, results.lithologies["orthogneiss+mg_water"].melt_frac_total_vol_cum, color=:purple, label="OG + H₂O (MG)")
    lines!(ax, T, results.lithologies["metapelite"].melt_frac_total_vol_cum, color=:black, label="Metapelite")
    lines!(ax, T, results.lithologies["metagraywacke"].melt_frac_total_vol_cum, color=:blue, label="Metagraywacke")
    lines!(ax, T, results.lithologies["orthogneiss"].melt_frac_total_vol_cum, color=:pink, label="Orthogneiss")
    
    # Solidus lines
    vlines!(ax, results.lithologies["metapelite"].T_solidus, linestyle=:dash, color=:black, linewidth=1.5)
    vlines!(ax, results.lithologies["metagraywacke"].T_solidus, linestyle=:dash, color=:blue, linewidth=1.5)
    vlines!(ax, results.lithologies["orthogneiss"].T_solidus, linestyle=:dash, color=:red, linewidth=1.5)
    
    axislegend(ax, position=:lt)
    
    if !isempty(save_path)
        save(save_path, fig)
    end
    
    return fig
end

"""
    get_phase_color(phase::String)

Get color for a phase, using default if not in color map.
"""
function get_phase_color(phase::String)
    if haskey(PHASE_COLOR_MAP, phase)
        return PHASE_COLOR_MAP[phase]
    else
        # Default colors for unknown phases
        default_colors = [:lightcoral, :lightskyblue, :lightgreen, :lightsalmon, :plum]
        return default_colors[mod1(hash(phase), length(default_colors))]
    end
end

"""
    plot_mode_boxes(results::SimulationResults; save_path::String="")

Create mode box diagrams for all lithologies.
"""
function plot_mode_boxes(results::SimulationResults; save_path::String="")
    fig = Figure(size=(1800, 1400))
    
    lithology_keys = ["metapelite", "metagraywacke", "orthogneiss+mp_water", "orthogneiss+mg_water"]
    titles = ["Metapelite", "Metagraywacke", "Orthogneiss + H₂O (MP)", "Orthogneiss + H₂O (MG)"]
    positions = [(1,1), (1,2), (2,1), (2,2)]
    
    T = results.T_array
    n_steps = length(T)
    all_phases = String[]
    
    for (key, title, pos) in zip(lithology_keys, titles, positions)
        lith = results.lithologies[key]
        ax = Axis(fig[pos...], xlabel="T (°C)", ylabel="Phase Fraction (molar)", title="Mode Box - $title")
        xlims!(ax, minimum(T), maximum(T))
        ylims!(ax, 0.0, 1.0)
        
        # Collect phases (excluding garnet)
        phase_names = unique(vcat([lith.out[i].ph for i in 1:n_steps]...))
        phase_names = filter(x -> x != "g", phase_names)
        append!(all_phases, phase_names)
        
        # Calculate phase fractions
        n_phases = length(phase_names)
        phase_fractions = zeros(n_steps, n_phases)
        
        for i in 1:n_steps
            for (j, phase) in enumerate(lith.out[i].ph)
                if phase != "g"
                    phase_idx = findfirst(==(phase), phase_names)
                    if !isnothing(phase_idx)
                        phase_fractions[i, phase_idx] = lith.out[i].ph_frac[j] * lith.effect_fract_total[i]
                    end
                end
            end
        end
        
        # Plot bands
        cumulative = zeros(n_steps)
        for j in 1:n_phases
            lower = copy(cumulative)
            cumulative .+= phase_fractions[:, j]
            band!(ax, T, lower, cumulative, color=get_phase_color(phase_names[j]))
        end
        
        # Add fractionated garnet
        lower_g = copy(cumulative)
        cumulative_g = cumulative .+ lith.g_frac_total_cum
        band!(ax, T, lower_g, cumulative_g, color=:indianred)
        
        # Add fractionated H2O
        lower_h2o = copy(cumulative_g)
        cumulative_h2o = cumulative_g .+ lith.H2O_frac_total_cum
        band!(ax, T, lower_h2o, cumulative_h2o, color=:lightblue)
        
        # Solidus lines
        vlines!(ax, results.lithologies["orthogneiss"].T_solidus, linestyle=:dash, color=:red, linewidth=1.5)
    end
    
    # Create legend
    unique_phases = unique(all_phases)
    legend_elements = LegendElement[PolyElement(color=get_phase_color(p)) for p in sort(unique_phases)]
    legend_labels = sort(unique_phases)
    
    push!(legend_elements, PolyElement(color=:indianred))
    push!(legend_labels, "g (fractionated)")
    push!(legend_elements, PolyElement(color=:lightblue))
    push!(legend_labels, "H₂O (fractionated)")
    push!(legend_elements, LineElement(color=:red, linestyle=:dash))
    push!(legend_labels, "OG solidus")
    
    Legend(fig[1:2, 3], legend_elements, legend_labels)
    
    if !isempty(save_path)
        save(save_path, fig)
    end
    
    return fig
end

"""
    plot_results(results::SimulationResults; output_dir::String=".")

Generate and save all plots.
"""
function plot_results(results::SimulationResults; output_dir::String=".")
    fig1 = plot_water_release(results; save_path=joinpath(output_dir, "H2O_frac.svg"))
    fig2 = plot_melt_fractions(results; save_path=joinpath(output_dir, "Liq_frac.svg"))
    fig3 = plot_mode_boxes(results; save_path=joinpath(output_dir, "Mode_Boxes.svg"))
    
    display(fig1)
    display(fig2)
    display(fig3)
    
    return (fig1, fig2, fig3)
end

# ============================================================================
# Data Export Functions
# ============================================================================

"""
    save_results(results::SimulationResults, output_dir::String=".")

Save simulation results to CSV files.
"""
function save_results(results::SimulationResults, output_dir::String=".")
    # Create summary DataFrame
    df = DataFrame(
        T = results.T_array,
        P = results.P_array
    )
    
    for (name, lith) in results.lithologies
        short_name = replace(name, "+" => "_", " " => "_")
        df[!, "$(short_name)_H2O_cum"] = lith.H2O_frac_total_cum
        df[!, "$(short_name)_g_cum"] = lith.g_frac_total_cum
        df[!, "$(short_name)_melt_vol"] = lith.melt_frac_total_vol_cum
    end
    
    CSV.write(joinpath(output_dir, "simulation_results.csv"), df)
    
    println("Results saved to $(joinpath(output_dir, "simulation_results.csv"))")
    
    return df
end

end # module
