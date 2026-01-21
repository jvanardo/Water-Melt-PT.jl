#=
    Example: Water-Melt-PT Simulation

    This script demonstrates how to use the WaterMeltPT package to model
    water release and melt production during crustal metamorphism.

    Prerequisites:
    1. Install Julia 1.9 or higher
    2. Install dependencies: ] add MAGEMin_C CairoMakie CSV DataFrames ProgressMeter Reexport
    3. Or activate this project: ] activate . ; instantiate
=#

# # Load the package
# using Pkg
# Pkg.activate(dirname(@__DIR__))  # Activate the package directory
# Pkg.instantiate()                 # Install dependencies if needed

using WaterMeltPT

# ============================================================================
# Example 1: Run simulation with default settings
# ============================================================================

println("=" ^ 60)
println("Example 1: Default simulation")
println("=" ^ 60)

# Run with all defaults (Tmin=600°C, Tmax=750°C, Pmin=14kbar, Pmax=9.5kbar)
results = run_water_melt_simulation()

# Generate and save plots
plot_results(results; output_dir=".")

# Save results to CSV
save_results(results, ".")

# ============================================================================
# Example 2: Custom P-T path configuration
# ============================================================================

println("\n" * "=" ^ 60)
println("Example 2: Custom P-T path")
println("=" ^ 60)

# Define a custom P-T path
custom_config = PTPathConfig(
    Tmin = 550.0,      # Starting temperature (°C)
    Tmax = 800.0,      # Ending temperature (°C)
    T_step = 2.0,      # Temperature increment (°C)
    Pmin = 12.0,       # Starting pressure (kbar)
    Pmax = 8.0,        # Ending pressure (kbar)
    g_factor = 0.5     # Fractionate 50% of garnet
)

results_custom = run_water_melt_simulation(custom_config)

# Access specific results
mp_results = results_custom.lithologies["metapelite"]
println("\nMetapelite solidus: $(mp_results.T_solidus)°C")
println("Final cumulative water released: $(round(mp_results.H2O_frac_total_cum[end]*100, digits=2))%")
println("Final cumulative garnet: $(round(mp_results.g_frac_total_cum[end]*100, digits=2))%")

# ============================================================================
# Example 3: Custom bulk compositions
# ============================================================================

println("\n" * "=" ^ 60)
println("Example 3: Custom bulk compositions")
println("=" ^ 60)

# Define your own bulk compositions
my_compositions = Dict(
    "metapelite" => BulkComposition(
        "my_metapelite",
        # SiO2, Al2O3, CaO, MgO, FeO, Fe2O3, K2O, Na2O, TiO2, MnO, H2O
        [58.0, 19.0, 0.5, 2.5, 7.0, 0.0, 4.0, 1.0, 0.9, 0.0, 7.0],
        "Custom composition"
    ),
    "metagraywacke" => BulkComposition(
        "my_metagraywacke",
        [71.0, 12.5, 2.0, 2.5, 5.0, 0.0, 2.5, 3.0, 0.6, 0.0, 7.0],
        "Custom composition"
    ),
    "orthogneiss" => BulkComposition(
        "my_orthogneiss",
        [72.0, 15.0, 1.2, 0.8, 2.5, 0.0, 4.5, 3.0, 0.3, 0.0, 7.0],
        "Custom composition"
    )
)

results_custom_comp = run_water_melt_simulation(
    PTPathConfig();  # Use default P-T path
    compositions = my_compositions,
    mp_og_ratio = 1.5,   # 1.5:1 ratio of metapelite to orthogneiss
    mg_og_ratio = 2.0    # 2:1 ratio of metagraywacke to orthogneiss
)

# ============================================================================
# Example 4: Access detailed results
# ============================================================================

println("\n" * "=" ^ 60)
println("Example 4: Accessing detailed results")
println("=" ^ 60)

# Get the default simulation results
results = run_water_melt_simulation(verbose=false)

# Access temperature and pressure arrays
T = results.T_array
P = results.P_array
println("P-T path: $(length(T)) steps from $(T[1])°C to $(T[end])°C")

# Access individual lithology results
og_mp = results.lithologies["orthogneiss+mp_water"]

# Find temperature where water-fluxed OG produces more melt than dry metapelite
mp = results.lithologies["metapelite"]
idx_crossover = findfirst(mp.melt_frac_total_vol_cum .> og_mp.melt_frac_total_vol_cum)

if !isnothing(idx_crossover)
    println("Metapelite exceeds OG+water melt at T = $(T[idx_crossover])°C")
end

# Access phase assemblage at a specific temperature step
step_idx = 50  # Step 50
if step_idx <= length(mp.out)
    phases = mp.out[step_idx].ph
    fractions = mp.out[step_idx].ph_frac
    println("\nPhase assemblage at T=$(T[step_idx])°C:")
    for (ph, frac) in zip(phases, fractions)
        println("  $ph: $(round(frac*100, digits=1))%")
    end
end

println("\n" * "=" ^ 60)
println("Examples complete!")
println("=" ^ 60)
