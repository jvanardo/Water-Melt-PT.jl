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

"""
Define the parameters use to caluclate the P-T path

g_factor: fraction of garnet used for the fractionation (between 0 and 1, default = 1.0)

Tmin: minimum temperature of the P-T path (°C, default = 600.0)
Tmax: maximum temperature of the P-T path (°C, default = 750.0)
Pmin: minimum pressure of the P-T path (GPa, default = 8.0)
Pmax: maximum pressure of the P-T path (GPa, default = 12.5)

og_fraction: fraction of orthogneiss compare to metapelite in the system considered (between 0 and 1, default = 0.25)

X_init_mp_wt: bulk composition of the metapelite in wt% (default = [59.22; 18.12; 0.60; 2.22; 6.32; 0.0; 3.65; 1.27; 0.84; 0.0; 7.5])
X_init_mg_wt: bulk composition of the metagraywacke in wt% (default = [70.48; 13.05; 1.68; 2.38; 4.86; 0.0; 2.43; 2.97; 0.70; 0.0; 7.5])
X_init_og_wt: bulk composition of the orthogneiss in wt% (default = [71.93; 14.91; 1.28; 0.79; 2.40; 0.0; 4.73; 2.95; 0.26; 0.0; 7.5])
"""

g_factor = 1.0

Tmin = 600.0
Tmax = 750.0
Pmin = 8.0
Pmax = 12.5

og_fraction = 0.25 # change this value to modify the fraction of orthogneiss compare to metapelite

# Define the bulk composi and lithologies
X_init_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "MnO"; "H2O"]
X_init_mp_wt = [59.22; 18.12; 0.60; 2.22; 6.32; 0.0; 3.65; 1.27; 0.84; 0.0; 7.5] #bulk compo of metapelites (Forshaw and Pattison, 2024), Fe2O3 and MnO = 0.0
X_init_mg_wt = [70.48; 13.05; 1.68; 2.38; 4.86; 0.0; 2.43; 2.97; 0.70; 0.0; 7.5] #bulk compo of metagraywackes (Villaros et al., 2018)
X_init_og_wt = [71.93; 14.91; 1.28; 0.79; 2.40; 0.0; 4.73; 2.95; 0.26; 0.0; 7.5] #bulk compo of orthogneiss (Weisbrod, 1970)


"""
Initialized MAGEMin_C, define the PT path and convert the bulk composi from wt% to mol% and normalize them.
Then, run the first point of the path to extract the excess water and remove it from the bulk composition for the next steps.
"""

T_array, P_array = define_PT_path(Tmin, Tmax, Pmin, Pmax)

data = Initialize_MAGEMin("mp", verbose=false);
Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "H2O"];
sys_in = "mol"

X_mp, X_mp_wt = prepare_bulk_composition(X_init_mp_wt, X_init_ox, P_array, T_array, data, sys_in)
X_mg, X_mg_wt = prepare_bulk_composition(X_init_mg_wt, X_init_ox, P_array, T_array, data, sys_in)
X_og, X_og_wt = prepare_bulk_composition(X_init_og_wt, X_init_ox, P_array, T_array, data, sys_in)

"""
Run the calculation along the PT path to obtain the amount of water and melt generated in metasedimentary lithologies, and in orthogneiss without fluid influx
"""

out_mp, g_frac_mp_total_cum, H2O_frac_total_mp, H2O_frac_mp_total_cum, melt_frac_mp, T_solidus_mp, effect_frac_mp_total = run_PT_path(T_array, P_array, data, X_mp, Xoxides, g_factor, sys_in, Tmin);
out_mg, g_frac_mg_total_cum, H2O_frac_total_mg, H2O_frac_mg_total_cum, melt_frac_mg, T_solidus_mg, effect_frac_mg_total = run_PT_path(T_array, P_array, data, X_mg, Xoxides, g_factor, sys_in, Tmin);
out_og, g_frac_og_total_cum, H2O_frac_total_og, H2O_frac_og_total_cum, melt_frac_og, T_solidus_og, effect_frac_og_total = run_PT_path(T_array, P_array, data, X_og, Xoxides, g_factor, sys_in, Tmin);

# Calculate the amount of water released between the two solidi
H2O_released_mp = zeros(3)
H2O_released_mg = zeros(3)
H2O_released_mp .= H2O_frac_mp_total_cum[Int(T_solidus_mp-Tmin),:] .- H2O_frac_mp_total_cum[Int(T_solidus_og-Tmin),:]
H2O_released_mg .= H2O_frac_mg_total_cum[Int(T_solidus_mg-Tmin),:] .- H2O_frac_mg_total_cum[Int(T_solidus_og-Tmin),:]


"""
Calculation of the volume of melt created with the water-influx
"""

additional_H2O_mp = calculate_additional_H2O(T_array, out_og, og_fraction, H2O_frac_total_mp)
additional_H2O_mg = calculate_additional_H2O(T_array, out_og, og_fraction, H2O_frac_total_mg)

X_og_melt, X_og_melt_wt = prepare_bulk_composition(X_init_og_wt, X_init_ox, P_array, T_array, data, sys_in)
X_og_mp_melt = deepcopy(X_og_melt)
X_og_mg_melt = deepcopy(X_og_melt)
out_og_mp_melt, g_frac_og_mp_melt_total_cum, H2O_frac_og_mp_melt_total, H2O_frac_og_mp_melt_total_cum, melt_frac_og_mp, T_solidus_og_mp_melt, effect_frac_og_mp_melt_total = run_PT_path(T_array, P_array, data, X_og_mp_melt, Xoxides, g_factor, sys_in, Tmin, additional_H2O = additional_H2O_mp);
out_og_mg_melt, g_frac_og_mg_melt_total_cum, H2O_frac_og_mg_melt_total, H2O_frac_og_mg_melt_total_cum, melt_frac_og_mg, T_solidus_og_mg_melt, effect_frac_og_mg_melt_total = run_PT_path(T_array, P_array, data, X_og_mg_melt, Xoxides, g_factor, sys_in, Tmin, additional_H2O = additional_H2O_mg);

Finalize_MAGEMin(data)

# search for the temperature where metasedimentary rocks become more melted than orthogneiss+water
idx_mp_exceeds = findfirst(melt_frac_mp[:,3] .> melt_frac_og_mp[:,3])
idx_mg_exceeds = findfirst(melt_frac_mg[:,3] .> melt_frac_og_mg[:,3])

# If MP never exceeds OG+MP, set to last index (750°C)
if isnothing(idx_mp_exceeds)
    idx_mp_exceeds = length(T_array)
end
if isnothing(idx_mg_exceeds)
    idx_mg_exceeds = length(T_array)
end

T_og_mp = T_array[idx_mp_exceeds] - T_solidus_og
T_og_mg = T_array[idx_mg_exceeds] - T_solidus_og

# Print results
println("H2O released mp: $(H2O_released_mp[1]*100) mol%, $(H2O_released_mp[2]*100) vol%, $(H2O_released_mp[3]*100) wt%")
println("H2O released mg: $(H2O_released_mg[1]*100) mol%, $(H2O_released_mg[2]*100) vol%, $(H2O_released_mg[3]*100) wt%")
println("$T_solidus_mp °C metapelite solidus temperature")
println("$T_solidus_mg °C metagraywacke solidus temperature")
println("$T_solidus_og °C orthogneiss solidus temperature")
println("T og+mp melt = $(T_array[idx_mp_exceeds]) °C")
println("T og+mg melt = $(T_array[idx_mg_exceeds]) °C")
println("melt_fraction_og_mp = $(melt_frac_og_mp[idx_mp_exceeds,3]*100)")
println("melt_fraction_og_mg = $(melt_frac_og_mg[idx_mg_exceeds,3]*100)")
println("melt_fraction_og_mg at T_solidus_mp = $(melt_frac_og_mg[Int(T_solidus_mp-Tmin),3]*100) %")
println("melt_fraction_og_mp at T_solidus_mg = $(melt_frac_og_mp[Int(T_solidus_mg-Tmin),3]*100) %")


"""
Plot the results

fig_pl shows the evolution of the Xan content in plagioclase along PT path
fig_melt_comp shows the evolution of oxides content in the melt along the PT path for the three lithologies, with and without fluid transfer
fig1 shows the amount of water released in wt% along the PT path for the metasedimentary rocks
fig2 shows the amount of melt generated in vol% along the PT path for the metasedimentary rocks and the orthogneiss with and without fluid transfer
fig3 shows boxplots for the three lithologies
"""

fig_pl = plot_fig_pl(T_array, out_mp, out_mg, out_og, T_solidus_mp, T_solidus_mg, T_solidus_og)
fig_melt_comp = plot_fig_melt_comp(T_array, out_mp, out_mg, out_og, out_og_mp_melt, out_og_mg_melt)
fig1 = plot_fig1(T_array, H2O_frac_mp_total_cum, H2O_frac_mg_total_cum, T_solidus_mp, T_solidus_mg, T_solidus_og, H2O_released_mp, H2O_released_mg)
fig2 = plot_fig2(T_array, melt_frac_og_mp[:,3], melt_frac_og_mg[:,3], melt_frac_mp[:,3], melt_frac_mg[:,3], melt_frac_og[:,3], T_solidus_mp, T_solidus_mg, T_solidus_og)
fig3 = plot_fig3(T_array, out_mp, out_mg, out_og, out_og_mp_melt, out_og_mg_melt, effect_frac_mp_total, effect_frac_mg_total, effect_frac_og_mp_melt_total, effect_frac_og_mg_melt_total, g_frac_mp_total_cum, g_frac_mg_total_cum, g_frac_og_total_cum, H2O_frac_mp_total_cum, H2O_frac_mg_total_cum, H2O_frac_og_total_cum, melt_frac_mp, melt_frac_mg, melt_frac_og, T_solidus_mp, T_solidus_mg, T_solidus_og, g_frac_og_mp_melt_total_cum, g_frac_og_mg_melt_total_cum, H2O_frac_og_mp_melt_total_cum, H2O_frac_og_mg_melt_total_cum, melt_frac_og_mp, melt_frac_og_mg)
save("H2O_frac_PT.svg", fig1)
save("An_fraction_PT.svg", fig_pl)
save("Melt_comp_PT.svg", fig_melt_comp)
save("Liq_frac_PT.svg",  fig2)
save("Mode_Boxes_4Panel_PT.svg", fig3)
