"""
    Water-Melt-PT

A Julia package for modeling water release and melt production during 
metamorphism using MAGEMin thermodynamic calculations.


"""

using Reexport: @reexport
@reexport using MAGEMin_C
using CSV
using DataFrames
using ProgressMeter
using CairoMakie


data = Initialize_MAGEMin("mp", verbose=false);
Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "H2O"];
sys_in = "mol"


# Define the fraction of garnet used for the fractionation
g_factor = 1.0

# Define the P-T path
Tmin = 600.0
Tmax = 750.0
T_array = collect(Tmin:1.0:Tmax)
Pmin = 14.0
Pmax = 9.5
P_resolution = (Pmax - Pmin) / (length(T_array) - 1)
if P_resolution == 0.0
    P_array = Pmin * (ones(length(T_array)))
else
    P_array = collect(Pmin:P_resolution:Pmax)
end


# Define the bulk compositions and lithologies
X_init_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "MnO"; "H2O"]
X_init_mp_wt = [59.22; 18.12; 0.60; 2.22; 6.32; 0.0; 3.65; 1.27; 0.84; 0.0; 7.5] #bulk compo of metapelites (Forshaw and Pattison, 2024), Fe2O3 and MnO = 0.0 
X_init_mg_wt = [70.48; 13.05; 1.68; 2.38; 4.86; 0.0; 2.43; 2.97; 0.70; 0.0; 7.5] #bulk compo of metagraywackes (Villaros et al., 2018)
X_init_og_wt = [71.93; 14.91; 1.28; 0.79; 2.40; 0.0; 4.73; 2.95; 0.26; 0.0; 7.5] #bulk compo of orthogneiss (Weisbrod, 1970)
# Concert wt% to mol%
X_init_mp, ox = convertBulk4MAGEMin(X_init_mp_wt, X_init_ox, "wt", "mp")
X_init_mg, ox = convertBulk4MAGEMin(X_init_mg_wt, X_init_ox, "wt", "mp")
X_init_og, ox = convertBulk4MAGEMin(X_init_og_wt, X_init_ox, "wt", "mp")
X_init_mp = X_init_mp ./ sum(X_init_mp) # normalize to 1
X_init_mg = X_init_mg ./ sum(X_init_mg) # normalize to 1
X_init_og = X_init_og ./ sum(X_init_og) # normalize to 1
X_mp = copy(X_init_mp) # copy the initial composition for metapelite
X_mg = copy(X_init_mg) # copy the initial composition for metagraywacke
X_og = copy(X_init_og) # copy the initial composition for orthogneiss

# Define water composition to extract the excess water
H2O_comp = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 1.0]

# Calculation for the first P-T conditions
out_mp_init = single_point_minimization(P_array[1], T_array[1], data, X=X_init_mp, Xoxides=Xoxides, name_solvus=true, sys_in=sys_in)
out_mg_init = single_point_minimization(P_array[1], T_array[1], data, X=X_init_mg, Xoxides=Xoxides, name_solvus=true, sys_in=sys_in)
out_og_init = single_point_minimization(P_array[1], T_array[1], data, X=X_init_og, Xoxides=Xoxides, name_solvus=true, sys_in=sys_in)

# extract water-excess for the first point
if "H2O" in out_mp_init.ph
    local H2O_index = findfirst(==("H2O"), out_mp_init.ph)
    X_mp .= X_init_mp .- ((H2O_comp) .* out_mp_init.ph_frac[H2O_index])
else
    error("Initial composition must contain H2O to proceed with water-excess calculation.")
end

if "H2O" in out_mg_init.ph
    local H2O_index = findfirst(==("H2O"), out_mg_init.ph)
    X_mg = X_init_mg .- (H2O_comp) .* out_mg_init.ph_frac[H2O_index]
else
    error("Initial composition must contain H2O to proceed with water-excess calculation.")
end

if "H2O" in out_og_init.ph
    local H2O_index = findfirst(==("H2O"), out_og_init.ph)
    X_og = X_init_og .- (H2O_comp) .* out_og_init.ph_frac[H2O_index]
else
    error("Initial composition must contain H2O to proceed with water-excess calculation.")
end

X_mp = X_mp ./ sum(X_mp) # normalize to 1
X_mg = X_mg ./ sum(X_mg) # normalize to 1
X_og = X_og ./ sum(X_og) # normalize to 1

# create arrays to store results and initialize parameters
out_mp = Vector{MAGEMin_C.gmin_struct{Float64,Int64}}(undef, length(P_array))
out_mg = Vector{MAGEMin_C.gmin_struct{Float64,Int64}}(undef, length(P_array))
out_og = Vector{MAGEMin_C.gmin_struct{Float64,Int64}}(undef, length(P_array))

# Create variables to store results
# T_mp_solidus = Tmin
# T_mg_solidus = Tmin
# T_og_solidus = Tmin
g_frac_mp = zeros(length(P_array))
g_frac_mg = zeros(length(P_array))
g_frac_og = zeros(length(P_array))
g_frac_mp_total = zeros(length(P_array))
g_frac_mg_total = zeros(length(P_array))
g_frac_og_total = zeros(length(P_array))
g_frac_mp_vol = zeros(length(P_array))
g_frac_mg_vol = zeros(length(P_array))
g_frac_og_vol = zeros(length(P_array))
g_frac_mp_total_vol = zeros(length(P_array))
g_frac_mg_total_vol = zeros(length(P_array))
g_frac_og_total_vol = zeros(length(P_array))
H2O_frac_mp = zeros(length(P_array))
H2O_frac_mg = zeros(length(P_array))
H2O_frac_og = zeros(length(P_array))
H2O_frac_mp_total = zeros(length(P_array))
H2O_frac_mg_total = zeros(length(P_array))
H2O_frac_og_total = zeros(length(P_array))
H2O_frac_mp_vol = zeros(length(P_array))
H2O_frac_mg_vol = zeros(length(P_array))
H2O_frac_og_vol = zeros(length(P_array))
H2O_frac_mp_total_vol = zeros(length(P_array))
H2O_frac_mg_total_vol = zeros(length(P_array))
H2O_frac_og_total_vol = zeros(length(P_array))
Melt_frac_mp_total_vol_cum = zeros(length(P_array))
Melt_frac_mg_total_vol_cum = zeros(length(P_array))
Melt_frac_og_total_vol_cum = zeros(length(P_array))
effective_fractions_mp = zeros(length(P_array))
effective_fractions_mg = zeros(length(P_array))
effective_fractions_og = zeros(length(P_array))
effective_fractions_mp_vol = zeros(length(P_array))
effective_fractions_mg_vol = zeros(length(P_array))
effective_fractions_og_vol = zeros(length(P_array))

# Run the calculations along the path defined by P_array and T_array
@showprogress for step in eachindex(T_array)
    T = T_array[step]
    P = P_array[step]

    # Calculate the effective bulk composition after removing garnet and H2O fractions from the previous step
    if step > 1
        effective_fractions_mp[step] = 1 - g_frac_mp[step-1] - H2O_frac_mp[step-1]
        effective_fractions_mp_vol[step] = 1 - g_frac_mp_vol[step-1] - H2O_frac_mp_vol[step-1]
        effective_fractions_mg[step] = 1 - g_frac_mg[step-1] - H2O_frac_mg[step-1]
        effective_fractions_mg_vol[step] = 1 - g_frac_mg_vol[step-1] - H2O_frac_mg_vol[step-1]
        effective_fractions_og[step] = 1 - g_frac_og[step-1] - H2O_frac_og[step-1]
        effective_fractions_og_vol[step] = 1 - g_frac_og_vol[step-1] - H2O_frac_og_vol[step-1]
    else
        effective_fractions_mp[step] = 1
        effective_fractions_mp_vol[step] = 1
        effective_fractions_mg[step] = 1
        effective_fractions_mg_vol[step] = 1
        effective_fractions_og[step] = 1
        effective_fractions_og_vol[step] = 1
    end

    # Run the minimization for each lithology
    out_og[step] = deepcopy(single_point_minimization(P, T, data, X=X_og, Xoxides=Xoxides, name_solvus=true, sys_in=sys_in))
    out_mp[step] = deepcopy(single_point_minimization(P, T, data, X=X_mp, Xoxides=Xoxides, name_solvus=true, sys_in=sys_in))
    out_mg[step] = deepcopy(single_point_minimization(P, T, data, X=X_mg, Xoxides=Xoxides, name_solvus=true, sys_in=sys_in))

    # extract water-excess
    if "H2O" in out_mp[step].ph
        local H2O_index = findfirst(==("H2O"), out_mp[step].ph)
        H2O_frac_mp[step] = out_mp[step].ph_frac[H2O_index]
        H2O_frac_mp_vol[step] = out_mp[step].ph_frac_vol[H2O_index]
        X_mp .= X_mp .- ((H2O_comp) .* H2O_frac_mp[step])  # remove only the excess
    else
        H2O_frac_mp[step] = 0.0
        H2O_frac_mp_vol[step] = 0.0
    end

    if "H2O" in out_mg[step].ph
        local H2O_index = findfirst(==("H2O"), out_mg[step].ph)
        H2O_frac_mg[step] = out_mg[step].ph_frac[H2O_index]
        H2O_frac_mg_vol[step] = out_mg[step].ph_frac_vol[H2O_index]
        X_mg .= X_mg .- ((H2O_comp) .* H2O_frac_mg[step])  # remove only the excess
    else
        H2O_frac_mg[step] = 0.0
        H2O_frac_mg_vol[step] = 0.0
    end

    if "H2O" in out_og[step].ph
        local H2O_index = findfirst(==("H2O"), out_og[step].ph)
        H2O_frac_og[step] = out_og[step].ph_frac[H2O_index]
        H2O_frac_og_vol[step] = out_og[step].ph_frac_vol[H2O_index]
        X_og .= X_og .- ((H2O_comp) .* out_og[step].ph_frac[H2O_index])
    else
        H2O_frac_og[step] = 0.0
        H2O_frac_og_vol[step] = 0.0
    end

    # garnet fractionation
    if "g" in out_mp[step].ph
        local g_index = findfirst(==("g"), out_mp[step].ph)
        X_mp .= X_mp .- (out_mp[step].SS_vec[g_index].Comp) .* out_mp[step].ph_frac[g_index] .* g_factor
        g_frac_mp[step] = out_mp[step].ph_frac[g_index] .* g_factor
        g_frac_mp_vol[step] = out_mp[step].ph_frac_vol[g_index] .* g_factor
    else
        g_frac_mp[step] = 0.0
        g_frac_mp_vol[step] = 0.0
    end

    if "g" in out_mg[step].ph
        local g_index = findfirst(==("g"), out_mg[step].ph)
        X_mg .= X_mg .- (out_mg[step].SS_vec[g_index].Comp) .* out_mg[step].ph_frac[g_index] .* g_factor
        g_frac_mg[step] = out_mg[step].ph_frac[g_index] .* g_factor
        g_frac_mg_vol[step] = out_mg[step].ph_frac_vol[g_index] .* g_factor
    else
        g_frac_mg[step] = 0.0
        g_frac_mg_vol[step] = 0.0
    end

    if "g" in out_og[step].ph
        local g_index = findfirst(==("g"), out_og[step].ph)
        X_og .= X_og .- (out_og[step].SS_vec[g_index].Comp) .* out_og[step].ph_frac[g_index] .* g_factor
        g_frac_og[step] = out_og[step].ph_frac[g_index] .* g_factor
        g_frac_og_vol[step] = out_og[step].ph_frac_vol[g_index] .* g_factor
    else
        g_frac_og[step] = 0.0
        g_frac_og_vol[step] = 0.0
    end

    # # search for the solidus
    # if !in("liq", out_mp[step].ph)
    #     T_mp_solidus = T
    # end
    # if !in("liq", out_mg[step].ph)
    #     T_mg_solidus = T
    # end
    # if !in("liq", out_og[step].ph)
    #     T_og_solidus = T
    # end
end

T_mp_solidus = findfirst(x -> in("liq", out_mp[x].ph), eachindex(T_array)) + Tmin
T_mg_solidus = findfirst(x -> in("liq", out_mg[x].ph), eachindex(T_array)) + Tmin
T_og_solidus = findfirst(x -> in("liq", out_og[x].ph), eachindex(T_array)) + Tmin

# calculate the effective fractions at each step
effect_fract_mp_total = accumulate(*, effective_fractions_mp)
effect_fract_mg_total = accumulate(*, effective_fractions_mg)
effect_fract_og_total = accumulate(*, effective_fractions_og)
effect_fract_mp_total_vol = accumulate(*, effective_fractions_mp_vol)
effect_fract_mg_total_vol = accumulate(*, effective_fractions_mg_vol)
effect_fract_og_total_vol = accumulate(*, effective_fractions_og_vol)

# Calculate the garnet, H2O and melt fractions in function of effective fractions
for i in eachindex(P_array)
    g_frac_mp_total[i] = g_frac_mp[i] * effect_fract_mp_total[i]
    g_frac_mg_total[i] = g_frac_mg[i] * effect_fract_mg_total[i]
    g_frac_og_total[i] = g_frac_og[i] * effect_fract_og_total[i]
    g_frac_mp_total_vol[i] = g_frac_mp_vol[i] * effect_fract_mp_total_vol[i]
    g_frac_mg_total_vol[i] = g_frac_mg_vol[i] * effect_fract_mg_total_vol[i]
    g_frac_og_total_vol[i] = g_frac_og_vol[i] * effect_fract_og_total_vol[i]
    H2O_frac_mp_total[i] = H2O_frac_mp[i] * effect_fract_mp_total[i]
    H2O_frac_mg_total[i] = H2O_frac_mg[i] * effect_fract_mg_total[i]
    H2O_frac_og_total[i] = H2O_frac_og[i] * effect_fract_og_total[i]
    H2O_frac_mp_total_vol[i] = H2O_frac_mp_vol[i] * effect_fract_mp_total_vol[i]
    H2O_frac_mg_total_vol[i] = H2O_frac_mg_vol[i] * effect_fract_mg_total_vol[i]
    H2O_frac_og_total_vol[i] = H2O_frac_og_vol[i] * effect_fract_og_total_vol[i]
    if "liq" in out_mp[i].ph
        Melt_frac_mp_total_vol_cum[i] = out_mp[i].ph_frac_vol[findfirst(==("liq"), out_mp[i].ph)] * effect_fract_mp_total_vol[i]
    end
    if "liq" in out_mg[i].ph
        Melt_frac_mg_total_vol_cum[i] = out_mg[i].ph_frac_vol[findfirst(==("liq"), out_mg[i].ph)] * effect_fract_mg_total_vol[i]
    end
    if "liq" in out_og[i].ph
        Melt_frac_og_total_vol_cum[i] = out_og[i].ph_frac_vol[findfirst(==("liq"), out_og[i].ph)] * effect_fract_og_total_vol[i]
    end
end

# Cumulate the garnet and H2O fractions
g_frac_mp_total_cum = accumulate(+, g_frac_mp_total)
g_frac_mg_total_cum = accumulate(+, g_frac_mg_total)
g_frac_og_total_cum = accumulate(+, g_frac_og_total)
g_frac_mp_total_vol_cum = accumulate(+, g_frac_mp_total_vol)
g_frac_mg_total_vol_cum = accumulate(+, g_frac_mg_total_vol)
g_frac_og_total_vol_cum = accumulate(+, g_frac_og_total_vol)
H2O_frac_mp_total_cum = accumulate(+, H2O_frac_mp_total)
H2O_frac_mg_total_cum = accumulate(+, H2O_frac_mg_total)
H2O_frac_og_total_cum = accumulate(+, H2O_frac_og_total)
H2O_frac_mp_total_vol_cum = accumulate(+, H2O_frac_mp_total_vol)
H2O_frac_mg_total_vol_cum = accumulate(+, H2O_frac_mg_total_vol)
H2O_frac_og_total_vol_cum = accumulate(+, H2O_frac_og_total_vol)

# Calculate the amount of water released between the two solidi
H2O_mp_melt_og = H2O_frac_mp_total_cum[Int(T_mp_solidus-Tmin)] - H2O_frac_mp_total_cum[Int(T_og_solidus-Tmin)]
H2O_mg_melt_og = H2O_frac_mg_total_cum[Int(T_mg_solidus-Tmin)] - H2O_frac_mg_total_cum[Int(T_og_solidus-Tmin)]

# Print results
println("H2O released mp: $(H2O_mp_melt_og*100)")
println("H2O released mg: $(H2O_mg_melt_og*100)")
println("$T_mp_solidus °C metapelite solidus temperature")
println("$T_mg_solidus °C metagraywacke solidus temperature")
println("$T_og_solidus °C orthogneiss solidus temperature")

# Calculation of the volume of melt created with the water-influx
# Extract excess water at the solidus point
out_og_melTmin = single_point_minimization(P_array[1], T_array[1], data, X=X_init_og, Xoxides=Xoxides, name_solvus= true, sys_in=sys_in);
X_og_melt = copy(X_init_og)
if "H2O" in out_og_melTmin.ph
    local H2O_index = findfirst(==("H2O"), out_og_melTmin.ph)
    X_og_melt = X_init_og .- ((H2O_comp) .* out_og_melTmin.ph_frac[H2O_index])
end
X_og_melt = X_og_melt ./ sum(X_og_melt)
X_og_mp_melt = copy(X_og_melt)
X_og_mg_melt = copy(X_og_melt)

# Create variables to store results
g_frac_og_mp_melt = zeros(length(P_array))
g_frac_og_mp_melt_total = zeros(length(P_array))
g_frac_og_mp_melt_vol = zeros(length(P_array))
g_frac_og_mp_melt_total_vol = zeros(length(P_array))
H2O_frac_og_mp_melt = zeros(length(P_array))
H2O_frac_og_mp_melt_total = zeros(length(P_array))
H2O_frac_og_mp_melt_vol = zeros(length(P_array))
H2O_frac_og_mp_melt_total_vol = zeros(length(P_array))
Melt_frac_og_mp_melt_total_vol_cum = zeros(length(P_array))
effective_fractions_og_mp_melt = zeros(length(P_array))
effective_fractions_og_mp_melt_vol = zeros(length(P_array))

out_og_mp_melt = Vector{MAGEMin_C.gmin_struct{Float64, Int64}}(undef, length(P_array))

g_frac_og_mg_melt = zeros(length(P_array))
g_frac_og_mg_melt_total = zeros(length(P_array))
g_frac_og_mg_melt_vol = zeros(length(P_array))
g_frac_og_mg_melt_total_vol = zeros(length(P_array))
H2O_frac_og_mg_melt = zeros(length(P_array))
H2O_frac_og_mg_melt_total = zeros(length(P_array))
H2O_frac_og_mg_melt_vol = zeros(length(P_array))
H2O_frac_og_mg_melt_total_vol = zeros(length(P_array))
Melt_frac_og_mg_melt_total_vol_cum = zeros(length(P_array))
effective_fractions_og_mg_melt = zeros(length(P_array))
effective_fractions_og_mg_melt_vol = zeros(length(P_array))

out_og_mg_melt = Vector{MAGEMin_C.gmin_struct{Float64, Int64}}(undef, length(P_array))

# Calculation of og_mp melt generation along the path
@showprogress for step in eachindex(T_array)
    T = T_array[step]
    P = P_array[step]
    if "liq" in out_og[step].ph
        external_H2O = 1.0.*(H2O_frac_mp_total[step]) #change coefficient to modify the fraction of metapelite compare to orthogneiss
        external_H2O_vol = 1.0.*(H2O_frac_mp_total_vol[step])
    else
        external_H2O = 0.0
        external_H2O_vol = 0.0
    end
    X_og_mp_melt .= X_og_mp_melt .+ external_H2O .* H2O_comp

    # Calculate the effective bulk composition after removing garnet and H2O fractions from the previous step
    if step > 1
        effective_fractions_og_mp_melt[step] = 1 - g_frac_og_mp_melt[step-1] + external_H2O
        effective_fractions_og_mp_melt_vol[step] = 1 - g_frac_og_mp_melt_vol[step-1] + external_H2O_vol
    else
        effective_fractions_og_mp_melt[step] = 1
        effective_fractions_og_mp_melt_vol[step] = 1
    end

    # Run the minimization
    out_og_mp_melt[step] = deepcopy(single_point_minimization(P, T, data, X=X_og_mp_melt, Xoxides=Xoxides, name_solvus= true, sys_in=sys_in));

    # garnet fractionation
    if "g" in out_og_mp_melt[step].ph
        local g_index = findfirst(==("g"), out_og_mp_melt[step].ph)
        X_og_mp_melt .= X_og_mp_melt .- (out_og_mp_melt[step].SS_vec[g_index].Comp) .* out_og_mp_melt[step].ph_frac[g_index] .* g_factor
        g_frac_og_mp_melt[step] = out_og_mp_melt[step].ph_frac[g_index] .* g_factor
        g_frac_og_mp_melt_vol[step] = out_og_mp_melt[step].ph_frac_vol[g_index] .* g_factor
    else
        g_frac_og_mp_melt[step] = 0.0
        g_frac_og_mp_melt_vol[step] = 0.0
    end

    # extract water-excess
    if "H2O" in out_og_mp_melt[step].ph && step < (T_og_solidus - Tmin)
        local H2O_index = findfirst(==("H2O"), out_og_mp_melt[step].ph)
        H2O_frac_og_mp_melt[step] = out_og_mp_melt[step].ph_frac[H2O_index] 
        H2O_frac_og_mp_melt_vol[step] = out_og_mp_melt[step].ph_frac_vol[H2O_index]
        X_og_mp_melt .= X_og_mp_melt .- ((H2O_comp) .* out_og_mp_melt[step].ph_frac[H2O_index])
    else
        H2O_frac_og_mp_melt[step] = 0.0
        H2O_frac_og_mp_melt_vol[step] = 0.0
    end
end

# calculate the effective fractions at each step
effect_fract_og_mp_melt_total = accumulate(*, effective_fractions_og_mp_melt)
effect_fract_og_mp_melt_total_vol = accumulate(*, effective_fractions_og_mp_melt_vol)

# Calculate the garnet, H2O and melt fractions in function of effective fractions
for i in eachindex(P_array)
    g_frac_og_mp_melt_total[i] = g_frac_og_mp_melt[i] * effect_fract_og_mp_melt_total[i]
    g_frac_og_mp_melt_total_vol[i] = g_frac_og_mp_melt_vol[i] * effect_fract_og_mp_melt_total_vol[i]
    H2O_frac_og_mp_melt_total[i] = H2O_frac_og_mp_melt[i] * effect_fract_og_mp_melt_total[i]
    H2O_frac_og_mp_melt_total_vol[i] = H2O_frac_og_mp_melt_vol[i] * effect_fract_og_mp_melt_total_vol[i]
    if "liq" in out_og_mp_melt[i].ph
        Melt_frac_og_mp_melt_total_vol_cum[i] = out_og_mp_melt[i].ph_frac_vol[findfirst(==("liq"), out_og_mp_melt[i].ph)] * effect_fract_og_mp_melt_total_vol[i]
    end
end

# Cumulate the garnet and H2O fractions
g_frac_og_mp_melt_total_cum = accumulate(+, g_frac_og_mp_melt_total)
g_frac_og_mp_melt_total_vol_cum = accumulate(+, g_frac_og_mp_melt_total_vol)
H2O_frac_og_mp_melt_total_cum = accumulate(+, H2O_frac_og_mp_melt_total)
H2O_frac_og_mp_melt_total_vol_cum = accumulate(+, H2O_frac_og_mp_melt_total_vol)


# Calculation of og_mg melt generation along the path
@showprogress for step in eachindex(T_array)
    T = T_array[step]
    P = P_array[step]
    if "liq" in out_og[step].ph
        external_H2O = 3.0.*(H2O_frac_mg_total[step])
        external_H2O_vol = 3.0.*(H2O_frac_mg_total_vol[step])
    else
        external_H2O = 0.0
        external_H2O_vol = 0.0
    end
    X_og_mg_melt .= X_og_mg_melt .+ external_H2O .* H2O_comp

    if step > 1
        effective_fractions_og_mg_melt[step] = 1 - g_frac_og_mg_melt[step-1] + external_H2O
        effective_fractions_og_mg_melt_vol[step] = 1 - g_frac_og_mg_melt_vol[step-1] + external_H2O_vol
    else
        effective_fractions_og_mg_melt[step] = 1
        effective_fractions_og_mg_melt_vol[step] = 1
    end

    out_og_mg_melt[step] = deepcopy(single_point_minimization(P, T, data, X=X_og_mg_melt, Xoxides=Xoxides, name_solvus= true, sys_in=sys_in));

    # garnet fractionation
    if "g" in out_og_mg_melt[step].ph
        local g_index = findfirst(==("g"), out_og_mg_melt[step].ph)
        X_og_mg_melt .= X_og_mg_melt .- (out_og_mg_melt[step].SS_vec[g_index].Comp) .* out_og_mg_melt[step].ph_frac[g_index] .* g_factor
        g_frac_og_mg_melt[step] = out_og_mg_melt[step].ph_frac[g_index] .* g_factor
        g_frac_og_mg_melt_vol[step] = out_og_mg_melt[step].ph_frac_vol[g_index] .* g_factor
    else
        g_frac_og_mg_melt[step] = 0.0
        g_frac_og_mg_melt_vol[step] = 0.0
    end

    # extract water-excess
    if "H2O" in out_og_mg_melt[step].ph && step < (T_og_solidus - Tmin)
        local H2O_index = findfirst(==("H2O"), out_og_mg_melt[step].ph)
        H2O_frac_og_mg_melt[step] = out_og_mg_melt[step].ph_frac[H2O_index]
        H2O_frac_og_mg_melt_vol[step] = out_og_mg_melt[step].ph_frac_vol[H2O_index]
        X_og_mg_melt .= X_og_mg_melt .- ((H2O_comp) .* out_og_mg_melt[step].ph_frac[H2O_index])
    else
        H2O_frac_og_mg_melt[step] = 0.0
        H2O_frac_og_mg_melt_vol[step] = 0.0
    end
end

# calculate the effective fractions at each step
effect_fract_og_mg_melt_total = accumulate(*, effective_fractions_og_mg_melt)
effect_fract_og_mg_melt_total_vol = accumulate(*, effective_fractions_og_mg_melt_vol)

# Calculate the garnet, H2O and melt fractions in function of effective fractions
for i in eachindex(P_array)
    g_frac_og_mg_melt_total[i] = g_frac_og_mg_melt[i] * effect_fract_og_mg_melt_total[i]
    g_frac_og_mg_melt_total_vol[i] = g_frac_og_mg_melt_vol[i] * effect_fract_og_mg_melt_total_vol[i]
    H2O_frac_og_mg_melt_total[i] = H2O_frac_og_mg_melt[i] * effect_fract_og_mg_melt_total[i]
    H2O_frac_og_mg_melt_total_vol[i] = H2O_frac_og_mg_melt_vol[i] * effect_fract_og_mg_melt_total_vol[i]
    if "liq" in out_og_mg_melt[i].ph
        Melt_frac_og_mg_melt_total_vol_cum[i] = out_og_mg_melt[i].ph_frac_vol[findfirst(==("liq"), out_og_mg_melt[i].ph)] * effect_fract_og_mg_melt_total_vol[i]
    end
end

# Cumulate the garnet and H2O fractions
g_frac_og_mg_melt_total_cum = accumulate(+, g_frac_og_mg_melt_total)
g_frac_og_mg_melt_total_vol_cum = accumulate(+, g_frac_og_mg_melt_total_vol)
H2O_frac_og_mg_melt_total_cum = accumulate(+, H2O_frac_og_mg_melt_total)
H2O_frac_og_mg_melt_total_vol_cum = accumulate(+, H2O_frac_og_mg_melt_total_vol)

# search for the temperature where metasedimentary rocks become more melted than orthogneiss+water
idx_mp_exceeds = findfirst(Melt_frac_mp_total_vol_cum .> Melt_frac_og_mp_melt_total_vol_cum)
idx_mg_exceeds = findfirst(Melt_frac_mg_total_vol_cum .> Melt_frac_og_mg_melt_total_vol_cum)

# If MP never exceeds OG+MP, set to last index (750°C)
if isnothing(idx_mp_exceeds)
    idx_mp_exceeds = length(T_array)
end
if isnothing(idx_mg_exceeds)
    idx_mg_exceeds = length(T_array)
end

T_og_mp = T_array[idx_mp_exceeds] - T_og_solidus
T_og_mg = T_array[idx_mg_exceeds] - T_og_solidus

# Print results
println("T og+mp melt = $(T_array[idx_mp_exceeds]) °C")
println("T og+mg melt = $(T_array[idx_mg_exceeds]) °C")
println("melt_fraction_og_mp = $(Melt_frac_og_mp_melt_total_vol_cum[idx_mp_exceeds]*100)")
println("melt_fraction_og_mg = $(Melt_frac_og_mg_melt_total_vol_cum[idx_mg_exceeds]*100)")

println("melt composition og+mp at T_mp = $(T_array[idx_mp_exceeds]) °C = $(out_og_mp_melt[idx_mp_exceeds].SS_vec[findfirst(==("liq"), out_og_mp_melt[idx_mp_exceeds].ph)].Comp)")
println("melt composition og+mg at T_mg = $(T_array[idx_mg_exceeds]) °C = $(out_og_mg_melt[idx_mg_exceeds].SS_vec[findfirst(==("liq"), out_og_mg_melt[idx_mg_exceeds].ph)].Comp)")
println("melt composition mp at T = $(T_array[idx_mp_exceeds]) °C = $(out_mp[idx_mp_exceeds].SS_vec[findfirst(==("liq"), out_mp[idx_mp_exceeds].ph)].Comp)")
println("melt composition mg at T = $(T_array[idx_mg_exceeds]) °C = $(out_mg[idx_mg_exceeds].SS_vec[findfirst(==("liq"), out_mg[idx_mg_exceeds].ph)].Comp)")

println("melt composition og+mp at T = $(Int(T_mp_solidus-Tmin)) °C = $(out_og_mp_melt[Int(T_mp_solidus-Tmin)].SS_vec[findfirst(==("liq"), out_og_mp_melt[Int(T_mp_solidus-Tmin)].ph)].Comp)")
println("melt composition og+mg at T = $(Int(T_mg_solidus-Tmin)) °C = $(out_og_mg_melt[Int(T_mg_solidus-Tmin)].SS_vec[findfirst(==("liq"), out_og_mg_melt[Int(T_mg_solidus-Tmin)].ph)].Comp)")
println("melt composition og at T = $(Int(T_mp_solidus-Tmin)) °C = $(out_og[Int(T_mp_solidus-Tmin)].SS_vec[findfirst(==("liq"), out_og[Int(T_mp_solidus-Tmin)].ph)].Comp)")
println("melt composition og at T = $(Int(T_mg_solidus-Tmin)) °C = $(out_og[Int(T_mg_solidus-Tmin)].SS_vec[findfirst(==("liq"), out_og[Int(T_mg_solidus-Tmin)].ph)].Comp)")

# Plot the results

# Plot the water released
fig1 = Figure(size = (800, 700))
ax = Axis(fig1[1, 1], xlabel="T (°C)", ylabel="H2O_frac", title="H2O_frac vs T")
xlims!(ax, minimum(T_array), maximum(T_array))
ylims!(ax, 0.0, 0.05)
lines!(ax, T_array, H2O_frac_mp_total_cum, color=:black, label="metapelite")
lines!(ax, T_array, H2O_frac_mg_total_cum, color=:blue, label="metagraywacke")
# Add vertical dashed lines for solidus temperatures
vlines!(ax, T_mp_solidus, linestyle=:dash, color=:black, linewidth=1.5, label="T solidus MP")
vlines!(ax, T_mg_solidus, linestyle=:dash, color=:blue, linewidth=1.5, label="T solidus MG")
vlines!(ax, T_og_solidus, linestyle=:dash, color=:red, linewidth=1.5, label="T solidus OG")
# Add text annotations with the water values
text!(ax, T_og_solidus + 5, 0.045, text="H₂O MP = $(round(H2O_mp_melt_og, digits=4))\nH₂O MG = $(round(H2O_mg_melt_og, digits=4))\nH₂O Total = $(round(H2O_total_melt_og, digits=4))", 
        fontsize=12, align=(:left, :top))
axislegend(ax, position=:lt)
display(fig1)  # Display the plot


# Plot the melt fractions
fig2 = Figure(size = (800, 700))
ax = Axis(fig2[1, 1], xlabel="T (°C)", ylabel="Liq_frac", title="Liq_frac vs T")
xlims!(ax, minimum(T_array), maximum(T_array))
ylims!(ax, 0.0, 0.1)
lines!(ax, T_array, Melt_frac_og_mp_melt_total_vol_cum, color=:red, label="orthogneiss+water-mp")
lines!(ax, T_array, Melt_frac_og_mg_melt_total_vol_cum, color=:purple, label="orthogneiss+water-mg")
lines!(ax, T_array, Melt_frac_mp_total_vol_cum, color=:black, label="metapelite")
lines!(ax, T_array, Melt_frac_mg_total_vol_cum, color=:blue, label="metagraywacke")
lines!(ax, T_array, Melt_frac_og_total_vol_cum, color=:pink, label="orthogneiss")
# Add vertical dashed lines for solidus temperatures
vlines!(ax, T_mp_solidus, linestyle=:dash, color=:black, linewidth=1.5, label="T solidus MP")
vlines!(ax, T_mg_solidus, linestyle=:dash, color=:blue, linewidth=1.5, label="T solidus MG")
vlines!(ax, T_og_solidus, linestyle=:dash, color=:red, linewidth=1.5, label="T solidus OG")
axislegend(ax, position=:lt)
display(fig2)  # Display the plot


# Create mode box diagrams for each rock type
# Extract phase fractions for metapelite (excluding fractionated garnet)
n_steps = length(T_array)

# Collect all unique phases across all temperature steps (excluding garnet)
phase_names_mp = unique(vcat([out_mp[i].ph for i in 1:n_steps]...))
phase_names_mp = filter(x -> x != "g", phase_names_mp)  # Exclude garnet
n_phases_mp = length(phase_names_mp)
phase_fractions_mp = zeros(n_steps, n_phases_mp)

for i in 1:n_steps
    for (j, phase) in enumerate(out_mp[i].ph)
        if phase != "g"  # Skip garnet as it's fractionated
            phase_idx = findfirst(==(phase), phase_names_mp)
            if !isnothing(phase_idx)
                # Scale by effective fraction to show what remains after fractionation
                phase_fractions_mp[i, phase_idx] = out_mp[i].ph_frac[j] * effect_fract_mp_total[i]
            end
        end
    end
end

# Extract phase fractions for metagraywacke (excluding fractionated garnet)
# Collect all unique phases across all temperature steps (excluding garnet)
phase_names_mg = unique(vcat([out_mg[i].ph for i in 1:n_steps]...))
phase_names_mg = filter(x -> x != "g", phase_names_mg)  # Exclude garnet
n_phases_mg = length(phase_names_mg)
phase_fractions_mg = zeros(n_steps, n_phases_mg)

for i in 1:n_steps
    for (j, phase) in enumerate(out_mg[i].ph)
        if phase != "g"  # Skip garnet as it's fractionated
            phase_idx = findfirst(==(phase), phase_names_mg)
            if !isnothing(phase_idx)
                # Scale by effective fraction to show what remains after fractionation
                phase_fractions_mg[i, phase_idx] = out_mg[i].ph_frac[j] * effect_fract_mg_total[i]
            end
        end
    end
end

# Extract phase fractions for orthogneiss_water-mp (excluding fractionated garnet)
# Collect all unique phases across all temperature steps (excluding garnet)
phase_names_og_mp = unique(vcat([out_og_mp_melt[i].ph for i in 1:n_steps]...))
phase_names_og_mp = filter(x -> x != "g", phase_names_og_mp)  # Exclude garnet
n_phases_og_mp = length(phase_names_og_mp)
phase_fractions_og_mp = zeros(n_steps, n_phases_og_mp)

for i in 1:n_steps
    for (j, phase) in enumerate(out_og_mp_melt[i].ph)
        if phase != "g"  # Skip garnet as it's fractionated
            phase_idx = findfirst(==(phase), phase_names_og_mp)
            if !isnothing(phase_idx)
                # Scale by effective fraction to show what remains after fractionation
                phase_fractions_og_mp[i, phase_idx] = out_og_mp_melt[i].ph_frac[j] * effect_fract_og_mp_melt_total[i]
            end
        end
    end
end

# Extract phase fractions for orthogneiss_water-mg (excluding fractionated garnet)
# Collect all unique phases across all temperature steps (excluding garnet)
phase_names_og_mg = unique(vcat([out_og_mg_melt[i].ph for i in 1:n_steps]...))
phase_names_og_mg = filter(x -> x != "g", phase_names_og_mg)  # Exclude garnet
n_phases_og_mg = length(phase_names_og_mg)
phase_fractions_og_mg = zeros(n_steps, n_phases_og_mg)

for i in 1:n_steps
    for (j, phase) in enumerate(out_og_mg_melt[i].ph)
        if phase != "g"  # Skip garnet as it's fractionated
            phase_idx = findfirst(==(phase), phase_names_og_mg)
            if !isnothing(phase_idx)
                # Scale by effective fraction to show what remains after fractionation
                phase_fractions_og_mg[i, phase_idx] = out_og_mg_melt[i].ph_frac[j] * effect_fract_og_mg_melt_total[i]
            end
        end
    end
end

# Create a unified color map for all phases across all three lithologies
all_phases = unique(vcat(phase_names_mp, phase_names_mg, phase_names_og_mp, phase_names_og_mg))

# Manual color assignment for each phase, from Lanari and Tedeschi (2024)
# You can customize these colors for each specific phase name
phase_color_map = Dict(
    "liq"  => :red,                   # Melt/liquid
    "pl"   => (0.925, 0.863, 0.620),  # Plagioclase
    "q"    => (1.000, 1.000, 1.000),  # Quartz
    "bi"   => (0.463, 0.224, 0.129),  # Biotite
    "mu"   => (0.851, 0.780, 0.796),  # Muscovite
    "afs"  => (0.914, 0.722, 0.800),  # K-feldspar
    "ilm"  => (0.878, 0.365, 0.165),  # Ilmenite
    "sill" => (0.000, 0.561, 0.737),  # Sillimanite
    "and"  => (0.498, 0.702, 0.757),  # Andalusite
    "ky"   => (0.000, 0.357, 0.576),  # Kyanite
    "crd"  => (0.533, 0.408, 0.659),  # Cordierite
    "zo"   => (0.733, 0.706, 0.216),  # Zoisite
    "H2O"  => :lightblue,             # Water
    "g"    => (0.710, 0.173, 0.122),  # Garnet
    "st"   => (0.886, 0.667, 0.000),  # Staurolite
    "ttn"  => (0.745, 0.541, 0.263),  # Titanite
    "ru"   => (0.024, 0.259, 0.463),  # Rutile
    "chl"  => (0.506, 0.753, 0.443),  # Chlorite
    "cld"  => (0.369, 0.467, 0.380),  # Chloritoid 
    "opx"  => (0.882, 0.486, 0.565),  # Orthopyroxene
    "pat"  => (0.855, 0.514, 0.314)   # Paragonite
)

# For any phases not manually defined, assign default colors
base_colors = [:lightcoral, :lightskyblue, :lightgreen, :lightsalmon, :plum, :lightcyan, :khaki, :lightpink, :tan, :lightgray, :wheat, :peachpuff, :lavender, :thistle, :powderblue]
for (i, phase) in enumerate(all_phases)
    if !haskey(phase_color_map, phase)
        phase_color_map[phase] = base_colors[mod1(i, length(base_colors))]
    end
end

# Create fig3 with 4 mode boxes (2x2 layout)
fig3 = Figure(size = (1800, 1400))

# Panel 1: Metapelite
ax3_mp = Axis(fig3[1, 1], xlabel="T (°C)", ylabel="Phase Fraction (molar)", title="Mode Box - Metapelite")
xlims!(ax3_mp, minimum(T_array), maximum(T_array))
ylims!(ax3_mp, 0.0, 1.0)

cumulative_mp = zeros(n_steps)
for j in 1:n_phases_mp
    lower = copy(cumulative_mp)
    cumulative_mp .+= phase_fractions_mp[:, j]
    band!(ax3_mp, T_array, lower, cumulative_mp, color=phase_color_map[phase_names_mp[j]], label=phase_names_mp[j])
end

# Add garnet on top to show fractionated amount
lower_g_mp = copy(cumulative_mp)
cumulative_g_mp = cumulative_mp .+ g_frac_mp_total_cum
band!(ax3_mp, T_array, lower_g_mp, cumulative_g_mp, color=:indianred, label="g (fractionated)")

# Add H2O on top to show fractionated water
lower_h2o_mp = copy(cumulative_g_mp)
cumulative_h2o_mp = cumulative_g_mp .+ H2O_frac_mp_total_cum
band!(ax3_mp, T_array, lower_h2o_mp, cumulative_h2o_mp, color=:lightblue, label="H2O (fractionated)")

# Add vertical dashed lines for solidus temperatures
vlines!(ax3_mp, T_mp_solidus, linestyle=:dash, color=:black, linewidth=1.5)
vlines!(ax3_mp, T_og_solidus, linestyle=:dash, color=:red, linewidth=1.5)

# Panel 2: Metagraywacke
ax3_mg = Axis(fig3[1, 2], xlabel="T (°C)", ylabel="Phase Fraction (molar)", title="Mode Box - Metagraywacke")
xlims!(ax3_mg, minimum(T_array), maximum(T_array))
ylims!(ax3_mg, 0.0, 1.0)

cumulative_mg = zeros(n_steps)
for j in 1:n_phases_mg
    lower = copy(cumulative_mg)
    cumulative_mg .+= phase_fractions_mg[:, j]
    band!(ax3_mg, T_array, lower, cumulative_mg, color=phase_color_map[phase_names_mg[j]], label=phase_names_mg[j])
end

# Add garnet on top to show fractionated amount
lower_g_mg = copy(cumulative_mg)
cumulative_g_mg = cumulative_mg .+ g_frac_mg_total_cum
band!(ax3_mg, T_array, lower_g_mg, cumulative_g_mg, color=:indianred, label="g (fractionated)")

# Add H2O on top to show fractionated water
lower_h2o_mg = copy(cumulative_g_mg)
cumulative_h2o_mg = cumulative_g_mg .+ H2O_frac_mg_total_cum
band!(ax3_mg, T_array, lower_h2o_mg, cumulative_h2o_mg, color=:lightblue, label="H2O (fractionated)")

# Add vertical dashed lines for solidus temperatures
vlines!(ax3_mg, T_mg_solidus, linestyle=:dash, color=:blue, linewidth=1.5)
vlines!(ax3_mg, T_og_solidus, linestyle=:dash, color=:red, linewidth=1.5)

# Panel 3: Orthogneiss + Water from Metapelite
ax3_og_mp = Axis(fig3[2, 1], xlabel="T (°C)", ylabel="Phase Fraction (molar)", title="Mode Box - Orthogneiss + H₂O from MP")
xlims!(ax3_og_mp, minimum(T_array), maximum(T_array))
ylims!(ax3_og_mp, 0.0, 1.0)

cumulative_og_mp = zeros(n_steps)
for j in 1:n_phases_og_mp
    lower = copy(cumulative_og_mp)
    cumulative_og_mp .+= phase_fractions_og_mp[:, j]
    band!(ax3_og_mp, T_array, lower, cumulative_og_mp, color=phase_color_map[phase_names_og_mp[j]], label=phase_names_og_mp[j])
end

# Add garnet on top to show fractionated amount
lower_g_og_mp = copy(cumulative_og_mp)
cumulative_g_og_mp = cumulative_og_mp .+ g_frac_og_mp_melt_total_cum
band!(ax3_og_mp, T_array, lower_g_og_mp, cumulative_g_og_mp, color=:indianred, label="g (fractionated)")

# Add H2O on top to show fractionated water
lower_h2o_og_mp = copy(cumulative_g_og_mp)
cumulative_h2o_og_mp = cumulative_g_og_mp .+ H2O_frac_og_mp_melt_total_cum
band!(ax3_og_mp, T_array, lower_h2o_og_mp, cumulative_h2o_og_mp, color=:lightblue, label="H2O (fractionated)")

# Add vertical dashed lines for solidus temperatures
vlines!(ax3_og_mp, T_og_solidus, linestyle=:dash, color=:red, linewidth=1.5)
vlines!(ax3_og_mp, T_mp_solidus, linestyle=:dash, color=:black, linewidth=1.5)

# Panel 4: Orthogneiss + Water from Metagraywacke
ax3_og_mg = Axis(fig3[2, 2], xlabel="T (°C)", ylabel="Phase Fraction (molar)", title="Mode Box - Orthogneiss + H₂O from MG")
xlims!(ax3_og_mg, minimum(T_array), maximum(T_array))
ylims!(ax3_og_mg, 0.0, 1.0)

cumulative_og_mg = zeros(n_steps)
for j in 1:n_phases_og_mg
    lower = copy(cumulative_og_mg)
    cumulative_og_mg .+= phase_fractions_og_mg[:, j]
    band!(ax3_og_mg, T_array, lower, cumulative_og_mg, color=phase_color_map[phase_names_og_mg[j]], label=phase_names_og_mg[j])
end

# Add garnet on top to show fractionated amount
lower_g_og_mg = copy(cumulative_og_mg)
cumulative_g_og_mg = cumulative_og_mg .+ g_frac_og_mg_melt_total_cum
band!(ax3_og_mg, T_array, lower_g_og_mg, cumulative_g_og_mg, color=:indianred, label="g (fractionated)")

# Add H2O on top to show fractionated water
lower_h2o_og_mg = copy(cumulative_g_og_mg)
cumulative_h2o_og_mg = cumulative_g_og_mg .+ H2O_frac_og_mg_melt_total_cum
band!(ax3_og_mg, T_array, lower_h2o_og_mg, cumulative_h2o_og_mg, color=:lightblue, label="H2O (fractionated)")

# Add vertical dashed lines for solidus temperatures
vlines!(ax3_og_mg, T_og_solidus, linestyle=:dash, color=:red, linewidth=1.5)
vlines!(ax3_og_mg, T_mg_solidus, linestyle=:dash, color=:blue, linewidth=1.5)

# Create a common legend for all unique phases across all four subplots
legend_elements = []
legend_labels = []

# Add all unique phases with their colors
for phase in sort(all_phases)
    push!(legend_elements, PolyElement(color=phase_color_map[phase]))
    push!(legend_labels, phase)
end

# Add garnet and H2O
push!(legend_elements, PolyElement(color=:indianred))
push!(legend_labels, "g (fractionated)")
push!(legend_elements, PolyElement(color=:lightblue))
push!(legend_labels, "H2O (fractionated)")

# Add solidus lines
push!(legend_elements, LineElement(color=:black, linestyle=:dash))
push!(legend_labels, "T solidus MP")
push!(legend_elements, LineElement(color=:blue, linestyle=:dash))
push!(legend_labels, "T solidus MG")
push!(legend_elements, LineElement(color=:red, linestyle=:dash))
push!(legend_labels, "T solidus OG")

# Place legend spanning both rows on the right side
Legend(fig3[1:2, 3], legend_elements, legend_labels, framevisible=true)

display(fig3)

save("H2O_frac_PT$(i).svg", fig1)
save("Liq_frac_PT$(i).svg",  fig2)
save("Mode_Boxes_4Panel_PT$(i).svg", fig3)

Finalize_MAGEMin(data)