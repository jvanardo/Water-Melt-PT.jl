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
- Vielzeuf & Montel (1994)  - Metagraywacke bulk compositions
- Weisbrod, A. (1970) - Orthogneiss bulk compositions
- Lanari, P. and Tedeschi, M. (2024) - Phase color conventions

"""
module WaterMeltPT

using Reexport
@reexport using MAGEMin_C
using CSV
using DataFrames
using ProgressMeter
@reexport using CairoMakie

export prepare_bulk_composition, define_PT_path, extract_init_H2O
export run_PT_path, calculate_additional_H2O, calculate_effect_frac, calculate_melt_frac
export extract_H2O_PT, fractionate_Grt
export plot_fig_pl, plot_fig_melt_comp, plot_fig1, plot_fig2, plot_fig3, collect_phases, ax_boxplot

function prepare_bulk_composition(X_init_wt, X_init_ox, P_array, T_array, data, sys_in;
                                M_oxides = [60.08; 101.96; 56.08; 40.30; 71.85; 159.69; 94.20; 61.98; 79.87; 70.94; 18.02],
                                H2O_comp = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 1.0])
    X_init, ox = convertBulk4MAGEMin(X_init_wt, X_init_ox, "wt", "mp")
    X_init = X_init ./ sum(X_init) # normalize to 1
    out_init = single_point_minimization(P_array[1], T_array[1], data, X=X_init, Xoxides=ox, name_solvus=true, sys_in=sys_in)
    X = extract_init_H2O(out_init, X_init)
    X_wt = X .* M_oxides
    X_wt = 100 * X_wt ./ sum(X_wt)
    return X, X_wt
end

function define_PT_path(Tmin, Tmax, Pmin, Pmax)
    T_array = collect(Tmin:1.0:Tmax)
    P_resolution = (Pmax - Pmin) / (length(T_array) - 1)
    if P_resolution == 0.0
        P_array = Pmin * (ones(length(T_array)))
    else
        P_array = collect(Pmin:P_resolution:Pmax)
    end
    return T_array, P_array
end

function extract_init_H2O(out, X; H2O_comp = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 1.0])
    if "H2O" in out.ph
        local H2O_index = findfirst(==("H2O"), out.ph)
        X .= X .- ((H2O_comp) .* out.ph_frac[H2O_index])
        X = X ./ sum(X) # normalize to 1 after removing the excess water
        return X  # remove only the excess
    else
        error("Initial composition  of $out must contain H2O to proceed with water-excess calculation.")
    end
end

function run_PT_path(T_array, P_array, data, X, Xoxides, g_factor, sys_in, Tmin;
    additional_H2O = zeros(length(T_array), 3),
    H2O_comp = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 1.0])

    out = Vector{MAGEMin_C.gmin_struct{Float64,Int64}}(undef, length(P_array))
    # mol%, vol%, wt%
    g_frac = zeros(length(P_array), 3)
    g_frac_total = zeros(length(P_array), 3)
    g_frac_total_cum = zeros(length(P_array), 3)
    H2O_frac = zeros(length(P_array), 3)
    H2O_frac_total = zeros(length(P_array), 3)
    H2O_frac_total_cum = zeros(length(P_array), 3)
    melt_frac = zeros(length(P_array), 3)
    effect_frac = zeros(length(P_array), 3)
    effect_frac_total = zeros(length(P_array), 3)

    @showprogress for step in eachindex(T_array)
        T = T_array[step]
        P = P_array[step]

        if additional_H2O[step,1] !== 0.0
            X .= X .+ (additional_H2O[step, 1] .* H2O_comp)
        end

        # Calculate the effective bulk composition after removing garnet and H2O frac from the previous step
        if step > 1
            effect_frac[step,:] .= 1 .- g_frac[step-1,:] .- H2O_frac[step-1,:] .+ additional_H2O[step,:]
        else
            effect_frac[step,:] .= 1.0
        end

        # Run the minimization for each lithology
        out[step] = deepcopy(single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, name_solvus=true, sys_in=sys_in))

        # extract water-excess
        X, H2O_frac[step,:] = extract_H2O_PT(out[step], X)

        # garnet fractionation
        X, g_frac[step,:] = fractionate_Grt(out[step], X, g_factor)
    end

    effect_frac_total .= accumulate(.*, effect_frac, dims=1)

    # Calculate the garnet, H2O and melt frac in function of effective frac
    for i in eachindex(P_array)
        g_frac_total[i,:], H2O_frac_total[i,:] = calculate_effect_frac(g_frac[i,:], H2O_frac[i,:], effect_frac_total[i,:])
        melt_frac[i,:] = calculate_melt_frac(out[i], effect_frac_total[i,:])
    end

    # Cumulate the garnet and H2O frac
    g_frac_total_cum .= accumulate(+, g_frac_total, dims=1)
    H2O_frac_total_cum .= accumulate(+, H2O_frac_total, dims=1)
    T_solidus = findfirst(x -> in("liq", out[x].ph), eachindex(T_array)) + (Tmin - 2)
    return out, g_frac_total_cum, H2O_frac_total, H2O_frac_total_cum, melt_frac, T_solidus, effect_frac_total
end

function calculate_additional_H2O(T_array, out_og, og_fraction, H2O_frac_total)
    metased_fraction = 1 - og_fraction
    ratio = metased_fraction / og_fraction
    additional_H2O = zeros(length(T_array), 3)
    for step in eachindex(T_array)
        if "liq" in out_og[step].ph
            additional_H2O[step,:] = ratio.*(H2O_frac_total[step,:])
        else
            additional_H2O[step,:] = zeros(3)
        end
    end
    return additional_H2O
end

function calculate_effect_frac(g_frac, H2O_frac, effect_frac)
    g_frac_tot = zeros(3)
    g_frac_tot .= g_frac .* effect_frac
    H2O_frac_tot = zeros(3)
    H2O_frac_tot .= H2O_frac .* effect_frac
    return g_frac_tot, H2O_frac_tot
end

function calculate_melt_frac(out, effect_frac)
    melt_fract = zeros(3)
     if "liq" in out.ph
        melt_fract[1] = out.ph_frac[findfirst(==("liq"), out.ph)] * effect_frac[1]
        melt_fract[2] = out.ph_frac_vol[findfirst(==("liq"), out.ph)] * effect_frac[2]
        melt_fract[3] = out.ph_frac_wt[findfirst(==("liq"), out.ph)] * effect_frac[3]
    end
    return melt_fract
end

function extract_H2O_PT(out, X; H2O_comp = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 1.0])
    H2O_frac = zeros(3)
    if "H2O" in out.ph
        H2O_index = findfirst(==("H2O"), out.ph)
        H2O_frac[1] = out.ph_frac[H2O_index]
        H2O_frac[2] = out.ph_frac_vol[H2O_index]
        H2O_frac[3] = out.ph_frac_wt[H2O_index]
        X .= X .- ((H2O_comp) .* H2O_frac[1])  # remove only the excess
    end
    return X, H2O_frac
end

function fractionate_Grt(out, X, g_factor)
    g_frac = zeros(3)
    if "g" in out.ph
        g_index = findfirst(==("g"), out.ph)
        X .= X .- (out.SS_vec[g_index].Comp) .* out.ph_frac[g_index] .* g_factor
        g_frac[1] = out.ph_frac[g_index] .* g_factor
        g_frac[2] = out.ph_frac_vol[g_index] .* g_factor
        g_frac[3] = out.ph_frac_wt[g_index] .* g_factor
    end
    return X, g_frac
end

function plot_fig_pl(T_array, out_mp, out_mg, out_og, T_solidus_mp, T_solidus_mg, T_solidus_og)
    Pl_mp = zeros(length(T_array))
    Pl_mg = zeros(length(T_array))
    Pl_og = zeros(length(T_array))
    for step in eachindex(T_array)
        if "pl" in out_mp[step].ph
            idx_an_mp = findfirst(==("an"), out_mp[step].SS_vec[findfirst(==("pl"), out_mp[step].ph)].emNames)
            Pl_mp[step] = out_mp[step].SS_vec[findfirst(==("pl"), out_mp[step].ph)].emFrac[idx_an_mp]
        else
            Pl_mp[step] = NaN
        end
        if "pl" in out_mg[step].ph
            idx_an_mg = findfirst(==("an"), out_mg[step].SS_vec[findfirst(==("pl"), out_mg[step].ph)].emNames)
            Pl_mg[step] = out_mg[step].SS_vec[findfirst(==("pl"), out_mg[step].ph)].emFrac[idx_an_mg]
        else
            Pl_mg[step] = NaN
        end
        if "pl" in out_og[step].ph
            idx_an_og = findfirst(==("an"), out_og[step].SS_vec[findfirst(==("pl"), out_og[step].ph)].emNames)
            Pl_og[step] = out_og[step].SS_vec[findfirst(==("pl"), out_og[step].ph)].emFrac[idx_an_og]
        else
            Pl_og[step] = NaN
        end
    end

    fig_pl = Figure(size = (800, 600), fontsize = 30)
    ax_pl = Axis(fig_pl[1, 1], xlabel="T (°C)", ylabel="Xan")
    xlims!(ax_pl, minimum(T_array), maximum(T_array))
    lines!(ax_pl, T_array, Pl_mp, color=:black, label="metapelite")
    lines!(ax_pl, T_array, Pl_mg, color=:gray, label="metagraywacke")
    lines!(ax_pl, T_array, Pl_og, color=:red, label="orthogneiss")
    vlines!(ax_pl, T_solidus_mp, linestyle=:dash, color=:black, linewidth=1.5, label="T solidus MP")
    vlines!(ax_pl, T_solidus_mg, linestyle=:dash, color=:gray, linewidth=1.5, label="T solidus MG")
    vlines!(ax_pl, T_solidus_og, linestyle=:dash, color=:red, linewidth=1.5, label="T solidus OG")
    display(fig_pl)  # Display the plot
    return fig_pl
end

function plot_fig_melt_comp(T_array, out_mp, out_mg, out_og, out_og_mp_melt, out_og_mg_melt)
    mp_melt_comp = Matrix{Float64}(undef, 11, length(T_array))
    mg_melt_comp = Matrix{Float64}(undef, 11, length(T_array))
    og_melt_comp = Matrix{Float64}(undef, 11, length(T_array))
    og_mp_melt_comp = Matrix{Float64}(undef, 11, length(T_array))
    og_mg_melt_comp = Matrix{Float64}(undef, 11, length(T_array))
    for step in eachindex(T_array)
        if "liq" in out_mp[step].ph
            mp_melt_comp[:,step] = out_mp[step].SS_vec[findfirst(==("liq"), out_mp[step].ph)].Comp_wt .* 100
        else
            mp_melt_comp[:,step] = [NaN for _ in 1:11]
        end
        if "liq" in out_mg[step].ph
            mg_melt_comp[:,step] = out_mg[step].SS_vec[findfirst(==("liq"), out_mg[step].ph)].Comp_wt .* 100
        else
            mg_melt_comp[:,step] = [NaN for _ in 1:11]
        end
        if "liq" in out_og[step].ph
            og_melt_comp[:,step] = out_og[step].SS_vec[findfirst(==("liq"), out_og[step].ph)].Comp_wt .* 100
        else
            og_melt_comp[:,step] = [NaN for _ in 1:11]
        end
        if "liq" in out_og_mp_melt[step].ph
            og_mp_melt_comp[:,step] = out_og_mp_melt[step].SS_vec[findfirst(==("liq"), out_og_mp_melt[step].ph)].Comp_wt .* 100
        else
            og_mp_melt_comp[:,step] = [NaN for _ in 1:11]
        end
        if "liq" in out_og_mg_melt[step].ph
            og_mg_melt_comp[:,step] = out_og_mg_melt[step].SS_vec[findfirst(==("liq"), out_og_mg_melt[step].ph)].Comp_wt .* 100
        else
            og_mg_melt_comp[:,step] = [NaN for _ in 1:11]
        end
    end

    fig_melt_comp = Figure(size = (1600, 2400), fontsize = 40)
    ax_h2O = Axis(fig_melt_comp[1, 1], xlabel="T (°C)", ylabel="H2O wt% in melt")
    xlims!(ax_h2O, minimum(T_array), maximum(T_array))
    lines!(ax_h2O, T_array, mp_melt_comp[11,:], color=:black, label="metapelite")
    lines!(ax_h2O, T_array, mg_melt_comp[11,:], color=:gray, label="metagraywacke")
    lines!(ax_h2O, T_array, og_melt_comp[11,:], color=:pink, label="orthogneiss without fluid transfer")
    lines!(ax_h2O, T_array, og_mp_melt_comp[11,:], color=:purple, label="orthogneiss with fluid transfer from metapelite")
    lines!(ax_h2O, T_array, og_mg_melt_comp[11,:], color=:red, label="orthogneiss with fluid transfer from metagraywacke")
    # axislegend(ax_h2O, position=:lb)
    ax_SiO2 = Axis(fig_melt_comp[1, 2], xlabel="T (°C)", ylabel="SiO2 wt% in melt")
    xlims!(ax_SiO2, minimum(T_array), maximum(T_array))
    lines!(ax_SiO2, T_array, mp_melt_comp[1,:], color=:black, label="metapelite")
    lines!(ax_SiO2, T_array, mg_melt_comp[1,:], color=:gray, label="metagraywacke")
    lines!(ax_SiO2, T_array, og_melt_comp[1,:], color=:pink, label="orthogneiss without fluid transfer")
    lines!(ax_SiO2, T_array, og_mp_melt_comp[1,:], color=:purple, label="orthogneiss with fluid transfer from metapelite")
    lines!(ax_SiO2, T_array, og_mg_melt_comp[1,:], color=:red, label="orthogneiss with fluid transfer from metagraywacke")
    ax_al2O3 = Axis(fig_melt_comp[2, 1], xlabel="T (°C)", ylabel="Al2O3 wt% in melt")
    xlims!(ax_al2O3, minimum(T_array), maximum(T_array))
    lines!(ax_al2O3, T_array, mp_melt_comp[2,:], color=:black, label="metapelite")
    lines!(ax_al2O3, T_array, mg_melt_comp[2,:], color=:gray, label="metagraywacke")
    lines!(ax_al2O3, T_array, og_melt_comp[2,:], color=:pink, label="orthogneiss without fluid transfer")
    lines!(ax_al2O3, T_array, og_mp_melt_comp[2,:], color=:purple, label="orthogneiss with fluid transfer from metapelite")
    lines!(ax_al2O3, T_array, og_mg_melt_comp[2,:], color=:red, label="orthogneiss with fluid transfer from metagraywacke")
    ax_CaO = Axis(fig_melt_comp[2, 2], xlabel="T (°C)", ylabel="CaO wt% in melt")
    xlims!(ax_CaO, minimum(T_array), maximum(T_array))
    lines!(ax_CaO, T_array, mp_melt_comp[3,:], color=:black, label="metapelite")
    lines!(ax_CaO, T_array, mg_melt_comp[3,:], color=:gray, label="metagraywacke")
    lines!(ax_CaO, T_array, og_melt_comp[3,:], color=:pink, label="orthogneiss without fluid transfer")
    lines!(ax_CaO, T_array, og_mp_melt_comp[3,:], color=:purple, label="orthogneiss with fluid transfer from metapelite")
    lines!(ax_CaO, T_array, og_mg_melt_comp[3,:], color=:red, label="orthogneiss with fluid transfer from metagraywacke")
    ax_MgO = Axis(fig_melt_comp[3, 1], xlabel="T (°C)", ylabel="MgO wt% in melt")
    xlims!(ax_MgO, minimum(T_array), maximum(T_array))
    lines!(ax_MgO, T_array, mp_melt_comp[4,:], color=:black, label="metapelite")
    lines!(ax_MgO, T_array, mg_melt_comp[4,:], color=:gray, label="metagraywacke")
    lines!(ax_MgO, T_array, og_melt_comp[4,:], color=:pink, label="orthogneiss without fluid transfer")
    lines!(ax_MgO, T_array, og_mp_melt_comp[4,:], color=:purple, label="orthogneiss with fluid transfer from metapelite")
    lines!(ax_MgO, T_array, og_mg_melt_comp[4,:], color=:red, label="orthogneiss with fluid transfer from metagraywacke")
    ax_FeO = Axis(fig_melt_comp[3, 2], xlabel="T (°C)", ylabel="FeO wt% in melt")
    xlims!(ax_FeO, minimum(T_array), maximum(T_array))
    lines!(ax_FeO, T_array, mp_melt_comp[5,:], color=:black, label="metapelite")
    lines!(ax_FeO, T_array, mg_melt_comp[5,:], color=:gray, label="metagraywacke")
    lines!(ax_FeO, T_array, og_melt_comp[5,:], color=:pink, label="orthogneiss without fluid transfer")
    lines!(ax_FeO, T_array, og_mp_melt_comp[5,:], color=:purple, label="orthogneiss with fluid transfer from metapelite")
    lines!(ax_FeO, T_array, og_mg_melt_comp[5,:], color=:red, label="orthogneiss with fluid transfer from metagraywacke")
    ax_K2O = Axis(fig_melt_comp[4, 1], xlabel="T (°C)", ylabel="K2O wt% in melt")
    xlims!(ax_K2O, minimum(T_array), maximum(T_array))
    lines!(ax_K2O, T_array, mp_melt_comp[6,:], color=:black, label="metapelite")
    lines!(ax_K2O, T_array, mg_melt_comp[6,:], color=:gray, label="metagraywacke")
    lines!(ax_K2O, T_array, og_melt_comp[6,:], color=:pink, label="orthogneiss without fluid transfer")
    lines!(ax_K2O, T_array, og_mp_melt_comp[6,:], color=:purple, label="orthogneiss with fluid transfer from metapelite")
    lines!(ax_K2O, T_array, og_mg_melt_comp[6,:], color=:red, label="orthogneiss with fluid transfer from metagraywacke")
    ax_Na2O = Axis(fig_melt_comp[4, 2], xlabel="T (°C)", ylabel="Na2O wt% in melt")
    xlims!(ax_Na2O, minimum(T_array), maximum(T_array))
    lines!(ax_Na2O, T_array, mp_melt_comp[7,:], color=:black, label="metapelite")
    lines!(ax_Na2O, T_array, mg_melt_comp[7,:], color=:gray, label="metagraywacke")
    lines!(ax_Na2O, T_array, og_melt_comp[7,:], color=:pink, label="orthogneiss without fluid transfer")
    lines!(ax_Na2O, T_array, og_mp_melt_comp[7,:], color=:purple, label="orthogneiss with fluid transfer from metapelite")
    lines!(ax_Na2O, T_array, og_mg_melt_comp[7,:], color=:red, label="orthogneiss with fluid transfer from metagraywacke")
    display(fig_melt_comp)
    return fig_melt_comp
end

function plot_fig1(T_array, H2O_frac_mp_total_cum, H2O_frac_mg_total_cum, T_solidus_mp, T_solidus_mg, T_solidus_og, H2O_released_mp, H2O_released_mg)
    fig1 = Figure(size = (800, 700))
    ax = Axis(fig1[1, 1], xlabel="T (°C)", ylabel="H2O_frac", title="H2O_frac vs T")
    xlims!(ax, minimum(T_array), maximum(T_array))
    ylims!(ax, 0.0, (H2O_frac_mp_total_cum[end,3]*100 + 0.1))
    lines!(ax, T_array, H2O_frac_mp_total_cum[:,3]*100, color=:black, label="metapelite")
    lines!(ax, T_array, H2O_frac_mg_total_cum[:,3]*100, color=:gray, label="metagraywacke")
    # Add vertical dashed lines for solidus temperatures
    vlines!(ax, T_solidus_mp, linestyle=:dash, color=:black, linewidth=1.5, label="T solidus MP")
    vlines!(ax, T_solidus_mg, linestyle=:dash, color=:gray, linewidth=1.5, label="T solidus MG")
    vlines!(ax, T_solidus_og, linestyle=:dash, color=:red, linewidth=1.5, label="T solidus OG")
    # Add text annotation with the water values
    text!(ax, 0.95, 0.05, text="H₂O MP = $(round(H2O_released_mp[3]*100, digits=4)) wt%\nH₂O MG = $(round(H2O_released_mg[3]*100, digits=4)) wt%",
            fontsize=12, align=(:right, :bottom), space=:relative)
    axislegend(ax, position=:lt)
    display(fig1)  # Display the plot
    return fig1
end

function plot_fig2(T_array, melt_frac_og_mp, melt_frac_og_mg, melt_frac_mp, melt_frac_mg, melt_frac_og, T_solidus_mp, T_solidus_mg, T_solidus_og)
    fig2 = Figure(size = (800, 700))
    ax = Axis(fig2[1, 1], xlabel="T (°C)", ylabel="Liq_frac", title="Liq_frac vs T")
    xlims!(ax, minimum(T_array), maximum(T_array))
    ylims!(ax, 0.0, 0.1)
    lines!(ax, T_array, melt_frac_og_mp, color=:red, label="orthogneiss+water-mp")
    lines!(ax, T_array, melt_frac_og_mg, color=:purple, label="orthogneiss+water-mg")
    lines!(ax, T_array, melt_frac_mp, color=:black, label="metapelite")
    lines!(ax, T_array, melt_frac_mg, color=:blue, label="metagraywacke")
    lines!(ax, T_array, melt_frac_og, color=:pink, label="orthogneiss")
    # Add vertical dashed lines for solidus temperatures
    vlines!(ax, T_solidus_mp, linestyle=:dash, color=:black, linewidth=1.5, label="T solidus MP")
    vlines!(ax, T_solidus_mg, linestyle=:dash, color=:blue, linewidth=1.5, label="T solidus MG")
    vlines!(ax, T_solidus_og, linestyle=:dash, color=:red, linewidth=1.5, label="T solidus OG")
    axislegend(ax, position=:lt)
    display(fig2)  # Display the plot
    return fig2
end

function collect_phases(n_steps, out, effect_frac_total)
    phase_names = unique(vcat([out[i].ph for i in 1:n_steps]...))
    phase_names = filter(x -> x != "g", phase_names)  # Exclude garnet
    n_phases = length(phase_names)
    phase_frac = zeros(n_steps, n_phases)

    for i in 1:n_steps
        for (j, phase) in enumerate(out[i].ph)
            if phase != "g"  # Skip garnet as it's fractionated
                phase_idx = findfirst(==(phase), phase_names)
                if !isnothing(phase_idx)
                    # Scale by effective fraction to show what remains after fractionation
                    phase_frac[i, phase_idx] = out[i].ph_frac[j] * effect_frac_total[i,1]
                end
            end
        end
    end

    return phase_names, phase_frac, n_phases
end

function ax_boxplot(ax, T_array, n_steps, phase_frac, g_frac_total_cum, H2O_frac_total_cum, phase_names, n_phases, T_solidus, T_og_solidus, phase_color_map)
    xlims!(ax, minimum(T_array), maximum(T_array))
    ylims!(ax, 0.0, 1.0)

    cumulative = zeros(n_steps)
    for j in 1:n_phases
        lower = copy(cumulative)
        cumulative .+= phase_frac[:, j]
        band!(ax, T_array, lower, cumulative, color=phase_color_map[phase_names[j]], label=phase_names[j])
    end

    # Add garnet on top to show fractionated amount
    lower_g = copy(cumulative)
    cumulative_g = cumulative .+ g_frac_total_cum
    band!(ax, T_array, lower_g, cumulative_g, color=(0.710, 0.173, 0.122), label="g (fractionated)")

    # Add H2O on top to show fractionated water
    lower_h2o = copy(cumulative_g)
    cumulative_h2o = cumulative_g .+ H2O_frac_total_cum
    band!(ax, T_array, lower_h2o, cumulative_h2o, color=:lightblue, label="H2O (fractionated)")

    # Add vertical dashed lines for solidus temperatures
    vlines!(ax, T_solidus, linestyle=:dash, color=:black, linewidth=1.5)
    vlines!(ax, T_og_solidus, linestyle=:dash, color=:red, linewidth=1.5)
    return ax
end

function plot_fig3(T_array, out_mp, out_mg, out_og, out_og_mp_melt, out_og_mg_melt, effect_frac_mp_total, effect_frac_mg_total, effect_frac_og_mp_melt_total, effect_frac_og_mg_melt_total, g_frac_mp_total_cum, g_frac_mg_total_cum, g_frac_og_total_cum, H2O_frac_mp_total_cum, H2O_frac_mg_total_cum, H2O_frac_og_total_cum, melt_frac_mp, melt_frac_mg, melt_frac_og, T_solidus_mp, T_solidus_mg, T_solidus_og, g_frac_og_mp_melt_total_cum, g_frac_og_mg_melt_total_cum, H2O_frac_og_mp_melt_total_cum, H2O_frac_og_mg_melt_total_cum, melt_frac_og_mp, melt_frac_og_mg)
    n_steps = length(T_array)

    phase_names_mp, phase_frac_mp, n_phases_mp = collect_phases(n_steps, out_mp, effect_frac_mp_total)
    phase_names_mg, phase_frac_mg, n_phases_mg = collect_phases(n_steps, out_mg, effect_frac_mg_total)
    phase_names_og_mp, phase_frac_og_mp, n_phases_og_mp = collect_phases(n_steps, out_og_mp_melt, effect_frac_og_mp_melt_total)
    phase_names_og_mg, phase_frac_og_mg, n_phases_og_mg = collect_phases(n_steps, out_og_mg_melt, effect_frac_og_mg_melt_total)

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
    ax3_mp = ax_boxplot(ax3_mp, T_array, n_steps, phase_frac_mp, g_frac_mp_total_cum[:,1], H2O_frac_mp_total_cum[:,1], phase_names_mp, n_phases_mp, T_solidus_mp, T_solidus_og, phase_color_map)

    # Panel 2: Metagraywacke
    ax3_mg = Axis(fig3[1, 2], xlabel="T (°C)", ylabel="Phase Fraction (molar)", title="Mode Box - Metagraywacke")
    ax3_mg = ax_boxplot(ax3_mg, T_array, n_steps, phase_frac_mg, g_frac_mg_total_cum[:,1], H2O_frac_mg_total_cum[:,1], phase_names_mg, n_phases_mg, T_solidus_mg, T_solidus_og, phase_color_map)

    # Panel 3: Orthogneiss + Water from Metapelite
    ax3_og_mp = Axis(fig3[2, 1], xlabel="T (°C)", ylabel="Phase Fraction (molar)", title="Mode Box - Orthogneiss + H₂O from MP")
    ax3_og_mp = ax_boxplot(ax3_og_mp, T_array, n_steps, phase_frac_og_mp, g_frac_og_mp_melt_total_cum[:,1], H2O_frac_og_mp_melt_total_cum[:,1], phase_names_og_mp, n_phases_og_mp, T_solidus_mp, T_solidus_og, phase_color_map)

    # Panel 4: Orthogneiss + Water from Metagraywacke
    ax3_og_mg = Axis(fig3[2, 2], xlabel="T (°C)", ylabel="Phase Fraction (molar)", title="Mode Box - Orthogneiss + H₂O from MG")
    ax3_og_mg = ax_boxplot(ax3_og_mg, T_array, n_steps, phase_frac_og_mg, g_frac_og_mg_melt_total_cum[:,1], H2O_frac_og_mg_melt_total_cum[:,1], phase_names_og_mg, n_phases_og_mg, T_solidus_mg, T_solidus_og, phase_color_map)

    # Create a common legend for all unique phases across all four subplots
    legend_elements = []
    legend_labels = []

    # Add all unique phases with their colors
    for phase in sort(all_phases)
        push!(legend_elements, PolyElement(color=phase_color_map[phase]))
        push!(legend_labels, phase)
    end

    # Add garnet and H2O
    push!(legend_elements, PolyElement(color=(0.710, 0.173, 0.122)))
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
    return fig3
end

end # module
