#=
CQDBasics.jl

This module wraps CQDBase.jl and CQDDataAnalysisBase.jl.
Author: Xukun Lin
Update: 10/24/2024
=#
module CQDBasics

using Reexport

include("CQDBase.jl")
include("CQDDataAnalysisBase.jl")

@reexport using .CQDBase
@reexport using .CQDDataAnalysisBase

include("QuantumMechanics.jl")

export get_quantum_mechanics_results

"""
    see_magnetic_fields(experiment::Experiment, current::Real, magnetic_field_computation_method::String, sigmoid_field::Union{String, Tuple{<:Real, <:Real}})

Visualize the z component of the magnetic fields.
"""
function see_magnetic_fields(experiment::Experiment, current::Real, magnetic_field_computation_method::String, sigmoid_field::Union{String, Tuple{<:Real, <:Real}})
    n = 2^9
    Bz = zeros(n)
    t_list = collect(range(experiment.time_span[1], experiment.time_span[2], n))
    for (i, t) ∈ enumerate(t_list)
        Bz[i] = get_external_magnetic_fields(t, current, experiment, magnetic_field_computation_method, sigmoid_field)[3]
    end
    Plots.default(fontfamily="Computer Modern", tickfont=10, linewidth=1.5, framestyle=:box, legendfontsize=9)
    field_plot = plot()
    plot!(field_plot, t_list * 1e6, Bz, label = "\$B_z\$")
    vspan!(collect(experiment.time_span) * 1e6, label = "", linecolor=:grey, fillcolor=:grey, alpha=0.2)
    xlabel!("\$t\$ [us]"); ylabel!("\$B_z\$ [T]"); title!("\$B_z\$")
    display(field_plot)
end

export see_magnetic_fields

end