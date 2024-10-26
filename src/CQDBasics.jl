"""
This module wraps CQDBase.jl, CQDDataAnalysisBase.jl, and QuantumMechanics.jl.

Author: Xukun Lin

Update: 10/25/2024

Required package (in addition to those required by wrapped modules): "Reexport".
"""
module CQDBasics

using Reexport

include("CQDBase.jl")
include("CQDDataAnalysisBase.jl")

@reexport using .CQDBase
@reexport using .CQDDataAnalysisBase

include("QuantumMechanics.jl")

export get_quantum_mechanics_results

"""
    see_magnetic_fields(plot_range::Tuple{<:Real, <:Real}, current::Real, zₐ::Real, Bᵣ::Vector{<:Real}, magnetic_field_computation_method::String, sigmoid_field::Union{String, Tuple{<:Real, <:Real}})

Visualize the z component of the magnetic fields.
"""
function see_magnetic_fields(plot_range::Tuple{<:Real, <:Real}, current::Real, zₐ::Real, Bᵣ::Vector{<:Real}, magnetic_field_computation_method::String, sigmoid_field::Union{String, Tuple{<:Real, <:Real}})
    n = 2^9
    Bz = zeros(n)
    y_list = collect(range(plot_range[1], plot_range[2], n))
    for (i, y) ∈ enumerate(y_list)
        Bz[i] = get_external_magnetic_fields(y, current, zₐ, Bᵣ, magnetic_field_computation_method, sigmoid_field)[3]
    end
    Plots.default(fontfamily="Computer Modern", tickfont=10, linewidth=1.5, framestyle=:box, legendfontsize=9)
    field_plot = plot()
    plot!(field_plot, y_list * 1e3, Bz, label = "\$B_z\$")
    # vspan!(collect(experiment.time_span) * 1e6, label = "", linecolor=:grey, fillcolor=:grey, alpha=0.2)
    xlabel!("\$y\$ [mm]"); ylabel!("\$B_z\$ [T]"); title!("\$B_z\$")
    display(field_plot)
end

export see_magnetic_fields

end