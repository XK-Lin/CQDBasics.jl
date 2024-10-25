using LinearAlgebra, Statistics, Plots, DifferentialEquations, ODEInterfaceDiffEq, WignerD, QuantumOptics

const ħ = 6.62607015e-34/(2π)
const γₑ = -1.76085963023e11
const γₙ = 1.25e7
const a_hfs = 230.8598601e6
const I = 3/2
const S = 1/2

"""
    get_quantum_mechanics_results(experiment::Experiment, initial_state::Union{String, Float64}, magnetic_field_computation_method::String, sigmoid_field::Union{String, Tuple{<:Real, <:Real}}, show_plot::Bool)

Calculate the quantum mechanical predictions of the spin flip probabilities.
"""
function get_quantum_mechanics_results(experiment::Experiment, initial_state::Union{String, Float64}, magnetic_field_computation_method::String, sigmoid_field::Union{String, Tuple{<:Real, <:Real}}, show_plot::Bool)
    βₑ = SpinBasis(Rational(S))
    σx, σy, σz = Matrix.((sigmax(βₑ).data, sigmay(βₑ).data, sigmaz(βₑ).data))
    βₙ = SpinBasis(Rational(I))
    τx, τy, τz = Matrix.((sigmax(βₙ).data, sigmay(βₙ).data, sigmaz(βₙ).data))
    nS, nI = trunc.(Int64, (2 * S + 1, 2 * I + 1))
    Iₑ = Matrix{ComplexF64}(LinearAlgebra.I, nS, nS)
    Iₙ = Matrix{ComplexF64}(LinearAlgebra.I, nI, nI)
    Sₑx, Sₑy, Sₑz = kron.((σx, σy, σz), fill(Iₙ, 3))
    Sₙx, Sₙy, Sₙz = kron.(fill(Iₑ, 3), (τx, τy, τz))
    σ_int = (Sₑx * Sₙx + Sₑy * Sₙy + Sₑz * Sₙz) / 4
    A_hfs = a_hfs * ħ * 2π
    H_int = A_hfs * σ_int
    function get_initial_ρ₀(state::Union{String, Float64})
        if state ∈ ("up", -1/2)
            ρₑ₀ = [0 0; 0 1]
        else
            ρₑ₀ = [1 0; 0 0]
        end
        return kron(ρₑ₀, Iₙ / 4)
    end
    function get_hamiltonian(B)
        Hₑ = -γₑ * ħ / 2 * kron(sum(B .* (σx, σy, σz)), Iₙ)
        Hₙ = -γₙ * ħ / 2 * kron(Iₑ, sum(B .* (τx, τy, τz)))
        return Hₑ + Hₙ + H_int
    end
    function vonNeumann!(du, u, p, t)
        # p = (current, experiment, magnetic_field_computation_method, sigmoid_field)
        # u[:, 1:nS * nI]: real(kron(ρₙ, ρₑ))
        # u[:, nS * nI + 1:2 * nS * nI]: imag(kron(ρₙ, ρₑ))
        B = get_external_magnetic_fields(t, p...)
        H = get_hamiltonian(B)
        ρ = u[:, 1:nS * nI] + 1im * u[:, nS * nI + 1:2 * nS * nI]
        dρ = (H * ρ - ρ * H) / (1im * ħ)
        du .= hcat(real.(dρ), imag.(dρ))
    end
    function get_eigenstates(t, current, experiment, magnetic_field_computation_method, sigmoid_field)
        eigenstates = eigvecs(get_hamiltonian(get_external_magnetic_fields(t, current, experiment, magnetic_field_computation_method, sigmoid_field)))
        return eigenstates ./ norm.(eachcol(eigenstates))
    end
    flip_probabilities = []
    for i ∈ eachindex(experiment.currents)
        v = get_eigenstates(experiment.time_span[1], experiment.currents[i], experiment, magnetic_field_computation_method, sigmoid_field)
        ρ₀ = get_initial_ρ₀(initial_state)
        ρ₀ = v * ρ₀ * v'
        u₀ = hcat(real.(ρ₀), imag.(ρ₀))
        p = (experiment.currents[i], experiment, magnetic_field_computation_method, sigmoid_field)
        prob = ODEProblem(vonNeumann!, u₀, experiment.time_span, p)
        @time solution = solve(prob, radau5(), reltol=1e-6, abstol=1e-6, dtmin=1e-30, force_dtmin=true, maxiters=1e14, saveat=2e-8, dt=1e-30)
        ρend = solution.u[end][:, 1:nS * nI] + 1im * solution.u[end][:, nS * nI + 1:2 * nS * nI]
        v = get_eigenstates(solution.t[end], experiment.currents[i], experiment, magnetic_field_computation_method, sigmoid_field)
        ρend = v' * ρend * v
        diagonal = abs.(diag(ρend))
        if initial_state ∈ ("up", -1/2)
            p_i = sum(diagonal[1:trunc(Int, nS * nI / 2)])
        else
            p_i = sum(diagonal[trunc(Int, nS * nI / 2 + 1):end])
        end
        push!(flip_probabilities, p_i)
    end
    if show_plot
        Plots.default(fontfamily="Computer Modern", tickfont=10, linewidth=1.5, framestyle=:box, legendfontsize=9)
        qm_plot = plot()
        start_index = experiment.currents[1] ≈ 0 ? 2 : 1
        plot!(qm_plot, abs.(experiment.currents[start_index:end]), flip_probabilities[start_index:end], marker=(:circle, 6), label="QM Solution for " * experiment.name, legend=:best, dpi=300, minorgrid=true, xscale=:log10)
        xlabel!("Current [A]"); ylabel!("Flip Probability"); title!("Quantum Mechanics Simulation")
        display(qm_plot)
    end
    return flip_probabilities
end