using Alert, Dates, ProgressMeter, CQDBasics

experiment = Experiment("Alex 156")

file_dir = pwd()

atom_number = 150
field_type = "exact"
initial_μₑ = "down"

kᵢ_list = [0, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4]
θₙ_is_fixed_list = [true, false]
branching_condition_list = ["B₀ dominant", "Bₑ dominant"]
average_method_list = [("ABC", 1/16), ("BCA", 1/16), "off"]
initial_μₙ_list = ["HS", "IHS", "Iso"]
sigmoid_field_list = ["off", (0.1, 2e-2)]
BₙBₑ_ratio_list = [(1.3, 1.3), (1.35, 1.35), (1.4, 1.4), (1.45, 1.45), (1.5, 1.5)]
total = length(kᵢ_list) * length(θₙ_is_fixed_list) * length(branching_condition_list) * length(average_method_list) * length(initial_μₙ_list) * length(sigmoid_field_list) * length(BₙBₑ_ratio_list)
p = Progress(total; desc="Running simulations (total $total)...")

for i ∈ eachindex(average_method_list)
    for j ∈ eachindex(branching_condition_list)
        for k ∈ eachindex(θₙ_is_fixed_list)
            for l ∈ eachindex(initial_μₙ_list)
                for m ∈ eachindex(sigmoid_field_list)
                    for n ∈ eachindex(BₙBₑ_ratio_list)
                        for o ∈ eachindex(kᵢ_list)
                            start_time = now()
                            simulation = Simulation(
                                "BE", # type
                                atom_number, # number of atoms
                                field_type, # field
                                initial_μₑ, # initial μₑ
                                initial_μₙ_list[l], # initial μₙ
                                "radau5", # solver
                                θₙ_is_fixed_list[k], # θₙ is fixed
                                branching_condition_list[j], # branching condition
                                "CQD", # BₙBₑ strength
                                BₙBₑ_ratio_list[n], # Bₙ Bₑ ratio
                                kᵢ_list[o], # kᵢ
                                average_method_list[i], # average method
                                sigmoid_field_list[m], # sigmoid field
                                "qm", # R2 comparison
                                "simulation info" # flies to save
                            )
                            raw_data, θₑ_plot, θₙ_plot, θₑθₙ_plot = simulate(experiment, simulation, false)
                            results = Results(experiment, simulation, raw_data, θₑ_plot, θₙ_plot, θₑθₙ_plot)
                            save_results(experiment, simulation, results, start_time, file_dir)
                            # cp(@__FILE__, joinpath(file_dir, Dates.format(start_time, "yyyy-mm-dd_HH-MM-SS-sss"), "Code.jl"))
                            next!(p)
                        end
                    end
                end
            end
        end
    end
end

alert("Simulation finished!")