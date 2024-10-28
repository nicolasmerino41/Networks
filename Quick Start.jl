PC = "MM-1"
include("Set-up.jl")
include("Build_Adjacency_Dictionary.jl")
include("HerpsVsBirmmals.jl")

iberian_interact_NA = Int.(deserialize("DFs\\iberian_interact_NA.jls"))
# Remove self-loops
for i in axes(iberian_interact_NA, 1), j in axes(iberian_interact_NA, 2)
    if i ==j && iberian_interact_NA[i, j] == 1
        iberian_interact_NA[i, j] = 0
    end
end
body_mass_vector = CSV.File("DFs\\body_mass_vector.csv") |> DataFrame
body_mass_vector = body_mass_vector[:, 1] # Turn it into a vector

# Build the adjacency dictionary
adjacency_dict = build_numbered_adjacency_dict(iberian_interact_NA)
# Build adjacency dictionary with names (NOT REALLY GOOD, not even know if it's correct)
adjacency_dict_names = Dict{String, Vector{String}}()
species_names = names(iberian_interact_NA, 1)
for (predator_index, prey_indices) in adjacency_dict
    predator_name = species_names[predator_index]
    prey_names = [species_names[prey_index] for prey_index in prey_indices]
    adjacency_dict_names[predator_name] = prey_names
end

# QUICK EXAMPLE OF THE MODEL BUILDING
fw = Foodweb(adjacency_dict)

m = END.Model(sp, fw, body_masses)
m1 = default_model(fw, BodyMass(body_mass_vector))
m1.species_names = names(iberian_interact_NA)
B0 = rand(256) # Vector of initial biomasses.
t = 100
@time sol = simulate(m1, B0, t)

# PLOT THE SIMULATION RESULTS
time = sol.t
Plots.plot(
    time,
    total_biomass(sol);
    xlabel = "Time",
    ylabel = "Observable",
    label = "Total biomass",
)
Plots.plot!(time, richness(sol); label = "Richness")
Plots.plot!(time, shannon_diversity(sol); label = "Shannon diversity")
