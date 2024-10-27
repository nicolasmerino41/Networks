PC = "nicol"
include("Set-up.jl")
include("Build_Adjacency_Dictionary.jl")

iberian_interact_NA = Int.(deserialize("Objects\\iberian_interact_NA.jls"))
body_mass_vector = CSV.File("DFs\\body_mass_vector.csv") |> DataFrame
body_mass_vector = body_mass_vector[:, 1]
# Build the adjacency dictionary
adjacency_dict = build_numbered_adjacency_dict(iberian_interact_NA)

fw = Foodweb(adjacency_dict)

m1 = default_model(fw, BodyMass(body_mass_vector))
B0 = rand(256) # Vector of initial biomasses.
t = 2
@time sol = simulate(m1, B0, t)