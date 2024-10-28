PC = "MM-1"
include("Set-up.jl")
include("Build_Adjacency_Dictionary.jl")
include("HerpsVsBirmmals.jl")
include("Provisional One-click-code.jl")
##### IMPORTING SOME DATA #####
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

##### SET UP SOME OBJECTS #####
herb_carv_vector = []
carv_herb_vector = []
for i in 1:256
    if any(iberian_interact_NA[i, :] .!= 0)
        herb_carv_vector = push!(herb_carv_vector, 0)
        carv_herb_vector = push!(carv_herb_vector, 1)
    elseif all(iberian_interact_NA[i, :] .== 0)
        herb_carv_vector = push!(herb_carv_vector, 1)
        carv_herb_vector = push!(carv_herb_vector, 0)
    end
end
# Identify herbivores and predators based on your data:
herbivores = Bool.(herb_carv_vector)
predators  = Bool.(carv_herb_vector)

# Total number of species
N = length(body_mass_vector)

##### METABOLIC CLASS VECTOR #####
metabolic_class = Vector{Symbol}(undef, N)
# Assign metabolic classes
for i in 1:N
    if herbivores[i]
        metabolic_class[i] = :producer
    elseif !herbivores[i]
        metabolic_class[i] = :invertebrate
    end
end

##### BODY MASSES #####
body_masses = BodyMass(body_mass_vector)

##### BUILDING THE MODEL #####
model = END.Model()
fw = Foodweb(adjacency_dict)
model += fw
model += body_masses

# Add metabolic classes
model += MetabolicClass(metabolic_class)
model += Metabolism(:Miele2019)

# Set mortality rates to 0.01 for all species
model += Mortality(0.01)

# Define logistic growth for herbivores
NPP = 150000000.0   # Total Net Primary Productivity
num_species = length(herbivores)
num_herbivores = length(herbivores[herbivores .== true])
Ki_value = NPP

ri = [herbivores[i] ? rand() : 0.0 for i in 1:num_species]
Ki = [herbivores[i] ? Ki_value : 0.0 for i in 1:num_species]

# Create competition matrix among herbivores (absolute competition)
producers_competition = zeros(num_species, num_species)
for i in axes(iberian_interact_NA, 1), j in axes(iberian_interact_NA, 2)
    if herbivores[i] == 1 && herbivores[j] == 1
        producers_competition[i, j] = 1.0
    end
end

# Add logistic growth component with competition matrix
g_competition = END.LogisticGrowth(
    r = ri,
    K = Ki,
    producers_competition = producers_competition
)
model += g_competition

# Add functional response for predators
model += END.ClassicResponse(h = 2.0)

w = NPP / num_herbivores

# Set initial conditions IF NEEDEED
# B0 = [herbivores[i] ? w : w/10 for i in 1:num_species]

callback = extinction_callback(model, body_mass_vector; verbose = true)
@time d = simulate(model, B0, 10; callback)
###### PLOTTING INDIVIDUAL SPECIES ######
MK.plot(d)

##### PLOTTING BIOMASS, RICHNESS, AND SHANNON #####
time = d.t
Plots.plot(
    time,
    total_biomass(d);
    xlabel = "Time",
    ylabel = "Observable",
    label = "Total biomass",
)
Plots.plot(time, richness(d); label = "Richness")
Plots.plot(time, shannon_diversity(d); label = "Shannon diversity")
MK.plot()