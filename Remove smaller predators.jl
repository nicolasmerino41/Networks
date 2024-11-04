function filter_prey_by_body_mass_ratio(predator_prey_dict::Dict{Int64, Set{Int64}}, body_masses::Vector{Float64}, ratio_threshold::Float64)
    # Dictionary to store the filtered prey sets
    filtered_dict = Dict{Int64, Set{Int64}}()
    
    for (predator, preys) in predator_prey_dict
        # Get the predator's body mass
        predator_mass = body_masses[predator]
        
        # Initialize a new set for filtered prey
        filtered_preys = Set{Int64}()
        
        # Filter prey based on body mass ratio
        for prey in preys
            if (predator_mass / body_masses[prey]) >= ratio_threshold
                push!(filtered_preys, prey)
            end
        end
        
        # Only add to the dictionary if there are any preys left
        if !isempty(filtered_preys)
            filtered_dict[predator] = filtered_preys
        end
    end
    
    return filtered_dict
end
filtered_adj = filter_prey_by_body_mass_ratio(adj, non_zero_body_masses, 1.0)

#############################################################################################
#############################################################################################
######################### USEFUL VARIABLES ###########################################
subcommunity = Int.(subcommunity_raster[cell...])

adj = build_numbered_adjacency_dict(subcommunity)
fw = Foodweb(filtered_adj)

non_zero_body_masses = body_mass_vector[raster_with_abundances_with_B0[cell...].a .> 0]

herb_carv_vector = []
carv_herb_vector = []
for i in 1:length(non_zero_body_masses)
    if any(subcommunity[i, :] .!= 0)
        herb_carv_vector = push!(herb_carv_vector, 0)
        carv_herb_vector = push!(carv_herb_vector, 1)
    elseif all(subcommunity[i, :] .== 0)
        herb_carv_vector = push!(herb_carv_vector, 1)
        carv_herb_vector = push!(carv_herb_vector, 0)
    end
end

# Identify herbivores and predators based on your data:
herbivores = Bool.(herb_carv_vector)
predators  = Bool.(carv_herb_vector)

# Total number of species
num_species = length(non_zero_body_masses)

##### METABOLIC CLASS VECTOR #####
metabolic_class = Vector{Symbol}(undef, num_species)
# Assign metabolic classes
for i in 1:num_species
    if herbivores[i]
        metabolic_class[i] = :producer
    elseif !herbivores[i]
        metabolic_class[i] = :ectotherm
    end
end

homogeneous_metabolic_class = 0.314 .* (non_zero_body_masses.^-0.25)
homogeneous_metabolic_class_with_no_herbivores = 0.314 .* (non_zero_body_masses.^-0.25) 
homogeneous_metabolic_class_with_no_herbivores[herbivores] .= 0.0
# Define logistic growth for herbivores
Ki_value = npp_raster[cell...]   # Total Net Primary Productivity
Ki_value_unity = 1.0
num_herbivores = length(herbivores[herbivores .== true])
num_predators = num_species - num_herbivores

ri = [herbivores[i] ? rand() : 0.0 for i in 1:num_species]
Ki = [herbivores[i] ? Ki_value/num_herbivores : 0.0 for i in 1:num_species]

# Create competition matrix among herbivores (absolute competition)
producers_competition = zeros(num_species, num_species)
for i in axes(subcommunity, 1), j in axes(subcommunity, 2)
    if herbivores[i] == 1 && herbivores[j] == 1 && i != j
        producers_competition[i, j] = 0.0
    elseif herbivores[i] == 1 && herbivores[j] == 1 && i == j
        producers_competition[i, j] = 1.0
    end
end
consumers_efficiency = zeros(num_species, num_species)
for i in axes(subcommunity, 1), j in axes(subcommunity, 2)
    if subcommunity[i, j] != 0 
        consumers_efficiency[i, j] = 1.0
    end
end

# Add logistic growth component with competition matrix
g_competition = END.LogisticGrowth(
    r = ri,
    K = Ki,
    producers_competition = ProducersCompetitionFromDiagonal(1.0)
)
######################################################################################
########################## SET UP MODEL STEP BY STEP #################################
######################################################################################
begin
    model = END.Model()
    model += fw
    model += BodyMass(non_zero_body_masses)
    model += MetabolicClass(:all_ectotherms)
    model += Mortality(0.0)
    model += Metabolism(homogeneous_metabolic_class_with_no_herbivores)
    model += g_competition
    interference_matrix = zeros(num_species, num_species)
    for i in axes(subcommunity, 1), j in axes(subcommunity, 2)
        if subcommunity[i, j] != 0
            interference_matrix[i, j] = 0.7
        end
    end
    model += END.BioenergeticResponse(c = IntraspecificInterference(1.0))
    # model += END.LinearResponse()
    # model += HillExponent(2.0)
end

w = Ki_value / num_herbivores

# Set initial conditions IF NEEDEED
B0 = [herbivores[i] ? w : w/10 for i in 1:num_species]

callback = extinction_callback(model, non_zero_body_masses; verbose = true)

@time a = simulate(model, B0, 10000; callback)

MK.plot(a)
richness(a[end])
sum(a[end])
println(length(model.metabolic_classes[model.metabolic_classes .!= :ectotherm]))