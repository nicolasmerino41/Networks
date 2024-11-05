#### I'M NOT SURE raster_with_abundances_with_B0 IS NECCESSARY ####
raster_with_abundances_with_B0 = deepcopy(raster_with_abundances)
for cell in idx
    abundance_vector = raster_with_abundances[cell...].a
    binary_vector = SVector{256, Float64}([!iszero(abundance_vector[i]) ? 1.0 : 0.0 for i in 1:length(abundance_vector)])
    raster_with_abundances_with_B0[cell...] = MyStructs256(binary_vector)
end

# simulation_matrix = Matrix(undef, 125, 76)
# @Threads.threads for cell in idx[1:8]
cell = idx[2000]    
######################### USEFUL VARIABLES ###########################################
subcommunity = Int.(subcommunity_raster[cell...])

adj = build_numbered_adjacency_dict(subcommunity)
fw = Foodweb(adj)
    
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

ri = [herbivores[i] ? 1.0 : 0.0 for i in 1:num_species]
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
    model += END.BioenergeticResponse(e = consumers_efficiency)
    # model += END.LinearResponse()
    # model += HillExponent(2.0)
end

w = Ki_value / num_herbivores

# Set initial conditions IF NEEDEED
B0 = [herbivores[i] ? w : w/1000 for i in 1:num_species]

callback = extinction_callback(model, non_zero_body_masses; verbose = true)

@time a = simulate(model, B0, 100; callback)

MK.plot(a)
richness(a[end])
sum(a[end])
###########################################################################################
########################## USE DEFAULT MODEL + MODIFICATIONS ##############################
model_default_bio = default_model(
    fw,
    BodyMass(non_zero_body_masses),
    MetabolicClass(metabolic_class),
    Metabolism(:Miele2019),
    END.LogisticGrowth(r = ri, K = Ki, producers_competition = producers_competition),
    Mortality(0.0),
    HillExponent(2.0),
    EfficiencyFromRawValues(1.0),
)

w = Ki_value / num_herbivores
B0 = [herbivores[i] ? w : w/10 for i in 1:num_species]

callback = extinction_callback(model_default_bio, non_zero_body_masses; verbose = true)

@time d = simulate(model, B0, 100; callback)

MK.plot(d)
###########################################################################################
########################### NICHE MODEL FOODWEB ###########################################
fw1 = Foodweb(:niche; S = 75, C = 0.1)
non_zero_body_masses_for_niche = Float64[]
for i in model_niche.metabolic_classes
    if i == :producer
        push!(non_zero_body_masses_for_niche, 0.1)
    elseif i == :ectotherm
        push!(non_zero_body_masses_for_niche, 1.0)
    end
end
begin
    model_niche = END.Model()
    model_niche += fw1
    model_niche += MetabolicClass(:all_ectotherms)
    model_niche += BodyMass(non_zero_body_masses_for_niche)
    model_niche += Mortality(:Miele2019)
    model_niche += MetabolismFromRawValues(homogeneous_metabolic_class)
    model_niche += END.LogisticGrowth()
    interference_matrix = zeros(num_species, num_species)
    for i in axes(subcommunity, 1), j in axes(subcommunity, 2)
        if subcommunity[i, j] != 0
            interference_matrix[i, j] = rand()
        end
    end
    model_niche += END.BioenergeticResponse()
    # model_niche += END.LinearResponse()
    # model_niche += HillExponent(2.0)
end

B0 = [rand() for _ in 1:num_species]

callback = extinction_callback(model_niche, 10e-6; verbose = true)

@time a = simulate(model_niche, B0, 100; callback)

MK.plot(a)
richness(a[end])
total_biomass(a[end])
