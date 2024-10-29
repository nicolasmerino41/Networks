#### I'M NOT SURE raster_with_abundances_with_B0 IS NECCESSARY ####
raster_with_abundances_with_B0 = deepcopy(raster_with_abundances)
for cell in idx
    abundance_vector = raster_with_abundances[cell...].a
    binary_vector = SVector{256, Float64}([!iszero(abundance_vector[i]) ? 1.0 : 0.0 for i in 1:length(abundance_vector)])
    raster_with_abundances_with_B0[cell...] = MyStructs256(binary_vector)
end

simulation_matrix = Matrix(undef, 125, 76)
# @Threads.threads for cell in idx[1:8]
    cell = idx[1]    
    model = END.Model()
    subcommunity = Int.(subcommunity_raster[cell...])
    adj = build_numbered_adjacency_dict(subcommunity)
    fw = Foodweb(adj)
    non_zero_body_masses = body_mass_vector[raster_with_abundances_with_B0[cell...].a .> 0]
    model += fw
    model += BodyMass(non_zero_body_masses)
    model += Mortality(0.01)
    
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
    model += MetabolicClass(metabolic_class)
    model += Metabolism(:Miele2019)

    # Define logistic growth for herbivores
    Ki_value = npp_raster[cell...]   # Total Net Primary Productivity
    # Ki_value = 1.0
    num_herbivores = length(herbivores[herbivores .== true])
    
    ri = [herbivores[i] ? rand() : 0.0 for i in 1:num_species]
    Ki = [herbivores[i] ? Ki_value : 0.0 for i in 1:num_species]

    # Create competition matrix among herbivores (absolute competition)
    producers_competition = zeros(num_species, num_species)
    for i in axes(subcommunity, 1), j in axes(subcommunity, 2)
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
    metabolism_vector = 0.314 .* (non_zero_body_masses.^-0.25)
    # Add functional response for predators
    model += END.ClassicResponse(h = 2.0)
    model += END.LinearResponse()
    model += HillExponent(2.0)

    w = Ki_value / num_herbivores

    # Set initial conditions IF NEEDEED
    B0 = [herbivores[i] ? w : w/10 for i in 1:num_species]
    model_default_classic = default_model(
        fw,
        BodyMass(non_zero_body_masses), 
        MetabolicClass(metabolic_class), 
        Metabolism(:Miele2019),
        hill_exponent = 3.0, 
        ClassicResponse(h = 2.0),
        END.LogisticGrowth(r = ri, K = Ki, producers_competition = producers_competition),
    )
    # Define custom allometry with specific parameters for producers
    custom_allometry = Allometry(
        producer = (a = 0.0138, b = -1/4),        # Set a and b for producers
        invertebrate = (a = 0.314, b = -1/4),     # Default values for invertebrates
        ectotherm = (a = 0.88, b = -1/4)          # Default values for ectotherms
    )
    model_default_bio = default_model(
        fw,
        BodyMass(non_zero_body_masses),
        MetabolicClass(metabolic_class),
        Metabolism(:Miele2019),
        END.LogisticGrowth(r = ri, K = Ki, producers_competition = producers_competition),
        Mortality(0.0),
        HillExponent(2.0)
    )

    callback = extinction_callback(model_default_bio, non_zero_body_masses; verbose = true)
    @time d = simulate(model, B0, 100; callback)
    MK.plot(d)
    simulation_matrix[cell...] = d
    index = findfirst(x -> x == 1.5799795f8, d[end])
    # Assume d is the output of the simulation, with d.t as time points and d.u as solution vectors for each time step
num_species = length(d.u[1])  # Assuming each d.u[i] has 54 species
num_timesteps = length(d.t)

# Initialize a vector of vectors to store time series for each species
species_timeseries = [Float32[] for _ in 1:num_species]
p = deepcopy(d)
# Populate species_timeseries with data from each timestep
for i in 1:num_timesteps
    for j in 1:num_species
        push!(species_timeseries[j], Float32(d.u[i][j]))
        p.u[i][j] = d.u[i][j]
    end
end

# species_timeseries now holds the time series data for each species
# species_timeseries[1] contains the time series for species 1, species_timeseries[2] for species 2, etc.

using WGLMakie

# Assuming `d.t` holds the time points and `species_timeseries` contains time series data for each species
time_points = d.t

# Plot each species' time series
fig = Figure()
ax = fig[1, 1] = Axis(fig, xlabel = "Time (t)", ylabel = "Abundance")

# Plot each species in a loop
for species_idx in 1:length(p)
    abundance_data = species_timeseries[species_idx]
    plot!(ax, time_points, abundance_data, label = "Species $species_idx")
end

# Add a legend to identify each species
axislegend(ax)

# Display the figure
fig 

MK.plot(p)