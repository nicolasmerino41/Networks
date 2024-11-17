#### I'M NOT SURE raster_with_abundances_with_B0 IS NECCESSARY ####
raster_with_abundances_with_B0 = deepcopy(raster_with_abundances)
for cell in idx
    abundance_vector = raster_with_abundances[cell...].a
    binary_vector = SVector{256, Float64}([!iszero(abundance_vector[i]) ? 1.0 : 0.0 for i in 1:length(abundance_vector)])
    raster_with_abundances_with_B0[cell...] = MyStructs256(binary_vector)
end

# simulation_matrix = Matrix(undef, 125, 76)
# @Threads.threads for cell in idx[1:8]
cell = sample(idx, 1)[1]    
######################### USEFUL VARIABLES ###########################################
IM = Bool.(Int.(subcommunity_raster[cell...]))

# Count herbivores and predators from the interaction matrix
function count_trophic_groups(interaction_matrix)
    num_prey = size(interaction_matrix, 2)
    num_herbivores = count(all(x -> x == 0, interaction_matrix[:, :], dims=2))  # Rows with all zeros
    num_predators = num_prey - num_herbivores  # Remaining rows are predators
    return num_herbivores, num_predators
end

# Herbivores struct definition
mutable struct Herbivores
    m::AbstractFloat        # Mortality rate
    H0::AbstractFloat       # Characteristic density (H_i^0)
    H_init::AbstractFloat   # Initial abundance (H_i(0))
    g::AbstractFloat        # Growth rate (to be calculated)
end

# Constructor for Herbivores
Herbivores(; m::AbstractFloat, H0::AbstractFloat, H_init::AbstractFloat, g::AbstractFloat=0.0) = Herbivores(m, H0, H_init, g)

# Create list of herbivores
function create_herbivores_list(interaction_matrix; 
                                m_mean::AbstractFloat=0.1, m_sd::AbstractFloat=0.02,
                                H0_mean::AbstractFloat, H0_sd::AbstractFloat=H0_mean/10)
    # Identify herbivores based on interaction matrix (columns with all zeros)
    num_herbivores = count(all(x -> x == 0, interaction_matrix[:, :], dims=1))
    herbivores_list = []

    # Generate herbivore parameters
    for i in 1:num_herbivores
        m = max(0.01, rand(Normal(m_mean, m_sd)))         # Mortality rate
        H0 = max(1.0, rand(Normal(H0_mean, H0_sd)))       # Characteristic density
        H_init = H0                                       # Initial abundance
        push!(herbivores_list, Herbivores(m=m, H0=H0, H_init=H_init))
    end
    return herbivores_list
end

# Predator struct definition
mutable struct Predator
    m::AbstractFloat        # Mortality rate
    a::AbstractFloat        # Attack rate
    h::AbstractFloat        # Handling time
    e::AbstractFloat        # Conversion efficiency
    P_init::AbstractFloat   # Initial abundance
end

# Constructor for Predator
Predator(; m::AbstractFloat, a::AbstractFloat, h::AbstractFloat, e::AbstractFloat, P_init::AbstractFloat) = Predator(m, a, h, e, P_init)

# Create list of predators
function create_predator_list(interaction_matrix; 
                               m_mean::AbstractFloat=0.1, m_sd::AbstractFloat=0.02,
                               a_mean::AbstractFloat=0.001, a_sd::AbstractFloat=0.0001,
                               h_mean::AbstractFloat=0.1, h_sd::AbstractFloat=0.01,
                               e_mean::AbstractFloat=0.1, e_sd::AbstractFloat=0.01,
                               P_init_mean::AbstractFloat=5.0, P_init_sd::AbstractFloat=1.0)
    # Identify predators based on interaction matrix (rows with at least one prey)
    num_predators = count(any(x -> x == 1, interaction_matrix, dims=2))
    predator_list = []

    # Generate predator parameters
    for _ in 1:num_predators
        m = max(0.01, rand(Normal(m_mean, m_sd)))         # Mortality rate
        a = max(0.0001, rand(Normal(a_mean, a_sd)))       # Attack rate
        h = max(0.01, rand(Normal(h_mean, h_sd)))         # Handling time
        e = max(0.01, rand(Normal(e_mean, e_sd)))         # Conversion efficiency
        P_init = max(0.1, rand(Normal(P_init_mean, P_init_sd)))  # Initial abundance
        push!(predator_list, Predator(m=m, a=a, h=h, e=e, P_init=P_init))
    end
    return predator_list
end

# Calculate growth rates for herbivores
function calculate_growth_rates(herbivores_list, NPP, mu)
    S_star = length(herbivores_list)
    # Calculate Fi for each herbivore
    F_list = [sp.H0 * sp.m for sp in herbivores_list]
    # Calculate the numerator of the competition term
    competition_numerator = 1 + mu * (S_star - 1)
    # Calculate gi for each herbivore
    for (i, sp) in enumerate(herbivores_list)
        Fi = F_list[i]
        sp.g = sp.m * sqrt((competition_numerator / S_star) * (NPP / Fi))
    end
end

# Ecosystem dynamics function
function ecosystem_dynamics!(du, u, p, t)
    herbivores_list, beta_matrix, predator_list, IM, extinction_threshold = p
    S_star = length(herbivores_list)
    num_predators = length(predator_list)
    H = u[1:S_star]  # Herbivore densities
    P = u[S_star+1:end]  # Predator densities
    du_H = zeros(S_star)
    du_P = zeros(num_predators)

    # Herbivore dynamics
    for i in 1:S_star
        sp = herbivores_list[i]
        # Apply extinction threshold
        if H[i] <= extinction_threshold
            H[i] = 0.0  # Set density to zero if below threshold
            continue     # Skip dynamics for extinct species
        end

        # Compute competition term
        competition = sum(beta_matrix[i, :] .* H) / sp.H0
        # Compute predation term
        predation = 0.0
        for k in 1:num_predators
            if IM[k, i]
                pred = predator_list[k]
                f_ki = (pred.a * H[i]) / (1 + pred.a * pred.h * H[i])  # Functional response
                predation += P[k] * f_ki
            end
        end
        # Compute derivative for herbivores
        du_H[i] = H[i] * sp.m * ((sp.g / sp.m) - 1 - competition) - predation
    end

    # Predator dynamics
    for k in 1:num_predators
        pred = predator_list[k]
        # Apply extinction threshold
        if P[k] <= extinction_threshold
            P[k] = 0.0  # Set density to zero if below threshold
            continue     # Skip dynamics for extinct species
        end

        ingestion = 0.0
        for i in 1:S_star
            if IM[k, i]
                f_ki = (pred.a * H[i]) / (1 + pred.a * pred.h * H[i])  # Functional response
                ingestion += f_ki
            end
        end
        # Compute derivative for predators
        du_P[k] = P[k] * (pred.e * ingestion - pred.m)
    end

    # Assign derivatives
    du[1:S_star] = du_H
    du[S_star+1:end] = du_P
end

# Main simulation function with extinction callback and plotting
function run_real_community_simulation(NPP_raster, cell; extinction_threshold=1.0, legend=false, npp = false, mu = 0.5, m_mean_pred = 0.1, plot = true)
    
    interaction_matrix = Bool.(Int.(subcommunity_raster[cell...]))
    
    NPP = Float64(NPP_raster[cell...])
    num_herbivores, num_predators = count_trophic_groups(interaction_matrix)
    if !isnothing(npp)
        NPP = Float64(npp)
    end
    println("NPP: ", typeof(NPP))
    # mu = 0.1
    H0_mean_aprox = NPP / num_herbivores
    
    herbivores_list = create_herbivores_list(interaction_matrix, H0_mean=H0_mean_aprox)
    predator_list = create_predator_list(interaction_matrix; m_mean = m_mean_pred)
    calculate_growth_rates(herbivores_list, NPP, mu)

    S_star = length(herbivores_list)
    beta_matrix = fill(mu, S_star, S_star)
    for i in 1:S_star
        beta_matrix[i, i] = 1.0
    end

    H_init_values = Float64[sp.H_init for sp in herbivores_list]
    P_init_values = Float64[pred.P_init for pred in predator_list]
    u_init = vcat(H_init_values, P_init_values)

    tspan = (0.0, 200.0)
    p = (herbivores_list, beta_matrix, predator_list, interaction_matrix, extinction_threshold)
    prob = ODEProblem(ecosystem_dynamics!, u_init, tspan, p)

    # Solve the ODE
    sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)

    # Extract time series data
    times = sol.t
    herbivore_data = sol[1:length(herbivores_list), :]  # Herbivore dynamics
    predator_data = sol[length(herbivores_list)+1:end, :]  # Predator dynamics

    # Extinction analysis
    extinct_herbivores = count(herbivore_data[:, end] .<= 1.0)
    extinct_predators = count(predator_data[:, end] .<= 1.0)
    alive_herbivores = num_herbivores - extinct_herbivores
    alive_predators = num_predators - extinct_predators

    # Final biomass
    final_herbivore_biomass = sum(herbivore_data[:, end])
    final_predator_biomass = sum(predator_data[:, end])

    println("Simulation Results:")
    println("Final Herbivore Biomass: ", final_herbivore_biomass)
    println("Final Predator Biomass: ", final_predator_biomass)
    println("$alive_herbivores/$num_herbivores herbivore species were present at the end of the simulation.")
    println("$alive_predators/$num_predators predator species were present at the end of the simulation.")

    if plot
    # Plotting dynamics
    fig = MK.Figure(; size = (600, 500))

    ax = MK.Axis(fig[1, 1], xlabel="Time", ylabel="Density",
                 title="Herbivore and Predator Dynamics Over Time")

    # Plot herbivore dynamics with solid lines
    for i in 1:length(herbivores_list)
        MK.lines!(ax, times, herbivore_data[i, :], label="Herbivore $(i)")
    end

    # Plot predator dynamics with dashed lines
    for k in 1:length(predator_list)
        MK.lines!(ax, times, predator_data[k, :], label="Predator $(k)", linestyle=:dash)
    end

    # Add a legend
    if legend
        MK.axislegend(ax; position=:rt)
    end

    # Display the figure
    display(fig)
    end
end

run_real_community_simulation(npp_raster, cell; plot = true, mu = 0.1, m_mean_pred = 0.1)