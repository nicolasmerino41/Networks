###### CONCLUSION: IT LOOKS LIKE THE BODY MASS RATIO IS USUALLY HIGHER THAN 1, WHICH IS GOOD ######
raster_with_names[cell...]
non_zero_body_masses
adj
function average_body_mass_ratios(predator_prey_dict, body_masses)
    # Dictionary to store the average body mass ratio for each predator
    average_ratios = Dict{Int, Float64}()
    
    for (predator, preys) in predator_prey_dict
        # Get the predator's body mass
        predator_mass = body_masses[predator]
        
        # Get the body masses for all prey and compute the mass ratio
        prey_masses = [body_masses[prey] for prey in preys]
        ratios = predator_mass ./ prey_masses
        
        # Calculate and store the average ratio for the current predator
        average_ratios[predator] = mean(ratios)
    end
    
    return average_ratios
end
avg_ratios =average_body_mass_ratios(adj, non_zero_body_masses)

function vector_body_mass_ratios(predator_prey_dict, body_masses)
    # Dictionary to store the average body mass ratio for each predator
    average_ratios = Dict{Int, Vector{Float64}}()
    
    for (predator, preys) in predator_prey_dict
        # Get the predator's body mass
        predator_mass = body_masses[predator]
        
        # Get the body masses for all prey and compute the mass ratio
        prey_masses = [body_masses[prey] for prey in preys]
        ratios = predator_mass ./ prey_masses
        
        # Calculate and store the average ratio for the current predator
        average_ratios[predator] = ratios
    end
    
    return average_ratios
end
vector_ratios = vector_body_mass_ratios(adj, non_zero_body_masses)

# Extract keys and values from the dictionary
predators = collect(keys(avg_ratios))
average_ratios = collect(values(avg_ratios))
avg_log_ratios = log.(average_ratios)
vector_ratios = collect(values(vector_ratios))
log_vector_ratios = deepcopy(vector_ratios)
for i in 1:length(vector_ratios)
    log_vector_ratios[i] = log.(vector_ratios[i])
    println(log_vector_ratios[i])
end
# Plot the data as a scatter plot
PL.scatter(predators, log_vector_ratios,
    xlabel = "Predator ID",
    ylabel = "Average Body Mass Ratio with Prey",
    title = "Predator-Prey Body Mass Ratios",
    legend = false,
    marker = (:circle, 1),
    color = :blue
)
# Add a horizontal red line at y = 1
hline!([0], color = :red, linestyle = :dash, label = "Ratio = 0")
display(plot)


PL.boxplot(predators, vector_ratios)
####################################################################################
####################################################################################
####################### LOOP FOR 20 CELLS SIMULTANEOUSLY ###########################
############################## AVERAGE OF ALL POINTS ###############################
# Prepare the figure with a grid layout for 20 scatter plots (4x5)
using CairoMakie
fig = Figure(resolution = (1200, 800))
grid = fig[1:4, 1:5]  # A 4x5 grid to fit 20 plots
for (i, cell) in enumerate(sample(idx, 20, replace = false))
    subcommunity = Int.(subcommunity_raster[cell...])
    adj = build_numbered_adjacency_dict(subcommunity)
    non_zero_body_masses = body_mass_vector[raster_with_abundances_with_B0[cell...].a .> 0]

    # Calculate the body mass ratios
    avg_ratios = average_body_mass_ratios(adj, non_zero_body_masses)
    predators = collect(keys(avg_ratios))
    average_ratios = collect(values(avg_ratios))
    avg_log_ratios = log.(average_ratios)

    # Sort predators by their body mass
    sorted_predators = sort(predators, by = p -> non_zero_body_masses[p])
    sorted_avg_log_ratios = [avg_log_ratios[predators .== p][1] for p in sorted_predators]
    
    # Calculate row and column in the 4x5 grid
    row = div(i - 1, 5) + 1
    col = mod(i - 1, 5) + 1

    # Create each scatter plot in the grid at the calculated row and column
    ax = Axis(grid[row, col], xlabel = "Body Mass Ordered Predator ID", ylabel = "Log Body Mass Ratio", 
              title = "Predator-Prey Ratios (Sample $i)", titlesize = 10)
    
    # Scatter plot for the current cell, using sorted predator order
    MK.scatter!(ax, sorted_predators, sorted_avg_log_ratios, color = :blue, markersize = 5)
    
    # Add a horizontal red line at y = 0 (log ratio = 0 implies a ratio of 1)
    hlines!(ax, [0], color = :red, linestyle = :dash)
end

fig  # Display the figure
################################################################################
############################ ALL POINTS ########################################
fig = Figure(resolution = (1200, 800))
grid = fig[1:4, 1:5]  # A 4x5 grid to fit 20 plots

for (i, cell) in enumerate(sample(idx, 20, replace = false))
    subcommunity = Int.(subcommunity_raster[cell...])
    adj = build_numbered_adjacency_dict(subcommunity)
    non_zero_body_masses = body_mass_vector[raster_with_abundances_with_B0[cell...].a .> 0]
    
    # Calculate vector ratios for each predator-prey relationship
    vector_ratios = vector_body_mass_ratios(adj, non_zero_body_masses)
    vector_log_ratios = Dict(predator => log.(ratios) for (predator, ratios) in vector_ratios)
    
    # Calculate average ratios and sort predators by body mass
    avg_ratios = average_body_mass_ratios(adj, non_zero_body_masses)
    sorted_predators = sort(collect(keys(avg_ratios)), by = p -> non_zero_body_masses[p])
    
    # Sort the avg_log_ratios according to the sorted_predators
    avg_log_ratios = log.(map(p -> avg_ratios[p], sorted_predators))
    
    # Calculate row and column in the 4x5 grid
    row = div(i - 1, 5) + 1
    col = mod(i - 1, 5) + 1

    # Create each scatter plot in the grid at the calculated row and column
    ax = Axis(grid[row, col], xlabel = "Body Mass Ordered Predator ID", ylabel = "Log Body Mass Ratio", 
              title = "Predator-Prey Ratios (Sample $i)", titlesize = 10)
    
    # Plot lines and scatter points for each predator, sorted by body mass
    for predator in sorted_predators
        ratios = vector_log_ratios[predator]
        MK.lines!(ax, fill(predator, length(ratios)), ratios, color = :blue)
        MK.scatter!(ax, fill(predator, length(ratios)), ratios, color = :blue, markersize = 5)
    end

    # Add a horizontal red line at y = 0 (log ratio = 0 implies a ratio of 1)
    hlines!(ax, [0], color = :red, linestyle = :dash)
end
fig  # Display the figure
################################################################################
################################################################################


   

