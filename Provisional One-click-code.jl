#################################################################################################
########################### LET'S GO ###############################
####################################################################
####################################################################
####################################################################
###################### remove_variables ############################
function remove_variable(var_name::Symbol)
    if isdefined(Main, var_name)
        eval(:(global $var_name = nothing))
        println("Variable $var_name removed from the environment.")
    else
        println("Variable $var_name is not defined.")
    end
end
######################## Size_selection kernel ########################
#######################################################################
function size_selection_kernel(predator_mass, prey_mass, sd, beta)
    intensity = exp(-((log(float(predator_mass)) / (beta * float(prey_mass)))^2.0) / (2.0 * float(sd)^2.0))
    return float(intensity)
end
beta = float(3)

gbif_sizes = CSV.read(joinpath(meta_path, "Lists\\gbif_sizes.csv"), DataFrame)[:, 2:end]

###################### ZERO DIAGONAL ##################################
######################################################################
function zero_out_diagonal!(matrix)
    n = size(matrix, 1)
    for i in 1:n
        matrix[i, i] = 0
    end
    return matrix
end
#################### FILL DIAGONAL ###################################
######################################################################
# If not, it needs to be defined as follows:
function fill_diagonal!(mat, val)
    for i in 1:min(size(mat)...)
        mat[i, i] = val
    end
end
#################### get_neighbors ###################################
######################################################################
# Helper function to get the neighbors
function get_neighbors(matrix, row, col)
    neighbors = []
    rows, cols = size(matrix)
    
    for r in max(row-1, 1):min(row+1, rows)
        for c in max(col-1, 1):min(col+1, cols)
            if (r != row || c != col) && !isnan(matrix[r, c])
                push!(neighbors, matrix[r, c])
            end
        end
    end
    return neighbors
end
#################### when_NA ###################################
######################################################################
function when_NA(array_output)
    for time in 1:length(array_output)
        value = 0.0
        for index in idx
            if any(isinf, array_output[time].state[index].a) || any(isnan, array_output[time].state[index].a)
              println("Time:", time, " Index:", index)
              value += 1
            end
        end
        if value != 0
            println("In time ", time,", ", value, " NA's where generated for the first time")
            return
        end
    end
end
#################### abundance_over_time##############################
######################################################################
function abundance_over_time(abundances)
    # Convert data to a matrix and transpose it
    transposed_abundances = hcat(abundances...)
    # Define timesteps
    timesteps = 1:size(transposed_abundances, 2)
    # Create the figure and axis
    fig = Figure()
    ax = Axis(fig[1, 1], title = "Abundance Over Time", xlabel = "Timesteps", ylabel = "Abundance") 
    # Plot each individual's abundance
    for i in 1:size(transposed_abundances, 1)
        lines!(ax, timesteps, transposed_abundances[i, :], label = "Individual $i") 
    end
    display(fig)
end
#################### count_zeros_ones ##############################
######################################################################
function count_zeros_ones(DA::DimArray{Vector{Float64}, 2}, idx::Vector{CartesianIndex{2}})
    zero_count = 0
    one_count = 0
    
    for index in idx
        current_cell = DA[index]
        zero_count += count(x -> x == 0.0, current_cell)
        one_count += count(x -> x == 1.0, current_cell)
    end
    
    return zero_count, one_count
end

################### map_plot ######################
###################################################
function map_plot(plot::AbstractArray; type = nothing, lambda_DA = nothing, palette = nothing, legend = false, show_grid = true, flip = true, title = nothing, kw...)
    
    if isa(plot, DimArray)
        plot = Matrix(plot)
    end
    
    fig = Figure()
    ax = Axis(fig[1, 1])

    # Set the colormap
    pal = isnothing(palette) ? :thermal : palette

    # Initialize the variable for the plotting object
    plt_obj = nothing

    # Plot according to the specified type
    if type == "image"
        if isnothing(lambda_DA)
            @error("No lambda_DA being used, choose one.")
        end
        plt_obj = image!(ax, plot, lambda_DA; colormap = pal, kw...)
    elseif type == "plot"
        plt_obj = plot!(ax, plot; kw...)
    else
        plt_obj = heatmap!(ax, plot; colormap = pal, kw...)
    end

    # Reverse the y-axis if needed
    ax.yreversed = flip

    # Add legend if requested
    if legend
        non_na_values = filter(!isnan, plot)
        # Ensure that non_na_values is not empty to prevent errors
        if !isempty(non_na_values)
            # Remove colormap and limits from the Colorbar call
            Colorbar(fig[1, 2], plt_obj)
        else
            @warn "No valid data for Colorbar."
        end
    end

    # Add title if provided
    if !isnothing(title)
        ax.title = title
    end

    hidexdecorations!(ax; grid = show_grid)
    hideydecorations!(ax; grid = show_grid)

    fig
end

# Define the custom colormap
custom_palette = cgrad([colorant"black", colorant"blue", colorant"yellow", colorant"green", colorant"red"], [0.0, 0.000000000001, 0.33, 0.66, 1.0]);
# map_plot(DA_richness; palette = custom_palette, rev = true, type = "plot")

####################################################################################################
web = CSV.read(joinpath(meta_path, "Metaweb_data\\TetraEU_pairwise_interactions.csv"), DataFrame)

web = DataFrame(predator = web.sourceTaxonName, prey = web.targetTaxonName)

web.predator = string.(web.predator)
web.prey = string.(web.prey)
web = web[:, [:predator, :prey]]

unique_predators = unique(web.predator)
unique_preys = unique(web.prey)

x = vcat(unique_predators, unique_preys)
unique_species = unique(x)

# Read the CSV file
diets = CSV.File(joinpath(meta_path, "Metaweb_data\\TetraEU_generic_diet.csv")) |> DataFrame
diets = hcat(diets.sourceTaxonName, diets.targetGenericItemName)

Amph = CSV.read(joinpath(meta_path, "Atlas_data/DB_Amphibians_IP.txt"), delim='\t', DataFrame)
Bird = CSV.read(joinpath(meta_path, "Atlas_data/DB_Birds_IP.txt"), delim='\t', DataFrame)
Mamm = CSV.read(joinpath(meta_path, "Atlas_data/DB_Mammals_IP.txt"), delim='\t', DataFrame)
Rept = CSV.read(joinpath(meta_path, "Atlas_data/DB_Reptiles_IP.txt"), delim='\t', DataFrame)
# data = load("Abundance lists ATLAS\\eq_dens_5928cells.RData")
amphibian_names = names(Amph)[2:end]
reptile_names = names(Rept)[2:end]
mammal_names = names(Mamm)[2:end]
bird_names = names(Bird)[2:end]
herps_names = append!(deepcopy(amphibian_names), deepcopy( reptile_names))
birmmals_names = append!(deepcopy(mammal_names), deepcopy(bird_names))
spain_fauna = append!(deepcopy(herps_names), deepcopy(birmmals_names)) 
# bird_names_in_unique_species = length(intersect(bird_names, unique_species_in_web))
# println("The number of bird names found in unique_species_in_web is $bird_names_in_unique_species")
# mammal_names_in_unique_species = length(intersect(mammal_names, unique_species_in_web))
# println("The number of mammal names found in unique_species_in_web is $mammal_names_in_unique_species")
# reptile_names_in_unique_species = length(intersect(reptile_names, unique_species_in_web))
# println("The number of reptile names found in unique_species_in_web is $reptile_names_in_unique_species")
# amphibian_names_in_unique_species = length(intersect(amphibian_names, unique_species_in_web))
# println("The number of amphibian names found in unique_species_in_web is $amphibian_names_in_unique_species")

# Merge the `web` DataFrame with `spain_fauna` using inner join on 'predator' from `web` and 'species' from `spain_fauna`
merged_web = innerjoin(web, DataFrame(species=spain_fauna), on=(:predator => :species))
# Filter the merged DataFrame based on the prey column to include only species found in Spain
merged_web = merged_web[in.(merged_web.prey, Ref(spain_fauna)), :]
# Obtaining the species names that are at least predator/prey
unique_species_in_web = unique(vcat(merged_web.predator, merged_web.prey))
println("There are ", length(unique_species_in_web), " unique species in the food web")

# Initializing an empty matrix with zeros for the Iberian species interaction
n = length(unique_species_in_web)
iberian_interact_matrix = zeros(Int, n, n)
iberian_interact_matrix = NamedArray(iberian_interact_matrix, (unique_species_in_web, unique_species_in_web))
# Ordering the matrix by Apmhibians, Reptiles, Mammals, Birds
spain_names = filter(name -> name in unique_species_in_web, names(hcat(Amph[:, Not(:UTMCODE)], Rept[:, Not(:UTMCODE)], Mamm[:, Not(:UTMCODE)], Bird[:, Not(:UTMCODE)])))
iberian_interact_matrix = iberian_interact_matrix[:, spain_names]
iberian_interact_matrix = iberian_interact_matrix[spain_names, :]

## Creating a mapping from species names to matrix indices
species_to_index = Dict(zip(spain_names, 1:n))
species_names = collect(keys(species_to_index))
#Filling the matrix with 1s where there are predator-prey interactions
for i in 1:nrow(merged_web)
    iberian_interact_matrix[merged_web.predator[i], merged_web.prey[i]] = 1
end
iberian_interact_matrix = float(iberian_interact_matrix)
# Count the amount of 1s in the iberian_interact_matrix
interaction_count = sum(iberian_interact_matrix .== 1)
println("The number of 1s in the iberian_interact_matrix is $interaction_count")

# Turn iberian_interact_matrix into a DataFrame
# Convert the iberian_interact_matrix into a DataFrame with appropriate column and row names
iberian_interact_df = DataFrame(iberian_interact_matrix, species_names)

################### EFFICIENT MATRIX FRAMEWORK #####################
####################################################################
####################################################################
####################################################################
# Load a DataFrame from a serialized file ('.jls' format).
iberian_interact_df = deserialize(joinpath(simbio_path, "Objects\\iberian_interact_df.jls"))
# Convert the DataFrame to a matrix for easier manipulation.
iberian_interact_matrix = iberian_interact_df |> Matrix
# Convert the modified matrix back to a DataFrame, preserving the original column names.
iberian_interact_df = DataFrame(iberian_interact_matrix, names(iberian_interact_df))
# Create a NamedArray from the matrix, using the DataFrame's column names for both dimensions.
iberian_interact_NA = NamedArray(
    iberian_interact_matrix, 
    (names(iberian_interact_df), names(iberian_interact_df)),
    ("Species", "Species")
)
iberian_interact_NA = iberian_interact_NA[spain_names, spain_names]
iberian_interact_NA