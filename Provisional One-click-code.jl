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

##################### CONVERT ARGUMENTS #######################
function Makie.convert_arguments(t::Type{<:Makie.Heatmap}, A::AbstractArray{<:MyStructs256, 2})
    scalars = map(mystruct -> mystruct.b, A)
    return Makie.convert_arguments(t, scalars)
end