###### BUILD ADJACENCY DICTIONARY WITH NUMBERS ######
function build_numbered_adjacency_dict(adj_matrix::NamedArray)
    # Initialize an empty dictionary where each predator maps to a set of preys
    adjacency_dict = Dict{Int, Set{Int}}()

    # Iterate over rows (predators) and columns (preys)
    for predator in axes(adj_matrix, 1)
        for prey in axes(adj_matrix, 2)
            # Check if there is a link (1 in the matrix)
            if adj_matrix[predator, prey] == 1
                # Initialize prey set if this predator has no entries yet
                if !haskey(adjacency_dict, predator)
                    adjacency_dict[predator] = Set()
                end
                # Add the prey to the predator's set of preys
                push!(adjacency_dict[predator], prey)
            end
        end
    end

    return adjacency_dict
end

###### BUILD ADJACENCY DICTIONARY WITH NAMES ######
function build_named_adjacency_dict(adj_matrix::NamedArray)
    # Initialize an empty dictionary where each predator maps to a set of prey species names
    adjacency_dict = Dict{String, Set{String}}()

    # Retrieve row and column names (assuming they are the same in this adjacency matrix)
    row_names = map(String, names(adj_matrix, 1))  # Convert predator names to strings
    col_names = map(String, names(adj_matrix, 2))  # Convert prey names to strings

    # Iterate over rows (predators) and columns (preys) by index
    for (predator_idx, predator_name) in enumerate(row_names)
        for (prey_idx, prey_name) in enumerate(col_names)
            # Check if there is a link (1 in the matrix)
            if adj_matrix[predator_idx, prey_idx] == 1
                # Initialize prey set if this predator has no entries yet
                if !haskey(adjacency_dict, predator_name)
                    adjacency_dict[predator_name] = Set()
                end
                # Add the prey to the predator's set of preys
                push!(adjacency_dict[predator_name], prey_name)
            end
        end
    end

    return adjacency_dict
end
