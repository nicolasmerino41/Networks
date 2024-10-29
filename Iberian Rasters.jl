# good_richness_raster = RS.Raster(joinpath(simbio_path, "Rasters\\good_richness_raster.tif"), lazy = true) 
raster_with_abundances = deserialize(joinpath(simbio_path, "Objects\\raster_with_abundances.jls"))
npp_raster = deserialize(joinpath(simbio_path, "Objects\\npp_raster.jls"))
# MK.heatmap(raster_with_abundances) # THIS DOES NOT WORK FOR SOME REASON
####### RATSER WITH NAMES #######
spain_names = names(iberian_interact_NA, 1)
function create_raster_with_names(raster_with_abundance, spain_names)
    # Initialize a new raster with the same dimensions, storing vectors of strings
    raster_with_names = Raster(
        map(cell -> get_present_species(cell, spain_names), raster_with_abundance.data),
        dims(raster_with_abundance)
    )
    return raster_with_names
end

function get_present_species(cell::MyStructs256, spain_names::Vector{String})
    # Identify indices where abundance > 0
    present_species_indices = findall(>(0), cell.a)
    
    # Map these indices to species names
    present_species_names = spain_names[present_species_indices]
    return present_species_names
end

# Usage example
# Assuming raster_with_abundance is your input Raster with MyStructs and spain_names is your names vector
raster_with_names = create_raster_with_names(raster_with_abundances, spain_names)
##########################################################################################
##########################################################################################
############################## idx vector ################################################
function get_nonzero_coordinates(raster_with_abundance)
    # Initialize an empty vector to store coordinates with non-zero presences
    idx = []

    # Iterate over each cell in the raster and check if it has non-zero presences
    for coord in CartesianIndices(raster_with_abundance.data)
        cell = raster_with_abundance.data[coord]
        if !iszero(cell)
            push!(idx, Tuple(coord))
        end
    end

    return idx
end
idx = get_nonzero_coordinates(raster_with_abundances)
##########################################################################################
##########################################################################################
############################### SUB COMMUNITY RASTER ################################################
function create_subcommunity_raster(raster_with_abundance, idx, iberian_interact_NA)
    # Initialize a new raster with the same dimensions to store submatrices
    subcommunity_raster = Raster(
        fill(Matrix{Float64}(undef, 0, 0), dims(raster_with_abundance)) # Placeholder empty matrices
    )

    # Populate subcommunity raster for each coordinate with non-zero presence
    for coord in idx
        # Get the abundance vector for the current cell
        vector = raster_with_abundance[coord...].a

        # Find indices of species with non-zero abundance
        non_zero_indices = findall(x -> x > 0, vector)

        # Create the submatrix from iberian_interact_NA using non-zero indices
        submatrix = iberian_interact_NA[non_zero_indices, non_zero_indices]

        # Store the submatrix in the new raster
        subcommunity_raster[coord...] = submatrix
    end

    return subcommunity_raster
end

# Assuming raster_with_abundance, idx, and iberian_interact_NA are already defined
subcommunity_raster = create_subcommunity_raster(raster_with_abundances, idx, iberian_interact_NA)