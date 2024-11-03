i = 1000
    cell = idx[i]    
    ######################### USEFUL VARIABLES ###########################################
    subcommunity = Int.(subcommunity_raster[cell...])
    subcommunity_NA = NamedArray(Int.(subcommunity_raster[cell...]), (raster_with_names[cell...], raster_with_names[cell...]))
    subcommunity_NA_float = NamedArray(Float64.(subcommunity_raster[cell...]), (raster_with_names[cell...], raster_with_names[cell...]))

    full_IM_NA = turn_adj_into_inter(subcommunity_NA_float, 0.1, 1.0, 1.0, 1.0)
    
    num_species = size(subcommunity_NA)[1]
    include("HerpsVsBirmmals.jl")
    ##### RULES #####
    function GLV(state::MyStructs256, k_DA::MyStructs256)
        return MyStructs256(
            SVector{256, Float64}(
                state.a + (state.a .* (k_DA.a - state.a) + ((full_IM * state.a) .* state.a) + ((full_comp * state.a) .* state.a)) 
            )
        )
    end

    biotic_GLV = Cell{Tuple{:state, :k_DA}, :state}() do data, (state, k_DA), I
        return MyStructs256(SVector{256, Float64}(max.(0.0, GLV(state, k_DA).a)))
    end

    pepe = (;
        state = MyStructs256(SVector{256, Float64}(rand([1.0, 10.0, 20.0, 100.0], num_species)))
    )

    array_output = ResultOutput(
        pepe; tspan = 1:2
    )

    @time resol = sim!(array_output, Ruleset(biotic_GLV))

