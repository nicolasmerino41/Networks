grid = [MyStructs256(SVector{256, Float64}(rand([1.0, 10.0, 20.0, 100.0], num_species))) for i in 1:2, j in 1:2]

pepe = (;
    state = grid
)

rule = Cell{:state, :state}() do data, state, I

    B0 = Vector(state.a)
    
    d = simulate(model, B0, 1; callback)
    
    return MyStructs256(SVector{256, Float64}(d[end]))
end

array_output = ResultOutput(
    pepe; tspan = 1:2
)
@time resol = sim!(array_output, Ruleset(rule; proc = ThreadedCPU()))

makie_output = MakieOutput(pepe, tspan = 1:3;
    fps = 10, ruleset = Ruleset(rule)
    ) do (; layout, frame)

    # Setup the keys and titles for each plot
    plot_keys = [:biomass]
    titles = ["Biomass"]

    # Create axes for each plot and customize them
    axes = [Axis(layout[1,1], title = titles[1])]

    # Apply the same color map and limits to all plots, ensure axes are hidden, and set titles
    for (ax, key, title) in zip(axes, plot_keys, titles)
        if key == :biomass
            Makie.heatmap!(ax, frame[:state]; interpolate=false)
        end
        
    end
end