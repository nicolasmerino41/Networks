# Parameter options
Ki_values = [Ki_value, Ki_value_unity]

ri_values = [rand, () -> 1.0, () -> 0.1]  # Using functions to generate values


interference_values = [1.0, 0.5, 0.0]

metabolism_values = [1.0, :Miele2019, 0.1]

mortality_values = [0.0, :Miele2019]

producers_competition_functions = [
    ProducersCompetition(1.0),
    ProducersCompetition(0.0),
    ProducersCompetition(0.1),
    ProducersCompetitionFromDiagonal(1.0),
    ProducersCompetitionFromDiagonal(0.0),
    ProducersCompetitionFromDiagonal(0.1)
]

efficiency_values = [:Miele2019, 1.0]

# Array to store results
results = []

# Optionally, collect parameter combinations for reference
param_combinations = []

for Ki_value_option in Ki_values
    for ri_function in ri_values
        for interference_value in interference_values
            for metabolism_value in metabolism_values
                for mortality_value in mortality_values
                    for producers_competition in producers_competition_functions
                        for efficiency_value in efficiency_values
                            # Now, within this loop, set up the model with the current parameters
                            println("start loop something")
                            # Prepare the parameters
                            Ki_value_current = Ki_value_option
                            
                            # Generate ri values for herbivores
                            ri = [herbivores[i] ? ri_function() : 0.0 for i in 1:num_species]
                            
                            # Set Ki values
                            Ki = [herbivores[i] ? Ki_value_current / num_herbivores : 0.0 for i in 1:num_species]
                            
                            # Set up metabolic class
                            if metabolism_value == :Miele2019
                                metabolic_class = :Miele2019
                            else
                                metabolic_class = metabolism_value
                            end
                            
                            # Set up mortality
                            if mortality_value == :Miele2019
                                mortality = Mortality(:Miele2019)
                            else
                                mortality = Mortality(mortality_value)
                            end
                            
                            # Set up efficiency
                            if efficiency_value == :Miele2019
                                efficiency = Efficiency(:Miele2019)
                            else
                                efficiency = Efficiency(efficiency_value)
                            end
                            
                            # Set up interference
                            interference_matrix = zeros(num_species, num_species)
                            for i in axes(subcommunity, 1), j in axes(subcommunity, 2)
                                if subcommunity[i, j] != 0
                                    interference_matrix[i, j] = interference_value
                                end
                            end
                            
                            # Set up bioenergetic response
                            bioenergetic_response = END.BioenergeticResponse(
                                c = IntraspecificInterference(interference_value),
                                e = efficiency
                            )
                            
                            # Set up logistic growth component
                            g_competition = END.LogisticGrowth(
                                r = ri,
                                K = Ki,
                                producers_competition = producers_competition
                            )
                            
                            # Build the model
                            model = END.Model()
                            model += fw
                            model += BodyMass(non_zero_body_masses)
                            model += MetabolicClass(:all_ectotherms)
                            model += mortality
                            model += Metabolism(metabolic_class)
                            model += g_competition
                            model += bioenergetic_response
                            
                            # Initial conditions
                            w = Ki_value_current / num_herbivores
                            B0 = [herbivores[i] ? w : w / 10 for i in 1:num_species]
                            
                            # Run the simulation
                            try
                                a = simulate(model, B0, 1000)
                                final_richness = richness(a[end])
                                final_abundances = a[end]
                                
                                # Save the results
                                push!(results, (
                                    Ki_value = Ki_value_option,
                                    ri_function = ri_function,
                                    interference_value = interference_value,
                                    metabolism_value = metabolism_value,
                                    mortality_value = mortality_value,
                                    producers_competition = producers_competition,
                                    efficiency_value = efficiency_value,
                                    richness = final_richness,
                                    abundances = final_abundances
                                ))
                            catch e
                                # Handle any errors during simulation
                                println("Simulation failed for parameters:")
                                println("Ki_value: $Ki_value_option")
                                println("ri_function: $ri_function")
                                println("interference_value: $interference_value")
                                println("metabolism_value: $metabolism_value")
                                println("mortality_value: $mortality_value")
                                println("producers_competition: $producers_competition")
                                println("efficiency_value: $efficiency_value")
                                println("Error: $e")
                                # Optionally, record the failure
                            end
                            
                        end
                    end
                end
            end
        end
    end
end

