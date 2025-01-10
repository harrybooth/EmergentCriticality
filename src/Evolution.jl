function evolutionary_algorithm(initial_network::BooleanNetwork,n_generations::Int,k_max::Int,μ::Float64,duplication_rate::Int,β,dev_access)
    # Initialize population

    network = deepcopy(initial_network)

    original_attractors = find_all_unique_attractors(network)

    sensitivity = calculate_network_sensitivity(network)

    all_s = [sensitivity]
    all_a = [length(original_attractors)]
    all_f = [alpha_h_fitness(original_attractors)]
    all_n = [network.n_nodes]

    for generation in 1:n_generations

        if generation % 100000 == 0
            println("Generation $generation")
        end

        # Apply duplication-divergence occasionally

        if generation % duplication_rate == 0
            mutated_network = duplication_event(network,β,k_max)
            
            if ! dev_access
                mutated_attractors = find_all_unique_attractors(mutated_network)
            else
                mutated_attractors = attractor_search(mutated_network,original_attractors,0.2,20)
            end

            mutated_attractors_reduced = [Set(Bool.(v[1:end-1]) for v in a) for a in mutated_attractors]

            acc = all([oa ∈ new_attractors_reduced for oa in old_attractors])

            new_attractors = findall([!(oa ∈ old_attractors) for oa in mutated_attractors_reduced])

            if acc
                if length(new_attractors) > 0 
                    αh_mutant = alpha_h_fitness(mutated_attractors[new_attractors])
                    αh_og = alpha_h_fitness(original_attractors)

                    if αh_mutant >= αh_og
                        original_attractors = mutated_attractors
                        network =  mutated_network
                        sensitivity = calculate_network_sensitivity(network)
                    end
                end
            end

        else

            mutated_network = apply_mutations(network,μ,β,k_max)

            if ! dev_access
                mutated_attractors = find_all_unique_attractors(mutated_network)
            else
                mutated_attractors = attractor_search(mutated_network,original_attractors,0.2,20)
            end

            acc =  all([a ∈ mutated_attractors for a in original_attractors])

            new_attractors =  findall([!(a ∈ original_attractors) for a in mutated_attractors])

            if acc
                if length(new_attractors) > 0 
                    αh_mutant = alpha_h_fitness(mutated_attractors[new_attractors])
                    αh_og = alpha_h_fitness(original_attractors)

                    if αh_mutant >= αh_og
                        original_attractors = mutated_attractors
                        network =  mutated_network
                        sensitivity = calculate_network_sensitivity(network)
                    end
                end
            end
        end

        push!(all_s,sensitivity)
        push!(all_a,length(original_attractors))
        push!(all_f,alpha_h_fitness(original_attractors))
        push!(all_n,network.n_nodes)

    end

    return network,all_s,all_a,all_f
end