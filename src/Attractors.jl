# Simulate network dynamics to find attractors
function find_all_attractors(network::BooleanNetwork)
    state_space = IterTools.product(fill(0:1, network.n_nodes)...)
    attractors = Dict{BitArray, Any}()
    for state in state_space
        seen_states = Set{BitArray}()
        current_state = BitArray(state)
    
        while !(current_state in seen_states) & !(current_state in keys(attractors))
            push!(seen_states, copy(current_state))
            current_state = step_network(network, current_state)
        end

        if current_state in keys(attractors)
            attractors[BitArray(state)] = attractors[current_state]
        else
            attractor_state = Set{BitArray}()
            
            while !(current_state ∈ attractor_state)
                push!(attractor_state, copy(current_state))
                current_state = step_network(network, current_state)
            end

            attractors[BitArray(state)] = attractor_state  # limit cycle
        end
    
    end
    return attractors
end

function find_all_unique_attractors(network::BooleanNetwork)
    attractors = find_all_attractors(network)
    unique(values(attractors))
end

function attractor_search(network::BooleanNetwork,start_attractors,bit_flip_p,trials)

    attractors = Dict{BitArray, Any}()

    for sa in start_attractors

        cycle_length = length(sa)

        for t in 1:trials

            if cycle_length == 1
                state = [i for i in sa][1]
            else
                state = [i for i in sa][rand(1:cycle_length)]
            end

            perturbed_state = copy(state)

            n_nodes = length(perturbed_state)

            n_perturbations = Int(ceil(n_nodes*bit_flip_p))

            perturb_id = sample(1:n_nodes,n_perturbations,replace = false)

            perturbed_state[perturb_id] .= .!(perturbed_state[perturb_id])

            if n_nodes != network.n_nodes
                perturbed_state = vcat(perturbed_state,rand([true,false],network.n_nodes - n_nodes))
            end

            seen_states = Set{BitArray}()

            current_state = BitArray(perturbed_state)
        
            while !(current_state in seen_states) & !(current_state in keys(attractors))
                push!(seen_states, copy(current_state))
                current_state = step_network(network, current_state)
            end

            if current_state in keys(attractors)
                attractors[BitArray(perturbed_state)] = attractors[current_state]
            else
                attractor_state = Set{BitArray}()
                
                while !(current_state ∈ attractor_state)
                    push!(attractor_state, copy(current_state))
                    current_state = step_network(network, current_state)
                end

                attractors[BitArray(perturbed_state)] = attractor_state  # limit cycle
            end

        end
    end

    unique(values(attractors))
end