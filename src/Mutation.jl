
function get_rows_where_state_on(truth_table::Vector{BitVector},row_id::Int)

    filter(x->Bool(digits(x, base=2,pad = Int(log2(length(truth_table))))[row_id]),0:length(truth_table)-1) .+ 1
end

function get_rows_where_state_off(truth_table::Vector{BitVector},row_id::Int)

    filter(x->!Bool(digits(x, base=2,pad = Int(log2(length(truth_table))))[row_id]),0:length(truth_table)-1) .+ 1
end

function get_rows_where_state_on(n_truth_values::Int,row_id::Int)

    filter(x->Bool(digits(x, base=2,pad = Int(log2(n_truth_values)))[row_id]),0:n_truth_values-1) .+ 1
end

function get_rows_where_state_off(n_truth_values::Int,row_id::Int)

    filter(x->!Bool(digits(x, base=2,pad = Int(log2(n_truth_values)))[row_id]),0:n_truth_values-1) .+ 1
end

# Mutation in regulatory regions (add/remove binding sites)

function mutate_regulatory_region(new_network::BooleanNetwork,k_max::Int64,target_node)

    if rand() < 0.5
        # Add a new binding site
        existing_regulators = new_network.all_regulators[target_node]
        new_regulator = rand(1:new_network.n_nodes)

        if new_regulator ∈ existing_regulators
            state_id = findall(existing_regulators .== new_regulator)[1]
            edit_boolean_id = get_rows_where_state_on(length(new_network.boolean_functions[target_node]),state_id)
            new_network.boolean_functions[target_node][edit_boolean_id] .= rand(0:1,length(edit_boolean_id))
                       
        elseif length(existing_regulators) < k_max
            new_network.adjacency_matrix[new_regulator, target_node] = true
            new_network.all_regulators[target_node] = vcat(existing_regulators,[new_regulator]) # append to end. So relevant truth table values are last half of new table 
            # Adjust the truth table for the added regulator
            new_network.boolean_functions[target_node] = vcat(new_network.boolean_functions[target_node],rand(0:1,length(new_network.boolean_functions[target_node])))
        else
            nothing
        end
    else
        # Remove an existing binding site

        existing_regulators = copy(new_network.all_regulators[target_node])

        if !(length(existing_regulators) < 2)
            regulator_to_remove = rand(existing_regulators)
            state_id = findall(existing_regulators .== regulator_to_remove)[1]
            new_network.adjacency_matrix[regulator_to_remove, target_node] = false
            new_network.all_regulators[target_node] = [r for r in existing_regulators if r != regulator_to_remove]
            # Adjust the truth table for the removed regulator
            keep_boolean_id = get_rows_where_state_off(length(new_network.boolean_functions[target_node]),state_id)
            new_network.boolean_functions[target_node] = copy(new_network.boolean_functions[target_node][keep_boolean_id])
        end
    end

    return new_network
end

# Mutation in coding regions (adjust targets based on out-degree and β)

function mutate_coding_region(new_network::BooleanNetwork, β,k_max::Int,target_node)

    βv = rand(β)

    # Determine the number of targets affected based on out-degree and β
    out_degree = sum(new_network.adjacency_matrix[target_node, :])  # Number of targets
    n_targets = ceil(Int,  βv * out_degree)

    for _ in 1:n_targets

        if (out_degree == 0) || rand() < 0.5

            print("1")

            an = sample(1:new_network.n_nodes,1,replace = false)[1]

            # Add a new binding site
            existing_regulators = new_network.all_regulators[an]
            new_regulator = target_node
    
            if new_regulator ∈ existing_regulators
                state_id = findall(existing_regulators .== new_regulator)[1]
                edit_boolean_id = get_rows_where_state_on(length(new_network.boolean_functions[an]),state_id)
                new_network.boolean_functions[an][edit_boolean_id] .= rand(0:1,length(edit_boolean_id))
                        
            elseif length(existing_regulators) < k_max

                new_network.adjacency_matrix[new_regulator, an] = true
                new_network.all_regulators[an] = vcat(new_network.all_regulators[an],[new_regulator]) # append to end. So relevant truth table values are last half of new table 
                # Adjust the truth table for the added regulator
                new_network.boolean_functions[an] = vcat(new_network.boolean_functions[an],rand(0:1,length(new_network.boolean_functions[an])))
            else
                nothing
            end
        else

            print("2")
            print(out_degree)
            # Remove an existing binding site

            an = sample(findall(new_network.adjacency_matrix[target_node, :]),1,replace = false)[1]

            existing_regulators = copy(new_network.all_regulators[an])

            if !(length(existing_regulators) < 2)
                regulator_to_remove = target_node

                print("3")

                state_id = findall(existing_regulators .== regulator_to_remove)[1]

                print("4")

                new_network.adjacency_matrix[regulator_to_remove, an] = false
                print("5")
                new_network.all_regulators[an] = [r for r in existing_regulators if r != regulator_to_remove]
                # Adjust the truth table for the removed regulator
                print("6")
                keep_boolean_id = get_rows_where_state_off(length(new_network.boolean_functions[an]),state_id)
                print("7")
                new_network.boolean_functions[an] = copy(new_network.boolean_functions[an][keep_boolean_id])
            end
        end

    end

    return new_network
end

function duplication_event(network::BooleanNetwork,β,k_max::Int)

    new_network = deepcopy(network)

    node_to_duplicate = rand(1:new_network.n_nodes)

    # Add a new node
    new_network.n_nodes += 1
    new_row = falses(1, new_network.n_nodes - 1)
    new_col = falses(new_network.n_nodes, 1)
    new_network.adjacency_matrix = vcat(new_network.adjacency_matrix,new_row)
    new_network.adjacency_matrix = hcat(new_network.adjacency_matrix,new_col)

    # Copy regulatory inputs from the duplicated node
    new_network.adjacency_matrix[:, new_network.n_nodes] .= new_network.adjacency_matrix[:, node_to_duplicate]

    # Copy targets of the duplicated node
    new_network.adjacency_matrix[new_network.n_nodes, :] .= new_network.adjacency_matrix[node_to_duplicate, :]

    new_network.all_regulators = vcat(new_network.all_regulators, [new_network.all_regulators[node_to_duplicate]])

    for i in 1:new_network.n_nodes
        new_network.all_regulators[i] = findall(new_network.adjacency_matrix[:,i])
    end # need to check this is right

    # Copy logic of the duplicated node
    new_network.boolean_functions = vcat(new_network.boolean_functions,[new_network.boolean_functions[node_to_duplicate]])

    if rand() < 0.5
        new_network = mutate_regulatory_region(new_network,k_max,new_network.n_nodes)
    else
        new_network = mutate_coding_region(new_network, β,k_max,new_network.n_nodes)
    end

    new_network
end

function apply_mutations(network::BooleanNetwork, μ::Float64,β,k_max::Int)

    new_network = deepcopy(network)

    for target_node in 1:new_network.n_nodes
        if rand() < μ
            if rand() < 0.5
                new_network = mutate_regulatory_region(new_network,k_max,target_node)
            else
                new_network = mutate_coding_region(new_network, β,k_max,target_node)
            end
        end
    end

    return new_network
end