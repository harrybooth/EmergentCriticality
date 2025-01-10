using Random
using StatsBase
using IterTools
using Distributions
using CairoMakie

bit2number(bit_v) = sum([v*2^(i-1) for (i,v) in enumerate(bit_v)])

# Define a Boolean Network
mutable struct BooleanNetwork
    n_nodes::Int                # Number of nodes
    adjacency_matrix::BitArray  # Adjacency matrix (who regulates whom)
    all_regulators :: Vector{Vector{Int}}
    boolean_functions::Vector{Vector{Bool}}  # Boolean functions for each node
end

# Initialize a random Boolean Network
function initialize_network(n_nodes::Int, k::Int)
    # Adjacency matrix: who regulates whom
    adjacency_matrix = falses(n_nodes, n_nodes)
    all_regulators = Vector{Vector{Int}}(undef, n_nodes)
    for i in 1:n_nodes
        regulators = sample(1:n_nodes, k; replace = false)
        adjacency_matrix[regulators, i] .= true
        all_regulators[i] = sort(regulators)
    end

    # Generate random Boolean functions for each node
    boolean_functions = Vector{Vector{Bool}}(undef, n_nodes)
    
    for i in 1:n_nodes
        n_inputs = sum(adjacency_matrix[:, i])
        boolean_functions[i] = rand(Bool, 2^n_inputs)
    end

    return BooleanNetwork(n_nodes, adjacency_matrix,all_regulators, boolean_functions)
end

# Update a single node based on its Boolean function
function update_node(network::BooleanNetwork, state::BitArray, node::Int)
    
    if isempty(network.all_regulators[node])
        return state[node] # No regulators, keep current state
    end

    # Calculate the input index for the Boolean function
    input_bits = state[network.all_regulators[node]]
    input_index = bit2number(input_bits) + 1
    
    # Apply the Boolean function
    return network.boolean_functions[node][input_index]
end

# Evolve the state of the network for one time step
function step_network(network::BooleanNetwork, state::BitArray)
    next_state = copy(state)
    for i in 1:network.n_nodes
        next_state[i] = update_node(network, state, i)
    end
    return next_state
end

# Run the Boolean network for multiple steps
function run_network(network::BooleanNetwork, initial_state::BitArray, steps::Int)
    state = copy(initial_state)
    states = [copy(state)]
    for _ in 1:steps
        state = step_network(network, state)
        push!(states, copy(state))
    end
    return states
end

function calculate_network_sensitivity(network::BooleanNetwork)
    bias_bf = mean.(network.boolean_functions)
    K_bf = length.(network.all_regulators)

    AvS_bf = 2 .* K_bf .* bias_bf .* (1 .- bias_bf)

    mean(AvS_bf)
end