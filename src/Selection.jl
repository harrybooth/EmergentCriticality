function alpha_h(attractor)

    g0 = [mean(.!(state)) for state in attractor]
    g1 = [mean(state) for state in attractor]

    α = mean(1/2 .* (1 .- abs.(g0 .- g1)))

end


function alpha_v(attractor,n_nodes)

    g0 = [mean(.!([state[i] for state in attractor])) for i in 1:n_nodes]
    g1 = [mean([state[i] for state in attractor]) for i in 1:n_nodes]

    α = mean(1/2 .* (1 .- abs.(g0 .- g1)))

end

alpha_h_fitness(all_attractors) = mean([alpha_h(attractor) for attractor in all_attractors])
alpha_v_fitness(all_attractors,n_nodes) = mean([alpha_v(attractor,n_nodes) for attractor in all_attractors])
    