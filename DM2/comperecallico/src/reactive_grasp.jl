include("codeDM2.jl")
include("grasp.jl")

function reactive_grasp(cost, liaisons_contraintes, liaisons_variables, vector_α, nb_iter, N_α)
	ensemble_z_max = Vector{Int64}(undef, nb_iter)
	nb_α = length(vector_α)
	probas_α = 1 / nb_α .* ones(Float64, nb_α)
	probas_α_cumulees = zeros(Float64, nb_α)
	q = zeros(Float64, nb_α) # notation cours
	compteur_chaque_α = zeros(Int64, nb_α)
	z_best::Int64, z_worst::Int64, z_average = 0, typemax(Int64), Vector{Float64}(undef, nb_α)
	x_max = zeros(Int64, length(cost))
	for i in 1:nb_iter
		s_prob = 0.0
		for j in 1:nb_α
			s_prob += probas_α[j]
			probas_α_cumulees[j] = s_prob
		end
		probas_α_cumulees[nb_α] = 1.0 # arrondis
		if i % N_α == 0
			a = rand()
			indice_α = findfirst(x -> x >= a, probas_α_cumulees)
			x, z_glouton = greedy_randomized_construction(cost, liaisons_contraintes, liaisons_variables, vector_α[indice_α])
			x, z = simple_descent(cost, liaisons_contraintes, liaisons_variables, x, z_glouton)
			if z > z_best
				z_best = z
				x_max = x
			end
			if z < z_worst
				z_worst = z
			end
			z_average[indice_α] = (compteur_chaque_α[indice_α] * z_average[indice_α] + z) / (compteur_chaque_α[indice_α] + 1)
			compteur_chaque_α[indice_α] += 1
		end
			for k in 1:nb_α
				q[k] = (z_average[k] - z_worst) / (z_best - z_worst)
			end
			probas_α = q / sum(q)
		ensemble_z_max[i] = z_best
	end
	return x_max, z_best, probas_α, ensemble_z_max
end
