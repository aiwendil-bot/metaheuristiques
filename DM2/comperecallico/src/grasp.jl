function grasp(C,liaisons_contraintes,liaisons_variables,α,nb_iter,target=0)
	ensemble_z_max = Vector{Int64}(undef, nb_iter)
	compteur::Int64 = 0
	z_max::Int64 = 0
	x = zeros(Int,length(C))
	t1::Float64,t2::Float64,t::Float64,check::Bool = 0.0,0.0,0.0,false
	t1=time()
	while compteur < nb_iter # time() - t1 < tlimit

		x, z = greedy_randomized_construction(C,liaisons_contraintes,liaisons_variables,α)
		x_amelio, z_amelio = simple_descent(C,liaisons_contraintes,liaisons_variables,x,z)
		if z_amelio > target && check == false
			t2=time()
			t = t2-t1
			check = true
		end
		if z_amelio > z_max
			z_max = z_amelio
			x = x_amelio
		end
	compteur += 1
	ensemble_z_max[compteur] = z_max
	end
	return x,z_max,ensemble_z_max,t
end
function grasp_v2(C,liaisons_contraintes,liaisons_variables,α)
	compteur::Int64 = 0
	x_cons,z_cons = greedy_randomized_construction(C,liaisons_contraintes,liaisons_variables,α)
	x_amelio,z_amelio = simple_descent(C,liaisons_contraintes,liaisons_variables,x_cons,z_cons)
	return x_cons,z_cons,x_amelio,z_amelio
end

function greedy_randomized_construction(cost, liaisons_contraintes,liaisons_variables, α)

	dim_matrix = (length(liaisons_contraintes),length(liaisons_variables))
	x_0 = zeros(Int,dim_matrix[2]) #on initie la solution avec des 0

	sous_ensembles_restants = ones(Int,dim_matrix[1]) # si 0 alors acheteur vu
	variables_restantes = ones(Int,dim_matrix[2]) # si 0 alors variable utilisé

	while sum(variables_restantes) > 0 #tant qu'il reste des variables non utilisées
		evaluations = utilities(cost,liaisons_contraintes,liaisons_variables,variables_restantes,sous_ensembles_restants)
		limit::Float64 = minimum(evaluations[1]) + α*(maximum(evaluations[1])-minimum(evaluations[1]))
		rcl = findall(evaluations[1] .>= limit-0.000001)

		rnd = rand(rcl)

		indice_choisi = (evaluations[2])[rnd]

		x_0[indice_choisi] = 1
		variables_restantes[indice_choisi] = 0

		#dans cette boucle, on exclut les sous_ensembles qui sont en conflit avec la solution mise à jour
		for i in liaisons_variables[indice_choisi]
			if sous_ensembles_restants[i] != 0
				for j in liaisons_contraintes[i]
					if j != indice_choisi && variables_restantes[j] == 1
						variables_restantes[j] = 0
					end
				end
			end
			sous_ensembles_restants[i] = 0
		end

	end
	return x_0,dot(x_0,cost)
end
