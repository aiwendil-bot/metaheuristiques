include("loadSPP.jl")
include("getfname.jl")

using LinearAlgebra

cd("$(homedir())/github/metaheuristiques")

fname = "DM1/Data/pb_100rnd0100.dat"
cost, matrix = loadSPP(fname)

#calcul les évaluations retourne l'indice de la meilleure
function indice_meilleure_evaluation(cost,matrix,subset_restant,cost_restant)
	evaluations = zeros(Float64,length(cost))
	for i in 1:length(cost)
		somme = 0
		for j in 1:length(matrix[:,i])
			if matrix[j,i] == 1
				somme = somme + cost_restant[j]
			end
		end
		if somme != 0
			evaluations[i] = cost[i] / somme
		end
	end
	max,indice_max = 0, 0
	for k in 1:size(subset_restant,1)
		if evaluations[k]*subset_restant[k] > max
			indice_max = k
			max = evaluations[k]
		end
	end
	return indice_max
end

function greedy_intelligent(cost, matrix)

	dim_matrix = size(matrix)
	x_0 = zeros(Int,dim_matrix[2]) #on initie la solution avec des 0

	sous_ensembles_restants = ones(Int,dim_matrix[2]) # si 0 alors sous ensemble vu
	cost_restants = ones(Int,dim_matrix[1]) # si 0 alors coût utilisé

	while sum(sous_ensembles_restants) > 0 #tant qu'il reste des sous ensembles non vus

		#on calcule les nouvelles évaluations et on prend la meilleure:
		indice_max = indice_meilleure_evaluation(cost,matrix,sous_ensembles_restants,cost_restants)

		#le cas indice_max =0 peut survenir s'il ne reste que des sous_ensembles avec des coûts égaux à 0, auquel cas on les prend
		if indice_max == 0
			restant = findall(x->x==1,sous_ensembles_restants)
			for k in restant
				x_0[k] = 1
				sous_ensembles_restants[k] = 0
			end
		else
			x_0[indice_max] = 1
			sous_ensembles_restants[indice_max] = 0
			#dans cette boucle, on exclut les sous_ensembles qui sont en conflit avec la solution mise à jour
			for i in 1:length(matrix[:,indice_max])
				if matrix[i,indice_max] == 1
					if cost_restants[i] != 0
						for k in 1:length(matrix[i,:])
							if matrix[i,k] == 1
								if k != indice_max && sous_ensembles_restants[k] == 1
									x_0[k] = 0
									sous_ensembles_restants[k] = 0
								end
							end
						end
						cost_restants[i] = 0
					end
				end
			end
		end
	end
		return x_0
end

println(greedy_intelligent(cost, matrix))
dot(cost,greedy_intelligent(cost,matrix))
