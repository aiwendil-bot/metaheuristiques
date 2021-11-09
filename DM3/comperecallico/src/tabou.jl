include("codeDM3.jl")

#tabou avec voisinage 1-1 echange
function tabou(cost::Vector{Int},liaisons_contraintes,liaisons_variables,timelim::Int,taille_tabou::Int)
	nb_variables = size(cost,1)
	n = 1
	xcons, zcons = greedy_construction(cost, liaisons_contraintes, liaisons_variables)
	xamelio, zamelio = deepest_descent(cost,liaisons_contraintes,liaisons_variables,xcons,zcons)
	xmax = copy(xamelio)
	zmax = zamelio
	forbidden_movs = zeros(Int,nb_variables,nb_variables)

	temps = time()
	while (time() < temps + timelim)
		current_sol = 0
		ligne_matrice_tabou= -1
		colonne_matrice_tabou = -1
		for j in 1:nb_variables
			if xamelio[j] == 0
				for kc in liaisons_variables[j]
					for kv in liaisons_contraintes[kc]
						if xamelio[kv]==1 && forbidden_movs[j,kv] < n && j!=kv
							xamelio[j] = 1
							xamelio[kv] = 0
							if  zamelio + cost[j] - cost[kv] >= current_sol && est_admissible(xamelio,liaisons_contraintes, liaisons_variables,j)
								current_sol = zamelio + cost[j] - cost[kv]
								ligne_matrice_tabou = j
								colonne_matrice_tabou = kv
							end
							xamelio[j] = 0
							xamelio[kv] = 1
						end
					end
				end
			end
		end
		n += 1

		if ligne_matrice_tabou!= -1 && colonne_matrice_tabou != -1
			xamelio[ligne_matrice_tabou] = 1
			xamelio[colonne_matrice_tabou] = 0
			zamelio = current_sol
			forbidden_movs[ligne_matrice_tabou,colonne_matrice_tabou] = forbidden_movs[colonne_matrice_tabou,ligne_matrice_tabou] = n + taille_tabou
		end
		if current_sol > zmax
			xmax = copy(xamelio)
			zmax = zamelio
		end
		println(forbidden_movs)
	end
	return xmax,zmax
end
