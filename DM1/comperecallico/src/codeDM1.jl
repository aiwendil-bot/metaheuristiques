# CONSTRUCTION GLOUTONNE

using LinearAlgebra

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

function greedy_construction(cost, matrix)

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
		return x_0,dot(x_0,cost)
end

# AMELIORATION PAR RECHERCHE PROFONDE (VOISINAGES : 01-EXCHANGE, 11-EXCHANGE, 21-EXCHANGE)

function est_admissible(x,matrix,i)
    for compteur in 1:length(matrix[:,i])
		if matrix[compteur,i] == 1
        	if dot(matrix[compteur,:],x) >1
            	return false
        	end
    	end
	end
    return true
end

function kp_exchange01(cost, matrix, x0, z)
    x = copy(x0)
    max::Int64 = copy(z)
    for i in 1:length(x)
        if x0[i] == 0
            x[i] = 1
            if est_admissible(x,matrix,i)
                if z + cost[i] > max
                    max = z + cost[i]
                    x0 = copy(x)
                    println("x")
                end
            end
            x[i] =0
        end
    end
    return x0,max
end

function kp_exchange11(cost, matrix, x0, z)
    x = copy(x0)
    max::Int64 = copy(z)
    for i in 1:length(x)
        if x0[i] == 0
            for j in 1:length(x)
                if x0[j] == 1
                    x[i],x[j] = 1,0
                    if est_admissible(x,matrix,i)
                        if z + cost[i] - cost[j] > max
                            max = z + cost[i] - cost[j]
                            x0 = copy(x)
                            print("x")
                        end
                    end
                    x[i],x[j] = 0,1
                end
            end
        end
    end
    return x0,max
end

function kp_exchange21(cost, matrix, x0, z)
    x = copy(x0)
    max::Int64 = copy(z)
    for i in 1:length(x)
        if x0[i] == 0
			for j in 1:length(x)
				if x0[j]==1
					for k in 1:length(x)
						if x0[k] == 1 && k != j
							x[i], x[j], x[k] = 1,0,0
							if est_admissible(x,matrix,i)
								if z + cost[i] - cost[j] - cost[k] > max
									max = z + cost[i] - cost[j] - cost[k]
									x0 = copy(x)
									print("x")
								end
							end
							x[i],x[j],x[k] = 0,1,1
						end
					end
				end
			end
		end
	end
    return x0,max
end

function deepest_descent(cost,matrix,x0,zInit)
	x,z = kp_exchange21(cost, matrix, x0, zInit)
	x,z = kp_exchange11(cost, matrix, x, z)
	x,z = kp_exchange01(cost, matrix, x, z)
	return x,z
end
