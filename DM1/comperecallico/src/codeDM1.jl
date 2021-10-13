# CONSTRUCTION GLOUTONNE

using LinearAlgebra

#calcul les évaluations retourne l'indice de la meilleure
function indice_meilleure_evaluation(cost,liaisons_contraintes,liaisons_variables,variables_restantes,sous_ensembles_restants)
	evaluations = zeros(Float64,length(cost))
	for i in 1:length(cost)
		somme::Int64 = 0
		for j in 1:length(liaisons_variables[i])
			somme += sous_ensembles_restants[liaisons_variables[i][j]]
		end
		if somme != 0
			evaluations[i] = cost[i] / somme
		end
	end
	max::Float64,indice_max::Int64 = 0, 0
	evaluations = evaluations .* variables_restantes #on ne prend plus en compte que les variables restantes
	for i in 1:length(evaluations)
		if evaluations[i] > max
			indice_max = i
			max = evaluations[i]
		end
	end
	return indice_max
end

#construction gloutonne d'une solution
function greedy_construction(cost, liaisons_contraintes,liaisons_variables)

	x_0 = zeros(Int,length(cost)) #on initie la solution avec des 0

	sous_ensembles_restants = ones(Int,length(liaisons_contraintes)) # si 0 alors acheteur vu
	variables_restantes = ones(Int,length(liaisons_variables)) # si 0 alors variable utilisé

	while sum(variables_restantes) > 0 #tant qu'il reste des variables non utilisées

		#on calcule les nouvelles évaluations et on prend la meilleure:
		indice_max::Int64 = indice_meilleure_evaluation(cost,liaisons_contraintes,liaisons_variables,variables_restantes,sous_ensembles_restants)
		#le cas indice_max = 0 peut survenir s'il ne reste que des coûts restant à 0. auquel cas on les prend tous
		if indice_max == 0
			liste_variables_restantes = findall(x->x==1,variables_restantes)
			for i in liste_variables_restantes
				x_0[i] = 1
				variables_restantes[i] = 0
			end
		else

			x_0[indice_max] = 1
			variables_restantes[indice_max] = 0

			#dans cette boucle, on exclut les sous_ensembles qui sont en conflit avec la solution mise à jour
			for i in liaisons_variables[indice_max]
				if sous_ensembles_restants[i] != 0
					for j in liaisons_contraintes[i]
						if j != indice_max && variables_restantes[j] == 1
							variables_restantes[j] = 0
						end
					end
				end
				sous_ensembles_restants[i] = 0
			end
		end
	end
	return x_0,dot(x_0,cost)
end

# AMELIORATION PAR SIMPLE DESCENTE (VOISINAGES : 01-EXCHANGE, 11-EXCHANGE, 21-EXCHANGE)

function est_admissible(x,matrix,i)::Bool
    for compteur in 1:length(matrix[:,i])
		if matrix[compteur,i] == 1
        	if dot(matrix[compteur,:],x) > 1
            	return false
        	end
    	end
	end
    return true
end

function utilities(cost,liaisons_contraintes,liaisons_variables,variables_restantes,sous_ensembles_restants)
	evaluation = zeros(Float64,length(cost))

	for i in 1:length(cost)
		somme::Int64 = 0
		for j in 1:length(liaisons_variables[i])
				somme = somme + sous_ensembles_restants[liaisons_variables[i][j]]
		end
		if somme != 0
			evaluation[i] = cost[i] / somme
		end
	end
	evaluation = evaluation[variables_restantes .== 1] #on ne prend plus en compte que les variables restantes
	indices = findall(variables_restantes .== 1)

	return evaluation, indices
end

#construction gloutonne d'une solution


#println(greedy_randomized_construction(C,A,0.1))

# AMELIORATION PAR SIMPLE DESCENTE (VOISINAGES : 01-EXCHANGE, 11-EXCHANGE, 21-EXCHANGE)

function est_admissible(x,liaisons_contraintes,liaisons_variables,i)::Bool
    for compteur in 1:length(liaisons_variables[i])
        if sum(x[liaisons_contraintes[liaisons_variables[i][compteur]]]) > 1
        	return false
    	end
	end
    return true
end

function kp_exchange01_simple(cost, liaisons_contraintes,liaisons_variables, x0, z)
    x = copy(x0)
    max::Int64 = copy(z)
    for i in 1:length(x)
        if x0[i] == 0
            x[i] = 1
            if est_admissible(x,liaisons_contraintes,liaisons_variables,i)
                if z + cost[i] > max
                    max = z + cost[i]
                    x0 = copy(x)
					return x0,max
                end
            end
            x[i] =0
        end
    end
    return x0,max
end

function kp_exchange11_simple(cost, liaisons_contraintes,liaisons_variables, x0, z)
    x = copy(x0)
    max::Int64 = copy(z)
    for i in 1:length(x)
        if x0[i] == 0
            for j in 1:length(x)
                if x0[j] == 1
                    x[i],x[j] = 1,0
                    if est_admissible(x,liaisons_contraintes,liaisons_variables,i)
                        if z + cost[i] - cost[j] > max
                            max = z + cost[i] - cost[j]
                            x0 = copy(x)
							return x0,max
                        end
                    end
                    x[i],x[j] = 0,1
                end
            end
        end
    end
    return x0,max
end


function kp_exchange21_simple(cost, liaisons_contraintes,liaisons_variables, x0, z)
    x = copy(x0)
    max::Int64 = copy(z)
    for i in 1:length(x)
        if x0[i] == 0
			for j in 1:length(x)
				if x0[j]==1
					for k in (j+1):(length(x)-1)
						if x0[k] == 1 && k != j
							x[i], x[j], x[k] = 1,0,0
							if est_admissible(x,liaisons_contraintes,liaisons_variables,i)
								if z + cost[i] - cost[j] - cost[k] > max
									max = z + cost[i] - cost[j] - cost[k]
									x0 = copy(x)
									return x0,max
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

function simple_descent(cost,liaisons_contraintes,liaisons_variables,x0,zInit)
	x,z = kp_exchange21_simple(cost, liaisons_contraintes,liaisons_variables, x0,zInit)
	x,z = kp_exchange11_simple(cost, liaisons_contraintes,liaisons_variables, x, z)
	x,z = kp_exchange01_simple(cost, liaisons_contraintes,liaisons_variables, x, z)
	return x,z
end

# AMELIORATION PAR PLUS PROFONDE DESCENTE (VOISINAGES : 01-EXCHANGE, 11-EXCHANGE, 21-EXCHANGE)

function kp_exchange01_profond(cost, liaisons_contraintes,liaisons_variables, x0, z)
    x = copy(x0)
    max::Int64 = copy(z)
    for i in 1:length(x)
        if x0[i] == 0
            x[i] = 1
            if est_admissible(x,liaisons_contraintes,liaisons_variables,i)
                if z + cost[i] > max
                    max = z + cost[i]
                    x0 = copy(x)
                end
            end
            x[i] =0
        end
    end
    return x0,max
end

function kp_exchange11_profond(cost, liaisons_contraintes,liaisons_variables, x0, z)
    x = copy(x0)
    max::Int64 = copy(z)
    for i in 1:length(x)
        if x0[i] == 0
            for j in 1:length(x)
                if x0[j] == 1
                    x[i],x[j] = 1,0
                    if est_admissible(x,liaisons_contraintes,liaisons_variables,i)
                        if z + cost[i] - cost[j] > max
                            max = z + cost[i] - cost[j]
                            x0 = copy(x)
                        end
                    end
                    x[i],x[j] = 0,1
                end
            end
        end
    end
    return x0,max
end

function kp_exchange21_profond(cost, liaisons_contraintes,liaisons_variables, x0, z)
    x = copy(x0)
    max::Int64 = copy(z)
    for i in 1:length(x)
        if x0[i] == 0
			for j in 1:length(x)
				if x0[j]==1
					for k in (j+1):(length(x)-1)
						if x0[k] == 1 && k != j
							x[i], x[j], x[k] = 1,0,0
							if est_admissible(x,liaisons_contraintes,liaisons_variables,i)
								if z + cost[i] - cost[j] - cost[k] > max
									max = z + cost[i] - cost[j] - cost[k]
									x0 = copy(x)
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

function deepest_descent(cost,liaisons_contraintes,liaisons_variables,x0,zInit)
	x,z = kp_exchange21_profond(cost, liaisons_contraintes,liaisons_variables, x0, zInit)
	x,z = kp_exchange11_profond(cost, liaisons_contraintes,liaisons_variables, x, z)
	x,z = kp_exchange01_profond(cost, liaisons_contraintes,liaisons_variables, x, z)
	return x,z
end
