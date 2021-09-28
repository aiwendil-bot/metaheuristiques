include("../../../libSPP/librarySPP.jl")

target = "../../../datatest"
C, A = loadSPP(string(target,"/","pb_100rnd0100.dat"))

# CONSTRUCTION GLOUTONNE

using LinearAlgebra

#calcul les évaluations retourne l'indice de la meilleure
function utilities(cost,matrix,variables_restantes,sous_ensembles_restants)
	evaluations = zeros(Float64,length(cost))
	for i in 1:length(cost)
		somme::Int64 = 0
		for j in 1:length(matrix[:,i])
			if matrix[j,i] == 1
				somme = somme + sous_ensembles_restants[j]
			end
		end
		if somme != 0
			evaluations[i] = cost[i] / somme
		end
	end
	evaluations = evaluations[variables_restantes .== 1] #on ne prend plus en compte que les variables restantes
	indices = findall(variables_restantes .== 1)

	return evaluations,indices
end

#construction gloutonne d'une solution
function greedy_randomized_construction(cost, matrix, α)

	dim_matrix = size(matrix)
	x_0 = zeros(Int,dim_matrix[2]) #on initie la solution avec des 0

	sous_ensembles_restants = ones(Int,dim_matrix[1]) # si 0 alors acheteur vu
	variables_restantes = ones(Int,dim_matrix[2]) # si 0 alors variable utilisé

	while sum(variables_restantes) > 0 #tant qu'il reste des variables non utilisées
		evaluations = utilities(cost,matrix,variables_restantes,sous_ensembles_restants)
		limit::Float64 = minimum(evaluations[1]) + α*(maximum(evaluations[1])-minimum(evaluations[1]))
		rcl = findall(evaluations[1] .>= limit-0.000001)
		println(evaluations)
		println(rcl)
		println("limite :",limit)
		rnd = rand(rcl)


		println(rnd)
		indice_choisi = (evaluations[2])[rnd]



		x_0[indice_choisi] = 1
		variables_restantes[indice_choisi] = 0

		#dans cette boucle, on exclut les sous_ensembles qui sont en conflit avec la solution mise à jour
		for i in 1:length(matrix[:,indice_choisi])
			if matrix[i,indice_choisi] == 1
				if sous_ensembles_restants[i] != 0
					for j in 1:length(matrix[i,:])
						if matrix[i,j] == 1
							if j != indice_choisi && variables_restantes[j] == 1
								variables_restantes[j] = 0
							end
						end
					end
					sous_ensembles_restants[i] = 0
				end
			end
		end
	end
	return x_0,dot(x_0,cost)
end

println(greedy_randomized_construction(C,A,0))
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

function kp_exchange01_simple(cost, matrix, x0, z)
    x = copy(x0)
    max::Int64 = copy(z)
    for i in 1:length(x)
        if x0[i] == 0
            x[i] = 1
            if est_admissible(x,matrix,i)
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

function kp_exchange11_simple(cost, matrix, x0, z)
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


function kp_exchange21_simple(cost, matrix, x0, z)
    x = copy(x0)
    max::Int64 = copy(z)
    for i in 1:length(x)
        if x0[i] == 0
			for j in 1:length(x)
				if x0[j]==1
					for k in (j+1):(length(x)-1)
						if x0[k] == 1 && k != j
							x[i], x[j], x[k] = 1,0,0
							if est_admissible(x,matrix,i)
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

function simple_descent(cost,matrix,x0,zInit)
	x,z = kp_exchange21_simple(cost, matrix, x0,zInit)
	x,z = kp_exchange11_simple(cost, matrix, x, z)
	x,z = kp_exchange01_simple(cost, matrix, x, z)
	return x,z
end

# AMELIORATION PAR PLUS PROFONDE DESCENTE (VOISINAGES : 01-EXCHANGE, 11-EXCHANGE, 21-EXCHANGE)

function kp_exchange01_profond(cost, matrix, x0, z)
    x = copy(x0)
    max::Int64 = copy(z)
    for i in 1:length(x)
        if x0[i] == 0
            x[i] = 1
            if est_admissible(x,matrix,i)
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

function kp_exchange11_profond(cost, matrix, x0, z)
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
                        end
                    end
                    x[i],x[j] = 0,1
                end
            end
        end
    end
    return x0,max
end

function kp_exchange21_profond(cost, matrix, x0, z)
    x = copy(x0)
    max::Int64 = copy(z)
    for i in 1:length(x)
        if x0[i] == 0
			for j in 1:length(x)
				if x0[j]==1
					for k in (j+1):(length(x)-1)
						if x0[k] == 1 && k != j
							x[i], x[j], x[k] = 1,0,0
							if est_admissible(x,matrix,i)
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

function deepest_descent(cost,matrix,x0,zInit)
	x,z = kp_exchange21_profond(cost, matrix, x0, zInit)
	x,z = kp_exchange11_profond(cost, matrix, x, z)
	x,z = kp_exchange01_profond(cost, matrix, x, z)
	return x,z
end
