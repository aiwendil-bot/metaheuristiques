include("../../../libSPP/librarySPP.jl")
include("/Users/nicolascompere/GitHub/metaheuristiques/DM1/comperecallico/src/codeDM1.jl")

target = "Data"
C, A = loadSPP(string(target,"/","pb_500rnd1600.dat"))

# CONSTRUCTION GLOUTONNE

using LinearAlgebra

function path_relinking(cost,matrix, x_s, x_t, i_max)
	symdiff = xor.(x_s, x_t)
	findall_symdiff = findall(x->x==1, symdiff) #indices des éléments qui ne sont pas dans l'intersection
	z_max = max(dot(cost,x_s),dot(cost,x_t))
	z_init = z_max
	x_max = dot(x_s,cost) > dot(x_t,cost) ? x_s : x_t
	x = copy(x_s)
	i=1
	println(length(findall_symdiff))
	while length(findall_symdiff) > 0 && i < i_max # ?
		valeurs_flips = Dict{Int64, Int64}()
		sizehint!(valeurs_flips, length(findall_symdiff))
		# indice_flip = rand(1:length(findall_symdiff)) # Passer de random à une évaluation
		for i in 1:length(findall_symdiff)
			x[findall_symdiff[i]] = x[findall_symdiff[i]] == 0 ? 1 : 0
			if est_admissible(cost, matrix, findall_symdiff[i])
				z = dot(cost,x)
				valeurs_flips[findall_symdiff[i]] = z
			end
			x[findall_symdiff[i]] = x[findall_symdiff[i]] == 0 ? 1 : 0
		end
		max_key, max_val = 0, 0
		for i in keys(valeurs_flips)
			if valeurs_flips[i] > max_val
				max_key = i
				max_val = valeurs_flips[i]
			end
		end
		return max_key, max_val # Reprendre ici
		#=if est_admissible(x,matrix,findall_symdiff[indice_flip])
			if z > z_init #solution "prometteuse" => autre critère ? z > 0.9 * z_max etc
				x, z = deepest_descent(cost,matrix,x,z)
			end
			symdiff = xor.(x, x_t)
			findall_symdiff = findall(x->x==1, symdiff)
			if z > z_max
				z_max = z
				x_max = copy(x)
			end
		else
			
			x[findall_symdiff[indice_flip]] = x[findall_symdiff[indice_flip]] == 0 ? 1 : 0
		end
		=#
		i=i+1
	end
	if i == i_max
		println("lien non effectué")
	end
	return x_max,z_max
end

function grasp(C,A,α,nb_iter)
	compteur::Int64 = 0
	z_max::Int64 = 0
	x = zeros(Int,length(C))
	while compteur < nb_iter
		x = greedy_randomized_construction(C,A,α)
		x_amelio = simple_descent(C,A,z,dot(C,x))
		z = dot(x_amelio,C)
		if z > z_max
			z_max = z
			x = x_amelio
		end
	i = i + 1
	end
	return x,z_max
end

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

		rnd = rand(rcl)

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

#println(greedy_randomized_construction(C,A,0.5))
x_1 =[0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1]
x_2 = [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]

println(path_relinking(C,A,x_1,x_2,100000))
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
