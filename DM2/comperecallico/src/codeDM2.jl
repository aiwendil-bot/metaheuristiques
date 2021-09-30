include("../../../libSPP/librarySPP.jl")
#include("/Users/nicolascompere/GitHub/metaheuristiques/DM1/comperecallico/src/codeDM1.jl")

target = "../../../Data"
C, A = loadSPP(string(target,"/","pb_500rnd0500.dat"))
println("pb_100rnd0300.dat")
# CONSTRUCTION GLOUTONNE

using LinearAlgebra
using Random
using Plots

function find_min_key(d)

	minkey = undef
    minval = -1

    for key in keys(d)
        if d[key] < minval
            minkey = key
            minval = d[key]
        end
    end

    return minkey,minval
end

function find_max_key(d)
	maxkey = undef
    maxval = -1

    for key in keys(d)
        if d[key] > maxval
            maxkey = key
            maxval = d[key]
        end
    end

    return maxkey,maxval
end



function grasppr(C,A,α,nb_iter,max_elite)
	ensemble_z_max = Vector{Int64}(undef,nb_iter)
	pool_elite = Dict()
	sizehint!(pool_elite,max_elite)
	x_max = Vector{Int64}(undef,length(C))
	z_max::Int64 = 0
	for i in 1:nb_iter
		x,z_greedy = greedy_randomized_construction(C,A,α)
		x,z = simple_descent(C,A,x,z_greedy)
		if i >= 2
			keys_elite = collect(keys(pool_elite))
			keys_elite = shuffle(keys_elite)
			y = keys_elite[1:rand(1:length(keys_elite))]
			for elite in y
				x_t = dot(C,x) > dot(C,elite) ? x : elite
				x_s = dot(C,x) < dot(C,elite) ? x : elite
				x_p, z_p = path_relinking(C, A, x_s, x_t)
				if length(pool_elite) < max_elite
					pool_elite[x_p] = z_p
				else
					minkey,minval = find_min_key(pool_elite)
					if z_p > minval
						delete!(pool_elite,minkey)
						pool_elite[x_p] = z_p
					end
				end
				if z_p > z_max
					z_max = z_p
					x_max = x_p
				end
			end
		else
			pool_elite[x] = z
		end
		ensemble_z_max[i] = find_max_key(pool_elite)[2]
	end
	return find_max_key(pool_elite),ensemble_z_max
end


function path_relinking(cost,matrix, x_s, x_t)
	symdiff = xor.(x_s, x_t)
	findall_symdiff = findall(x->x==1, symdiff) #indices des éléments qui ne sont pas dans l'intersection
	z_max = max(dot(cost,x_s),dot(cost,x_t))
	z_init = z_max
	x_max = dot(x_s,cost) > dot(x_t,cost) ? x_s : x_t
	x = copy(x_s)
	z_s = dot(cost,x_s)
	i=1
	while length(findall_symdiff) > 0

		valeurs_flips = Dict{Int64, Int64}()
		sizehint!(valeurs_flips, length(findall_symdiff))
		# indice_flip = rand(1:length(findall_symdiff)) # Passer de random à une évaluation
		for i in 1:length(findall_symdiff)
			x[findall_symdiff[i]] = x[findall_symdiff[i]] == 0 ? 1 : 0
			z = x[findall_symdiff[i]] == 1 ? z_s + cost[findall_symdiff[i]] : z_s - cost[findall_symdiff[i]]
			valeurs_flips[findall_symdiff[i]] = z
			x[findall_symdiff[i]] = x[findall_symdiff[i]] == 0 ? 1 : 0
		end
		max_key, max_val = 0, 0
		for i in keys(valeurs_flips)
			if valeurs_flips[i] > max_val
				max_key = i
				max_val = valeurs_flips[i]
			end
		end
		#if max_key != 0
		x[max_key] = x[max_key] == 0 ? 1 : 0
		#end
		symdiff = xor.(x, x_t)
		findall_symdiff = findall(x->x==1, symdiff)
		if isAdmissible(cost,matrix,x)
			if max_val > z_init #solution "prometteuse" => autre critère ? z > 0.9 * z_max etc
				x, max_val = deepest_descent(cost,matrix,x,max_val)
			end

			if max_val > z_max
				z_max = max_val
				x_max = x
			end
		end
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

#println(greedy_randomized_construction(C,A,0.1))

#
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
x = 1:100
y=grasppr(C,A,0.7,100,5)
println(y[1])
plot(x,y[2])
