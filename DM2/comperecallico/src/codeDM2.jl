include("../../../libSPP/librarySPP.jl")
#include("/Users/nicolascompere/GitHub/metaheuristiques/DM1/comperecallico/src/codeDM1.jl")

target = "../../../Data"
C, A = loadSPP(string(target,"/","pb_100rnd0100.dat"))
println("pb_100rnd0300.dat")
# CONSTRUCTION GLOUTONNE

using LinearAlgebra
using Random
using Plots
using Distributed

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
#=
x = 1:200
t=@elapsed y=grasppr(C,A,0.6,200,5)
print(t)
println(y[1])
plot(x,y[2])
=#
