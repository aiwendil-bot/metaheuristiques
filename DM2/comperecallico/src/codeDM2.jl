using LinearAlgebra
using Random
using Distributed

include("grasp.jl")


#contient : fonction d'√©valuation, simple descente, plus profonde descente


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

#calcul les √©valuations retourne l'indice de la meilleure

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
#=

=#
