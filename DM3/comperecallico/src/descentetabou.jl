function istabou(forbidden_sols,sol)::Bool
	return findfirst(x -> x == sol, forbidden_sols) != nothing
end

function kp_exchange01_profond_tabou(cost, liaisons_contraintes,liaisons_variables, x0, z,forbidden_sols)
	x = copy(x0)
	max::Int64 = copy(z)
	for i in 1:length(x)
		if x0[i] == 0
			x[i] = 1
			if !istabou(forbidden_sols,x) &&  est_admissible(x,liaisons_contraintes,liaisons_variables,i)
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

function kp_exchange11_profond_tabou(cost, liaisons_contraintes,liaisons_variables, x0, z,forbidden_sols)
	x = copy(x0)
	max::Int64 = copy(z)
	for i in 1:length(x)
		if x0[i] == 0
			for j in 1:length(x)
				if x0[j] == 1
					x[i],x[j] = 1,0
					if !istabou(forbidden_sols,x) && est_admissible(x,liaisons_contraintes,liaisons_variables,i)
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

function kp_exchange21_profond_tabou(cost, liaisons_contraintes,liaisons_variables, x0, z,forbidden_sols)
	x = copy(x0)
	max::Int64 = copy(z)
	for i in 1:length(x)
		if x0[i] == 0
			for j in 1:length(x)
				if x0[j]==1
					for k in (j+1):(length(x)-1)
						if x0[k] == 1 && k != j
							x[i], x[j], x[k] = 1,0,0
							if !istabou(forbidden_sols,x) && est_admissible(x,liaisons_contraintes,liaisons_variables,i)
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

function deepest_descent_tabou(cost::Vector{Int64},liaisons_contraintes,liaisons_variables,x0::Vector{Int64},zInit::Int64,forbidden_sols::Vector{Vector{Int64}})
	x,z = kp_exchange21_profond_tabou(cost, liaisons_contraintes,liaisons_variables, x0, zInit,forbidden_sols)
	x,z = kp_exchange11_profond_tabou(cost, liaisons_contraintes,liaisons_variables, x, z,forbidden_sols)
	x,z = kp_exchange01_profond_tabou(cost, liaisons_contraintes,liaisons_variables, x, z,forbidden_sols)
	return x,z
end
