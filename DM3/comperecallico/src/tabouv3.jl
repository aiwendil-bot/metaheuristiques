include("codeDM3.jl")
include("descentetabou.jl")

function tabouv3(cost::Vector{Int},liaisons_contraintes,liaisons_variables,timelim::Int,taille_tabou::Int)
	nb_variables::Int64 = size(cost,1)
	n::Int64 = 1
	xcons::Vector{Int64}, zcons::Int64 = greedy_construction(cost, liaisons_contraintes, liaisons_variables)
	sol::Vector{Int64}, zsol::Int64 = deepest_descent(cost,liaisons_contraintes,liaisons_variables,xcons,zcons)
	xmax = copy(sol)
	zmax = zsol
	forbidden_sols = Vector{Vector{Int64}}(undef,taille_tabou)
	for i in 1:taille_tabou
		forbidden_sols[i] = zeros(length(cost))
	end
	temps = time()

	while (time() < temps + timelim)


		sol, z_current = deepest_descent_tabou(cost,liaisons_contraintes,liaisons_variables,sol,zsol,forbidden_sols)
		if z_current > zmax
			z_max = z_current
			xmax = sol
		end
		#màj tabou
		forbidden_sols[n] = sol
		if n < taille_tabou
			n += 1
		else
			n = 1
		end
	end
	return xmax, zmax
end
