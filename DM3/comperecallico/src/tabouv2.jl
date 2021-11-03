include("codeDM3.jl")


function tabouv2(cost::Vector{Int},liaisons_contraintes,liaisons_variables,timelim::Int,taille_tabou::Int,taille_echantillonage::Int)
	nb_variables = size(cost,1)
	n = 1
	xcons, zcons = greedy_construction(cost, liaisons_contraintes, liaisons_variables)
	xamelio, zamelio = deepest_descent(cost,liaisons_contraintes,liaisons_variables,xcons,zcons)
	xmax = copy(xamelio)
	init = copy(xmax)
	zmax = zamelio
	forbidden_movs = Vector{Any}(undef,taille_tabou)
	for i in 1:taille_tabou
		forbidden_movs[i] = (0,0,0)
	end
	max_xechantillon, max_zechantillon, max_mov = zeros(length(cost)), 0, (0,0,0)
	temps = time()

	while (time() < temps + timelim)
		# échantillonage
		current = copy(xamelio)
		echantillonage = Vector{Any}(undef,taille_echantillonage) #definir selon taille instance
		i = 1
		while i <= taille_echantillonage
			k,j,l = rand(1:length(cost)),rand(1:length(cost)),rand(1:length(cost))
			if findfirst(x->x==(k,j,l),forbidden_movs) == nothing && current[k] == 0 && current[j] == 1 && current[l] == 1
				echantillonage[i] = (k,j,l)
				i += 1
			end
		end
		max_xechantillon, max_zechantillon, max_mov = zeros(length(cost)), 0, (0,0,0)
		for i in 1:taille_echantillonage
			current[echantillonage[i][1]],current[echantillonage[i][2]],current[echantillonage[i][3]] = 1,0,0
			if dot(current,cost) > max_zechantillon
				max_zechantillon = dot(current,cost)
				max_xechantillon = current
				max_mov = echantillonage[i]
			end
			current[echantillonage[i][1]],current[echantillonage[i][2]],current[echantillonage[i][3]] = 0,1,1
		end
		if max_zechantillon == 0
			println("echantillon non améliorant")
		else
			forbidden_movs[n] = max_mov
			if n < taille_tabou
				n += 1
			else
				n = 1
			end
		end
	end
	return max_zechantillon, max_xechantillon
end
